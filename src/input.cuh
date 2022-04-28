#ifndef __MADNALAB_INPUT__
#define __MADNALAB_INPUT__

std::vector<std::string>& getBasePairTypeComponents(std::string& bpType,
                                                    std::shared_ptr<uammd::System> sys){
    
    static bool init = false;
    static std::map<std::string,std::vector<std::string>> bpTypesList;
    
    if(!init){
        bpTypesList["All"] = {"S","P","A","C","G","T"};
        bpTypesList["S"]   = {"S"};
        init = true;
    }

    if(bpTypesList.count(bpType) == 0){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Requested base pair type (%s) than does not exist",bpType.c_str());
    } 

    return bpTypesList[bpType];

}

struct setInfo{
    int setSize;
    thrust::host_vector<int> set2id;
};

template<class infoType>
std::map<int,std::vector<typename infoType::info>> getInfoList(std::shared_ptr<uammd::System>    sys,
                                                               std::string name,std::string filePath){
    std::ifstream inFile(filePath);

    std::map<int,std::vector<typename infoType::info>> infoList;
    std::string line;

    while(std::getline(inFile,line)){
        std::stringstream ss(line);
        
        typename infoType::info infoBuffer = infoType::getInfo(ss);

        infoList[infoBuffer.simId].push_back(infoBuffer);
    }

    return infoList;
}

template<class infoType>
setInfo getSet2id(std::shared_ptr<uammd::System>       sys,
                  std::shared_ptr<uammd::ParticleData>  pd,
                  std::map<int,std::shared_ptr<uammd::ParticleGroup>>& simGroups,
                  std::shared_ptr<ff::Topology> top,
                  std::map<int,std::vector<typename infoType::info>>& infoList,
                  int infoIndex,
                  int bpIndex,
                  std::string model){

    int setSize=-1;
    thrust::host_vector<int> set2id;

    for(auto& info : infoList){

        if(info.second.size() > infoIndex){
            std::vector<int> bp_ids =  basePair2id(pd,top,sys,
                                                   simGroups[info.second[infoIndex].simId],
                                                   getBasePairTypeComponents(info.second[infoIndex].type,sys),
                                                   info.second[infoIndex].bp[bpIndex],
                                                   model);

            if(bp_ids.size() == setSize or setSize < 0){
                setSize = bp_ids.size();
            } else {
                sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                                  "Set size has to be equal for all base pairs");
            }
            
            set2id.insert(set2id.end(),bp_ids.begin(),bp_ids.end());
        }
    }

    setInfo stI;

    stI.setSize = setSize;
    stI.set2id  = set2id;

    return stI;
}

#endif
