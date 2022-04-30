#ifndef __BASE_PAIR__
#define __BASE_PAIR__

int2 basePair2res(int basePairPosition,
                  int nBasis, 
                  std::string model,
                  std::shared_ptr<uammd::System> sys){
    
    if(nBasis%2!=0){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                "The total number of basis should be even, but the value is: %i",
                nBasis);
    }

    int nBasisPairs = nBasis/2;

    if(basePairPosition == 0){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid pair base position, it must be different than 0");
    }

    if(abs(basePairPosition) > nBasisPairs ){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid pair base position,"
                                          "its absolute value can not be larger "
                                          "than the total number of base pairs (%i)," 
                                          "but the value is: %i",
                                           nBasisPairs,basePairPosition);
    }

    int2 res;

    if(model == "MADna" or model == "MADnaFast"){
            
        if(basePairPosition > 0){
            res.x = basePairPosition-1;
            res.y = nBasis-1-res.x;
        } else {
            res.x = nBasisPairs+basePairPosition;
            res.y = nBasis-1-res.x;
        }

    } else if (model == "WLC"){
        
        if(basePairPosition > 0){
            res.x = basePairPosition-1;
            res.y = res.x;
        } else {
            res.x = nBasisPairs+basePairPosition;
            res.y = res.x;
        }

    } else {
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid selected model: %s",model.c_str());
    }

    return res;
}


std::vector<int> basePair2id(std::shared_ptr<uammd::ParticleData> pd,
                             std::shared_ptr<ff::Topology> top,
                             std::shared_ptr<uammd::System> sys,
                             std::shared_ptr<uammd::ParticleGroup> simIdGroup,
                             std::vector<std::string>& types,
                             int basePairPosition,
                             std::string model){

    std::vector<int> idSet;

    auto pos   = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
        
    auto id    = pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
    auto res   = pd->getResId(uammd::access::location::cpu,uammd::access::mode::read);
        
    auto groupIndex = simIdGroup->getIndexIterator(uammd::access::location::cpu);
        
    //Compute nBasis
    int nBasis=1;
    int basis_prev = res[groupIndex[0]];
    for(int i=0;i<simIdGroup->getNumberParticles();i++){
        int index = groupIndex[i];
        
        int basis = res[index];

        if(basis != basis_prev){
            basis_prev=basis;
            nBasis++;
        }
    }

    int2 b1;
    if(model == "MADna" or model == "MADnaFast"){
        b1 = basePair2res(basePairPosition,
                          nBasis, 
                          model,
                          sys);
    } else if (model == "WLC"){
        b1 = basePair2res(basePairPosition,
                          nBasis*2, //Fow WLC nBasis == nBasisPairs
                          model,
                          sys);
    } else {
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid selected model: %s",model.c_str());
    }

    int b11 = b1.x;
    int b12 = b1.y;

    for(int i=0;i<simIdGroup->getNumberParticles();i++){
        
        int index = groupIndex[i];
    
        int type = pos[index].w;

        int pid   = id[index];
        int basis = res[index];

        std::string typeName = top->getTypes()->getTypeParameters(type).name;

        if(std::find(types.begin(),types.end(),typeName) != types.end()){
            if        (basis == b11){
                idSet.push_back(pid);
            } else if (basis == b12){
                idSet.push_back(pid);
            }                 
        }
    }

    return idSet;
}

std::vector<int> basePair2id(std::shared_ptr<uammd::ParticleData> pd,
                             std::shared_ptr<ff::Topology> top,
                             std::shared_ptr<uammd::System> sys,
                             std::shared_ptr<uammd::ParticleGroup> simIdGroup,
                             std::vector<std::string>& types,
                             std::vector<int>& basePairPosition,
                             std::string model){
    
    std::vector<int> idSet;

    for(int bpp : basePairPosition){
        auto idSetBuffer = basePair2id(pd,top,sys,
                                       simIdGroup,
                                       types,bpp,
                                       model);

        idSet.insert(idSet.end(),idSetBuffer.begin(),idSetBuffer.end());
    }
    
    return idSet;
}

#endif
