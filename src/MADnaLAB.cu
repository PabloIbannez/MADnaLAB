#include <UAMMDstructured.cuh>

using namespace uammd::structured;

using ff  = forceField::EMPTY_KCALMOL_A;

using ffMADna = forceField::MADna;
using ffWLC   = forceField::KratkyPorodModel::KratkyPorodModel<forceField::ForceFieldBase<ff::Units,Types::BASIC>>;

using SIM = Simulation<ff,SteepestDescent,LangevinNVT::BBK>;

#include "basePair.cuh"
#include "interactors.cuh"
#include "measures.cuh"
#include "input.cuh"

int main(int argc, char** argv){

    auto sys = std::make_shared<uammd::System>();
    
    uammd::InputFile in(argv[1]);

    bool loadFromBackup = false;
    if(in.getOption("loadFromBackup",uammd::InputFile::Optional)){
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Loading from backup...");
        loadFromBackup = true;
    }
    std::map<std::string,std::string> backupInfo;
    
    ullint seed;
    if(in.getOption("seed",uammd::InputFile::Optional)){
        seed = std::atoi(in.getOption("seed",uammd::InputFile::Optional).str().c_str()); 
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Reading seed from input file, current seed:%lli",seed);
    } else {
        seed = 0xf31337Bada55D00dULL^time(NULL);
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Generating seed using time, current seed:%lli",seed);
    }
    sys->rng().setSeed(seed);

    std::string simulationSetFolder = in.getOption("simulationSetFolder",uammd::InputFile::Required).str();
    sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                     "Simulation folder: %s",simulationSetFolder.c_str());

    std::shared_ptr<SIM> sim = std::make_shared<SIM>(sys,in);

    if(loadFromBackup){
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Loading particle data from backup...");
        
        sim->loadParticleBuffer(simulationSetFolder+"/backup.back");
        sim->updateParticleData(true);
    
        int initStep;
        in.getOption("initStep",uammd::InputFile::Required)>>initStep;

        sim->setStep(initStep);
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Starting from step:%i",initStep);
    }

    auto top = sim->getTopology();
    auto pd  = sim->getParticleData();
    auto pg  = sim->getParticleGroup();
    
    std::string model = in.getOption("modelName",uammd::InputFile::Required).str();
    {
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Selected model: %s",model.c_str());

               if(model == "MADna"){
            std::shared_ptr<ffMADna> madna = std::make_shared<ffMADna>(sys,pd,pg,in);
            sim->addInteractor(madna);
        } else if(model == "WLC"){
            std::shared_ptr<ffWLC> wlc = std::make_shared<ffWLC>(sys,pd,pg,in);
            sim->addInteractor(wlc);
        } else {
            sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                              "Selected model: %s, is not implemented",model.c_str());
        }
    }
        
    std::map<int,std::shared_ptr<uammd::ParticleGroup>> simGroups;
    {
        std::set<int> simList;
        {
            auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);

            fori(0,pd->getNumParticles()){
                simList.emplace(simId[i]);
            }
        }
        
        for(const int& s : simList){
            selectors::simulationId selector(s);

            auto pgs = std::make_shared<uammd::ParticleGroup>(selector,
                                                              pd,
                                                              sys,
                                                              "simId_"+std::to_string(s));
            simGroups[s]=pgs;
            
        }
    }
        
    std::map<int,std::string> simId2folder;
    {
        std::string outPutFolders = in.getOption("outPutFolders",uammd::InputFile::Required).str();
        
        std::ifstream outPutFoldersFile(outPutFolders);

        std::string line;
        while(std::getline(outPutFoldersFile,line)){
            std::stringstream ss(line);

            int simId;
            std::string path;

            ss >> simId >> path;

            simId2folder[simId] = path;
        }
    }

    //Trajectory output
    {
        WriteStep<ff::Units>::Parameters paramBase = WriteStep<ff::Units>::inputFileToParam(in);

        int interval = std::stoi(in.getOption("nStepsWriteInterval",uammd::InputFile::Required).str());

        for(auto pg : simGroups){    
            
            auto groupIndex  = pg.second->getIndexIterator(uammd::access::location::cpu);
            
            auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
            int  s = simId[groupIndex[0]];
            
            WriteStep<ff::Units>::Parameters param = paramBase;
            param.outPutFilePath = simId2folder[s]+"/"+paramBase.outPutFilePath+"_"+std::to_string(s);
            
            if(loadFromBackup){
                if(paramBase.outPutFormat == "lammpstrj"){
                    int initStep;
                    uammd::real dt;
                    
                    in.getOption("initStep",uammd::InputFile::Required)>>initStep;
                    in.getOption("dt",uammd::InputFile::Required)>>dt;
                    
                    uammd::real backupTime = dt*initStep;
                    
                    std::string prevTrajName = param.outPutFilePath + "."+
                                               paramBase.outPutFormat;
                    std::ifstream prevTraj(prevTrajName);
                    
                    std::string tmpTrajName = simId2folder[s]+"/tmpTraj."+paramBase.outPutFormat;
                    std::ofstream tmpTraj(tmpTrajName);

                    std::string line;
                    while (std::getline(prevTraj, line))
                    {
                        std::istringstream ss(line);

                        if(line.find("TIMESTEP") != std::string::npos){
                            std::getline(prevTraj, line);
                            uammd::real fileTime = std::atof(line.c_str());

                            if(fileTime>backupTime){break;}
                            else{
                                tmpTraj << "ITEM: TIMESTEP" << std::endl;
                            }
                        }
                        tmpTraj << line << std::endl;
                    }

                    std::filesystem::copy(tmpTrajName, prevTrajName,
                                          std::filesystem::copy_options::overwrite_existing);
                    std::filesystem::remove(tmpTrajName);

                } else {
                    sys->log<uammd::System::ERROR>("[MADnaLAB] "
                                                   "Current output format:%s, can not be fixed for backup. "
                                                   "Some frames of the previous simulations could be present in the output file !!!",
                                                   paramBase.outPutFormat.c_str());
                }
            }

            param.append = loadFromBackup;
            std::shared_ptr<WriteStep<ff::Units>> wStep = std::make_shared<WriteStep<ff::Units>>(sys,
                                                                                                 pd,
                                                                                                 pg.second,
                                                                                                 interval,
                                                                                                 param);

            wStep->setPBC(false);
            sim->addSimulationStep(wStep);
    
        }
    }

    //Backup
    std::shared_ptr<WriteStep<ff::Units>> backupStep;
    {
        WriteStep<ff::Units>::Parameters param;
        param.outPutFormat = "back";

        int interval = std::stoi(in.getOption("nStepsBackupInterval",uammd::InputFile::Required).str());

        auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);
        
        param.outPutFilePath = simulationSetFolder+"/backup";
        
        backupStep = std::make_shared<WriteStep<ff::Units>>(sys,
                                                            pd,
                                                            pg,
                                                            interval,
                                                            param);

        backupStep->setPBC(false);
        sim->addSimulationStep(backupStep);

        backupStep->tryApplyStep(sim->getStep(),0,true);
    }

    //External force beteew centers of mass
    std::string name = "externalForceBtwCOM";
    bool externalForceBtwCOMActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(externalForceBtwCOMActive){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());

        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalForceBtwCOM{
            
            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real force;
            };
        
            static info getInfo(std::stringstream& ss){

                int n1;
                int n2;

                info infoBuffer;
                infoBuffer.bp.resize(2);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n1;
                infoBuffer.bp[0].resize(n1);
                for(uint i=0;i<n1;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> n2;
                infoBuffer.bp[1].resize(n2);
                for(uint i=0;i<n2;i++){ss >> infoBuffer.bp[1][i];}

                ss >> infoBuffer.force;

                return infoBuffer;
            }
        };
        
        auto infoList = getInfoList<externalForceBtwCOM>(sys,name,filePath);

        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
            
            int nSets = 0;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    nSets++;
                }
            }

            setInfo sI1 = getSet2id<externalForceBtwCOM>(sys,pd,simGroups,top,infoList,i,0,model);
            setInfo sI2 = getSet2id<externalForceBtwCOM>(sys,pd,simGroups,top,infoList,i,1,model);

            thrust::host_vector<uammd::real> F;

            for(auto& info : infoList){
                if(info.second.size() > i){
                    F.push_back(ff::Units::TO_INTERNAL_FORCE*info.second[i].force);
                }
            }
            
            std::shared_ptr<Interactor::Sets::ConstantForceBtwCOM> externalForceBtwCOMInteractor
            = std::make_shared<Interactor::Sets::ConstantForceBtwCOM>(sys,pd,pg,
                                                                      sI1.setSize,sI2.setSize,
                                                                      nSets,
                                                                      sI1.set2id,sI2.set2id,
                                                                      F,
                                                                      Interactor::Sets::ConstantForceBtwCOM::Parameters());
             sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                              "Added \"%s\" with setSize1: %i, setSize2: %i and nSets: %i",
                                               name.c_str(),sI1.setSize,sI2.setSize,nSets);
          
            sim->addInteractor(externalForceBtwCOMInteractor);
        }
    } 
    
    //External force
    name = "externalForce";
    bool externalForceActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(externalForceActive){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalForce{
            
            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real3 force;
            };
        
            static info getInfo(std::stringstream& ss){
                
                int n;

                info infoBuffer;
                infoBuffer.bp.resize(1);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n;
                infoBuffer.bp[0].resize(n);
                for(uint i=0;i<n;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> infoBuffer.force.x >> 
                      infoBuffer.force.y >> 
                      infoBuffer.force.z;

                return infoBuffer;
            }
        };
        
        auto infoList = getInfoList<externalForce>(sys,name,filePath);

        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
            
            int nSets = 0;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    nSets++;
                }
            }

            setInfo sI = getSet2id<externalForce>(sys,pd,simGroups,top,infoList,i,0,model);

            thrust::host_vector<uammd::real3> F;

            for(auto& info : infoList){
                if(info.second.size() > i){
                    F.push_back(ff::Units::TO_INTERNAL_FORCE*info.second[i].force);
                }
            }
            
            std::shared_ptr<Interactor::Sets::ExternalForceOverCOM> externalForceInteractor
            = std::make_shared<Interactor::Sets::ExternalForceOverCOM>(sys,pd,pg,
                                                                       sI.setSize,
                                                                       nSets,
                                                                       sI.set2id,
                                                                       F,
                                                                       Interactor::Sets::ExternalForceOverCOM::Parameters());
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Added \"%s\" with setSize: %i and nSets: %i",
                                              name.c_str(),sI.setSize,nSets);
            
            sim->addInteractor(externalForceInteractor);
        }
    }

    //External torque
    name = "externalTorque";
    bool externalTorqueActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(externalTorqueActive){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalTorque{
            
            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real3 torque;
            };
        
            static info getInfo(std::stringstream& ss){
                
                int n;

                info infoBuffer;
                infoBuffer.bp.resize(1);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n;
                infoBuffer.bp[0].resize(n);
                for(uint i=0;i<n;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> infoBuffer.torque.x >> 
                      infoBuffer.torque.y >> 
                      infoBuffer.torque.z;

                return infoBuffer;
            }
        };
        
        auto infoList = getInfoList<externalTorque>(sys,name,filePath);

        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
            
            int nSets = 0;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    nSets++;
                }
            }

            setInfo sI = getSet2id<externalTorque>(sys,pd,simGroups,top,infoList,i,0,model);

            thrust::host_vector<uammd::real3> T;

            for(auto& info : infoList){
                if(info.second.size() > i){
                    T.push_back(ff::Units::TO_INTERNAL_FORCE*info.second[i].torque);
                }
            }
            
            std::shared_ptr<Interactor::Sets::ExternalTorqueOverCOM> externalTorqueInteractor
            = std::make_shared<Interactor::Sets::ExternalTorqueOverCOM>(sys,pd,pg,
                                                                        sI.setSize,
                                                                        nSets,
                                                                        sI.set2id,
                                                                        T,
                                                                        Interactor::Sets::ExternalTorqueOverCOM::Parameters());
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Added \"%s\" with setSize: %i and nSets: %i",
                                              name.c_str(),sI.setSize,nSets);
            
            sim->addInteractor(externalTorqueInteractor);
        }
    }
    
    //Constraints

    //Constraint distance beteew centers of mass
    name = "constraintsDistanceBtwCOM";
    bool constraintsDistanceBtwCOMActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(constraintsDistanceBtwCOMActive){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsDistanceBtwCOM{
            
            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real r0;
                uammd::real K;
            };
        
            static info getInfo(std::stringstream& ss){

                int n1;
                int n2;

                info infoBuffer;
                infoBuffer.bp.resize(2);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n1;
                infoBuffer.bp[0].resize(n1);
                for(uint i=0;i<n1;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> n2;
                infoBuffer.bp[1].resize(n2);
                for(uint i=0;i<n2;i++){ss >> infoBuffer.bp[1][i];}

                ss >> infoBuffer.r0;
                ss >> infoBuffer.K;

                return infoBuffer;
            }
        };
        
        auto infoList = getInfoList<constraintsDistanceBtwCOM>(sys,name,filePath);

        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
            
            int nSets = 0;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    nSets++;
                }
            }

            setInfo sI1 = getSet2id<constraintsDistanceBtwCOM>(sys,pd,simGroups,top,infoList,i,0,model);
            setInfo sI2 = getSet2id<constraintsDistanceBtwCOM>(sys,pd,simGroups,top,infoList,i,1,model);

            thrust::host_vector<uammd::real> r0;
            thrust::host_vector<uammd::real> K;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    r0.push_back(info.second[i].r0);
                    K.push_back(info.second[i].K);
                }
            }

            std::shared_ptr<Interactor::Sets::HarmonicBondBtwCOM> constraintsDistanceBtwCOMInteractor
            = std::make_shared<Interactor::Sets::HarmonicBondBtwCOM>(sys,pd,pg,
                                                                      sI1.setSize,sI2.setSize,
                                                                      nSets,
                                                                      sI1.set2id,sI2.set2id,
                                                                      r0,K,
                                                                      Interactor::Sets::HarmonicBondBtwCOM::Parameters());
             sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                              "Added \"%s\" with setSize1: %i, setSize2: %i and nSets: %i",
                                               name.c_str(),sI1.setSize,sI2.setSize,nSets);
          
            sim->addInteractor(constraintsDistanceBtwCOMInteractor);
        }
    }

    //Constraint position of center of mass
    name = "constraintsPositionOfCOM";
    bool constraintsPositionOfCOMActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(constraintsPositionOfCOMActive){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsPositionOfCOM{
           

            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real3 fixedPoint;
                uammd::real3 K;
            };
        
            static info getInfo(std::stringstream& ss){
                
                int n;

                info infoBuffer;
                infoBuffer.bp.resize(1);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n;
                infoBuffer.bp[0].resize(n);
                for(uint i=0;i<n;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> infoBuffer.fixedPoint.x >> 
                      infoBuffer.fixedPoint.y >> 
                      infoBuffer.fixedPoint.z;

                ss >> infoBuffer.K.x >> 
                      infoBuffer.K.y >> 
                      infoBuffer.K.z;

                return infoBuffer;
            }
        };
        
        auto infoList = getInfoList<constraintsPositionOfCOM>(sys,name,filePath);

        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
            
            int nSets = 0;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    nSets++;
                }
            }

            setInfo sI = getSet2id<constraintsPositionOfCOM>(sys,pd,simGroups,top,infoList,i,0,model);
            
            thrust::host_vector<uammd::real3> fixedPoint;
            thrust::host_vector<uammd::real3> K;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    fixedPoint.push_back(info.second[i].fixedPoint);
                    K.push_back(info.second[i].K);
                }
            }

            std::shared_ptr<Interactor::Sets::HarmonicFixedCOM> constraintsPositionOfCOMInteractor 
            = std::make_shared<Interactor::Sets::HarmonicFixedCOM>(sys,pd,pg,
                                                                       sI.setSize,
                                                                       nSets,
                                                                       sI.set2id,
                                                                       fixedPoint,
                                                                       K,
                                                                       Interactor::Sets::HarmonicFixedCOM::Parameters());
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Added \"%s\" with setSize: %i and nSets: %i",
                                              name.c_str(),sI.setSize,nSets);
            
            sim->addInteractor(constraintsPositionOfCOMInteractor);
        }
    } 
    
    //Constraint position of beads
    name = "constraintsPositionOfBeads";
    bool constraintsPositionOfBeadsActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    if(constraintsPositionOfBeadsActive){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsPositionOfBeads{
            
            struct info{
                int simId;
                std::string type;
                std::vector<std::vector<int>> bp;
                uammd::real3 K;
            };
        
            static info getInfo(std::stringstream& ss){
                int n;

                info infoBuffer;
                infoBuffer.bp.resize(1);

                ss >> infoBuffer.simId >>  
                      infoBuffer.type;

                ss >> n;
                infoBuffer.bp[0].resize(n);
                for(uint i=0;i<n;i++){ss >> infoBuffer.bp[0][i];}
                
                ss >> infoBuffer.K.x >> 
                      infoBuffer.K.y >> 
                      infoBuffer.K.z;

                return infoBuffer;
            }
        };
    
        using FixedType                 = Potentials::Bond1::HarmonicConst_r0;
        
        using InteractorHarmonicFixed   = Interactor::BondedInteractor<FixedType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<FixedType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromVector<FixedType>>;

        FixedType::Parameters fixedParameters;

        fixedParameters.r0 = {0.0,0.0,0.0};

        std::shared_ptr<FixedType> fb = std::make_shared<FixedType>(pd,
                                                                    fixedParameters);
        
        typename InteractorHarmonicFixed::Parameters interactorFixedParameters;
        
        interactorFixedParameters.bondName = "FIXED";
                
        //Load fixed

        std::shared_ptr<std::vector<FixedType::Bond>> fixedVector = std::make_shared<std::vector<FixedType::Bond>>();
        
        const int* id2index = pd->getIdOrderedIndices(uammd::access::location::cpu);
        auto pos    = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
        auto simId  = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);

        /////////////////////////////////////////////////////////////////
            
        auto infoList = getInfoList<constraintsPositionOfBeads>(sys,name,filePath);
        
        uint maxEntrySize=0;
        for(auto& info : infoList){
            maxEntrySize = std::max(maxEntrySize,uint(info.second.size()));
        }

        for(int i=0;i<maxEntrySize;i++){
        
            setInfo sI = getSet2id<constraintsPositionOfBeads>(sys,pd,simGroups,top,infoList,i,0,model);
            
            std::map<int,uammd::real3> simId2K;
            for(auto& info : infoList){
                if(info.second.size() > i){
                    simId2K[info.second[i].simId] = info.second[i].K;
                }
            }
                
            if(!loadFromBackup){
                for(int id : sI.set2id){
                    int index = id2index[id];

                    FixedType::Bond bi;

                    bi.i = id;
                    bi.bondInfo.K = simId2K[simId[index]];
                    bi.bondInfo.pos = uammd::make_real3(pos[index]);

                    fixedVector->push_back(bi);
                }
            } else {
                for(int id : sI.set2id){
                    int index = id2index[id];

                    FixedType::Bond bi;

                    bi.i = id;
                    bi.bondInfo.K = simId2K[simId[index]];

                    in.getOption(std::to_string(id)+"fixedPos",uammd::InputFile::Required) >> bi.bondInfo.pos.x 
                                                                                           >> bi.bondInfo.pos.y 
                                                                                           >> bi.bondInfo.pos.z;
                    
                    fixedVector->push_back(bi);
                }
            }
            
            for(auto fi : *fixedVector){
                backupInfo[std::to_string(fi.i)+"fixedPos"]= std::to_string(fi.bondInfo.pos.x) + " " +
                                                             std::to_string(fi.bondInfo.pos.y) + " " +
                                                             std::to_string(fi.bondInfo.pos.z);
            }

        }

        std::shared_ptr<InteractorHarmonicFixed> constraintsPositionOfBeadsInteractor = std::make_shared<InteractorHarmonicFixed>(sys, pd, pg,
                                                                                                                                  fixedVector, fb,
                                                                                                                                  interactorFixedParameters);
        
        sim->addInteractor(constraintsPositionOfBeadsInteractor);
    } 
    
    //Boundaries

    //Zplates
    name = "boundaryZPlates";
    bool boundaryZPlatesActive = bool(in.getOption(name+"Active",uammd::InputFile::Optional));
    
    std::shared_ptr<CompressiblePlates> boundaryZPlatesInteractor; //Exposed for backup 
    if(boundaryZPlatesActive){
        
        uammd::real initialPlatesSeparation;
        uammd::real finalPlatesSeparation;
        
        uammd::real compressionVelocity;
        
        in.getOption("initialPlatesSeparation",uammd::InputFile::Required) >> initialPlatesSeparation;
        in.getOption("finalPlatesSeparation",uammd::InputFile::Required) >> finalPlatesSeparation;
        
        in.getOption("compressionVelocity",uammd::InputFile::Required) >> compressionVelocity;
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Initial plates separation: %f",
                                          initialPlatesSeparation);
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Final plates separation: %f",
                                          finalPlatesSeparation);
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Compression velocity: %f",
                                          compressionVelocity);

        CompressiblePlates::Parameters par;

        par.initialPlatesSeparation = initialPlatesSeparation;
        par.endPlatesSeparation = finalPlatesSeparation;
        
        par.compressionVelocity = compressionVelocity/ff::Units::TO_INTERNAL_TIME;
        
        boundaryZPlatesInteractor = std::make_shared<CompressiblePlates>(sys,
                                                                         pd,
                                                                         pg,
                                                                         par);
        
        sim->addInteractor(boundaryZPlatesInteractor);

        if(loadFromBackup){

            uammd::real topPlatePos; 
            uammd::real bottomPlatePos; 

            in.getOption("topPlatePosition",uammd::InputFile::Required) >> topPlatePos;
            in.getOption("bottomPlatePosition",uammd::InputFile::Required) >> bottomPlatePos;
            
            boundaryZPlatesInteractor->setTopPlatePosition(topPlatePos);
            boundaryZPlatesInteractor->setBottomPlatePosition(bottomPlatePos);
        }
    }

    //Measures
    bool measuresActive = bool(in.getOption("measuresActive",uammd::InputFile::Optional));
    if(measuresActive){

        int nStepsMeasure;
        in.getOption("nStepsMeasure",uammd::InputFile::Required) >> nStepsMeasure;
        
        std::stringstream measuresList_stream = std::stringstream(in.getOption("measuresList",uammd::InputFile::Required).str());
        
        std::vector<std::string> measuresList;
        std::string measure;
        while(measuresList_stream >> measure){
            measuresList.push_back(measure);
        }
        
        if(loadFromBackup){
            for(auto pg : simGroups){    
            
                auto groupIndex  = pg.second->getIndexIterator(uammd::access::location::cpu);
                
                auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
                int  s = simId[groupIndex[0]];

                int initStep;
                in.getOption("initStep",uammd::InputFile::Required)>>initStep;
                
                std::string prevMeasureName = simId2folder[s]+"/measures_"+std::to_string(s)+".dat";
                std::ifstream prevMeasure(prevMeasureName);
                
                std::string tmpMeasureName = simId2folder[s]+"/tmpMeasure.dat";
                std::ofstream tmpMeasure(tmpMeasureName);

                std::string line;
                while (std::getline(prevMeasure, line))
                {
                    std::istringstream ss(line);

                    int measureStep;
                    ss >> measureStep;

                    if(measureStep > initStep){break;}

                    tmpMeasure << line << std::endl;
                }

                std::filesystem::copy(tmpMeasureName, prevMeasureName,
                                      std::filesystem::copy_options::overwrite_existing);
                std::filesystem::remove(tmpMeasureName);

            }
        }

        std::shared_ptr<MeasuresList<SIM>> mStep = std::make_shared<MeasuresList<SIM>>(sys,
                                                                                       pd,
                                                                                       pg,
                                                                                       nStepsMeasure,
                                                                                       simGroups,simId2folder,
                                                                                       measuresList,
                                                                                       sim,
                                                                                       loadFromBackup);

        sim->addSimulationStep(mStep);
    }

    //Remove loadFromBackup from input file
    if(loadFromBackup){
        std::filesystem::copy(simulationSetFolder+"/options.back", argv[1],
                              std::filesystem::copy_options::overwrite_existing);
    }
        
    try{
        if(!loadFromBackup){
            sim->run();
        } else {
            sim->getIntegrator()->init();
            sim->run(false);
        }
    } catch (uammd::exception &e) {
        
        sys->log<uammd::System::ERROR>("  [MADnaLAB] "
                                       "Error detected at simulation step: %u",sim->getStep());

        int errorCode = 1;

        if(loadFromBackup){
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Previous backup detected");

            int initStep;
            in.getOption("initStep",uammd::InputFile::Required)>>initStep;

            if(initStep == backupStep->getLastStepApplied()){
                sys->log<uammd::System::WARNING>("[MADnaLAB] "
                                                 "Previus and current backup points match: %i",initStep);
                errorCode = 2;
            } else {
                sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                                 "Previus and current backup points are different. "
                                                 "Number of steps between backup points: %i",backupStep->getLastStepApplied()-initStep
                                                 );
            }
        }

        in.~InputFile();
        
        std::filesystem::copy(argv[1], simulationSetFolder+"/options.back",
                              std::filesystem::copy_options::overwrite_existing);
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Writting backup info. Last backup at step: %u",backupStep->getLastStepApplied());

        std::ofstream inOut;
        inOut.open(argv[1], std::ios_base::app);

        //Write info to input
        inOut << "loadFromBackup" << std::endl;
        inOut << "initStep " << backupStep->getLastStepApplied() << std::endl;

        for(auto bckInf : backupInfo){
            inOut << bckInf.first << " " << bckInf.second << std::endl;
        }
        
        if(boundaryZPlatesActive){
            inOut << "topPlatePosition "    << boundaryZPlatesInteractor->getTopPlatePosition() << std::endl;
            inOut << "bottomPlatePosition " << boundaryZPlatesInteractor->getBottomPlatePosition() << std::endl;
        }
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Exit with error code %i",errorCode);

        std::exit(errorCode);
    }

    return EXIT_SUCCESS;
}


