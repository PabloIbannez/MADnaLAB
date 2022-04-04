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
    
    //ullint seed = 0xf31337Bada55D00dULL^time(NULL);
    ullint seed = 0xf31337Bada55D00dULL;
    sys->rng().setSeed(seed);
    
    uammd::InputFile in(argv[1]);

    std::shared_ptr<SIM> sim = std::make_shared<SIM>(sys,in);

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
            
            std::shared_ptr<WriteStep<ff::Units>> wStep = std::make_shared<WriteStep<ff::Units>>(sys,
                                                                                                 pd,
                                                                                                 pg.second,
                                                                                                 interval,
                                                                                                 param);

            wStep->setPBC(false);
            sim->addSimulationStep(wStep);
        }
    }

    //External force beteew centers of mass
    std::string name = "externalForceBtwCOM";
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());

        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalForceBtwCOM{
            
            struct info{
                int simId;
                std::string type;
                int bp[2];
                uammd::real force;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.bp[1] >> 
                      infoBuffer.force;

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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalForce{
            
            struct info{
                int simId;
                std::string type;
                int bp[1];
                uammd::real3 force;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.force.x >> 
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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Detected %s ...",name.c_str());
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct externalTorque{
            
            struct info{
                int simId;
                std::string type;
                int bp[1];
                uammd::real3 torque;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.torque.x >> 
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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsDistanceBtwCOM{
            
            struct info{
                int simId;
                std::string type;
                int bp[2];
                uammd::real r0;
                uammd::real K;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.bp[1] >> 
                      infoBuffer.r0 >> 
                      infoBuffer.K;

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

            setInfo sI1 = getSet2id<constraintsDistanceBtwCOM>(sys,pd,simGroups,top,infoList,i,0,model);;
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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsPositionOfCOM{
            
            struct info{
                int simId;
                std::string type;
                int bp[1];
                uammd::real3 fixedPoint;
                uammd::real3 K;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.fixedPoint.x >> infoBuffer.fixedPoint.y >> infoBuffer.fixedPoint.z >>
                      infoBuffer.K.x >> infoBuffer.K.y >> infoBuffer.K.z ;

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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
        std::string filePath;
        in.getOption(name+"FilePath",uammd::InputFile::Required) >> filePath;
        
        struct constraintsPositionOfBeads{
            
            struct info{
                int simId;
                std::string type;
                int bp[1];
                uammd::real3 K;
            };
        
            static info getInfo(std::stringstream& ss){
                info infoBuffer;

                ss >> infoBuffer.simId >>  
                      infoBuffer.type  >>
                      infoBuffer.bp[0] >> 
                      infoBuffer.K.x >> infoBuffer.K.y >> infoBuffer.K.z ;

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
            
            for(int id : sI.set2id){
                int index = id2index[id];

                FixedType::Bond bi;

                bi.i = id;
                bi.bondInfo.K = simId2K[simId[index]];
                bi.bondInfo.pos = uammd::make_real3(pos[index]);

                fixedVector->push_back(bi);
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
    if(in.getOption(name+"Active",uammd::InputFile::Optional)){
        
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
        
        std::shared_ptr<CompressiblePlates> boundaryZPlatesInteractor = std::make_shared<CompressiblePlates>(sys,
                                                                                                             pd,
                                                                                                             pg,
                                                                                                             par);
        
        sim->addInteractor(boundaryZPlatesInteractor);

    }

    //Measures
    if(in.getOption("measuresActive",uammd::InputFile::Optional)){

        int nStepsMeasure;
        in.getOption("nStepsMeasure",uammd::InputFile::Required) >> nStepsMeasure;
        
        std::stringstream measuresList_stream = std::stringstream(in.getOption("measuresList",uammd::InputFile::Required).str());
        
        std::vector<std::string> measuresList;
        std::string measure;
        while(measuresList_stream >> measure){
            measuresList.push_back(measure);
        }
        
        std::shared_ptr<MeasuresList<SIM>> mStep = std::make_shared<MeasuresList<SIM>>(sys,
                                                                                       pd,
                                                                                       pg,
                                                                                       nStepsMeasure,
                                                                                       simGroups,simId2folder,
                                                                                       measuresList,
                                                                                       sim);

        sim->addSimulationStep(mStep);
    }

    sim->run();

    return EXIT_SUCCESS;
}


