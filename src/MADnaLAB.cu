#include <UAMMDstructured.cuh>

using namespace uammd::structured;

using ff = forceField::MADna;
using SIM = Simulation<ff,SteepestDescent,LangevinNVT::BBK>;

int2 getPairsIndex(int pairBasePosition,
                   int nBasis, 
                   std::shared_ptr<uammd::System> sys){
    
    if(nBasis%2!=0){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                "The total number of basis should be even, but the value is: %i",
                nBasis);
    }

    int nBasisPairs = nBasis/2;

    if(pairBasePosition == 0){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid pair base position, it must be different than 0");
    }

    if(abs(pairBasePosition) > nBasisPairs ){
        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                          "Invalid pair base position,"
                                          "its absolute value can not be larger "
                                          "than the total number of base pairs (%i)," 
                                          "but the value is: %i",
                                           nBasisPairs,pairBasePosition);
    }

    int2 pairsIndex;
            
    if(pairBasePosition > 0){
        pairsIndex.x = pairBasePosition-1;
        pairsIndex.y = nBasis-1-pairsIndex.x;
    } else {
        pairsIndex.x = nBasisPairs+pairBasePosition;
        pairsIndex.y = nBasis-1-pairsIndex.x;
    }

    return pairsIndex;
}

thrust::host_vector<int> getBasisPairsSet(int pairBasePosition,
                                          std::shared_ptr<uammd::ParticleData> pd,
                                          std::shared_ptr<ff::Topology> top,
                                          std::vector<std::shared_ptr<uammd::ParticleGroup>>& simGroups,
                                          std::shared_ptr<uammd::System> sys){
        
    thrust::host_vector<int> set;
    
    auto pos   = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
        
    auto id    = pd->getId(uammd::access::location::cpu,uammd::access::mode::read);
    auto res   = pd->getResId(uammd::access::location::cpu,uammd::access::mode::read);

    for(auto pg : simGroups){
        
        auto groupIndex = pg->getIndexIterator(uammd::access::location::cpu);
        
        int len=1;
        int basis_prev = res[groupIndex[0]];
        for(int i=0;i<pg->getNumberParticles();i++){
            int index = groupIndex[i];
            
            int basis = res[index];

            if(basis != basis_prev){
                basis_prev=basis;
                len++;
            }
        }

        int2 b1 = getPairsIndex(pairBasePosition,
                                len, 
                                sys);
        int b11 = b1.x;
        int b12 = b1.y;

        for(int i=0;i<pg->getNumberParticles();i++){
            
            int index = groupIndex[i];
        
            int type = pos[index].w;

            int pid   = id[index];
            int basis = res[index];

            if(top->getTypes()->getTypeParameters(type).name == "S"){
                if        (basis == b11){
                    set.push_back(pid);
                } else if (basis == b12){
                    set.push_back(pid);
                }                 
            }
        }
    }

    return set;
}
    
class CompressiblePlates: public Interactor::Plates{

    private:
            
        uammd::real initialPlatesSeparation;
        uammd::real endPlatesSeparation;
            
        uammd::real compressionVelocity;

    public:

        struct Parameters : public Plates::Parameters {

            uammd::real initialPlatesSeparation;
            uammd::real endPlatesSeparation;
            
            uammd::real compressionVelocity;
        
        };

        CompressiblePlates(std::shared_ptr<uammd::System>       sys,
                           std::shared_ptr<uammd::ParticleData>  pd,
                           std::shared_ptr<uammd::ParticleGroup> pg,
                           Parameters par):initialPlatesSeparation(par.initialPlatesSeparation),
                                           endPlatesSeparation(par.endPlatesSeparation),
                                           compressionVelocity(par.compressionVelocity),
                                           Plates(sys,pd,pg,par){
            
            this->setTopPlatePosition(initialPlatesSeparation/2.0);
            this->setBottomPlatePosition(-initialPlatesSeparation/2.0);
        }


        void updateSimulationTime(uammd::real simulationTime) override {
            
            uammd::real platesSep = abs(this->getTopPlatePosition()-this->getBottomPlatePosition());

            if(platesSep > endPlatesSeparation){
                this->setTopPlatePosition(double(initialPlatesSeparation/2.0)-double(simulationTime*compressionVelocity));
                this->setBottomPlatePosition(double(-initialPlatesSeparation/2.0)+double(simulationTime*compressionVelocity));
            }
        }
};

int main(int argc, char** argv){

    auto sys = std::make_shared<uammd::System>();
    
    //ullint seed = 0xf31337Bada55D00dULL^time(NULL);
    ullint seed = 0xf31337Bada55D00dULL;
    sys->rng().setSeed(seed);
    
    uammd::InputFile in("options.dat");

    SIM sim(sys,in);

    auto top = sim.getTopology();
    auto pd  = sim.getParticleData();
    auto pg  = sim.getParticleGroup();
        
    std::vector<std::shared_ptr<uammd::ParticleGroup>> simGroups;
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

            auto pg = std::make_shared<uammd::ParticleGroup>(selector,
                                                             pd,
                                                             sys,
                                                             "simId_"+std::to_string(s));
            simGroups.push_back(pg);
            
        }
    }

    {
        std::ofstream list(in.getOption("outPutFilePath",uammd::InputFile::Required).str()+".list");
        
        WriteStep::Parameters paramBase = WriteStep::inputFileToParam(in);

        int interval = std::stoi(in.getOption("nStepsWriteInterval",uammd::InputFile::Required).str());
        
        for(auto pg : simGroups){    
            
            auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);

            auto pos   = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
            auto chain = pd->getChainId(uammd::access::location::cpu,uammd::access::mode::read);
            
            auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
            int  s = simId[groupIndex[0]];

            for(int i=0;i<pg->getNumberParticles();i++){
            
                int index = groupIndex[i];
            
                int type = pos[index].w;

                auto typeName = top->getTypes()->getTypeParameters(type).name;

                if(typeName != "S" and 
                   typeName != "P" and
                   chain[index] == 0){
                    list << typeName;
                }
            }

            WriteStep::Parameters param = paramBase;
            param.outPutFilePath = paramBase.outPutFilePath+"_"+std::to_string(s);
            
            list << " " << s << " " << param.outPutFilePath << std::endl;
            
            std::shared_ptr<WriteStep> wStep = std::make_shared<WriteStep>(sys,
                                                                           pd,
                                                                           pg,
                                                                           interval,
                                                                           param);

            wStep->setPBC(false);
            sim.addSimulationStep(wStep);
        }
    }

    //Pulling
    if(in.getOption("pullingActive",uammd::InputFile::Optional)){
        
        uammd::real pullingForce;
        in.getOption("pullingForce",uammd::InputFile::Required) >> pullingForce;

        int pullingBasePair1;
        in.getOption("pullingBasePair1",uammd::InputFile::Required) >> pullingBasePair1;

        int pullingBasePair2;
        in.getOption("pullingBasePair2",uammd::InputFile::Required) >> pullingBasePair2;
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Pulling force: %f",
                                          pullingForce);
        
        thrust::host_vector<int> set1 = getBasisPairsSet(pullingBasePair1,
                                                         pd,top,
                                                         simGroups,
                                                         sys);
        
        thrust::host_vector<int> set2 = getBasisPairsSet(pullingBasePair2,
                                                         pd,top,
                                                         simGroups,
                                                         sys);
        
        std::shared_ptr<Interactor::ConstantForceCOMCopies> pullingInteractor = std::make_shared<Interactor::ConstantForceCOMCopies>(sys,
                                                                                                                                     pd,pg,
                                                                                                                                     2,2,
                                                                                                                                     set1,set2,
                                                                                                                                     simGroups.size(),
                                                                                                                                     Interactor::ConstantForceCOMCopies::Parameters());
            
        pullingInteractor->setState(ff::Units::TO_INTERNAL_FORCE*pullingForce);

        sim.addInteractor(pullingInteractor);

    } 
    
    //Pulling external force
    if(in.getOption("pullingExternalForceActive",uammd::InputFile::Optional)){
        
        uammd::real3 pullingExternalForce;
        in.getOption("pullingExternalForce",uammd::InputFile::Required) >> pullingExternalForce.x 
                                                                        >> pullingExternalForce.y
                                                                        >> pullingExternalForce.z;

        int pullingExternalBasePair1;
        in.getOption("pullingExternalBasePair1",uammd::InputFile::Required) >> pullingExternalBasePair1;

        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "Pulling external force: %f,%f,%f",
                                          pullingExternalForce.x,
                                          pullingExternalForce.y,
                                          pullingExternalForce.z);
            
        thrust::host_vector<int> set1 = getBasisPairsSet(pullingExternalBasePair1,
                                                         pd,top,
                                                         simGroups,
                                                         sys);
        
        std::shared_ptr<Interactor::ExternalCOMCopies> pullingExternalInteractor = std::make_shared<Interactor::ExternalCOMCopies>(sys,
                                                                                                                                   pd,pg,
                                                                                                                                   2,
                                                                                                                                   set1,
                                                                                                                                   simGroups.size(),
                                                                                                                                   Interactor::ExternalCOMCopies::Parameters());
            
        pullingExternalInteractor->setState(ff::Units::TO_INTERNAL_FORCE*pullingExternalForce);

        sim.addInteractor(pullingExternalInteractor);
    } 
    
    //Constraints 

    //fixed distance between COM of two particles (Harmonic)
    if(in.getOption("harmonicCOMdistanceActive",uammd::InputFile::Optional)){
        
        std::stringstream constrainedBasePair1_stream = std::stringstream(in.getOption("constrainedBasePair1",uammd::InputFile::Required).str());
        std::stringstream constrainedBasePair2_stream = std::stringstream(in.getOption("constrainedBasePair2",uammd::InputFile::Required).str());
        std::stringstream COMdistance_stream          = std::stringstream(in.getOption("COMdistance",uammd::InputFile::Required).str());
        std::stringstream COMk_stream                 = std::stringstream(in.getOption("COMk",uammd::InputFile::Required).str());
        
        int constrainedBasePair1;
        int constrainedBasePair2;
        
        while(constrainedBasePair1_stream >> constrainedBasePair1){
              constrainedBasePair2_stream >> constrainedBasePair2;
            
            uammd::real COMdistance;
            uammd::real COMk;

            COMdistance_stream >> COMdistance;
            COMk_stream >> COMk;
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Base pairs constrained: %i, %i",
                                              constrainedBasePair1,
                                              constrainedBasePair2);

            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Constrained distance: %f",
                                              COMdistance);
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Constrained K: %f",
                                              COMk);
                
            thrust::host_vector<int> set1 = getBasisPairsSet(constrainedBasePair1,
                                                             pd,top,
                                                             simGroups,
                                                             sys);
            
            thrust::host_vector<int> set2 = getBasisPairsSet(constrainedBasePair2,
                                                             pd,top,
                                                             simGroups,
                                                             sys);

            Interactor::HarmonicCOMCopies::Parameters param;
            param.K = COMk;
            
            std::shared_ptr<Interactor::HarmonicCOMCopies> harmonicCOMInteractor = std::make_shared<Interactor::HarmonicCOMCopies>(sys,
                                                                                                                                   pd,pg,
                                                                                                                                   2,2,
                                                                                                                                   set1,set2,
                                                                                                                                   simGroups.size(),
                                                                                                                                   param);
                
            harmonicCOMInteractor->setState(COMdistance);

            sim.addInteractor(harmonicCOMInteractor);

        }
    } 
    
    //COM fixed (Harmonic)
    if(in.getOption("constraintHarmonicCOMfixedActive",uammd::InputFile::Optional)){
        
        std::stringstream fixedBasePair1_stream = std::stringstream(in.getOption("fixedBasePair1",uammd::InputFile::Required).str());
        std::stringstream fixedPoint_stream     = std::stringstream(in.getOption("fixedPoint",uammd::InputFile::Required).str());
        std::stringstream Kfixed_stream         = std::stringstream(in.getOption("Kfixed",uammd::InputFile::Required).str());
        
        int fixedBasePair1;
        while(fixedBasePair1_stream >> fixedBasePair1){
            
            uammd::real3 fixedPoint;
            uammd::real3 Kfixed;

            fixedPoint_stream >> fixedPoint.x >> fixedPoint.y >> fixedPoint.z;
            Kfixed_stream >> Kfixed.x >> Kfixed.y >> Kfixed.z;
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Fixed point on base pair: %i",
                                              fixedBasePair1);

            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Fixed point: %f,%f,%f",
                                              fixedPoint.x,
                                              fixedPoint.y,
                                              fixedPoint.z);
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Kfixed: %f,%f,%f",
                                              Kfixed.x,
                                              Kfixed.y,
                                              Kfixed.z);
                
            thrust::host_vector<int> set1 = getBasisPairsSet(fixedBasePair1,
                                                             pd,top,
                                                             simGroups,
                                                             sys);

            Interactor::HarmonicFixedCOMCopies::Parameters param;
            param.K = Kfixed;
            
            std::shared_ptr<Interactor::HarmonicFixedCOMCopies> harmonicFixedInteractor = std::make_shared<Interactor::HarmonicFixedCOMCopies>(sys,
                                                                                                                                               pd,pg,
                                                                                                                                               2,
                                                                                                                                               set1,
                                                                                                                                               simGroups.size(),
                                                                                                                                               param);
                
            harmonicFixedInteractor->setState(fixedPoint);

            sim.addInteractor(harmonicFixedInteractor);

        }

    } 
    
    //Boundaries

    //Zplates
    if(in.getOption("boundaryZPlatesActive",uammd::InputFile::Optional)){
        
        uammd::real initialPlatesSeparation;
        uammd::real endPlatesSeparation;
        
        uammd::real compressionVelocity;
        
        in.getOption("initialPlatesSeparation",uammd::InputFile::Required) >> initialPlatesSeparation;
        in.getOption("endPlatesSeparation",uammd::InputFile::Required) >> endPlatesSeparation;
        
        in.getOption("compressionVelocity",uammd::InputFile::Required) >> compressionVelocity;

        CompressiblePlates::Parameters par;

        par.initialPlatesSeparation = initialPlatesSeparation;
        par.endPlatesSeparation = endPlatesSeparation;
        
        par.compressionVelocity = compressionVelocity/ff::Units::TO_INTERNAL_TIME;
        
        std::shared_ptr<CompressiblePlates> compressiblePlates = std::make_shared<CompressiblePlates>(sys,
                                                                                                      pd,
                                                                                                      pg,
                                                                                                      par);
        
        sim.addInteractor(compressiblePlates);

    }

    sim.run();

    return EXIT_SUCCESS;
}


