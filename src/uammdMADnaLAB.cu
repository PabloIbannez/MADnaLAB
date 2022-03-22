#include <UAMMDstructured.cuh>

using namespace uammd::structured;

using ff = forceField::MADna;
using SIM = Simulation<ff,SteepestDescent,LangevinNVT::BBK>;

//Aux functions

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
                                          std::vector<std::string>& types,
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

            std::string typeName = top->getTypes()->getTypeParameters(type).name;

            if(std::find(types.begin(),types.end(),typeName) != types.end()){
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

//Additional interactos
    
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

//Measures

template<class SimulationType>
class MeasuresList: public SimulationStep{

        uammd::real simulationTime;

        std::map<int,std::ofstream> measuresFiles;

        std::vector<std::shared_ptr<uammd::ParticleGroup>> simGroups;
        std::vector<std::string> measuresList;
        
        std::shared_ptr<SimulationType> sim;

    public:
        
        MeasuresList(std::shared_ptr<uammd::System>       sys,
                     std::shared_ptr<uammd::ParticleData>  pd,
                     std::shared_ptr<uammd::ParticleGroup> pg,
                     int interval,
                     std::vector<std::shared_ptr<uammd::ParticleGroup>> simGroups,
                     std::vector<std::string> measuresList,
                     std::shared_ptr<SimulationType> sim):SimulationStep(sys,pd,pg,"Measures",interval),
                                                          simGroups(simGroups),
                                                          measuresList(measuresList),
                                                          sim(sim),
                                                          simulationTime(0.0){
            for(auto pg : simGroups){    
                
                auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);

                auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
                int  s = simId[groupIndex[0]];

                measuresFiles.emplace(s,std::ofstream("measures_"+std::to_string(s)+".dat"));
            
                measuresFiles[s] << "# step time ";

                for(std::string m : measuresList){
                    if        (m == "temperature"){
                        measuresFiles[s] << "temperature ";
                    } else if (m == "energy"){
                        for(std::string pot : sim->getForceField()->getComponentsList()){
                            measuresFiles[s] << pot << " ";
                        }
                        for(auto& inter : sim->getIntegrator()->getInteractors()){
                            std::string interactorName = inter->getName();

                            interactorName.erase(std::remove_if(std::begin(interactorName),std::end(interactorName),
                                                                [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), std::end(interactorName));

                            measuresFiles[s] << interactorName << " ";
                        }
                    } else {
                        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                                          "Measure, %s, is not implemented",
                                                          m.c_str());
                    }
                }
                    
                measuresFiles[s] << std::endl;
            }
        }

        void init(cudaStream_t st) override{}

        void applyStep(int step, cudaStream_t st) override{
            
            for(auto pg : simGroups){    
                
                auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);
                
                int  s;
                {
                    auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
                    s = simId[groupIndex[0]];
                }

                measuresFiles[s] << step << " " << simulationTime*SimulationType::ForceField::Units::FROM_INTERNAL_TIME << " ";

                for(std::string m : measuresList){
                    if        (m == "temperature"){
                        uammd::real kE = Measures::totalKineticEnergy(this->sys,
                                                                      this->pd,
                                                                      this->pg,st);
                        
                        int N = this->pg->getNumberParticles();
                        uammd::real T = uammd::real(2.0/(3.0*N*SimulationType::ForceField::Units::KBOLTZ))*kE;
                        measuresFiles[s] << T << " ";

                    } else if (m == "energy"){
                        for(std::string pot : sim->getForceField()->getComponentsList()){

                            {
                                auto energy = pd->getEnergy(uammd::access::location::gpu, uammd::access::mode::write);     
                                thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), uammd::real(0));
                            }
            
                            uammd::Interactor::Computables comp;
                            comp.energy = true;

                            sim->getForceField()->sum(pot,comp,st);

                            uammd::real E = Measures::totalPotentialEnergy(this->sys,
                                                                           this->pd,
                                                                           this->pg,st);

                            measuresFiles[s] << E << " ";
                        }
                        for(auto& inter : sim->getIntegrator()->getInteractors()){
                            
                            {
                                auto energy = pd->getEnergy(uammd::access::location::gpu, uammd::access::mode::write);     
                                thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), uammd::real(0));
                            }
            
                            uammd::Interactor::Computables comp;
                            comp.energy = true;

                            inter->sum(comp,st);

                            uammd::real E = Measures::totalPotentialEnergy(this->sys,
                                                                           this->pd,
                                                                           this->pg,st);

                            measuresFiles[s] << E << " ";
                        }
                    } else {
                        sys->log<uammd::System::CRITICAL>("[MADnaLAB] "
                                                          "Measure, %s, is not implemented",
                                                          m.c_str());
                    }
                }
                    
                measuresFiles[s] << std::endl;
                
            }
        }
        
        void updateSimulationTime(uammd::real newSimulationTime) override {
            simulationTime = newSimulationTime;
        }


};

int main(int argc, char** argv){

    std::vector<std::string> allTypes  = {"S","P","A","C","G","T"};
    std::vector<std::string> sugarType = {"S"};

    auto sys = std::make_shared<uammd::System>();
    
    //ullint seed = 0xf31337Bada55D00dULL^time(NULL);
    ullint seed = 0xf31337Bada55D00dULL;
    sys->rng().setSeed(seed);
    
    uammd::InputFile in("options.dat");

    std::shared_ptr<SIM> sim = std::make_shared<SIM>(sys,in);

    auto top = sim->getTopology();
    auto pd  = sim->getParticleData();
    auto pg  = sim->getParticleGroup();
        
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

        std::ofstream list("simulations.list");
        
        for(auto pg : simGroups){    
            
            auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);

            auto pos   = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);
            auto chain = pd->getChainId(uammd::access::location::cpu,uammd::access::mode::read);
            
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
            
            auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
            int  s = simId[groupIndex[0]];
            
            list << " " << s << std::endl;
        }

    }

    {
        WriteStep<ff::Units>::Parameters paramBase = WriteStep<ff::Units>::inputFileToParam(in);

        int interval = std::stoi(in.getOption("nStepsWriteInterval",uammd::InputFile::Required).str());
        
        for(auto pg : simGroups){    
            
            auto groupIndex  = pg->getIndexIterator(uammd::access::location::cpu);
            
            auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
            int  s = simId[groupIndex[0]];
            
            WriteStep<ff::Units>::Parameters param = paramBase;
            param.outPutFilePath = paramBase.outPutFilePath+"_"+std::to_string(s);
            
            std::shared_ptr<WriteStep<ff::Units>> wStep = std::make_shared<WriteStep<ff::Units>>(sys,
                                                                                                 pd,
                                                                                                 pg,
                                                                                                 interval,
                                                                                                 param);

            wStep->setPBC(false);
            sim->addSimulationStep(wStep);
        }
    }

    //Pulling
    if(in.getOption("externalForceBtwCOMActive",uammd::InputFile::Optional)){
        
        uammd::real forceBtwCOMforce;
        in.getOption("forceBtwCOMforce",uammd::InputFile::Required) >> forceBtwCOMforce;

        int forceBtwCOMbasePair1;
        int forceBtwCOMbasePair2;
        
        in.getOption("forceBtwCOMbasePair1",uammd::InputFile::Required) >> forceBtwCOMbasePair1;
        in.getOption("forceBtwCOMbasePair2",uammd::InputFile::Required) >> forceBtwCOMbasePair2;
        
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "External force between center of mass: %f",
                                          forceBtwCOMforce);
        
        thrust::host_vector<int> set1 = getBasisPairsSet(forceBtwCOMbasePair1,
                                                         pd,top,
                                                         simGroups,
                                                         sugarType,
                                                         sys);
        
        thrust::host_vector<int> set2 = getBasisPairsSet(forceBtwCOMbasePair2,
                                                         pd,top,
                                                         simGroups,
                                                         sugarType,
                                                         sys);
        
        std::shared_ptr<Interactor::ConstantForceCOMCopies> externalForceBtwCOM 
        = std::make_shared<Interactor::ConstantForceCOMCopies>(sys,
                                                               pd,pg,
                                                               2,2,
                                                               set1,set2,
                                                               simGroups.size(),
                                                               Interactor::ConstantForceCOMCopies::Parameters());
            
        externalForceBtwCOM->setState(ff::Units::TO_INTERNAL_FORCE*forceBtwCOMforce);

        sim->addInteractor(externalForceBtwCOM);

    } 
    
    //Pulling external force
    if(in.getOption("externalForceActive",uammd::InputFile::Optional)){
        
        uammd::real3 forceforce;
        in.getOption("forceforce",uammd::InputFile::Required) >> forceforce.x 
                                                              >> forceforce.y
                                                              >> forceforce.z;

        int forcebasePair;
        in.getOption("forcebasePair",uammd::InputFile::Required) >> forcebasePair;

        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "External force on pair: %i",
                                          forcebasePair);
            
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "External force: %f,%f,%f",
                                          forceforce.x,
                                          forceforce.y,
                                          forceforce.z);
        
        thrust::host_vector<int> set1 = getBasisPairsSet(forcebasePair,
                                                         pd,top,
                                                         simGroups,
                                                         sugarType,
                                                         sys);
        
        std::shared_ptr<Interactor::ExternalCOMCopies> externalForceInteractor = std::make_shared<Interactor::ExternalCOMCopies>(sys,
                                                                                                                                 pd,pg,
                                                                                                                                 2,
                                                                                                                                 set1,
                                                                                                                                 simGroups.size(),
                                                                                                                                 Interactor::ExternalCOMCopies::Parameters());
            
        externalForceInteractor->setState(ff::Units::TO_INTERNAL_FORCE*forceforce);

        sim->addInteractor(externalForceInteractor);
    } 
    
    //external torque
    if(in.getOption("externalTorqueActive",uammd::InputFile::Optional)){
        
        uammd::real3 torquetorque;
        in.getOption("torquetorque",uammd::InputFile::Required) >> torquetorque.x 
                                                                >> torquetorque.y
                                                                >> torquetorque.z;

        int torquebasePair;
        in.getOption("torquebasePair",uammd::InputFile::Required) >> torquebasePair;

        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "External torque on pair: %i",
                                          torquebasePair);
            
        sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                         "External torque: %f,%f,%f",
                                          torquetorque.x,
                                          torquetorque.y,
                                          torquetorque.z);
        
        thrust::host_vector<int> set1 = getBasisPairsSet(torquebasePair,
                                                         pd,top,
                                                         simGroups,
                                                         allTypes,
                                                         sys);
        
        std::shared_ptr<Interactor::ExternalTorqueCOMCopies> externalTorqueInteractor = std::make_shared<Interactor::ExternalTorqueCOMCopies>(sys,
                                                                                                                                              pd,pg,
                                                                                                                                              2,
                                                                                                                                              set1,
                                                                                                                                              simGroups.size(),
                                                                                                                                              Interactor::ExternalTorqueCOMCopies::Parameters());
            
        externalTorqueInteractor->setState(ff::Units::TO_INTERNAL_FORCE*torquetorque);

        sim->addInteractor(externalTorqueInteractor);
    } 
    
    //Constraints 

    //fixed distance between COM of two particles (Harmonic)
    if(in.getOption("constraintDistanceBtwCOMActive",uammd::InputFile::Optional)){
        
        std::stringstream distanceBtwCOMbasePair1_stream = std::stringstream(in.getOption("distanceBtwCOMbasePair1",uammd::InputFile::Required).str());
        std::stringstream distanceBtwCOMbasePair2_stream = std::stringstream(in.getOption("distanceBtwCOMbasePair2",uammd::InputFile::Required).str());
        std::stringstream distanceBtwCOMdistance_stream  = std::stringstream(in.getOption("distanceBtwCOMdistance",uammd::InputFile::Required).str());
        std::stringstream distanceBtwCOMK_stream         = std::stringstream(in.getOption("distanceBtwCOMK",uammd::InputFile::Required).str());
        
        int distanceBtwCOMbasePair1;
        int distanceBtwCOMbasePair2;
        
        while(distanceBtwCOMbasePair1_stream >> distanceBtwCOMbasePair1){
            
            distanceBtwCOMbasePair2_stream >> distanceBtwCOMbasePair2;
            
            uammd::real distanceBtwCOMdistance;
            uammd::real distanceBtwCOMK;

            distanceBtwCOMdistance_stream >> distanceBtwCOMdistance;
            distanceBtwCOMK_stream >> distanceBtwCOMK;
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Base pairs constrained: %i, %i",
                                              distanceBtwCOMbasePair1,
                                              distanceBtwCOMbasePair2);

            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Constrained distance: %f",
                                              distanceBtwCOMdistance);
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Constrained K: %f",
                                              distanceBtwCOMK);
                
            thrust::host_vector<int> set1 = getBasisPairsSet(distanceBtwCOMbasePair1,
                                                             pd,top,
                                                             simGroups,
                                                             sugarType,
                                                             sys);
            
            thrust::host_vector<int> set2 = getBasisPairsSet(distanceBtwCOMbasePair2,
                                                             pd,top,
                                                             simGroups,
                                                             sugarType,
                                                             sys);

            Interactor::HarmonicCOMCopies::Parameters param;
            param.K = distanceBtwCOMK;
            
            std::shared_ptr<Interactor::HarmonicCOMCopies> constraintDistanceBtwCOMInteractor 
            = std::make_shared<Interactor::HarmonicCOMCopies>(sys,
                                                              pd,pg,
                                                              2,2,
                                                              set1,set2,
                                                              simGroups.size(),
                                                              param);
                
            constraintDistanceBtwCOMInteractor->setState(distanceBtwCOMdistance);

            sim->addInteractor(constraintDistanceBtwCOMInteractor);

        }
    } 
    
    //COM fixed (Harmonic)
    if(in.getOption("constraintPositionOfCOMActive",uammd::InputFile::Optional)){
        
        std::stringstream positionOfCOMbasePair_stream = std::stringstream(in.getOption("positionOfCOMbasePair",uammd::InputFile::Required).str());
        std::stringstream positionOfCOMposition_stream     = std::stringstream(in.getOption("positionOfCOMposition",uammd::InputFile::Required).str());
        std::stringstream positionOfCOMK_stream         = std::stringstream(in.getOption("positionOfCOMK",uammd::InputFile::Required).str());
        
        int positionOfCOMbasePair;
        while(positionOfCOMbasePair_stream >> positionOfCOMbasePair){
            
            uammd::real3 positionOfCOMposition;
            uammd::real3 positionOfCOMK;

            positionOfCOMposition_stream >> positionOfCOMposition.x >> positionOfCOMposition.y >> positionOfCOMposition.z;
            positionOfCOMK_stream >> positionOfCOMK.x >> positionOfCOMK.y >> positionOfCOMK.z;
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Base pair COM constrained: %i",
                                              positionOfCOMbasePair);

            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "COM position: %f,%f,%f",
                                              positionOfCOMposition.x,
                                              positionOfCOMposition.y,
                                              positionOfCOMposition.z);
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "K: %f,%f,%f",
                                              positionOfCOMK.x,
                                              positionOfCOMK.y,
                                              positionOfCOMK.z);
                
            thrust::host_vector<int> set1 = getBasisPairsSet(positionOfCOMbasePair,
                                                             pd,top,
                                                             simGroups,
                                                             sugarType,
                                                             sys);

            Interactor::HarmonicFixedCOMCopies::Parameters param;
            param.K = positionOfCOMK;
            
            std::shared_ptr<Interactor::HarmonicFixedCOMCopies> constraintPositionOfCOMInteractor 
            = std::make_shared<Interactor::HarmonicFixedCOMCopies>(sys,
                                                                   pd,pg,
                                                                   2,
                                                                   set1,
                                                                   simGroups.size(),
                                                                   param);
                
            constraintPositionOfCOMInteractor->setState(positionOfCOMposition);

            sim->addInteractor(constraintPositionOfCOMInteractor);

        }
    } 
    
    //fixed beads pos (Harmonic)
    if(in.getOption("constraintPositionOfBeadsActive",uammd::InputFile::Optional)){
        
        std::stringstream positionOfBeadsbasePair_stream = std::stringstream(in.getOption("positionOfBeadsbasePair",uammd::InputFile::Required).str());
        std::stringstream positionOfBeadsK_stream         = std::stringstream(in.getOption("positionOfBeadsK",uammd::InputFile::Required).str());
        
        int positionOfBeadsbasePair;
        while(positionOfBeadsbasePair_stream >> positionOfBeadsbasePair){
            
            uammd::real3 positionOfBeadsK;

            positionOfBeadsK_stream >> positionOfBeadsK.x >> positionOfBeadsK.y >> positionOfBeadsK.z;
            
            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "Fixed point on base pair: %i",
                                              positionOfBeadsbasePair);

            sys->log<uammd::System::MESSAGE>("[MADnaLAB] "
                                             "positionOfBeadsK: %f,%f,%f",
                                              positionOfBeadsK.x,
                                              positionOfBeadsK.y,
                                              positionOfBeadsK.z);
                
            thrust::host_vector<int> set1 = getBasisPairsSet(positionOfBeadsbasePair,
                                                             pd,top,
                                                             simGroups,
                                                             allTypes,
                                                             sys);

            using FixedType                 = Potentials::Bond1::HarmonicConst_K_r0;
            using InteractorHarmonicFixed   = Interactor::BondedInteractor<FixedType,
                                                                           Interactor::BondedInteractor_ns::BondProcessor<FixedType>,
                                                                           Interactor::BondedInteractor_ns::BondReaderFromVector<FixedType>>;

            FixedType::Parameters fixedParameters;

            fixedParameters.K  = positionOfBeadsK;
            fixedParameters.r0 = {0.0,0.0,0.0};

            std::shared_ptr<FixedType> fb = std::make_shared<FixedType>(pd,
                                                                        fixedParameters);
            
            typename InteractorHarmonicFixed::Parameters interactorFixedParameters;
            
            interactorFixedParameters.bondName = "FIXED";
                    
            //Load fixed
            std::shared_ptr<std::vector<FixedType::Bond>> fixedVector = std::make_shared<std::vector<FixedType::Bond>>();
            
            const int* id2index = pd->getIdOrderedIndices(uammd::access::location::cpu);
            auto pos   = pd->getPos(uammd::access::location::cpu,uammd::access::mode::read);

            for(int id : set1){
                int index = id2index[id];

                FixedType::Bond bi;

                bi.i = id;
                bi.bondInfo.pos = uammd::make_real3(pos[index]);

                fixedVector->push_back(bi);
            }

            std::shared_ptr<InteractorHarmonicFixed> fixed = std::make_shared<InteractorHarmonicFixed>(sys, pd, pg,
                                                                                                       fixedVector, fb,
                                                                                                       interactorFixedParameters);
            
            sim->addInteractor(fixed);

        }
    } 
    
    //Boundaries

    //Zplates
    if(in.getOption("boundaryZPlatesActive",uammd::InputFile::Optional)){
        
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
        
        std::shared_ptr<CompressiblePlates> compressiblePlates = std::make_shared<CompressiblePlates>(sys,
                                                                                                      pd,
                                                                                                      pg,
                                                                                                      par);
        
        sim->addInteractor(compressiblePlates);

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
                                                                                       simGroups,measuresList,
                                                                                       sim);

        sim->addSimulationStep(mStep);


    }

    sim->run();

    return EXIT_SUCCESS;
}


