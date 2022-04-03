#ifndef __MADNALAB_MEASURES__
#define __MADNALAB_MEASURES__

//Measures
template<class SimulationType>
class MeasuresList: public SimulationStep{

        uammd::real simulationTime;

        std::map<int,std::ofstream> measuresFiles;

        std::map<int,std::shared_ptr<uammd::ParticleGroup>> simGroups;
        std::map<int,std::string> simId2folder;
        std::vector<std::string> measuresList;
        
        std::shared_ptr<SimulationType> sim;

    public:
        
        MeasuresList(std::shared_ptr<uammd::System>       sys,
                     std::shared_ptr<uammd::ParticleData>  pd,
                     std::shared_ptr<uammd::ParticleGroup> pg,
                     int interval,
                     std::map<int,std::shared_ptr<uammd::ParticleGroup>>& simGroups,
                     std::map<int,std::string> simId2folder,
                     std::vector<std::string> measuresList,
                     std::shared_ptr<SimulationType> sim):SimulationStep(sys,pd,pg,"Measures",interval),
                                                          simGroups(simGroups),
                                                          simId2folder(simId2folder),
                                                          measuresList(measuresList),
                                                          sim(sim),
                                                          simulationTime(0.0){
            for(auto pgs : simGroups){    
                
                auto groupIndex  = pgs.second->getIndexIterator(uammd::access::location::cpu);

                auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);
                int  s = simId[groupIndex[0]];

                measuresFiles.emplace(s,std::ofstream(simId2folder[s]+"/measures_"+std::to_string(s)+".dat"));
            
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
            
            for(auto pgs : simGroups){    
                
                auto groupIndex  = pgs.second->getIndexIterator(uammd::access::location::cpu);
                
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
                                                                      pgs.second,st);
                        
                        int N = pgs.second->getNumberParticles();
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
                                                                           pgs.second,st);

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
                                                                           pgs.second,st);

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

#endif
