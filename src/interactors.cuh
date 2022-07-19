#ifndef __MADNALAB_INTERACTORS__
#define __MADNALAB_INTERACTORS__

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

        CompressiblePlates(std::shared_ptr<uammd::ParticleGroup> pg,
                           Parameters par):initialPlatesSeparation(par.initialPlatesSeparation),
                                           endPlatesSeparation(par.endPlatesSeparation),
                                           compressionVelocity(par.compressionVelocity),
                                           Plates(pg,par){
            
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


#endif
