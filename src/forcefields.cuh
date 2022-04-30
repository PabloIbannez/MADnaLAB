#ifndef __MADNA_FF__
#define __MADNA_FF__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace WormLikeChain{

    template<class Base_ >
    class WormLikeChain : public Base_{
        
        protected:

            using Base = Base_;
            
            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::KratkyPorod;
            
            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            
            std::shared_ptr<InteractorBondType>   bonds;
            std::shared_ptr<InteractorAngleType>  angles;
        
        public:

            WormLikeChain(std::shared_ptr<System>        sys,
                          std::shared_ptr<ParticleData>  pd,
                          std::shared_ptr<ParticleGroup> pg,
                          InputFile&                     in):Base(sys,pd,pg,in){
                
                //Add bonds
                BondType::Parameters bondParameters;

                std::shared_ptr<BondType> bH_PBC = std::make_shared<BondType>(this->pd,
                                                                              bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, bH_PBC,
                                                             interactorBondParameters);

                //Add angles
                AngleType::Parameters angleParameters;

                std::shared_ptr<AngleType> aKP_PBC = std::make_shared<AngleType>(this->pd,
                                                                                 angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                               this->top, aKP_PBC,
                                                               interactorAngleParameters);
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
                angles->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
                angles->updateBox(box);
            }
    
    };

}

using WLC = WormLikeChain::WormLikeChain<
            ForceFieldBase<UnitsSystem::KCALMOL_A,
                           Types::BASIC>>;

}}}

namespace uammd{
namespace structured{ 
namespace forceField{
namespace MechanicallyAccurateDNA{
    
    template<class Base_ >
    class MechanicallyAccurateDNABonded : public Base_{

            using Base = Base_;
        
        protected:

            using BondType     = Potentials::Bond2::Harmonic;
            using AngleType    = Potentials::Bond3::HarmonicAngular;
            using DihedralType = Potentials::Bond4::Dihedral;

            using InteractorBondType   = Interactor::BondedInteractor<BondType,
                                                                      Interactor::BondedInteractor_ns::BondProcessor<BondType>,
                                                                      Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondType>>;
            using InteractorAngleType   = Interactor::BondedInteractor<AngleType,
                                                                       Interactor::BondedInteractor_ns::BondProcessor<AngleType>,
                                                                       Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,AngleType>>;
            using InteractorDihedralType   = Interactor::BondedInteractor<DihedralType,
                                                                          Interactor::BondedInteractor_ns::BondProcessor<DihedralType>,
                                                                          Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,DihedralType>>;
            
            std::vector<std::string> componentsList = {"bonds","angles","dihedrals"};

            std::shared_ptr<InteractorBondType>     bonds;
            std::shared_ptr<InteractorAngleType>    angles;
            std::shared_ptr<InteractorDihedralType> dihedrals;
            
        public:
        
            MechanicallyAccurateDNABonded(std::shared_ptr<System>        sys,
                                          std::shared_ptr<ParticleData>  pd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          InputFile&                     in):Base(sys,pd,pg,in){

                Base::name = "MechanicallyAccurateDNABonded";
                
                if(!std::is_same<typename Base::Units,
                                 UnitsSystem::KCALMOL_A>::value){
                    sys->log<System::CRITICAL>("[%s] Mechanically Accurate DNA force field is parametrized in the %s units system,"
                                               "but %s units system is provied",
                                                Base::name.c_str(),
                                                UnitsSystem::KCALMOL_A::NAME.c_str(),
                                                Base::Units::NAME.c_str());
                }

                //Add bonds
                BondType::Parameters bondParameters;

                std::shared_ptr<BondType> b_PBC = std::make_shared<BondType>(this->pd,
                                                                             bondParameters);
                
                typename InteractorBondType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS";

                bonds = std::make_shared<InteractorBondType>(this->sys, this->pd, this->pg,
                                                             this->top, b_PBC,
                                                             interactorBondParameters);
                
                //Add angles
                AngleType::Parameters angleParameters;

                std::shared_ptr<AngleType> a_PBC = std::make_shared<AngleType>(this->pd,
                                                                               angleParameters);
                
                typename InteractorAngleType::Parameters interactorAngleParameters;
                
                interactorAngleParameters.bondName = "ANGLES";

                angles = std::make_shared<InteractorAngleType>(this->sys, this->pd, this->pg,
                                                               this->top, a_PBC,
                                                               interactorAngleParameters);
                
                //Add dihedrals
                DihedralType::Parameters dihedralParameters;

                std::shared_ptr<DihedralType> d_PBC = std::make_shared<DihedralType>(this->pd,
                                                                                     dihedralParameters);
                
                typename InteractorDihedralType::Parameters interactorDihedralParameters;
                
                interactorDihedralParameters.bondName = "DIHEDRALS";

                dihedrals = std::make_shared<InteractorDihedralType>(this->sys, this->pd, this->pg,
                                                                     this->top, d_PBC,
                                                                     interactorDihedralParameters);
            }


            std::vector<std::string> getComponentsList(){return componentsList;}
            
            void sum(std::string component,Computables comp,cudaStream_t st) {
                if(std::find(componentsList.begin(),componentsList.end(),component) != componentsList.end()){
                    if       (component=="bonds"){
                        bonds->sum(comp,st);
                    } else if(component=="angles"){
                        angles->sum(comp,st);
                    } else if(component=="dihedrals"){
                        dihedrals->sum(comp,st);
                    } 
                } else {
                    this->sys->template log<System::CRITICAL>("[%s] Requested potential %s to sum. "
                                                                "But %s is not present in the force field",
                                                                Base::name.c_str(),
                                                                component.c_str(),component.c_str());
                }
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds->sum(comp,st);
                angles->sum(comp,st);
                dihedrals->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds->updateBox(box);
                angles->updateBox(box);
                dihedrals->updateBox(box);
            }
    
    };

    template<class Base_ >
    class MechanicallyAccurateDNA : public MechanicallyAccurateDNABonded<Base_>{

            using Base = MechanicallyAccurateDNABonded<Base_>;
        
        protected:

            const real epsilon = 1.0;
            
            using DHType    = Potentials::UnBound::DebyeHuckel<typename Base::Units>;
            using WCAType  = Potentials::UnBound::WCA;
            
            using InteractorDHType   = Interactor::PairInteractor<DHType,typename Base::NeighbourList>;
            using InteractorWCAType = Interactor::PairInteractor<WCAType,typename Base::NeighbourList>;

            std::shared_ptr<InteractorDHType>  dh;
            std::shared_ptr<InteractorWCAType> wca;
        
        
        protected:
            
            real cutOffDstDH;
            real cutOffDstWCA;

            real dielectricConstant;
            real debyeLength;

        public:
        
            MechanicallyAccurateDNA(std::shared_ptr<System>        sys,
                                    std::shared_ptr<ParticleData>  pd,
                                    std::shared_ptr<ParticleGroup> pg,
                                    InputFile&                     in):Base(sys,pd,pg,in),
                                                                       cutOffDstDH(std::stof(in.getOption("cutOffDstDH",InputFile::Required).str())),
                                                                       cutOffDstWCA(std::stof(in.getOption("cutOffDstWCA",InputFile::Required).str())),
                                                                       dielectricConstant(std::stof(in.getOption("dielectricConstant",InputFile::Required).str())),
                                                                       debyeLength(std::stof(in.getOption("debyeLength",InputFile::Required).str())){
                
                Base::name = "MechanicallyAccurateDNA";

                Base::componentsList.push_back("dh");
                Base::componentsList.push_back("wca");

                if(cutOffDstDH >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[%s] cutOffDstDH (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                 Base::name.c_str(),
                                                 cutOffDstDH,this->nl->getCutOffVerlet());
                }
                
                if(cutOffDstWCA >= this->nl->getCutOffVerlet()){
                    sys->log<System::CRITICAL>("[%s] cutOffDstWCA (%f) "
                                                 "has to be smaller than VerletListDst (%f)",
                                                 Base::name.c_str(),
                                                 cutOffDstWCA,this->nl->getCutOffVerlet());
                }
                
                this->sys->template log<System::MESSAGE>("[%s] "
                                                         "Parameter cutOffDstDH added: %f",
                                                          Base::name.c_str(),cutOffDstDH);
                this->sys->template log<System::MESSAGE>("[%s] "
                                                         "Parameter cutOffDstWCA added: %f",
                                                          Base::name.c_str(),cutOffDstWCA);
                
                this->sys->template log<System::MESSAGE>("[%s] "
                                                         "Parameter dielectricConstant added: %f",
                                                          Base::name.c_str(),dielectricConstant);
                this->sys->template log<System::MESSAGE>("[%s] "
                                                         "Parameter debyeLength added: %f",
                                                          Base::name.c_str(),debyeLength);
                
                this->nl->setCutOff(std::max(this->nl->getCutOff(),
                                             std::max(cutOffDstDH,cutOffDstWCA)));
                
                //Add dh
                typename DHType::Parameters dhPotentialParam;
                
                dhPotentialParam.dielectricConstant = dielectricConstant;
                dhPotentialParam.debyeLength        = debyeLength;
                
                dhPotentialParam.cutOff = cutOffDstDH;
                
                std::shared_ptr<DHType> potDH = std::make_shared<DHType>(this->pd,dhPotentialParam);

                typename InteractorDHType::Parameters interactorDHParameters;

                interactorDHParameters.name = "DH";
                interactorDHParameters.pot  = potDH;
                interactorDHParameters.nl   = this->nl;
                interactorDHParameters.conditionInteractionName = "charged";

                dh = std::make_shared<InteractorDHType>(this->sys,this->pd,this->pg,
                                                        interactorDHParameters);
                
                //Add wca
                typename WCAType::Parameters wcaPotentialParam;
                
                wcaPotentialParam.epsilon = epsilon;
                
                wcaPotentialParam.cutOff = cutOffDstWCA;
                
                std::shared_ptr<WCAType> potWCA = std::make_shared<WCAType>(this->pd,wcaPotentialParam);

                typename InteractorWCAType::Parameters interactorWCAParameters;

                interactorWCAParameters.name = "WCA";
                interactorWCAParameters.pot  = potWCA;
                interactorWCAParameters.nl   = this->nl;
                interactorWCAParameters.conditionInteractionName = "nonExcluded";

                wca = std::make_shared<InteractorWCAType>(this->sys,this->pd,this->pg,
                                                          interactorWCAParameters);
            }


            void sum(std::string component,Computables comp,cudaStream_t st) {
                if(std::find(Base::componentsList.begin(),Base::componentsList.end(),component) != Base::componentsList.end()){
                           if(component=="dh"){
                        dh->sum(comp,st);
                    } else if(component=="wca"){
                        wca->sum(comp,st);
                    } else {
                        Base::sum(component,comp,st);
                    }
                } else {
                    this->sys->template log<System::CRITICAL>("[%s] Requested potential %s to sum. "
                                                                "But %s is not present in the force field",
                                                                Base::name.c_str(),
                                                                component.c_str(),component.c_str());
                }
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                dh->sum(comp,st);
                wca->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                dh->updateBox(box);
                wca->updateBox(box);
            }
    
    };
    
    template<class Base_ >
    class MechanicallyAccurateDNAFast : public MechanicallyAccurateDNABonded<Base_>{

            using Base = MechanicallyAccurateDNABonded<Base_>;
        
        protected:

            const real epsilon = 1.0;
            
            using BondDHType = Potentials::Bond2::DebyeHuckel<typename Base::Units>;
            
            using InteractorBondDHType   = Interactor::BondedInteractor<BondDHType,
                                                                        Interactor::BondedInteractor_ns::BondProcessor<BondDHType>,
                                                                        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,BondDHType>>;
            

            std::shared_ptr<InteractorBondDHType> bonds_dh;
        

        public:
        
            MechanicallyAccurateDNAFast(std::shared_ptr<System>        sys,
                                        std::shared_ptr<ParticleData>  pd,
                                        std::shared_ptr<ParticleGroup> pg,
                                        InputFile&                     in):Base(sys,pd,pg,in){
                
                Base::name = "MechanicallyAccurateDNAFast";

                Base::componentsList.push_back("bonds_dh");
                
                //Add bonds DH
                typename BondDHType::Parameters bondParameters;

                std::shared_ptr<BondDHType> bdh_PBC = std::make_shared<BondDHType>(this->pd,
                                                                                   bondParameters);
                
                typename InteractorBondDHType::Parameters interactorBondParameters;
                
                interactorBondParameters.bondName = "BONDS_DH";

                bonds_dh = std::make_shared<InteractorBondDHType>(this->sys, this->pd, this->pg,
                                                                  this->top, bdh_PBC,
                                                                  interactorBondParameters);
            }


            void sum(std::string component,Computables comp,cudaStream_t st) {
                if(std::find(Base::componentsList.begin(),Base::componentsList.end(),component) != Base::componentsList.end()){
                           if(component=="bonds_dh"){
                        bonds_dh->sum(comp,st);
                    } else {
                        Base::sum(component,comp,st);
                    }
                } else {
                    this->sys->template log<System::CRITICAL>("[%s] Requested potential %s to sum. "
                                                                "But %s is not present in the force field",
                                                                Base::name.c_str(),
                                                                component.c_str(),component.c_str());
                }
            }
            
            void sum(Computables comp,cudaStream_t st) override {
                Base::sum(comp,st);
                bonds_dh->sum(comp,st);
            }
            
            void updateBox(Box box){
                Base::updateBox(box);
                bonds_dh->updateBox(box);
            }
    
    };
}

using MADna = MechanicallyAccurateDNA::MechanicallyAccurateDNA<
              ForceFieldNeighbourBase<UnitsSystem::KCALMOL_A,
                                      Types::BASIC,
                                      conditions::chargedExcluded>>;

using MADnaFast = MechanicallyAccurateDNA::MechanicallyAccurateDNAFast<
                  ForceFieldBase<UnitsSystem::KCALMOL_A,
                                 Types::BASIC>>;

}}}



#endif
