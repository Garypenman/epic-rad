#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order
*/
#include "ElectroIonReaction.h"
#include "ePICUtilities.h"
#include "ReactionUtilities.h"

namespace rad{
  namespace config{
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;
    using rad::reaction::util::DeduceColumnVectorType;
    using rad::reaction::util::ColType;
    using rad::epic::UndoAfterBurn;

    //! Class definition

    class ePICReaction : public ElectroIonReaction {
      
    private:
 
    public:

      ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,fileNameGlob,columns} {

      }
     ePICReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,filenames,columns} {

      }

      /**
       * Only alias ReconstructedParticles columns
       */ 
       void AliasColumns(Bool_t IsEnd=kTRUE){
	AddType(Rec());
	setBranchAlias("ReconstructedParticles.momentum.x",Rec()+"px");
	setBranchAlias("ReconstructedParticles.momentum.y",Rec()+"py");
	setBranchAlias("ReconstructedParticles.momentum.z",Rec()+"pz");
	setBranchAlias("ReconstructedParticles.mass",Rec()+"m");
	setBranchAlias("ReconstructedParticles.PDG",Rec()+"pid");

	reaction::util::CountParticles(this,Rec());
	
	
	if(IsEnd){
	  reaction::util::RedefineFundamentalAliases(this);

	  ApplyAfterBurner(Rec());
	  ApplyAfterBurnerOnBeams();

	  DefineBeamElectron();
	  DefineBeamIon();

	}
       }
      
      /**
       * Alias ReconstructedParticles and MCParticle columns
       */ 
     void AliasColumnsAndMC(Bool_t IsEnd=kTRUE){
	AliasColumns(kFALSE);
	AliasColumnsMC(kFALSE);
	//Matching reconstructed to truth :
	// 1) Find final state truth tru_genStat == 1 =>rec_match_id 0,,N
	// 2) Map from old to new : final_match_id
	//    final_match_id = sizeof(MCParticles)
	//             value = order in new arrays or -1 if not included
	// 3) Create new simID : final_match_id[simID[]]
	//    converts to position in new truth array
	// 4) Create new arrays : sizeof(tru_final_state)
	//    tru_ entries = all filled as truth
	//    rec_ entries = filled if that particle was reconstructed
	
	//remove all but true generated beam+final state particles
	//rec_match_id : 0,1,2,...N=number final particles 
	Define(Rec()+"match_id",[](const ROOT::RVecI& stat){

	  auto filtstat = stat[stat==1];
	  auto id = helpers::Enumerate<uint>(filtstat.size());
	  return id;//[filtstat==1];
	},{Truth()+"genStat"}); //just need tru gen status to get N

	//ReconstructedParticleAssociations.simID is going to be deprecated
	//old files only have that, so here check if it exists, if not use new version
	std::string simID = "ReconstructedParticleAssociations.simID";
	if( rad::config::OriginalColumnExists(simID,CurrFrame()) == false ){
	  simID = "_ReconstructedParticleAssociations_sim.index";
	  //Also need to Define match_id seperately due to tyoe chunge from uint to int
	  //make an mc_match branch cut on actual generated particles (no secondaries)
	  //Points rec array to tru array. rec_array has no beam particles, so can ignore
	  Define(Truth()+"match_id",[](const ROOT::RVecI& stat,const ROOT::RVecI& simID,const ROOT::RVecU& finalID){
	    const auto n = finalID.size(); //mcparticles stat==1
	    ROOT::RVecI final_match_id(n,-1);
	    for(uint i=0;i<n;++i){
	      if(i>=simID.size())break;
	      //if this truth particle was reconstructed add its new id
	      // if(rad::helpers::Contains(simID,finalID[i]))
	      //final_match_id[finalID[i]]=i;

	      if(rad::helpers::Contains(finalID,static_cast<UInt_t>(simID[i]))){
		//final_match_id[finalID[i]]=simID[i]-2;
		final_match_id[i]=rad::helpers::findIndex(finalID,simID[i]);
	      }
	    }
	    
	    
	    ROOT::RVecU tru_match_id =final_match_id[final_match_id!=-1]; //Filter valid ids
	    //std::cout<<"tru_match_id "<<simID<<" "<<simID.size()<<" "<<finalID<<" "<<finalID.size()<<" "<<tru_match_id<<" "<<tru_match_id.size()<<std::endl;
	    return tru_match_id;
	    
	  },{Truth()+"genStat",simID,Truth()+"final_id"});//simID points from rec to tru
	  
	}
	else{
	  //make an mc_match branch cut on actual generated particles (no secondaries)
	  //Points rec array to tru array. rec_array has no beam particles, so can ignore
	  Define(Truth()+"match_id",[](const ROOT::RVecI& stat,const ROOT::RVecU& simID,const ROOT::RVecU& finalID){
	    const auto n = finalID.size(); //mcparticles stat==1
	    ROOT::RVecI final_match_id(n,-1);
	    for(uint i=0;i<n;++i){
	      if(i>=simID.size())break;
	      //if this truth particle was reconstructed add its new id
	      // if(rad::helpers::Contains(simID,finalID[i]))
	      //final_match_id[finalID[i]]=i;

	      if(rad::helpers::Contains(finalID,simID[i])){
		//final_match_id[finalID[i]]=simID[i]-2;
		final_match_id[i]=rad::helpers::findIndex(finalID,simID[i]);
	      }
	    }	  
	    
	    ROOT::RVecU tru_match_id =final_match_id[final_match_id!=-1]; //Filter valid ids
	    //std::cout<<"tru_match_id "<<simID<<" "<<simID.size()<<" "<<finalID<<" "<<finalID.size()<<" "<<tru_match_id<<" "<<tru_match_id.size()<<std::endl;
	    return tru_match_id;
	    
	  },{Truth()+"genStat",simID,Truth()+"final_id"});//simID points from rec to tru, tru_final_id enumerates final state particles
	}

	//make an branch with size of number of generator particles (status 1 or 4)
	//used to truncate tru arrays
	//Define("tru_n","rad::helpers::Count(tru_genStat,1)+rad::helpers::Count(tru_genStat,4)");
	Define(Truth()+"n",Form("rad::helpers::Count(%sgenStat,1)",Truth().data()) );
	
	
	if(IsEnd){
	  reaction::util::RedefineFundamentalAliases(this);


	  ApplyAfterBurner(Truth());
	  ApplyAfterBurner(Rec());
	  ApplyAfterBurnerOnBeams();

	  DefineBeamElectron();
	  DefineBeamIon();

	  //after undo afterburn
	  //AddAdditionalComponents();
	  //needed to make sure tru and rec are defined at the same time
	  //and therefore contain same elements.
	  //	  Filter([](const ROOT::RVecF& th,const ROOT::RVecF& pmag,const ROOT::RVecF& ph,const ROOT::RVecF& eta,const ROOT::RVecF& rth,const ROOT::RVecF& rpmag,const ROOT::RVecF& rph,const ROOT::RVecF& reta){return true;},{Truth()+"pmag",Truth()+"theta",Truth()+"phi",Truth()+"eta",Rec()+"pmag",Rec()+"theta",Rec()+"phi",Rec()+"eta"});
	}
     }
      /**
       * Only alias MCParticle columns
       */ 
      void AliasColumnsMC(Bool_t IsEnd=kTRUE){
	AddType(Truth());
	setBranchAlias("MCParticles.momentum.x",Truth()+"px");
	setBranchAlias("MCParticles.momentum.y",Truth()+"py");
	setBranchAlias("MCParticles.momentum.z",Truth()+"pz");
	setBranchAlias("MCParticles.mass",Truth()+"m");
	setBranchAlias("MCParticles.PDG",Truth()+"pid");
	setBranchAlias("MCParticles.generatorStatus",Truth()+"genStat");
    
	Define(Truth()+"final_id",[](const ROOT::RVecI& stat){
	  //map full array to final state only array  
	  auto indices = helpers::Enumerate<uint>(stat.size());
	  return indices[stat==1];

	},{Truth()+"genStat"});//simID points from rec to tru


	reaction::util::CountParticles(this,Truth());

	if(IsEnd){
	  reaction::util::RedefineFundamentalAliases(this);


	  ApplyAfterBurner(Truth());
	  ApplyAfterBurnerOnBeams();
	  
	  DefineBeamElectron();
	  DefineBeamIon();

	}
      }
      /**
       * Alias the columns and rearrange entries 
       * according to ReconstructedParticleAssociations
       * this reorders reconstructed to match mc
       * Note, we should consider changing to the opposite
       * As this is slower and produces larger output trees
       * than matching to reconstructed order.
       * The advantage this way is that the particles are
       * well defined.
       */

      //this function is very ugly at the moment. This is due to template type requiring
      //explicit type when call the Re* functions.
      void AliasColumnsAndMatchWithMC(Bool_t IsEnd=kTRUE){

	/*Remake this function, probably using string define
	  currently OK for rec, but need to truncate mc columns to tru_n entries*/
	//for tru : RedefineViaAlias(alias,"rad::helpers::Truncate(tru_n))
	//for rec : RedefineViaAlias(alias,Form("helpers::Reorder(%s,rec_match_id,tru_match_id,tru_n)",alias.data());
	
	AliasColumnsAndMC(kFALSE);


	//need to do it here or these calculations are done before synching
	//which will not work as different number of elements
	
	if(IsEnd){
	  reaction::util::RedefineFundamentalAliases(this);

	  //use true masses for each rec track
	  //i.e. ignore PID
	  //Must truncate to make sure return array is same size as in array
	  //RDF may add beams to tru_m before calling this giving a size mismatch
	  RedefineViaAlias(Rec()+"m",[](const ROOT::RVecF& recm,const ROOT::RVecD& trum){return helpers::Truncate(ROOT::RVecF(trum),recm.size());},{Rec()+"m",Truth()+"m"});

	  ApplyAfterBurner(Truth());
	  ApplyAfterBurner(Rec());
	  ApplyAfterBurnerOnBeams();
	  
	  DefineBeamElectron();
	  DefineBeamIon();
	
	  _truthMatched =true;
	}
      }//AliasColumnsAndMatchWithMC

      
      void PostParticles() override{
	//auto finalNames=FinalParticleNames();
	// if(CheckForType(Rec())){

	//   for(const auto& name:finalNames){
	//     Redefine(Rec()+name,Form("rad::indice::InvalidateIndices(%s)",(Rec()+name).data()));
	//   }

	// }
	//once particle are added to vectors
	//we can calculate additional components
	AddAdditionalComponents();

	if(IsTruthMatched()){
	  //needed to make sure tru and rec are defined at the same time
	  //and therefore contain same elements.
	  //need a function string that uses all component vectors and always return true
	  Filter(Form("bool(%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty() +1); ",(Truth()+"pmag").data(),(Truth()+"theta").data(),(Truth()+"phi").data(),(Truth()+"eta").data(),(Rec()+"pmag").data(),(Rec()+"theta").data(),(Rec()+"phi").data(),(Rec()+"eta").data()),"truthmatch");
 
	  //add resolution functions
	  reaction::util::ResolutionFraction(this,"pmag");
	  reaction::util::Resolution(this,"theta");
	  reaction::util::Resolution(this,"phi");
	  reaction::util::Resolution(this,"eta");

	}
      }
      
      void AddAdditionalComponents(){
	//and add some additional columns
	DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
	DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
	DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
	DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));


    }
      
      template<typename T> 
      void RedefineFundamental( const string& name ){
	
	auto contains = [](const std::string&  s1,const std::string& s2){
	  return (s1.find(s2) != std::string::npos);
	};
	
	if(contains(name,"rec") ){
	  RedefineViaAlias(name,helpers::Reorder<T,UInt_t,UInt_t>,{name.data(),Rec()+"match_id",Truth()+"match_id",Truth()+"n"});
	  }
	else if(contains(name,"tru") ){
	  RedefineViaAlias(name,helpers::Rearrange<T,UInt_t>,{name.data(),Truth()+"final_id"});
	  
	}
	
      }
      

      void SetBeamsFromMC(Long64_t nrows=100){
	_useBeamsFromMC=true;
	auto nthreads  = ROOT::GetThreadPoolSize();
	//Range only works in single thread mode
	if(nthreads) ROOT::DisableImplicitMT();

	auto tempframe = GetFileNames().size()==0 ?
	  ROOT::RDataFrame{GetTreeName(),GetFileName()} :
	  ROOT::RDataFrame{GetTreeName(),GetFileNames()};
	
	auto beamdf = tempframe.Range(nrows).Define("emean","MCParticles.momentum.z[MCBeamElectrons_objIdx.index[0]]").Define("pzmean","MCParticles.momentum.z[MCBeamProtons_objIdx.index[0]]").Define("pxmean","MCParticles.momentum.x[MCBeamProtons_objIdx.index[0]]");
	auto pze  = beamdf.Mean("emean");
	auto pzp  = beamdf.Mean("pzmean");
	auto pxp  = beamdf.Mean("pxmean");

	setBeamElectron(0,0,*pze);//afterburned number
	setBeamIon(*pxp,0,*pzp);//afterburned number

	if(nthreads) ROOT::EnableImplicitMT(nthreads);
      }

      void ApplyAfterBurner(std::string type){
	  //Undo the afterburner procedure
	  //here we just account for crossing angle
	  //just need to redefine 1 component. Other 2 have been updated
	  //need to redefine at least one to make sure this is called before
	  //any of the components
	if(DeduceColumnVectorType(this,type + "px")==ColType::Float &&
	   DeduceColumnVectorType(this,type + "m")==ColType::Double){

	  RedefineViaAlias(type+"px",
			   UndoAfterBurn<float,double>{_p4ion_beam,_p4el_beam,-0.025} ,
			   {type+"px",type+"py",type+"pz",type+"m"});
	}
	else if(DeduceColumnVectorType(this,type + "px")==ColType::Double &&
		DeduceColumnVectorType(this,type + "m")==ColType::Double){
	  
	  RedefineViaAlias(type+"px",
			   UndoAfterBurn<double,double>{_p4ion_beam,_p4el_beam,-0.025} ,
			   {type+"px",type+"py",type+"pz",type+"m"});
	}
  	else if(DeduceColumnVectorType(this,type + "px")==ColType::Float &&
		DeduceColumnVectorType(this,type + "m")==ColType::Float){
	  
	  RedefineViaAlias(type+"px",
			   UndoAfterBurn<float,float>{_p4ion_beam,_p4el_beam,-0.025} ,
			   {type+"px",type+"py",type+"pz",type+"m"});
	}
 	else if(DeduceColumnVectorType(this,type + "px")==ColType::Double &&
		DeduceColumnVectorType(this,type + "m")==ColType::Float){
	  
	  RedefineViaAlias(type+"px",
			   UndoAfterBurn<double,float>{_p4ion_beam,_p4el_beam,-0.025} ,
			   {type+"px",type+"py",type+"pz",type+"m"});
	}
	else{
	  std::cerr <<"Error ePICReaction::ApplyAfterBurner momentum and mass types not valid "<<CurrFrame().GetColumnType(type + "px")<<" "<<CurrFrame().GetColumnType(type + "m")<<std::endl;
	  exit(0);
	}
      }
      void ApplyAfterBurnerOnBeams(){

	//undo afterburn on beam components
	auto beams_px = ROOT::RVecD{_p4el_beam.X(),_p4ion_beam.X()};
	auto beams_py = ROOT::RVecD{_p4el_beam.Y(),_p4ion_beam.Y()};
	auto beams_pz = ROOT::RVecD{_p4el_beam.Z(),_p4ion_beam.Z()};
	auto beams_m = ROOT::RVecD{_p4el_beam.M(),_p4ion_beam.M()};
	cout<<"Pre Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	rad::epic::UndoAfterBurn<double,double> undoAB_DD{_p4ion_beam,_p4el_beam,-0.025};
	undoAB_DD(beams_px,beams_py,beams_pz,beams_m);
	_p4el_beam.SetCoordinates(beams_px[0],beams_py[0],beams_pz[0],beams_m[0]);
	_p4ion_beam.SetCoordinates(beams_px[1],beams_py[1],beams_pz[1],beams_m[1]);
	
	cout<<"Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	
      }

      bool IsTruthMatched()const {return _truthMatched;}

    private:
      
      bool _truthMatched =false;

    };//ePICReaction

    
  }//config
}//rad
