#pragma once

#include "ParticleCreator.h"
#include "DefineNames.h"

namespace rad{
  namespace epic{
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;

    ///\brief Add scattered e- from tagger
    ///to the particle 4-vector lists
    /// p4 is the beam particle vector
    /// iafter is to keep the prototype of other particle adding
    /// functions which depend on other particles
    template<typename Tp, typename Tm>
      int Particle(const RVec<Tp> &tpx, const  RVec<Tp> &tpy, const  RVec<Tp> &tpz, const Tm &tmass, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tp> &m,const RVecI& iafter){
      
      //std::cout<<"ParticleFixedBeam "<< px.size()<<m<<m.size()<<std::endl;
      UInt_t entry = 0;
      auto idx = px.size();
      if(tpx.empty()==false){
	//add new components
	px.push_back(tpx[entry]);
	py.push_back(tpy[entry]);
	pz.push_back(tpz[entry]);
	m.push_back(tmass);
	//m.push_back(0.00051099900);
      }
      else{
	px.push_back(0.);
	py.push_back(0.);
	pz.push_back(0.);
	m.push_back(0.);
	//m.push_back(0.00051099900);
      }
      return idx;
    }
    ///\brief Place scattered e- from tagger
    ///to the particle 4-vector lists
    ///synched with the tru_ scattered e-
    /// iafter is to keep the prototype of other particle adding
    /// functions which depend on other particles
    template<typename Tp, typename Tm,typename Tmatch>
    int ParticleMCMatched(const double threshold,const int idx,const RVec<Tp> &tpx,const  RVec<Tp> &tpy,const  RVec<Tp> &tpz, const Tm &tmass, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tp> &m,RVec<Tmatch>& imatch){
      
      //check for valid mc particle
      if(idx==-1) return -1;

      UInt_t entry = 0;
      if(tpx.empty()==false){
	//check threshold
	if((tpx[entry]*tpx[entry]+tpy[entry]*tpy[entry]+tpz[entry]*tpz[entry])<threshold*threshold){
	  return -1;
	}

	//add new components
  	px[idx]=tpx[entry];
	py[idx]=tpy[entry];
	pz[idx]=tpz[entry];
	m[idx] = tmass;
	//Add to Truth()+"match_id";
	imatch.push_back(idx);
	return idx;
      }
      else{	
	return -1;
      }
    }
    
    class ePICParticleCreator : public rad::config::ParticleCreator{
      
    public:
      
    ePICParticleCreator() : mcmatched_forward_proton(false) {};
    ePICParticleCreator(rad::config::ConfigReaction& cr) : rad::config::ParticleCreator{cr}{}, mcmatched_forward_proton(false);
      
      //Detector functions (not mc matched)
      void LowQ2Electron() {
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.x","tagger_px");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.y","tagger_py");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.z","tagger_pz");
	
	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::Particle(tagger_px,tagger_py,tagger_pz,0.00051099900,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());

      }
      
      void RomanPotProton(const std::string name="pprime") {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	Reaction()->Define(Rec()+name,Form("rad::epic::Particle(rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,%s)",Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+name);
      }
      
      void B0Proton(const std::string name="pprime"){
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.x","B0_px");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.y","B0_py");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.z","B0_pz");
	
	Reaction()->Define(Rec()+"B0proton",Form("rad::epic::Particle(B0_px,B0_py,B0_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+"B0proton");
      }
      
      void ZDCLambda(const std::string name="lambda"){
	Reaction()->setBranchAlias("ReconstructedFarForwardZDCLambdas.momentum.x","ZDC_px");
	Reaction()->setBranchAlias("ReconstructedFarForwardZDCLambdas.momentum.y","ZDC_py");
	Reaction()->setBranchAlias("ReconstructedFarForwardZDCLambdas.momentum.z","ZDC_pz");
	
	Reaction()->Define(Rec()+"ZDClambda",Form("rad::epic::Particle(ZDC_px,ZDC_py,ZDC_pz,1.1115683,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+"ZDClambda");

      }
      
      //to do - swap above functions to call this one
      //and check for bugs
      void DetectorParticle(const std::string name, const std::string branchname, const std::string aliasname, const double mass){
	std::string smass = to_string(mass);
	
	Reaction()->setBranchAlias(branchname+".momentum.x",aliasname+"_px");
	Reaction()->setBranchAlias(branchname+".momentum.y",aliasname+"_py");
	Reaction()->setBranchAlias(branchname+".momentum.z",aliasname+"_pz");
	
	Reaction()->Define(Rec()+aliasname,Form("rad::epic::Particle(%s_px,%s_py,%s_pz,%s,%spx,%spy,%spz,%sm,{0}",aliasname.data(),aliasname.data(),aliasname.data(),smass.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+aliasname);
	
      }


      /////////////////////
      //MC Matched Functions
      void MCMatchedParticle(const std::string name, const std::string matchname, const std::string branchname, const std::string aliasname, const double mass, const double thresh){
	
	std::string smass = to_string(mass);
	std::string sthresh = to_string(thresh);
	
	Reaction()->setBranchAlias(branchname+".momentum.x",aliasname+"_px");
	Reaction()->setBranchAlias(branchname+".momentum.y",aliasname+"_py");
	Reaction()->setBranchAlias(branchname+".momentum.z",aliasname+"_pz");
	
	//to do - string utility to replace form
	//make function
	Reaction()->Define(Rec()+name,Form("rad::epic::ParticleMCMatched(%s,%s,%s_px,%s_py,%s_pz,%s,%spx,%spy,%spz,%sm,%s)",sthresh.data(),matchname.data(),aliasname.data(),aliasname.data(),aliasname.data(),smass.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+name);
      }
      
      void MCMatchedLowQ2Electron() {
	MCMatchedParticle(rad::names::ScatEle(), rad::names::ScatEle(), "TaggerTrackerTracks", "tagger", rad::constant::M_ele(), 0.1);
      }
      
      void MCMatchedRomanPotProton(const std::string matchname="pprime") {
	MCMatchedParticle("RPproton", matchname, "ForwardRomanPotRecParticles", "rp", rad::constant::M_pro(), 10);
      }
      
      void MCMatchedB0Proton(const std::string matchname="pprime") {
	MCMatchedParticle("B0proton", matchname, "ReconstructedTruthSeededChargedParticles", "B0", rad::constant::M_pro(), 10);
      }
      
      void MCMatchedZDCLambda(const std::string matchname="lambda"){
	MCMatchedParticle("ZDClambda", matchname, "ReconstructedFarForwardZDCLambdas", "ZDC", rad::constant::M_Lambda(), 10);
      }
      
      void MCMatchedFarForwardProton(const std::string name="pprime") {
	MCMatchedRomanPotProton(name);
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.x","B0_px");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.y","B0_py");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.z","B0_pz");
	
	//to do - make this work with the new MCMatchedParticle intermediate function style above
	//if RP exists use that
	//if not consider B0 candidates
	//Note B0 threshold = 10
	Reaction()->Define(Rec()+"B0proton",Form("if(rec_RPproton==-1) return rad::epic::ParticleMCMatched(10,%s,B0_px,B0_py,B0_pz,0.93827208943,%spx,%spy,%spz,%sm,%s); return -1;",name.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+"B0proton");
      }
      
      //probably redundant
      //as this approach wont/cant work per event
      bool getMCMatchedForwardProtonStatus(){
	return mcmatched_forward_proton;
      }
      
    private:
      bool mcmatched_forward_proton;
  
    };//end class
  }//end epic NS
}//end rad NS
