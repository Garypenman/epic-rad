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
      int Particle(const RVec<Tp> &tpx,const  RVec<Tp> &tpy,const  RVec<Tp> &tpz, const Tm &tmass, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tp> &m,const RVecI& iafter){
      
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
      //add new components

      UInt_t entry = 0;
      if(tpx.empty()==false){
	//check threshold
	if((tpx[entry]*tpx[entry]+tpy[entry]*tpy[entry]+tpz[entry]*tpz[entry])<threshold*threshold){
	  return -1;
	}
	//if(tpz[0]<0)std::cout<<"ParticleMCMatched "<<m<<" "<<tpz<<" "<<idx<<" "<<(tpx[0]*tpx[0]+tpy[0]*tpy[0]+tpz[0]*tpz[0])<<" "<<tmass<<std::endl;
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
      
      ePICParticleCreator() = default;
    ePICParticleCreator(rad::config::ConfigReaction& cr):rad::config::ParticleCreator{cr}{};
      
      //////////////////////////////////////////////////////////////////
      void LowQ2Electron(/*const string& name,const string& p4name*/) {
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.x","tagger_px");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.y","tagger_py");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.z","tagger_pz");
	//empty parts string as not dependent on others

	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::Particle(tagger_px,tagger_py,tagger_pz,0.00051099900,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	//g. penman 21.05 this redefines scat_ele and doesnt work with mcmatching
	//copy rec_scat_ele to scat_ele
	//Reaction()->Define(rad::names::ScatEle(),Rec()+rad::names::ScatEle());
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());

      }
 
      void MCMatchedLowQ2Electron() {
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.x","tagger_px");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.y","tagger_py");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.z","tagger_pz");
	
	//cant find a way to use this yet
	//float electron_mass=0.00051099900;

	//Note threshold = 0.1
	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::ParticleMCMatched(0.1,%s,tagger_px,tagger_py,tagger_pz,0.00051099900,%spx,%spy,%spz,%sm,%s)",rad::names::ScatEle().data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());
	
      }
      
      void RomanPotProton(const std::string name="pprime") {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	//cant find a way to use this yet
	//float proton_mass=0.93827208943;
	Reaction()->Define(Rec()+name,Form("rad::epic::Particle(rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,%s)",Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	//copy rec_pprime to pprime?
	//Reaction()->Define(name,Rec()+name);
	Reaction()->AddParticleName(Rec()+name);
      }
      void MCMatchedRomanPotProton(const std::string name="pprime") {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	//Note threshold = 10
	Reaction()->Define(Rec()+"RPproton",Form("rad::epic::ParticleMCMatched(10,%s,rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,%s)",name.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+"RPproton");
      }
      void MCMatchedB0Proton(const std::string name="pprime") {
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.x","B0_px");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.y","B0_py");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.z","B0_pz");
	
	//Note threshold = 10
	Reaction()->Define(Rec()+"B0proton",Form("rad::epic::ParticleMCMatched(10,%s,B0_px,B0_py,B0_pz,0.93827208943,%spx,%spy,%spz,%sm,%s)",name.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+"B0proton");
      }
      
      void MCMatchedFarForwardProton(const std::string name="pprime") {
	MCMatchedRomanPotProton(name);
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.x","B0_px");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.y","B0_py");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.z","B0_pz");

	//if RP exists use that
	//if not consider B0 candidates
	//Note B0 threshold = 10
	Reaction()->Define(Rec()+"B0proton",Form("if(rec_RPproton==-1) return rad::epic::ParticleMCMatched(10,%s,B0_px,B0_py,B0_pz,0.93827208943,%spx,%spy,%spz,%sm,%s); return -1;",name.data(),Rec().data(),Rec().data(),Rec().data(),Rec().data(),(Truth()+"match_id").data()));
	Reaction()->AddParticleName(Rec()+"B0proton");
     }
   };
    
  }
}
