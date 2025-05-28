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
      auto idx = px.size();
      if(tpx.empty()==false){
	//add new components
	px.push_back(tpx[0]);
	py.push_back(tpy[0]);
	pz.push_back(tpz[0]);
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
    template<typename Tp, typename Tm>
      int ParticleMCMatched(const int idx,const RVec<Tp> &tpx,const  RVec<Tp> &tpy,const  RVec<Tp> &tpz, const Tm &tmass, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tp> &m,const RVecI& iafter){
      
      //add new components
      if(tpx.empty()==false){
	std::cout << "Old momenta         " <<idx<<" "<<pz[idx]<<" "<<m[idx]<<" "<<std::endl; 
	px[idx]=tpx[0];
	py[idx]=tpy[0];
	pz[idx]=tpz[0];
	m[idx] = tmass;
	//m[idx] = 0.00051099900;
	
	std::cout<<"ParticleMCMatched "<<idx<<" "<<pz[idx]<<" "<<m[idx]<<" \n"<<std::endl;
	return idx;
      }
      else{
	//std::cout << "Container Empty         " <<idx<<" "<<pz[idx]<<" "<<m[idx]<<" \n"<<std::endl; 
	//px[idx]=0;
	//py[idx]=0;
	//pz[idx]=0;
	//m[idx] = tmass;
	
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
	
	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::ParticleMCMatched(%s,tagger_px,tagger_py,tagger_pz,0.00051099900,%spx,%spy,%spz,%sm,{0})",rad::names::ScatEle().data(),Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());
	
      }
      
      void RomanPotProton() {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	//cant find a way to use this yet
	//float proton_mass=0.93827208943;
	Reaction()->Define(Rec()+"pprime",Form("rad::epic::Particle(rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	//copy rec_pprime to pprime?
	//Reaction()->Define("pprime",Rec()+"pprime");
	Reaction()->AddParticleName(Rec()+"pprime");
      }
      void MCMatchedRomanPotProton() {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	Reaction()->Define(Rec()+"RPproton",Form("rad::epic::ParticleMCMatched(%s,rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})","pprime",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+"RPproton");
      }
      void MCMatchedB0Proton() {
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.x","B0_px");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.y","B0_py");
	Reaction()->setBranchAlias("ReconstructedTruthSeededChargedParticles.momentum.z","B0_pz");
	
	Reaction()->Define(Rec()+"B0proton",Form("rad::epic::ParticleMCMatched(%s,B0_px,B0_py,B0_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})","pprime",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+"B0proton");
      }
      
    };
    
  }
}
