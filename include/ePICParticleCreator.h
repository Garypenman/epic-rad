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
    int ParticleLowQ2Electron(const RVec<Tp> &tpx,const  RVec<Tp> &tpy,const  RVec<Tp> &tpz, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){

      //std::cout<<"ParticleFixedBeam "<< px.size()<<m<<m.size()<<std::endl;
      auto idx = px.size();
      if(tpx.empty()==false){
	//add new components
	px.push_back(tpx[0]);
	py.push_back(tpy[0]);
	pz.push_back(tpz[0]);
	m.push_back(0.00051099900);
      }
      else{
	px.push_back(0.);
	py.push_back(0.);
	pz.push_back(0.);
	m.push_back(0.00051099900);
    
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
	px[idx]=tpx[0];
	py[idx]=tpy[0];
	pz[idx]=tpz[0];
	m[idx] = tmass;
	//m[idx] = 0.00051099900;
	
	//std::cout<<"ParticleMCMatched "<<idx<<" "<<pz[idx]<<" "<<m[idx]<<" "<<std::endl;
	return idx;
      }
      return -1;
    }
    
    class ePICParticleCreator : public rad::config::ParticleCreator{
      
    public:
      
      ePICParticleCreator() = default;
    ePICParticleCreator(rad::config::ConfigReaction& cr):rad::config::ParticleCreator{cr}{};

      //////////////////////////////////////////////////////////////////
      void LowQ2Electron(const string& name,const string& p4name) {
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.x","tagger_px");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.y","tagger_py");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.z","tagger_pz");
	//empty parts string as not dependent on others

	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::ParticleLowQ2Electron(tagger_px,tagger_py,tagger_pz,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	//copy rec_scat_ele to scat_ele
	Reaction()->Define(rad::names::ScatEle(),Rec()+rad::names::ScatEle());
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
      
      void MCMatchedRomanPotProton() {
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.x","rp_px");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.y","rp_py");
	Reaction()->setBranchAlias("ForwardRomanPotRecParticles.momentum.z","rp_pz");
	
	//cant find a way to use this yet
	//float proton_mass=0.93827208943;
	//is rad::names::Baryons() the correct name to use here? 
	//no not in current form because its a string_view. Need to change how its initialised in rad/include/DefineNames.h
	//or use different approach. Will hard code a string for  now and consult derek later.
	//Reaction()->Define(Rec()+rad::names::Baryons(),Form("rad::epic::ParticleMCMatched(%s,rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})",rad::names::Baryons().data(),Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	//Reaction()->AddParticleName(Rec()+rad::names::Baryons());
	//Reaction()->Define(Rec()+rad::names::Baryons(),Form("rad::epic::ParticleMCMatched(%s,rp_px,rp_py,rp_pz,0.93827208943,%spx,%spy,%spz,%sm,{0})",rad::names::Baryons().data(),Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	//Reaction()->AddParticleName(Rec()+rad::names::Baryons());
	std::cout << Rec() << " " << Rec().data() << std::endl;
      }
      
    };
    
  }
}
