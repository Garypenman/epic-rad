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
      int ParticleMCMatchedLowQ2Electron(const int idx,const RVec<Tp> &tpx,const  RVec<Tp> &tpy,const  RVec<Tp> &tpz, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz,RVec<Tp> &trupz, RVec<Tm> &m,const RVecI& iafter){

      //add new components
      if(tpx.empty()==false){
	px[idx]=tpx[0];
	py[idx]=tpy[0];
	pz[idx]=tpz[0];
	m[idx] = 0.00051099900;
	
	//std::cout<<"ParticleMCMatchedLowQ2Electron "<<idx<<" "<<pz[idx]<<" "<<trupz[idx]<<" "<<tpz[0]<<std::endl;
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

	Reaction()->Define(Rec()+"tagger",Form("rad::epic::ParticleLowQ2Electron(tagger_px,tagger_py,tagger_pz,%spx,%spy,%spz,%sm,{0})",Rec().data(),Rec().data(),Rec().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());

      }
 
      void MCMatchedLowQ2Electron() {
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.x","tagger_px");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.y","tagger_py");
	Reaction()->setBranchAlias("TaggerTrackerTracks.momentum.z","tagger_pz");

	Reaction()->Define(Rec()+rad::names::ScatEle(),Form("rad::epic::ParticleMCMatchedLowQ2Electron(%s,tagger_px,tagger_py,tagger_pz,%spx,%spy,%spz,%spz,%sm,{0})",rad::names::ScatEle().data(),Rec().data(),Rec().data(),Rec().data(),Truth().data(),Rec().data()));
	Reaction()->AddParticleName(Rec()+rad::names::ScatEle());
	
      }
       
    };
    
  }
}
