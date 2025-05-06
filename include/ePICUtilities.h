#pragma once
#include "Beams.h"
#include <Math/VectorUtil.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>

//inline  std::array<double,4>  rad::beams::InitBotComponents() {return {0,0,99.9339,0.938272};}
//inline  std::array<double,4>  rad::beams::InitTopComponents() {return {0,0,-10.007,0.000510999};}


namespace rad{
  namespace epic{

    using PxPyPzEVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
    using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;
    using ROOT::Math::RotationX;
    using ROOT::Math::RotationY;
    using ROOT::RVecD;
    using ROOT::RVecF;
    
 
    
    //UndoAfterBurn undo{xangle};
    class UndoAfterBurn 
    {
    public:
      UndoAfterBurn(PxPyPzMVector p_beam,PxPyPzMVector e_beam,Float_t angle=-0.025):_crossAngle{angle}{
	//calculate and store current boosts and rotations
	RotsAndBoosts(p_beam,e_beam);
      };

      //can't template as foreach requires compiletime knowledge
      template<typename Tp, typename Tm>
      // void operator()(RVec<Tp> &px,RVec<Tp> &py,RVec<Tp> &pz, const RVec<Tm> &m) const
      bool operator()(RVec<Tp> &px,RVec<Tp> &py,RVec<Tp> &pz, const RVec<Tm> &m) const
      {
	
	//apply to all particles
	auto n_parts = m.size();
	for(uint i=0;i<n_parts;++i){
	  undoAfterburn(i,px,py,pz,m);
	}
	return true;
	//return px;
      }
      //void  operator()(RVecF &px,RVecF &py,RVecF &pz, const RVecF &m) const
      //void operator()(unsigned int slot, RVecF &px,RVecF &py,RVecF &pz, const RVecD//  &m) const
      // {
	
      // 	//apply to all particles
      // 	auto n_parts = m.size();
      // 	for(auto i=0;i<n_parts;++i){
      // 	  undoAfterburn(i,px,py,pz,m);
      // 	}
      // 	//return px;
      // }

    private:
      void RotsAndBoosts(PxPyPzMVector p_beam,PxPyPzMVector e_beam);

      template<typename Tp, typename Tm>
      void undoAfterburn(uint idx,RVec<Tp> &px,RVec<Tp> &py,RVec<Tp> &pz, const RVec<Tm> &m) const;
      //  void undoAfterburn(uint idx,RVecF &px,RVecF &py,RVecF &pz, const RVecD &m) const;

      // Objects for undoing afterburn boost
      Float_t _crossAngle{-0.025}; // Crossing angle in radians
      RotationX _rotAboutX;
      RotationY _rotAboutY;
      MomVector _vBoostToCoM;
      MomVector _vBoostToHoF;

    };

    
    //----------------------------------------------------
    //----------------------------------------------------
    //            UNDO AFTERBURNER PROCEDURE
    //----------------------------------------------------
    //----------------------------------------------------
    
    // Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
    // USE BEAM VECTORS
    void UndoAfterBurn::RotsAndBoosts(PxPyPzMVector p_beam,PxPyPzMVector e_beam){
       //We need MCParticle beams ?
      //   auto p_beam = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      //auto e_beam = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
    
      // auto p_beam = rad::beams::InitialBotVector();
      //auto e_beam = rad::beams::InitialTopVector();
      
      // Holding vectors for beam - undoing crossing angle ONLY
      //PxPyPzMVector p_beam(_crossAngle*p.E(), 0., p.E(), p.M());
      //PxPyPzMVector e_beam(0., 0., -k.E(), k.M());
      p_beam.SetCoordinates(_crossAngle*p_beam.E(), 0., p_beam.E(), p_beam.M());
      e_beam.SetCoordinates(0., 0., -e_beam.E(), e_beam.M());
     
      // Define boost vector to CoM frame
      auto CoM_boost = p_beam+e_beam;
      //vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
      _vBoostToCoM =  CoM_boost.BoostToCM();
      
      // Apply boost to beam vectors
      p_beam = boost(p_beam, _vBoostToCoM);
      e_beam = boost(e_beam, _vBoostToCoM);
  
      // Calculate rotation angles and create rotation objects
      auto rotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
      auto rotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

      _rotAboutY = ROOT::Math::RotationY(rotY);
      _rotAboutX = ROOT::Math::RotationX(rotX);

      // Apply rotation to beam vectors
      p_beam = _rotAboutY(p_beam);
      p_beam = _rotAboutX(p_beam);
      e_beam = _rotAboutY(e_beam);
      e_beam = _rotAboutX(e_beam);

      // Define boost vector back to head-on frame
      PxPyPzEVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
      _vBoostToHoF = -HoF_boost.BoostToCM();

      
      //vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

      // Apply boost back to head on frame to beam vectors
      p_beam = boost(p_beam, _vBoostToHoF);
      e_beam = boost(e_beam, _vBoostToHoF);
      //      std::cout<<"UndoAfterBurn::RotsAndBoosts() "<<p_beam<<" "<<e_beam<<std::endl;
      // Make changes to input vectors
      //p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
      // k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
    }
    //////////////////////////////////////////////////////
    // Undo afterburn procedure only
    
    //void UndoAfterBurn::undoAfterburn(uint idx,RVecF &px,RVecF &py,RVecF&pz, const RVecD &m) const{
    template<typename Tp, typename Tm>
    void UndoAfterBurn::undoAfterburn(uint idx,RVec<Tp> &px,RVec<Tp> &py,RVec<Tp> &pz, const RVec<Tm> &m) const{
      //   std::cout<<"undoAfterburn in "<<idx<<" px "<<px[idx]<<" m "<<m[idx]<<std::endl;
      auto a = FourVector(idx,px,py,pz,m);
      //std::cout<<"undoAfterburn in "<<a<<_rotAboutY<<_rotAboutX<<" "<<_vBoostToCoM<<_vBoostToHoF<<std::endl;
      // Undo AB procedure for single vector, a^{mu}
      a = boost(a, _vBoostToCoM); // BOOST TO COM FRAME
      a = _rotAboutY(a);          // ROTATE TO Z-AXIS
      a = _rotAboutX(a);          // ROTATE TO Z-AXIS
      a = boost(a, _vBoostToHoF); // BOOST BACK TO HEAD ON FRAME
      
      //apply to original components
      px[idx]=a.X();
      py[idx]=a.Y();
      pz[idx]=a.z();
    }

  }
}
