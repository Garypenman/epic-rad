#pragma once
#include "Beams.h"
#include <Math/VectorUtil.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>

namespace rad {
  namespace epic {

    using PxPyPzEVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
    using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;
    using ROOT::Math::RotationX;
    using ROOT::Math::RotationY;
    using ROOT::RVecD;
    using ROOT::RVecF;

    /**
     * @brief Class to undo afterburner transformations on particle momenta.
     * 
     * This class computes and stores the required boosts and rotations to 
     * undo the crossing angle and put the event in the desired reference frame.
     */
    template<typename Tp, typename Tm>
    class UndoAfterBurn {
      using momeType_t = Tp;
      using massType_t = Tm;

    public:
      /**
       * @brief Constructor. Initializes the undo parameters using beam vectors and crossing angle.
       */
      UndoAfterBurn(PxPyPzMVector p_beam, PxPyPzMVector e_beam, Float_t angle = -0.025)
        : _crossAngle{angle} {
        RotsAndBoosts(p_beam, e_beam);
      }

      /**
       * @brief Applies the undo transformation to all particles in the event.
       * 
       * @param px Vector of x components of momentum.
       * @param py Vector of y components of momentum.
       * @param pz Vector of z components of momentum.
       * @param m Vector of masses.
       * @return Transformed px vector (for chaining if desired).
       */
      RVec<momeType_t> operator()(RVec<momeType_t>& px, RVec<momeType_t>& py, RVec<momeType_t>& pz, const RVec<massType_t>& m) const {
        auto n_parts = m.size();
        for (size_t i = 0; i < n_parts; ++i) {
          undoAfterburn(i, px, py, pz, m);
        }
        return px;
      }

    private:
      Float_t _crossAngle{-0.025};          // Crossing angle in radians
      RotationX _rotAboutX;                 // Rotation about X axis
      RotationY _rotAboutY;                 // Rotation about Y axis
      MomVector _vBoostToCoM;               // Boost vector to CoM frame
      MomVector _vBoostToHoF;               // Boost vector back to head-on frame

      /**
       * @brief Calculates and stores the required boosts and rotations based on beam vectors.
       * 
       * This must be called at construction for each new event to set up the transformations.
       * Requires beam 4-vectors as input
       */
      void RotsAndBoosts(PxPyPzMVector p_beam, PxPyPzMVector e_beam);
      

      /**
       * @brief Undo afterburner transformation for a single particle.
       * 
       * @param idx Index of the particle in the vectors.
       * @param px Vector of x momentum.
       * @param py Vector of y momentum.
       * @param pz Vector of z momentum.
       * @param m Vector of masses.
       */
      void undoAfterburn(size_t idx, RVec<momeType_t>& px, RVec<momeType_t>& py, RVec<momeType_t>& pz, const RVec<massType_t>& m) const {
        auto a = FourVector(idx, px, py, pz, m);
        // Step 1: Boost to CoM frame
        a = boost(a, _vBoostToCoM);
        // Step 2: Rotate to align with Z axis
        a = _rotAboutY(a);
        a = _rotAboutX(a);
        // Step 3: Boost back to head-on frame
        a = boost(a, _vBoostToHoF);

        // Update momentum components
        px[idx] = a.X();
        py[idx] = a.Y();
        pz[idx] = a.Z();
      }
    };


    template<typename Tp, typename Tm>
    inline void UndoAfterBurn<Tp, Tm>::RotsAndBoosts(PxPyPzMVector p_beam, PxPyPzMVector e_beam) {
        // Set beam coordinates with crossing angle and energy
        p_beam.SetCoordinates(_crossAngle * p_beam.E(), 0., p_beam.E(), p_beam.M());
        e_beam.SetCoordinates(0., 0., -e_beam.E(), e_beam.M());

        // Calculate boost to center-of-mass (CoM) frame
        auto CoM_boost = p_beam + e_beam;
        _vBoostToCoM = CoM_boost.BoostToCM();

        // Boost beams to CoM frame
        p_beam = boost(p_beam, _vBoostToCoM);
        e_beam = boost(e_beam, _vBoostToCoM);

        // Calculate required rotation angles
        auto rotY = -1.0 * TMath::ATan2(p_beam.X(), p_beam.Z());
        auto rotX = 1.0 * TMath::ATan2(p_beam.Y(), p_beam.Z());

        // Set up rotation objects
        _rotAboutY = ROOT::Math::RotationY(rotY);
        _rotAboutX = ROOT::Math::RotationX(rotX);

        // Rotate beams to align with axes
        p_beam = _rotAboutY(p_beam);
        p_beam = _rotAboutX(p_beam);
        e_beam = _rotAboutY(e_beam);
        e_beam = _rotAboutX(e_beam);

        // Calculate boost back to head-on frame
        PxPyPzEVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
        _vBoostToHoF = -HoF_boost.BoostToCM();

        // Apply final boost to beams (not strictly necessary for all workflows)
        p_beam = boost(p_beam, _vBoostToHoF);
        e_beam = boost(e_beam, _vBoostToHoF);
     }
  } // namespace epic
} // namespace rad
