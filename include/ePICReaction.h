#pragma once

#include "ElectroIonReaction.h"
#include "ePICParticleCreator.h"
#include "ePICUtilities.h"
#include "ReactionUtilities.h"
#include "StringUtilities.h"

/**
 * @file ePICReaction.h
 * @brief Defines the rad::config::ePICReaction class for configuring ePIC ROOT files.
 */

namespace rad {
  namespace config {
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;
    using rad::reaction::util::DeduceColumnVectorType;
    using rad::reaction::util::ColType;
    using rad::epic::UndoAfterBurn;
    using rad::utils::createFunctionCallString;
 
    /**
     * @class ePICReaction
     * @brief Configures ePIC ROOT files for hadronic final state data analysis.
     *
     * Derived from ElectroIonReaction, this class provides methods for setting up
     * data analysis and calculations for ePIC files with fixed particle order.
     */
    class ePICReaction : public ElectroIonReaction {
    private:

    public:
      /**
       * @brief Construct an ePICReaction for a set of files by glob.
       * @param treeName Name of the ROOT tree.
       * @param fileNameGlob Glob pattern for input files.
       * @param columns Optional list of columns to load.
       */
      ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns = {});

      /**
       * @brief Construct an ePICReaction from an explicit list of files.
       * @param treeName Name of the ROOT tree.
       * @param filenames Vector of input file names.
       * @param columns Optional list of columns to load.
       */
      ePICReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns = {});

      /**
       * @brief Alias only ReconstructedParticles columns.
       * @param IsEnd Whether to finalize aliasing and apply afterburners.
       */
      void AliasColumns(Bool_t IsEnd = kTRUE);

      /**
       * @brief Alias both ReconstructedParticles and MCParticle columns and perform matching.
       * @param IsEnd Whether to finalize aliasing and apply afterburners.
       */
      void AliasColumnsAndMC(Bool_t IsEnd = kTRUE);

      /**
       * @brief Alias only MCParticle columns.
       * @param IsEnd Whether to finalize aliasing and apply afterburners.
       */
      void AliasColumnsMC(Bool_t IsEnd = kTRUE);

      /**
       * @brief Alias and reorder columns so reconstructed particles match MC order.
       * @param IsEnd Whether to finalize aliasing and apply afterburners.
       */
      void AliasColumnsAndMatchWithMC(Bool_t IsEnd = kTRUE);

      /**
       * @brief Final processing after particle columns are defined.
       *
       * Adds additional components and, if truth-matching is enabled, defines resolution columns.
       */
      void PostParticles() override;

      /**
       * @brief Add additional kinematic components (phi, theta, eta, pmag) for all types.
       */
      void AddAdditionalComponents();

      /**
       * @brief Redefine fundamental columns to match truth or reconstructed ordering.
       * @tparam T Column type.
       * @param name Name of the column to redefine.
       */
      template<typename T>
      void RedefineFundamental(const std::string& name);

      /**
       * @brief Set beam momenta from MCParticles for a sample of rows.
       * @param nrows Number of rows to use for mean calculation.
       */
      void SetBeamsFromMC(Long64_t nrows = 100);
     /**
       * @brief Set beam momenta from MCParticles for a sample of rows.
       * @param iel index of electron beam
       * @param iion index of ion beam.
       * @param nrows Number of rows to use for mean calculation.
       */
      void SetBeamsFromMC(UInt_t iel,UInt_t iion,Long64_t nrows = 100);
      /**
       * @brief Undo the afterburner procedure for a given particle type.
       * @param type Particle type prefix, e.g., "rec" or "tru".
       */
      void ApplyAfterBurner(std::string type);

      /**
       * @brief Undo the afterburner procedure for the beam particles.
       */
      void ApplyAfterBurnerOnBeams();

      /**
       * @brief Check if columns have been matched to truth particles.
       * @return True if truth matching has been performed.
       */
      bool IsTruthMatched() const;

    private:
      bool _truthMatched = false; ///< True if truth matching has been performed.
    };


    // ---------------- Inline function definitions ----------------

    inline ePICReaction::ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns)
      : ElectroIonReaction{treeName, fileNameGlob, columns} {
    }

    inline ePICReaction::ePICReaction(const std::string_view treeName, const std::vector<std::string>& filenames, const ROOT::RDF::ColumnNames_t& columns)
      : ElectroIonReaction{treeName, filenames, columns} {
    }

    inline void ePICReaction::AliasColumns(Bool_t IsEnd) {
      AddType(Rec());
      setBranchAlias("ReconstructedParticles.momentum.x", Rec() + "px");
      setBranchAlias("ReconstructedParticles.momentum.y", Rec() + "py");
      setBranchAlias("ReconstructedParticles.momentum.z", Rec() + "pz");
      setBranchAlias("ReconstructedParticles.mass", Rec() + "m");
      setBranchAlias("ReconstructedParticles.PDG", Rec() + "pid");

      reaction::util::CountParticles(this, Rec());
      Define(Rec() + "n", Form("%sm.size()", Rec().data()));

      if (IsEnd) {
        reaction::util::RedefineFundamentalAliases(this);
        ApplyAfterBurner(Rec());
        ApplyAfterBurnerOnBeams();
        DefineBeamElectron();
        DefineBeamIon();
      }
    }

    inline void ePICReaction::AliasColumnsAndMC(Bool_t IsEnd) {
      AliasColumns(kFALSE);
      AliasColumnsMC(kFALSE);

      // Remove all but true generated beam+final state particles and perform matching.
      Define(Rec() + "match_id", [](const ROOT::RVecI& stat) {
        auto filtstat = stat[stat == 1];
        auto id = helpers::Enumerate<uint>(filtstat.size());
        return id;
      }, {Truth() + "genStat"});

      // Determine simID column for matching reconstructed and truth particles.
      std::string simID = "ReconstructedParticleAssociations.simID";
      if (OriginalColumnExists(simID) == false) {
        simID = "_ReconstructedParticleAssociations_sim.index";
        Define(Truth() + "match_id", [](const ROOT::RVecI& stat, const ROOT::RVecI& simID, const ROOT::RVecU& finalID) {
	  const auto n = finalID.size();
	  ROOT::RVecI final_match_id(n, -1);
	  for (uint i = 0; i < n; ++i) {
	    if (i >= simID.size()) break;
	    if (rad::helpers::Contains(finalID, static_cast<UInt_t>(simID[i]))) {
	      final_match_id[i] = rad::helpers::findIndex(finalID, simID[i]);
	    }
	  }
	  ROOT::RVecU tru_match_id = final_match_id[final_match_id != -1];
	  return tru_match_id;
        }, {Truth() + "genStat", simID, Truth() + "final_id"});
      } else {
        Define(Truth() + "match_id", [](const ROOT::RVecI& stat, const ROOT::RVecU& simID, const ROOT::RVecU& finalID) {
	  const auto n = finalID.size();
	  ROOT::RVecI final_match_id(n, -1);
	  for (uint i = 0; i < n; ++i) {
	    if (i >= simID.size()) break;
	    if (rad::helpers::Contains(finalID, simID[i])) {
	      final_match_id[i] = rad::helpers::findIndex(finalID, simID[i]);
	    }
	  }
	  ROOT::RVecU tru_match_id = final_match_id[final_match_id != -1];
	  return tru_match_id;
        }, {Truth() + "genStat", simID, Truth() + "final_id"});
      }

      // Define the number of truth particles.
      Define(Truth() + "n", createFunctionCallString("rad::helpers::Count", Truth()+"genStat","1"));

      if (IsEnd) {
        reaction::util::RedefineFundamentalAliases(this);
        ApplyAfterBurner(Truth());
        ApplyAfterBurner(Rec());
        ApplyAfterBurnerOnBeams();
        DefineBeamElectron();
        DefineBeamIon();
      }
    }

    inline void ePICReaction::AliasColumnsMC(Bool_t IsEnd) {
      AddType(Truth());
      setBranchAlias("MCParticles.momentum.x", Truth() + "px");
      setBranchAlias("MCParticles.momentum.y", Truth() + "py");
      setBranchAlias("MCParticles.momentum.z", Truth() + "pz");
      setBranchAlias("MCParticles.mass", Truth() + "m");
      setBranchAlias("MCParticles.PDG", Truth() + "pid");
      setBranchAlias("MCParticles.generatorStatus", Truth() + "genStat");

      // Map full MC array to final state-only array.
      Define(Truth() + "final_id", [](const ROOT::RVecI& stat) {
        auto indices = helpers::Enumerate<uint>(stat.size());
        return indices[stat == 1];
      }, {Truth() + "genStat"});

      reaction::util::CountParticles(this, Truth());

      if (IsEnd) {
        reaction::util::RedefineFundamentalAliases(this);
        ApplyAfterBurner(Truth());
        ApplyAfterBurnerOnBeams();
        DefineBeamElectron();
        DefineBeamIon();
      }
    }

    inline void ePICReaction::AliasColumnsAndMatchWithMC(Bool_t IsEnd) {
      AliasColumnsAndMC(kFALSE);

      if (IsEnd) {
        reaction::util::RedefineFundamentalAliases(this);
        RedefineViaAlias(Rec() + "m", [](const ROOT::RVecF& recm, const ROOT::RVecD& trum) { return helpers::Truncate(ROOT::RVecF(trum), recm.size()); }, {Rec() + "m", Truth() + "m"});
        ApplyAfterBurner(Truth());
        ApplyAfterBurner(Rec());
        ApplyAfterBurnerOnBeams();
        DefineBeamElectron();
        DefineBeamIon();
        _truthMatched = true;
      }
    }

    inline void ePICReaction::PostParticles() {
      AddAdditionalComponents();

      if (IsTruthMatched()) {
        // Ensure tru and rec columns are defined together and add resolution columns.
        Filter(Form("bool((%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty()+%s.empty())*0 +1); ",
                    (Truth() + "pmag").data(),
                    (Truth() + "theta").data(),
                    (Truth() + "phi").data(),
                    (Truth() + "eta").data(),
                    (Rec() + "pmag").data(),
                    (Rec() + "theta").data(),
                    (Rec() + "phi").data(),
                    (Rec() + "eta").data()),"truthmatch");
	
        reaction::util::ResolutionFraction(this, "pmag");
        reaction::util::Resolution(this, "theta");
        reaction::util::Resolution(this, "phi");
        reaction::util::Resolution(this, "eta");
      }
    }

    inline void ePICReaction::AddAdditionalComponents() {
      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
      DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
    }

    template<typename T>
    inline void ePICReaction::RedefineFundamental(const std::string& name) {
      auto contains = [](const std::string& s1, const std::string& s2) {
        return (s1.find(s2) != std::string::npos);
      };
      if (contains(name, "rec")) {
        RedefineViaAlias(name, helpers::Reorder<T, UInt_t, UInt_t>, {name.data(), Rec() + "match_id", Truth() + "match_id", Truth() + "n"});
      } else if (contains(name, "tru")) {
        RedefineViaAlias(name, helpers::Rearrange<T, UInt_t>, {name.data(), Truth() + "final_id"});
      }
    }

    inline void ePICReaction::SetBeamsFromMC(Long64_t nrows) {
      _useBeamsFromMC = true;
      auto nthreads = ROOT::GetThreadPoolSize();
      if (nthreads) ROOT::DisableImplicitMT();

      auto tempframe = GetFileNames().size() == 0 ?
        ROOT::RDataFrame{GetTreeName(), GetFileName()} :
        ROOT::RDataFrame{GetTreeName(), GetFileNames()};

      auto beamdf = tempframe.Range(nrows)
        .Define("emean", "MCParticles.momentum.z[MCBeamElectrons_objIdx.index[0]]")
        .Define("pzmean", "MCParticles.momentum.z[MCBeamProtons_objIdx.index[0]]")
        .Define("pxmean", "MCParticles.momentum.x[MCBeamProtons_objIdx.index[0]]");
      auto pze = beamdf.Mean("emean");
      auto pzp = beamdf.Mean("pzmean");
      auto pxp = beamdf.Mean("pxmean");

      setBeamElectron(0, 0, *pze);
      setBeamIon(*pxp, 0, *pzp);

      if (nthreads) ROOT::EnableImplicitMT(nthreads);
    }
    inline void ePICReaction::SetBeamsFromMC(UInt_t iel,UInt_t iion,Long64_t nrows) {
      _useBeamsFromMC = true;
      auto nthreads = ROOT::GetThreadPoolSize();
      if (nthreads) ROOT::DisableImplicitMT();

      auto tempframe = GetFileNames().size() == 0 ?
        ROOT::RDataFrame{GetTreeName(), GetFileName()} :
        ROOT::RDataFrame{GetTreeName(), GetFileNames()};

      auto beamdf = tempframe.Range(nrows)
        .Define("emean", Form("MCParticles.momentum.z[%d]",iel))
        .Define("pzmean", Form("MCParticles.momentum.z[%d]",iion))
        .Define("pxmean", Form("MCParticles.momentum.x[%d]",iion));
      auto pze = beamdf.Mean("emean");
      auto pzp = beamdf.Mean("pzmean");
      auto pxp = beamdf.Mean("pxmean");

      setBeamElectron(0, 0, *pze);
      setBeamIon(*pxp, 0, *pzp);

      if (nthreads) ROOT::EnableImplicitMT(nthreads);
    }

    inline void ePICReaction::ApplyAfterBurner(std::string type) {
      if (DeduceColumnVectorType(this, type + "px") == ColType::Float &&
	  DeduceColumnVectorType(this, type + "m") == ColType::Double) {
        RedefineViaAlias(type + "px",
			 UndoAfterBurn<float, double>{_p4ion_beam, _p4el_beam, -0.025},
			 {type + "px", type + "py", type + "pz", type + "m"});
      } else if (DeduceColumnVectorType(this, type + "px") == ColType::Double &&
		 DeduceColumnVectorType(this, type + "m") == ColType::Double) {
        RedefineViaAlias(type + "px",
			 UndoAfterBurn<double, double>{_p4ion_beam, _p4el_beam, -0.025},
			 {type + "px", type + "py", type + "pz", type + "m"});
      } else if (DeduceColumnVectorType(this, type + "px") == ColType::Float &&
		 DeduceColumnVectorType(this, type + "m") == ColType::Float) {
        RedefineViaAlias(type + "px",
			 UndoAfterBurn<float, float>{_p4ion_beam, _p4el_beam, -0.025},
			 {type + "px", type + "py", type + "pz", type + "m"});
      } else if (DeduceColumnVectorType(this, type + "px") == ColType::Double &&
		 DeduceColumnVectorType(this, type + "m") == ColType::Float) {
        RedefineViaAlias(type + "px",
			 UndoAfterBurn<double, float>{_p4ion_beam, _p4el_beam, -0.025},
			 {type + "px", type + "py", type + "pz", type + "m"});
      } else {
        std::cerr << "Error ePICReaction::ApplyAfterBurner momentum and mass types not valid "
                  << CurrFrame().GetColumnType(type + "px") << " " << CurrFrame().GetColumnType(type + "m") << std::endl;
        exit(0);
      }
    }

    inline void ePICReaction::ApplyAfterBurnerOnBeams() {
      auto beams_px = ROOT::RVecD{_p4el_beam.X(), _p4ion_beam.X()};
      auto beams_py = ROOT::RVecD{_p4el_beam.Y(), _p4ion_beam.Y()};
      auto beams_pz = ROOT::RVecD{_p4el_beam.Z(), _p4ion_beam.Z()};
      auto beams_m = ROOT::RVecD{_p4el_beam.M(), _p4ion_beam.M()};
      std::cout << "Pre Undo afterburn head on beam 4-vectors : " << _p4el_beam << " " << _p4ion_beam << std::endl;
      rad::epic::UndoAfterBurn<double, double> undoAB_DD{_p4ion_beam, _p4el_beam, -0.025};
      undoAB_DD(beams_px, beams_py, beams_pz, beams_m);
      _p4el_beam.SetCoordinates(beams_px[0], beams_py[0], beams_pz[0], beams_m[0]);
      _p4ion_beam.SetCoordinates(beams_px[1], beams_py[1], beams_pz[1], beams_m[1]);
      std::cout << "Undo afterburn head on beam 4-vectors : " << _p4el_beam << " " << _p4ion_beam << std::endl;
    }

    inline bool ePICReaction::IsTruthMatched() const { return _truthMatched; }

  } // namespace config
} // namespace rad
