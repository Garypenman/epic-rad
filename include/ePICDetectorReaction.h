// ePICDetectorReaction.h
#pragma once

//! Derived class for configuring ePIC root files with associations to detector and track info

/*!
  ePICDetectorReaction helps set up data analysis and calculations
  for specific hadronic final states using ePIC files.
  This class configures files with a fixed particle order, and
  uses podio to link particles to their full detector reconstruction.
*/

#include "ePICReaction.h"
#include "PodioMetadata.h"
#include "PodioRelated.h"
#include "RVecHelpers.h"
#include <TFile.h>
#include <TTree.h>

namespace rad{
  namespace config{
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;
    using rad::utils::createFunctionCallString;
    using rad::utils::replaceAll;
    using RVecS = ROOT::RVec<std::string>;

    //! ePICDetectorReaction: Handles association of detector objects and tracks to particles for ePIC analyses
    class ePICDetectorReaction : public ePICReaction {

    public:

      /**
       * Constructor: Open files matching a glob pattern, set up tree and columns.
       * Also loads associated podio metadata for collection IDs.
       */
      ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} );

      /**
       * Constructor: Accepts vector of file names. Sets up tree, columns, and podio metadata.
       */
      ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} );

      /**
       * Associates clusters data collections to corresponding particles.
       * @param types   List of cluster collection names.
       * @param members Data fields within each cluster collection to associate.
       */
      void AssociateClusters(const RVecS& types,const RVecS& members ){
        AssociateObjects("clusters",types,members);
      }

      /**
       * Associates tracks data collections to corresponding particles.
       * @param types   List of track collection names.
       * @param members Data fields within each track collection to associate.
       */
      void AssociateTracks(const RVecS& types,const RVecS& members ){
        AssociateObjects("tracks",types,members);
      }

      /**
       * General association function for any detector object (e.g. clusters, tracks).
       * Creates columns matched to reconstructed particles.
       * @param object  Type of object (e.g. "clusters", "tracks").
       * @param types   List of collection names for this object type.
       * @param members Data field names to extract from each collection.
       */
      void AssociateObjects(const std::string& object,const RVecS& types,const RVecS& members);


      /**
       * Getter for podio metadata
       */
      rad::podio::PodioMetadata &PodioInfo(){return _podio_meta;}
    private:
      // Holds podio metadata for the currently loaded file(s)
      rad::podio::PodioMetadata _podio_meta;
    };

    ///////////////////////////////////////////////
    /// Inline function implementations
    ///////////////////////////////////////////////

    inline ePICDetectorReaction::ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t& columns )
      : ePICReaction{treeName,fileNameGlob,columns} // initialize base class
    {
      // Open podio metadata table using TChain and copy metadata for the first file
      TChain chain("podio_metadata");
      chain.Add(fileNameGlob.data());

      _podio_meta = rad::podio::PodioMetadata(chain.GetListOfFiles()->At(0)->GetTitle());
    }

    inline ePICDetectorReaction::ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t& columns )
      : ePICReaction{treeName,filenames,columns} // initialize base class
    {
      // Copy podio metadata from the first file in the vector
      _podio_meta = rad::podio::PodioMetadata(filenames[0].data());
    }

    inline void ePICDetectorReaction::AssociateObjects(const std::string& object, const RVecS& types, const RVecS& members) {
      // Exit if no types provided
      if (types.empty()) return;

      ROOT::RVecU collIndices;
      collIndices.reserve(types.size());

      // Collect valid collection names and their IDs, warn if not found in metadata
      std::vector<std::string> valid_types;
      for (const auto& assoc_name : types) {
        if (!_podio_meta.Exists(assoc_name)) {
          std::cerr << "Warning: ePICDetectorReaction::AssociateObjects - No detector object "
                    << assoc_name << " in podio_metadata." << std::endl;
          continue;
        }
        collIndices.push_back(_podio_meta.CollectionIDFor(assoc_name));
        valid_types.push_back(assoc_name);
      }
      if (valid_types.empty()) return;

      // Set up function to map event collection IDs to local indices for this object type
      auto local_collIDs = rad::podio::ConvertCollectionId(collIndices);
      auto coll_IdxsName = object + "_idxs";
      Define(coll_IdxsName, local_collIDs, {"_ReconstructedParticles_" + object + ".collectionID"});

      // For each data member, define a vector for each association and create the corresponding columns
      for (const auto& member : members) {
        std::ostringstream memberNames;
        for (const auto& assoc_name : valid_types) {
          memberNames << assoc_name << '.' << member << ',';
        }
        std::string memberNamesStr = memberNames.str();
        if (!memberNamesStr.empty()) memberNamesStr.pop_back(); // remove last comma

        // Create column name (replace '.' with '_')
        auto collListName = object + member + DoNotWriteTag();
        replaceAll(collListName, ".", "_");

        // Get type of this member from the first valid collection
        auto member_type = CurrFrame().GetColumnType(valid_types[0] + '.' + member);
        Define(collListName, Form("ROOT::RVec<%s>{%s}", member_type.data(), memberNamesStr.data()));

        // Set up transformation for the object collection
        auto rec_index = "_ReconstructedParticles_" + object + ".index";
        std::string rec_begin = Form("ReconstructedParticles.%s_begin", object.data());
        std::string rec_end = Form("ReconstructedParticles.%s_end", object.data());
        auto func_rec_object = createFunctionCallString("rad::podio::combineOneToMany",
                                                        collListName, coll_IdxsName, rec_index,
                                                        rec_begin, rec_end);
	
        std::string rec_object = Rec() + object + "_" + member;
        replaceAll(rec_object, ".", "_");

        Define(rec_object, func_rec_object);

        // If truth/reco matching exists, reorder to align with truth particles
        if (rad::config::ColumnExists(Truth() + "match_id", CurrFrame())) {
          auto func_rec_match = createFunctionCallString("rad::helpers::Reorder",
                                                         rec_object, Rec() + "match_id", Truth() + "match_id", Truth() + "n");
          RedefineExpr(rec_object, func_rec_match);
        }
      }
    }

  } // namespace config
} // namespace rad
