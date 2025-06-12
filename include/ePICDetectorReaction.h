#pragma once

//!  Derived class to configure ePIC root files with associations made to detector and track info

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order. It also uses podio to link particles to their full detector reconstruction
*/
#include "ePICReaction.h"
#include "PodioMetadata.h"
#include "PodioRelated.h"
#include "RVecHelpers.h"
#include "StringUtilities.h"
#include <TFile.h>
#include <TTree.h>


namespace rad{
  namespace config{
    using rad::names::data_type::Rec;
    using rad::names::data_type::Truth;
    using rad::utils::createFunctionCallString;
    using rad::utils::replaceAll;

    //! Class definition

    class ePICDetectorReaction : public ePICReaction {

      
    public:

      ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} );
      ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} );

      
      /**
       * Configure association of clusters data to particles
       */
      void AssociateClusters(const std::vector<std::string>& types,const std::vector<std::string>& members ){
	AssociateObjects("clusters",types,members);
      }
     /**
       * Configure association of tracks data to particles
       */
      void AssociateTracks(const std::vector<std::string>& types,const std::vector<std::string>& members ){
	AssociateObjects("tracks",types,members);
      }

      /**
       * Configure given associations to create objects which may 
       * be used to create columns matched to reconstructed particles
       */
      void AssociateObjects(const string& object,const std::vector<std::string>& types,const std::vector<std::string>& members);


      
    private:
      rad::podio::PodioMetadata _podio_meta;
      
    };
    ///////////////////////////////////////////////
    /// Function definitions
    ///////////////////////////////////////////////
    
    ePICDetectorReaction::ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ) : ePICReaction{treeName,fileNameGlob,columns} {
	//open file and copy podio table
	TChain chain("podio_metadata");
	chain.Add(fileNameGlob.data());

	_podio_meta = rad::podio::PodioMetadata(chain.GetListOfFiles()->At(0)->GetTitle());
      }
    ePICDetectorReaction::ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ) : ePICReaction{treeName,filenames,columns} {
      //open file and copy podio table
      _podio_meta = rad::podio::PodioMetadata(filenames[0].data());
      
    }

    
    void ePICDetectorReaction::AssociateObjects(const std::string& object, const std::vector<std::string>& types, const std::vector<std::string>& members) {
      if (types.empty()) return;

      ROOT::RVecU collIndices;
      collIndices.reserve(types.size());

      // Collect valid association names and their collection IDs
      std::vector<std::string> valid_types;
      for (const auto& assoc_name : types) {
        if (!_podio_meta.Exists(assoc_name)) {
	  //#ifdef RAD_DEBUG
	  std::cerr << "Warning : ePICDetectorReaction::AssociateObjects, no detector object "
		    << assoc_name << " in podio_metadata." << std::endl;
	  //#endif
	  continue;
        }
        collIndices.push_back(_podio_meta.CollectionIDFor(assoc_name));
        valid_types.push_back(assoc_name);
	// #ifdef RAD_DEBUG
        //std::cout << "AssociateObjects " << object << " " << assoc_name << " " << collIndices.back() << std::endl;
        //#endif
      }
      if (valid_types.empty()) return;

      // Create functor to convert event collectionIDs to local indices
      auto local_collIDs = rad::podio::ConvertCollectionId(collIndices);
      auto coll_IdxsName = object + "_idxs";
      Define(coll_IdxsName, local_collIDs, {"_ReconstructedParticles_" + object + ".collectionID"});
    
      // Loop over data members and Define a vector for each association
      for (const auto& member : members) {
        std::ostringstream memberNames;
        for (const auto& assoc_name : valid_types) {
	  memberNames << assoc_name << '.' << member << ',';
        }
        std::string memberNamesStr = memberNames.str();
        if (!memberNamesStr.empty()) memberNamesStr.pop_back(); // remove last comma

        auto collListName = object + member + DoNotWriteTag();
        replaceAll(collListName, ".", "_");
        auto member_type = CurrFrame().GetColumnType(valid_types[0] + '.' + member);
        Define(collListName, Form("ROOT::RVec<%s>{%s}", member_type.data(), memberNamesStr.data()));

        // Transformation for object collection
        auto rec_index = "_ReconstructedParticles_" + object + ".index";
        std::string rec_begin = Form("ReconstructedParticles.%s_begin", object.data());
        std::string rec_end = Form("ReconstructedParticles.%s_end", object.data());
        auto func_rec_object = createFunctionCallString("rad::podio::combineOneToMany",
                                                        collListName, coll_IdxsName, rec_index,
                                                        rec_begin, rec_end);
        std::string rec_object = Rec() + object + "_" + member;
        replaceAll(rec_object, ".", "_");

        Define(rec_object, func_rec_object);

        // Reorder to match truth and rec if required
        if (rad::config::ColumnExists(Truth() + "match_id", CurrFrame())) {
	  auto func_rec_match = createFunctionCallString("rad::helpers::Reorder",
							 rec_object, Rec() + "match_id", Truth() + "match_id", Truth() + "n");
	  RedefineExpr(rec_object, func_rec_match);
        }
      }
    
    }



    
  }
}
