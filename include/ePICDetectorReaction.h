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

      ePICDetectorReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ePICReaction{treeName,fileNameGlob,columns} {
	//open file and copy podio table
	TChain chain("podio_metadata");
	chain.Add(fileNameGlob.data());

	_podio_meta = rad::podio::PodioMetadata(chain.GetListOfFiles()->At(0)->GetTitle());
      }
    ePICDetectorReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ePICReaction{treeName,filenames,columns} {
      //open file and copy podio table
      _podio_meta = rad::podio::PodioMetadata(filenames[0].data());

      }

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
      void AssociateObjects(const string& object,const std::vector<std::string>& types,const std::vector<std::string>& members){
	if(types.empty()==true) return;
	ROOT::RVecU collIndices; //create indices for the given collection types
      
	for(const auto& assoc_name:types){ //loop over the specified detector associations
	  if(_podio_meta.Exists(assoc_name)==false){
	    std::cerr<<"Warning : ePICDetectorReaction::AssociateObjects, no detector object "
		     <<assoc_name<<" in podio_metadata. "<<std::endl;
	    continue;
	  }
	  
	  //save the collectionID value in indices
	  //this allows us to define a local index for the collection
	  collIndices.push_back(_podio_meta.CollectionIDFor(assoc_name));
	  std::cout<<"AssociateObjects "<<object<<" "<<assoc_name<<" "<<collIndices.back()<<std::endl;
	}
	
	//create functor to convert event collectionIDs to local indices
	auto local_collIDs = rad::podio::ConvertCollectionId(collIndices);
	//define collection indices for this object association
	auto coll_IdxsName = object+"_idxs";
	Define(coll_IdxsName , local_collIDs, {"_ReconstructedParticles_"+object+".collectionID"});
	
	//loop over given data members and Define a vector for each association 
	for(const auto& member:members){//loop over the data member we are interested in
	  
	  std::string memberNames ;
	  for(const auto& assoc_name:types){ //loop over the specified detector associations
	    
	    memberNames+=assoc_name; //association name  e.g. CentralCKFTracks
	    memberNames+="."+member; //plus specific data member required e.g. momentum.x
	    memberNames+=","; //seperate each member by a comma e.g. CentralCKFTracks.momentum.x,
	    
	  }
	  memberNames.pop_back();//remove last ,
	  
	  //Define a list of corresponding collection names
	  //convert to indices of local collections
	  auto collListName = object+member+DoNotWriteTag();
	  replaceAll(collListName,".","_");
	  auto member_type = CurrFrame().GetColumnType(types[0]+'.'+member);
	  Define(collListName,Form("ROOT::RVec<%s>{%s}",member_type.data(),memberNames.data()));

	  //Here we define the transformation from the object collection (e.g. EcalEndcapNClusters.energy to a vector ordered in ReconstructedParticles
	  //memberName e.g. {{EcalEndcapPClusters.energy,EcalEndcapNClusters.energy,...}, {*Tracks*.energy,...} }
	  //coll_IdxsName e.g. {1}, i.e. use  EcalEndcapNClusters.energy //local version of _ReconstructedParticles_clusters.collectionID
	  // and put this value at "_ReconstructedParticles_"+object+".index" -> "_ReconstructedParticles_clusters.index

	  //create strings for the Define
	  auto rec_index = string("_ReconstructedParticles_")+object+".index"; //object indices
	  std::string rec_begin = Form("ReconstructedParticles.%s_begin",object.data());
	  std::string rec_end = Form("ReconstructedParticles.%s_end",object.data());
	  auto func_rec_object =createFunctionCallString("rad::podio::combineOneToMany" ,
							 collListName, coll_IdxsName, rec_index,
							 rec_begin ,rec_end); 
	  std::string rec_object = Rec()+object+"_"+member;
	  replaceAll(rec_object,".","_");
	  
	  Define(rec_object , func_rec_object);

	  //reorder to match truth and rec if required
	  if(rad::config::ColumnExists(Truth()+"match_id",CurrFrame()) == true ){
	    auto func_rec_match =createFunctionCallString("rad::helpers::Reorder" ,
							  rec_object, Rec()+"match_id",Truth()+"match_id",Truth()+"n"); 
	    RedefineExpr(rec_object, func_rec_match );
	  }
	}
	
      }


      
    private:
      rad::podio::PodioMetadata _podio_meta;
      
      //  podio::CollectionIDTable _idTable; //to be cloned from podio_metadata:events___idTable
      // podio::CollectionIDTable *_idTable=nullptr; //to be cloned from podio_metadata:events___idTable, it has a deleted constructor so cannot use copy contructor
      //std::shared_ptr<podio::CollectionIDTable> _idTable; //to be cloned from podio_metadata:events___idTable
      /**
       * Local index for particle associated collection IDs
       * should be in range 0 - N declared trackers
       */
      /* ROOT::RVecU _clustersCollectionIdxs; */
      /* ROOT::RVecU _tracksCollectionIdxs; */
      /* ROOT::RVecU _particleIDsCollectionIdxs; */
      /* ROOT::RVecU _startVertexCollectionIdxs; */
      /* ROOT::RVecU _particleIDUsedCollectionIdxs; */
      
    };
  }
}
