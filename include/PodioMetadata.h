#pragma once

#include "StringUtilities.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

namespace podio{
  namespace root_utils{
    // Struct definition from podio::root_utils
    struct CollectionWriteInfo {
      uint32_t collectionID{static_cast<uint32_t>(-1)};
      std::string dataType{};
      bool isSubset{false};
      unsigned int schemaVersion{0};
      std::string name{};
      std::string storageType{};
    };
    
  }
}

#pragma link C++ class podio::root_utils::CollectionWriteInfo+;
#pragma link C++ class vector<podio::root_utils::CollectionWriteInfo>+;

namespace rad{
  namespace podio{
    using RVecS = ROOT::RVec<std::string>;
    using RVecU  = ROOT::RVecU;
    
    class PodioMetadata {

    public:
      PodioMetadata() = default;
      PodioMetadata(const std::string& filename, const std::string& treename="podio_metadata"){

	std::cout<<"PodioMetadata file : "  <<filename<<" "<<treename<<std::endl;
	/* ROOT::RDataFrame df(treename,filename); */
	/* //Get the collectionId and names column data */
	/* auto ids_ptr = df.Take<RVecU>("m_collectionIDs"); */
	/* auto names_ptr = df.Take<RVecS>("m_names"); */

	/* //Take the vectors of the first event */
	/* //and keep them */
	/* _collectionIDs = (*(ids_ptr ))[0]; */
	/* _names = (*(names_ptr ))[0]; */
	
	/* for(uint i=0;i<_names.size();++i){ */
	/*   cout<<_names[i]<<"\t collection id = "<<_collectionIDs[i]<<endl; */
	/* } */
	
	// g. penman 17.11.25
	//new podio version fix
        
	// Open file
        auto file = std::unique_ptr<TFile>(TFile::Open(filename.c_str()));
        if (!file || file->IsZombie()) {
          throw std::runtime_error("Failed to open file: " + filename);
        }

        // Get metadata tree
        TTree *tree = file->Get<TTree>(treename.c_str());
        if (!tree) {
          throw std::runtime_error("Metadata tree not found: " + treename);
        }

        // Read CollectionWriteInfo vector
        std::vector<::podio::root_utils::CollectionWriteInfo> *info = nullptr;
        tree->SetBranchAddress("events___CollectionTypeInfo", &info);
        tree->GetEntry(0);

        if (!info || info->empty()) {
          throw std::runtime_error("No CollectionTypeInfo found in metadata.");
        }

        // Fill internal vectors
        for (const auto &entry : *info) {
          _collectionIDs.push_back(entry.collectionID);
          _names.push_back(entry.name);
        }
      }

      bool Exists(const std::string name){
	if(std::find(_names.begin(),_names.end(),name)==_names.end()) return false;
	return true;
      }

      UInt_t CollectionIDFor(const std::string& name){
	//no checks made, use Exists first
	auto it = std::find(_names.begin(),_names.end(),name);
	UInt_t index =  it - _names.begin();
	return _collectionIDs[index];
      }

      RVecS FilterNames(std::string sub_string){
	return rad::utils::filterStrings(_names,sub_string);
      }
    private:
      
      RVecU _collectionIDs;
      RVecS _names;
      
    };

  }
}
