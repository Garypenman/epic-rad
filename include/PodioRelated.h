#pragma once

#pragma link C++ class ROOT::VecOps::RVec<float>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>+;

/**
 *associations tell us how additional information are synchronised with a base collection
 * Base+object+.collectionID  e.g. _ReconstructedParticles_clusters.collectionID
 * Base+object+.index         e.g. _ReconstructedParticles_clusters.index
 */
#include "Constants.h"
#include <ROOT/RVec.hxx>
#include <map>        // For std::map (if considering sorted unique IDs)
#include <unordered_map> // For std::unordered_map (for efficient ID-to-index mapping)

namespace rad{
  namespace podio{

    /**
     * @brief An optimized class to convert a vector of raw collection IDs to a
     * vector of local indices based on a pre-defined mapping.
     *
     * This class is designed to be a callable object (a functor) that can be
     * used, for example, in a data processing pipeline (like ROOT's RDataFrame
     * transformations) to map large sets of collection IDs to their corresponding
     * integer indices, handling cases where an ID is not found.
     *
     * This version uses an `std::unordered_map` internally for highly efficient
     * ID-to-index lookups (average O(1) complexity).
     */
    class ConvertCollectionId {
    public:
      /**
       * @brief Constructor for convertCollectionId.
       * It builds an internal hash map to store the mapping from collection IDs
       * to their local indices. This map is built once during construction
       * for fast subsequent lookups.
       *
       * @param collIDs A ROOT::RVecU containing the original (unique) collection IDs
       * that define the mapping. The index of each ID in this vector
       * will be its corresponding local index.
       */
      ConvertCollectionId(const ROOT::RVecU& collIDs) {
        // Clear any previous map content (important if object is reused or default constructed)
        _idToIndexMap.clear();
        // Pre-allocate map memory based on the size of collIDs for better performance
        // and to reduce rehashes.
        _idToIndexMap.reserve(collIDs.size());

        // Populate the hash map: iterate through the input IDs and store
        // each ID as a key, with its vector index as the value.
        // This is an O(M) operation on average, where M is collIDs.size().
        for (size_t i = 0; i < collIDs.size(); ++i) {
	  _idToIndexMap[collIDs[i]] = static_cast<int>(i);
        }
      }

      /**
       * @brief Overloaded function call operator to perform the ID-to-index conversion.
       * This makes the class object callable like a function.
       *
       * This method uses the pre-built `_idToIndexMap` for very fast lookups.
       *
       * @param collID A ROOT::RVecU containing the collection IDs for which local
       * indices need to be found. These are the IDs from the actual
       * data being processed.
       * @return A ROOT::RVecI where each element is the local index corresponding
       * to the input collection ID, or -1 if the ID is not found in the
       * pre-defined mapping (`_idToIndexMap`).
       */
      ROOT::RVecI operator()(const ROOT::RVecU& collID) {
        // Create a new RVecI to store the resulting local indices.
        // It's sized to match the input vector of IDs for efficient population.
        ROOT::RVecI localID(collID.size(),rad::constant::InvalidEntry<int>());
        uint i = 0; // Initialize a counter for the output vector index.

        // Iterate through each collection ID in the input vector.
        // Each lookup in the unordered_map is, on average, an O(1) operation.
        // Thus, the entire loop is O(N) on average, where N is collID.size().
        for (const auto id : collID) {
	  // Attempt to find the current 'id' in the hash map.
	  auto it = _idToIndexMap.find(id);

	  // Check if the iterator returned by find() points to the end of the map.
	  // If it does, the ID was not found in our pre-defined mapping.
	  if (it != _idToIndexMap.end()) {
	    // If the ID is found, 'it->second' contains its corresponding local index.
	    localID[i] = it->second;
	  }//  else {
	  //   // If not found, assign -1 to indicate an unknown ID.
	  //   localID[i] = -1;
	  // }
	  ++i; // Move to the next position in the output vector.
        }
        return localID; // Return the vector of local indices.
      }

    private:
      // Member variable: An unordered_map to store the mapping from original
      // collection ID (unsigned int) to its corresponding local index (int).
      // This provides fast average O(1) lookup times.
      std::unordered_map<unsigned int, int> _idToIndexMap;
      // The original _origIds vector is no longer needed as the map stores the mapping.
    };

  
    
    /**
     * @brief Take all given collections and combine
     * The combined vector is ordered by the given association which defines inputs local_collIds and order
     * @param collections list of all values in all given collections : vector of collection's with vectors of values
     * @param local_collIds codes for collections in local scope, converted from podio collectionIDs
     * @param order of values from original collection in result
     * @return vector of values ordered as per association
     */
    template <typename Tval, typename Tind1, typename Tind2>
      ROOT::RVec<Tval> combineCollections(const ROOT::RVec<ROOT::RVec<Tval>>& collections,const ROOT::RVec<Tind1>& local_collIds,const ROOT::RVec<Tind2>& order){

      auto Nelements = local_collIds.size(); //number of collections to combine

      //note idxs.size() - number of objects(e.g.clusters) with rec_particle association
      ROOT::RVec<Tval> result(Nelements,rad::constant::InvalidEntry<Tval>());
      //loop over ReconstructedParticles and get the associated data
      for(uint i=0;i<Nelements;++i){
	result[i] = (rad::constant::IsInvalidEntry(local_collIds[i]) ) ?
	  rad::constant::InvalidEntry<Tval>() : collections[local_collIds[i]][order[i]];
      }
      // std::cout<<"CombineCollections "<<collections<<" "<<local_collIds<<" "<<order<<" "<<result<<std::endl;
      return result;
    };

    /**
     * @brief take in a base collection with one-to-many relation;
     * its begin and end points; and a vector of values to arrange
     * loop through its limits and associate the corresponing values
     * from object_vals
     * @param rec_begin indices of first entry on object_vals
     * @param rec_begin indices of last entry on object_vals
     * @param object_vals vector of values this relation points to.
     *        This may be created by CombineCollections   
     * @return for each entry in our base collection a vector of these object values
     */
    template <typename Tval, typename Tlims>
      ROOT::RVec<ROOT::RVec<Tval>> vecOfOneToMany(const ROOT::RVec<Tlims>& rec_begin,const ROOT::RVec<Tlims>& rec_end,const ROOT::RVec<Tval>& object_vals){
      
      auto n_rec = rec_begin.size(); //loop size
      ROOT::RVec<ROOT::RVec<Tval>> vals; //results vector
      for(size_t i = 0;i < n_rec; ++i){
	//for each entry add a vector of values
	//with range given by begin and end
	auto start_iterator = object_vals.begin() + rec_begin[i];
	auto end_iterator = object_vals.begin() + rec_end[i];
	vals.emplace_back(start_iterator,end_iterator);
      }
      // cout<< "VectorValsByLimits "<<vals<< endl;
      return vals;
    }

    /**
     * @brief Combine the 2 steps for making the associated data structure
     */
    template <typename Tval, typename Tind1, typename Tind2, typename Tlims>
    ROOT::RVec<ROOT::RVec<Tval>> combineOneToMany(const ROOT::RVec<ROOT::RVec<Tval>>& collections,const ROOT::RVec<Tind1>& local_collIds,const ROOT::RVec<Tind2>& order,const ROOT::RVec<Tlims>& rec_begin,const ROOT::RVec<Tlims>& rec_end){

      auto object_vals = combineCollections(collections,local_collIds,order);
      return vecOfOneToMany(rec_begin,rec_end,object_vals);
      
    }
    
  }//podio
}//rad
