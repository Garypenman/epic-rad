#pragma once

#pragma link C++ class ROOT::VecOps::RVec<float>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>+;

/**
 * Associations describe how additional information is synchronized with a base collection.
 * For example:
 *   - Base+object+.collectionID  (e.g., _ReconstructedParticles_clusters.collectionID)
 *   - Base+object+.index         (e.g., _ReconstructedParticles_clusters.index)
 */
#include "Constants.h"
#include <ROOT/RVec.hxx>
#include <map>           // For std::map (sorted, unique IDs)
#include <unordered_map> // For std::unordered_map (fast ID-to-index mapping)

namespace rad{
  namespace podio{

    /**
     * @brief Functor to map a vector of collection IDs to their corresponding local indices.
     *
     * This class builds an efficient lookup (hash map) from collection IDs to indices,
     * enabling fast translation of arbitrary ID lists to index lists, typically for use
     * in data processing workflows (e.g., ROOT RDataFrame).
     */
    class ConvertCollectionId {
    public:
      /**
       * @brief Constructor. Builds the lookup map from collection IDs to their local indices.
       *
       * @param collIDs Vector of unique collection IDs. Each ID's position is its local index.
       */
      ConvertCollectionId(const ROOT::RVecU& collIDs) {
        _idToIndexMap.clear();
        _idToIndexMap.reserve(collIDs.size());
        for (size_t i = 0; i < collIDs.size(); ++i) {
          _idToIndexMap[collIDs[i]] = static_cast<int>(i);
        }
      }

      /**
       * @brief Map input collection IDs to their local indices using the pre-built lookup.
       *
       * @param collID Vector of (possibly non-unique) collection IDs to map.
       * @return Vector of local indices, or rad::constant::InvalidEntry<int>() if an ID is not found.
       */
      ROOT::RVecI operator()(const ROOT::RVecU& collID) const {
	ROOT::RVecI localID(collID.size(), rad::constant::InvalidEntry<int>());
	for (size_t i = 0; i < collID.size(); ++i) {
	  auto it = _idToIndexMap.find(collID[i]);
	  if (it != _idToIndexMap.end()) {
	    localID[i] = it->second;
	  }
	}
	return localID;
      }

    private:
      // Maps original collection IDs (unsigned int) to local indices (int) for fast lookup.
      std::unordered_map<unsigned int, int> _idToIndexMap;
    };

    /**
     * @brief Merge multiple collections into a single vector, ordering by association.
     *
     * @param collections  Vector of collections (each a vector of values).
     * @param local_collIds  Vector of local collection IDs (converted from global IDs).
     * @param order  Vector specifying which item to take from each collection.
     * @return Vector of combined values, each selected according to the association.
     */
    template <typename Tval, typename Tind1, typename Tind2>
      ROOT::RVec<Tval> combineCollections(const ROOT::RVec<ROOT::RVec<Tval>>& collections,
                                          const ROOT::RVec<Tind1>& local_collIds,
                                          const ROOT::RVec<Tind2>& order)
    {
      // std::cout<<"combineCollections "<<collections<<" "<<local_collIds<<" "<<order<<std::endl;
      auto Nelements = local_collIds.size();
      ROOT::RVec<Tval> result(Nelements, rad::constant::InvalidEntry<Tval>());
      for (uint i = 0; i < Nelements; ++i) {
        result[i] = (rad::constant::IsInvalidEntry(local_collIds[i]))
                    ? rad::constant::InvalidEntry<Tval>()
                    : collections[local_collIds[i]][order[i]];
      }
      return result;
    };

    /**
     * @brief For a one-to-many relationship, group values for each base entry.
     *
     * Given begin and end indices for each base object, slice out the associated values
     * from a flat vector and group into sub-vectors.
     *
     * @param rec_begin  Start indices for each group in object_vals.
     * @param rec_end    End indices (exclusive) for each group in object_vals.
     * @param object_vals Flat vector of all associated values (to be grouped).
     * @return Vector of vectors: one vector of values per base entry.
     */
    template <typename Tval, typename Tlims>
      ROOT::RVec<ROOT::RVec<Tval>> vecOfOneToMany(const ROOT::RVec<Tlims>& rec_begin,
                                                  const ROOT::RVec<Tlims>& rec_end,
                                                  const ROOT::RVec<Tval>& object_vals)
    {
      auto n_rec = rec_begin.size();
      //std::cout<< "vecOfOneToMany " << rec_begin<<" "<<rec_end<<" "<<object_vals<<std::endl;
      //vals vector has 1 entry for every element in rec_begin (e.g. for every reconstructed particle)
      ROOT::RVec<ROOT::RVec<Tval>> vals;
      vals.reserve(n_rec);
      for (size_t i = 0; i < n_rec; ++i) {
        auto start_iterator = object_vals.begin() + rec_begin[i];
        auto end_iterator = object_vals.begin() + rec_end[i];
        vals.emplace_back(start_iterator, end_iterator);
      }
      return vals;
    }

    /**
     * @brief Combine collection merging and one-to-many grouping into a single step.
     *
     * @param collections     Vector of collections (each a vector of values).
     * @param local_collIds   Vector of local collection IDs.
     * @param order           Vector specifying which item to take from each collection.
     * @param rec_begin       Start indices for each one-to-many group.
     * @param rec_end         End indices for each one-to-many group.
     * @return Nested vector: grouped values for each base entry.
     */
    template <typename Tval, typename Tind1, typename Tind2, typename Tlims>
    ROOT::RVec<ROOT::RVec<Tval>> combineOneToMany(const ROOT::RVec<ROOT::RVec<Tval>>& collections,
                                                  const ROOT::RVec<Tind1>& local_collIds,
                                                  const ROOT::RVec<Tind2>& order,
                                                  const ROOT::RVec<Tlims>& rec_begin,
                                                  const ROOT::RVec<Tlims>& rec_end)
    {
 
      auto object_vals = combineCollections(collections, local_collIds, order);
      return vecOfOneToMany(rec_begin, rec_end, object_vals);
    }

  } // namespace podio
} // namespace rad
