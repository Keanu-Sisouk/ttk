/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::PersistenceDiagramDictDecoding
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %PersistenceDiagramDictDecoding class that computes
/// for each vertex of a triangulation the average scalar value of itself and
/// its direct neighbors.
///
/// \b Related \b publication: \n
/// 'PersistenceDiagramDictDecoding'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <PersistenceDiagramClustering.h>
#include <Wrapper.h>

namespace ttk {
  using DiagramTuple = std::tuple<
    /** Vertex Id of low pair element */
    ttk::SimplexId,
    /** Critical Type of low pair element */
    ttk::CriticalType,
    /** Vertex Id of high pair element */
    ttk::SimplexId,
    /** Critical Type of high pair element */
    ttk::CriticalType,
    /** Pair persistence value */
    double,
    /** Pair type */
    ttk::SimplexId,
    /** Pair birth */
    double,
    /** Low pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float,
    /** Pair death */
    double,
    /** High pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float>;

  using Diagram = std::vector<DiagramTuple>;
  using Matrice = std::vector<std::vector<double>>;
  using MatchingTuple = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;
  using VectorMatchingTuple = std::vector<MatchingTuple>;

  /**
   * The PersistenceDiagramDictDecoding class provides methods to compute for
   * each vertex of a triangulation the average scalar value of itself and its
   * direct neighbors.
   */
  class PersistenceDiagramDictDecoding : virtual public Debug {

  public:
    PersistenceDiagramDictDecoding() {
      this->setDebugMsgPrefix("PersistenceDiagramDictDecoding");
    }

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    void execute(std::vector<Diagram> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 std::vector<Diagram> &Barycenters) const;

  }; // PersistenceDiagramDictDecoding class

} // namespace ttk
