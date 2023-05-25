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
#include <PersistenceDiagramDictEncoding.h>
#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {
  using Matrice = std::vector<std::vector<double>>;
  using VectorMatchingTuple = std::vector<MatchingType>;

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

    void execute(std::vector<DiagramType> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 std::vector<DiagramType> &Barycenters) const;

  protected:
    bool ProgBarycenter{false};

    void computeAtomsCoordinates(
      std::vector<ttk::DiagramType> &atoms,
      const std::vector<std::vector<double>> &vectorWeights,
      std::vector<std::pair<double, double>> &coords,
      std::vector<std::pair<double, double>> &true_coords,
      std::vector<double> &xVector,
      std::vector<double> &yVector,
      const double spacing,
      const double max_persistence,
      const size_t nAtoms) const;

  }; // PersistenceDiagramDictDecoding class

} // namespace ttk
