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
#include <PersistenceDiagramDictionary.h>
#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>

#include <DimensionReduction.h>

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
    enum class BACKEND { MDS = 0, DICTIONARY = 1 };

    PersistenceDiagramDictDecoding() {
      this->setDebugMsgPrefix("PersistenceDiagramDictDecoding");
    }

    void execute(std::vector<DiagramType> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 std::vector<DiagramType> &Barycenters) const;

  protected:
    BACKEND ProjMet{BACKEND::MDS};
    bool ProgBarycenter{false};

    void computeAtomsCoordinates(
      std::vector<ttk::DiagramType> &atoms,
      const std::vector<std::vector<double>> &vectorWeights,
      std::vector<std::array<double, 3>> &coords,
      std::vector<std::array<double, 3>> &trueCoords,
      std::vector<double> &xVector,
      std::vector<double> &yVector,
      std::vector<double> &zVector,
      const double spacing,
      const double maxPersistence,
      const size_t nAtoms) const;

  }; // PersistenceDiagramDictDecoding class

} // namespace ttk
