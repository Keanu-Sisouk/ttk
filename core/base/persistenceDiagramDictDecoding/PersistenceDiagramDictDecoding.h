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
  }; // PersistenceDiagramDictDecoding class

} // namespace ttk
