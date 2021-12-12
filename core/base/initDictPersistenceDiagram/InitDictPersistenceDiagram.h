#pragma once

//#include <PDClustering.h>
//#include <PersistenceDiagramBarycenter.h>
#include <Wrapper.h>
#include <algorithm>
#include <array>
#include <PersistenceDiagramAuction.h>

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
  class InitFarBorderDict : public Debug {

  public:
    InitFarBorderDict() {
      this->setDebugMsgPrefix("InitFarBorderDict");
    };


    void execute(std::vector<Diagram> &DictDiagrams,
                      const std::vector<Diagram> &datas,
                      int nbAtoms);


    // void executeAtoms(std::vector<Diagram> &DictDiagrams);

    // inline void setNbAtoms(const int nbAtoms) {
    // NbAtoms = nbAtoms;
    //}
    // inline void setDos(const bool min, const bool sad, const bool max) {
    //   do_min_ = min;
    //   do_sad_ = sad;
    //   do_max_ = max;
    // }
  protected:
    void
      setBidderDiagrams(const size_t nInputs,
                        std::vector<Diagram> &inputDiagrams,
                        std::vector<BidderDiagram<double>> &bidder_diags) const;

    double computeDistance(const BidderDiagram<double> &D1,
                           const BidderDiagram<double> &D2) const;


    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda{0};
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
    bool do_min_{true}, do_sad_{true}, do_max_{true};
  };
} // namespace ttk
