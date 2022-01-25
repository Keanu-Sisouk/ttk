#pragma once

//#include <PDClustering.h>
//#include <PersistenceDiagramBarycenter.h>
#include <PersistenceDiagramAuction.h>
#include <Wrapper.h>
#include <algorithm>
#include <array>

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
  using Matrix = std::vector<std::vector<double>>;
  using MatchingTuple = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;
  class InitFarBorderDict : public Debug {

  public:
    InitFarBorderDict() {
      this->setDebugMsgPrefix("InitFarBorderDict");
    };

    void execute(std::vector<Diagram> &DictDiagrams,
                 const std::vector<Diagram> &datas,
                 const int nbAtoms,
                 bool do_min_,
                 bool do_sad_,
                 bool do_max_);

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
    double Lambda{0};
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
    // bool do_min_{true}, do_sad_{true}, do_max_{true};
  };
} // namespace ttk
