#pragma once

//#include <PDClustering.h>
//#include <PersistenceDiagramBarycenter.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramDistanceMatrix.h>
#include <PersistenceDiagramUtils.h>
#include <Wrapper.h>

#include <algorithm>
#include <array>

namespace ttk {
  using Matrix = std::vector<std::vector<double>>;

  class InitFarBorderDict : public Debug {

  public:
    InitFarBorderDict() {
      this->setDebugMsgPrefix("InitFarBorderDict");
    };

    void execute(std::vector<ttk::DiagramType> &DictDiagrams,
                 const std::vector<ttk::DiagramType> &datas,
                 const int nbAtoms,
                 bool do_min_,
                 bool do_sad_,
                 bool do_max_);

  protected:
    void setBidderDiagrams(const size_t nInputs,
                           std::vector<ttk::DiagramType> &inputDiagrams,
                           std::vector<BidderDiagram> &bidder_diags) const;

    double computeDistance(const BidderDiagram &D1,
                           const BidderDiagram &D2) const;

    int getNextIndex(const Matrix &distMatrix,
                     const std::vector<int> &indices) const;

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    double Lambda{0};
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
    // bool do_min_{true}, do_sad_{true}, do_max_{true};
  };
} // namespace ttk
