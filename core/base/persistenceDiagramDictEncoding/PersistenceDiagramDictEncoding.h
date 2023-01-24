/// \ingroup base
/// \class ttk::PersistenceDiagramDictEncoding
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering

#pragma once

#include "PersistenceDiagramUtils.h"
#include <array>
#include <limits>

#include <ConstrainedGradientDescent.h>
#include <InitDictBorder.h>
#include <InitDictRandomly.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramClustering.h>
#include <Wrapper.h>

namespace ttk {
  class PersistenceDiagramDictEncoding : virtual public Debug {

  public:
    PersistenceDiagramDictEncoding() {
      this->setDebugMsgPrefix("PersistenceDiagramDictEncoding");
    }

    void execute(std::vector<ttk::DiagramType> &intermediateDiagrams,
                 const std::vector<ttk::DiagramType> &intermediateAtoms,
                 std::vector<ttk::DiagramType> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 const std::array<size_t, 2> &nInputs,
                 const int seed,
                 const int numAtom,
                 std::vector<double> &loss_tab,
                 std::vector<double> &timers,
                 std::vector<double> &true_loss_tab,
                 std::vector<std::vector<double>> &allLosses,
                 double percent_);

    void method(const std::vector<ttk::DiagramType> &intermediateDiagrams,
                const std::vector<ttk::DiagramType> &intermediateAtoms,
                std::vector<ttk::DiagramType> &dictDiagrams,
                std::vector<std::vector<double>> &vectorWeights,
                const std::array<size_t, 2> &nInputs,
                const int seed,
                const int numAtom,
                std::vector<double> &loss_tab,
                std::vector<double> &true_loss_tab,
                std::vector<double> &timers,
                std::vector<std::vector<double>> &allLosses,
                std::vector<std::vector<double>> &histoVectorWeights,
                std::vector<ttk::DiagramType> &histoDictDiagrams,
                bool preWeightOpt,
                double acc,
                std::vector<BidderDiagram> &true_bidder_diagram_min,
                std::vector<BidderDiagram> &true_bidder_diagram_sad,
                std::vector<BidderDiagram> &true_bidder_diagram_max,
                Timer &tm_method,
                double percent_,
                bool do_compression);

    enum class BACKEND {
      BORDER_INIT = 0,
      RANDOM_INIT = 1,
      FIRST_DIAGS = 2,
      INPUT_ATOMS = 3,
      GREEDY_INIT = 4
    };

    inline void setWasserstein(const int data) {
      Wasserstein = data;
    }
    inline void setDos(const bool min, const bool sad, const bool max) {
      do_min_ = min;
      do_sad_ = sad;
      do_max_ = max;
    }
    inline void setAlpha(const double alpha) {
      Alpha = alpha;
    }
    inline void setLambda(const double lambda) {
      Lambda = lambda;
    }
    inline void setDeltaLim(const double deltaLim) {
      DeltaLim = deltaLim;
    }
    inline void setMaxNumberOfPairs(const size_t data) {
      MaxNumberOfPairs = data;
    }
    inline void setMinPersistence(const double data) {
      MinPersistence = data;
    }
    inline void setConstraint(const int data) {
      if(data == 0) {
        this->Constraint = ConstraintType::FULL_DIAGRAMS;
      } else if(data == 1) {
        this->Constraint = ConstraintType::NUMBER_PAIRS;
      } else if(data == 2) {
        this->Constraint = ConstraintType::ABSOLUTE_PERSISTENCE;
      } else if(data == 3) {
        this->Constraint = ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG;
      } else if(data == 4) {
        this->Constraint = ConstraintType::RELATIVE_PERSISTENCE_GLOBAL;
      }
    }

  protected:
    BACKEND BackEnd{BACKEND::BORDER_INIT};
    double distVect(const std::vector<double> &vec1,
                    const std::vector<double> &vec2) const;

    double getMostPersistent(
      const std::vector<BidderDiagram> &bidder_diags) const;
    double computeDistance(const BidderDiagram &D1,
                           const BidderDiagram &D2,
                           std::vector<ttk::MatchingType> &matching) const;

    void computeGradientWeights(
      std::vector<double> &gradWeights,
      std::vector<Matrix> &hessianList,
      const std::vector<ttk::DiagramType> &dictDiagrams,
      const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsMin,
      const std::vector<ttk::MatchingType> &matchingsMax,
      const std::vector<ttk::MatchingType> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad,
      const bool do_optimizeAtoms) const;

    void computeExplicitSolutionWeights(
      std::vector<double> &weights,
      const std::vector<ttk::DiagramType> &dictDiagrams,
      const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsMin,
      const std::vector<ttk::MatchingType> &matchingsMax,
      const std::vector<ttk::MatchingType> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad,
      const bool do_optimizeAtoms) const;

    void computeGradientAtoms(
      std::vector<Matrix> &gradsAtoms,
      const std::vector<double> &weights,
      const ttk::DiagramType &Barycenter,
      const ttk::DiagramType &newData,
      const std::vector<ttk::MatchingType> &matchingsMin,
      const std::vector<ttk::MatchingType> &matchingsMax,
      const std::vector<ttk::MatchingType> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad,
      std::vector<int> &checker,
      std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
      ttk::DiagramType &infoToAdd,
      int nbDiags,
      bool do_DimReduct) const;
    
    void computeExplicitSolution(
        std::vector<ttk::DiagramType> &DictDiagrams,
        const std::vector<std::vector<double>> &vectorWeights,
        const std::vector<std::vector<Matrix>> &gradsAtomsList,
        const std::vector<std::vector<std::vector<ttk::MatchingType>>> &allMatchings,
        const std::vector<ttk::DiagramType> &Barycenters,
        const int nb_points,
        const std::vector<std::vector<int>> &checkerAtomsList,
        std::vector<std::vector<std::vector<int>>> &allProjForDiag,
        std::vector<ttk::DiagramType> &allFeaturesToAdd,
        std::vector<std::vector<std::array<double,2 >>> &allProjLocations,
        std::vector<std::vector<std::vector<double>>> &allVectorForProjContrib,
        std::vector<std::vector<std::vector<std::array<double, 2>>>> &allPairToAddGradList,
        std::vector<ttk::DiagramType> &allInfoToAdd);

    // A modifier
    void
      setBidderDiagrams(const size_t nInputs,
                        std::vector<ttk::DiagramType> &inputDiagrams,
                        std::vector<BidderDiagram> &bidder_diags) const;

    // A modifier
    void enrichCurrentBidderDiagrams(
      const std::vector<BidderDiagram> &bidder_diags,
      std::vector<BidderDiagram> &current_bidder_diags,
      const std::vector<double> &maxDiagPersistence) const;

    int InitDictionary(std::vector<ttk::DiagramType> &dictDiagrams,
                       const std::vector<ttk::DiagramType> &datas,
                       const std::vector<ttk::DiagramType> &inputAtoms,
                       const int nbAtom,
                       bool do_min_,
                       bool do_sad_,
                       bool do_max_,
                       int seed,
                       double percent_);

    void gettingBidderDiagrams(
      const std::vector<ttk::DiagramType> &intermediateDiagrams,
      std::vector<BidderDiagram> &bidder_diagrams_min,
      std::vector<BidderDiagram> &bidder_diagrams_sad,
      std::vector<BidderDiagram> &bidder_diagrams_max);

    void controlAtomsSize(
      const std::vector<ttk::DiagramType> &intermediateDiagrams,
      std::vector<ttk::DiagramType> &dictDiagrams
    ) const;

    void controlAtomsSize2(
      const std::vector<ttk::DiagramType> &intermediateDiagrams,
      std::vector<ttk::DiagramType> &dictDiagrasm
    ) const;

    double getMaxPers(const ttk::DiagramType &data);

    int Wasserstein{2};
    double Alpha{1.0};
    double DeltaLim{0.01};
    // lambda : 0<=lambda<=1
    // parametrizes the point used for the physical (critical) coordinates of
    // the persistence paired lambda = 1 : extremum (min if pair min-sad, max if
    // pair sad-max) lambda = 0 : saddle (bad stability) lambda = 1/2 : middle
    // of the 2 critical points of the pair
    double Lambda;
    size_t MaxNumberOfPairs{20};
    double MinPersistence{0.1};
    double CompressionFactor{1.5};
    bool do_min_{true}, do_sad_{true}, do_max_{true};

    int maxLag2;

    int MaxEpoch;
    bool MaxEigenValue{true};
    bool OptimizeWeights{true};
    bool OptimizeAtoms{true};

    bool CreationFeatures{true};
    bool explicitSolWeights{false};
    bool explicitSol{false};
    bool Fusion{false};
    bool ProgBarycenter{false};

    bool sortedForTest{false};
    bool ProgApproach{false};
    bool StopCondition{true};

    bool CompressionMode{false};
    bool DimReductMode{false};

    enum class ConstraintType {
      FULL_DIAGRAMS,
      NUMBER_PAIRS,
      ABSOLUTE_PERSISTENCE,
      RELATIVE_PERSISTENCE_PER_DIAG,
      RELATIVE_PERSISTENCE_GLOBAL,
    };
    ConstraintType Constraint{ConstraintType::RELATIVE_PERSISTENCE_GLOBAL};
  };
} // namespace ttk
