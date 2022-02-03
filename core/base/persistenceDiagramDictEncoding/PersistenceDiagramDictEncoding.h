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

#include <array>
#include <limits>

#include <InitDictBorder.h>
#include <InitDictRandomly.h>
#include <ConstrainedGradientDescent.h>
#include <PersistenceDiagramAuction.h>
#include <PersistenceDiagramClustering.h>
#include <Wrapper.h>

namespace ttk {
  class PersistenceDiagramDictEncoding : virtual public Debug {

  public:
    PersistenceDiagramDictEncoding() {
      this->setDebugMsgPrefix("PersistenceDiagramDictEncoding");
    }

    void execute(const std::vector<Diagram> &intermediateDiagrams,
                 std::vector<Diagram> &dictDiagrams,
                 std::vector<std::vector<double>> &vectorWeights,
                 const std::array<size_t, 2> &nInputs,
                 const int seed,
                 const int numAtom);

    enum class BACKEND{BORDER_INIT = 0 , RANDOM_INIT = 1 , FIRST_DIAGS = 2};

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
      const std::vector<BidderDiagram<double>> &bidder_diags) const;
    double computeDistance(const BidderDiagram<double> &D1,
                           const BidderDiagram<double> &D2,
                           std::vector<MatchingTuple> &matching) const;

    void computeGradientWeights(
      std::vector<double> &gradWeights,
      std::vector<Matrix> &hessianList,
      const std::vector<Diagram> &dictDiagrams,
      const std::vector<std::vector<MatchingTuple>> &matchingsAtoms,
      const Diagram &Barycenter,
      const Diagram &newData,
      const std::vector<MatchingTuple> &matchingsMin,
      const std::vector<MatchingTuple> &matchingsMax,
      const std::vector<MatchingTuple> &matchingsSad,
      const std::vector<size_t> &indexBaryMin,
      const std::vector<size_t> &indexBaryMax,
      const std::vector<size_t> &indexBarySad,
      const std::vector<size_t> &indexDataMin,
      const std::vector<size_t> &indexDataMax,
      const std::vector<size_t> &indexDataSad) const;

    void computeGradientAtoms(std::vector<Matrix> &gradsAtoms,
                              const std::vector<double> &weights,
                              const Diagram &Barycenter,
                              const Diagram &newData,
                              const std::vector<MatchingTuple> &matchingsMin,
                              const std::vector<MatchingTuple> &matchingsMax,
                              const std::vector<MatchingTuple> &matchingsSad,
                              const std::vector<size_t> &indexBaryMin,
                              const std::vector<size_t> &indexBaryMax,
                              const std::vector<size_t> &indexBarySad,
                              const std::vector<size_t> &indexDataMin,
                              const std::vector<size_t> &indexDataMax,
                              const std::vector<size_t> &indexDataSad,
                              std::vector<int> &checker) const;

    // A modifier
    void
      setBidderDiagrams(const size_t nInputs,
                        std::vector<Diagram> &inputDiagrams,
                        std::vector<BidderDiagram<double>> &bidder_diags) const;

    // A modifier
    void enrichCurrentBidderDiagrams(
      const std::vector<BidderDiagram<double>> &bidder_diags,
      std::vector<BidderDiagram<double>> &current_bidder_diags,
      const std::vector<double> &maxDiagPersistence) const;


    int InitDictionary(std::vector<ttk::Diagram> &dictDiagrams,
                       const std::vector<ttk::Diagram> &datas,
                       const int nbAtom,
                       bool do_min_,
                       bool do_sad_,
                       bool do_max_,
                       int seed);

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
    bool do_min_{true}, do_sad_{true}, do_max_{true};

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
