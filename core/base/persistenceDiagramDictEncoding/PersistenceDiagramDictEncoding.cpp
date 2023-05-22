#include "PersistenceDiagramUtils.h"
#include <algorithm>
#include <cmath>
#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif

#include <PersistenceDiagramDictEncoding.h>



using namespace ttk;


void PersistenceDiagramDictEncoding::execute(
  std::vector<ttk::DiagramType> &intermediateDiagrams,
  const std::vector<ttk::DiagramType> &intermediateAtoms,
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const int seed,
  const int numAtom,
  std::vector<double> &loss_tab,
  std::vector<double> &timers,
  std::vector<double> &true_loss_tab,
  std::vector<std::vector<double>> &allLosses,
  double percent) {

  if(!ProgApproach_) {

    for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
      std::cout << "SIZE OF DIAG " << i
                << " IS: " << intermediateDiagrams[i].size() << std::endl;
      if(sortedForTest_) {
        auto &diag = intermediateDiagrams[i];
        std::sort(diag.begin(), diag.end(),
                  [](const PersistencePair &t1, const PersistencePair &t2) {
                    return (t1.death.sfValue - t1.birth.sfValue)
                           > (t2.death.sfValue - t2.birth.sfValue);
                  });
      }
    }

    bool do_compression = false;
    if(CompressionMode_) {
      do_compression = true;
    }

    std::vector<BidderDiagram> trueBidderDiagramMin(intermediateDiagrams.size());
    std::vector<BidderDiagram> trueBidderDiagramSad(intermediateDiagrams.size());
    std::vector<BidderDiagram> trueBidderDiagramMax(intermediateDiagrams.size());

    std::vector<std::vector<double>> histoVectorWeights(
      intermediateDiagrams.size());
    std::vector<ttk::DiagramType> histoDictDiagrams(numAtom);
    this->maxLag2_ = 5;

    Timer tm_method{};
    Timer tm_init{};
    bool preWeightOpt = false;
    initDictionary(dictDiagrams, intermediateDiagrams, intermediateAtoms,
                   numAtom, this->do_min_, this->do_sad_, this->do_max_, seed, percent);
    this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
                   threadNumber_, debug::LineMode::NEW);
    if(true){
    // if(CompressionMode && !CreationFeatures){
      controlAtomsSize(intermediateDiagrams, dictDiagrams);
    }
    method(intermediateDiagrams, dictDiagrams, vectorWeights, seed,
           numAtom, loss_tab, true_loss_tab, timers, allLosses,
           histoVectorWeights, histoDictDiagrams, preWeightOpt, 0.01,
           trueBidderDiagramMin, trueBidderDiagramSad, trueBidderDiagramMax,
           tm_method, percent, do_compression);
  } else {

    bool do_compression = false;
    for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
      auto &diag = intermediateDiagrams[i];
      std::sort(diag.begin(), diag.end(),
                [](const PersistencePair &t1, const PersistencePair &t2) {
                  return (t1.death.sfValue - t1.birth.sfValue)
                         > (t2.death.sfValue - t2.birth.sfValue);
                });
    }

    std::vector<BidderDiagram> trueBidderDiagramMin{};
    std::vector<BidderDiagram> trueBidderDiagramSad{};
    std::vector<BidderDiagram> trueBidderDiagramMax{};

    gettingBidderDiagrams(intermediateDiagrams, trueBidderDiagramMin,
                          trueBidderDiagramSad, trueBidderDiagramMax);

    // std::vector<double> percentages{0.2 , 0.15 , 0.1 , 0.05};
    // std::vector<double> percentages{0.8 , 0.6 , 0.5, 0.4, 0.3 , 0.2};
    // std::vector<double> percentages{0.3 , 0.2 , 0.1 , 0.05, 0.01};
    int start = 20;
    double stop = percent;
    std::vector<double> percentages;
    for(int value = start; value > stop; value -= 5) {
      percentages.push_back(static_cast<double>(value) / 100.);
    }
    percentages.push_back(static_cast<double>(percent) / 100.);
    for(size_t k = 0; k < percentages.size(); ++k) {
      std::cout << "PERCENT " << percentages[k] << std::endl;
    }
    // std::vector<double> percentages{0.4};
    std::vector<std::vector<double>> histoVectorWeights(
      intermediateDiagrams.size());
    std::vector<ttk::DiagramType> histoDictDiagrams(numAtom);
    std::vector<ttk::DiagramType> dataTemp(intermediateDiagrams.size());
    bool preWeightOpt = true;
    Timer tm_method{};
    for(size_t j = 0; j < 1; ++j) {
      double percentage = percentages[j];
      // std::vector<Diagram> dataTemp(intermediateDiagrams.size());
      // for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      // auto diag = intermediateDiagrams[i];
      // dataTemp[i] = diag;
      //}

      this->maxLag2_ = 0;

      for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
        auto &diag = intermediateDiagrams[i];
        auto &t = diag[0];
        // double maxPers = t.death.sfValue - t.birth.sfValue;
        double maxPers = getMaxPers(diag);
        std::cout << "MAX PERS" << maxPers << std::endl;
        auto &diagTemp = dataTemp[i];
        dataTemp[i].push_back(t);
        for(size_t p = 1; p < diag.size(); ++p) {
          auto &t2 = diag[p];
          if(percentage * maxPers <= (t2.death.sfValue - t2.birth.sfValue)) {
            dataTemp[i].push_back(t2);
          } else {
            continue;
          }
        }
        // double maxPers = t.death.sfValue - t.birth.sfValue;
        // diag.erase(std::remove_if(diag.begin() , diag.end(), [maxPers,
        // percentage](DiagramTuple &t){ return (t.death.sfValue -
        // t.birth.sfValue) < percentage*maxPers;}), diag.end());
        std::cout << "SIZE OF DIAG " << i << " IS: " << diagTemp.size()
                  << std::endl;
      }

      Timer tm_init{};
      initDictionary(dictDiagrams, dataTemp, intermediateAtoms, numAtom,
                     this->do_min_, this->do_sad_, this->do_max_, seed, percent);
      this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
                     threadNumber_, debug::LineMode::NEW);

      method(dataTemp, dictDiagrams, vectorWeights, seed, numAtom,
             loss_tab, true_loss_tab, timers, allLosses, histoVectorWeights,
             histoDictDiagrams, preWeightOpt, 0.01, trueBidderDiagramMin,
             trueBidderDiagramSad, trueBidderDiagramMax, tm_method, percent,
             do_compression);
    }

    // int min_pairs_to_add = 0;
    std::vector<int> sizeCheck(intermediateDiagrams.size(), 0);
    int sum = 0;
    for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
      sum += sizeCheck[i];
    }
    // std::vector<int> newPercentages{1, 5 , 25, 50 , 100};
    // std::vector<int> newPercentages(10 , 10);
    // int q = 1;
    // while(sum != static_cast<int>(intermediateDiagrams.size())){
    // for(size_t j = 1 ; j < newPercentages.size() ; ++j){
    // preWeightOpt = false;
    int counter = 0;
    for(size_t j = 1; j < percentages.size(); ++j) {
      if(j == percentages.size() - static_cast<size_t>(1) && CompressionMode_) {
        do_compression = true;
      }
      double percentage = percentages[j];
      double previousPerc = percentages[j - 1];
      if(j < percentages.size() - 1) {
        this->maxLag2_ = 0;
      } else {
        this->maxLag2_ = 5;
      }

      // std::cout << "PERCENTAGE " << percentage
      //  std::vector<Diagram> dataTemp(intermediateDiagrams.size());
      //  for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      //  auto diag = intermediateDiagrams[i];
      //  dataTemp[i] = diag;
      // }
      for(size_t i = 0; i < intermediateDiagrams.size(); ++i) {
        // if(sizeCheck[i] == 1){
        // continue;
        //}
        auto &diag = intermediateDiagrams[i];
        auto &diagTemp = dataTemp[i];
        // int m = diag.size();

        // std::cout << "SIZE ORIGINAL " << m << std::endl;
        int n = diagTemp.size();
        // int counter = 0;
        // auto &lastTuple = diagTemp[n - 1];
        // int max_pairs_to_add = 10 + m*10/100;
        // int max_pairs_to_add = 5 + m * percentage / 100;
        // std::cout << "MAX PAIRS TO ADD " << max_pairs_to_add << std::endl;
        // std::cout << "TRIPLE WHAT" << m*q/100 << std::endl;
        // std::cout << "WHAT " << static_cast<int> (m*(q/100)) << std::endl;
        // double previousPers = lastTuple.death.sfValue - lastTuple.birth.sfValue;
        // auto &t = diag[0];
        double maxPers = getMaxPers(diag);
        // double maxPers = t.death.sfValue - t.birth.sfValue;
        // auto &diagTemp = dataTemp[i];
        // dataTemp[i].push_back(t);
        for(size_t p = 0; p < diag.size(); ++p) {
          auto &t2 = diag[p];

          if(percentage * maxPers <= (t2.death.sfValue - t2.birth.sfValue)
             && (t2.death.sfValue - t2.birth.sfValue)
                  < previousPerc * maxPers) {
            dataTemp[i].push_back(t2);
            counter += 1;
          }

          // if((counter <= max_pairs_to_add) && ((t2.death.sfValue -
          // t2.birth.sfValue) < previousPers)) { dataTemp[i].push_back(t2);
          // counter += 1;
          //}
        }

        // double maxPers = t.death.sfValue - t.birth.sfValue;
        // diag.erase(std::remove_if(diag.begin() , diag.end(), [maxPers,
        // percentage](DiagramTuple &t){ return (t.death.sfValue -
        // t.birth.sfValue) < percentage*maxPers;}), diag.end());
        std::cout << "SIZE OF DIAG " << i << " IS: " << diagTemp.size()
                  << std::endl;
      }
      if(counter == 0) {
        continue;
      }
      method(dataTemp, dictDiagrams, vectorWeights, seed, numAtom,
             loss_tab, true_loss_tab, timers, allLosses, histoVectorWeights,
             histoDictDiagrams, preWeightOpt, 0.01, trueBidderDiagramMin,
             trueBidderDiagramSad, trueBidderDiagramMax, tm_method, percent,
             do_compression);
      // sum = 0;
      // for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      // sum += sizeCheck[i];
      //}

      // q += 20;
    }

    this->printMsg(
      "Total time", 1.0, tm_method.getElapsedTime(), threadNumber_);
  }
}

void PersistenceDiagramDictEncoding::method(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
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
  std::vector<BidderDiagram> &trueBidderDiagramMin,
  std::vector<BidderDiagram> &trueBidderDiagramSad,
  std::vector<BidderDiagram> &trueBidderDiagramMax,
  Timer &tm_method,
  double percent,
  bool do_compression) {

  Timer tm{};
  double tm_part = 0.;

  bool do_optimizeAtoms = false;
  bool do_optimizeWeights = false;

  if(OptimizeWeights_) {
    printMsg("Weight Optimization activated");
    do_optimizeWeights = true;
  } else {
    printWrn("Weight Optimization desactivated");
  }
  if(OptimizeAtoms_) {
    do_optimizeAtoms = true;
    printMsg("Atom Optimization activated");
  } else {
    printWrn("Atom Optimization desactivated");
  }

  const auto nDiags = intermediateDiagrams.size();

  if(do_min_ && do_sad_ && do_max_) {
    this->printMsg("Processing all critical pairs types");
  } else if(do_min_) {
    this->printMsg("Processing only MIN-SAD pairs");
  } else if(do_sad_) {
    this->printMsg("Processing only SAD-SAD pairs");
  } else if(do_max_) {
    this->printMsg("Processing only SAD-MAX pairs");
  }

  // inputDiagrams = newDatas here
  // tracking the original indices
  std::vector<ttk::DiagramType> inputDiagramsMin(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsSad(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsMax(nDiags);

  std::vector<BidderDiagram> bidderDiagramsMin{};
  std::vector<BidderDiagram> bidderDiagramsSad{};
  std::vector<BidderDiagram> bidderDiagramsMax{};

  std::vector<std::vector<size_t>> originIndexDatasMin(nDiags);
  std::vector<std::vector<size_t>> originIndexDatasSad(nDiags);
  std::vector<std::vector<size_t>> originIndexDatasMax(nDiags);

  // std::vector<BidderDiagram> current_bidderDiagramsMin{};
  // std::vector<BidderDiagram> current_bidderDiagramsSad{};
  // std::vector<BidderDiagram> current_bidderDiagramsMax{};

  // Store the persistence of the global min-max pair
  // std::vector<double> maxDiagPersistence(nDiags);

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; i++) {
    const ttk::DiagramType &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const ttk::PersistencePair &t = CTDiagram[j];
      const ttk::CriticalType nt1 = t.birth.type;
      const ttk::CriticalType nt2 = t.death.type;
      const double pers = t.persistence();
      // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMin[i].emplace_back(t);
          originIndexDatasMin[i].push_back(j);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            originIndexDatasMax[i].push_back(j);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
            originIndexDatasMin[i].push_back(j);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
            originIndexDatasSad[i].push_back(j);
          }
        }
      }
    }
  }

  if(this->do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidderDiagramsMin);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidderDiagramsSad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidderDiagramsMax);
  }

  switch(this->Constraint) {
    case ConstraintType::FULL_DIAGRAMS:
      this->printMsg("Using all diagram pairs");
      break;
    case ConstraintType::NUMBER_PAIRS:
      this->printMsg("Using the " + std::to_string(this->MaxNumberOfPairs)
                     + " most persistent pairs");
      break;
    case ConstraintType::ABSOLUTE_PERSISTENCE: {
      std::stringstream pers{};
      pers << std::fixed << std::setprecision(2) << this->MinPersistence_;
      this->printMsg("Using diagram pairs above a persistence threshold of "
                     + pers.str());
    } break;
    case ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence_)))
        + "% most persistent pairs of every diagram");
      break;
    case ConstraintType::RELATIVE_PERSISTENCE_GLOBAL:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence_)))
        + "% most persistent pairs of all diagrams");
      break;
  }

  // std::vector<std::vector<double>> distMat{};

  std::vector<ttk::DiagramType> barycentersList(nDiags);
  std::vector<std::vector<std::vector<ttk::MatchingType>>> allMatchingsAtoms(
    nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasMin(nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasSad(nDiags);
  std::vector<std::vector<ttk::MatchingType>> matchingsDatasMax(nDiags);
  ConstrainedGradientDescent gradActor;
  double loss;
  double true_loss = 0.;
  // double loss1;
  // int epoch = 1;
  // std::vector<double> loss_tab;
  int lag = 0;
  int lag2 = 0;
  int lag3 = 0;
  int lagLimit = 10;
  int MIN_EPOCH = 20;
  int MAX_EPOCH = MaxEpoch_;
  bool cond = true;
  int epoch = 0;
  int nbEpochPrevious = loss_tab.size();
  std::vector<size_t> initSizes(dictDiagrams.size());
  for(size_t j = 0; j < dictDiagrams.size(); ++j) {
    initSizes[j] = dictDiagrams[j].size();
  }

  double factEquiv = static_cast<double>(numAtom);
  // double factEquiv = 1.;
  // double step = 1. / (sqrt(factEquiv) * 1e1);
  double step = 1. / (2. * 2. * factEquiv);
  gradActor.setStep(factEquiv);

  // std::vector<ttk::DiagramType> histoDictDiagrams(dictDiagrams.size());
  // BUFFERS
  std::vector<std::vector<int>> bufferHistoAllEpochLife(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferHistoAllBoolLife(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckUnderDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> bufferCheckAboveGlobal(dictDiagrams.size());

  std::vector<std::vector<int>> histoAllEpochLife(dictDiagrams.size());
  std::vector<std::vector<bool>> histoAllBoolLife(dictDiagrams.size());
  std::vector<std::vector<bool>> checkUnderDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkAboveGlobal(dictDiagrams.size());
  // std::vector<std::vector<double>> histoVectorWeights(nDiags);
  //  std::vector<double> allLosses(nDiags , 0.);
  std::vector<double> allLossesAtEpoch(nDiags, 0.);
  std::vector<double> trueAllLossesAtEpoch(nDiags, 0.);

  // bool condition = true;
  while(epoch < MAX_EPOCH && cond) {
    // for(int epoch = 1; epoch < MAX_EPOCH; ++epoch) {

    loss = 0.;
    true_loss = 0.;
    // auto vectorWeightsOld = vectorWeights;
    // this->printMsg("Epoch: " + std::to_string(epoch));
    Timer tm_it{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    //////////////////////////////WEIGHTS///////////////////////////////////
    for(size_t i = 0; i < nDiags; ++i) {
      auto &barycenter = barycentersList[i];
      std::vector<double> &weight = vectorWeights[i];
      std::vector<std::vector<ttk::MatchingType>> &matchings
        = allMatchingsAtoms[i];
      computeWeightedBarycenter(
        dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter_);
    }
    this->printMsg(
      "Computed 1st Barycenters for epoch " + std::to_string(epoch),
      epoch / static_cast<double>(MAX_EPOCH), tm_it.getElapsedTime(),
      threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);
    tm_part += static_cast<double>(tm_it.getElapsedTime());
   
    //   "====================BARYCENTER FINISHED======================");

    // tracking the original indices
    std::vector<ttk::DiagramType> barycentersListMin(nDiags);
    std::vector<ttk::DiagramType> barycentersListSad(nDiags);
    std::vector<ttk::DiagramType> barycentersListMax(nDiags);

    std::vector<BidderDiagram> bidderBarycentersListMin{};
    std::vector<BidderDiagram> bidderBarycentersListSad{};
    std::vector<BidderDiagram> bidderBarycentersListMax{};

    std::vector<std::vector<size_t>> originIndexBarysMin(nDiags);
    std::vector<std::vector<size_t>> originIndexBarysSad(nDiags);
    std::vector<std::vector<size_t>> originIndexBarysMax(nDiags);

    computeAllDistances(barycentersList, nDiags, barycentersListMin, barycentersListSad, barycentersListMax,
      bidderBarycentersListMin, bidderBarycentersListSad, bidderBarycentersListMax,
      originIndexBarysMin, originIndexBarysSad, originIndexBarysMax,
      bidderDiagramsMin, bidderDiagramsMax, bidderDiagramsSad,
      matchingsDatasMin, matchingsDatasMax, matchingsDatasSad, trueBidderDiagramMin,
      trueBidderDiagramSad, trueBidderDiagramMax, allLossesAtEpoch,
      trueAllLossesAtEpoch, true);

    for(size_t p = 0; p < nDiags; ++p) {
      loss += allLossesAtEpoch[p];
    }

    for(size_t p = 0; p < nDiags; ++p) {
      true_loss += trueAllLossesAtEpoch[p];
    }

    for(size_t p = 0; p < nDiags; ++p) {
      if(!ProgApproach_) {
        allLosses[p].push_back(allLossesAtEpoch[p]);
      } else {
        allLosses[p].push_back(trueAllLossesAtEpoch[p]);
      }
    }

    loss_tab.push_back(loss);

    printMsg(
      " Epoch " + std::to_string(epoch) + ", loss = " + std::to_string(loss), 1,
      threadNumber_, ttk::debug::LineMode::REPLACE);
    true_loss_tab.push_back(true_loss);
    timers.push_back(tm_method.getElapsedTime());

    // if(OptimizeWeights && OptimizeAtoms){
    if(preWeightOpt && OptimizeAtoms_) {
      if(epoch < 10) {
        do_optimizeAtoms = false;
      } else {
        do_optimizeAtoms = true;
      }
    }

    if(loss < 1e-7) {
      cond = false;
      do_optimizeWeights = false;
      do_optimizeAtoms = false;
    }

    double mini = *std::min_element(
      loss_tab.begin() + nbEpochPrevious, loss_tab.end() - 1);
    if(loss <= mini) {
      for(size_t p = 0; p < dictDiagrams.size(); ++p) {
        const auto atom = dictDiagrams[p];
        histoDictDiagrams[p] = atom;

        const auto &histoEpochAtom = histoAllEpochLife[p];
        const auto &histoBoolAtom = histoAllBoolLife[p];
        const auto &boolUnderDiag = checkUnderDiag[p];
        const auto &boolDiag = checkDiag[p];
        const auto &boolAboveGlobal = checkAboveGlobal[p];

        bufferHistoAllEpochLife[p] = histoEpochAtom;
        bufferHistoAllBoolLife[p] = histoBoolAtom;
        bufferCheckUnderDiag[p] = boolUnderDiag;
        bufferCheckDiag[p] = boolDiag;
        bufferCheckAboveGlobal[p] = boolAboveGlobal;

      }
      for(size_t p = 0; p < nDiags; ++p) {
        const auto weights = vectorWeights[p];
        histoVectorWeights[p] = weights;
      }
      lag = 0;
      lag3 = 0;
    } else {
      if(epoch > MIN_EPOCH) {
        lag += 1;
      } else {
        lag = 0;
      }
    }


    if((epoch > MIN_EPOCH) && (loss_tab[epoch + nbEpochPrevious] / loss_tab[epoch + nbEpochPrevious - 1] > 0.99)) {
      if(loss_tab[epoch + nbEpochPrevious] < loss_tab[epoch + nbEpochPrevious - 1]) {
        // if(true){
        if(lag2 == this->maxLag2_) {
          // lag = 0;
          this->printMsg("Loss not decreasing enough");
          if(StopCondition_) {
            for(size_t p = 0; p < dictDiagrams.size(); ++p) {
              const auto atom = histoDictDiagrams[p];
              dictDiagrams[p] = atom;
            }
            for(size_t p = 0; p < nDiags; ++p) {
              const auto weights = histoVectorWeights[p];
              vectorWeights[p] = weights;
            }
            do_optimizeWeights = false;
            do_optimizeAtoms = false;
            cond = false;
          }
        } else {
          lag2 += 1;
        }
      }
    } else {
      lag2 = 0;
    }

    if(epoch > MIN_EPOCH && lag > lagLimit) {

      if(StopCondition_) {
        for(size_t p = 0; p < dictDiagrams.size(); ++p) {
          const auto atom = histoDictDiagrams[p];
          dictDiagrams[p] = atom;
        }
        for(size_t p = 0; p < nDiags; ++p) {
          const auto weights = histoVectorWeights[p];
          vectorWeights[p] = weights;
        }
        this->printMsg("Minimum not passed");
        // if(StopCondition){
        do_optimizeWeights = false;
        do_optimizeAtoms = false;

        cond = false;
      }
    }

    if(cond  && (epoch > 0) &&  (loss_tab[epoch + nbEpochPrevious] > 2. * mini)) {
      lag3 += 1;
      lag2 = 0;


      if(lag3 > 2) {
        this->printMsg("Loss increasing too much, reducing step and recompute Barycenters and matchings");
        for(size_t p = 0; p < dictDiagrams.size(); ++p) {
          const auto atom = histoDictDiagrams[p];
          dictDiagrams[p] = atom;

          const auto &bufferHistoEpochAtom = bufferHistoAllEpochLife[p];
          const auto &bufferHistoBoolAtom = bufferHistoAllBoolLife[p];
          const auto &bufferBoolUnderDiag = bufferCheckUnderDiag[p];
          const auto &bufferBoolDiag = bufferCheckDiag[p];
          const auto &bufferBoolAboveGlobal = bufferCheckAboveGlobal[p];

          histoAllEpochLife[p] = bufferHistoEpochAtom;
          histoAllBoolLife[p] = bufferHistoBoolAtom;
          checkUnderDiag[p] = bufferBoolUnderDiag;
          checkDiag[p] = bufferBoolDiag;
          checkAboveGlobal[p] = bufferBoolAboveGlobal;
        }


        for(size_t p = 0; p < nDiags; ++p) {
          const auto weights = histoVectorWeights[p];
          vectorWeights[p] = weights;
        }
        gradActor.reduceStep();
        step = step/2.;
        lag3 = 0;


        barycentersList.clear();
        barycentersList.resize(nDiags);
        allMatchingsAtoms.clear();
        allMatchingsAtoms.resize(nDiags);

        barycentersListMin.clear();
        barycentersListSad.clear();
        barycentersListMax.clear();

        barycentersListMin.resize(nDiags);
        barycentersListSad.resize(nDiags);
        barycentersListMax.resize(nDiags);

        bidderBarycentersListMin.clear();
        bidderBarycentersListSad.clear();
        bidderBarycentersListMax.clear();

        originIndexBarysMin.clear();
        originIndexBarysSad.clear();
        originIndexBarysMax.clear();

        originIndexBarysMin.resize(nDiags);
        originIndexBarysSad.resize(nDiags);
        originIndexBarysMax.resize(nDiags);

        matchingsDatasMin.clear();
        matchingsDatasSad.clear();
        matchingsDatasMax.clear();

        matchingsDatasMin.resize(nDiags);
        matchingsDatasSad.resize(nDiags);
        matchingsDatasMax.resize(nDiags);


#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      //////////////////////////////WEIGHTS///////////////////////////////////
        for(size_t i = 0; i < nDiags; ++i) {
          auto &barycenter = barycentersList[i];
          std::vector<double> &weight = vectorWeights[i];
          std::vector<std::vector<ttk::MatchingType>> &matchings
            = allMatchingsAtoms[i];
          computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings,
                                    *this, ProgBarycenter_);
        }

        computeAllDistances(barycentersList, nDiags, barycentersListMin, barycentersListSad, barycentersListMax,
        bidderBarycentersListMin, bidderBarycentersListSad, bidderBarycentersListMax,
        originIndexBarysMin, originIndexBarysSad, originIndexBarysMax,
        bidderDiagramsMin, bidderDiagramsMax, bidderDiagramsSad,
        matchingsDatasMin, matchingsDatasMax, matchingsDatasSad, trueBidderDiagramMin,
        trueBidderDiagramSad, trueBidderDiagramMax, allLossesAtEpoch,
        trueAllLossesAtEpoch, false);
      }
    }

    std::vector<std::vector<Matrix>> allHessianLists(nDiags);
    std::vector<std::vector<double>> gradWeightsList(nDiags);
    Timer tm_opt1{};

    // WEIGHT OPTIMIZATION
    if(do_optimizeWeights) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &gradWeights = gradWeightsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const auto &Barycenter = barycentersList[i];
        const auto &Data = intermediateDiagrams[i];
        std::vector<Matrix> &hessianList = allHessianLists[i];
        const std::vector<ttk::MatchingType> &matchingsMin
          = matchingsDatasMin[i];
        const std::vector<ttk::MatchingType> &matchingsMax
          = matchingsDatasMax[i];
        const std::vector<ttk::MatchingType> &matchingsSad
          = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = originIndexBarysMin[i];
        const std::vector<size_t> &indexBarySad = originIndexBarysSad[i];
        const std::vector<size_t> &indexBaryMax = originIndexBarysMax[i];
        const std::vector<size_t> &indexDataMin = originIndexDatasMin[i];
        const std::vector<size_t> &indexDataSad = originIndexDatasSad[i];
        const std::vector<size_t> &indexDataMax = originIndexDatasMax[i];
        std::vector<double> &weights = vectorWeights[i];
        computeGradientWeights(
          gradWeights, hessianList, dictDiagrams, matchingsAtoms, Barycenter,
          Data, matchingsMin, matchingsMax, matchingsSad, indexBaryMin,
          indexBaryMax, indexBarySad, indexDataMin, indexDataMax, indexDataSad,
          do_optimizeAtoms);
        int nb_points = Barycenter.size();
        gradActor.executeWeightsProjected(
          hessianList, weights, gradWeights, epoch, nb_points, MaxEigenValue_);
        }

      this->printMsg("Computed 1st opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt1.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    // std::vector<double> &weight = vectorWeights[0];

    for(size_t p = 0; p < nDiags; ++p) {
      allLossesAtEpoch[p] = 0.;
      trueAllLossesAtEpoch[p] = 0.;
    }
    barycentersList.clear();
    barycentersList.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);
    ////////////////////////////////ATOM////////////////////////////////////////

    // this->printMsg(
    // "========================ATOM NOW=============================");

    if(do_optimizeAtoms) {
      Timer tm_it2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &barycenter = barycentersList[i];
        std::vector<double> &weight = vectorWeights[i];
        // double sum = 0.;
        // for(int q = 0; q < weight.size(); ++q) {
        //   sum += weight[q];
        //   std::cout << weight[q] << std::endl;
        // }
        // std::cout << "sum: " << sum << std::endl;
        // this->printMsg(std::to_string(sum_temp));
        std::vector<std::vector<ttk::MatchingType>> &matchings
          = allMatchingsAtoms[i];
        computeWeightedBarycenter(
          dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter_);
        // for(int i = 0; i < barycenter.size(); ++i) {
        // ttk::PersistencePair &t = barycenter[i];
        // std::cout << "Pair: " << t.birth.sfValue << ", " << t.death.sfValue
        //          << std::endl;
        // }
      }
      this->printMsg(
        "Computed 2nd Barycenters for epoch " + std::to_string(epoch),
        epoch / static_cast<double>(MAX_EPOCH), tm_it2.getElapsedTime(),
        threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);
      tm_part += static_cast<double>(tm_it2.getElapsedTime());

      barycentersListMin.clear();
      barycentersListSad.clear();
      barycentersListMax.clear();

      barycentersListMin.resize(nDiags);
      barycentersListSad.resize(nDiags);
      barycentersListMax.resize(nDiags);

      bidderBarycentersListMin.clear();
      bidderBarycentersListSad.clear();
      bidderBarycentersListMax.clear();

      originIndexBarysMin.clear();
      originIndexBarysSad.clear();
      originIndexBarysMax.clear();

      originIndexBarysMin.resize(nDiags);
      originIndexBarysSad.resize(nDiags);
      originIndexBarysMax.resize(nDiags);

      matchingsDatasMin.clear();
      matchingsDatasSad.clear();
      matchingsDatasMax.clear();

      matchingsDatasMin.resize(nDiags);
      matchingsDatasSad.resize(nDiags);
      matchingsDatasMax.resize(nDiags);
      // std::vector<BidderDiagram> bidderBarycentersListMin{};
      // std::vector<BidderDiagram> bidderBarycentersListSad{};
      // std::vector<BidderDiagram> bidderBarycentersListMax{};

      computeAllDistances(barycentersList, nDiags, barycentersListMin, barycentersListSad, barycentersListMax,
        bidderBarycentersListMin, bidderBarycentersListSad, bidderBarycentersListMax,
        originIndexBarysMin, originIndexBarysSad, originIndexBarysMax,
        bidderDiagramsMin, bidderDiagramsMax, bidderDiagramsSad,
        matchingsDatasMin, matchingsDatasMax, matchingsDatasSad, trueBidderDiagramMin,
        trueBidderDiagramSad, trueBidderDiagramMax, allLossesAtEpoch,
        trueAllLossesAtEpoch, false);

      // if(do_optimizeAtoms) {
      std::vector<std::vector<std::vector<std::array<double, 2>>>>
        allPairToAddToGradList(nDiags);
      std::vector<ttk::DiagramType> allInfoToAdd(nDiags);
      std::vector<std::vector<Matrix>> gradsAtomsList(nDiags);
      std::vector<std::vector<int>> checkerAtomsList(nDiags);
      // for(size_t i = 0 ; i < nDiags ; ++i){
      //   auto &checkerAtoms = checkerAtomsList[i];
      //   checkerAtoms.resize()
      // }

      bool do_DimReduct = false;
      if(DimReductMode_ && numAtom <= 3) {
        do_DimReduct = true;
      } else {
        do_DimReduct = false;
      }

      Timer tm_opt2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &gradsAtoms = gradsAtomsList[i];
        auto &checkerAtoms = checkerAtomsList[i];
        // const auto &matchingsAtoms = allMatchingsAtoms[i];
        const auto &Barycenter = barycentersList[i];
        const auto &Data = intermediateDiagrams[i];
        // std::vector<Matrix> &gradsAtoms = gradsAtomsList[i];
        const std::vector<ttk::MatchingType> &matchingsMin
          = matchingsDatasMin[i];
        const std::vector<ttk::MatchingType> &matchingsMax
          = matchingsDatasMax[i];
        const std::vector<ttk::MatchingType> &matchingsSad
          = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = originIndexBarysMin[i];
        const std::vector<size_t> &indexBarySad = originIndexBarysSad[i];
        const std::vector<size_t> &indexBaryMax = originIndexBarysMax[i];
        const std::vector<size_t> &indexDataMin = originIndexDatasMin[i];
        const std::vector<size_t> &indexDataSad = originIndexDatasSad[i];
        const std::vector<size_t> &indexDataMax = originIndexDatasMax[i];
        const std::vector<double> &weights = vectorWeights[i];
        // int nb_points = barycentersList[i].size();
        // std::vector<int> checkerAtoms(Barycenter.size(), 0);
        // std::cout << "DIAG: " << i << " =====================" << "\n";
        computeGradientAtoms(
          gradsAtoms, weights, Barycenter, Data, matchingsMin, matchingsMax,
          matchingsSad, indexBaryMin, indexBaryMax, indexBarySad, indexDataMin,
          indexDataMax, indexDataSad, checkerAtoms, pairToAddGradList,
          infoToAdd, static_cast<int>(nDiags), do_DimReduct);
        // std::cout << "=========================================" << "\n";
        //  gradActor.executeAtoms(dictDiagrams, matchingsAtoms, Barycenter,
        //                         gradsAtoms, nb_points, checkerAtoms, epoch);
      }

      std::vector<double> maxiDeath(numAtom);
      std::vector<double> minBirth(numAtom);
      for(int j = 0; j < numAtom; ++j) {
        auto &atom = dictDiagrams[j];
        auto &temp = atom[0];
        maxiDeath[j] = temp.death.sfValue;
        minBirth[j] = temp.birth.sfValue;
      }

      std::vector<std::vector<std::vector<int>>> allProjectionsList(nDiags);
      std::vector<ttk::DiagramType> allFeaturesToAdd(nDiags);
      std::vector<std::vector<std::array<double, 2>>> allProjLocations(nDiags);
      std::vector<std::vector<std::vector<double>>>
        allVectorForProjContributions(nDiags);
      std::vector<ttk::DiagramType> allTrueFeaturesToAdd(numAtom);
      std::vector<std::vector<std::vector<int>>> allTrueProj(numAtom);
      std::vector<std::vector<std::array<double, 2>>> allTrueProjLoc(numAtom);
      // std::vector<std::vector<int>> allAtomIndices(nDiags);

      for(size_t i = 0; i < nDiags; ++i) {
        // std::cout << " OPTIM DIAG: " << i << std::endl;
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &projForDiag = allProjectionsList[i];
        auto &vectorForProjContrib = allVectorForProjContributions[i];
        auto &featuresToAdd = allFeaturesToAdd[i];
        auto &projLocations = allProjLocations[i];
        auto &gradsAtoms = gradsAtomsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const auto &Barycenter = barycentersList[i];
        const auto &checkerAtoms = checkerAtomsList[i];
        int nb_points = barycentersList[i].size();
        gradActor.executeAtoms(
          dictDiagrams, matchingsAtoms, Barycenter, gradsAtoms, nb_points,
          checkerAtoms, epoch, projForDiag, featuresToAdd, projLocations,
          vectorForProjContrib, pairToAddGradList, infoToAdd);
        // std::cout << " OPTIM DIAG: " << i << std::endl;
        // std::cout << "=================================================="
        // << std::endl;
      }

      if(CreationFeatures_) {

        // std::cout << "CREATING FEATURES" << std::endl;
        // double factEquiv = sqrt(static_cast<double>(numAtom));

        for(size_t i = 0; i < nDiags; ++i) {
          auto &projForDiag = allProjectionsList[i];
          auto &featuresToAdd = allFeaturesToAdd[i];
          auto &projLocations = allProjLocations[i];
          auto &vectorForProjContrib = allVectorForProjContributions[i];
          for(size_t j = 0; j < projForDiag.size(); ++j) {
            auto &t = featuresToAdd[j];
            std::array<double, 2> &pair = projLocations[j];
            std::vector<double> &vectorContrib = vectorForProjContrib[j];
            std::vector<int> &projAndIndex = projForDiag[j];
            std::vector<int> proj(numAtom);
            for(int m = 0; m < numAtom; ++m) {
              proj[m] = projAndIndex[m];
            }
            int atomIndex = static_cast<int>(projAndIndex[numAtom]);
            // auto it = std::find(allTrueProj[atomIndex].begin() ,
            // allTrueProj[atomIndex].end() , proj); bool ralph = it !=
            // allTrueProj[atomIndex].end();
            bool lenNull = allTrueProj[atomIndex].size() == 0;
            if(lenNull) {
              const CriticalType c1 = t.birth.type;
              const CriticalType c2 = t.death.type;
              const SimplexId idTemp = t.dim;
              pair[0] = pair[0] - step * vectorContrib[0];
              pair[1] = pair[1] - step * vectorContrib[1];
              if(pair[0] > pair[1]) {
                pair[1] = pair[0];
              }
              if(pair[0] < minBirth[atomIndex]) {
                continue;
              }
              if(pair[1] - pair[0] < 1e-7){
                continue;
              }

              PersistencePair newPair{CriticalVertex{0, c1, pair[0], {}},
                                      CriticalVertex{0, c2, pair[1], {}},
                                      idTemp, true};
              allTrueProj[atomIndex].push_back(proj);
              allTrueFeaturesToAdd[atomIndex].push_back(newPair);
              allTrueProjLoc[atomIndex].push_back(pair);
            } else {
              // auto it = std::find(allTrueProj[atomIndex].begin() ,
              // allTrueProj[atomIndex].end() , proj); bool ralph = it ==
              // allTrueProj[atomIndex].end();
              bool ralph = true;
              size_t index = 0;
              if(Fusion_) {
                for(size_t n = 0; n < allTrueProj[atomIndex].size(); ++n) {
                  auto &projStocked = allTrueProj[atomIndex][n];
                  auto &projLocStocked = allTrueProj[atomIndex][n];
                  double distance
                    = sqrt(pow((pair[0] - projLocStocked[0]), 2)
                          + pow((pair[1] - projLocStocked[1]), 2));
                  if(proj == projStocked && distance < 1e-3) {
                    ralph = false;
                    index = n;
                    break;
                  }
                }
              }
              if(ralph) {

                const CriticalType c1 = t.birth.type;
                const CriticalType c2 = t.death.type;
                const SimplexId idTemp = t.dim;
                pair[0] = pair[0] - step * vectorContrib[0];
                pair[1] = pair[1] - step * vectorContrib[1];
                if(pair[0] > pair[1]) {
                  pair[1] = pair[0];
                }
                if(pair[0] < minBirth[atomIndex]) {
                  continue;
                }
                if(pair[1] - pair[0] < 1e-7 ){
                  continue;
                }


                PersistencePair newPair{CriticalVertex{0, c1, pair[0], {}},
                                        CriticalVertex{0, c2, pair[1], {}},
                                        idTemp, true};
                allTrueProj[atomIndex].push_back(proj);
                allTrueFeaturesToAdd[atomIndex].push_back(newPair);
                allTrueProjLoc[atomIndex].push_back(pair);

              } else {
                // auto index = std::distance(allTrueProj[atomIndex].begin() ,
                // it);
                auto &tReal = allTrueFeaturesToAdd[atomIndex][index];
                tReal.birth.sfValue
                  = tReal.birth.sfValue - step * vectorContrib[0];
                tReal.death.sfValue
                  = tReal.death.sfValue - step * vectorContrib[1];
                if(tReal.birth.sfValue > tReal.death.sfValue) {
                  tReal.death.sfValue = tReal.birth.sfValue;
                }
              }
            }
          }
        }

        if(!do_compression){
          // if (CreationFeatures){
          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto initSize = initSizes[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            if(histoEpochAtom.size() > 0) {
              for(size_t j = 0; j < histoEpochAtom.size(); ++j) {
                auto &t = atom[initSize + j];
                histoEpochAtom[j] += 1;
                histoBoolAtom[j] = t.death.sfValue - t.birth.sfValue < 0.1*(percent/100.)*maxiDeath[i];
                boolDiag[j] = t.death.sfValue - t.birth.sfValue < 1e-6;
                boolUnderDiag[j] = t.death.sfValue < t.birth.sfValue;
                boolAboveGlobal[j] = t.birth.sfValue > maxiDeath[i];
              }
            }
          }

          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &trueFeaturesToAdd = allTrueFeaturesToAdd[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            for(size_t j = 0; j < trueFeaturesToAdd.size(); ++j) {
              auto &t = trueFeaturesToAdd[j];
              atom.push_back(t);
              histoEpochAtom.push_back(0);
              histoBoolAtom.push_back(t.death.sfValue - t.birth.sfValue < 0.1*(percent/100.)*maxiDeath[i]);
              boolDiag.push_back(t.death.sfValue - t.birth.sfValue < 1e-6);
              boolUnderDiag.push_back(t.death.sfValue < t.birth.sfValue);
              boolAboveGlobal.push_back(t.birth.sfValue > maxiDeath[i]);
            }
          }

          std::vector<std::vector<size_t>> allIndicesToDelete(numAtom);
          for(int i = 0; i < numAtom; ++i) {
            auto &indicesAtomToDelete = allIndicesToDelete[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            for(size_t j = 0; j < histoEpochAtom.size(); ++j) {
              if(boolUnderDiag[j] || boolDiag[j] || boolAboveGlobal[j]
                || (histoEpochAtom[j] > 5 && histoBoolAtom[j])) {
                // if(boolUnderDiag[j] || (histoEpochAtom[j] > 2 &&
                // histoBoolAtom[j])){
                indicesAtomToDelete.push_back(j);
              }
            }
          }

          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &histoEpochAtom = histoAllEpochLife[i];
            auto &histoBoolAtom = histoAllBoolLife[i];
            auto &indicesAtomToDelete = allIndicesToDelete[i];
            auto &boolUnderDiag = checkUnderDiag[i];
            auto &boolDiag = checkDiag[i];
            auto &boolAboveGlobal = checkAboveGlobal[i];
            auto initSize = initSizes[i];
            // size_t initAtomSize = atom.size() -histoEpochAtom.size();
            if(static_cast<int>(indicesAtomToDelete.size()) > 0) {
              for(int j = static_cast<int>(indicesAtomToDelete.size()) - 1;
                  j >= 0; j--) {
                atom.erase(atom.begin() + initSize + indicesAtomToDelete[j]);
                histoEpochAtom.erase(histoEpochAtom.begin()
                                    + indicesAtomToDelete[j]);
                histoBoolAtom.erase(histoBoolAtom.begin()
                                    + indicesAtomToDelete[j]);
                boolUnderDiag.erase(boolUnderDiag.begin()
                                    + indicesAtomToDelete[j]);
                boolDiag.erase(boolDiag.begin() + indicesAtomToDelete[j]);
                boolAboveGlobal.erase(boolAboveGlobal.begin()
                                      + indicesAtomToDelete[j]);
              }
            } else {
              continue;
            }
          }
        } else {
          for(int i = 0; i < numAtom; ++i) {
            auto &atom = dictDiagrams[i];
            auto &trueFeaturesToAdd = allTrueFeaturesToAdd[i];
            for(size_t j = 0; j < trueFeaturesToAdd.size(); ++j) {
              auto &t = trueFeaturesToAdd[j];
              atom.push_back(t);
            }
          }
          controlAtomsSize(intermediateDiagrams, dictDiagrams);
        }
      } else {
        if(do_compression){
          if(epoch > -1){
              controlAtomsSize(intermediateDiagrams, dictDiagrams);
          }
        }
      }

      for(int i = 0; i < numAtom; ++i) {
        auto &atom = dictDiagrams[i];
        auto &globalPair = atom[0];
        for(size_t j = 0; j < atom.size(); ++j) {
          auto &t = atom[j];

          if(t.death.sfValue > globalPair.death.sfValue) {
            t.death.sfValue = globalPair.death.sfValue;
          }

          if(t.birth.sfValue > t.death.sfValue) {
            t.death.sfValue = t.birth.sfValue;
          }

          if(t.birth.sfValue < globalPair.birth.sfValue) {
            t.birth.sfValue = globalPair.birth.sfValue;
          }


        }
      }
      this->printMsg("Computed 2nd opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt2.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);

      // ATOM OPTIMIZATION
    }
    epoch += 1;
    // this->printMsg("=====================================================");
    barycentersList.clear();
    barycentersList.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);

  }
  printMsg(
    " Epoch " + std::to_string(epoch) + ", loss = " + std::to_string(loss), 1,
    threadNumber_);

  printMsg("Loss returned "
           + std::to_string(*std::min_element(
             loss_tab.begin() + nbEpochPrevious, loss_tab.end()))
           + " at Epoch "
           + std::to_string(std::min_element(loss_tab.begin() + nbEpochPrevious,
                                             loss_tab.end())
                            - loss_tab.begin()));

  for(size_t p = 0; p < dictDiagrams.size(); ++p) {
    auto atom = histoDictDiagrams[p];
    dictDiagrams[p] = atom;
  }
  for(size_t p = 0; p < nDiags; ++p) {
    auto weights = histoVectorWeights[p];
    vectorWeights[p] = weights;
  }

  // this->printMsg("time spent computing barycenter" +
  // std::to_string(tm_part));
  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

}

double PersistenceDiagramDictEncoding::distVect(
  const std::vector<double> &vec1, const std::vector<double> &vec2) const {

  double dist = 0.;
  for(size_t i = 0; i < vec1.size(); ++i) {
    dist = dist + (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
  }
  return std::sqrt(dist);
}

double PersistenceDiagramDictEncoding::getMostPersistent(
  const std::vector<BidderDiagram> &bidder_diags) const {

  double max_persistence = 0;

  for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
    for(size_t j = 0; j < bidder_diags[i].size(); ++j) {
      const double persistence = bidder_diags[i].at(j).getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }

  return max_persistence;
}

double PersistenceDiagramDictEncoding::computeDistance(
  const BidderDiagram &D1,
  const BidderDiagram &D2,
  std::vector<ttk::MatchingType> &matching) const {

  GoodDiagram D2_bis{};
  for(size_t i = 0; i < D2.size(); i++) {
    const Bidder &b = D2.at(i);
    Good g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_);
    g.setPrice(0);
    D2_bis.emplace_back(g);
  }

  PersistenceDiagramAuction auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(D1, D2_bis);
  double loss = auction.run(matching);
  return loss;
}

void PersistenceDiagramDictEncoding::computeGradientWeights(
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
  const bool do_optimizeAtoms) const {

  // initialization
  std::vector<std::vector<std::array<double, 2>>> grad_list(Barycenter.size());
  std::vector<std::vector<std::array<double, 2>>> pairToAddGradList;
  for(size_t i = 0; i < grad_list.size(); ++i) {
    grad_list[i].resize(matchingsAtoms.size());
  }
  // std::vector<ttk::MatchingType> matching;
  std::vector<std::array<double, 2>> directions(Barycenter.size());
  std::vector<std::array<double, 2>> data_assigned(Barycenter.size());
  // std::vector<double> gradient(dictDiagrams.size(), 0.);
  gradWeights.resize(dictDiagrams.size());
  for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    gradWeights[i] = 0.;
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  for(size_t j = 0; j < Barycenter.size(); ++j) {
    checker[j].resize(matchingsAtoms.size());
  }
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<int> tracker2(Barycenter.size(), 0);

  // this->printMsg("error?2");
  // computing gradients
  for(size_t i = 0; i < matchingsAtoms.size(); ++i) {
    for(size_t j = 0; j < matchingsAtoms[i].size(); ++j) {
      const ttk::MatchingType &t = matchingsAtoms[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      // if(Id2 < 0) {
      if(Id2 < 0 || static_cast<int>(grad_list.size()) <= Id2
         || static_cast<int>(dictDiagrams[i].size()) <= Id1) {
        continue;
      } else if(Id1 < 0) {
        // this->printMsg("========DIAGONAL=========");
        const PersistencePair &t3 = Barycenter[Id2];
        auto &point = grad_list[Id2][i];
        const double birth_barycenter = t3.birth.sfValue;
        const double death_barycenter = t3.death.sfValue;
        // std::cout << "Barycenter Pair:" << birth_barycenter << " "
        //           << death_barycenter << std::endl;
        const double birth_death_atom
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        point[0] = birth_death_atom;
        point[1] = birth_death_atom;
        // std::cout << "Proj coordinates" << birth_death_atom << std::endl;
        // checker[Id2].push_back(i);
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      } else {
        // this->printMsg("====UPDATE GRADLIST========");
        const PersistencePair &t2 = dictDiagrams[i][Id1];
        auto &point = grad_list[Id2][i];
        const double birth_atom = t2.birth.sfValue;
        const double death_atom = t2.death.sfValue;
        point[0] = birth_atom;
        point[1] = death_atom;
        // checker[Id2].push_back(i);
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      }
    }
  }

  // this->printMsg("===================PASSED============================");
  // computeDistance(newDataBidder, barycenterBidder, matching);

  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsMin, indexBaryMin, 
    indexDataMin, pairToAddGradList, directions, data_assigned, tracker2, do_optimizeAtoms);


  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsMax, indexBaryMax, 
    indexDataMax, pairToAddGradList, directions, data_assigned, tracker2, do_optimizeAtoms);


  computeDirectionsGradWeight(matchingsAtoms, Barycenter, newData, matchingsSad, indexBarySad, 
    indexDataSad, pairToAddGradList, directions, data_assigned, tracker2, do_optimizeAtoms);

  std::vector<int> temp(pairToAddGradList.size(), 1);
  grad_list.insert(
    grad_list.end(), pairToAddGradList.begin(), pairToAddGradList.end());
  tracker2.insert(tracker2.end(), temp.begin(), temp.end());
  tracker.insert(tracker.end(), temp.begin(), temp.end());
  std::vector<int> temp2(matchingsAtoms.size());
  for(size_t j = 0; j < matchingsAtoms.size(); ++j) {
    temp2.push_back(static_cast<int>(j));
  }
  for(size_t j = 0; j < pairToAddGradList.size(); ++j) {
    checker.push_back(temp2);
  }

  // this->printMsg("======================PASSED2==========================");
  for(size_t i = 0; i < grad_list.size(); ++i) {
    const auto &data_point = data_assigned[i];
    for(size_t j = 0; j < checker[i].size(); ++j) {
      auto &point = grad_list[i][checker[i][j]];
      point[0] -= data_point[0];
      point[1] -= data_point[1];
    }
  }

  for(size_t i = 0; i < grad_list.size(); ++i) {
    // for(auto it = std::begin(checker); it != std::end(checker); ++it) {
    if(tracker[i] == 0 || tracker2[i] == 0) {
      // this->printMsg("Not tracked");
      continue;
    } else {
      for(size_t j = 0; j < checker[i].size(); ++j) {
        const auto &point = grad_list[i][checker[i][j]];
        // const double birth = t.birth.sfValue;
        // const double death = t.death.sfValue;
        const auto &direction = directions[i];
        // gradient[j] += -2 * (birth * direction[0] + death * direction[1]);
        gradWeights[checker[i][j]]
          += -2 * (point[0] * direction[0] + point[1] * direction[1]);
      }
    }
  }
  // this->printMsg("======================PASSED2==========================");
  hessianList.resize(grad_list.size());
  for(size_t i = 0; i < grad_list.size(); ++i) {
    Matrix &hessian = hessianList[i];
    hessian.resize(checker[i].size());
    for(size_t j = 0; j < checker[i].size(); ++j) {
      auto &line = hessian[j];
      // hessian[j].resize(checker[i].size());
      line.resize(checker[i].size());
      const auto &point = grad_list[i][checker[i][j]];
      for(size_t q = 0; q < checker[i].size(); ++q) {
        const auto &point_temp = grad_list[i][checker[i][q]];
        line[q] = point[0] * point_temp[0] + point[1] * point_temp[1];
        // this->printMsg("======================COMPUTING==========================");
      }
    }
  }
  // this->printMsg("======================PASSED3==========================");

  // return gradient;
}

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

void PersistenceDiagramDictEncoding::computeGradientAtoms(
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
  std::vector<PersistencePair> &infoToAdd,
  int nbDiags,
  bool do_DimReduct) const {

  // std::vector<ttk::MatchingType> matching;
  gradsAtoms.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    gradsAtoms[i].resize(weights.size());
  }

  checker.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    checker[i] = 0;
  }
  std::vector<std::vector<double>> directions(Barycenter.size());

  computeDirectionsGradAtoms( gradsAtoms, Barycenter, weights, newData,
  matchingsMin, indexBaryMin, indexDataMin, 
  pairToAddGradList, directions, checker, 
  infoToAdd, do_DimReduct);

  computeDirectionsGradAtoms( gradsAtoms, Barycenter, weights, newData,
  matchingsMax, indexBaryMax, indexDataMax, 
  pairToAddGradList, directions, checker, 
  infoToAdd, do_DimReduct);

computeDirectionsGradAtoms( gradsAtoms, Barycenter, weights, newData,
  matchingsSad, indexBarySad, indexDataSad, 
  pairToAddGradList, directions, checker, 
  infoToAdd, do_DimReduct);


  for(size_t i = 0; i < Barycenter.size(); ++i) {
    if(checker[i] == 0) {
      continue;
    } else {
      for(size_t j = 0; j < weights.size(); ++j) {
        std::vector<double> temp(2);
        const std::vector<double> &direction = directions[i];
        temp[0] = -2 * weights[j] * direction[0];
        temp[1] = -2 * weights[j] * direction[1];
        gradsAtoms[i][j] = temp;
      }
    }
  }
}


// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

void PersistenceDiagramDictEncoding::setBidderDiagrams(
  const size_t nInputs,
  std::vector<ttk::DiagramType> &inputDiagrams,
  std::vector<BidderDiagram> &bidder_diags) const {

  bidder_diags.resize(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP

  for(size_t i = 0; i < nInputs; i++) {
    auto &diag = inputDiagrams[i];
    auto &bidders = bidder_diags[i];

    for(size_t j = 0; j < diag.size(); j++) {
      // Add bidder to bidders
      Bidder b(diag[j], j, this->Lambda);
      b.setPositionInAuction(bidders.size());
      bidders.emplace_back(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printMsg("Diagonal point in diagram " + std::to_string(i) + "!",
                       ttk::debug::Priority::DETAIL);
      }
    }
  }
}

void PersistenceDiagramDictEncoding::enrichCurrentBidderDiagrams(
  const std::vector<BidderDiagram> &bidder_diags,
  std::vector<BidderDiagram> &current_bidder_diags,
  const std::vector<double> &maxDiagPersistence) const {

  current_bidder_diags.resize(bidder_diags.size());
  const auto nInputs = current_bidder_diags.size();
  const auto maxPersistence
    = *std::max_element(maxDiagPersistence.begin(), maxDiagPersistence.end());

  if(this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL) {
    for(size_t i = 0; i < nInputs; ++i) {
      for(size_t j = 0; j < bidder_diags[i].size(); ++j) {
        auto b = bidder_diags[i].at(j);

        if( // filter out pairs below absolute persistence threshold
          (this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
           && b.getPersistence() > this->MinPersistence_)
          || // filter out pairs below persistence threshold relative to
          // the most persistent pair *of each diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
           && b.getPersistence()
                > this->MinPersistence_ * maxDiagPersistence[i])
          || // filter out pairs below persistence threshold relative to the
             // most persistence pair *in all diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL
           && b.getPersistence() > this->MinPersistence_ * maxPersistence)) {
          b.id_ = current_bidder_diags[i].size();
          b.setPositionInAuction(current_bidder_diags[i].size());
          current_bidder_diags[i].emplace_back(b);
        }
      }
    }
    return;
  }

  const double prev_min_persistence = 2.0 * getMostPersistent(bidder_diags);
  double new_min_persistence = 0.0;

  // 1. Get size of the largest current diagram, deduce the maximal number
  // of points to append
  size_t max_diagram_size = 0;
  for(const auto &diag : current_bidder_diags) {
    max_diagram_size
      = std::max(static_cast<size_t>(diag.size()), max_diagram_size);
  }
  size_t max_points_to_add = std::max(
    this->MaxNumberOfPairs, this->MaxNumberOfPairs + max_diagram_size / 10);
  // 2. Get which points can be added, deduce the new minimal persistence
  std::vector<std::vector<int>> candidates_to_be_added(nInputs);
  std::vector<std::vector<size_t>> idx(nInputs);

  for(size_t i = 0; i < nInputs; i++) {
    double local_min_persistence = std::numeric_limits<double>::min();
    std::vector<double> persistences;
    for(size_t j = 0; j < bidder_diags[i].size(); j++) {
      Bidder b = bidder_diags[i].at(j);
      double persistence = b.getPersistence();
      if(persistence >= 0.0 && persistence <= prev_min_persistence) {
        candidates_to_be_added[i].emplace_back(j);
        idx[i].emplace_back(idx[i].size());
        persistences.emplace_back(persistence);
      }
    }
    const auto cmp = [&persistences](const size_t a, const size_t b) {
      return ((persistences[a] > persistences[b])
              || ((persistences[a] == persistences[b]) && (a > b)));
    };
    std::sort(idx[i].begin(), idx[i].end(), cmp);
    const auto size = candidates_to_be_added[i].size();
    if(size >= max_points_to_add) {
      double last_persistence_added
        = persistences[idx[i][max_points_to_add - 1]];
      if(last_persistence_added > local_min_persistence) {
        local_min_persistence = last_persistence_added;
      }
    }
    if(i == 0) {
      new_min_persistence = local_min_persistence;
    } else {
      if(local_min_persistence < new_min_persistence) {
        new_min_persistence = local_min_persistence;
      }
    }
    // 3. Add the points to the current diagrams
    const auto s = candidates_to_be_added[i].size();
    for(size_t j = 0; j < std::min(max_points_to_add, s); j++) {
      Bidder b = bidder_diags[i].at(candidates_to_be_added[i][idx[i][j]]);
      const double persistence = b.getPersistence();
      if(persistence >= new_min_persistence) {
        b.id_ = current_bidder_diags[i].size();
        b.setPositionInAuction(current_bidder_diags[i].size());
        current_bidder_diags[i].emplace_back(b);
      }
    }
  }
}

int PersistenceDiagramDictEncoding::initDictionary(
  std::vector<ttk::DiagramType> &dictDiagrams,
  const std::vector<ttk::DiagramType> &datas,
  const std::vector<ttk::DiagramType> &inputAtoms,
  const int nbAtom,
  bool do_min,
  bool do_sad,
  bool do_max,
  int seed,
  double percent) {
  switch(this->BackEnd) {
    case BACKEND::INPUT_ATOMS: {
      if(static_cast<int>(inputAtoms.size()) != nbAtom) {
        printErr("number of atoms don't match, going with "
                 + std::to_string(inputAtoms.size()) + " atoms");
      }
      for(size_t i = 0; i < inputAtoms.size(); i++) {
        const auto &t = inputAtoms[i];
        double maxPers = getMaxPers(t);

        dictDiagrams.push_back(t);
        dictDiagrams[i].erase(
          std::remove_if(dictDiagrams[i].begin(), dictDiagrams[i].end(),
                         [maxPers](ttk::PersistencePair &p) {
                           return (p.death.sfValue - p.birth.sfValue)
                                  < 0. * maxPers;
                         }),
          dictDiagrams[i].end());
      }
      break;
    }

    case BACKEND::BORDER_INIT: {
      ttk::InitFarBorderDict initializer;
      initializer.setThreadNumber(this->threadNumber_);
      initializer.execute(dictDiagrams, datas, nbAtom, do_min, do_sad, do_max);
      break;
    }

    case BACKEND::RANDOM_INIT:
      ttk::InitRandomDict{}.execute(dictDiagrams, datas, nbAtom, seed); // TODO
      break;

    case BACKEND::FIRST_DIAGS: {
      for(int i = 0; i < nbAtom; ++i) {
        const auto &t = datas[i];
        dictDiagrams.push_back(t);
      }
      break;
    }
    case BACKEND::GREEDY_INIT: {
      dictDiagrams.resize(datas.size());
      for(size_t i = 0; i < datas.size(); ++i) {
        dictDiagrams[i] = datas[i];
      }
      // const std::array<size_t, 2> nInputsUseless{100, 100};
      while(static_cast<int>(dictDiagrams.size()) != nbAtom) {
        // std::vector<ttk::DiagramType> dictTemp;
        std::vector<double> allEnergy(dictDiagrams.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t j = 0; j < dictDiagrams.size(); ++j) {
          std::vector<ttk::DiagramType> dictTemp;
          std::vector<ttk::DiagramType> dataAlone;
          std::vector<double> lossTabTemp;
          std::vector<double> timersTemp;
          std::vector<double> trueLossTabTemp;
          std::vector<std::vector<double>> allLossesTemp(1);
          Timer tm_temp{};

          for(size_t p = 0; p < dictDiagrams.size(); ++p) {
            if(p != j) {
              dictTemp.push_back(dictDiagrams[p]);
            } else {
              dataAlone.push_back(dictDiagrams[p]);
            }
          }
          std::vector<std::vector<double>> weightsTemp(1);
          std::vector<double> weights(
            dictTemp.size(), 1. / (static_cast<int>(dictTemp.size()) * 1.));
          weightsTemp[0] = weights;
          std::vector<std::vector<double>> histoVectorWeights(1);
          std::vector<ttk::DiagramType> histoDictDiagrams(dictTemp.size());
          std::vector<BidderDiagram> bidderTempMin(dataAlone.size());
          std::vector<BidderDiagram> bidderTempMax(dataAlone.size());
          std::vector<BidderDiagram> bidderTempSad(dataAlone.size());
          bool do_compression = false;
          this->method(dataAlone, dictTemp, weightsTemp, seed,
                       static_cast<int>(dictTemp.size()), lossTabTemp,
                       trueLossTabTemp, timersTemp, allLossesTemp,
                       histoVectorWeights, histoDictDiagrams, false, 0.01,
                       bidderTempMin, bidderTempSad, bidderTempMax, tm_temp,
                       percent, do_compression);
          double min_loss
            = *std::min_element(lossTabTemp.begin(), lossTabTemp.end());
          allEnergy[j] = min_loss;
        }
        int index = std::min_element(allEnergy.begin(), allEnergy.end())
                    - allEnergy.begin();
        dictDiagrams.erase(dictDiagrams.begin() + index);
      }
      break;
    }
    default:
      break;
  }
  return 0;
}

void PersistenceDiagramDictEncoding::gettingBidderDiagrams(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<BidderDiagram> &bidderDiagramsMin,
  std::vector<BidderDiagram> &bidderDiagramsSad,
  std::vector<BidderDiagram> &bidderDiagramsMax) {

  size_t nDiags = intermediateDiagrams.size();
  // double distance = 0.;

  std::vector<ttk::DiagramType> inputDiagramsMin(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsSad(nDiags);
  std::vector<ttk::DiagramType> inputDiagramsMax(nDiags);

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; i++) {
    const auto &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const PersistencePair &t = CTDiagram[j];
      const ttk::CriticalType nt1 = t.birth.type;
      const ttk::CriticalType nt2 = t.death.type;
      const double pers = t.persistence();
      // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMin[i].emplace_back(t);
          // originIndexDatasMax[i].push_back(j);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            // originIndexDatasMax[i].push_back(j);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
            // originIndexDatasMin[i].push_back(j);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
            // originIndexDatasSad[i].push_back(j);
          }
        }
      }
    }
  }

  if(this->do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidderDiagramsMin);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidderDiagramsSad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidderDiagramsMax);
  }

  // return distance;
}

double PersistenceDiagramDictEncoding::getMaxPers(const ttk::DiagramType &data) {
  double maxPers = 0.;
  for(size_t j = 0; j < data.size(); ++j) {
    auto &t = data[j];
    double pers = t.death.sfValue - t.birth.sfValue;
    if(pers > maxPers) {
      maxPers = pers;
    }
  }

  return maxPers;
}

void PersistenceDiagramDictEncoding::controlAtomsSize(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<ttk::DiagramType> &dictDiagrams) const{

  size_t m = dictDiagrams.size();
  int globalSize = 0;

  for(size_t j = 0; j < intermediateDiagrams.size() ; ++j){
    auto &data = intermediateDiagrams[j];
    globalSize += static_cast<int>(data.size());
  }

  int dictSize = 0;

  for(size_t j = 0; j < m ; ++j){
    auto &atom = dictDiagrams[j];
    dictSize += static_cast<int>(atom.size());
  }

  if(static_cast<double>(dictSize) > (1./this->CompressionFactor)*static_cast<double>(globalSize)){
    double factor = (1./this->CompressionFactor)*(static_cast<double>(globalSize)/static_cast<double>(dictSize));
    std::vector<std::vector<double>> tempDictPersistencePairs(m);

    for(size_t j = 0; j < m ; ++j){
      auto &temp = tempDictPersistencePairs[j];
      auto &atom = dictDiagrams[j];
      temp.resize(atom.size());
      for(size_t p = 0; p < atom.size(); ++p){
        auto &t = atom[p];
        temp[p] = t.persistence();
      }
    }

    for(size_t j = 0; j < m; ++j) {
      auto &temp = tempDictPersistencePairs[j];
      std::sort(temp.begin(), temp.end(), std::greater<double>());
    }

    std::vector<double> persThresholds(m);
    for(size_t j = 0; j < m ; ++j){
      auto &temp = tempDictPersistencePairs[j];
      int index = static_cast<int>(floor(factor*static_cast<double>(temp.size())));
      persThresholds[j] = temp[index];
    }

    for(size_t j = 0; j < m ; ++j){
      double persThreshold = persThresholds[j];
      auto &atom = dictDiagrams[j];
      atom.erase(
        std::remove_if(atom.begin(),
                      atom.end(),
                      [persThreshold](ttk::PersistencePair &t) {
                        return (t.death.sfValue - t.birth.sfValue)
                                < persThreshold;
                      }),
        atom.end());
    }
  }
}

void PersistenceDiagramDictEncoding::controlAtomsSize2(
  const std::vector<ttk::DiagramType> &intermediateDiagrams,
  std::vector<ttk::DiagramType> &dictDiagrams) const{
  size_t m = dictDiagrams.size();
  int globalSize = 0;

  for(size_t j = 0; j < intermediateDiagrams.size() ; ++j){
    auto &data = intermediateDiagrams[j];
    globalSize += static_cast<int>(data.size());
  }

  int dictSize = 0;

  for(size_t j = 0; j < m ; ++j){
    auto &atom = dictDiagrams[j];
    dictSize += static_cast<int>(atom.size());
  }

  if(static_cast<double>(dictSize) > (1./this->CompressionFactor)*static_cast<double>(globalSize)){
    double factor = (1./this->CompressionFactor)*static_cast<double>(globalSize)/static_cast<double>(dictSize);
    std::vector<double> tempDictPersistencePairs;
    for(size_t j = 0; j < m ; ++j){
      auto &atom = dictDiagrams[j];
      for(size_t p = 0; p < atom.size(); ++p){
        auto &t = atom[p];
        tempDictPersistencePairs.emplace_back(t.persistence());
      }
    }

    std::sort(tempDictPersistencePairs.begin(), tempDictPersistencePairs.end(), std::greater<double>());
    int index = static_cast<int>(floor(factor*static_cast<double>(tempDictPersistencePairs.size())));
    double persThreshold = tempDictPersistencePairs[index];

    for(size_t j = 0; j < m ; ++j){
      auto &atom = dictDiagrams[j];
      atom.erase(
        std::remove_if(atom.begin(),
                      atom.end(),
                      [persThreshold](ttk::PersistencePair &t) {
                        return (t.death.sfValue - t.birth.sfValue)
                                < persThreshold;
                      }),
        atom.end());
    }
  }

}

void PersistenceDiagramDictEncoding::computeDirectionsGradWeight(
  const std::vector<std::vector<ttk::MatchingType>> &matchingsAtoms,
  const ttk::DiagramType &Barycenter,
  const ttk::DiagramType &newData,
  const std::vector<ttk::MatchingType> &matchingsCritType,
  const std::vector<size_t> &indexBaryCritType,
  const std::vector<size_t> &indexDataCritType,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  std::vector<std::array<double, 2>> &directions,
  std::vector<std::array<double, 2>> &data_assigned,
  std::vector<int> &tracker2,
  const bool do_optimizeAtoms) const {

  size_t m = matchingsCritType.size();
  int k = 0;
  for(size_t i = 0; i < m; ++i) {
    const ttk::MatchingType &t = matchingsCritType[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    if(Id2 < 0) {
      k += 1;

      if(Id1 < 0) {
        continue;
      } else {
        // if(false){
        if(do_optimizeAtoms && CreationFeatures_ && ProgApproach_) {
          const PersistencePair &t2 = newData[indexDataCritType[Id1]];
          const double birth_data = t2.birth.sfValue;
          const double death_data = t2.death.sfValue;
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::array<double, 2> direction;
          direction[0] = birth_data - birth_death_barycenter;
          direction[1] = death_data - birth_death_barycenter;
          /* std::vector<std::vector<double>> temp3(weights.size()); */
          /* std::vector<double> temp2(2); */
          /* for(size_t j = 0; j < weights.size(); ++j) { */
          /*   temp2[0] += -2 * weights[j] * direction[0]; */
          /*   temp2[1] += -2 * weights[j] * direction[1]; */
          /*   temp3[j] = temp2; */
          /* } */

          std::vector<std::array<double, 2>> newPairs(matchingsAtoms.size());
          for(size_t j = 0; j < matchingsAtoms.size(); ++j) {
            std::array<double, 2> pair{
              birth_death_barycenter, birth_death_barycenter};
            newPairs[j] = pair;
          }
          pairToAddGradList.push_back(newPairs);
          data_assigned.push_back({birth_data, death_data});
          directions.push_back(direction);
        } else {
          continue;
        }
      }
      // this->printMsg("k = " + std::to_string(k));
    } else {
      // this->printMsg("==Here?==");
      // this->printMsg(std::to_string(static_cast<int>(indexBaryCritType[Id2])));
      const PersistencePair &t3 = Barycenter[indexBaryCritType[Id2]];
      const double birth_barycenter = t3.birth.sfValue;
      const double death_barycenter = t3.death.sfValue;
      auto &direction = directions[indexBaryCritType[Id2]];
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
        data_assigned[indexBaryCritType[Id2]] = {birth_death_data, birth_death_data};

      } else {
        // checker[indexBaryCritType[Id2]].push_back(indexDataCritType[Id1]);
        const PersistencePair &t2 = newData[indexDataCritType[Id1]];
        const double birth_data = t2.birth.sfValue;
        const double death_data = t2.death.sfValue;
        // direction[0] = t2.birth.sfValue - t3.birth.sfValue;
        // direction[1] = t2.death.sfValue - t3.death.sfValue;
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        data_assigned[indexBaryCritType[Id2]] = {birth_data, death_data};
        // directions[Id2].push_back(direction);
      }
      tracker2[indexBaryCritType[Id2]] = 1;
    }
  }
}


void PersistenceDiagramDictEncoding::computeDirectionsGradAtoms(
  std::vector<Matrix> &gradsAtoms,
  const ttk::DiagramType &Barycenter,
  const std::vector<double> &weights,
  const ttk::DiagramType &newData,
  const std::vector<ttk::MatchingType> &matchingsCritType,
  const std::vector<size_t> &indexBaryCritType,
  const std::vector<size_t> &indexDataCritType,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  std::vector<std::vector<double>> &directions,
  std::vector<int> &checker,
  std::vector<PersistencePair> &infoToAdd,
  const bool do_DimReduct) const{


  for(size_t i = 0; i < matchingsCritType.size(); ++i) {
    const ttk::MatchingType &t = matchingsCritType[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    // std::cout << Id1 << ""

    int maxWeights = std::max_element(weights.begin(), weights.end()) - weights.begin();
    int k = 0;
    if(Id2 < 0) {
      k += 1;
      if(Id1 < 0) {
        continue;
      } else {
        if(CreationFeatures_) {
          const PersistencePair &t2 = newData[indexDataCritType[Id1]];
          const double birth_data = t2.birth.sfValue;
          const double death_data = t2.death.sfValue;
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::vector<double> direction(2);
          direction[0] = birth_data - birth_death_barycenter;
          direction[1] = death_data - birth_death_barycenter;
          std::vector<std::vector<double>> temp3(weights.size());
          std::vector<double> temp2(2);
          if(CompressionMode_ && !do_DimReduct) {
            for(size_t j = 0; j < weights.size(); ++j) {
              if(j == static_cast<size_t>(maxWeights)){
                temp2[0] = -2 * weights[j] * direction[0];
                temp2[1] = -2 * weights[j] * direction[1];
              } else {
                temp2[0] = 0.;
                temp2[1] = 0.;
              }
              temp3[j] = temp2;
            }
          } else {
            for(size_t j = 0; j < weights.size(); ++j) {
              temp2[0] = -2 * weights[j] * direction[0];
              temp2[1] = -2 * weights[j] * direction[1];
              temp3[j] = temp2;
            }
          }
          gradsAtoms.push_back(temp3);
          checker.push_back(1);
          std::vector<std::array<double, 2>> newPairs(weights.size());
          for(size_t j = 0; j < weights.size(); ++j) {
            std::array<double, 2> pair{
              birth_death_barycenter, birth_death_barycenter};
            newPairs[j] = pair;
          }
          pairToAddGradList.push_back(newPairs);
          infoToAdd.push_back(t2);

        } else {
          continue;
        }
      }
      // this->printMsg("k = " + std::to_string(k));
    } else {
      const PersistencePair &t3 = Barycenter[indexBaryCritType[Id2]];
      const double birth_barycenter = t3.birth.sfValue;
      const double death_barycenter = t3.death.sfValue;
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;

      } else {
        const PersistencePair &t2 = newData[indexDataCritType[Id1]];
        const double birth_data = t2.birth.sfValue;
        const double death_data = t2.death.sfValue;
        // direction[0] = t2.birth.sfValue - t3.birth.sfValue;
        // direction[1] = t2.death.sfValue - t3.death.sfValue;
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
      }
      directions[indexBaryCritType[Id2]] = std::move(direction);
      checker[indexBaryCritType[Id2]] = 1;
    }
  }


}



void PersistenceDiagramDictEncoding::computeAllDistances(
  std::vector<ttk::DiagramType> &barycentersList,
  const size_t nDiags,
  std::vector<ttk::DiagramType> &barycentersListMin,
  std::vector<ttk::DiagramType> &barycentersListSad,
  std::vector<ttk::DiagramType> &barycentersListMax,
  std::vector<BidderDiagram> &bidderBarycentersListMin,
  std::vector<BidderDiagram> &bidderBarycentersListSad,
  std::vector<BidderDiagram> &bidderBarycentersListMax,
  std::vector<std::vector<size_t>> &originIndexBarysMin,
  std::vector<std::vector<size_t>> &originIndexBarysSad,
  std::vector<std::vector<size_t>> &originIndexBarysMax,
  std::vector<BidderDiagram> &bidderDiagramsMin,
  std::vector<BidderDiagram> &bidderDiagramsMax,
  std::vector<BidderDiagram> &bidderDiagramsSad,
  std::vector<std::vector<ttk::MatchingType>> &matchingsDatasMin,
  std::vector<std::vector<ttk::MatchingType>> &matchingsDatasMax,
  std::vector<std::vector<ttk::MatchingType>> &matchingsDatasSad,
  std::vector<BidderDiagram> &trueBidderDiagramMin,
  std::vector<BidderDiagram> &trueBidderDiagramSad,
  std::vector<BidderDiagram> &trueBidderDiagramMax,
  std::vector<double> &allLossesAtEpoch,
  std::vector<double> &trueAllLossesAtEpoch,
  bool firstDistComputation) const{
    // setting BidderDiagram barycentersList
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; i++) {
      const auto &barycenter = barycentersList[i];

      for(size_t j = 0; j < barycenter.size(); ++j) {
        const ttk::PersistencePair &t = barycenter[j];
        const ttk::CriticalType nt1 = t.birth.type;
        const ttk::CriticalType nt2 = t.death.type;
        const double pers = t.persistence();
        // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

        if(pers > 0) {
          if(nt1 == CriticalType::Local_minimum
             && nt2 == CriticalType::Local_maximum) {
            barycentersListMin[i].emplace_back(t);
            originIndexBarysMin[i].push_back(j);
          } else {
            if(nt1 == CriticalType::Local_maximum
               || nt2 == CriticalType::Local_maximum) {
              barycentersListMax[i].emplace_back(t);
              originIndexBarysMax[i].push_back(j);
            }
            if(nt1 == CriticalType::Local_minimum
               || nt2 == CriticalType::Local_minimum) {
              barycentersListMin[i].emplace_back(t);
              originIndexBarysMin[i].push_back(j);
            }
            if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
               || (nt1 == CriticalType::Saddle2
                   && nt2 == CriticalType::Saddle1)) {
              barycentersListSad[i].emplace_back(t);
              originIndexBarysSad[i].push_back(j);
            }
          }
        }
      }
    }
    if(this->do_min_) {
      setBidderDiagrams(nDiags, barycentersListMin, bidderBarycentersListMin);
    }
    if(this->do_sad_) {
      setBidderDiagrams(nDiags, barycentersListSad, bidderBarycentersListSad);
    }
    if(this->do_max_) {
      setBidderDiagrams(nDiags, barycentersListMax, bidderBarycentersListMax);
    }


    // Compute distance and matchings
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; ++i) {
      std::vector<ttk::MatchingType> matching_min;
      std::vector<ttk::MatchingType> matching_sad;
      std::vector<ttk::MatchingType> matching_max;

      std::vector<ttk::MatchingType> matching_min_temp;
      std::vector<ttk::MatchingType> matching_sad_temp;
      std::vector<ttk::MatchingType> matching_max_temp;

      if(this->do_min_) {
        auto &barycentermin = bidderBarycentersListMin[i];
        auto &datamin = bidderDiagramsMin[i];
        auto &truedatamin = trueBidderDiagramMin[i];
        size_t sizeMin = truedatamin.size();

        //#ifdef TTK_ENABLE_OPENMP
        //#pragma omp atomic update
        //#endif // TTK_ENABLE_OPENMP
        if(firstDistComputation){
          allLossesAtEpoch[i]
            += computeDistance(datamin, barycentermin, matching_min);

          if((ProgApproach_) && (sizeMin != 0)) {
            trueAllLossesAtEpoch[i]
              += computeDistance(truedatamin, barycentermin, matching_min_temp);
          }
        } else {
          computeDistance(datamin, barycentermin, matching_min);
        }
      }
      if(this->do_max_) {
        auto &barycentermax = bidderBarycentersListMax[i];
        auto &datamax = bidderDiagramsMax[i];
        auto &truedatamax = trueBidderDiagramMax[i];
        size_t sizeMax = truedatamax.size();
        //#ifdef TTK_ENABLE_OPENMP
        //#pragma omp atomic update
        //#endif // TTK_ENABLE_OPENMP
        if(firstDistComputation){
          allLossesAtEpoch[i]
            += computeDistance(datamax, barycentermax, matching_max);

          if((ProgApproach_) && (sizeMax != 0)) {
            trueAllLossesAtEpoch[i]
              += computeDistance(truedatamax, barycentermax, matching_max_temp);
          }
        } else {
          computeDistance(datamax, barycentermax, matching_max);
        }
      }
      if(this->do_sad_) {
        auto &barycentersListad = bidderBarycentersListSad[i];
        auto &datasad = bidderDiagramsSad[i];
        auto &truedatasad = trueBidderDiagramSad[i];
        size_t sizeSad = truedatasad.size();
        //#ifdef TTK_ENABLE_OPENMP
        //#pragma omp atomic update
        //#endif // TTK_ENABLE_OPENMP
        if(firstDistComputation){
          allLossesAtEpoch[i]
            += computeDistance(datasad, barycentersListad, matching_sad);

          if((ProgApproach_) && (sizeSad != 0)) {
            trueAllLossesAtEpoch[i]
              += computeDistance(truedatasad, barycentersListad, matching_sad_temp);
          }
        } else {
          computeDistance(datasad, barycentersListad, matching_sad);
        }
      }
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }

  }