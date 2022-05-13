#include <algorithm>
#include <math.h>

#include <PersistenceDiagramDictEncoding.h>

using namespace ttk;

static bool testNeg(DiagramTuple &t) {
  double birth = std::get<6>(t);
  double death = std::get<10>(t);
  return death < birth;
}


void PersistenceDiagramDictEncoding::execute(
  std::vector<Diagram> &intermediateDiagrams,
  const std::vector<Diagram> &intermediateAtoms,
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const std::array<size_t, 2> &nInputs,
  const int seed,
  const int numAtom,
  std::vector<double> &loss_tab,
  std::vector<double> &true_loss_tab,
  std::vector<std::vector<double>> &allLosses,
  int percent_) {

  
  if(!ProgApproach){
  

    for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      std::cout << "SIZE OF DIAG " << i << " IS: " << intermediateDiagrams[i].size() << std::endl;
      if(sortedForTest){
        auto &diag = intermediateDiagrams[i];
        std::sort(diag.begin(), diag.end() , [](DiagramTuple &t1 , DiagramTuple &t2){
            return (std::get<10>(t1) - std::get<6>(t1)) > (std::get<10>(t2) - std::get<6>(t2));});
      }
    }

    
    std::vector<BidderDiagram<double>> bidder_diagram_min(intermediateDiagrams.size());
    std::vector<BidderDiagram<double>> bidder_diagram_sad(intermediateDiagrams.size());
    std::vector<BidderDiagram<double>> bidder_diagram_max(intermediateDiagrams.size());

    std::vector<std::vector<double>> histoVectorWeights(intermediateDiagrams.size());
    std::vector<Diagram> histoDictDiagrams(numAtom);
    this->maxLag2 = 10;
    Timer tm_init{};
    bool preWeightOpt = false;
    InitDictionary(dictDiagrams, intermediateDiagrams, intermediateAtoms, numAtom, this->do_min_, this->do_sad_, this->do_max_, seed);
    this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(), threadNumber_, debug::LineMode::NEW);
    method(intermediateDiagrams, intermediateAtoms, dictDiagrams, vectorWeights, nInputs, seed, numAtom, loss_tab, true_loss_tab, allLosses, histoVectorWeights, histoDictDiagrams, preWeightOpt, 0.01, bidder_diagram_min, bidder_diagram_sad, bidder_diagram_max);
  } else {
    for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      auto &diag = intermediateDiagrams[i];
      std::sort(diag.begin(), diag.end() , [](DiagramTuple &t1 , DiagramTuple &t2){
          return (std::get<10>(t1) - std::get<6>(t1)) > (std::get<10>(t2) - std::get<6>(t2));});
      
    }

    std::vector<BidderDiagram<double>> bidder_diagram_min{};
    std::vector<BidderDiagram<double>> bidder_diagram_sad{};
    std::vector<BidderDiagram<double>> bidder_diagram_max{};

    gettingBidderDiagrams(intermediateDiagrams, bidder_diagram_min , bidder_diagram_sad, bidder_diagram_max);

    //std::vector<double> percentages{0.2 , 0.15 , 0.1 , 0.05};
    //std::vector<double> percentages{0.8 , 0.6 , 0.5, 0.4, 0.3 , 0.2};
    //std::vector<double> percentages{0.3 , 0.2 , 0.1 , 0.05, 0.01};
    int start = 40;
    int stop = percent_;
    std::vector<double> percentages;
    for(int value = start ; value > stop ; value -=5){
      percentages.push_back(static_cast<double>(value)/100.);
    }
    percentages.push_back(static_cast<double>(percent_)/100.);
    for(size_t k = 0 ; k < percentages.size() ; ++k){
      std::cout << "PERCENT " << percentages[k] << std::endl;
    }
    //std::vector<double> percentages{0.4};
    std::vector<std::vector<double>> histoVectorWeights(intermediateDiagrams.size());
    std::vector<Diagram> histoDictDiagrams(numAtom);
    std::vector<Diagram> dataTemp(intermediateDiagrams.size());
    bool preWeightOpt = true;
    for(size_t j = 0 ; j < 1 ; ++j){
      double percentage = percentages[j];
      // std::vector<Diagram> dataTemp(intermediateDiagrams.size());
      // for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      // auto diag = intermediateDiagrams[i];
      // dataTemp[i] = diag;
      //}
      
      this->maxLag2 = 10;

      for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
        auto &diag = intermediateDiagrams[i];
        auto &t = diag[0];
        //double max_pers = std::get<10>(t) - std::get<6>(t);
        double max_pers = getMaxPers(diag);
        std::cout << "MAX PERS" << max_pers << std::endl;
        auto &diagTemp = dataTemp[i];
        dataTemp[i].push_back(t);
        for(size_t p = 1; p < diag.size(); ++p) {
          auto &t2 = diag[p];
          if(percentage * max_pers  <= (std::get<10>(t2) - std::get<6>(t2))) {
            dataTemp[i].push_back(t2);
          } else {
            continue;
          }
        }
        // double max_pers = std::get<10>(t) - std::get<6>(t);
        // diag.erase(std::remove_if(diag.begin() , diag.end(), [max_pers,
        // percentage](DiagramTuple &t){ return (std::get<10>(t) -
        // std::get<6>(t)) < percentage*max_pers;}), diag.end());
        std::cout << "SIZE OF DIAG " << i << " IS: " << diagTemp.size()
                  << std::endl;
      }
      
      Timer tm_init{};
      InitDictionary(dictDiagrams, dataTemp, intermediateAtoms, numAtom, this->do_min_, this->do_sad_, this->do_max_, seed);
      this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(), threadNumber_, debug::LineMode::NEW);

      method(dataTemp, intermediateAtoms, dictDiagrams, vectorWeights, nInputs, seed, numAtom, loss_tab, true_loss_tab, allLosses, histoVectorWeights, histoDictDiagrams, preWeightOpt, 0.01, bidder_diagram_min, bidder_diagram_sad, bidder_diagram_max);
      
    }



    int min_pairs_to_add = 0;
    std::vector<int> sizeCheck(intermediateDiagrams.size(), 0);
    int sum = 0;
    for(size_t i = 0; i < intermediateDiagrams.size() ; ++i){
      sum += sizeCheck[i];
    }
    //std::vector<int> newPercentages{1, 5 , 25, 50 , 100};
    //std::vector<int> newPercentages(10 , 10);
    int q = 1;
    //while(sum != static_cast<int>(intermediateDiagrams.size())){
    //for(size_t j = 1 ; j < newPercentages.size() ; ++j){
    //preWeightOpt = false;
    int counter = 0;
    for(size_t j = 1 ; j < percentages.size() ; ++j){
      double percentage = percentages[j];
      double previousPerc = percentages[j-1];
      if( j < percentages.size() - 1){
        this->maxLag2 = 10;
      } else {
        this->maxLag2 = 10;
      }


      //std::cout << "PERCENTAGE " << percentage
      // std::vector<Diagram> dataTemp(intermediateDiagrams.size());
      // for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
      // auto diag = intermediateDiagrams[i];
      // dataTemp[i] = diag;
      //}
      for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
        //if(sizeCheck[i] == 1){
          //continue;
        //}
        auto &diag = intermediateDiagrams[i];
        auto &diagTemp = dataTemp[i];
        int m = diag.size();

        //std::cout << "SIZE ORIGINAL " << m << std::endl;
        int n = diagTemp.size();
        //int counter = 0;
        auto &lastTuple = diagTemp[n - 1];
        //int max_pairs_to_add = 10 + m*10/100;
        int max_pairs_to_add = 5 + m*percentage/100;
        //std::cout << "MAX PAIRS TO ADD " << max_pairs_to_add << std::endl;
        //std::cout << "TRIPLE WHAT" << m*q/100 << std::endl;
        //std::cout << "WHAT " << static_cast<int> (m*(q/100)) << std::endl;
        double previousPers = std::get<10>(lastTuple) - std::get<6>(lastTuple);
        auto &t = diag[0];
        double max_pers = getMaxPers(diag);
        //double max_pers = std::get<10>(t) - std::get<6>(t);
        //auto &diagTemp = dataTemp[i];
        //dataTemp[i].push_back(t);
        for(size_t p = 0; p < diag.size(); ++p) {
          auto &t2 = diag[p];

          if(percentage * max_pers <= (std::get<10>(t2) - std::get<6>(t2))
             && (std::get<10>(t2) - std::get<6>(t2)) < previousPerc * max_pers) {
            dataTemp[i].push_back(t2);
            counter += 1;
          }
          
          //if((counter <= max_pairs_to_add) && ((std::get<10>(t2) - std::get<6>(t2)) < previousPers)) {
            //dataTemp[i].push_back(t2);
            //counter += 1;
          //}
        }
        
        // double max_pers = std::get<10>(t) - std::get<6>(t);
        // diag.erase(std::remove_if(diag.begin() , diag.end(), [max_pers,
        // percentage](DiagramTuple &t){ return (std::get<10>(t) -
        // std::get<6>(t)) < percentage*max_pers;}), diag.end());
        std::cout << "SIZE OF DIAG " << i << " IS: " << diagTemp.size()
                  << std::endl;
      }
      if(counter == 0){
        continue;
      }
      method(dataTemp, intermediateAtoms, dictDiagrams, vectorWeights, nInputs, seed, numAtom, loss_tab, true_loss_tab, allLosses, histoVectorWeights, histoDictDiagrams, preWeightOpt, 0.01, bidder_diagram_min, bidder_diagram_sad, bidder_diagram_max );
      //sum = 0;
      //for(size_t i = 0 ; i < intermediateDiagrams.size() ; ++i){
        //sum += sizeCheck[i];
      //}
      
      //q += 20;
    }
  }
}

void PersistenceDiagramDictEncoding::method(
  const std::vector<Diagram> &intermediateDiagrams,
  const std::vector<Diagram> &intermediateAtoms,
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const std::array<size_t, 2> &nInputs,
  const int seed,
  const int numAtom,
  std::vector<double> &loss_tab,
  std::vector<double> &true_loss_tab,
  std::vector<std::vector<double>> &allLosses,
  std::vector<std::vector<double>> &histoVectorWeights,
  std::vector<Diagram> &histoDictDiagrams,
  bool preWeightOpt,
  double acc,
  std::vector<BidderDiagram<double>> &true_bidder_diagram_min,
  std::vector<BidderDiagram<double>> &true_bidder_diagram_sad,
  std::vector<BidderDiagram<double>> &true_bidder_diagram_max) {

  Timer tm{};
  double tm_part = 0.;

  bool do_optimizeAtoms = false;
  bool do_optimizeWeights = false;

  if(OptimizeWeights) {
    printMsg("Weight Optimization activated");
    do_optimizeWeights = true;
  } else {
    printWrn("Weight Optimization desactivated");
  }
  if(OptimizeAtoms) {
    do_optimizeAtoms = true;
    printMsg("Atom Optimization activated");
  } else {
    printWrn("Atom Optimization desactivated");
  }

  //Timer tm_init{};
  //InitDictionary(dictDiagrams, intermediateDiagrams, intermediateAtoms, numAtom, this->do_min_,
  //               this->do_sad_, this->do_max_, seed);
  //this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
  //               threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);

  // for(size_t i = 0; i < dictDiagrams.size(); ++i) {
  //   std::cout << "SIZE: " << dictDiagrams[i].size();
  // }

  // for(size_t i = 0; i < dictDiagrams.size(); ++i) {
  //   std::cout << "Atom " << i << std::endl;
  //   for(size_t j = 0; j < dictDiagrams[i].size(); ++j) {
  //     DiagramTuple &t = dictDiagrams[i][j];
  //     std::cout << "Pair atoms: " << std::get<6>(t) << ", " <<
  //     std::get<10>(t)
  //               << std::endl;
  //   }
  // }

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
  std::vector<Diagram> inputDiagramsMin(nDiags);
  std::vector<Diagram> inputDiagramsSad(nDiags);
  std::vector<Diagram> inputDiagramsMax(nDiags);

  std::vector<BidderDiagram<double>> bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> bidder_diagrams_max{};

  std::vector<std::vector<size_t>> origin_index_datasMin(nDiags);
  std::vector<std::vector<size_t>> origin_index_datasSad(nDiags);
  std::vector<std::vector<size_t>> origin_index_datasMax(nDiags);

  // std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
  // std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
  // std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

  // Store the persistence of the global min-max pair
  // std::vector<double> maxDiagPersistence(nDiags);

  // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nDiags; i++) {
    const Diagram &CTDiagram = intermediateDiagrams[i];

    for(size_t j = 0; j < CTDiagram.size(); ++j) {
      const DiagramTuple &t = CTDiagram[j];
      const ttk::CriticalType nt1 = std::get<1>(t);
      const ttk::CriticalType nt2 = std::get<3>(t);
      const double pers = std::get<4>(t);
      // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax[i].emplace_back(t);
          origin_index_datasMax[i].push_back(j);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            origin_index_datasMax[i].push_back(j);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
            origin_index_datasMin[i].push_back(j);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
            origin_index_datasSad[i].push_back(j);
          }
        }
      }
    }
  }

  if(this->do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);
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
      pers << std::fixed << std::setprecision(2) << this->MinPersistence;
      this->printMsg("Using diagram pairs above a persistence threshold of "
                     + pers.str());
    } break;
    case ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence)))
        + "% most persistent pairs of every diagram");
      break;
    case ConstraintType::RELATIVE_PERSISTENCE_GLOBAL:
      this->printMsg(
        "Using the "
        + std::to_string(static_cast<int>(100 * (1 - this->MinPersistence)))
        + "% most persistent pairs of all diagrams");
      break;
  }

  // std::vector<std::vector<double>> distMat{};

  std::vector<Diagram> Barycenters(nDiags);
  std::vector<std::vector<std::vector<MatchingTuple>>> allMatchingsAtoms(
    nDiags);
  std::vector<std::vector<MatchingTuple>> matchingsDatasMin(nDiags);
  std::vector<std::vector<MatchingTuple>> matchingsDatasSad(nDiags);
  std::vector<std::vector<MatchingTuple>> matchingsDatasMax(nDiags);
  ConstrainedGradientDescent gradActor;
  double loss;
  double true_loss = 0.;
  // double loss1;
  // int epoch = 1;
  // std::vector<double> loss_tab;
  int lag = 0;
  int lag2 = 0;
  int lagLimit = 10;
  int MIN_EPOCH = 5;
  int MAX_EPOCH = MaxEpoch;
  bool cond = true;
  int epoch = 0;
  int nbEpochPrevious = loss_tab.size();
  std::vector<size_t> initSizes(dictDiagrams.size());
  for(size_t j = 0 ; j < dictDiagrams.size() ; ++j){
    initSizes[j] = dictDiagrams[j].size();

  }

  //std::vector<Diagram> histoDictDiagrams(dictDiagrams.size());
  std::vector<std::vector<int>> histoAllEpochLife(dictDiagrams.size());
  std::vector<std::vector<bool>> histoAllBoolLife(dictDiagrams.size());
  std::vector<std::vector<bool>> checkUnderDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkDiag(dictDiagrams.size());
  std::vector<std::vector<bool>> checkAboveGlobal(dictDiagrams.size());
  //std::vector<std::vector<double>> histoVectorWeights(nDiags);
  // std::vector<double> allLosses(nDiags , 0.);
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
    for(int i = 0; i < nDiags; ++i) {
      Diagram &barycenter = Barycenters[i];
      std::vector<double> &weight = vectorWeights[i];
      // double sum = 0.;
      // for(int q = 0; q < weight.size(); ++q) {
      //   sum += weight[q];
      //   std::cout << weight[q] << std::endl;
      // }
      // std::cout << "sum: " << sum << std::endl;
      // std::cout << "Poids: " << weight[0] << weight[1] << weight[2]
      //           << std::endl;
      // std::cout << "================================================="
      //          << std::endl;
      std::vector<std::vector<MatchingTuple>> &matchings = allMatchingsAtoms[i];
      computeWeightedBarycenter(
        dictDiagrams, weight, barycenter, matchings, ProgBarycenter);
      // std::cout << "Barycenter" << i << std::endl;
      // for(int j = 0; j < barycenter.size(); ++j) {
      //   DiagramTuple &t = barycenter[j];
      //   std::cout << "Pair: " << std::get<6>(t) << ", " << std::get<10>(t)
      //             << std::endl;
      // 
    }
    this->printMsg(
     "Computed 1st Barycenters for epoch " + std::to_string(epoch),
      epoch / static_cast<double>(MAX_EPOCH), tm_it.getElapsedTime(),
      threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);
    tm_part += static_cast<double>(tm_it.getElapsedTime());
    // this->printMsg(
    //   "====================BARYCENTER FINISHED======================");

    // tracking the original indices
    std::vector<Diagram> BarycentersMin(nDiags);
    std::vector<Diagram> BarycentersSad(nDiags);
    std::vector<Diagram> BarycentersMax(nDiags);

    std::vector<BidderDiagram<double>> bidder_barycenters_min{};
    std::vector<BidderDiagram<double>> bidder_barycenters_sad{};
    std::vector<BidderDiagram<double>> bidder_barycenters_max{};

    std::vector<std::vector<size_t>> origin_index_barysMin(nDiags);
    std::vector<std::vector<size_t>> origin_index_barysSad(nDiags);
    std::vector<std::vector<size_t>> origin_index_barysMax(nDiags);

    // setting BidderDiagram Barycenters
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; i++) {
      const Diagram &barycenter = Barycenters[i];

      for(size_t j = 0; j < barycenter.size(); ++j) {
        const DiagramTuple &t = barycenter[j];
        const ttk::CriticalType nt1 = std::get<1>(t);
        const ttk::CriticalType nt2 = std::get<3>(t);
        const double pers = std::get<4>(t);
        // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

        if(pers > 0) {
          if(nt1 == CriticalType::Local_minimum
             && nt2 == CriticalType::Local_maximum) {
            BarycentersMax[i].emplace_back(t);
            origin_index_barysMax[i].push_back(j);
          } else {
            if(nt1 == CriticalType::Local_maximum
               || nt2 == CriticalType::Local_maximum) {
              BarycentersMax[i].emplace_back(t);
              origin_index_barysMax[i].push_back(j);
            }
            if(nt1 == CriticalType::Local_minimum
               || nt2 == CriticalType::Local_minimum) {
              BarycentersMin[i].emplace_back(t);
              origin_index_barysMin[i].push_back(j);
            }
            if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
               || (nt1 == CriticalType::Saddle2
                   && nt2 == CriticalType::Saddle1)) {
              BarycentersSad[i].emplace_back(t);
              origin_index_barysSad[i].push_back(j);
            }
          }
        }
      }
    }
    if(this->do_min_) {
      setBidderDiagrams(nDiags, BarycentersMin, bidder_barycenters_min);
    }
    if(this->do_sad_) {
      setBidderDiagrams(nDiags, BarycentersSad, bidder_barycenters_sad);
    }
    if(this->do_max_) {
      setBidderDiagrams(nDiags, BarycentersMax, bidder_barycenters_max);
    }

    // Compute distance and matchings
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; ++i) {
      std::vector<MatchingTuple> matching_min;
      std::vector<MatchingTuple> matching_sad;
      std::vector<MatchingTuple> matching_max;

      
      std::vector<MatchingTuple> matching_min_temp;
      std::vector<MatchingTuple> matching_sad_temp;
      std::vector<MatchingTuple> matching_max_temp;
      
      if(this->do_min_) {
        auto &barycentermin = bidder_barycenters_min[i];
        auto &datamin = bidder_diagrams_min[i];
        auto &truedatamin = true_bidder_diagram_min[i];
        size_t s = truedatamin.size();

//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
        allLossesAtEpoch[i]
          += computeDistance(datamin, barycentermin, matching_min);

        if((ProgApproach) && (s != 0)){
          trueAllLossesAtEpoch[i] += computeDistance(truedatamin, barycentermin, matching_min_temp);
        }
      }
      if(this->do_max_) {
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];
        auto &truedatamax = true_bidder_diagram_max[i];
        size_t s = truedatamax.size();
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
        allLossesAtEpoch[i]
          += computeDistance(datamax, barycentermax, matching_max);
        
        if((ProgApproach) && (s != 0)){
          trueAllLossesAtEpoch[i] += computeDistance(truedatamax, barycentermax, matching_max_temp);
        }
      }
      if(this->do_sad_) {
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];
        auto &truedatasad = true_bidder_diagram_sad[i];
        size_t s = truedatasad.size();
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
        allLossesAtEpoch[i]
          += computeDistance(datasad, barycentersad, matching_sad);
        
        if((ProgApproach) && (s != 0)){
          trueAllLossesAtEpoch[i] += computeDistance(truedatasad, barycentersad, matching_sad_temp);
        }
      }
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }

    for(size_t p = 0 ; p < nDiags ; ++p){
      loss += allLossesAtEpoch[p];
    }

    for(size_t p = 0 ; p < nDiags ; ++p){
      true_loss += trueAllLossesAtEpoch[p];
    }


    for(size_t p = 0 ; p < nDiags ; ++p){

      allLosses[p].push_back(allLossesAtEpoch[p]);

    }


    loss_tab.push_back(loss);
    true_loss_tab.push_back(true_loss);

    
    if(preWeightOpt){
      if(epoch < 5){
        do_optimizeAtoms = false;
      } else {
        do_optimizeAtoms = true;
      }
    }


    double mini = *std::min_element(loss_tab.begin() + nbEpochPrevious, loss_tab.end() - 1);
    if(loss <= mini) {
      for(size_t p = 0; p < dictDiagrams.size(); ++p) {
        const auto atom = dictDiagrams[p];
        histoDictDiagrams[p] = atom;
      }
      for(size_t p = 0; p < nDiags; ++p) {
        const auto weights = vectorWeights[p];
        histoVectorWeights[p] = weights;
      }
      lag = 0;
    } else {
      if(epoch > MIN_EPOCH) {
        lag += 1;
      } else {
        lag = 0;
      }
    }


    // this->printMsg("LAG" + std::to_string(lag));
    // std::cout << "LAG" << lag << std::endl;
    if((epoch > MIN_EPOCH)
       && (loss_tab[epoch] / loss_tab[epoch - 1] > 0.99)) {
      if(loss_tab[epoch] < loss_tab[epoch - 1]) {
      //if(true){  
        if (lag2 == maxLag2){
          // lag = 0;
          this->printMsg("Loss not decreasing enough");
          if(StopCondition){
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
          lag2 +=1; 
        }
      }
    } else {
      lag2 = 0;
    }

    if(epoch > MIN_EPOCH && lag > lagLimit) {
      
      if(StopCondition){
        for(size_t p = 0; p < dictDiagrams.size(); ++p) {
          const auto atom = histoDictDiagrams[p];
          dictDiagrams[p] = atom;
        }
        for(size_t p = 0; p < nDiags; ++p) {
          const auto weights = histoVectorWeights[p];
          vectorWeights[p] = weights;
        }
        this->printMsg("Minimum not passed");
      //if(StopCondition){
        do_optimizeWeights = false;
        do_optimizeAtoms = false;

        cond = false;
      }
    }

    // if(epoch == 1) {
    //   loss1 = loss;
    // }

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
        const Diagram &Barycenter = Barycenters[i];
        const Diagram &Data = intermediateDiagrams[i];
        std::vector<Matrix> &hessianList = allHessianLists[i];
        const std::vector<MatchingTuple> &matchingsMin = matchingsDatasMin[i];
        const std::vector<MatchingTuple> &matchingsMax = matchingsDatasMax[i];
        const std::vector<MatchingTuple> &matchingsSad = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = origin_index_barysMin[i];
        const std::vector<size_t> &indexBarySad = origin_index_barysSad[i];
        const std::vector<size_t> &indexBaryMax = origin_index_barysMax[i];
        const std::vector<size_t> &indexDataMin = origin_index_datasMin[i];
        const std::vector<size_t> &indexDataSad = origin_index_datasSad[i];
        const std::vector<size_t> &indexDataMax = origin_index_datasMax[i];
        std::vector<double> &weights = vectorWeights[i];
        computeGradientWeights(
          gradWeights, hessianList, dictDiagrams, matchingsAtoms, Barycenter,
          Data, matchingsMin, matchingsMax, matchingsSad, indexBaryMin,
          indexBaryMax, indexBarySad, indexDataMin, indexDataMax, indexDataSad);
        int nb_points = Barycenter.size();
        gradActor.executeWeightsProjected(
          hessianList, weights, gradWeights, epoch, nb_points, MaxEigenValue);
      }

      this->printMsg("Computed 1st opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt1.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    


    std::vector<double> &weight = vectorWeights[0];


    for(size_t p = 0 ; p < nDiags ; ++p){
      allLossesAtEpoch[p] = 0.;
      trueAllLossesAtEpoch[p] = 0.;
    }
    Barycenters.clear();
    Barycenters.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);
    ////////////////////////////////ATOM////////////////////////////////////////

    // for(size_t k = 0 ; k < vectorWeights.size() ; ++k){
    //   std::vector<double> &weight = vectorWeights[k];
    //   std::cout << "Weights" << k << std::endl;
    //   for(size_t p = 0 ; p < weight.size() ; ++p){
    //     std::cout << weight[p] << std::endl;
    //   }
    // }

    // this->printMsg(
    // "========================ATOM NOW=============================");
    
    if(do_optimizeAtoms) {
      Timer tm_it2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(int i = 0; i < nDiags; ++i) {
        Diagram &barycenter = Barycenters[i];
        std::vector<double> &weight = vectorWeights[i];
      // double sum = 0.;
      // for(int q = 0; q < weight.size(); ++q) {
      //   sum += weight[q];
      //   std::cout << weight[q] << std::endl;
      // }
      // std::cout << "sum: " << sum << std::endl;
      // this->printMsg(std::to_string(sum_temp));
        std::vector<std::vector<MatchingTuple>> &matchings = allMatchingsAtoms[i];
        computeWeightedBarycenter(
          dictDiagrams, weight, barycenter, matchings, ProgBarycenter);
      // for(int i = 0; i < barycenter.size(); ++i) {
      // DiagramTuple &t = barycenter[i];
      // std::cout << "Pair: " << std::get<6>(t) << ", " << std::get<10>(t)
      //          << std::endl;
      // }
      }
      this->printMsg(
        "Computed 2nd Barycenters for epoch " + std::to_string(epoch),
        epoch / static_cast<double>(MAX_EPOCH), tm_it2.getElapsedTime(),
        threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);
      tm_part += static_cast<double>(tm_it2.getElapsedTime());

      BarycentersMin.clear();
      BarycentersSad.clear();
      BarycentersMax.clear();

      BarycentersMin.resize(nDiags);
      BarycentersSad.resize(nDiags);
      BarycentersMax.resize(nDiags);

      bidder_barycenters_min.clear();
      bidder_barycenters_sad.clear();
      bidder_barycenters_max.clear();

      origin_index_barysMin.clear();
      origin_index_barysSad.clear();
      origin_index_barysMax.clear();

      origin_index_barysMin.resize(nDiags);
      origin_index_barysSad.resize(nDiags);
      origin_index_barysMax.resize(nDiags);

      matchingsDatasMin.clear();
      matchingsDatasSad.clear();
      matchingsDatasMax.clear();

      matchingsDatasMin.resize(nDiags);
      matchingsDatasSad.resize(nDiags);
      matchingsDatasMax.resize(nDiags);
    // std::vector<BidderDiagram<double>> bidder_barycenters_min{};
    // std::vector<BidderDiagram<double>> bidder_barycenters_sad{};
    // std::vector<BidderDiagram<double>> bidder_barycenters_max{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; i++) {
        const Diagram &barycenter = Barycenters[i];

        for(size_t j = 0; j < barycenter.size(); ++j) {
          const DiagramTuple &t = barycenter[j];
          const ttk::CriticalType nt1 = std::get<1>(t);
          const ttk::CriticalType nt2 = std::get<3>(t);
          const double pers = std::get<4>(t);
        // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

          if(pers > 0) {
            if(nt1 == CriticalType::Local_minimum
              && nt2 == CriticalType::Local_maximum) {
              BarycentersMax[i].emplace_back(t);
              origin_index_barysMax[i].push_back(j);
            } else {
              if(nt1 == CriticalType::Local_maximum
                || nt2 == CriticalType::Local_maximum) {
                BarycentersMax[i].emplace_back(t);
                origin_index_barysMax[i].push_back(j);
              }
              if(nt1 == CriticalType::Local_minimum
                || nt2 == CriticalType::Local_minimum) {
                BarycentersMin[i].emplace_back(t);
                origin_index_barysMin[i].push_back(j);
              }
              if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
                || (nt1 == CriticalType::Saddle2
                    && nt2 == CriticalType::Saddle1)) {
                BarycentersSad[i].emplace_back(t);
                origin_index_barysSad[i].push_back(j);
              }
            }
          }
        }
      }
      if(this->do_min_) {
        setBidderDiagrams(nDiags, BarycentersMin, bidder_barycenters_min);
      }
      if(this->do_sad_) {
        setBidderDiagrams(nDiags, BarycentersSad, bidder_barycenters_sad);
      }
      if(this->do_max_) {
        setBidderDiagrams(nDiags, BarycentersMax, bidder_barycenters_max);
      }
    //double temp2 = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        std::vector<MatchingTuple> matching_min;
        std::vector<MatchingTuple> matching_sad;
        std::vector<MatchingTuple> matching_max;
        if(this->do_min_) {
          auto &barycentermin = bidder_barycenters_min[i];
          auto &datamin = bidder_diagrams_min[i];
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
          computeDistance(datamin, barycentermin, matching_min);
        }
        if(this->do_max_) {
          auto &barycentermax = bidder_barycenters_max[i];
          auto &datamax = bidder_diagrams_max[i];

//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
          computeDistance(datamax, barycentermax, matching_max);
        }
        if(this->do_sad_) {
          auto &barycentersad = bidder_barycenters_sad[i];
          auto &datasad = bidder_diagrams_sad[i];

//#ifdef TTK_ENABLE_OPENMP
//#pragma omp atomic update
//#endif // TTK_ENABLE_OPENMP
          computeDistance(datasad, barycentersad, matching_sad);
        }
        matchingsDatasMin[i] = std::move(matching_min);
        matchingsDatasSad[i] = std::move(matching_sad);
        matchingsDatasMax[i] = std::move(matching_max);
      }


    // this->printMsg("====================NOW ATOM
    // UPDATE======================"); ATOM OPTIMIZATION

    

    //if(do_optimizeAtoms) {
      std::vector<std::vector<std::vector<std::array<double, 2>>>>
        allPairToAddToGradList(nDiags);
      std::vector<std::vector<DiagramTuple>> allInfoToAdd(nDiags);
      std::vector<std::vector<Matrix>> gradsAtomsList(nDiags);
      std::vector<std::vector<int>> checkerAtomsList(nDiags);
      // for(size_t i = 0 ; i < nDiags ; ++i){
      //   auto &checkerAtoms = checkerAtomsList[i];
      //   checkerAtoms.resize()
      // }
      Timer tm_opt2{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < nDiags; ++i) {
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &gradsAtoms = gradsAtomsList[i];
        auto &checkerAtoms = checkerAtomsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const Diagram &Barycenter = Barycenters[i];
        const Diagram &Data = intermediateDiagrams[i];
        // std::vector<Matrix> &gradsAtoms = gradsAtomsList[i];
        const std::vector<MatchingTuple> &matchingsMin = matchingsDatasMin[i];
        const std::vector<MatchingTuple> &matchingsMax = matchingsDatasMax[i];
        const std::vector<MatchingTuple> &matchingsSad = matchingsDatasSad[i];
        const std::vector<size_t> &indexBaryMin = origin_index_barysMin[i];
        const std::vector<size_t> &indexBarySad = origin_index_barysSad[i];
        const std::vector<size_t> &indexBaryMax = origin_index_barysMax[i];
        const std::vector<size_t> &indexDataMin = origin_index_datasMin[i];
        const std::vector<size_t> &indexDataSad = origin_index_datasSad[i];
        const std::vector<size_t> &indexDataMax = origin_index_datasMax[i];
        const std::vector<double> &weights = vectorWeights[i];
        int nb_points = Barycenters[i].size();
        // std::vector<int> checkerAtoms(Barycenter.size(), 0);
        computeGradientAtoms(gradsAtoms, weights, Barycenter, Data,
                             matchingsMin, matchingsMax, matchingsSad,
                             indexBaryMin, indexBaryMax, indexBarySad,
                             indexDataMin, indexDataMax, indexDataSad,
                             checkerAtoms, pairToAddGradList, infoToAdd);
        // gradActor.executeAtoms(dictDiagrams, matchingsAtoms, Barycenter,
        //                        gradsAtoms, nb_points, checkerAtoms, epoch);
      }
      std::vector<double> maxiDeath(numAtom);
      std::vector<double> minBirth(numAtom);
      for(int j = 0; j < numAtom; ++j) {
        auto &atom = dictDiagrams[j];
        auto &temp = atom[0];
        maxiDeath[j] = std::get<10>(temp);
        minBirth[j] = std::get<6>(temp);
      }

      std::vector<std::vector<std::vector<int>>> allProjectionsList(nDiags);
      std::vector<std::vector<DiagramTuple>> allFeaturesToAdd(nDiags);
      std::vector<std::vector<std::array<double, 2>>> allProjLocations(nDiags);
      std::vector<std::vector<std::vector<double>>>
        allVectorForProjContributions(nDiags);
      std::vector<std::vector<DiagramTuple>> allTrueFeaturesToAdd(numAtom);
      std::vector<std::vector<std::vector<int>>> allTrueProj(numAtom);
      std::vector<std::vector<std::array<double, 2>>> allTrueProjLoc(numAtom);
      // std::vector<std::vector<int>> allAtomIndices(nDiags);

      for(size_t i = 0; i < nDiags; ++i) {
        //std::cout << " OPTIM DIAG: " << i << std::endl;
        auto &pairToAddGradList = allPairToAddToGradList[i];
        auto &infoToAdd = allInfoToAdd[i];
        auto &projForDiag = allProjectionsList[i];
        auto &vectorForProjContrib = allVectorForProjContributions[i];
        auto &featuresToAdd = allFeaturesToAdd[i];
        auto &projLocations = allProjLocations[i];
        auto &gradsAtoms = gradsAtomsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const Diagram &Barycenter = Barycenters[i];
        const auto &checkerAtoms = checkerAtomsList[i];
        int nb_points = Barycenters[i].size();
        gradActor.executeAtoms(
          dictDiagrams, matchingsAtoms, Barycenter, gradsAtoms, nb_points,
          checkerAtoms, epoch, projForDiag, featuresToAdd, projLocations,
          vectorForProjContrib, pairToAddGradList, infoToAdd);
        //std::cout << " OPTIM DIAG: " << i << std::endl;
        //std::cout << "==================================================" << std::endl;
      }
      
      if (CreationFeatures){
        // std::cout << "CREATING FEATURES" << std::endl;
        // double factEquiv = sqrt(static_cast<double>(numAtom));
        double factEquiv = static_cast<double>(numAtom);
        //double factEquiv = 1.;
        //double step = 1. / (sqrt(factEquiv) * 1e1);
        double step = 1. / (2. *  2. * factEquiv );
        for(size_t i = 0 ; i < nDiags ; ++i){
          auto &projForDiag = allProjectionsList[i];
          auto &featuresToAdd = allFeaturesToAdd[i];
          auto &projLocations = allProjLocations[i];
          auto &vectorForProjContrib = allVectorForProjContributions[i];
          for(size_t j = 0 ; j < projForDiag.size() ; ++j){
            DiagramTuple &t = featuresToAdd[j];
            std::array<double, 2> &pair = projLocations[j];
            std::vector<double> &vectorContrib = vectorForProjContrib[j];
            std::vector<int> &projAndIndex = projForDiag[j];
            std::vector<int> proj(numAtom);
            for(int m = 0; m < numAtom; ++m) {
              proj[m] = projAndIndex[m];
            }
            int atomIndex = static_cast<int>(projAndIndex[numAtom]);
            //auto it = std::find(allTrueProj[atomIndex].begin() , allTrueProj[atomIndex].end() , proj);
            //bool ralph = it != allTrueProj[atomIndex].end();
            bool lenNull = allTrueProj[atomIndex].size() == 0;
            if (lenNull){
              const CriticalType c1 = std::get<1>(t);
              const CriticalType c2 = std::get<3>(t);
              const SimplexId idTemp = std::get<5>(t);
              pair[0] = pair[0] - step * vectorContrib[0];
              pair[1] = pair[1] - step * vectorContrib[1];
              if(pair[0] > pair[1]) {
                pair[1] = pair[0];
              }
              if(pair[0] < minBirth[atomIndex]){
                continue;
              }
              DiagramTuple newPair{0,       c1,      0,  c2, pair[1] - pair[0],
                                   idTemp,  pair[0], 0., 0., 0.,
                                   pair[1], 0.,      0., 0.};
              allTrueProj[atomIndex].push_back(proj);
              allTrueFeaturesToAdd[atomIndex].push_back(newPair);
              allTrueProjLoc[atomIndex].push_back(pair);
            } else {
              // auto it = std::find(allTrueProj[atomIndex].begin() ,
              // allTrueProj[atomIndex].end() , proj); bool ralph = it ==
              // allTrueProj[atomIndex].end();
              bool ralph = true;
              size_t index = 0;
              if(Fusion) {
                for(size_t n = 0; n < allTrueProj[atomIndex].size(); ++n) {
                  auto &projStocked = allTrueProj[atomIndex][n];
                  auto &projLocStocked = allTrueProj[atomIndex][n];
                  double distance = sqrt(pow((pair[0] - projLocStocked[0]), 2)
                                         + pow((pair[1] - projLocStocked[1]), 2));
                  if(proj == projStocked && distance < 1e-3) {
                    ralph = false;
                    index = n;
                    break;
                  }
                }
              }
              if (ralph){

                const CriticalType c1 = std::get<1>(t);
                const CriticalType c2 = std::get<3>(t);
                const SimplexId idTemp = std::get<5>(t);
                pair[0] = pair[0] - step * vectorContrib[0];
                pair[1] = pair[1] - step * vectorContrib[1];
                if(pair[0] > pair[1]) {
                  pair[1] = pair[0];
                }
                if(pair[0] < minBirth[atomIndex]){
                  continue;
                }

                DiagramTuple newPair{0,       c1,      0,  c2, pair[1] - pair[0],
                                     idTemp,  pair[0], 0., 0., 0.,
                                     pair[1], 0.,      0., 0.};
                allTrueProj[atomIndex].push_back(proj);
                allTrueFeaturesToAdd[atomIndex].push_back(newPair);
                allTrueProjLoc[atomIndex].push_back(pair);

              } else {
                // auto index = std::distance(allTrueProj[atomIndex].begin() ,
                // it);
                auto &tReal = allTrueFeaturesToAdd[atomIndex][index];
                std::get<6>(tReal) = std::get<6>(tReal) - step * vectorContrib[0];
                std::get<10>(tReal)
                  = std::get<10>(tReal) - step * vectorContrib[1];
                if(std::get<6>(tReal) > std::get<10>(tReal)) {
                  std::get<10>(tReal) = std::get<6>(tReal);
                }
              }
            }
          }
        }
      
        //if (CreationFeatures){
        for(int i = 0 ; i < numAtom ; ++i){
          auto &atom = dictDiagrams[i];
          auto &histoEpochAtom = histoAllEpochLife[i];
          auto &histoBoolAtom = histoAllBoolLife[i];
          auto &boolUnderDiag = checkUnderDiag[i];
          auto &boolDiag = checkDiag[i];
          auto initSize = initSizes[i];
          auto &boolAboveGlobal = checkAboveGlobal[i];
          if(histoEpochAtom.size() > 0){
            for(size_t j = 0; j < histoEpochAtom.size(); ++j) {
              auto &t = atom[initSize + j];
              histoEpochAtom[j] += 1;
              histoBoolAtom[j] = std::get<10>(t) - std::get<6>(t) < 1e-3;
              boolDiag[j] = std::get<10>(t) - std::get<6>(t) < 1e-6;
              boolUnderDiag[j] = std::get<10>(t) < std::get<6>(t);
              boolAboveGlobal[j] = std::get<6>(t) > maxiDeath[i];
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
          for(size_t j = 0 ; j < trueFeaturesToAdd.size() ; ++j){
            auto &t = trueFeaturesToAdd[j];
            atom.push_back(t);
            histoEpochAtom.push_back(0);
            histoBoolAtom.push_back(std::get<10>(t) - std::get<6>(t) < 1e-3);
            boolDiag.push_back(std::get<10>(t) - std::get<6>(t) < 1e-6);
            boolUnderDiag.push_back(std::get<10>(t) < std::get<6>(t));
            boolAboveGlobal.push_back(std::get<6>(t) > maxiDeath[i]);
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
               || (histoEpochAtom[j] > 1000 && histoBoolAtom[j])) {
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
          //size_t initAtomSize = atom.size() -histoEpochAtom.size();
          if(static_cast<int>(indicesAtomToDelete.size()) > 0) {
            for(int j = static_cast<int>(indicesAtomToDelete.size()) - 1; j >= 0;
                j--) {
              atom.erase(atom.begin() + initSize  + indicesAtomToDelete[j]);
              histoEpochAtom.erase(histoEpochAtom.begin() + indicesAtomToDelete[j]);
              histoBoolAtom.erase(histoBoolAtom.begin() + indicesAtomToDelete[j]);
              boolUnderDiag.erase(boolUnderDiag.begin() + indicesAtomToDelete[j]);
              boolDiag.erase(boolDiag.begin() + indicesAtomToDelete[j]);
              boolAboveGlobal.erase(boolAboveGlobal.begin()
                                    + indicesAtomToDelete[j]);
            }
          } else {
            continue;
          }
        }
      }

      for(int i = 0; i < numAtom; ++i){
        auto &atom = dictDiagrams[i];
        for(size_t j = 0 ; j < atom.size() ; ++j){
          auto &t = atom[j];
          if(std::get<6>(t) > std::get<10>(t)){
            std::get<10>(t) = std::get<6>(t);
          }
        }
      }
      this->printMsg("Computed 2nd opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt2.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);


      // ATOM OPTIMIZATION
    }
    epoch +=1;
    // this->printMsg("=====================================================");
    Barycenters.clear();
    Barycenters.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);

    //for(size_t p = 0 ; p < dictDiagrams.size() ; ++p){
      //auto &atom = dictDiagrams[p];
      //std::cout << "ATOM: " + std::to_string(p) << std::endl;
      //for(size_t k = 0 ; k < atom.size() ; ++k){
        //auto &t = atom[k];
        //std::cout << "PAIR: " + std::to_string(std::get<6>(t)) + " AND " + std::to_string(std::get<10>(t)) << std::endl;
      //}
    //}



    printMsg(" Epoch "+std::to_string(epoch)+", loss = "+std::to_string(loss), 1, threadNumber_, ttk::debug::LineMode::REPLACE);
  }
  printMsg(" Epoch "+std::to_string(epoch)+", loss = "+std::to_string(loss), 1, threadNumber_);

  printMsg("Loss returned "
           + std::to_string(*std::min_element(
             loss_tab.begin() + nbEpochPrevious, loss_tab.end()))
           + " at Epoch "
           + std::to_string(std::min_element(loss_tab.begin() + nbEpochPrevious,
                                             loss_tab.end())
                            - loss_tab.begin()));

  // this->printMsg("Epoch" + std::to_string(epoch) + "==================");
  // this->printMsg("loss1 " + std::to_string(loss1) + "=================");
  // this->printMsg("loss " + std::to_string(loss) + "===================");

  for(size_t p = 0; p < dictDiagrams.size(); ++p) {
    auto atom = histoDictDiagrams[p];
    dictDiagrams[p] = atom;
  }
  for(size_t p = 0; p < nDiags; ++p) {
    auto weights = histoVectorWeights[p];
    vectorWeights[p] = weights;
  }

  // this->printMsg("time spent computing barycenter" + std::to_string(tm_part));
  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  // for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    // std::cout << "Atom " << i << std::endl;
    // for(size_t j = 0; j < dictDiagrams[i].size(); ++j) {
    //   DiagramTuple &t = dictDiagrams[i][j];
    //   std::cout << "Pair atoms: " << std::get<6>(t) << ", " << std::get<10>(t)
    //             << " and " << (0. <= std::get<6>(t)) << " and "
    //             << (std::get<10>(t) >= std::get<6>(t)) << std::endl;
    // }
  // }
}

double PersistenceDiagramDictEncoding::distVect(
  const std::vector<double> &vec1, const std::vector<double> &vec2) const {

  double dist = 0.;
  for(int i = 0; i < vec1.size(); ++i) {
    dist = dist + (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
  }
  return std::sqrt(dist);
}

double PersistenceDiagramDictEncoding::getMostPersistent(
  const std::vector<BidderDiagram<double>> &bidder_diags) const {

  double max_persistence = 0;

  for(unsigned int i = 0; i < bidder_diags.size(); ++i) {
    for(int j = 0; j < bidder_diags[i].size(); ++j) {
      const double persistence = bidder_diags[i].get(j).getPersistence();
      if(persistence > max_persistence) {
        max_persistence = persistence;
      }
    }
  }

  return max_persistence;
}

double PersistenceDiagramDictEncoding::computeDistance(
  const BidderDiagram<double> &D1,
  const BidderDiagram<double> &D2,
  std::vector<MatchingTuple> &matching) const {

  GoodDiagram<double> D2_bis{};
  for(int i = 0; i < D2.size(); i++) {
    const Bidder<double> &b = D2.get(i);
    Good<double> g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_x_, b.coords_y_, b.coords_z_);
    g.setPrice(0);
    D2_bis.addGood(g);
  }

  PersistenceDiagramAuction<double> auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(&D1, &D2_bis);
  double loss;
  loss = auction.run(&matching);
  return loss;
}

void PersistenceDiagramDictEncoding::computeGradientWeights(
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
  const std::vector<size_t> &indexDataSad) const {

  // initialization
  std::vector<std::vector<std::array<double, 2>>> grad_list(Barycenter.size());
  std::vector<std::vector<std::array<double, 2>>> pairToAddGradList;
  for(int i = 0; i < grad_list.size(); ++i) {
    grad_list[i].resize(matchingsAtoms.size());
  }
  // std::vector<MatchingTuple> matching;
  std::vector<std::array<double, 2>> directions(Barycenter.size());
  std::vector<std::array<double, 2>> data_assigned(Barycenter.size());
  // std::vector<double> gradient(dictDiagrams.size(), 0.);
  gradWeights.resize(dictDiagrams.size());
  for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    gradWeights[i] = 0.;
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  for(size_t j = 0 ; j < Barycenter.size() ; ++j){
    checker[j].resize(matchingsAtoms.size());
  }
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<int> tracker2(Barycenter.size(), 0);

  // this->printMsg("error?2");
  // computing gradients
  for(int i = 0; i < matchingsAtoms.size(); ++i) {
    // this->printMsg("Atom " + std::to_string(i));
    // this->printMsg("======================= atom size: "
    //             + std::to_string(static_cast<int>(dictDiagrams[i].size()))
    //             + ", and nb matchings: "
    //             + std::to_string(static_cast<int>(matchingsAtoms[i].size()))
    //             + ", and bary size: "
    //             + std::to_string(static_cast<int>(Barycenter.size()))
    //             + "=====================================");
    for(int j = 0; j < matchingsAtoms[i].size(); ++j) {
      const MatchingTuple &t = matchingsAtoms[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      // if(Id2 < 0) {
      if(Id2 < 0 || static_cast<int>(grad_list.size() <= Id2)
         || static_cast<int>(dictDiagrams[i].size()) <= Id1) {
        continue;
      } else if(Id1 < 0) {
        // this->printMsg("========DIAGONAL=========");
        const DiagramTuple &t3 = Barycenter[Id2];
        auto &point = grad_list[Id2][i];
        const double birth_barycenter = std::get<6>(t3);
        const double death_barycenter = std::get<10>(t3);
        // std::cout << "Barycenter Pair:" << birth_barycenter << " "
        //           << death_barycenter << std::endl;
        const double birth_death_atom
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        point[0] = birth_death_atom;
        point[1] = birth_death_atom;
        // std::cout << "Proj coordinates" << birth_death_atom << std::endl;
        //checker[Id2].push_back(i);
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      } else {
        // this->printMsg("====UPDATE GRADLIST========");
        const DiagramTuple &t2 = dictDiagrams[i][Id1];
        auto &point = grad_list[Id2][i];
        const double birth_atom = std::get<6>(t2);
        const double death_atom = std::get<10>(t2);
        point[0] = birth_atom;
        point[1] = death_atom;
        //checker[Id2].push_back(i);
        checker[Id2][i] = i;
        tracker[Id2] = 1;
      }
    }
  }

  // this->printMsg("===================PASSED============================");
  // computeDistance(newDataBidder, barycenterBidder, matching);

  int k = 0;
  for(int i = 0; i < matchingsMin.size(); ++i) {
    const MatchingTuple &t = matchingsMin[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);

    if(Id2 < 0) {
      k += 1;
      
      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataMin[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
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
      // this->printMsg(std::to_string(static_cast<int>(indexBaryMin[Id2])));
      const DiagramTuple &t3 = Barycenter[indexBaryMin[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      auto &direction = directions[indexBaryMin[Id2]];
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
        data_assigned[indexBaryMin[Id2]] = {birth_death_data, birth_death_data};

      } else {
        // checker[indexBaryMin[Id2]].push_back(indexDataMin[Id1]);
        const DiagramTuple &t2 = newData[indexDataMin[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        data_assigned[indexBaryMin[Id2]] = {birth_data, death_data};
        // directions[Id2].push_back(direction);
      }
      tracker2[indexBaryMin[Id2]] = 1;
    }
  }

  for(int i = 0; i < matchingsMax.size(); ++i) {
    const MatchingTuple &t = matchingsMax[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);

    if(Id2 < 0) {
      k += 1;
      
      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataMax[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::array<double , 2> direction;
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
      // this->printMsg(std::to_string(static_cast<int>(indexBaryMax[Id2])));
      const DiagramTuple &t3 = Barycenter[indexBaryMax[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      auto &direction = directions[indexBaryMax[Id2]];
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
        data_assigned[indexBaryMax[Id2]] = {birth_death_data, birth_death_data};
      } else {
        // checker[indexBaryMax[Id2]].push_back(indexDataMax[Id1]);
        const DiagramTuple &t2 = newData[indexDataMax[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        data_assigned[indexBaryMax[Id2]] = {birth_data, death_data};
        // directions[Id2].push_back(direction);
      }
      tracker2[indexBaryMax[Id2]] = 1;
    }
  }

  for(int i = 0; i < matchingsSad.size(); ++i) {
    const MatchingTuple &t = matchingsSad[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);

    if(Id2 < 0) {
      k += 1;
      
      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataSad[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
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
      // this->printMsg(std::to_string(static_cast<int>(indexBarySad[Id2])));
      const DiagramTuple &t3 = Barycenter[indexBarySad[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      auto &direction = directions[indexBarySad[Id2]];
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
        data_assigned[indexBarySad[Id2]] = {birth_death_data, birth_death_data};
      } else {
        // checker[indexBarySad[Id2]].push_back(indexDataSad[Id1]);
        const DiagramTuple &t2 = newData[indexDataSad[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        data_assigned[indexBarySad[Id2]] = {birth_data, death_data};
        // directions[Id2].push_back(direction);
      }
      tracker2[indexBarySad[Id2]] = 1;
    }
  }

  std::vector<int> temp(pairToAddGradList.size(), 1);
  grad_list.insert(grad_list.end(), pairToAddGradList.begin() , pairToAddGradList.end());
  tracker2.insert(tracker2.end(), temp.begin() , temp.end());
  tracker.insert(tracker.end() , temp.begin() , temp.end());
  std::vector<int> temp2(matchingsAtoms.size());
  for(size_t j = 0 ; j < matchingsAtoms.size() ; ++j){
    temp2.push_back(static_cast<int>(j));
  }
  for(size_t j = 0 ; j < pairToAddGradList.size() ; ++j){
    checker.push_back(temp2);
  }

  // this->printMsg("======================PASSED2==========================");
  for(int i = 0; i < grad_list.size(); ++i) {
    const auto &data_point = data_assigned[i];
    for(int j = 0; j < checker[i].size(); ++j) {
      auto &point = grad_list[i][checker[i][j]];
      point[0] -= data_point[0];
      point[1] -= data_point[1];
    }
  }

  // this->printMsg("error?3");
  for(int i = 0; i < grad_list.size(); ++i) {
    // for(auto it = std::begin(checker); it != std::end(checker); ++it) {
    if(tracker[i] == 0 || tracker2[i] == 0) {
      //this->printMsg("Not tracked");
      continue;
    } else {
      for(int j = 0; j < checker[i].size(); ++j) {
        // for(int j = 0; j < grad_list[i].size(); ++j) {
        // this->printMsg("======GRADIENT INDEX " +
        // std::to_string(checker[i][j])
        //               + "=======");
        const auto &point = grad_list[i][checker[i][j]];
        // const double birth = std::get<6>(t);
        // const double death = std::get<10>(t);
        const auto &direction = directions[i];
        // gradient[j] += -2 * (birth * direction[0] + death * direction[1]);
        gradWeights[checker[i][j]]
          += -2 * (point[0] * direction[0] + point[1] * direction[1]);
      }
    }
  }
  // this->printMsg("======================PASSED2==========================");
  hessianList.resize(grad_list.size());
  for(int i = 0; i < grad_list.size(); ++i) {
    Matrix &hessian = hessianList[i];
    hessian.resize(checker[i].size());
    for(int j = 0; j < checker[i].size(); ++j) {
      auto &line = hessian[j];
      // hessian[j].resize(checker[i].size());
      line.resize(checker[i].size());
      const auto &point = grad_list[i][checker[i][j]];
      for(int q = 0; q < checker[i].size(); ++q) {
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
  std::vector<int> &checker,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  std::vector<DiagramTuple> &infoToAdd) const {

  // std::vector<MatchingTuple> matching;
  gradsAtoms.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    gradsAtoms[i].resize(weights.size());
  }

  checker.resize(Barycenter.size());
  for(size_t i = 0; i < Barycenter.size(); ++i) {
    checker[i] = 0;
  }
  std::vector<std::vector<double>> directions(Barycenter.size());
  // std::vector<int> checker(Barycenter.size(), 0);
  // computeDistance(newDataBidder, barycenterBidder, matching);

  int k = 0;

  for(size_t i = 0; i < matchingsMin.size(); ++i) {
    const MatchingTuple &t = matchingsMin[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    // std::cout << Id1 << ""

    if(Id2 < 0) {
      k += 1;
      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataMin[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::vector<double> direction(2);
          direction[0] = birth_data - birth_death_barycenter;
          direction[1] = death_data - birth_death_barycenter;
          std::vector<std::vector<double>> temp3(weights.size());
          std::vector<double> temp2(2);
          for(size_t j = 0; j < weights.size(); ++j) {
            temp2[0] += -2 * weights[j] * direction[0];
            temp2[1] += -2 * weights[j] * direction[1];
            temp3[j] = temp2;
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
      const DiagramTuple &t3 = Barycenter[indexBaryMin[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        const DiagramTuple &t2 = newData[indexDataMin[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        // directions[Id2].push_back(direction);
      }
      directions[indexBaryMin[Id2]] = std::move(direction);
      checker[indexBaryMin[Id2]] = 1;
    }
  }

  for(size_t i = 0; i < matchingsMax.size(); ++i) {
    const MatchingTuple &t = matchingsMax[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);

    if(Id2 < 0) {
      k += 1;

      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataMax[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::vector<double> direction(2);
          direction[0] = birth_data - birth_death_barycenter;
          direction[1] = death_data - birth_death_barycenter;
          std::vector<std::vector<double>> temp3(weights.size());
          std::vector<double> temp2(2);
          for(size_t j = 0; j < weights.size(); ++j) {
            temp2[0] += -2 * weights[j] * direction[0];
            temp2[1] += -2 * weights[j] * direction[1];
            temp3[j] = temp2;
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
      // this->printMsg("==Here?==");
      // this->printMsg(std::to_string(static_cast<int>(indexBaryMax[Id2])));
      const DiagramTuple &t3 = Barycenter[indexBaryMax[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        const DiagramTuple &t2 = newData[indexDataMax[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        // directions[Id2].push_back(direction);
      }
      directions[indexBaryMax[Id2]] = std::move(direction);
      checker[indexBaryMax[Id2]] = 1;
    }
  }

  for(size_t i = 0; i < matchingsSad.size(); ++i) {
    const MatchingTuple &t = matchingsSad[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);

    if(Id2 < 0) {
      k += 1;
      if (Id1 < 0){
        continue;
      } else {
        if(CreationFeatures) {
          const DiagramTuple &t2 = newData[indexDataSad[Id1]];
          const double birth_data = std::get<6>(t2);
          const double death_data = std::get<10>(t2);
          const double birth_death_barycenter
            = birth_data + (death_data - birth_data) / 2.;
          std::vector<double> direction(2);
          direction[0] = birth_data - birth_death_barycenter;
          direction[1] = death_data - birth_death_barycenter;
          std::vector<std::vector<double>> temp3(weights.size());
          std::vector<double> temp2(2);
          for(size_t j = 0; j < weights.size(); ++j) {
            temp2[0] += -2 * weights[j] * direction[0];
            temp2[1] += -2 * weights[j] * direction[1];
            temp3[j] = temp2;
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
      // this->printMsg("==Here?==");
      // this->printMsg(std::to_string(static_cast<int>(indexBarySad[Id2])));
      const DiagramTuple &t3 = Barycenter[indexBarySad[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        const DiagramTuple &t2 = newData[indexDataSad[Id1]];
        const double birth_data = std::get<6>(t2);
        const double death_data = std::get<10>(t2);
        // direction[0] = std::get<6>(t2) - std::get<6>(t3);
        // direction[1] = std::get<10>(t2) - std::get<10>(t3);
        direction[0] = birth_data - birth_barycenter;
        direction[1] = death_data - death_barycenter;
        // directions[Id2].push_back(direction);
      }
      directions[indexBarySad[Id2]] = std::move(direction);
      checker[indexBarySad[Id2]] = 1;
    }
  }

  for(size_t i = 0; i < Barycenter.size(); ++i) {
    if(checker[i] == 0) {
      // this->printMsg("NOT CHECKED");
      // printf("NOT CHECKED")
      // std::cout << "NOT CHECKED" << std::endl;
      // continue;
      // this->printMsg("bouh");
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
  // return gradsLists;
}

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

void PersistenceDiagramDictEncoding::setBidderDiagrams(
  const size_t nInputs,
  std::vector<Diagram> &inputDiagrams,
  std::vector<BidderDiagram<double>> &bidder_diags) const {

  bidder_diags.resize(nInputs);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP

  for(size_t i = 0; i < nInputs; i++) {
    auto &diag = inputDiagrams[i];
    auto &bidders = bidder_diags[i];

    for(size_t j = 0; j < diag.size(); j++) {
      // Add bidder to bidders
      Bidder<double> b(diag[j], j, this->Lambda);
      b.setPositionInAuction(bidders.size());
      bidders.addBidder(b);
      if(b.isDiagonal() || b.x_ == b.y_) {
        this->printMsg("Diagonal point in diagram " + std::to_string(i) + "!",
                       ttk::debug::Priority::DETAIL);
      }
    }
  }
}

void PersistenceDiagramDictEncoding::enrichCurrentBidderDiagrams(
  const std::vector<BidderDiagram<double>> &bidder_diags,
  std::vector<BidderDiagram<double>> &current_bidder_diags,
  const std::vector<double> &maxDiagPersistence) const {

  current_bidder_diags.resize(bidder_diags.size());
  const auto nInputs = current_bidder_diags.size();
  const auto maxPersistence
    = *std::max_element(maxDiagPersistence.begin(), maxDiagPersistence.end());

  if(this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
     || this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL) {
    for(size_t i = 0; i < nInputs; ++i) {
      for(int j = 0; j < bidder_diags[i].size(); ++j) {
        auto b = bidder_diags[i].get(j);

        if( // filter out pairs below absolute persistence threshold
          (this->Constraint == ConstraintType::ABSOLUTE_PERSISTENCE
           && b.getPersistence() > this->MinPersistence)
          || // filter out pairs below persistence threshold relative to
          // the most persistent pair *of each diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
           && b.getPersistence() > this->MinPersistence * maxDiagPersistence[i])
          || // filter out pairs below persistence threshold relative to the
             // most persistence pair *in all diagrams*
          (this->Constraint == ConstraintType::RELATIVE_PERSISTENCE_GLOBAL
           && b.getPersistence() > this->MinPersistence * maxPersistence)) {
          b.id_ = current_bidder_diags[i].size();
          b.setPositionInAuction(current_bidder_diags[i].size());
          current_bidder_diags[i].addBidder(b);
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
    for(int j = 0; j < bidder_diags[i].size(); j++) {
      Bidder<double> b = bidder_diags[i].get(j);
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
      Bidder<double> b
        = bidder_diags[i].get(candidates_to_be_added[i][idx[i][j]]);
      const double persistence = b.getPersistence();
      if(persistence >= new_min_persistence) {
        b.id_ = current_bidder_diags[i].size();
        b.setPositionInAuction(current_bidder_diags[i].size());
        current_bidder_diags[i].addBidder(b);
      }
    }
  }
}

int PersistenceDiagramDictEncoding::InitDictionary(
  std::vector<ttk::Diagram> &dictDiagrams,
  const std::vector<ttk::Diagram> &datas,
  const std::vector<ttk::Diagram> &inputAtoms,
  const int nbAtom,
  bool do_min,
  bool do_sad,
  bool do_max,
  int seed) {
  switch(this->BackEnd) {
    case BACKEND::INPUT_ATOMS: {
      if(static_cast<int>(inputAtoms.size())!=nbAtom){
        printErr("number of atoms don't match, going with " + std::to_string(inputAtoms.size())+" atoms");
      }
      for(size_t i=0; i<inputAtoms.size(); i++){
        const auto &t = inputAtoms[i];
        dictDiagrams.push_back(t);
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
      for(size_t i = 0 ; i < datas.size() ; ++i){
        dictDiagrams[i] = datas[i];
      }
      const std::array<size_t , 2> nInputsUseless{100 , 100};
      while(static_cast<int>(dictDiagrams.size()) != nbAtom){
        //std::vector<Diagram> dictTemp;
        std::vector<double> allEnergy(dictDiagrams.size());
        
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t j = 0 ; j < dictDiagrams.size() ; ++j){
          std::vector<Diagram> dictTemp;
          std::vector<Diagram> dataAlone;
          std::vector<double> lossTabTemp;
          std::vector<double> trueLossTabTemp;
          std::vector<std::vector<double>> allLossesTemp(1);
           
          for(size_t p = 0 ; p < dictDiagrams.size() ; ++p){
            if(p != j){
              dictTemp.push_back(dictDiagrams[p]);
            } else {
              dataAlone.push_back(dictDiagrams[p]);
            }
          }
          std::vector<std::vector<double>> weightsTemp(1);
          std::vector<double> weights(dictTemp.size() , 1./ (static_cast<int>(dictTemp.size())* 1.));
          weightsTemp[0] = weights;
          std::vector<std::vector<double>> histoVectorWeights(1);
          std::vector<Diagram> histoDictDiagrams(dictTemp.size());
          std::vector<BidderDiagram<double>> bidderTempMin(dataAlone.size());
          std::vector<BidderDiagram<double>> bidderTempMax(dataAlone.size());
          std::vector<BidderDiagram<double>> bidderTempSad(dataAlone.size());
          this->method(dataAlone, inputAtoms, dictTemp, weightsTemp, nInputsUseless, seed, static_cast<int>(dictTemp.size()), lossTabTemp, trueLossTabTemp, allLossesTemp, histoVectorWeights, histoDictDiagrams, false, 0.01, bidderTempMin, bidderTempSad, bidderTempMax);
          double min_loss = *std::min_element(lossTabTemp.begin() , lossTabTemp.end());
          allEnergy[j] = min_loss;
        }
        int index = std::min_element(allEnergy.begin() , allEnergy.end()) - allEnergy.begin();
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
    const std::vector<ttk::Diagram> &intermediateDiagrams,
    std::vector<BidderDiagram<double>> &bidder_diagrams_min,
    std::vector<BidderDiagram<double>> &bidder_diagrams_sad,
    std::vector<BidderDiagram<double>> &bidder_diagrams_max){

    size_t nDiags = intermediateDiagrams.size();
    //double distance = 0.;
    
    std::vector<Diagram> inputDiagramsMin(nDiags);
    std::vector<Diagram> inputDiagramsSad(nDiags);
    std::vector<Diagram> inputDiagramsMax(nDiags);

    //std::vector<BidderDiagram<double>> bidder_diagrams_min{};
    //std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
    //std::vector<BidderDiagram<double>> bidder_diagrams_max{};

    //std::vector<std::vector<size_t>> origin_index_datasMin(nDiags);
    //std::vector<std::vector<size_t>> origin_index_datasSad(nDiags);
    //std::vector<std::vector<size_t>> origin_index_datasMax(nDiags);

    // std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
    // std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
    // std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

    // Store the persistence of the global min-max pair
    // std::vector<double> maxDiagPersistence(nDiags);

    // Create diagrams for min, saddle and max persistence pairs
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; i++) {
      const Diagram &CTDiagram = intermediateDiagrams[i];

      for(size_t j = 0; j < CTDiagram.size(); ++j) {
        const DiagramTuple &t = CTDiagram[j];
        const ttk::CriticalType nt1 = std::get<1>(t);
        const ttk::CriticalType nt2 = std::get<3>(t);
        const double pers = std::get<4>(t);
        // maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

        if(pers > 0) {
          if(nt1 == CriticalType::Local_minimum
            && nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            //origin_index_datasMax[i].push_back(j);
          } else {
            if(nt1 == CriticalType::Local_maximum
              || nt2 == CriticalType::Local_maximum) {
              inputDiagramsMax[i].emplace_back(t);
              //origin_index_datasMax[i].push_back(j);
            }
            if(nt1 == CriticalType::Local_minimum
              || nt2 == CriticalType::Local_minimum) {
              inputDiagramsMin[i].emplace_back(t);
              //origin_index_datasMin[i].push_back(j);
            }
            if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
              || (nt1 == CriticalType::Saddle2
                  && nt2 == CriticalType::Saddle1)) {
              inputDiagramsSad[i].emplace_back(t);
              //origin_index_datasSad[i].push_back(j);
            }
          }
        }
      }
    }

    if(this->do_min_) {
      setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
    }
    if(this->do_sad_) {
      setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
    }
    if(this->do_max_) {
      setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);
    }

    

    //return distance;
}

double PersistenceDiagramDictEncoding::getMaxPers(const Diagram &data){
  double max_pers = 0.;
  for(size_t j = 0 ; j < data.size() ; ++j){
    auto &t = data[j];
    double pers = std::get<10>(t) - std::get<6>(t);
    if(pers > max_pers){
      max_pers = pers;
    }
  }

  return max_pers;
}


