#include <algorithm>
#include <math.h>

#include <PersistenceDiagramDictEncoding.h>

using namespace ttk;

void PersistenceDiagramDictEncoding::execute(
  const std::vector<Diagram> &intermediateDiagrams,
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const std::array<size_t, 2> &nInputs,
  const int seed,
  const int numAtom) {

  Timer tm{};
  double tm_part = 0.;

  if(OptimizeWeights) {
    printMsg("Weight Optimization activated");
  } else {
    printWrn("Weight Optimization desactivated");
  }
  if(OptimizeAtoms) {
    printMsg("Atom Optimization activated");
  } else {
    printWrn("Atom Optimization desactivated");
  }

  Timer tm_init{};
  InitDictionary(dictDiagrams, intermediateDiagrams, numAtom, this->do_min_,
                 this->do_sad_, this->do_max_, seed);
  this->printMsg("Initialization computed ", 1, tm_init.getElapsedTime(),
                 threadNumber_, debug::LineMode::NEW, debug::Priority::DETAIL);

  for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    std::cout << "SIZE: " << dictDiagrams[i].size();
  }

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
  // double loss1;
  // int epoch = 1;
  std::vector<double> loss_tab;
  int lag = 0;
  int lagLimit = 50;
  int MAX_EPOCH = 500;
  bool cond = true;
  int epoch = 0;
  std::vector<Diagram> histoDictDiagrams(dictDiagrams.size());
  std::vector<std::vector<double>> histoVectorWeights(nDiags);
  std::vector<double> allLosses(nDiags , 0.);
  std::ofstream myFile("/home/keanu/ttk-data/weightsTimeLine2.csv");
  for(int j = 0 ; j < numAtom ; ++j){
    myFile << "weight" + std::to_string(j+1);
    if(j != numAtom - 1) myFile << ",";
  }
  myFile << "\n";


  std::ofstream lossHisto("/home/keanu/Bureau/python_trash/loss_histo.csv");
  lossHisto << "loss";
  lossHisto << "\n";

  std::ofstream allLossesEnd("/home/keanu/ttk-data/all_losses.csv");
  for(size_t j = 0 ; j < nDiags ; ++j){
    allLossesEnd << "diag" + std::to_string(j+1);
    if(j != nDiags - 1) allLossesEnd << ",";
  }
  allLossesEnd << "\n";

  // bool condition = true;
  while(cond && epoch < MAX_EPOCH) {
    // for(int epoch = 1; epoch < MAX_EPOCH; ++epoch) {
    
    loss = 0.;
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
      computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
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
      if(this->do_min_) {
        auto &barycentermin = bidder_barycenters_min[i];
        auto &datamin = bidder_diagrams_min[i];


        allLosses[i] += computeDistance(datamin, barycentermin, matching_min);
        
      }
      if(this->do_max_) {
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];


        allLosses[i] += computeDistance(datamax, barycentermax, matching_max);
        
      }
      if(this->do_sad_) {
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];


        allLosses[i] += computeDistance(datasad, barycentersad, matching_sad);
        
      }
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }

    for(size_t p = 0 ; p < nDiags ; ++p){
      loss+=allLosses[p];
    }


    for(size_t p = 0 ; p < nDiags ; ++p){
      allLossesEnd << allLosses[p];
      if(p != nDiags - 1) allLossesEnd << ","; 
    }
    allLossesEnd << "\n";

    loss_tab.push_back(loss);
    lossHisto << loss;
    lossHisto << "\n";

    double mini = *std::min_element(loss_tab.begin(), loss_tab.end() - 1);
    if(loss <= mini) {
      for(size_t p = 0; p < dictDiagrams.size(); ++p) {
        const auto &atom = dictDiagrams[p];
        histoDictDiagrams[p] = atom;
      }
      for(size_t p = 0; p < nDiags; ++p) {
        const auto &weights = vectorWeights[p];
        histoVectorWeights[p] = weights;
      }
      lag = 0;
    } else {
      lag += 1;
    }

    this->printMsg("LAG" + std::to_string(lag));
    // std::cout << "LAG" << lag << std::endl;
    if((epoch > 1) && (loss_tab[epoch] / loss_tab[epoch - 1] > 0.9999)) {
      if(loss_tab[epoch] < loss_tab[epoch - 1]) {
        this->printMsg("Loss not decreasing enough");
        OptimizeWeights = 0;
        OptimizeAtoms = 0;
        cond = false;
      }
    }

    if(lag > lagLimit) {
      // for(size_t p = 0; p < dictDiagrams.size(); ++p) {
      //   const auto &atom = histoDictDiagrams[p];
      //   dictDiagrams[p] = atom;
      // }
      // for(size_t p = 0; p < nDiags; ++p) {
      //   const auto &weights = histoVectorWeights[p];
      //   vectorWeights[p] = weights;
      // }
      this->printMsg("Minimum not passed");
      OptimizeWeights = 0;
      OptimizeAtoms = 0;
      cond = false;
    }

    // if(epoch == 1) {
    //   loss1 = loss;
    // }

    std::vector<std::vector<Matrix>> allHessianLists(nDiags);
    std::vector<std::vector<double>> gradWeightsList(nDiags);
    Timer tm_opt1{};

    // WEIGHT OPTIMIZATION
    if(OptimizeWeights) {
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
          hessianList, weights, gradWeights, epoch, nb_points);
      }

      this->printMsg("Computed 1st opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt1.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
    }

    std::vector<double> &weight = vectorWeights[0];
    for(int j = 0 ; j < numAtom ; ++j){
      myFile << weight.at(j);
      if(j!= numAtom - 1) myFile << ",";
    }
    myFile << "\n";

    for(size_t p = 0 ; p < nDiags ; ++p){
      allLosses[p] = 0.;
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
      computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
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
    double temp2 = 0;

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
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif // TTK_ENABLE_OPENMP
        temp2 += computeDistance(datamin, barycentermin, matching_min);
      }
      if(this->do_max_) {
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif // TTK_ENABLE_OPENMP
        temp2 += computeDistance(datamax, barycentermax, matching_max);
      }
      if(this->do_sad_) {
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif // TTK_ENABLE_OPENMP
        temp2 += computeDistance(datasad, barycentersad, matching_sad);
      }
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }

    // this->printMsg("====================NOW ATOM
    // UPDATE======================"); ATOM OPTIMIZATION
    if(OptimizeAtoms) {

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
        computeGradientAtoms(
          gradsAtoms, weights, Barycenter, Data, matchingsMin, matchingsMax,
          matchingsSad, indexBaryMin, indexBaryMax, indexBarySad, indexDataMin,
          indexDataMax, indexDataSad, checkerAtoms);
        // gradActor.executeAtoms(dictDiagrams, matchingsAtoms, Barycenter,
        //                        gradsAtoms, nb_points, checkerAtoms, epoch);
      }

      for(size_t i = 0; i < nDiags; ++i) {
        auto &gradsAtoms = gradsAtomsList[i];
        const auto &matchingsAtoms = allMatchingsAtoms[i];
        const Diagram &Barycenter = Barycenters[i];
        const auto &checkerAtoms = checkerAtomsList[i];
        int nb_points = Barycenters[i].size();

        gradActor.executeAtoms(dictDiagrams, matchingsAtoms, Barycenter,
                               gradsAtoms, nb_points, checkerAtoms, epoch);
      }

      this->printMsg("Computed 2nd opt for epoch " + std::to_string(epoch),
                     epoch / static_cast<double>(MAX_EPOCH),
                     tm_opt2.getElapsedTime(), threadNumber_,
                     debug::LineMode::NEW, debug::Priority::DETAIL);
      // ATOM OPTIMIZATION
    }

    // this->printMsg("=====================================================");
    Barycenters.clear();
    Barycenters.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);

    epoch += 1;
  }

  myFile.close();
  // this->printMsg("Epoch" + std::to_string(epoch) + "==================");
  // this->printMsg("loss1 " + std::to_string(loss1) + "=================");
  // this->printMsg("loss " + std::to_string(loss) + "===================");

  for(size_t p = 0; p < dictDiagrams.size(); ++p) {
    const auto &atom = histoDictDiagrams[p];
    dictDiagrams[p] = atom;
  }
  for(size_t p = 0; p < nDiags; ++p) {
    const auto &weights = histoVectorWeights[p];
    vectorWeights[p] = weights;
  }

  this->printMsg("time spent computing barycenter" + std::to_string(tm_part));
  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  for(size_t i = 0; i < dictDiagrams.size(); ++i) {
    std::cout << "Atom " << i << std::endl;
    for(size_t j = 0; j < dictDiagrams[i].size(); ++j) {
      DiagramTuple &t = dictDiagrams[i][j];
      std::cout << "Pair atoms: " << std::get<6>(t) << ", " << std::get<10>(t)
                << " and " << (0. <= std::get<6>(t)) << " and "
                << (std::get<10>(t) >= std::get<6>(t)) << std::endl;
    }
  }
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
        checker[Id2].push_back(i);
        tracker[Id2] = 1;
      } else {
        // this->printMsg("====UPDATE GRADLIST========");
        const DiagramTuple &t2 = dictDiagrams[i][Id1];
        auto &point = grad_list[Id2][i];
        const double birth_atom = std::get<6>(t2);
        const double death_atom = std::get<10>(t2);
        point[0] = birth_atom;
        point[1] = death_atom;
        checker[Id2].push_back(i);
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

  // this->printMsg("======================PASSED2==========================");
  for(int i = 0; i < Barycenter.size(); ++i) {
    const auto &data_point = data_assigned[i];
    for(int j = 0; j < checker[i].size(); ++j) {
      auto &point = grad_list[i][checker[i][j]];
      point[0] -= data_point[0];
      point[1] -= data_point[1];
    }
  }

  // this->printMsg("error?3");
  for(int i = 0; i < Barycenter.size(); ++i) {
    // for(auto it = std::begin(checker); it != std::end(checker); ++it) {
    if(tracker[i] == 0 || tracker2[i] == 0) {
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
  hessianList.resize(Barycenter.size());
  for(int i = 0; i < Barycenter.size(); ++i) {
    Matrix &hessian = hessianList[i];
    hessian.resize(checker[i].size());
    for(int j = 0; j < checker[i].size(); ++j) {
      auto &line = hessian[j];
      // hessian[j].resize(checker[i].size());
      line.resize(checker[i].size());
      const auto &point = grad_list[i][checker[i][j]];
      for(int k = 0; k < checker[i].size(); ++k) {
        const auto &point_temp = grad_list[i][checker[i][k]];
        line[k] = point[0] * point_temp[0] + point[1] * point_temp[1];
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
  std::vector<int> &checker) const {

  // std::vector<MatchingTuple> matching;
  gradsAtoms.resize(Barycenter.size());
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
      // this->printMsg("k = " + std::to_string(k));
    } else {
      const DiagramTuple &t3 = Barycenter[indexBaryMin[Id2]];
      const double birth_barycenter = std::get<6>(t3);
      const double death_barycenter = std::get<10>(t3);
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
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
    } else {
      for(size_t j = 0; j < weights.size(); ++j) {
        std::vector<double> temp(2);
        const std::vector<double> &direction = directions[i];
        temp[0] = -2 * weights[j] * direction[0];
        temp[1] = -2 * weights[j] * direction[1];
        gradsAtoms[i].push_back(temp);
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
  const int nbAtom,
  bool do_min,
  bool do_sad,
  bool do_max,
  int seed) {
  switch(this->BackEnd) {
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
    default:
      break;
  }
  return 0;
}
