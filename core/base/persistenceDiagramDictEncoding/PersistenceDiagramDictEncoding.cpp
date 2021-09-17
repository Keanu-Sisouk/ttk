#include <algorithm>

#include <PersistenceDiagramDictEncoding.h>

using namespace ttk;

void PersistenceDiagramDictEncoding::execute(
  const std::vector<Diagram> &intermediateDiagrams,
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const std::array<size_t, 2> &nInputs) const {

  Timer tm{};

  // for(size_t i = 0; i < dictDiagrams.size(); ++i) {
  // std::cout << "Atom " << i << std::endl;
  // for(size_t j = 0; j < dictDiagrams[i].size(); ++j) {
  // DiagramTuple &t = dictDiagrams[i][j];
  // std::cout << "Pair atoms: " << std::get<6>(t) << ", " << std::get<10>(t)
  //          << std::endl;
  //}
  //}

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
  double loss1;
  for(int epoch = 1; epoch < 100; ++epoch) {
    int k = 0;
    loss = 0;
    // this->printMsg("Epoch: " + std::to_string(epoch));
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif // TTK_ENABLE_OPENMP
    //////////////////////////////WEIGHTS/////////////////////////////////////////
    for(int i = 0; i < nDiags; ++i) {
      Diagram &barycenter = Barycenters[i];
      std::vector<double> &weight = vectorWeights[i];
      //  std::cout << "Poids: " << weight[0] << weight[1] << weight[2]
      //          << std::endl;
      // std::cout << "================================================="
      //          << std::endl;
      std::vector<std::vector<MatchingTuple>> &matchings = allMatchingsAtoms[i];
      computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
      // for(int i = 0; i < barycenter.size(); ++i) {
      //  DiagramTuple &t = barycenter[i];
      // std::cout << "Pair: " << std::get<6>(t) << ", " << std::get<10>(t)
      //        << std::endl;
      //}
    }
    this->printMsg("=============================WE ARE HERE "
                   "================================");
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
    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
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

    // this->printMsg("Complete", 1.0, tm.getElapsedTime(),
    // this->threadNumber_);
    // std::vector<double> gradient = compute
    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; ++i) {
      std::vector<MatchingTuple> matching_min;
      std::vector<MatchingTuple> matching_sad;
      std::vector<MatchingTuple> matching_max;
      if(this->do_min_) {
        auto &barycentermin = bidder_barycenters_min[i];
        auto &datamin = bidder_diagrams_min[i];
        loss += computeDistance(datamin, barycentermin, matching_min);
      }
      if(this->do_max_) {
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];
        loss += computeDistance(datamax, barycentermax, matching_max);
      }
      if(this->do_sad_) {
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];
        loss += computeDistance(datasad, barycentersad, matching_sad);
      }
      // matching_sad.insert(matching_sad.end(),
      //                    std::make_move_iterator(matching_max.begin()),
      //                    std::make_move_iterator(matching_max.end()));
      // matching_min.insert(matching_min.end(),
      //                    std::make_move_iterator(matching_sad.begin()),
      //                    std::make_move_iterator(matching_sad.end()));
      // matchingsDatas[i] = matching_min.insert()
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }

    if(epoch == 1) {
      loss1 = loss;
    }
    this->printMsg("====================PRINT WEIGHT======================");
    for(int j = 0; j < vectorWeights[3].size(); ++j) {
      std::cout << vectorWeights[3][j] << std::endl;
    }
    this->printMsg("====================PRINT WEIGHT======================");

    std::vector<std::vector<double>> gradientsWeights(nDiags);
    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; ++i) {
      const std::vector<std::vector<MatchingTuple>> &matchingsAtoms
        = allMatchingsAtoms[i];
      const Diagram &Barycenter = Barycenters[i];
      const Diagram &Data = intermediateDiagrams[i];
      const std::vector<MatchingTuple> &matchingsMin = matchingsDatasMin[i];
      const std::vector<MatchingTuple> &matchingsMax = matchingsDatasMax[i];
      const std::vector<MatchingTuple> &matchingsSad = matchingsDatasSad[i];
      const std::vector<size_t> &indexBaryMin = origin_index_barysMin[i];
      const std::vector<size_t> &indexBarySad = origin_index_barysSad[i];
      const std::vector<size_t> &indexBaryMax = origin_index_barysMax[i];
      const std::vector<size_t> &indexDataMin = origin_index_datasMin[i];
      const std::vector<size_t> &indexDataSad = origin_index_datasSad[i];
      const std::vector<size_t> &indexDataMax = origin_index_datasMax[i];

      gradientsWeights[i] = std::move(computeGradientWeights(
        dictDiagrams, matchingsAtoms, Barycenter, Data, matchingsMin,
        matchingsMax, matchingsSad, indexBaryMin, indexBaryMax, indexBarySad,
        indexDataMin, indexDataMax, indexDataSad));
    }

    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
    for(int i = 0; i < nDiags; ++i) {
      // const std::vector<double> gradWeight = computeGradientWeights(
      // dictDiagrams, matchingsAtoms, Barycenter, Data, matchingsMin,
      // matchingsMax, matchingsSad, indexBaryMin, indexBaryMax, indexBarySad,
      // indexDataMin, indexDataMax, indexDataSad);
      int nb_points = Barycenters[i].size();
      std::vector<double> &weights = vectorWeights[i];
      std::vector<double> &gradWeights = gradientsWeights[i];
      gradActor.executeWeightsProjected(weights, gradWeights, epoch, nb_points);
    }

    Barycenters.clear();
    Barycenters.resize(nDiags);
    allMatchingsAtoms.clear();
    allMatchingsAtoms.resize(nDiags);
    ////////////////////////////////ATOM////////////////////////////////////////

//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(threadNumber_)
//#endif // TTK_ENABLE_OPENMP
    for(int i = 0; i < nDiags; ++i) {
      Diagram &barycenter = Barycenters[i];
      std::vector<double> &weight = vectorWeights[i];
      // this->printMsg(std::to_string(sum_temp));
      std::vector<std::vector<MatchingTuple>> &matchings = allMatchingsAtoms[i];
      computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
      // for(int i = 0; i < barycenter.size(); ++i) {
      //  DiagramTuple &t = barycenter[i];
      //  std::cout << "Pair: " << std::get<6>(t) << ", " << std::get<10>(t)
      //            << std::endl;
      //}
    }

    this->printMsg("====================NOW ATOM======================");

    // std::vector<Diagram> BarycentersMin(nDiags);
    // std::vector<Diagram> BarycentersSad(nDiags);
    // std::vector<Diagram> BarycentersMax(nDiags);
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

    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
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

    //#ifdef TTK_ENABLE_OPENMP
    //#pragma omp parallel for num_threads(threadNumber_)
    //#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < nDiags; ++i) {
      std::vector<MatchingTuple> matching_min;
      std::vector<MatchingTuple> matching_sad;
      std::vector<MatchingTuple> matching_max;
      if(this->do_min_) {
        auto &barycentermin = bidder_barycenters_min[i];
        auto &datamin = bidder_diagrams_min[i];
        temp2 += computeDistance(datamin, barycentermin, matching_min);
      }
      if(this->do_max_) {
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];
        temp2 += computeDistance(datamax, barycentermax, matching_max);
      }
      if(this->do_sad_) {
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];
        temp2 += computeDistance(datasad, barycentersad, matching_sad);
      }
      // matching_sad.insert(matching_sad.end(),
      //                    std::make_move_iterator(matching_max.begin()),
      //                    std::make_move_iterator(matching_max.end()));
      // matching_min.insert(matching_min.end(),
      //                    std::make_move_iterator(matching_sad.begin()),
      //                    std::make_move_iterator(matching_sad.end()));
      // matchingsDatas[i] = matching_min.insert()
      matchingsDatasMin[i] = std::move(matching_min);
      matchingsDatasSad[i] = std::move(matching_sad);
      matchingsDatasMax[i] = std::move(matching_max);
    }
    this->printMsg("====================NOW ATOM UPDATE======================");

    for(size_t i = 0; i < nDiags; ++i) {
      std::vector<std::vector<MatchingTuple>> &matchingsAtoms
        = allMatchingsAtoms[i];
      Diagram &Barycenter = Barycenters[i];
      const Diagram &Data = intermediateDiagrams[i];
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
      std::vector<int> checkerAtoms(Barycenter.size(), 0);
      std::vector<Matrice> gradsAtoms = computeGradientAtoms(
        weights, Barycenter, Data, matchingsMin, matchingsMax, matchingsSad,
        indexBaryMin, indexBaryMax, indexBarySad, indexDataMin, indexDataMax,
        indexDataSad, checkerAtoms);
      gradActor.executeAtoms(
        dictDiagrams, matchingsAtoms, Barycenter, gradsAtoms, nb_points , checkerAtoms);
    }

    for(size_t i = 0; i < dictDiagrams.size(); ++i) {
      std::cout << "Atom " << i << std::endl;
      for(size_t j = 0; j < dictDiagrams[i].size(); ++j) {
        DiagramTuple &t = dictDiagrams[i][j];
        std::cout << "Pair atoms: " << std::get<6>(t) << ", " << std::get<10>(t)
                  << std::endl;
      }
    }

    this->printMsg("=====================================================");
    for(size_t i = 0; i < allMatchingsAtoms[0].size(); ++i) {
      std::cout << "matchings atom" << i << std::endl;
      for(size_t j = 0; j < allMatchingsAtoms[0][i].size(); ++j) {
        MatchingTuple &t = allMatchingsAtoms[0][i][j];
        std::cout << "Matching: " << std::get<0>(t) << ", " << std::get<1>(t)
                  << std::endl;
      }
    }
    this->printMsg("=====================================================");
    for(size_t i = 0; i < allMatchingsAtoms[1].size(); ++i) {
      std::cout << "matchings atom" << i << std::endl;
      for(size_t j = 0; j < allMatchingsAtoms[1][i].size(); ++j) {
        MatchingTuple &t = allMatchingsAtoms[1][i][j];
        std::cout << "Matching: " << std::get<0>(t) << ", " << std::get<1>(t)
                  << std::endl;
      }
    }

    this->printMsg("=====================================================");
    for(size_t i = 0; i < allMatchingsAtoms[2].size(); ++i) {
      std::cout << "matchings atom" << i << std::endl;
      for(size_t j = 0; j < allMatchingsAtoms[2][i].size(); ++j) {
        MatchingTuple &t = allMatchingsAtoms[2][i][j];
        std::cout << "Matching: " << std::get<0>(t) << ", " << std::get<1>(t)
                  << std::endl;
      }
    }

    this->printMsg("=====================================================");
    for(size_t i = 0; i < allMatchingsAtoms[3].size(); ++i) {
      std::cout << "matchings atom" << i << std::endl;
      for(size_t j = 0; j < allMatchingsAtoms[3][i].size(); ++j) {
        MatchingTuple &t = allMatchingsAtoms[3][i][j];
        std::cout << "Matching: " << std::get<0>(t) << ", " << std::get<1>(t)
                  << std::endl;
      }
    }

    this->printMsg("=====================================================");

  } // return distMat;
  this->printMsg("loss1 " + std::to_string(loss1) + "=================");
  this->printMsg("loss " + std::to_string(loss) + "===================");

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);
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

std::vector<double> PersistenceDiagramDictEncoding::computeGradientWeights(
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
  std::vector<Matrice> grad_list(Barycenter.size());
  for(int i = 0; i < grad_list.size(); ++i) {
    grad_list[i].resize(matchingsAtoms.size());
  }
  // std::vector<MatchingTuple> matching;
  std::vector<std::vector<double>> directions(Barycenter.size());
  std::vector<double> gradient(dictDiagrams.size(), 0.);
  std::vector<std::vector<int>> checker(Barycenter.size());
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<int> tracker2(Barycenter.size(), 0);

  // this->printMsg("error?2");
  // computing gradients
  for(int i = 0; i < matchingsAtoms.size(); ++i) {
    this->printMsg(
      "======================="
      + std::to_string(static_cast<int>(dictDiagrams[i].size())) + " and "
      + std::to_string(static_cast<int>(matchingsAtoms[i].size())) + " and "
      + std::to_string(static_cast<int>(Barycenter.size()))
      + "=====================================");
    for(int j = 0; j < matchingsAtoms[i].size(); ++j) {
      const MatchingTuple &t = matchingsAtoms[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      //if(Id2 < 0){
      if(Id2 < 0 || Id2 >= grad_list.size() || Id1 >= dictDiagrams[i].size()) {
        continue;
      } else if(Id1 < 0) {
        const DiagramTuple &t3 = Barycenter[Id2];
        std::vector<double> point(2);
        const double birth_barycenter = std::get<6>(t3);
        const double death_barycenter = std::get<10>(t3);
        const double birth_death_atom
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        point[0] = birth_death_atom;
        point[1] = birth_death_atom;
        // this->printMsg(std::to_string(Id2));
        // grad_list[Id2].push_back(point);
        grad_list[Id2][i] = std::move(point);
        checker[Id2].push_back(i);
        tracker[Id2] = 1;
      } else {
        // this->printMsg("====UPDATE GRADLIST========");
        const DiagramTuple &t2 = dictDiagrams[i][Id1];
        std::vector<double> point(2);
        const double birth_atom = std::get<6>(t2);
        const double death_atom = std::get<10>(t2);
        point[0] = birth_atom;
        point[1] = death_atom;
        // this->printMsg(std::to_string(Id2));
        // grad_list[Id2].push_back(point);
        grad_list[Id2][i] = std::move(point);
        checker[Id2].push_back(i);
        tracker[Id2] = 1;
      }
    }
  }

  this->printMsg("===================PASSED============================");
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
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        // checker[indexBaryMin[Id2]].push_back(indexDataMin[Id1]);
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
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        // checker[indexBaryMax[Id2]].push_back(indexDataMax[Id1]);
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
      std::vector<double> direction(2);
      if(Id1 < 0) {
        const double birth_death_data
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2;
        direction[0] = birth_death_data - birth_barycenter;
        direction[1] = birth_death_data - death_barycenter;
      } else {
        // checker[indexBarySad[Id2]].push_back(indexDataSad[Id1]);
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
      tracker2[indexBarySad[Id2]] = 1;
    }
  }

  this->printMsg("======================PASSED2==========================");

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
        const std::vector<double> &point = grad_list[i][checker[i][j]];
        // const double birth = std::get<6>(t);
        // const double death = std::get<10>(t);
        const std::vector<double> &direction = directions[i];
        // gradient[j] += -2 * (birth * direction[0] + death * direction[1]);
        gradient[checker[i][j]]
          += -2 * (point[0] * direction[0] + point[1] * direction[1]);
      }
    }
  }
  return gradient;
}

std::vector<Matrice> PersistenceDiagramDictEncoding::computeGradientAtoms(
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
  std::vector<Matrice> gradsLists(Barycenter.size());
  std::vector<std::vector<double>> directions(Barycenter.size());
  //std::vector<int> checker(Barycenter.size(), 0);
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
      this->printMsg(std::to_string(static_cast<int>(indexBaryMin[Id2])));
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
      this->printMsg(std::to_string(static_cast<int>(indexBaryMax[Id2])));
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
      this->printMsg(std::to_string(static_cast<int>(indexBarySad[Id2])));
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

  for(int i = 0; i < Barycenter.size(); ++i) {
    if(checker[i] == 0) {
      continue;
    } else {
      for(int j = 0; j < weights.size(); ++j) {
        std::vector<double> temp(2);
        const std::vector<double> &direction = directions[i];
        temp[0] = -2 * weights[j] * direction[0];
        temp[1] = -2 * weights[j] * direction[1];
        gradsLists[i].push_back(temp);
      }
    }
  }
  return gradsLists;
}

void PersistenceDiagramDictEncoding::setBidderDiagrams(
  const size_t nInputs,
  std::vector<Diagram> &inputDiagrams,
  std::vector<BidderDiagram<double>> &bidder_diags) const {

  bidder_diags.resize(nInputs);

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
