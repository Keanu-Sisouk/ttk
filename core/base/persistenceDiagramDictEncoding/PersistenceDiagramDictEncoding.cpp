#include <algorithm>

#include <PersistenceDiagramDictEncoding.h>

using namespace ttk;

void PersistenceDiagramDictEncoding::execute(
  const std::vector<Diagram> &intermediateDiagrams,
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  const std::array<size_t, 2> &nInputs) const {

  Timer tm{};

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
  //std::vector<BidderDiagram<double>> current_bidder_diagrams_min{};
  //std::vector<BidderDiagram<double>> current_bidder_diagrams_sad{};
  //std::vector<BidderDiagram<double>> current_bidder_diagrams_max{};

  // Store the persistence of the global min-max pair
  //std::vector<double> maxDiagPersistence(nDiags);

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
      //maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

      if(pers > 0) {
        if(nt1 == CriticalType::Local_minimum
           && nt2 == CriticalType::Local_maximum) {
          inputDiagramsMax[i].emplace_back(t);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
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

  //std::vector<std::vector<double>> distMat{};

  std::vector<Diagram> Barycenters{};
  std::vector<std::vector<std::vector<MatchingTuple>>> allMatchingsAtoms;
  std::vector<std::vector<MatchingTuple>> matchingsDatas;
  ConstrainedGradientDescent gradActor ;

  for(int epoch = 0; epoch < 50; ++epoch) {
    int k = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(int i = 0 ; i < nDiags ; ++i){
      Diagram &barycenter = Barycenters[i];
      std::vector<double> &weight = vectorWeights[i];
      std::vector<std::vector<MatchingTuple>> &matchings = allMatchingsAtoms[i];
      computeWeightedBarycenter(dictDiagrams , weight , barycenter , matchings);
    }
    std::vector<Diagram> BarycentersMin(nDiags);
    std::vector<Diagram> BarycentersSad(nDiags);
    std::vector<Diagram> BarycentersMax(nDiags);

    std::vector<BidderDiagram<double>> bidder_barycenters_min{};
    std::vector<BidderDiagram<double>> bidder_barycenters_sad{};
    std::vector<BidderDiagram<double>> bidder_barycenters_max{};

//setting BidderDiagram Barycenters
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
            //maxDiagPersistence[i] = std::max(pers, maxDiagPersistence[i]);

            if(pers > 0) {
              if(nt1 == CriticalType::Local_minimum
                && nt2 == CriticalType::Local_maximum) {
                  BarycentersMax[i].emplace_back(t);
              } else {
              if(nt1 == CriticalType::Local_maximum
                || nt2 == CriticalType::Local_maximum) {
                BarycentersMax[i].emplace_back(t);
              }
              if(nt1 == CriticalType::Local_minimum
                || nt2 == CriticalType::Local_minimum) {
                BarycentersMin[i].emplace_back(t);
                }
              if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
                || (nt1 == CriticalType::Saddle2
                  && nt2 == CriticalType::Saddle1)) {
                BarycentersSad[i].emplace_back(t);
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


    this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);
  // std::vector<double> gradient = compute
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for (size_t i = 0 ; i < nDiags ; ++i){
      std::vector<MatchingTuple> matching_min;
      std::vector<MatchingTuple> matching_sad;
      std::vector<MatchingTuple> matching_max;
      if (this->do_min_){
        auto &barycentermin = bidder_barycenters_min[i];
        auto &datamin = bidder_diagrams_min[i];
        computeDistance(datamin , barycentermin , matching_min);
      }
      if (this->do_max_){
        auto &barycentermax = bidder_barycenters_max[i];
        auto &datamax = bidder_diagrams_max[i];
        computeDistance(datamax , barycentermax , matching_max);
      }
      if (this->do_sad_){
        auto &barycentersad = bidder_barycenters_sad[i];
        auto &datasad = bidder_diagrams_sad[i];
        computeDistance(datasad , barycentersad , matching_sad);
      }
      matching_sad.insert(matching_sad.end() , std::make_move_iterator(matching_max.begin()) , std::make_move_iterator(matching_max.end()));
      matching_min.insert(matching_min.end() , std::make_move_iterator(matching_sad.begin()) , std::make_move_iterator(matching_sad.end()));
    //matchingsDatas[i] = matching_min.insert()
      matchingsDatas[i] = matching_min;
    }
    for (size_t i = 0 ; i < nDiags ; ++i){
      std::vector<double> gradWeight = computeGradientWeights(dictDiagrams, allMatchingsAtoms[i], Barycenters[i], intermediateDiagrams[i] , matchingsDatas[i]);
      int nb_points = Barycenters[i].size();
      gradActor.executeWeightsProjected(vectorWeights[i] , gradWeight , epoch , nb_points);
    }


    
  }// return distMat;
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

void PersistenceDiagramDictEncoding::computeDistance(
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
  auction.run(&matching);
}

std::vector<double> PersistenceDiagramDictEncoding::computeGradientWeights(
  const std::vector<Diagram> &dictDiagrams,
  const std::vector<std::vector<MatchingTuple>> &matchingsAtoms,
  const Diagram &Barycenter,
  const Diagram &newData,
  const std::vector<MatchingTuple> &matchings) const {

  //initialization
  std::vector<std::vector<DiagramTuple>> grad_list(Barycenter.size());
  //std::vector<MatchingTuple> matching;
  std::vector<std::vector<double>> directions(Barycenter.size());
  std::vector<double> gradient(dictDiagrams.size(), 0);

  //computing gradients
  for(int i = 0; i < matchingsAtoms.size(); ++i) {
    for(int j = 0; j < matchingsAtoms[i].size(); ++j) {
      const MatchingTuple &t = matchingsAtoms[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      const DiagramTuple &t2 = dictDiagrams[i][Id1];
      grad_list[Id2].push_back(t2);
    }
  }

  //computeDistance(newDataBidder, barycenterBidder, matching);
  for(int i = 0; i < matchings.size(); ++i) {
    const MatchingTuple &t = matchings[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    const DiagramTuple &t2 = newData[Id1];
    const DiagramTuple &t3 = Barycenter[Id2];
    std::vector<double> direction(2);
    direction[0] = std::get<6>(t2) - std::get<6>(t3);
    direction[1] = std::get<10>(t2) - std::get<10>(t3);
    // directions[Id2].push_back(direction);
    directions[Id2] = direction;
  }
  for(int i = 0; i < Barycenter.size(); ++i) {
    for(int j = 0; j < dictDiagrams.size(); ++j) {
      const DiagramTuple &t = grad_list[i][j];
      const double birth = std::get<6>(t);
      const double death = std::get<10>(t);
      const std::vector<double> &direction = directions[i];
      gradient[j] += -2 * (birth * direction[0] + death * direction[1]);
    }
  }
  return gradient;
}

std::vector<Matrice> PersistenceDiagramDictEncoding::computeGradientAtoms(
  const std::vector<double> &weights,
  const Diagram &Barycenter,
  const Diagram &newData,
  const std::vector<MatchingTuple> &matchings) const {

  //std::vector<MatchingTuple> matching;
  std::vector<Matrice> gradsLists(Barycenter.size());
  std::vector<std::vector<double>> directions(Barycenter.size());
  //computeDistance(newDataBidder, barycenterBidder, matching);
  for(int i = 0; i < matchings.size(); ++i) {
    const MatchingTuple &t = matchings[i];
    // Id in newData
    const SimplexId Id1 = std::get<0>(t);
    // Id in barycenter
    const SimplexId Id2 = std::get<1>(t);
    const DiagramTuple &t2 = newData[Id1];
    const DiagramTuple &t3 = Barycenter[Id2];
    std::vector<double> direction(2);
    direction[0] = std::get<6>(t2) - std::get<6>(t3);
    direction[1] = std::get<10>(t2) - std::get<10>(t3);
    // directions[Id2].push_back(direction);
    directions[Id2] = direction;
  }
  for(int i = 0; i < Barycenter.size(); ++i) {
    for(int j = 0; j < weights.size(); ++j) {
      std::vector<double> temp(2);
      const std::vector<double> &direction = directions[i];
      temp[0] = -2 * weights[j] * direction[0];
      temp[1] = -2 * weights[j] * direction[1];
      gradsLists[i].push_back(temp);
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

void PersistenceDiagramDictEncoding::setBidderDiagram(
  Diagram &inputDiagram, BidderDiagram<double> &bidder_diag) const {

  // auto &diag = inputDiagrams[i];
  // auto &bidders = bidder_diags[i];

  for(size_t j = 0; j < inputDiagram.size(); j++) {
    // Add bidder to bidders
    Bidder<double> b(inputDiagram[j], j, this->Lambda);
    b.setPositionInAuction(bidder_diag.size());
    bidder_diag.addBidder(b);
    if(b.isDiagonal() || b.x_ == b.y_) {
      this->printMsg(
        "Diagonal point in diagram !", ttk::debug::Priority::DETAIL);
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
