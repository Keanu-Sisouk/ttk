#include <InitDictBorder.h>
// #include <Shuffle.h>
//
// #include <random>

using namespace ttk;

void InitFarBorderDict::execute(std::vector<Diagram> &DictDiagrams,
                                const std::vector<Diagram> &datas,
                                const int nbAtoms,
                                bool do_min_,
                                bool do_sad_,
                                bool do_max_) {

  const int nDiags = datas.size();

  // if(do_min_ && do_sad_ && do_max_) {
  //   this->printMsg("Processing all critical pairs types");
  // } else if(do_min_) {
  //   this->printMsg("Processing only MIN-SAD pairs");
  // } else if(do_sad_) {
  //   this->printMsg("Processing only SAD-SAD pairs");
  // } else if(do_max_) {
  //   this->printMsg("Processing only SAD-MAX pairs");
  // }

  // inputDiagrams = newDatas here
  // tracking the original indices
  std::vector<Diagram> inputDiagramsMin(nDiags);
  std::vector<Diagram> inputDiagramsSad(nDiags);
  std::vector<Diagram> inputDiagramsMax(nDiags);

  std::vector<BidderDiagram<double>> bidder_diagrams_min{};
  std::vector<BidderDiagram<double>> bidder_diagrams_sad{};
  std::vector<BidderDiagram<double>> bidder_diagrams_max{};

  for(int i = 0; i < nDiags; i++) {
    const Diagram &CTDiagram = datas[i];

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
          // origin_index_datasMax[i].push_back(j);
        } else {
          if(nt1 == CriticalType::Local_maximum
             || nt2 == CriticalType::Local_maximum) {
            inputDiagramsMax[i].emplace_back(t);
            // origin_index_datasMax[i].push_back(j);
          }
          if(nt1 == CriticalType::Local_minimum
             || nt2 == CriticalType::Local_minimum) {
            inputDiagramsMin[i].emplace_back(t);
            // origin_index_datasMin[i].push_back(j);
          }
          if((nt1 == CriticalType::Saddle1 && nt2 == CriticalType::Saddle2)
             || (nt1 == CriticalType::Saddle2
                 && nt2 == CriticalType::Saddle1)) {
            inputDiagramsSad[i].emplace_back(t);
            // origin_index_datasSad[i].push_back(j);
          }
        }
      }
    }
  }
  // setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
  // setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
  // setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);

  if(do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
  }
  if(do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
  }
  if(do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);
  }

  Matrix allDists(nDiags);
  for(int i = 0; i < nDiags; ++i) {
    auto &dists = allDists[i];
    for(int j = 0; j < nDiags; ++j) {
      dists.push_back(0.);
    }
  }
  std::vector<double> allDistsSummed(nDiags, 0.);
  for(int i = 0; i < nDiags; ++i) {
    std::vector<double> &dists = allDists[i];
    auto &datamin1 = bidder_diagrams_min[i];
    auto &datamax1 = bidder_diagrams_max[i];
    auto &datasad1 = bidder_diagrams_sad[i];
    for(int j = 0; j < nDiags; ++j) {
      auto &datamin2 = bidder_diagrams_min[j];
      auto &datamax2 = bidder_diagrams_max[j];
      auto &datasad2 = bidder_diagrams_sad[j];
      dists[j] += computeDistance(datamin1, datamin2);
      dists[j] += computeDistance(datamax1, datamax2);
      dists[j] += computeDistance(datasad1, datasad2);
      allDistsSummed[i] += dists[j];
    }
  }
  std::vector<int> indices;
  int Id1 = std::max_element(allDistsSummed.begin(), allDistsSummed.end())
            - allDistsSummed.begin();
  indices.push_back(Id1);
  // Diagram atom1 = datas[Id1];

  for(int i = 1; i < nbAtoms; ++i) {
    std::vector<double> distsToPtsSummed(nDiags, 0);
    for(int j = 0; j < nDiags; ++j) {
      if(std::find(indices.begin(), indices.end(), j) != indices.end()) {
        continue;
      } else {
        auto &dataMin = bidder_diagrams_min[j];
        auto &dataMax = bidder_diagrams_max[j];
        auto &dataSad = bidder_diagrams_sad[j];
        for(size_t k = 0; k < indices.size(); ++k) {
          auto &dataMinAtom = bidder_diagrams_min[indices[k]];
          auto &dataMaxAtom = bidder_diagrams_max[indices[k]];
          auto &dataSadAtom = bidder_diagrams_sad[indices[k]];
          double dist = computeDistance(dataMinAtom, dataMin)
                        + computeDistance(dataMaxAtom, dataMax)
                        + computeDistance(dataSadAtom, dataSad);
          distsToPtsSummed[j] += dist;
        }
      }
    }
    int newId
      = std::max_element(distsToPtsSummed.begin(), distsToPtsSummed.end())
        - distsToPtsSummed.begin();
    indices.push_back(newId);
  }

  DictDiagrams.resize(nbAtoms);
  for(int i = 0; i < nbAtoms; ++i) {
    const Diagram &atom = datas[indices[i]];
    DictDiagrams[i] = atom;
  }
}

void InitFarBorderDict::setBidderDiagrams(
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

double
  InitFarBorderDict::computeDistance(const BidderDiagram<double> &D1,
                                     const BidderDiagram<double> &D2) const {

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
  loss = auction.run();
  return loss;
}

// void InitRandomDict::execute(std::vector<Diagram> &DictDiagrams,
//                              const std::vector<Diagram> &datas,
//                              const int nbAtom,
//                              const int seed) {
//   int nDiags = datas.size();
//   DictDiagrams.resize(nbAtom);
//   std::vector<int> indices(nDiags);
//   std::iota(indices.begin(), indices.end(), 0);
//   std::mt19937 random_engine{};
//   random_engine.seed(seed);
//   ttk::shuffle(indices, random_engine);
//   for (int i = 0 ; i < nbAtom ; ++i){
//     const Diagram &atom = datas[indices[i]];
//     DictDiagrams[i] = atom;
//   }
// }
