#include <InitDictBorder.h>

using namespace ttk;

void InitFarBorderDict::execute(std::vector<ttk::DiagramType> &DictDiagrams,
                                const std::vector<ttk::DiagramType> &datas,
                                const int nbAtoms,
                                bool do_min_,
                                bool do_sad_,
                                bool do_max_) {

  const size_t nDiags = datas.size();

  PersistenceDiagramDistanceMatrix MatrixCalculator;
  std::array<size_t, 2> nInputs{nDiags, 0};
  MatrixCalculator.setDos(do_min_, do_sad_, do_max_);
  MatrixCalculator.setThreadNumber(this->threadNumber_);
  Matrix distMatrix = MatrixCalculator.execute(datas, nInputs);
  // std::vector<int> indices;
  std::vector<double> allDistsSummed(nDiags);
  for(size_t i = 0; i < nDiags; ++i) {
    const auto &line = distMatrix[i];
    for(size_t j = 0; j < nDiags; ++j) {
      allDistsSummed[j] += line[j];
    }
  }
  // int Id1 = std::max_elemebnt(allDistsSummed.begin() ; allDistsSummed.end())
  std::vector<int> indices;
  int Id1 = std::max_element(allDistsSummed.begin(), allDistsSummed.end())
            - allDistsSummed.begin();
  indices.push_back(Id1);

  for(int i = 1; i < nbAtoms; ++i) {
    indices.push_back(getNextIndex(distMatrix, indices));
  }

  std::vector<double> tempDistsSummed(nbAtoms);
  for(int i = 0; i < nbAtoms; ++i) {
    tempDistsSummed[i] = allDistsSummed[indices[i]];
    std::cout << "VALUE: " << tempDistsSummed[i] << std::endl;
  }

  int FirstId = std::min_element(tempDistsSummed.begin(), tempDistsSummed.end())
                - tempDistsSummed.begin();
  DictDiagrams.push_back(datas[indices[FirstId]]);
  for(int i = 0; i < nbAtoms; ++i) {
    if(i == FirstId) {
      continue;
    } else {
      DictDiagrams.push_back(datas[indices[i]]);
    }
  }
  // DictDiagrams.resize(nbAtoms);
  // for(int i = 0; i < nbAtoms; ++i) {
  // const Diagram &atom = datas[indices[i]];
  // std::cout << "INDICE :" << indices[i] << std::endl;
  // DictDiagrams[i] = atom;
  //}
}

void InitFarBorderDict::setBidderDiagrams(
  const size_t nInputs,
  std::vector<ttk::DiagramType> &inputDiagrams,
  std::vector<BidderDiagram> &bidder_diags) const {

  bidder_diags.resize(nInputs);

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

double InitFarBorderDict::computeDistance(const BidderDiagram &D1,
                                          const BidderDiagram &D2) const {

  GoodDiagram D2_bis{};
  for(size_t i = 0; i < D2.size(); i++) {
    const Bidder &b = D2.at(i);
    Good g(b.x_, b.y_, b.isDiagonal(), D2_bis.size());
    g.SetCriticalCoordinates(b.coords_[0], b.coords_[1], b.coords_[2]);
    g.setPrice(0);
    D2_bis.emplace_back(g);
  }

  PersistenceDiagramAuction auction(
    this->Wasserstein, this->Alpha, this->Lambda, this->DeltaLim, true);
  auction.BuildAuctionDiagrams(D1, D2_bis);
  double loss;
  loss = auction.run();
  return loss;
}

int InitFarBorderDict::getNextIndex(const Matrix &distMatrix,
                                    const std::vector<int> &indices) const {
  std::vector<double> allSumCumul(distMatrix.size(), 0.);
  for(size_t k = 0; k < indices.size(); ++k) {
    const auto &line = distMatrix[indices[k]];
    for(size_t i = 0; i < distMatrix.size(); ++i) {
      if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
        allSumCumul[i] += 0.;
      } else {
        allSumCumul[i] += line[i];
      }
    }
  }
  int newId = std::max_element(allSumCumul.begin(), allSumCumul.end())
              - allSumCumul.begin();
  return newId;
}
