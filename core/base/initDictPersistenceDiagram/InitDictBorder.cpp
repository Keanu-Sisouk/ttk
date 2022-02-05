#include <InitDictBorder.h>

using namespace ttk;

void InitFarBorderDict::execute(std::vector<Diagram> &DictDiagrams,
                                const std::vector<Diagram> &datas,
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

  DictDiagrams.resize(nbAtoms);
  for(int i = 0; i < nbAtoms; ++i) {
    const Diagram &atom = datas[indices[i]];
    DictDiagrams[i] = atom;
  }

  // this->printMsg("Initialisation time", 1,
  // tm1.getElapsedTime(),threadNumber_, debug::LineMode::NEW,
  // debug::Priority::DETAIL);
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
