#include <InitDictPersistenceDiagram.h>

using namespace ttk;

void InitFarBorderDict::execute(
  std::vector<Diagram> &DictDiagrams,
  const std::vector<Diagram> &datas,
  const int nbAtoms){

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

  // std::vector<std::vector<size_t>> origin_index_datasMin(nDiags);
  // std::vector<std::vector<size_t>> origin_index_datasSad(nDiags);
  // std::vector<std::vector<size_t>> origin_index_datasMax(nDiags);



  for(size_t i = 0; i < nDiags; i++) {
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


  if(this->do_min_) {
    setBidderDiagrams(nDiags, inputDiagramsMin, bidder_diagrams_min);
  }
  if(this->do_sad_) {
    setBidderDiagrams(nDiags, inputDiagramsSad, bidder_diagrams_sad);
  }
  if(this->do_max_) {
    setBidderDiagrams(nDiags, inputDiagramsMax, bidder_diagrams_max);
  }

  std::vector<std::vector<double>> allDists(nDiags);
  for(int i = 0 ; i < nDiags ; ++i){
    std::vector<double> &dists = allDists[i];
    for (int j = 0 ; j < nDiags ; ++j){
      dists.push_back(0.);
    }
  }
  std::vector<double> allDistsSummed(nDiags, 0.);
  for (int i = 0 ; i < nDiags ; ++i){
    std::vector<double> &dists = allDists[i];
    auto &datamin1 = bidder_diagrams_min[i];
    auto &datamax1 = bidder_diagrams_max[i];
    auto &datasad1 = bidder_diagrams_sad[i];
    for (int j = 0 ; j < nDiags ; ++j){
      auto &datamin2 = bidder_diagrams_min[j];
      auto &datamax2 = bidder_diagrams_max[j];
      auto &datasad2 = bidder_diagrams_sad[j];
      dists[j]+= computeDistance(datamin1 , datamin2);
      dists[j]+= computeDistance(datamax1 , datamax2);
      dists[j]+= computeDistance(datasad1 , datasad2);
      allDistsSummed[i]+= dists[j];
    }
  }
  if (nbAtoms == 2){
    int Id1 = std::max_element(allDistsSummed.begin() , allDistsSummed.end()) - allDistsSummed.begin();
    Diagram atom1 = datas[Id1];
    DictDiagrams.push_back(atom1);
    int Id2 = std::max_element(allDists[Id1].begin() , allDists[Id1].end()) - allDists[Id1].begin();
    Diagram atom2 = datas[Id2];
    DictDiagrams.push_back(atom2);
  }

  if(nbAtoms == 3){
    int Id1 = std::max_element(allDistsSummed.begin() , allDistsSummed.end()) - allDistsSummed.begin();
    Diagram atom1 = datas[Id1];
    DictDiagrams.push_back(atom1);
    int Id2 = std::max_element(allDists[Id1].begin() , allDists[Id1].end()) - allDists[Id1].begin();
    Diagram atom2 = datas[Id2];
    DictDiagrams.push_back(atom2);
    std::vector<double> distsTo2Pts(nDiags , 0.);
    for(int i = 0 ; i < nDiags ; ++i){
      if (i == Id1 || i == Id2){
        continue;
      } else {
        auto &dataMinAtom1 = bidder_diagrams_min[Id1];
        auto &dataMaxAtom1 = bidder_diagrams_max[Id1];
        auto &dataSadAtom1 = bidder_diagrams_sad[Id1];
        auto &dataMinAtom2 = bidder_diagrams_min[Id2];
        auto &dataMaxAtom2 = bidder_diagrams_max[Id2];
        auto &dataSadAtom2 = bidder_diagrams_sad[Id2];
        auto &dataMin = bidder_diagrams_min[i];
        auto &dataMax = bidder_diagrams_max[i];
        auto &dataSad = bidder_diagrams_sad[i];
        distsTo2Pts[i] = computeDistance(dataMinAtom1 , dataMin) +
          computeDistance(dataMaxAtom1 , dataMax) +
          computeDistance(dataSadAtom1 , dataSad) +
          computeDistance(dataMinAtom2 , dataMin) +
          computeDistance(dataMaxAtom2 , dataMax) +
          computeDistance(dataSadAtom2 , dataSad);
      }
    }
    int Id3 = std::max_element(distsTo2Pts.begin() , distsTo2Pts.end()) - distsTo2Pts.begin();
    Diagram atom3 = datas[Id3];
    DictDiagrams.push_back(atom3);
  }

  if(nbAtoms == 4){
    int Id1 = std::max_element(allDistsSummed.begin() , allDistsSummed.end()) - allDistsSummed.begin();
    Diagram atom1 = datas[Id1];
    DictDiagrams.push_back(atom1);
    int Id2 = std::max_element(allDists[Id1].begin() , allDists[Id1].end()) - allDists[Id1].begin();
    Diagram atom2 = datas[Id2];
    DictDiagrams.push_back(atom2);
    std::vector<double> distsTo2Pts(nDiags , 0.);
    for(int i = 0 ; i < nDiags ; ++i){
      if (i == Id1 || i == Id2){
        continue;
      } else {
        auto &dataMinAtom1 = bidder_diagrams_min[Id1];
        auto &dataMaxAtom1 = bidder_diagrams_max[Id1];
        auto &dataSadAtom1 = bidder_diagrams_sad[Id1];
        auto &dataMinAtom2 = bidder_diagrams_min[Id2];
        auto &dataMaxAtom2 = bidder_diagrams_max[Id2];
        auto &dataSadAtom2 = bidder_diagrams_sad[Id2];
        auto &dataMin = bidder_diagrams_min[i];
        auto &dataMax = bidder_diagrams_max[i];
        auto &dataSad = bidder_diagrams_sad[i];
        distsTo2Pts[i] = computeDistance(dataMinAtom1 , dataMin) +
          computeDistance(dataMaxAtom1 , dataMax) +
          computeDistance(dataSadAtom1 , dataSad) +
          computeDistance(dataMinAtom2 , dataMin) +
          computeDistance(dataMaxAtom2 , dataMax) +
          computeDistance(dataSadAtom2 , dataSad);
      }
    }
    int Id3 = std::max_element(distsTo2Pts.begin() , distsTo2Pts.end()) - distsTo2Pts.begin();
    Diagram atom3 = datas[Id3];
    DictDiagrams.push_back(atom3);
    int Id4 = std::max_element(allDists[Id3].begin() , allDists[Id3].end()) - allDists[Id3].begin();
    Diagram atom4 = datas[Id4];
    DictDiagrams.push_back(atom4);
  }

  // if(nbAtoms == 5){
  //   int Id1 = std::max_element(allDistsSummed.begin() , allDistsSummed.end()) - allDistsSummed.begin();
  //   Diagram atom1 = datas[Id1];
  //   DictDiagrams.push_back(atom1);
  //   int Id2 = std::max_element(allDists[Id1].begin() , allDists[Id1]end()) - allDists[Id1].begin();
  //   Diagram atom2 = datas[Id2];
  //   DictDiagrams.push_back(atom2);
  //   std::vector<double> distsTo2Pts(nDiags , 0.);
  //   for(int i = 0 ; i < nDiags ; ++i){
  //     if (i == Id1 || i == Id2){
  //       continue;
  //     } else {
  //       const Diagram &data = datas[i];
  //       distsTo2Pts[i] = computeDistance(atom1 , data) + computeDistance(atom2 , data);
  //     }
  //   }
  //   int Id3 = std::max_element(distsTo2Pts.begin() , distsTo2Pts.end()) - distsTo2Pts.begin();
  //   Diagram atom3 = datas[Id3];
  //   DictDiagrams.push_back(atom3);
  //   int Id4 = std::max_element(allDists[Id3].begin() , allDists[Id3].end()) - allDists[Id3].begin();
  //   Diagram atom4 = datas[Id4];
  //   DictDiagrams.push_back(atom4);
  //   std::vector<double> distsTo4pts(nDiags , 0.);
  //   for(int i = 0 ; i < nDiags ; ++i){
  //     if(i == Id1 || i == Id2 || i == Id3 || i == Id4){
  //       continue;
  //     } else {
  //       Diagram &data = datas[i];
  //       distsTo4Pts[i] = computeDistance(atom1 , data) + computeDistance(atom2 , data) + computeDistance(atom3 , data) + computeDistance(atom4 , data);
  //     }
  //   }
  //   int Id5 = std::max_element(distsTo4Pts.begin() , distsTo4Pts.end()) - distsTo4Pts.begin();
  //   Diagram &atom5 = datas[Id5];
  //   DictDiagrams.push_back(atom5);
  // }
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

double InitFarBorderDict::computeDistance(
  const BidderDiagram<double> &D1,
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
