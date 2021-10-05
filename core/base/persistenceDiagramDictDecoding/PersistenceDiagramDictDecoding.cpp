#include <PersistenceDiagramDictDecoding.h>
#include <algorithm>
using namespace ttk;

void PersistenceDiagramDictDecoding::execute(
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  std::vector<Diagram> &Barycenters) const {

  std::vector<std::vector<std::vector<MatchingTuple>>> AllMatchingsAtoms(
    Barycenters.size());
  for(int i = 0; i < Barycenters.size(); ++i) {
    Diagram &barycenter = Barycenters[i];
    std::vector<double> &weight = vectorWeights[i];
    //  std::cout << "Poids: " << weight[0] << weight[1] << weight[2]
    //          << std::endl;
    // std::cout << "================================================="
    //          << std::endl;
    std::vector<std::vector<MatchingTuple>> &matchings = AllMatchingsAtoms[i];
    computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
  }
}
