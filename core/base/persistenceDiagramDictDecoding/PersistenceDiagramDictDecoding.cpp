#include <PersistenceDiagramDictDecoding.h>
#include <algorithm>
using namespace ttk;

void PersistenceDiagramDictDecoding::execute(
  std::vector<Diagram> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  std::vector<Diagram> &Barycenters) const {

  std::vector<std::vector<VectorMatchingTuple>> AllMatchingsAtoms(
    Barycenters.size());
  for(size_t i = 0; i < Barycenters.size(); ++i) {
    auto &barycenter = Barycenters[i];
    auto &weight = vectorWeights[i];
    auto &matchings = AllMatchingsAtoms[i];
    computeWeightedBarycenter(dictDiagrams, weight, barycenter, matchings);
  }
}
