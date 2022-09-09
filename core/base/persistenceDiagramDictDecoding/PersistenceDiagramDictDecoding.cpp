#include <PersistenceDiagramClustering.h>
#include <PersistenceDiagramDictDecoding.h>

void ttk::PersistenceDiagramDictDecoding::execute(
  std::vector<ttk::DiagramType> &dictDiagrams,
  std::vector<std::vector<double>> &vectorWeights,
  std::vector<ttk::DiagramType> &Barycenters) const {

  Timer tm{};

  std::vector<std::vector<std::vector<MatchingType>>> AllMatchingsAtoms(
    Barycenters.size());

  for(size_t i = 0; i < Barycenters.size(); ++i) {
    auto &barycenter = Barycenters[i];
    auto &weight = vectorWeights[i];
    auto &matchings = AllMatchingsAtoms[i];
    computeWeightedBarycenter(
      dictDiagrams, weight, barycenter, matchings, *this, ProgBarycenter);
  }

  this->printMsg(
    "Computed barycenters", 1.0, tm.getElapsedTime(), this->threadNumber_);
}
