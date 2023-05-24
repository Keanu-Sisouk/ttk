#include "PersistenceDiagramUtils.h"
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


void ttk::PersistenceDiagramDictDecoding::computeAtomsCoordinates(
  const std::vector<ttk::DiagramType> &atoms,
  const std::vector<std::vector<double>> &vectorWeights,
  std::vector<std::pair<double, double>> &coords,
  std::vector<std::pair<double, double>> &true_coords,
  const double spacing,
  const double max_persistence,
  const size_t nAtoms) const{
  
  if(nAtoms == 2) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(2);
    const auto distMatrix = MatrixCalculator.execute(atoms, nInputs);
    coords[0].first = 0.;
    true_coords[0].first = 0.;
    coords[0].second = 0.;
    true_coords[0].first = 0.;
    coords[1].first = spacing * distMatrix[0][1];
    true_coords[1].first = distMatrix[0][1];
    
  } else if(nAtoms == 3) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix
      = MatrixCalculator.execute(atoms, nInputs);
    coords[0].first = 0.;
    true_coords[0].first = 0.;
    coords[0].second = 0.;
    true_coords[0].second = 0.;
    coords[1].first = spacing * distMatrix[0][1];
    true_coords[1].first = distMatrix[0][1];
    coords[1].second = 0.;
    true_coords[0].second = 0.;
    double distOpposed = distMatrix[2][1];
    double firstDist = distMatrix[0][1];
    double distAdja = distMatrix[0][2];
    double alpha = std::acos(
      (distOpposed * distOpposed - firstDist * firstDist - distAdja * distAdja)
      / (-2. * firstDist * distAdja));
    coords[2].first = spacing * distAdja * std::cos(alpha);
    true_coords[2].first = distAdja * std::cos(alpha);
    coords[2].second = spacing * distAdja * std::sin(alpha);
    true_coords[2].second = distAdja * std::sin(alpha);

  } else {
    for(size_t i = 0; i < nAtoms; ++i) {
      const auto angle
        = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nAtoms);
      true_coords[i].first = max_persistence * std::cos(angle);
      true_coords[i].second = max_persistence * std::sin(angle);
    }
  }
}