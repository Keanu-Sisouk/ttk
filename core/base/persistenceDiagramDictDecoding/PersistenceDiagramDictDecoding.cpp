#include "DimensionReduction.h"
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
  std::vector<ttk::DiagramType> &atoms,
  const std::vector<std::vector<double>> &vectorWeights,
  std::vector<std::pair<double, double>> &coords,
  std::vector<std::pair<double, double>> &true_coords,
  std::vector<double> &xVector,
  std::vector<double> &yVector,
  const double spacing,
  const double max_persistence,
  const size_t nAtoms) const {

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
    true_coords[1].second = 0.;
    
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
    true_coords[1].second = 0.;
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

  } else if(nAtoms == 4){
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix
      = MatrixCalculator.execute(atoms, nInputs);
    // coords[0].first = 0.;
    // true_coords[0].first = 0.;
    // coords[0].second = 0.;
    // true_coords[0].second = 0.;
    // coords[1].first = spacing * distMatrix[0][1];
    // true_coords[1].first = distMatrix[0][1];
    // coords[1].second = 0.;
    // true_coords[0].second = 0.;
    // double distOpposed = distMatrix[2][1];
    // double firstDist = distMatrix[0][1];
    // double distAdja = distMatrix[0][2];
    // double alpha = std::acos(
    //   (distOpposed * distOpposed - firstDist * firstDist - distAdja * distAdja)
    //   / (-2. * firstDist * distAdja));
    // coords[2].first = spacing * distAdja * std::cos(alpha);
    // true_coords[2].first = distAdja * std::cos(alpha);
    // coords[2].second = spacing * distAdja * std::sin(alpha);
    // true_coords[2].second = distAdja * std::sin(alpha);





    ttk::DimensionReduction DimProjector;
    DimProjector.setIsInputDistanceMatrix(true);
    // ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    // std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    // std::vector<std::vector<double>> distMatrix
      // = MatrixCalculator.execute(atoms, nInputs);
    int nRow = distMatrix.size();
    std::vector<double> matrixForProjector;
    for(int i = 0; i < nRow; ++i) {
      for(int j = 0; j < nRow; ++j) {
        matrixForProjector.push_back(distMatrix[j][i]);
      }
    }
    std::vector<std::vector<double>> coordsAtom;
    DimProjector.execute(coordsAtom, matrixForProjector, nRow, nRow);

    for(size_t i = 0; i < 2; ++i) {
      for(size_t j = 0; j < nAtoms; ++j) {
        if(i == 0) {
          true_coords[j].first = coordsAtom[0][j];
        } else {
          true_coords[j].second = coordsAtom[1][j];
        }
      }
    }
  
  
  
  } else {
    // std::vector<ttk::DiagramType> dictDiagrams;
    // std::vector<ttk::DiagramType> intermediateAtoms;
    // std::vector<double> lossTab;
    // std::vector<double> trueLossTab;
    // std::vector<double> timers;
    // std::vector<std::vector<double>> allLosses(nAtoms);
    // const int seed = 0;
    // const int m = 3;
    // std::vector<std::vector<double>> tempWeights(nAtoms);
    // for(size_t i = 0; i < tempWeights.size(); ++i) {
    //   std::vector<double> weights(m, 1. / (m * 1.));
    //   tempWeights[i] = std::move(weights);
    // }
    // ttk::PersistenceDiagramDictEncoding DictionaryEncoder;
    // DictionaryEncoder.setUseDimReduct(false);
    // DictionaryEncoder.setUseProgApproach(true);
    // DictionaryEncoder.execute(atoms, atoms, dictDiagrams, tempWeights, seed,
    // m,
    //                           lossTab, timers, trueLossTab, allLosses, 0.);
    // std::vector<std::pair<double, double>> tempCoords(3);
    // std::vector<std::pair<double, double>> temp_true_coords(3);
    // ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    // std::array<size_t, 2> nInputs{3, 0};
    // MatrixCalculator.setDos(true, true, true);
    // MatrixCalculator.setThreadNumber(3);
    // std::vector<std::vector<double>> distMatrix
    //   = MatrixCalculator.execute(dictDiagrams, nInputs);
    // tempCoords[0].first = 0.;
    // temp_true_coords[0].first = 0.;
    // tempCoords[0].second = 0.;
    // temp_true_coords[0].second = 0.;
    // tempCoords[1].first = spacing * distMatrix[0][1];
    // temp_true_coords[1].first = distMatrix[0][1];
    // tempCoords[1].second = 0.;
    // temp_true_coords[0].second = 0.;
    // double distOpposed = distMatrix[2][1];
    // double firstDist = distMatrix[0][1];
    // double distAdja = distMatrix[0][2];
    // double alpha = std::acos(
    //   (distOpposed * distOpposed - firstDist * firstDist - distAdja *
    //   distAdja) / (-2. * firstDist * distAdja));
    // tempCoords[2].first = spacing * distAdja * std::cos(alpha);
    // temp_true_coords[2].first = distAdja * std::cos(alpha);
    // tempCoords[2].second = spacing * distAdja * std::sin(alpha);
    // temp_true_coords[2].second = distAdja * std::sin(alpha);  
    // for(int i = 0; i < 2; ++i) {
    //   for(size_t j = 0; j < nAtoms; ++j) {
    //     double temp = 0.;
    //     for(int iAtom = 0; iAtom < 3; ++iAtom) {
    //       if(i == 0) {
    //         temp += tempWeights[j][iAtom] * temp_true_coords[iAtom].first;
    //         true_coords[j].first = temp;
    //       } else {
    //         temp += tempWeights[j][iAtom] * temp_true_coords[iAtom].second;
    //         true_coords[j].second = temp;
    //       }
    //     }
    //   }
    // }


    ttk::DimensionReduction DimProjector;
    DimProjector.setIsInputDistanceMatrix(true);
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix
      = MatrixCalculator.execute(atoms, nInputs);
    int nRow = distMatrix.size();
    std::vector<double> matrixForProjector;
    for(int i = 0; i < nRow; ++i) {
      for(int j = 0; j < nRow; ++j) {
        matrixForProjector.push_back(distMatrix[j][i]);
      }
    }
    std::vector<std::vector<double>> coordsAtom;
    DimProjector.execute(coordsAtom, matrixForProjector, nRow, nRow);

    for(size_t i = 0; i < 2; ++i) {
      for(size_t j = 0; j < nAtoms; ++j) {
        if(i == 0) {
          true_coords[j].first = coordsAtom[0][j];
        } else {
          true_coords[j].second = coordsAtom[1][j];
        }
      }
    }
  }
  size_t nDiags = vectorWeights.size();

  for(int i = 0; i < 2; ++i) {
    for(size_t j = 0; j < nDiags; ++j) {
      double temp = 0.;
      for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        if(i == 0) {
          temp += vectorWeights[j][iAtom] * true_coords[iAtom].first;
        } else {
          temp += vectorWeights[j][iAtom] * true_coords[iAtom].second;
        }
      }
      if(i == 0) {
        xVector[j] = temp;
      } else {
        yVector[j] = temp;
      }
    }
  }
}