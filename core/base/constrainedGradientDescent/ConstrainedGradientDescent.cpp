#include <ConstrainedGradientDescent.h>
#include <cmath>
#include <csignal>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif // TTK_ENABLE_EIGEN

using namespace ttk;

void ConstrainedGradientDescent::executeWeightsProjected(
  std::vector<Matrix> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points,
  bool MaxEigenValue) {
  gradientDescentWeights(
    hessianList, weights, grad, epoch, nb_points, MaxEigenValue);
  projectionOnSimplex(weights);
}

void ConstrainedGradientDescent::executeAtoms(
  std::vector<ttk::DiagramType> &DictDiagrams,
  const std::vector<std::vector<ttk::MatchingType>> &matchings,
  const ttk::DiagramType &Barycenter,
  const std::vector<Matrix> &gradsLists,
  const int nb_points,
  const std::vector<int> &checkerAtomsExt,
  int epoch,
  std::vector<std::vector<int>> &projForDiag,
  ttk::DiagramType &featuresToAdd,
  std::vector<std::array<double, 2>> &projLocations,
  std::vector<std::vector<double>> &vectorForProjContrib,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  ttk::DiagramType &infoToAdd) {
  gradientDescentAtoms(DictDiagrams, matchings, Barycenter, gradsLists,
                       nb_points, checkerAtomsExt, epoch, projForDiag,
                       featuresToAdd, projLocations, vectorForProjContrib,
                       pairToAddGradList, infoToAdd);
}

// simple projection on simplex, aka where a vector has positive elements and
// sum to 1.
void ConstrainedGradientDescent::projectionOnSimplex(
  std::vector<double> &weights) {

  int n = weights.size();
  std::vector<double> copy_temp = weights;
  std::sort(copy_temp.rbegin(), copy_temp.rend());
  // std::vector<double> u = std::sort(weights.begin(), weights.end(),
  // std::greater<double>());
  double K = 1.;
  double somme_u = copy_temp[0];
  double theta = (somme_u - 1.) / K;
  while(K < n && (somme_u + copy_temp[K] - 1.) / (K + 1.) < copy_temp[K]) {
    somme_u += copy_temp[K];
    K += 1.;
    theta = (somme_u - 1.) / K;
  }
  for(int i = 0; i < n; ++i) {
    weights[i] = std::max(weights[i] - theta, 0.);
  }

  double sum = 0.;
  for(int i = 0; i < n - 1; ++i) {
    weights[i] = trunc(weights[i] * 1e8) / 1e8;
    sum += weights[i];
  }
  weights[n - 1] = 1. - sum;
}

void ConstrainedGradientDescent::gradientDescentWeights(
  std::vector<Matrix> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points,
  bool MaxEigenValue) {

  double mini = *std::min_element(weights.begin(), weights.end());
  int n = weights.size();
  double stepWeight;
  double L = 0.;
  // std::cout << "STEP = " << step << std::endl;
#ifndef TTK_ENABLE_EIGEN
  MaxEigenValue = false;
#endif // TTK_ENABLE_EIGEN
  if(MaxEigenValue) {
#ifdef TTK_ENABLE_EIGEN
    for(size_t i = 0; i < hessianList.size(); ++i) {
      auto &hessian = hessianList[i];
      int m = hessian.size();
      Eigen::MatrixXd H(m, m);
      for(size_t j = 0; j < hessian.size(); ++j) {
        for(size_t k = 0; k < hessian.size(); ++k) {
          H(j, k) = hessian[j][k];
        }
      }
      Eigen::EigenSolver<Eigen::MatrixXd> es;
      es.compute(H, false);
      Eigen::VectorXcd eigvals = es.eigenvalues();
      L += eigvals.lpNorm<Eigen::Infinity>();
    }
#endif // TTK_ENABLE_EIGEN
  } else {
    for(size_t i = 0; i < hessianList.size(); ++i) {
      auto &hessian = hessianList[i];
      for(size_t k = 0; k < hessian.size(); ++k) {
        double diag = hessian[k][k];
        L += 1. * diag;
      }
    }
  }

  if(L > 0){
    stepWeight = 1. / L;
  } else {
    stepWeight = 0;
  }
  // std::cout << "STEP" << step << std::endl;

  for(int i = 0; i < n; ++i) {
    weights[i] = weights[i] - stepWeight * grad[i];
    // std::cout << "GRAD: " + std::to_string(grad[i]) << std::endl;
  }
}

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// METTRE TIMER POUR VOIR QUOI PARALELLISER
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

void ConstrainedGradientDescent::gradientDescentAtoms(
  std::vector<ttk::DiagramType> &DictDiagrams,
  const std::vector<std::vector<ttk::MatchingType>> &matchings,
  const ttk::DiagramType &Barycenter,
  const std::vector<Matrix> &gradsLists,
  const int nb_points,
  const std::vector<int> &checkerAtomsExt,
  int epoch,
  std::vector<std::vector<int>> &projForDiag,
  ttk::DiagramType &featuresToAdd,
  std::vector<std::array<double, 2>> &projLocations,
  std::vector<std::vector<double>> &vectorForProjContrib,
  std::vector<std::vector<std::array<double, 2>>> &pairToAddGradList,
  ttk::DiagramType &infoToAdd) {

  // Here vector of diagramTuple because it is not a persistence diagram per
  // say.
  // we get the right pairs to update for each barycenter pair.

  std::vector<double> miniBirth(matchings.size());
  for(size_t i = 0; i < matchings.size(); ++i) {
    auto &t = DictDiagrams[i][0];
    miniBirth[i] = t.birth.sfValue;
  }

  std::vector<std::vector<std::array<double, 2>>> grad_list(Barycenter.size());
  std::vector<std::vector<double>> projectionsBuffer(Barycenter.size());

  for(size_t i = 0; i < Barycenter.size(); ++i) {
    projectionsBuffer[i].resize(matchings.size());
  }

  for(size_t i = 0; i < grad_list.size(); ++i) {
    grad_list[i].resize(matchings.size());
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  for(size_t i = 0; i < checker.size(); ++i) {
    checker[i].resize(matchings.size());
  }
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<std::vector<int>> tracker_diagonal(Barycenter.size());
  std::vector<std::vector<int>> tracker_match(Barycenter.size());
  for(size_t i = 0; i < tracker_diagonal.size(); ++i) {
    tracker_diagonal[i].resize(matchings.size());
  }
  for(size_t i = 0; i < tracker_match.size(); ++i) {
    tracker_match[i].resize(matchings.size());
  }

  for(size_t i = 0; i < matchings.size(); ++i) {
    for(size_t j = 0; j < matchings[i].size(); ++j) {
      const auto &t = matchings[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      if(Id2 < 0 || static_cast<SimplexId>(grad_list.size()) <= Id2
         || static_cast<SimplexId>(DictDiagrams[i].size()) <= Id1) {
        continue;
      } else {
        if(Id1 < 0) {
          const auto &t3 = Barycenter[Id2];
          auto &point = grad_list[Id2][i];
          const double birth_barycenter = t3.birth.sfValue;
          const double death_barycenter = t3.death.sfValue;
          const double birth_death_atom
            = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
          point[0] = birth_death_atom;
          point[1] = birth_death_atom;

          // checker[Id2].push_back(i);
          checker[Id2][i] = i;
          if(checker[Id2].size() > matchings.size()) {
            std::raise(SIGINT);
          }
          tracker[Id2] = 1;
          // tracker_match[Id2].push_back(Id1);
          tracker_match[Id2][i] = Id1;
          if(static_cast<SimplexId>(DictDiagrams[i].size()) <= Id1) {
            std::cout << "ID1: " << Id1 << std::endl;
          }
          tracker_diagonal[Id2][i] = 1;
          projectionsBuffer[Id2][i] = birth_death_atom;

        } else {
          // this->printMsg("====UPDATE GRADLIST========");
          const auto &t2 = DictDiagrams[i][Id1];
          auto &point = grad_list[Id2][i];
          const double birth_atom = t2.birth.sfValue;
          const double death_atom = t2.death.sfValue;
          point[0] = birth_atom;
          point[1] = death_atom;
          // checker[Id2].push_back(i);
          checker[Id2][i] = i;
          if(checker[Id2].size() > matchings.size()) {
            std::raise(SIGINT);
          }
          tracker[Id2] = 1;
          // tracker_match[Id2].push_back(Id1);

          tracker_match[Id2][i] = Id1;

          tracker_diagonal[Id2][i] = 0;
          projectionsBuffer[Id2][i] = 0.;
        }
      }
    }
  }

  grad_list.insert(
    grad_list.end(), pairToAddGradList.begin(), pairToAddGradList.end());
  for(size_t j = 0; j < pairToAddGradList.size(); ++j) {
    std::vector<int> temp1(matchings.size(), 1);
    std::vector<int> temp2(matchings.size(), -1);
    std::vector<int> temp3(matchings.size());
    for(size_t l = 0; l < matchings.size(); ++l) {
      temp3[l] = static_cast<int>(l);
    }
    tracker_diagonal.push_back(temp1);
    tracker_match.push_back(temp1);
    checker.push_back(temp3);
    tracker.push_back(1);
  }

  // #ifdef TTK_ENABLE_OPENMP
  // #pragma omp parallel for num_threads(threadNumber_)
  // #endif // TTK_ENABLE_OPENMP

  for(size_t i = 0; i < grad_list.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {

      continue;
    } else {

      std::vector<double> pos(grad_list[i].size(), 0.);
      int k = 0;
      // int k = 1;
      // for(int j = 0; j < grad_list[i].size(); ++j) {
      for(size_t j = 0; j < checker[i].size(); ++j) {
        // DiagramTuple &t = grad_list[i][j];
        // std::vector<double> &t = grad_list[i][j];
        auto &t = grad_list[i][checker[i][j]];
        // double birth = std::get<6>(t);
        // double death = std::get<10>(t);
        double birth = t[0];
        double death = t[1];
        pos[j] = death - birth;
      }
      // if(true) {
      //  printf("==============BOOL VERIFIED==============");
      std::vector<bool> pos2(pos.size(), false);
      // std::vector<double> temp2(pos.size(), 0.);
      std::vector<double> temp2;
      // for(int p = 0; p < pos.size(); ++p) {
      for(size_t p = 0; p < checker[i].size(); ++p) {
        // DiagramTuple &t = grad_list[i][p];
        // std::vector<double> &t = grad_list[i][p];
        auto &t = grad_list[i][checker[i][p]];
        // double birth = std::get<6>(t);
        double birth = t[0];
        pos2[p] = birth == 0.;
        if(birth > 0.) {
          temp2.push_back(birth);
          // temp2[p] = birth;
        }
      }
      // std::vector<double> temp(pos.size(), 0.);
      /* std::vector<double> temp; */

      /* for(size_t p = 0; p < pos.size(); ++p) { */
      /*   double val = pos[p]; */
      /*   if(val > 0.) { */
      /*     temp.push_back(val); */
      /*     // temp[p] = val; */
      /*   } */
      /* } */
      // std::cout << "TEMP SIZE " <<temp.size() << std::endl;
      // double mini = *std::min_element(temp.begin(), temp.end());

      // double step;
      // double factEquiv = static_cast<double>(DictDiagrams.size());
      // double factEquiv = 1.;
      // step = 1. / (sqrt(factEquiv) * 1e1);
      // step = 1. / (2. * 2. * factEquiv);
      // double factEquiv = DictDiagrams.size();

      // setStep(factEquiv);

      for(size_t p = 0; p < checker[i].size(); ++p) {
        if(pos2[p]) {

          auto &t = grad_list[i][checker[i][p]];

          t[1] = t[1] - (this->stepAtom) * gradsLists[i][checker[i][p]][1];
        } else if(pos[p] < 1e-7) {

          // auto &t0 = grad_list[i][checker[i][p]];

          // t0[0] = t0[0] - step * gradsLists[i][checker[i][p]][0];
          // t0[1] = t0[1] - step * gradsLists[i][checker[i][p]][1];
          continue;
        } else {

          auto &t0 = grad_list[i][checker[i][p]];

          t0[0] = t0[0] - (this->stepAtom) * gradsLists[i][checker[i][p]][0];
          t0[1] = t0[1] - (this->stepAtom) * gradsLists[i][checker[i][p]][1];

          if(t0[0] > t0[1]) {
            t0[1] = t0[0];
          }

          // if (i == 0 && checker[i][p] == 0) std::cout << "PAIRE GLOBALE
          // AFTER: " << t0[0] << " and " << t0[1] << std::endl;
        }
      }
      // printf("PASSED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      //} else {
      //  continue;
      //}
    }
  }

  // printf("PASSED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
  int count = 0;
  for(size_t i = 0; i < checker.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {
      // this->printMsg("SAUT1");
      // printf("SAUT1");
      continue;
    } else {
      if(i < Barycenter.size()) {
        const auto &infos = Barycenter[i];
        for(size_t j = 0; j < checker[i].size(); ++j) {
          auto &tracker_temp = tracker_match[i][j];
          if(tracker_diagonal[i][j] == 1 || tracker_temp == -1
             || tracker_temp > 10000) {
            auto &index = checker[i][j];
            auto &t2 = grad_list[i][index];
            // if(t2[1] - t2[0] < 1e-7) {
            // continue;
            //} else {
            std::vector<int> projAndIndex(DictDiagrams.size() + 1);
            // double proj = projectionsBuffer[i][index];
            int atomIndex = static_cast<int>(index);
            // projAndIndex[0] = proj;
            // projAndIndex[1] = atomIndex;
            for(size_t m = 0; m < DictDiagrams.size(); ++m) {
              projAndIndex[m] = tracker_match[i][m];
            }
            // projAndIndex[DictDiagrams.size()] = atomIndex;
            // projForDiag.push_back(projAndIndex);
            // featuresToAdd.push_back(infos);
            // projLocations.push_back(grad_list[i][index]);
            // vectorForProjContrib.push_back(gradsLists[i][index]);
            //  DictDiagrams[index].push_back(newPair);
            // }

          } else {
            auto &index = checker[i][j];
            auto &t2 = grad_list[i][index];

            auto &t1 = DictDiagrams[index][tracker_temp];

            if(t2[0] < miniBirth[index]) {
              t1.birth.sfValue = miniBirth[index];
            } else {
              t1.birth.sfValue = t2[0];
            }
            if(t2[1] < t2[0]) {
              // count +=1;

              t1.birth.sfValue = t2[0];
              t1.death.sfValue = t2[0];
              //t1.persistence = t1.death.sfValue - t1.birth.sfValue;
              // std::get<10>(t1) = t2[0];
              // std::cout << "Under diag" << std::endl;
              //continue;
            } else {
              t1.death.sfValue = t2[1];
            }
          }
        }
      } else {
        const auto &infos = infoToAdd[static_cast<int>(i)
                                      - static_cast<int>(Barycenter.size())];

        for(size_t j = 0; j < checker[i].size(); ++j) {
          auto &tracker_temp = tracker_match[i][j];
          if(tracker_diagonal[i][j] == 1 || tracker_temp == -1
             || tracker_temp > 10000) {
            auto &index = checker[i][j];
            auto &t2 = grad_list[i][index];
            // if(t2[1] - t2[0] < 1e-7) {
            // continue;
            //} else {
            std::vector<int> projAndIndex(DictDiagrams.size() + 1);
            // double proj = projectionsBuffer[i][index];
            int atomIndex = static_cast<int>(index);
            // projAndIndex[0] = proj;
            // projAndIndex[1] = atomIndex;
            for(size_t m = 0; m < DictDiagrams.size(); ++m) {
              projAndIndex[m] = tracker_match[i][m];
            }
            projAndIndex[DictDiagrams.size()] = atomIndex;
            projForDiag.push_back(projAndIndex);
            featuresToAdd.push_back(infos);
            projLocations.push_back(grad_list[i][index]);
            vectorForProjContrib.push_back(gradsLists[i][index]);
            // DictDiagrams[index].push_back(newPair);
            //}

          } else {
            auto &index = checker[i][j];
            auto &t2 = grad_list[i][index];
            auto &t1 = DictDiagrams[index][tracker_temp];
            if(t2[0] < miniBirth[index]) {
              t1.birth.sfValue = miniBirth[index];
            } else {
              t1.birth.sfValue = t2[0];
            }
            if(t2[1] < t2[0]) {
              // count+=1;
              t1.birth.sfValue = t2[0];
              t1.death.sfValue = t2[0];
              //t1.persistence = t1.death.sfValue - t1.birth.sfValue;

              // std::cout << "Under diag" << std::endl;
              //continue;
            } else {
              t1.death.sfValue = t2[1];
            }

          }
        }
      }
    }
  }
  // std::cout << "COUNT OF UNDER DIAG " << count << std::endl;
}


void ConstrainedGradientDescent::setStep(double factEquiv){
  this->stepAtom = 1. / (2. * 2. * factEquiv);
}

void ConstrainedGradientDescent::reduceStep(){
  this->stepAtom = this->stepAtom/2.;
}