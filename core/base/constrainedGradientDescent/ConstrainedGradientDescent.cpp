#include <ConstrainedGradientDescent.h>
#include <math.h>

using namespace ttk;

void ConstrainedGradientDescent::executeWeightsProjected(
  std::vector<Matrice> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points) {
  // printf("===========WEIGHT UPDATE=============");
  gradientDescentWeights(hessianList, weights, grad, epoch, nb_points);
  projectionOnSimplex(weights);
}

void ConstrainedGradientDescent::executeAtoms(
  std::vector<Diagram> &DictDiagrams,
  const std::vector<std::vector<MatchingTuple>> &matchings,
  const Diagram &Barycenter,
  const std::vector<Matrice> &gradsLists,
  const int nb_points,
  const std::vector<int> &checkerAtomsExt,
  int epoch) {
  // this->printMsg("==========ATOM UPDATING=============");
  gradientDescentAtoms(DictDiagrams, matchings, Barycenter, gradsLists,
                       nb_points, checkerAtomsExt, epoch);
}

// simple projection on simplex, aka where a vector has positive elements and
// sum to 1.
void ConstrainedGradientDescent::projectionOnSimplex(
  std::vector<double> &weights) {

  int n = weights.size();
  std::vector<double> copy_temp(n);
  for(int i = 0; i < n; ++i) {
    copy_temp[i] = weights[i];
  }
  std::sort(copy_temp.begin(), copy_temp.end(), std::greater<double>());
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

  // double sum = 0.;
  // for(int i = 0; i < n - 1; ++i) {
  //   weights[i] = trunc(weights[i] * 1e6) / 1e6;
  //   sum += weights[i];
  // }
  // weights[n - 1] = 1. - sum;
}

void ConstrainedGradientDescent::gradientDescentWeights(
  std::vector<Matrice> &hessianList,
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points) {

  double mini = *std::min_element(weights.begin(), weights.end());
  int n = weights.size();
  // double norm_grad = 0.;
  // for(int i = 0; i < n; ++i) {
  //   norm_grad += grad[i] * grad[i];
  // }
  double step;
  double L = 0.;
  for(int i = 0 ; i < hessianList.size() ; ++i){
    for(int j = 0 ; j < hessianList[i].size() ; ++j){
      double diag = hessianList[i][j][j];
      L+= 2.*diag;
    }
  }
  step = 1./L;
  std::cout << "STEP = " << step << std::endl;


  // if(nb_points < 100) {
  //   step = mini / (5e3 * (epoch + 1.));
  //   // step = mini / 5e3 * (epoch + 1.);
  // } else if(100 <= nb_points < 700) {
  //   step = std::min(mini, 1.) / (norm_grad); // * pow(2, 10));
  //   // step = std::min(mini , 1.) / norm_grad;
  // } else {
  //   if(epoch < 10) {
  //     // float step = 1 / 2 * *(epoch + 1);
  //     step = std::min(mini, 1.) / pow(2, epoch + 1);
  //   } else {
  //     // float step = 1 / 2 * *11;
  //     step = std::min(mini, 1.) / pow(2, epoch + 1);
  //   }
  // }
  for(int i = 0; i < n; ++i) {
    weights[i] = weights[i] - mini * step * grad[i];
  }
}

void ConstrainedGradientDescent::gradientDescentAtoms(
  std::vector<Diagram> &DictDiagrams,
  const std::vector<std::vector<MatchingTuple>> &matchings,
  const Diagram &Barycenter,
  const std::vector<Matrice> &gradsLists,
  const int nb_points,
  const std::vector<int> &checkerAtomsExt,
  int epoch) {
  // Here vector of diagramTuple because it is not a persistence diagram per
  // say.
  // we get the right pairs to update for each barycenter pair.
  std::vector<Matrice> grad_list(Barycenter.size());

  for(int i = 0; i < grad_list.size(); ++i) {
    grad_list[i].resize(matchings.size());
  }

  std::vector<std::vector<int>> checker(Barycenter.size());
  std::vector<int> tracker(Barycenter.size(), 0);
  std::vector<std::vector<int>> tracker_diagonal(Barycenter.size());
  // std::vector<int> tracker2(Barycenter.size(), 0);
  std::vector<std::vector<int>> tracker_match(Barycenter.size());

  // for(int i = 0; i < matchings.size(); ++i) {
  //  for(int j = 0; j < matchings[i].size(); j++) {
  //    const MatchingTuple &t = matchings[i][j];
  //    const SimplexId Id1 = std::get<0>(t);
  //    const SimplexId Id2 = std::get<1>(t);
  //    DiagramTuple &t2 = DictDiagrams[i][Id1];
  //    grad_list[Id2].push_back(t2);
  //  }
  //}

  for(int i = 0; i < matchings.size(); ++i) {
    for(int j = 0; j < matchings[i].size(); ++j) {
      const MatchingTuple &t = matchings[i][j];
      // Id in atom
      const SimplexId Id1 = std::get<0>(t);
      // Id in barycenter
      const SimplexId Id2 = std::get<1>(t);
      // if(Id2 < 0) {
      if(Id2 < 0 || Id2 >= grad_list.size()
         || Id1 >= static_cast<int>(DictDiagrams[i].size())) {
        continue;
      } else if(Id1 < 0) {
        const DiagramTuple &t3 = Barycenter[Id2];
        std::vector<double> point(2);
        const double birth_barycenter = std::get<6>(t3);
        const double death_barycenter = std::get<10>(t3);
        const double birth_death_atom
          = birth_barycenter + (death_barycenter - birth_barycenter) / 2.;
        point[0] = birth_death_atom;
        point[1] = birth_death_atom;
        // grad_list[Id2].push_back(point);
        grad_list[Id2][i] = std::move(point);
        checker[Id2].push_back(i);
        tracker[Id2] = 1;
        tracker_match[Id2].push_back(-1);
        tracker_diagonal[Id2].push_back(1);

      } else {
        // this->printMsg("====UPDATE GRADLIST========");
        const DiagramTuple &t2 = DictDiagrams[i][Id1];
        std::vector<double> point(2);
        const double birth_atom = std::get<6>(t2);
        const double death_atom = std::get<10>(t2);
        point[0] = birth_atom;
        point[1] = death_atom;
        // grad_list[Id2].push_back(point);
        grad_list[Id2][i] = std::move(point);
        checker[Id2].push_back(i);
        tracker[Id2] = 1;
        tracker_match[Id2].push_back(Id1);
        tracker_diagonal[Id2].push_back(0);
      }
    }
  }

  // printf("=========ATOM GRADIENT STEP==============");
  std::cout << "GRAD_LIST SIZE = LOOP NB: " << grad_list.size() << std::endl;
  for(int i = 0; i < grad_list.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {
      // printf("SAUT!!!!!!!!");
      // std::cout << "SAUT!!!!!!" << std::endl;
      continue;
    } else {
      std::vector<double> pos(grad_list[i].size(), 0.);
      int k = 0;
      // for(int j = 0; j < grad_list[i].size(); ++j) {
      for(int j = 0; j < checker[i].size(); ++j) {
        // DiagramTuple &t = grad_list[i][j];
        // std::vector<double> &t = grad_list[i][j];
        std::vector<double> &t = grad_list[i][checker[i][j]];
        // double birth = std::get<6>(t);
        // double death = std::get<10>(t);
        double birth = t[0];
        double death = t[1];
        pos[j] = death - birth;
        //std::cout << "PERSISTENCE: " << pos[j] << std::endl;
        if(death - birth > 1e-10) {
          k += 1;
        }
      }
      if(k > 0) {

        // printf("==============BOOL VERIFIED==============");
        std::vector<bool> pos2(pos.size(), false);
        // std::vector<double> temp2(pos.size(), 0.);
        std::vector<double> temp2;
        // for(int p = 0; p < pos.size(); ++p) {
        for(int p = 0; p < checker[i].size(); ++p) {
          // DiagramTuple &t = grad_list[i][p];
          // std::vector<double> &t = grad_list[i][p];
          std::vector<double> &t = grad_list[i][checker[i][p]];
          // double birth = std::get<6>(t);
          double birth = t[0];
          pos2[p] = birth == 0;
          if(birth > 0) {
            temp2.push_back(birth);
            // temp2[p] = birth;
          }
        }
        // std::vector<double> temp(pos.size(), 0.);
        std::vector<double> temp;

        for(int p = 0; p < pos.size(); ++p) {
          double val = pos[p];
          if(val > 0.) {
            temp.push_back(val);
            // temp[p] = val;
          }
        }
        // std::cout << "TEMP SIZE " <<temp.size() << std::endl;
        double mini = *std::min_element(temp.begin(), temp.end());

        double step;
        if(temp2.size() == 0) {
          // step = std::min(1., mini) / (1e1 + 1.0 * epoch);
          step = std::min(1. , mini) / (1e1);
        } else {
          double mini2 = *std::min_element(temp2.begin(), temp2.end());
          // double maxi = std::max_element(pos.begin() ; pos.end());
          // step = std::min(std::min(1., mini), mini2) / (1e1 + 1.0 * epoch);
          step = std::min(std::min(1.,mini) , mini2) / (1e1);
        }

        std::cout << "STEP : " << step << std::endl;
        // std::cout << "STEP : " << step << std::endl;
        // double step = 1. / 1e1;
        // for(int p = 0; p < pos.size(); ++p) {
        for(int p = 0; p < checker[i].size(); ++p) {
          if(pos2[p]) {
            // if(pos2[checker[i][p]]) {
            // printf("==========ATOM UPDATING2=============");
            // DiagramTuple &t = grad_list[i][p];
            // std::vector<double> &t = grad_list[i][p];
            std::vector<double> &t = grad_list[i][checker[i][p]];
            // std::vector<double> &t =
            // grad_list[tracker_match[i][p]][checker[i][p]];
            // std::get<10>(t) = std::get<10>(t) - step * gradsLists[i][p][1];
            // t[1] = t[1] - step * gradsLists[i][p][1];
            t[1] = t[1] - step * gradsLists[i][checker[i][p]][1];
          } else if(pos[p] < 1e-17) {
            //} else if(pos[checker[i][p]] < 1e-17) {
            continue;
          } else {
            // printf("==========ATOM UPDATING3=============");
            // DiagramTuple &t = grad_list[i][p];
            // std::get<6>(t) = std::get<6>(t) - step * gradsLists[i][p][0];
            // std::get<10>(t) = std::get<10>(t) - step * gradsLists[i][p][1];
            // std::vector<double> &t = grad_list[i][p];
            // t[0] = t[0] - step * gradsLists[i][p][0];
            // t[1] = t[1] - step * gradsLists[i][p][1];
            std::vector<double> &t0 = grad_list[i][checker[i][p]];
            std::cout << "GRADS" << gradsLists[i][checker[i][p]][0] << ","
                      << gradsLists[i][checker[i][p]][1] << std::endl;
            t0[0] = t0[0] - step * gradsLists[i][checker[i][p]][0];
            t0[1] = t0[1] - step * gradsLists[i][checker[i][p]][1];
          }
        }
        // printf("PASSED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      } else {
        continue;
      }
    }
  }
  // printf("===============PASSED================");
  // for(int i = 0; i < DictDiagrams.size(); ++i) {
  // for(int j = 0; j < nb_points; ++j) {
  // DiagramTuple &t1 = DictDiagrams[i][j];
  // DiagramTuple &t2 = grad_list[j][i];
  // std::vector<double> &t2 = grad_list[j][i];
  // std::get<6>(t1) = std::get<6>(t2);
  // std::get<10>(t1) = std::get<10>(t2);
  // std::get<6>(t1) = t2[0];
  // std::get<10>(t1) = t2[1];
  //}
  //}

  // printf("PASSED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
  for(int i = 0; i < checker.size(); ++i) {
    if(tracker[i] == 0 || checkerAtomsExt[i] == 0) {
      // this->printMsg("SAUT1");
      // printf("SAUT1");
      // continue;
    } else {
      // bool test = false;
      for(int j = 0; j < checker[i].size(); ++j) {
        if(tracker_diagonal[i][j] == 1 || tracker_match[i][j] == -1) {
          // if(test){
          // this->printMsg("SAUT2");
          // printf("SAUT2");
          // continue;
        } else {
          std::vector<double> &t2 = grad_list[i][checker[i][j]];

          // if(t2[1] - t2[0] < 1e-17) {
          //   DictDiagrams[checker[i][j]].erase(
          //     DictDiagrams[checker[i][j]].begin() + tracker_match[i][j]);
          // } else {
          //   DiagramTuple &t1 =
          //   DictDiagrams[checker[i][j]][tracker_match[i][j]];
          //   // printf("ATOM" + std::to_string(checker[i][j]) " , PAIR " +
          //   // std::to_string(tracker_match[i][j]));
          //   // std::cout << "ATOM " << checker[i][j] << " , SIZE"
          //   //           << DictDiagrams[checker[i][j]].size() << std::endl;
          //   // std::cout << "ATOM " << checker[i][j] << " , PAIR"
          //   //           << tracker_match[i][j] << std::endl;
          //   std::get<6>(t1) = t2[0];
          //   std::get<10>(t1) = t2[1];
          // }
          if(t2[1] > t2[0]) {
            DiagramTuple &t1 = DictDiagrams[checker[i][j]][tracker_match[i][j]];
            // printf("ATOM" + std::to_string(checker[i][j]) " , PAIR " +
            // std::to_string(tracker_match[i][j]));
            // std::cout << "ATOM " << checker[i][j] << " , SIZE"
            //           << DictDiagrams[checker[i][j]].size() << std::endl;
            // std::cout << "ATOM " << checker[i][j] << " , PAIR"
            //           << tracker_match[i][j] << std::endl;
            std::get<6>(t1) = t2[0];
            std::get<10>(t1) = t2[1];
          } else {
            DictDiagrams[checker[i][j]].erase(
              DictDiagrams[checker[i][j]].begin() + tracker_match[i][j]);
          }
        }
      }
    }
  }
}
