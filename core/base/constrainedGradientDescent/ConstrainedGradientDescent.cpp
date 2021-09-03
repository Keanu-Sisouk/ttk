#include <ConstrainedGradientDescent.h>

using namespace ttk;

void ConstrainedGradientDescent::executeWeightsProjected(
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points) {
  gradientDescentWeights(weights, grad, epoch, nb_points);
  projectionOnSimplex(weights);
}

void ConstrainedGradientDescent::executeAtoms(
  std::vector<Diagram> &DictDiagrams,
  const std::vector<std::vector<MatchingTuple>> &matchings,
  const Diagram &Barycenter,
  const std::vector<Matrice> &gradsLists,
  const int nb_points) {
  gradientDescentAtoms(
    DictDiagrams, matchings, Barycenter, gradsLists, nb_points);
}


//simple projection on simplex, aka where a vector has positive elements and 
//sum to 1.
void ConstrainedGradientDescent::projectionOnSimplex(
  std::vector<double> &weights) {
  int n = weights.size();
  std::sort(weights.begin(), weights.end(), std::greater<double>());
  // std::vector<double> u = std::sort(weights.begin(), weights.end(),
  // std::greater<double>());
  double K = 1.;
  double somme_u = weights[0];
  double theta = (somme_u - 1.) / K;
  while(K < n && (somme_u + weights[K] - 1.) / (K + 1.) < weights[K]) {
    somme_u += weights[K];
    K += 1;
    theta = (somme_u - 1.) / K;
  }
  for(int i = 0; i < n; ++i) {
    weights[i] = std::max(weights[i] - theta, 0.);
  }
}

void ConstrainedGradientDescent::gradientDescentWeights(
  std::vector<double> &weights,
  const std::vector<double> &grad,
  const int epoch,
  const int nb_points) {
  double mini = *std::min_element(weights.begin(), weights.end());
  int n = weights.size();
  double norm_grad = 0;
  for(int i = 0; i < n; ++i) {
    norm_grad += weights[i] * weights[i];
  }
  float step = 0;
  if(nb_points < 100) {
    step = mini / (5e3 * (epoch + 1));
  } else if(100 <= nb_points < 700) {
    step = std::min(mini, 1.) / norm_grad;
  } else {
    if(epoch < 10) {
      // float step = 1 / 2 * *(epoch + 1);
      step = 1 / pow(2, epoch + 1);
    } else {
      // float step = 1 / 2 * *11;
      step = 1 / pow(2, 11);
    }
  }
  for(int i = 0; i < n; ++i) {
    weights[i] = weights[i] - step * grad[i];
  }
}

void ConstrainedGradientDescent::gradientDescentAtoms(
  std::vector<Diagram> &DictDiagrams,
  const std::vector<std::vector<MatchingTuple>> &matchings,
  const Diagram &Barycenter,
  const std::vector<Matrice> &gradsLists,
  const int nb_points) {
  // Here vector of diagramTuple because it is not a persistence diagram per
  // say.
  //we get the right paris to update for each barycenter pair.
  std::vector<std::vector<DiagramTuple>> grad_list(Barycenter.size());
  for(int i = 0; i < matchings.size(); ++i) {
    for(int j = 0; j < matchings[i].size(); j++) {
      const MatchingTuple &t = matchings[i][j];
      const SimplexId Id1 = std::get<0>(t);
      const SimplexId Id2 = std::get<1>(t);
      DiagramTuple &t2 = DictDiagrams[i][Id1];
      grad_list[Id2].push_back(t2);
    }
  }
  for(int i = 0; i < grad_list.size(); ++i) {
    std::vector<double> pos(grad_list[i].size(), 1);
    int k = 0;
    for(int j = 0; j < grad_list[i].size(); ++j) {
      DiagramTuple &t = grad_list[i][j];
      double birth = std::get<6>(t);
      double death = std::get<10>(t);
      pos[j] = death - birth;
      if(death - birth < 1e-10) {
        k += 1;
      }
    }
    if(k > 0) {
      std::vector<bool> pos2(pos.size(), false);
      std::vector<float> temp2;
      for(int k = 0; k < pos.size(); ++k) {
        DiagramTuple &t = grad_list[i][k];
        double birth = std::get<6>(t);
        pos2[k] = birth > 0;
        if(birth > 0) {
          temp2.push_back(birth);
        }
      }
      std::vector<double> temp;
      for(int k = 0; k < pos.size(); ++k) {
        double val = pos[k];
        if(val > 0) {
          temp.push_back(val);
        }
      }
      double mini = *std::min_element(temp.begin(), temp.end());
      double mini2 = *std::min_element(temp2.begin(), temp2.end());
      // double maxi = std::max_element(pos.begin() ; pos.end());
      double step = std::min(std::min(1., mini), mini2) / 1e1;
      for(int k = 0; k < pos.size(); ++k) {
        if(pos2[k]) {
          DiagramTuple &t = grad_list[i][k];
          std::get<10>(t) = std::get<10>(t) - step * gradsLists[i][k][1];
        } else if(pos[k] < 1e-17) { ;
        } else {
          DiagramTuple &t = grad_list[i][k];
          std::get<6>(t) = std::get<6>(t) - step * gradsLists[i][k][0];
          std::get<10>(t) = std::get<10>(t) - step * gradsLists[i][k][1];
        }
      }
    }
  }
}
