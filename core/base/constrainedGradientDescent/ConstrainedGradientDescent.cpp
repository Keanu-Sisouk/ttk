#include <ConstrainedGradientDescent.h>

using namespace ttk;




void ConstrainedGradientDescent::executeWeightsProjected(std::vector<double> &weights , const std::vector<double> &grad , const int epoch, const int nb_points){
  gradientDescentWeights(weights , grad , epoch , nb_points);
  projectionOnSimplex(weights);
}


void ConstrainedGradientDescent::projectionOnSimplex(std::vector<double> &weights){
  int n = weights.size();
  std::vector<double> u = sort(weights , weights + n , greater<double>());
  double K = 1.;
  double somme_u = u[0];
  double theta = (somme_u - 1.)/K;
  while(K < n && (somme_u + u[K] -1)/(K+1) < u[K]){
    somme_u +=u[K];
    K+=1;
    theta = (somme_u - a)/K;
  }
  for (int i = 0 ; i < n ; ++i){
    weights[i] = std::max(weights[i] - theta , 0);
  }
}


void ConstrainedGradientDescent::gradientDescentWeights(std::vector<double> &weights,const std::vector<double> &grad, const int epoch,const int nb_points){
  double mini = std::min_element(weights.begin() , weights.end());
  int n = weights.size();
  double norm_grad = 0;
  for(int i = 0 ; i < n ; ++i){
    norm_grad+= weights[i]**2;
  }
  if(nb_points < 100){
    float step = mini/(5e3*(epoch+1));
  }else if(100 <= nb_points < 700 ){
    float step = std::min(mini , 1.)/norm_grad;
  }else{
    if(epoch<10){
      float step = 1/2**(epoch+1);
    }else{
      float step = 1/2**11;
    }
  }
  for(int i = 0 ; i < n ; ++i){
    weights[i] = weights[i] - step*grad[i];
  }
}


void ConstrainedGradientDescent::gradientDescentAtoms(std::vector<Diagram> &DictDiagrams , const std::vector<std::vector<matchingTuple>> &matchings , const Diagram &Barycenter , const std::vector<Matrice> &gradLists , const int nb_points){
  std::vector<matrice> grad_list;
  //for(int i = 0 ; i< DictDiagrams.size() ; ++i){}
}
