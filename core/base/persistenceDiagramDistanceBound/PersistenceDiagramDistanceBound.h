/// \ingroup base
/// \class ttk::PersistenceDiagramDistanceBound
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK %persistenceDiagramDistanceBound processing package.
///
/// %PersistenceDiagramDistanceBound is a TTK processing package that takes a scalar field on the
/// input and produces a scalar field on the output.
///
/// \sa ttkPersistenceDiagramDistanceBound.cpp for a usage example.

#pragma once

// Standard.
#include "AbstractTriangulation.h"
#include "DataTypes.h"
#include <cmath>
#include <cstdlib>
#include <string>

// Base code.
#include <Debug.h>
#include <Geometry.h>
#include <Triangulation.h>

#ifdef TTK_ENABLE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif // TTK_ENABLE_EIGEN



namespace ttk {

  class PersistenceDiagramDistanceBound : virtual public Debug {
  public:
    PersistenceDiagramDistanceBound();

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
        triangulation->preconditionBoundaryVertices();
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    template <class dataType>
    double execute(const dataType *const inputData1,
                const dataType *const inputData2,
                const AbstractTriangulation *triangulation,
                const SimplexId vertexNumber);

    template <class dataType>
    int computeLinf(const dataType *const input1,
                    const dataType *const input2,
                    const SimplexId vertexNumber);

    template <class dataType>
    double computeLipsCst(const AbstractTriangulation *triangulation,
                          const dataType *const input);


    template <class dataType>
    int computeLp(const dataType *const input1,
                  const dataType *const input2,
                  const AbstractTriangulation *triangulation);

    inline double getResult() {
      return result;
    }

    inline void setPrintRes(const bool data) {
      this->printRes = data;
    }

    template <typename type>
    static type abs_diff(const type var1, const type var2) {
      return (var1 > var2) ? var1 - var2 : var2 - var1;
    }

  protected:
    double result{};
    bool printRes{true};
    double kDim = 2.;
    int wassersteinParam = 2;
  };
} // namespace ttk

// template functions
template <class dataType>
double ttk::PersistenceDiagramDistanceBound::execute(const dataType *const inputData1,
                            const dataType *const inputData2,
                            const AbstractTriangulation *triangulation,
                            const SimplexId vertexNumber) {

  Timer t;
  int status;
  double bound;
// Check variables consistency
#ifndef TTK_ENABLE_KAMIKAZE
  if(inputData1 == nullptr || inputData2 == nullptr) {
    return -1;
  }
#endif

  status = computeLp(inputData1, inputData2, triangulation);
  bound = pow(result, 1./wassersteinParam);
  // status = computeLinf(inputData1, inputData2, vertexNumber);

  // double lips1 = computeLipsCst(triangulation, inputData1);
  // double lips2 = computeLipsCst(triangulation, inputData2);

  // unsigned int dim = triangulation->getNumberOfVertices();
  // float p0[3], p1[3];
  // triangulation->getVertexPoint(0, p0[0], p0[1], p0[2]);
  // triangulation->getVertexPoint(static_cast<int>(dim) - 1, p1[0], p1[1], p1[2]);

  // float wide = p1[0] - p0[0];
  // float height = p1[1] - p0[1];

  // double lebesgueMeasure = wide*height;


  // if(lips1 > lips2){
  //   bound = pow((2. / M_PI) * lebesgueMeasure * pow(lips1, kDim), 1./static_cast<double>(wassersteinParam))
  //       *pow(result, 1. - kDim/static_cast<double>(wassersteinParam));
  // } else {
  //   bound = pow((2. / M_PI) * lebesgueMeasure * pow(lips2, kDim), 1./static_cast<double>(wassersteinParam))
  //       *pow(result, 1. - kDim/static_cast<double>(wassersteinParam));
  // }
  if(this->printRes) {

    this->printMsg("Stability bound: " + std::to_string(bound));
    this->printMsg(
      "Data-set processed", 1.0, t.getElapsedTime(), this->threadNumber_);

  }

  return bound;
}

template <class dataType>
int ttk::PersistenceDiagramDistanceBound::computeLinf(const dataType *const input1,
                                const dataType *const input2,
                                const ttk::SimplexId vertexNumber) {
  if(vertexNumber < 1)
    return 0;

  dataType maxValue = abs_diff<dataType>(input1[0], input2[0]);

// Compute difference for each point.
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) reduction(max : maxValue)
#endif
  for(ttk::SimplexId i = 1; i < vertexNumber; ++i) {
    const dataType iter = abs_diff<dataType>(input1[i], input2[i]);
    if(iter > maxValue)
      maxValue = iter;
  }

  // Affect result.
  result = (double)maxValue;


  return 0;
}


template <class dataType>
double ttk::PersistenceDiagramDistanceBound::computeLipsCst(
  const AbstractTriangulation *triangulation,
  const dataType *const input) {

  unsigned int dim = triangulation->getNumberOfVertices();

  double lips = 0.;
  std::vector<double> allNorm(dim);

#ifdef TTK_ENABLE_EIGEN

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(unsigned int i = 0; i < dim ; ++i){

    double iter;
    auto neighborNum = triangulation->getVertexNeighborNumber(i);

    float p0[3];
    triangulation->getVertexPoint(i, p0[0], p0[1], p0[2]);
    Eigen::MatrixXf V(neighborNum, 2);
    Eigen::VectorXf f(neighborNum);
    float value = static_cast<float>(input[i]);

    for(int j = 0; j < neighborNum ; ++j){

      SimplexId neighborId;
      triangulation->getVertexNeighbor(i, j, neighborId);
      float p1[3];
      triangulation->getVertexPoint(neighborId, p1[0], p1[1], p1[2]);
      float value2 =  static_cast<float>(input[neighborId]);

      V(j,0) = (float)p1[0] - (float)p0[0];
      V(j,1) = (float)p1[1] - (float)p0[1];
      f(j) = value2 - value;
    }

    Eigen::VectorXf grad = V.colPivHouseholderQr().solve(f);
    iter = grad.norm();

    allNorm[i] = iter;
  }

  lips = *std::max_element(allNorm.begin(), allNorm.end());
#endif //TTK_ENABLE_EIGEN
  return lips;
}

template <class dataType>
int ttk::PersistenceDiagramDistanceBound::computeLp(const dataType *const input1,
                                const dataType *const input2,
                                const AbstractTriangulation *triangulation) {

unsigned int dim = triangulation->getNumberOfCells();

dataType Lp_norm = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
  for(unsigned int i = 1; i < dim; ++i) {
    dataType temp1 = 0;
    dataType temp2 = 0;
    std::vector<SimplexId> verticesId{3};

    for(int j = 0; j < 3 ; ++j){
      triangulation->getTriangleVertex(i, j , verticesId[j]);
    }

    dataType iter1;
    dataType iter2;
    for(int j = 0; j < 3; ++j){
      iter1 = input1[verticesId[j]];
      iter2 = input2[verticesId[j]];
      if(temp1 < iter1){
        temp1 = iter1;
      }
      if(temp2 < iter2){
        temp2 = iter2;
      }
    }
    Lp_norm += pow(static_cast<double>(abs_diff<dataType>(temp1, temp2)), wassersteinParam);
  }

  Lp_norm = pow(Lp_norm , 1./wassersteinParam);
  result = (double)Lp_norm;
  return 0;
}