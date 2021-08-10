#pragma once

#include <algorithm>
#include <array>
#include <PersistenceDiagramBarycenter.h>
#include <PDClustering.h>
#include <Wrapper.h>


namespace ttk {
  using diagramTuple = std::tuple<
    /** Vertex Id of low pair element */
    ttk::SimplexId,
    /** Critical Type of low pair element */
    ttk::CriticalType,
    /** Vertex Id of high pair element */
    ttk::SimplexId,
    /** Critical Type of high pair element */
    ttk::CriticalType,
    /** Pair persistence value */
    double,
    /** Pair type */
    ttk::SimplexId,
    /** Pair birth */
    double,
    /** Low pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float,
    /** Pair death */
    double,
    /** High pair element 3D coordinates */
    // TODO use std::array<float, 3>
    float,
    float,
    float>;

  using Diagram = std::vector<diagramTuple>;
  using Matrice = std::vector<std::vector<double>>;


  class ConstrainedGradientDescent : public Debug{

  public:
    ConstrainedGradientDescent(){
      this->setDebugMsgPrefix("ConstrainedGradientDescent");
    };

    void executeWeightsProjected(std::vector<double> &weights , const std::vector<double> &grad , const int epoch,const int nb_points);
    //void executeAtoms(std::vector<Diagram> &DictDiagrams);

    inline void setNbAtoms(const int nbAtoms){
      NbAtoms = nbAtoms;
    }

  protected:
    void projectionOnSimplex(std::vector<double> &weights);

    void gradientDescentWeights(std::vector<double> &weights , const std::vector<double> &grad ,const int epoch, const int nb_points);

    void gradientDescentAtoms(std::vector<Diagram> &DictDiagrams ,const std::vector<std::vector<matchingTuple>> &matchings , const Diagram &Barycenter , const std::vector<Matrice> &gradLists , const int nb_points);

  };
}
