#pragma once

#include <Wrapper.h>
#include <algorithm>
#include <array>


namespace ttk {
  using DiagramTuple = std::tuple<
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

  using Diagram = std::vector<DiagramTuple>;
  using Matrix = std::vector<std::vector<double>>;
  using MatchingTuple = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;
  class InitRandomDict : public Debug {

  public:
    InitRandomDict() {
      this->setDebugMsgPrefix("InitRandomDict");
    };

    void execute(std::vector<Diagram> &DictDiagrams,
                 const std::vector<Diagram> &datas,
                 const int nbAtoms,
                 const int seed);

  protected:
  };
} // namespace ttk
