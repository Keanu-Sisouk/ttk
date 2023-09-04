/// \ingroup base
/// \class ttk::PersistenceDiagramDistanceBoundMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2020
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

#include <PersistenceDiagramDistanceBound.h>
#include <Wrapper.h>

#include <string>
#include <vector>

namespace ttk {
  class PersistenceDiagramDistanceBoundMatrix : virtual public Debug {
  public:
    PersistenceDiagramDistanceBoundMatrix();

    int preconditionTriangulation(AbstractTriangulation *triangulation) {
      // Pre-condition functions.
      if(triangulation) {
        triangulation->preconditionBoundaryVertices();
        triangulation->preconditionVertexNeighbors();
      }
      return 0;
    }

    template <typename TIn>
    int execute(std::vector<std::vector<double>> &output,
                const ttk::AbstractTriangulation *triangulation,
                const std::vector<const TIn *> &inputs,
                const size_t nPoints) const;

  protected:
  };
} // namespace ttk

template <typename TIn>
int ttk::PersistenceDiagramDistanceBoundMatrix::execute(std::vector<std::vector<double>> &output,
                                  const ttk::AbstractTriangulation *triangulation,
                                  const std::vector<const TIn *> &inputs,
                                  const size_t nPoints) const {

  const size_t nInputs = inputs.size();

  output.resize(nInputs);
  for(size_t i = 0; i < nInputs; ++i){
    output[i].resize(nInputs);
  }

  PersistenceDiagramDistanceBound worker{};
  worker.setThreadNumber(1);
  worker.setPrintRes(false);

  // compute matrix upper triangle
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) firstprivate(worker)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nInputs; ++i) {
    output[i][i] = 0;
    for(size_t j = i + 1; j < nInputs; ++j) {
      // call execute with nullptr output
      double bound = worker.execute(inputs[i], inputs[j], triangulation,nPoints);
      // store result
      output[i][j] = bound;
      output[j][i] = bound;
    }
  }

  return 0;
}

