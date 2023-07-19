/// \ingroup base
/// \class ttk::PersistenceDiagramDistanceMatrix
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramClustering
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

#include <array>

#include <Debug.h>
#include <PersistenceDiagramUtils.h>

namespace ttk {

  class PersistenceDiagramSlicedWasserstein : virtual public Debug {

  public:
    PersistenceDiagramSlicedWasserstein() {
      this->setDebugMsgPrefix("PersistenceDiagramDistanceMatrix");
    }

    double execute(const ttk::DiagramType &diag1,
                   const ttk::DiagramType &diag2,
                   int sampleNumber);


  protected:





    void projectionOnThetaLine(const ttk::DiagramType &diag,
                               const std::vector<std::array<double, 2>> &proj,
                               std::vector<std::array<double, 2>> &projOnTheta,
                               double theta);
    
    void augmentDiagram(const ttk::DiagramType &diag,
                        std::vector<std::array<double, 2>> &proj);

  };
} // namespace ttk
