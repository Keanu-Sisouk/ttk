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
                   int &sampleNumber,
                   int maxSampleNb);


    inline void setNbPoints(int data){
      nbPoints = data;
    }

    inline void setUseQuasiMC(bool data){
      useQuasiMC = data;
    }

    inline void setNbOfProjused(int data){
      nbOfProjUsed += data;
    }

  protected:


    void projectionOnThetaLine(const ttk::DiagramType &diag,
                               const std::vector<std::array<double, 2>> &proj,
                               std::vector<std::array<double, 2>> &projOnTheta,
                               std::vector<int> &originIndices,
                               std::vector<double> &scalarProd,
                               double theta,
                               bool vertical);
    
    void augmentDiagram(const ttk::DiagramType &diag,
                        std::vector<std::array<double, 2>> &proj);

    void slicedTransport(std::vector<ttk::MatchingType> &matchings,
                        const ttk::DiagramType &diag1,
                        const ttk::DiagramType &diag2,
                        int sampleNumber);

    void projectionOnThetaLine(const std::vector<std::array<double,2>> &limitMeasure,
                              std::vector<std::array<double, 2>> &projOnTheta,
                              double theta);

    void getMatchings(std::vector<MatchingType> &matchings,
                      const ttk::DiagramType &diag1,
                      const ttk::DiagramType &diag2,
                      const std::vector<std::array<double, 2>> &proj,
                      const std::vector<std::array<double, 2>> &limitMeasure);

    double computeNormGradient(const ttk::DiagramType &diag1,
                               const ttk::DiagramType &diag2,
                               const std::vector<std::array<double, 2>> &proj1,
                               const std::vector<std::array<double, 2>> &proj2,
                               const std::vector<int> &originIndices1,
                               const std::vector<int> &originIndices2,
                               const std::vector<double> &scalarProd1,
                               const std::vector<double> &scalarProd2,
                               double theta);

    double classicMonteCarlo(const ttk::DiagramType &diag1,
                             const ttk::DiagramType &diag2,
                             const std::vector<std::array<double, 2>> &proj1,
                             const std::vector<std::array<double, 2>> &proj2,
                             int &sampleNumber,
                             int maxSampleNb);

    double quasiMonteCarlo(const ttk::DiagramType &diag1,
                           const ttk::DiagramType &diag2,
                           const std::vector<std::array<double, 2>> &proj1,
                           const std::vector<std::array<double, 2>> &proj2,
                           int &sampleNumber,
                           int maxSampleNb);

    int EPOCH_MAX = 300;

    int nbPoints = 1;

    bool useQuasiMC = true;

    int nbOfProjUsed = 0;
  };
} // namespace ttk
