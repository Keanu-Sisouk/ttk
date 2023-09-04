/// \ingroup vtk
/// \class ttkPersistenceDiagramDistanceBoundMatrix
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief Computes a distance matrix using LDistance between several
/// input datasets with the same number of points
///
/// \sa PersistenceDiagramDistanceBoundMatrix
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

// VTK Module
#include <ttkPersistenceDiagramDistanceBoundMatrixModule.h>

// TTK code includes
#include <PersistenceDiagramDistanceBoundMatrix.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCEDIAGRAMDISTANCEBOUNDMATRIX_EXPORT ttkPersistenceDiagramDistanceBoundMatrix
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDistanceBoundMatrix {

public:
  static ttkPersistenceDiagramDistanceBoundMatrix *New();

  vtkTypeMacro(ttkPersistenceDiagramDistanceBoundMatrix, ttkAlgorithm);


protected:
  ttkPersistenceDiagramDistanceBoundMatrix();
  ~ttkPersistenceDiagramDistanceBoundMatrix() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  template <typename T>
  int dispatch(std::vector<std::vector<double>> &distanceMatrix,
               const ttk::AbstractTriangulation *triangulation,
               const std::vector<vtkDataSet *> &inputData,
               const size_t nPoints);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
