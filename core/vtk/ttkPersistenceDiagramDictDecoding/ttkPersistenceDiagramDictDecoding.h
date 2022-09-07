#pragma once

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>
// VTK Module
#include <ttkPersistenceDiagramDictDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <PersistenceDiagramDictDecoding.h>

class TTKPERSISTENCEDIAGRAMDICTDECODING_EXPORT ttkPersistenceDiagramDictDecoding
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDictDecoding {

private:
public:
  static ttkPersistenceDiagramDictDecoding *New();
  vtkTypeMacro(ttkPersistenceDiagramDictDecoding, ttkAlgorithm);

  vtkGetMacro(Spacing, double);
  vtkSetMacro(Spacing, double);

  vtkGetMacro(ShowAtoms, int);
  vtkSetMacro(ShowAtoms, int);

  vtkGetMacro(ProgBarycenter, int);
  vtkSetMacro(ProgBarycenter, int);

protected:
  ttkPersistenceDiagramDictDecoding();
  ~ttkPersistenceDiagramDictDecoding() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  void outputDiagrams(vtkMultiBlockDataSet *output,
                      vtkTable *output_coordinates,
                      const std::vector<ttk::DiagramType> &diags,
                      const std::vector<ttk::DiagramType> &atoms,
                      vtkTable *weights_vtk,
                      const std::vector<std::vector<double>> &weights,
                      const double spacing,
                      const double max_persistence) const;

  double getMaxPersistence(const ttk::DiagramType &diagram) const;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  double Spacing{};
  int ShowAtoms{1};
  bool ComputePoints{false};
};
