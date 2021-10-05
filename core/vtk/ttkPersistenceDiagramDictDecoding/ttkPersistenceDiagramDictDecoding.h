#pragma once

#include <vtkInformation.h>
#include <vtkInformationVector.h>
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

protected:
  ttkPersistenceDiagramDictDecoding();
  ~ttkPersistenceDiagramDictDecoding() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  double getPersistenceDiagram(ttk::Diagram &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  void diagramToVTU(vtkUnstructuredGrid *output,
                    const Diagram &diagram,
                    const double max_persistence) const;

  double getMaxPersistence(Diagram &diagram);

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
