/// \ingroup base
/// \class ttkPersistenceDiagramDictEncoding
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa PersistenceDiagramDictEncoding

#pragma once

// VTK includes
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSetGet.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceDiagramDictEncodingModule.h>

// ttk code includes

#include <PersistenceDiagramDictEncoding.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class TTKPERSISTENCEDIAGRAMDICTENCODING_EXPORT ttkPersistenceDiagramDictEncoding
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagramDictEncoding {

private:
  int atomNumber_{3};
  int seed_{0};
  double percent_{0};

public:
  // enum class BACKEND{BORDER_INIT = 0 , RANDOM_INIT = 1 , FIRST_DIAGS = 2};

  static ttkPersistenceDiagramDictEncoding *New();

  vtkTypeMacro(ttkPersistenceDiagramDictEncoding, ttkAlgorithm);

  void SetWassersteinMetric(const std::string &data) {
    Wasserstein = (data == "inf") ? -1 : stoi(data);
    Modified();
  }
  std::string GetWassersteinMetric() {
    return Wasserstein == -1 ? "inf" : std::to_string(Wasserstein);
  }

  void SetAntiAlpha(double data) {
    data = 1 - data;
    if(data > 0 && data <= 1) {
      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
  }

  vtkGetMacro(Alpha, double);

  vtkSetMacro(percent_, double);
  vtkGetMacro(percent_, double);

  vtkSetMacro(OptimizeWeights, int);
  vtkGetMacro(OptimizeWeights, int);

  vtkSetMacro(OptimizeAtoms, int);
  vtkGetMacro(OptimizeAtoms, int);

  vtkSetMacro(MaxEigenValue, int);
  vtkGetMacro(MaxEigenValue, int);

  vtkSetMacro(Fusion, int);
  vtkGetMacro(Fusion, int);

  vtkSetMacro(ProgBarycenter, int);
  vtkGetMacro(ProgBarycenter, int);

  vtkSetMacro(explicitSolWeights, int);
  vtkGetMacro(explicitSolWeights, int);

  vtkSetMacro(explicitSol, int);
  vtkGetMacro(explicitSol, int);

  vtkSetMacro(MaxEpoch, int);
  vtkGetMacro(MaxEpoch, int);

  vtkSetMacro(ProgApproach, int);
  vtkGetMacro(ProgApproach, int);

  vtkSetMacro(StopCondition, int);
  vtkGetMacro(StopCondition, int);

  vtkSetMacro(CompressionMode, int);
  vtkGetMacro(CompressionMode, int);

  vtkSetMacro(DimReductMode, int);
  vtkGetMacro(DimReductMode, int);

  vtkSetMacro(sortedForTest, int);
  vtkGetMacro(sortedForTest, int);

  vtkSetMacro(CreationFeatures, int);
  vtkGetMacro(CreationFeatures, int);

  vtkSetMacro(atomNumber_, int);
  vtkGetMacro(atomNumber_, int);

  vtkSetMacro(seed_, int);
  vtkGetMacro(seed_, int);

  vtkSetMacro(DeltaLim, double);
  vtkGetMacro(DeltaLim, double);

  vtkSetMacro(Lambda, double);
  vtkGetMacro(Lambda, double);

  ttkSetEnumMacro(BackEnd, BACKEND);
  vtkGetEnumMacro(BackEnd, BACKEND);

  void SetPairType(const int data) {
    switch(data) {
      case(0):
        this->setDos(true, false, false);
        break;
      case(1):
        this->setDos(false, true, false);
        break;
      case(2):
        this->setDos(false, false, true);
        break;
      default:
        this->setDos(true, true, true);
        break;
    }
    Modified();
  }
  int GetPairType() {
    if(do_min_ && do_sad_ && do_max_) {
      return -1;
    } else if(do_min_) {
      return 0;
    } else if(do_sad_) {
      return 1;
    } else if(do_max_) {
      return 2;
    }
    return -1;
  }

  void SetConstraint(const int arg_) {
    this->setConstraint(arg_);
    this->Modified();
  }
  int GetConstraint() {
    switch(this->Constraint) {
      case ConstraintType::FULL_DIAGRAMS:
        return 0;
      case ConstraintType::NUMBER_PAIRS:
        return 1;
      case ConstraintType::ABSOLUTE_PERSISTENCE:
        return 2;
      case ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG:
        return 3;
      case ConstraintType::RELATIVE_PERSISTENCE_GLOBAL:
        return 4;
    }
    return -1;
  }

  vtkSetMacro(MaxNumberOfPairs, unsigned int);
  vtkGetMacro(MaxNumberOfPairs, unsigned int);

  vtkSetMacro(MinPersistence, double);
  vtkGetMacro(MinPersistence, double);

  vtkSetMacro(CompressionFactor, double);
  vtkGetMacro(CompressionFactor, double);

protected:
  ttkPersistenceDiagramDictEncoding();
  ~ttkPersistenceDiagramDictEncoding() override = default;

  // BACKEND BackEnd{BACKEND::BORDER_INIT};
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
