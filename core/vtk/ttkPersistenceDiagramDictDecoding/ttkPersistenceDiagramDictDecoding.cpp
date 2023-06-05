#include <ttkPersistenceDiagramDictDecoding.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkVariantArray.h>

#include <ttkPersistenceDiagramUtils.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagramDictDecoding);

ttkPersistenceDiagramDictDecoding::ttkPersistenceDiagramDictDecoding() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

int ttkPersistenceDiagramDictDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

int ttkPersistenceDiagramDictDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

int ttkPersistenceDiagramDictDecoding::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto weights_vtk = vtkTable::GetData(inputVector[1]);

  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  if(blocks != nullptr) {
    int numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
    }
  }

  const size_t nDiags = inputDiagrams.size();

  std::vector<ttk::DiagramType> dictDiagrams(nDiags);
  for(size_t i = 0; i < nDiags; ++i) {
    auto &atom = dictDiagrams[i];
    const auto ret = VTUToDiagram(atom, inputDiagrams[i], *this);
    if(ret != 0) {
      this->printWrn("Could not read Persistence Diagram");
    }
  }

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  std::vector<vtkDataArray *> inputWeights;
  int numWeights = weights_vtk->GetNumberOfRows();
  for(int i = 0; i < weights_vtk->GetNumberOfColumns(); ++i) {
    std::cout << weights_vtk->GetColumnName(i) << "\n";
  }

  if(weights_vtk != nullptr) {
    inputWeights.resize(nDiags);
    for(size_t i = 0; i < nDiags; ++i) {
      std::string name{"Atom"};
      zeroPad(name, nDiags, i);
      inputWeights[i] = vtkDataArray::SafeDownCast(
        weights_vtk->GetColumnByName(name.c_str()));
    }
  }

  std::vector<std::vector<double>> vectorWeights(numWeights);
  for(int i = 0; i < numWeights; ++i) {
    std::vector<double> &t1 = vectorWeights[i];
    for(size_t j = 0; j < nDiags; ++j) {
      double weight = inputWeights[j]->GetTuple1(i);
      t1.push_back(weight);
    }
  }
  std::vector<ttk::DiagramType> Barycenters(numWeights);

  if(!ComputePoints) {
    this->execute(dictDiagrams, vectorWeights, Barycenters);
  }

  auto output_dgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_coordinates = vtkTable::GetData(outputVector, 1);
  output_dgm->SetNumberOfBlocks(numWeights);
  output_coordinates->SetNumberOfRows(numWeights);

  outputDiagrams(output_dgm, output_coordinates, Barycenters, dictDiagrams,
                 weights_vtk, vectorWeights, Spacing, 1);

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation

  // make a SHALLOW copy of the input
  // outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  // outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}

void ttkPersistenceDiagramDictDecoding::outputDiagrams(
  vtkMultiBlockDataSet *output,
  vtkTable *output_coordinates,
  const std::vector<ttk::DiagramType> &diags,
  std::vector<ttk::DiagramType> &atoms,
  vtkTable *weights_vtk,
  const std::vector<std::vector<double>> &weights,
  const double spacing,
  const double max_persistence) const {

  const auto nDiags = diags.size();
  const auto nAtoms = atoms.size();

  ttk::SimplexId n_existing_blocks = ShowAtoms ? nAtoms : 0;

  output->SetNumberOfBlocks(nDiags + n_existing_blocks);
  std::vector<std::array<double, 3>> coords(nAtoms);
  std::vector<std::array<double, 3>> true_coords(nAtoms);
  std::vector<double> xVector(nDiags);
  std::vector<double> yVector(nDiags);
  std::vector<double> zVector(nDiags, 0.);
  vtkNew<vtkDoubleArray> dummy{};

  computeAtomsCoordinates(atoms, weights, coords, true_coords, xVector, yVector,
                          zVector, spacing, max_persistence, nAtoms);

  if(nAtoms == 2) {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i][0];
        double Y = coords[i][1];
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else if(nAtoms == 3) {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i][0];
        double Y = coords[i][1];
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else {

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        const auto angle
          = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nAtoms);
        double X = spacing * max_persistence * std::cos(angle);
        double Y = spacing * max_persistence * std::sin(angle);
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      
      }
    }
  }

  int numDiags = diags.size();

  for(int i = 0; i < numDiags; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    if(!ComputePoints) {
      DiagramToVTU(vtu, diags[i], dummy, *this, 3, false);
    }

    double X = 0;
    double Y = 0;
    for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
      X += weights[i][iAtom] * coords[iAtom][0];
      Y += weights[i][iAtom] * coords[iAtom][1];
    }

    TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
    output->SetBlock(i + n_existing_blocks, vtu);
  }

  for(size_t i = 0; i < 3; ++i) {
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(nDiags);
    std::string name;

    if(i == 0) {
      name = "X";
    } else if(i == 1) {
      name = "Y";
    } else {
      name = "Z";
    }
    col->SetName(name.c_str());
    for(size_t j = 0; j < nDiags; ++j) {
      if(i == 0) {
        col->SetValue(j, xVector[j]);
      } else if(i == 1) {
        col->SetValue(j, yVector[j]);
      } else {
        col->SetValue(j, zVector[j]);
      }
    }
    col->Modified();
    output_coordinates->AddColumn(col);
  }

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };

  for(int i = 0; i < weights_vtk->GetNumberOfColumns(); ++i) {
    int test = 0;
    const auto array = weights_vtk->GetColumn(i);

    for(size_t j = 0; j < nDiags; ++j) {
      std::string name{"Atom"};
      zeroPad(name, nDiags, j);
      if(strcmp(name.c_str(), weights_vtk->GetColumnName(i)) == 0) {
        test += 1;
      }
    }
    if(test > 0) {
      continue;
    }
    output_coordinates->AddColumn(array);
  }

  for(size_t i = 0; i < nAtoms; ++i) {
    vtkNew<vtkVariantArray> row{};
    row->SetNumberOfValues(output_coordinates->GetNumberOfColumns());
    std::cout << "number of values: " << row->GetNumberOfValues() << std::endl;
    for(int j = 0; j < output_coordinates->GetNumberOfColumns(); ++j) {
      if(strcmp(output_coordinates->GetColumnName(j), "X") == 0) {
        row->SetValue(j, true_coords[i][0]);
      } else if(strcmp(output_coordinates->GetColumnName(j), "Y") == 0) {
        row->SetValue(j, true_coords[i][1]);
      } else if(strcmp(output_coordinates->GetColumnName(j), "Z") == 0) {
        row->SetValue(j, true_coords[i][2]);
      } else if(strcmp(output_coordinates->GetColumnName(j), "ClusterID")
                == 0) {
        row->SetValue(j, -1);
      } else {
        continue;
      }
    }
    row->Modified();
    output_coordinates->InsertNextRow(row);
  }
}

double ttkPersistenceDiagramDictDecoding::getMaxPersistence(
  const ttk::DiagramType &diagram) const {

  double max_persistence{0};
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &t = diagram[i];
    const double &pers = t.persistence();
    max_persistence = std::max(pers, max_persistence);
  }
  return max_persistence;
}
