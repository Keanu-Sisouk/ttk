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
#include <vtkVariantArray.h>
//#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
//#include <vtkTable.h>

#include <PersistenceDiagramDistanceMatrix.h>
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
    // info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
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

  // Number of input diagrams
  // int numInputs = 0;

  if(blocks != nullptr) {
    int numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
      // if(this->GetMTime() < input[i]->GetMTime()) {
      //  needUpdate_ = true;
      //}
    }
  }

  const size_t nDiags = inputDiagrams.size();

  std::vector<ttk::DiagramType> dictDiagrams(nDiags);
  for(size_t i = 0; i < nDiags; ++i) {
    auto &atom = dictDiagrams[i];
    // double max_dimension2 = getPersistenceDiagram(
    // atom, vtkUnstructuredGrid::SafeDownCast(inputDiagrams[i]));
    const auto ret = VTUToDiagram(atom, inputDiagrams[i], *this);
    // for(size_t k = 0; k < atom.size(); ++k) {
    //   DiagramTuple &t = atom[k];
    //   std::cout << "Pair atoms: " << std::get<6>(t) << ", " <<
    //   std::get<10>(t)
    //             << std::endl;
    // }
    if(ret != 0) {
      this->printWrn("Could not read Persistence Diagram");
      // return 0;
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

  std::cout << "PASSED !!!!" << std::endl;
  std::vector<vtkDataArray *> inputWeights;
  int numWeights = weights_vtk->GetNumberOfRows();
  for(int i = 0; i < weights_vtk->GetNumberOfColumns(); ++i) {
    std::cout << weights_vtk->GetColumnName(i) << "\n";
  }

  // this->printMsg(std::to_string(numWeights));
  if(weights_vtk != nullptr) {
    // int numWeights = weights_vtk->GetNumberOfColumns();
    // this->printMsg(std::to_string(numWeights));
    inputWeights.resize(nDiags);
    for(size_t i = 0; i < nDiags; ++i) {
      std::string name{"Atom"};
      zeroPad(name, nDiags, i);
      // std::cout << name << "\n";
      // const auto array = weights_vtk->GetColumnByName(name.c_str());
      // array->PrintSelf(std::cout, vtkIndent{});
      inputWeights[i] = vtkDataArray::SafeDownCast(
        weights_vtk->GetColumnByName(name.c_str()));
      // if(this->GetMTime() < input[i]->GetMTime()) {
      //  needUpdate_ = true;
      //}
    }
  }

  std::cout << "PASSED 2 !!!!!" << std::endl;

  // const int nWeights = ;
  std::vector<std::vector<double>> vectorWeights(numWeights);
  for(int i = 0; i < numWeights; ++i) {
    std::vector<double> &t1 = vectorWeights[i];
    // vtkDoubleArray &t2 = inputWeights[i];
    for(size_t j = 0; j < nDiags; ++j) {
      // double weight = t1[j];
      std::cout << "ICI???" << std::endl;
      double weight = inputWeights[j]->GetTuple1(i);
      t1.push_back(weight);
    }
  }
  std::cout << "PASSED 3!!!!!!" << std::endl;
  std::vector<ttk::DiagramType> Barycenters(numWeights);

  if(!ComputePoints) {
    this->execute(dictDiagrams, vectorWeights, Barycenters);
  }



  // this->printMsg("=====ICI?======");
  auto output_dgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_coordinates = vtkTable::GetData(outputVector, 1);
  output_dgm->SetNumberOfBlocks(numWeights);
  output_coordinates->SetNumberOfRows(numWeights);
  // this->printMsg(std::to_string(nWeights));
  // for(int i = 0; i < nWeights; ++i) {
  //   // vtkUnstructuredGrid temp =
  //   // vtkUnstructuredGrid::SafeDownCast(output_dgm->GetBlock(i));
  //   vtkNew<vtkUnstructuredGrid> vtu;
  //   ttk::Diagram &diagram = Barycenters[i];
  //   double max_persistence = getMaxPersistence(diagram);
  //   diagramToVTU(vtu, diagram, max_persistence);
  //   // this->printMsg("=====HERE?======");
  //   output_dgm->SetBlock(i, vtu);
  //   // this->printMsg("=====HERE2?=====");
  // }
  //
  // double max_persistence = getMaxPersistence(diagram);


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
  const std::vector<ttk::DiagramType> &atoms,
  vtkTable *weights_vtk,
  const std::vector<std::vector<double>> &weights,
  const double spacing,
  const double max_persistence) const {

  const auto nDiags = diags.size();
  const auto nAtoms = atoms.size();

  ttk::SimplexId n_existing_blocks = ShowAtoms ? nAtoms : 0;

  output->SetNumberOfBlocks(nDiags + n_existing_blocks);
  std::vector<std::pair<double, double>> coords(nAtoms);
  std::vector<std::pair<double, double>> true_coords(nAtoms);

  vtkNew<vtkDoubleArray> dummy{};



  if(nAtoms == 2) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(2);
    const auto distMatrix = MatrixCalculator.execute(atoms, nInputs);
    coords[0].first = 0.;
    true_coords[0].first = 0.;
    coords[0].second = 0.;
    true_coords[0].first = 0.;
    coords[1].first = spacing * distMatrix[0][1];
    true_coords[1].first = distMatrix[0][1];

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i].first;
        double Y = coords[i].second;
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else if(nAtoms == 3) {
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t, 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix
      = MatrixCalculator.execute(atoms, nInputs);
    coords[0].first = 0.;
    true_coords[0].first = 0.;
    coords[0].second = 0.;
    true_coords[0].second = 0.;
    coords[1].first = spacing * distMatrix[0][1];
    true_coords[1].first = distMatrix[0][1];
    coords[1].second = 0.;
    true_coords[0].second = 0.;
    double distOpposed = distMatrix[2][1];
    double firstDist = distMatrix[0][1];
    double distAdja = distMatrix[0][2];
    double alpha = std::acos(
      (distOpposed * distOpposed - firstDist * firstDist - distAdja * distAdja)
      / (-2. * firstDist * distAdja));
    coords[2].first = spacing * distAdja * std::cos(alpha);
    true_coords[2].first = distAdja * std::cos(alpha);
    coords[2].second = spacing * distAdja * std::sin(alpha);
    true_coords[2].second = distAdja * std::sin(alpha);

    if(ShowAtoms) {
      for(size_t i = 0; i < nAtoms; ++i) {
        double X = coords[i].first;
        double Y = coords[i].second;
        vtkNew<vtkUnstructuredGrid> vtu{};
        DiagramToVTU(vtu, atoms[i], dummy, *this, 3, false);
        TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
        output->SetBlock(i, vtu);
      }
    }

  } else {

    for(size_t i = 0; i < nAtoms; ++i) {
      const auto angle
        = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nAtoms);
      double X = spacing * max_persistence * std::cos(angle);
      double Y = spacing * max_persistence * std::sin(angle);
      coords[i].first = X;
      coords[i].second = Y;
      true_coords[i].first = max_persistence * std::cos(angle);
      true_coords[i].second = max_persistence * std::sin(angle);

      if(ShowAtoms) {
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


      std::cout << "SIZE OF BARY: " << i << " IS " << diags[i].size() << "\n";
      DiagramToVTU(vtu, diags[i], dummy, *this, 3, false);

    }

    // for(int j = 0 ; j < weights_vtk->GetNumberOfColumns() ; ++i){
    // auto array = weights_vtk->GetColumn(j);
    // std::string name{"ClusterID"};
    // if(strcmp(name.c_str(), weights_vtk->GetColumnName(j)) != 0){
    // continue;
    //}
    // vtkNew<vtkDataArray> clusterId{};
    // clusterId->SetName("ClusterID");
    // clusterId->SetNumberOfComponents(1);
    // clusterId->SetNumberOfTuples(1);
    // clusterId->SetNumberOfValues(1);
    // clusterId->SetTuple(0,i,array);
    // clusterId->Fill(array
    // array->GetTuples(i,i,clusterId);
    // vtu->GetFieldData()->AddArray(clusterId);
    //}

    double X = 0;
    double Y = 0;
    for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
      X += weights[i][iAtom] * coords[iAtom].first;
      Y += weights[i][iAtom] * coords[iAtom].second;
    }

    TranslateDiagram(vtu, std::array<double, 3>{X, Y, 0.0});
    output->SetBlock(i + n_existing_blocks, vtu);
  }

  for(size_t i = 0; i < 2; ++i) {
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(nDiags);
    std::string name;

    if(i == 0) {
      name = "X";
    } else {
      name = "Y";
    }
    // name.append(std::to_string(i));
    col->SetName(name.c_str());
    for(size_t j = 0; j < nDiags; ++j) {
      double temp = 0;
      for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        if(i == 0) {
          temp += weights[j][iAtom] * true_coords[iAtom].first;
        } else {
          temp += weights[j][iAtom] * true_coords[iAtom].second;
        }
      }
      col->SetValue(j, temp);
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
        row->SetValue(j, true_coords[i].first);
      } else if(strcmp(output_coordinates->GetColumnName(j), "Y") == 0) {
        row->SetValue(j, true_coords[i].second);
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
