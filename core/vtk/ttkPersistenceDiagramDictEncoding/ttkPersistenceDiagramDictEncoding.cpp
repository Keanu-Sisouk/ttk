#include <ttkMacros.h>
#include <ttkPersistenceDiagramDictEncoding.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkTable.h>

vtkStandardNewMacro(ttkPersistenceDiagramDictEncoding);

ttkPersistenceDiagramDictEncoding::ttkPersistenceDiagramDictEncoding() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(6);
}

int ttkPersistenceDiagramDictEncoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    // info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagramDictEncoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    /*info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);*/
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else if(port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else if(port == 3) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else if(port == 4) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else if(port == 5) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceDiagramDictEncoding::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  ttk::Memory m;

  // Get input data
  // std::vector<vtkUnstructuredGrid *> inputDiagrams;

  // auto nBlocks = inputVector[0]->GetNumberOfInformationObjects();
  // std::vector<vtkMultiBlockDataSet *> blocks(nBlocks);

  // if(nBlocks > 2) {
  //  this->printWrn("Only dealing with the first two MultiBlockDataSets");
  //  nBlocks = 2;
  //}

  // number of diagrams per input block
  std::array<size_t, 2> nInputs{0, 0};

  // for(int i = 0; i < nBlocks; ++i) {
  //  blocks[i] = vtkMultiBlockDataSet::GetData(inputVector[0], i);
  //  if(blocks[i] != nullptr) {
  //    nInputs[i] = blocks[i]->GetNumberOfBlocks();
  //    for(size_t j = 0; j < nInputs[i]; ++j) {
  //      inputDiagrams.emplace_back(
  //        vtkUnstructuredGrid::SafeDownCast(blocks[i]->GetBlock(j)));
  //    }
  //  }
  //}

  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto atomBlocks = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // Flat storage for diagrams extracted from blocks
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  std::vector<vtkUnstructuredGrid *> inputAtoms;

  // Number of input diagrams
  // int numInputs = 0;
  // printf("Atom number %d", numAtom);

  int numInputAtoms = 0;
  if(atomBlocks != nullptr) {
    numInputAtoms = atomBlocks->GetNumberOfBlocks();
    inputAtoms.resize(numInputAtoms);
    for(int i = 0; i < numInputAtoms; ++i) {
      inputAtoms[i]
        = vtkUnstructuredGrid::SafeDownCast(atomBlocks->GetBlock(i));
    }
  }
  if(BackEnd == BACKEND::INPUT_ATOMS) {
    atomNumber_ = numInputAtoms;
  }
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

  const int numAtom = this->GetatomNumber_();
  printMsg("Number of atoms: " + ttk::debug::output::YELLOW
           + ttk::debug::output::UNDERLINED + std::to_string(numAtom)
           + ttk::debug::output::ENDCOLOR + ttk::debug::output::ENDCOLOR);

  // total number of diagrams
  const int nDiags = inputDiagrams.size();

  // Sanity check
  for(const auto vtu : inputDiagrams) {
    if(vtu == nullptr) {
      this->printErr("Input diagrams are not all vtkUnstructuredGrid");
      return 0;
    }
  }

  // Set output
  // auto diagramsDistTable = vtkTable::GetData(outputVector);

  auto output_dgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_weights = vtkTable::GetData(outputVector, 1);
  auto output_loss = vtkTable::GetData(outputVector, 2);
  auto output_allLosses = vtkTable::GetData(outputVector, 3);
  auto true_output_loss = vtkTable::GetData(outputVector, 4);
  auto output_timers = vtkTable::GetData(outputVector, 5);
  // int numAtom = this->GetAtomNumber();
  output_dgm->SetNumberOfBlocks(numAtom);

  if(BackEnd == BACKEND::INPUT_ATOMS) {
    std::cout << "KONICHIWA" << std::endl;
    for(int i = 0; i < numAtom; ++i) {
      vtkNew<vtkUnstructuredGrid> vtu;
      vtu->DeepCopy(inputAtoms[i]);
      output_dgm->SetBlock(i, vtu);
    }
  } else {
    std::cout << "KOMBAWA" << std::endl;
    for(int i = 0; i < numAtom; ++i) {
      vtkNew<vtkUnstructuredGrid> vtu;
      vtu->DeepCopy(inputDiagrams[i]);
      output_dgm->SetBlock(i, vtu);
    }
  }
  std::cout << "BOUYASHAKA" << std::endl;

  std::vector<ttk::DiagramType> intermediateDiagrams(nDiags);
  std::vector<ttk::DiagramType> intermediateAtoms(numInputAtoms);

  double max_dimension_total = 0.0;
  for(int i = 0; i < nDiags; ++i) {

    const auto ret
      = VTUToDiagram(intermediateDiagrams[i], inputDiagrams[i], *this);

    double max_pers = this->getMaxPers(intermediateDiagrams[i]);
    double percentage = this->percent_;
    std::cout << "MAX PERS BEFORE FILTERING " << max_pers << std::endl;
    // if (max_dimension < this->getMaxPers(intermediateDiagrams[i])){
    // std::sort(intermediateDiagrams[i].begin(), intermediateDiagrams[i].end(),
    //           [](ttk::DiagramTuple &t1, ttk::DiagramTuple &t2) {
    //             return (std::get<10>(t1) - std::get<6>(t1))
    //                   > (std::get<10>(t2) - std::get<6>(t2));
    //           });
    // max_dimension =
    //}
    // ttk::DiagramTuple &temp = intermediateDiagrams[i][0];
    // double max_pers = std::get<10>(temp) - std::get<6>(temp);
    intermediateDiagrams[i].erase(
      std::remove_if(intermediateDiagrams[i].begin(),
                     intermediateDiagrams[i].end(),
                     [max_pers, percentage](ttk::PersistencePair &t) {
                       return (t.death.sfValue - t.birth.sfValue)
                              < (percentage / 100.) * max_pers;
                     }),
      intermediateDiagrams[i].end());

    auto &t = intermediateDiagrams[i][0];
    std::cout << "MAX PERS BEFORE ORDERING "
              << t.death.sfValue - t.birth.sfValue << std::endl;
    if(ret != 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
    if(max_dimension_total < max_pers) {
      max_dimension_total = max_pers;
    }
  }
  for(int i = 0; i < numInputAtoms; ++i) {
    const auto ret = VTUToDiagram(intermediateAtoms[i], inputAtoms[i], *this);
    if(ret != 0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
  }

  // double max_dimension_total = 0.0;
  // for(int i = 0; i < 1; ++i) {
  //   double max_dimension = getPersistenceDiagram(
  //     intermediateDiagrams[i], inputDiagrams[i + nDiags - 1]);
  //   if(max_dimension < 0.0) {
  //     this->printErr("Could not read Persistence Diagram");
  //     return 0;
  //   }
  //   if(max_dimension_total < max_dimension) {
  //     max_dimension_total = max_dimension;
  //   }
  // }

  std::vector<ttk::DiagramType> dictDiagrams;
  const int seed = this->Getseed_();

  // this->printMsg("==============COUCHE TTK=======================");

  std::vector<std::vector<double>> vectorWeights(nDiags);
  for(size_t i = 0; i < vectorWeights.size(); ++i) {
    // std::vector<double> weights{0.333, 0.333, 0.334};
    // std::vector<double> weights{1. / 3., 1. / 3., 1. / 3.};
    std::vector<double> weights(numAtom, 1. / (numAtom * 1.));
    vectorWeights[i] = std::move(weights);
  }

  std::vector<double> loss_tab;
  std::vector<double> true_loss_tab;
  std::vector<double> timers;
  std::vector<std::vector<double>> allLosses(nDiags);
  // const auto diagramsDistMat = this->execute(intermediateDiagrams,
  // dictDiagrams, vectorWeights,  nInputs);
  this->execute(intermediateDiagrams, intermediateAtoms, dictDiagrams,
                vectorWeights, nInputs, seed, numAtom, loss_tab, timers, true_loss_tab,
                allLosses, this->percent_);
  // zero-padd column name to keep Row Data columns ordered
  // this->printMsg("============WE ARE HERE 173 AFTER EXECUTE============");
  output_weights->SetNumberOfRows(nDiags);

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };
  // output_weights->SetNumberOfTuples(3);
  // this->printMsg("============WE ARE HERE 184 AFTER EXECUTE============");
  for(int i = 0; i < numAtom; ++i) {
    std::string name{"Atom"};
    zeroPad(name, numAtom, i);
    // name
    vtkNew<vtkDoubleArray> col{};
    // vtkDoubleArray *col=vtkDoubleArray::New();
    // col->SetNumberOfComponents(1);
    // col->SetNumberOfTuples(3);
    col->SetNumberOfValues(nDiags);
    col->SetName(name.c_str());
    for(int j = 0; j < nDiags; ++j) {
      col->SetValue(j, vectorWeights[j][i]);
    }
    col->Modified();
    // col->Modified();
    // printf("number of values %d", int(col->GetNumberOfValues()));
    output_weights->AddColumn(col);
  }

  vtkNew<vtkFieldData> fd{};
  fd->CopyStructure(inputDiagrams[0]->GetFieldData());
  fd->SetNumberOfTuples(nDiags);
  for(int i = 0; i < nDiags; ++i) {
    fd->SetTuple(i, 0, inputDiagrams[i]->GetFieldData());
  }

  for(int i = 0; i < fd->GetNumberOfArrays(); ++i) {
    output_weights->AddColumn(fd->GetAbstractArray(i));
  }

  vtkNew<vtkDoubleArray> colLoss{};
  colLoss->SetNumberOfValues(loss_tab.size());
  colLoss->SetName("Loss evolution");
  for(size_t j = 0; j < loss_tab.size(); ++j) {
    colLoss->SetValue(j, loss_tab[j]);
  }
  colLoss->Modified();
  output_loss->AddColumn(colLoss);

  vtkNew<vtkDoubleArray> trueColLoss{};
  trueColLoss->SetNumberOfValues(true_loss_tab.size());
  trueColLoss->SetName("True loss evolution");
  for(size_t j = 0; j < true_loss_tab.size(); ++j) {
    trueColLoss->SetValue(j, true_loss_tab[j]);
  }
  trueColLoss->Modified();
  true_output_loss->AddColumn(trueColLoss);

  for(int i = 0; i < nDiags; ++i) {
    std::vector<double> &loss = allLosses[i];
    std::string name{"Loss squared"};
    zeroPad(name, nDiags, i);
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(loss.size());
    col->SetName(name.c_str());
    for(size_t j = 0; j < loss.size(); ++j) {
      col->SetValue(j, loss[j]);
    }
    col->Modified();
    output_allLosses->AddColumn(col);
  }

  vtkNew<vtkDoubleArray> colTimers{};
  colTimers->SetNumberOfValues(timers.size());
  colTimers->SetName("Timers");
  for(size_t j = 0; j < timers.size() ; ++j){
    colTimers->SetValue(j , timers[j]);
  }
  colTimers->Modified();
  output_timers->AddColumn(colTimers);

  // this->printMsg("============WE ARE HERE 204 AFTER EXECUTE============");

  vtkNew<vtkFloatArray> dummy{};

  for(int i = 0; i < numAtom; ++i) {
    // vtkUnstructuredGrid temp =
    // vtkUnstructuredGrid::SafeDownCast(output_dgm->GetBlock(i));
    vtkNew<vtkUnstructuredGrid> vtu;
    ttk::DiagramType &diagram = dictDiagrams[i];
    DiagramToVTU(vtu, diagram, dummy, *this, 3, false);
    // this->printMsg("=====HERE?======");
    output_dgm->SetBlock(i, vtu);
    // this->printMsg("=====HERE2?=====");
  }
  // this->printMsg("========JUST BEFORE RETURN============");
  return 1;
}
