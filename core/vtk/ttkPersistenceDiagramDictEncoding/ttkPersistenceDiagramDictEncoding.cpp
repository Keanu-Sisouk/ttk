#include <ttkMacros.h>
#include <ttkPersistenceDiagramDictEncoding.h>
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
  SetNumberOfOutputPorts(4);
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

  int numInputAtoms =0;
  if(atomBlocks != nullptr) {
    numInputAtoms = atomBlocks->GetNumberOfBlocks();
    inputAtoms.resize(numInputAtoms);
    for(int i = 0; i < numInputAtoms; ++i) {
      inputAtoms[i] = vtkUnstructuredGrid::SafeDownCast(atomBlocks->GetBlock(i));
    }
  }
  if(BackEnd == BACKEND::INPUT_ATOMS){
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
  printMsg("Number of atoms: "+ttk::debug::output::YELLOW+ttk::debug::output::UNDERLINED+std::to_string(numAtom)+ttk::debug::output::ENDCOLOR+ttk::debug::output::ENDCOLOR);

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

  // int numAtom = this->GetAtomNumber();
  output_dgm->SetNumberOfBlocks(numAtom);

  for(int i = 0; i < numAtom; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu;
    vtu->DeepCopy(inputDiagrams[i]);
    output_dgm->SetBlock(i, vtu);
  }

  std::vector<ttk::Diagram> intermediateDiagrams(nDiags);
  std::vector<ttk::Diagram> intermediateAtoms(numInputAtoms);

  double max_dimension_total = 0.0;
  double percentage = static_cast<double>(this->percent_);
  for(int i = 0; i < nDiags; ++i) {
    double max_dimension
      = getPersistenceDiagram(intermediateDiagrams[i], inputDiagrams[i]);
    std::cout << "MAX PERS BEFORE FILTERING " << max_dimension << std::endl;

    intermediateDiagrams[i].erase(
      std::remove_if(intermediateDiagrams[i].begin(),
                     intermediateDiagrams[i].end(),
                     [max_dimension, percentage](ttk::DiagramTuple &t) {
                       return (std::get<10>(t) - std::get<6>(t))
                              < (percentage / 100.) * max_dimension;
                     }),
      intermediateDiagrams[i].end());

    auto &t = intermediateDiagrams[i][0];
    std::cout << "MAX PERS BEFORE ORDERING " << std::get<10>(t) - std::get<6>(t) << std::endl;
    //std::sort(intermediateDiagrams[i].begin(), intermediateDiagrams[i].end(),
    //          [](ttk::DiagramTuple &t1, ttk::DiagramTuple &t2) {
    //            return (std::get<10>(t1) - std::get<6>(t1))
    //                   < (std::get<10>(t2) - std::get<6>(t2));
    //          });
    if(max_dimension < 0.0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
    if(max_dimension_total < max_dimension) {
      max_dimension_total = max_dimension;
    }
  }
  for(int i = 0; i < numInputAtoms; ++i) {
    double max_dimension
      = getPersistenceDiagram(intermediateAtoms[i], inputAtoms[i]);
    if(max_dimension < 0.0) {
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

  std::vector<ttk::Diagram> dictDiagrams;
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
  std::vector<std::vector<double>> allLosses(nDiags);
  // const auto diagramsDistMat = this->execute(intermediateDiagrams,
  // dictDiagrams, vectorWeights,  nInputs);
  this->execute(intermediateDiagrams, intermediateAtoms, dictDiagrams,
                vectorWeights, nInputs, seed, numAtom, loss_tab, allLosses, this->percent_);
  // zero-padd column name to keep Row Data columns ordered
  // this->printMsg("============WE ARE HERE 173 AFTER EXECUTE============");
  output_weights->SetNumberOfRows(numAtom);

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };
  // output_weights->SetNumberOfTuples(3);
  // this->printMsg("============WE ARE HERE 184 AFTER EXECUTE============");
  for(int i = 0; i < nDiags; ++i) {
    std::string name{"weights"};
    zeroPad(name, nDiags, i);
    // name
    vtkNew<vtkDoubleArray> col{};
    // vtkDoubleArray *col=vtkDoubleArray::New();
    // col->SetNumberOfComponents(1);
    // col->SetNumberOfTuples(3);
    col->SetNumberOfValues(numAtom);
    col->SetName(name.c_str());
    for(int j = 0; j < numAtom; ++j) {
      col->SetValue(j, vectorWeights[i][j]);
    }
    col->Modified();
    // col->Modified();
    // printf("number of values %d", int(col->GetNumberOfValues()));
    output_weights->AddColumn(col);
  }

  vtkNew<vtkDoubleArray> colLoss{};
  colLoss->SetNumberOfValues(loss_tab.size());
  colLoss->SetName("Loss evolution");
  for(size_t j = 0; j < loss_tab.size(); ++j) {
    colLoss->SetValue(j, loss_tab[j]);
  }
  colLoss->Modified();
  output_loss->AddColumn(colLoss);

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
  // this->printMsg("============WE ARE HERE 204 AFTER EXECUTE============");

  for(int i = 0; i < numAtom; ++i) {
    // vtkUnstructuredGrid temp =
    // vtkUnstructuredGrid::SafeDownCast(output_dgm->GetBlock(i));
    vtkNew<vtkUnstructuredGrid> vtu;
    ttk::Diagram &diagram = dictDiagrams[i];
    double max_persistence = getMaxPersistence(diagram);
    diagramToVTU(vtu, diagram, max_persistence);
    // this->printMsg("=====HERE?======");
    output_dgm->SetBlock(i, vtu);
    // this->printMsg("=====HERE2?=====");
  }
  // this->printMsg("========JUST BEFORE RETURN============");
  return 1;
}

double ttkPersistenceDiagramDictEncoding::getPersistenceDiagram(
  ttk::Diagram &diagram, vtkUnstructuredGrid *CTPersistenceDiagram_) {

  const auto pd = CTPersistenceDiagram_->GetPointData();
  const auto cd = CTPersistenceDiagram_->GetCellData();
  const auto points = CTPersistenceDiagram_->GetPoints();

  if(pd == nullptr || cd == nullptr || points == nullptr) {
    this->printErr("Missing Diagram PointData, CellData or Points");
    return -1.0;
  }

  const auto vertexIdentifierScalars
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pd->GetArray("CriticalType"));
  const auto pairIdentifierScalars
    = vtkIntArray::SafeDownCast(cd->GetArray("PairIdentifier"));
  const auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cd->GetArray("PairType"));
  const auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cd->GetArray("Persistence"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Death"));
  const auto critCoordinates
    = vtkFloatArray::SafeDownCast(pd->GetArray("Coordinates"));

  const bool embed = birthScalars != nullptr && deathScalars != nullptr;

  if(!embed && critCoordinates == nullptr) {
    this->printErr("Malformed Persistence Diagram");
    return -2.0;
  }

  int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
      pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }

  // If diagram has the diagonal (we assume it is last)
  if(*pairIdentifierScalars->GetTuple(pairingsSize - 1) == -1)
    pairingsSize -= 1;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points) {
    this->printErr("Missing Persistence Diagram data array");
    return -3.0;
  }

  diagram.resize(pairingsSize);
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize; ++i) {
    // this->printMsg("=====" + std::to_string(i) + "=====DEBUT=====");
    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    std::array<double, 3> coordsBirth{}, coordsDeath{};

    const auto i0 = 2 * i;
    const auto i1 = 2 * i + 1;

    double birth, death;

    if(embed) {
      points->GetPoint(i0, coordsBirth.data());
      points->GetPoint(i1, coordsDeath.data());
      birth = birthScalars->GetValue(i0);
      death = deathScalars->GetValue(i1);
    } else {
      critCoordinates->GetTuple(i0, coordsBirth.data());
      critCoordinates->GetTuple(i1, coordsDeath.data());
      birth = points->GetPoint(i0)[0];
      death = points->GetPoint(i1)[1];
    }

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
      if(pairIdentifier == 0) {
        max_dimension = persistence;

        // diagram[0] = std::make_tuple(
        // vertexId1, ttk::CriticalType::Local_minimum, vertexId2,
        // ttk::CriticalType::Saddle1, persistence, pairType, birth,
        // coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
        // coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        diagram[0] = std::make_tuple(
          vertexId1, ttk::CriticalType::Saddle1, vertexId2,
          ttk::CriticalType::Local_maximum, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);

      } else {
        diagram[pairIdentifier] = std::make_tuple(
          vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
          persistence, pairType, birth, coordsBirth[0], coordsBirth[1],
          coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
      }
    }
    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        this->printWrn("Diagram pair identifiers must be compact (not exceed "
                       "the diagram size).");
      }
    }
    // this->printMsg("=====" + std::to_string(i) + "=====FIN=====");
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }

  return max_dimension;
}

// double ttkPersistenceDiagramDictEncoding::getPersistenceDiagram(
//   ttk::Diagram &diagram, vtkUnstructuredGrid *CTPersistenceDiagram_) {
//
//   const auto pd = CTPersistenceDiagram_->GetPointData();
//   const auto cd = CTPersistenceDiagram_->GetCellData();
//   const auto points = CTPersistenceDiagram_->GetPoints();
//
//   if(pd == nullptr || cd == nullptr || points == nullptr) {
//     this->printErr("Missing Diagram PointData, CellData or Points");
//     return -1.0;
//   }
//
//   const auto vertexIdentifierScalars
//     = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
//   const auto nodeTypeScalars
//     = vtkIntArray::SafeDownCast(pd->GetArray("CriticalType"));
//   const auto pairIdentifierScalars
//     = vtkIntArray::SafeDownCast(cd->GetArray("PairIdentifier"));
//   const auto extremumIndexScalars
//     = vtkIntArray::SafeDownCast(cd->GetArray("PairType"));
//   const auto persistenceScalars
//     = vtkDoubleArray::SafeDownCast(cd->GetArray("Persistence"));
//   const auto birthScalars =
//   vtkDoubleArray::SafeDownCast(pd->GetArray("Birth")); const auto
//   deathScalars = vtkDoubleArray::SafeDownCast(pd->GetArray("Death")); const
//   auto critCoordinates
//     = vtkFloatArray::SafeDownCast(pd->GetArray("Coordinates"));
//
//   const bool embed = birthScalars != nullptr && deathScalars != nullptr;
//
//   if(!embed && critCoordinates == nullptr) {
//     this->printErr("Malformed Persistence Diagram");
//     return -2.0;
//   }
//
//   int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
//   // FIX : no more missed pairs
//   for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
//     const float index_of_pair = pair_index;
//     if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
//       pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
//   }
//
//   // If diagram has the diagonal (we assume it is last)
//   if(*pairIdentifierScalars->GetTuple(pairingsSize - 1) == -1)
//     pairingsSize -= 1;
//
//   if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
//      || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
//      || !points) {
//     this->printErr("Missing Persistence Diagram data array");
//     return -3.0;
//   }
//
//   diagram.resize(pairingsSize + 1);
//   int nbNonCompact = 0;
//   double max_dimension = 0;
//
//   // skip diagonal cell (corresponding points already dealt with)
//   for(int i = 0; i < pairingsSize; ++i) {
//
//     int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
//     int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
//     int nodeType1 = nodeTypeScalars->GetValue(2 * i);
//     int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);
//
//     int pairIdentifier = pairIdentifierScalars->GetValue(i);
//     int pairType = extremumIndexScalars->GetValue(i);
//     double persistence = persistenceScalars->GetValue(i);
//
//     std::array<double, 3> coordsBirth{}, coordsDeath{};
//
//     const auto i0 = 2 * i;
//     const auto i1 = 2 * i + 1;
//
//     double birth, death;
//
//     if(embed) {
//       points->GetPoint(i0, coordsBirth.data());
//       points->GetPoint(i1, coordsDeath.data());
//       birth = birthScalars->GetValue(i0);
//       death = deathScalars->GetValue(i1);
//     } else {
//       critCoordinates->GetTuple(i0, coordsBirth.data());
//       critCoordinates->GetTuple(i1, coordsDeath.data());
//       birth = points->GetPoint(i0)[0];
//       death = points->GetPoint(i1)[1];
//     }
//
//     if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
//       if(pairIdentifier == 0) {
//         max_dimension = persistence;
//
//         diagram[0] = std::make_tuple(
//           vertexId1, ttk::CriticalType::Local_minimum, vertexId2,
//           ttk::CriticalType::Saddle1, persistence, pairType, birth,
//           coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
//           coordsDeath[0], coordsDeath[1], coordsDeath[2]);
//         diagram[pairingsSize] = std::make_tuple(
//           vertexId1, ttk::CriticalType::Saddle1, vertexId2,
//           ttk::CriticalType::Local_maximum, persistence, pairType, birth,
//           coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
//           coordsDeath[0], coordsDeath[1], coordsDeath[2]);
//
//       } else {
//         diagram[pairIdentifier] = std::make_tuple(
//           vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
//           persistence, pairType, birth, coordsBirth[0], coordsBirth[1],
//           coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
//           coordsDeath[2]);
//       }
//     }
//     if(pairIdentifier >= pairingsSize) {
//       nbNonCompact++;
//       if(nbNonCompact == 0) {
//         this->printWrn("Diagram pair identifiers must be compact (not exceed
//         the diagram size).");
//       }
//     }
//   }
//
//   if(nbNonCompact > 0) {
//     this->printWrn("Missed " + std::to_string(nbNonCompact)
//                    + " pairs due to non-compactness.");
//   }
//
//   return max_dimension;
// }

void ttkPersistenceDiagramDictEncoding::diagramToVTU(
  vtkUnstructuredGrid *output,
  const ttk::Diagram &diagram,
  const double max_persistence) const {

  const auto nPoints = 2 * diagram.size();
  if(nPoints == 0) {
    this->printWrn("Diagram with no points");
    return;
  }

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(nPoints);
  output->SetPoints(points);

  // point data
  vtkNew<vtkIntArray> critType{};
  critType->SetName("CriticalType");
  critType->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(critType);

  vtkNew<vtkFloatArray> coords{};
  coords->SetNumberOfComponents(3);
  coords->SetName("Coordinates");
  coords->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(coords);

  vtkNew<vtkDoubleArray> pointPers{};
  pointPers->SetName("Persistence");
  pointPers->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(pointPers);

  vtkNew<ttkSimplexIdTypeArray> vsf{};
  vsf->SetName(ttk::VertexScalarFieldName);
  vsf->SetNumberOfTuples(nPoints);
  output->GetPointData()->AddArray(vsf);

  // cell data
  vtkNew<vtkIntArray> pairId{};
  pairId->SetName("PairIdentifier");
  pairId->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairId);

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");
  pairType->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairType);

  vtkNew<vtkDoubleArray> pairPers{};
  pairPers->SetName("Persistence");
  pairPers->SetNumberOfTuples(diagram.size() + 1);
  output->GetCellData()->AddArray(pairPers);

  for(size_t j = 0; j < diagram.size(); ++j) {
    const auto &pair{diagram[j]};
    const auto birth{std::get<6>(pair)};
    const auto death{std::get<10>(pair)};
    const auto birtVertId{std::get<0>(pair)};
    const auto deathVertId{std::get<2>(pair)};
    const auto pType{std::get<5>(pair)};
    const auto birthType{std::get<1>(pair)};
    const auto deathType{std::get<3>(pair)};
    std::array<float, 3> coordsBirth{
      std::get<7>(pair), std::get<8>(pair), std::get<9>(pair)};
    std::array<float, 3> coordsDeath{
      std::get<11>(pair), std::get<12>(pair), std::get<13>(pair)};

    // cell data
    pairId->SetTuple1(j, j);
    pairType->SetTuple1(j, pType);
    pairPers->SetTuple1(j, death - birth);

    // point data
    coords->SetTuple(2 * j + 0, coordsBirth.data());
    coords->SetTuple(2 * j + 1, coordsDeath.data());
    pointPers->SetTuple1(2 * j + 0, death - birth);
    pointPers->SetTuple1(2 * j + 1, death - birth);
    critType->SetTuple1(2 * j + 0, static_cast<int>(birthType));
    critType->SetTuple1(2 * j + 1, static_cast<int>(deathType));
    vsf->SetTuple1(2 * j + 0, birtVertId);
    vsf->SetTuple1(2 * j + 1, deathVertId);

    points->SetPoint(2 * j + 0, birth, birth, 0);
    points->SetPoint(2 * j + 1, birth, death, 0);

    const std::array<vtkIdType, 2> ids{
      2 * static_cast<vtkIdType>(j) + 0,
      2 * static_cast<vtkIdType>(j) + 1,
    };
    output->InsertNextCell(VTK_LINE, 2, ids.data());
  }

  // add diagonal
  const auto minmax_birth = std::minmax_element(
    diagram.begin(), diagram.end(),
    [](const ttk::DiagramTuple &a, const ttk::DiagramTuple &b) {
      return std::get<6>(a) < std::get<6>(b);
    });
  const std::array<vtkIdType, 2> ids{
    2 * (minmax_birth.first - diagram.begin()),
    2 * (minmax_birth.second - diagram.begin()),
  };
  output->InsertNextCell(VTK_LINE, 2, ids.data());
  pairId->SetTuple1(diagram.size(), -1);
  pairType->SetTuple1(diagram.size(), -1);
  // use twice the max persistence of all input diagrams...
  pairPers->SetTuple1(diagram.size(), 2.0 * max_persistence);
}

double
  ttkPersistenceDiagramDictEncoding::getMaxPersistence(ttk::Diagram &diagram) {
  double max_persistence{0};
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &t = diagram[i];
    const double pers = std::get<4>(t);
    max_persistence = std::max(pers, max_persistence);
  }
  return max_persistence;
}
