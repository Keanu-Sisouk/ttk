#include <ttkPersistenceDiagramDictEncoding.h>

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
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

int ttkPersistenceDiagramDictEncoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    //info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
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
    } else if (port == 1){
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
  //std::vector<vtkUnstructuredGrid *> inputDiagrams;

  //auto nBlocks = inputVector[0]->GetNumberOfInformationObjects();
  //std::vector<vtkMultiBlockDataSet *> blocks(nBlocks);

  //if(nBlocks > 2) {
  //  this->printWrn("Only dealing with the first two MultiBlockDataSets");
  //  nBlocks = 2;
  //}

   //number of diagrams per input block
  std::array<size_t, 2> nInputs{0, 0};

  //for(int i = 0; i < nBlocks; ++i) {
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

  // Flat storage for diagrams extracted from blocks
  std::vector<vtkUnstructuredGrid *> inputDiagrams;

  // Number of input diagrams
  int numInputs = 0;
  int numAtom = this->GetatomNumber_();
  printf("Atom number %d" , numAtom);

  if(blocks != nullptr) {
    numInputs = blocks->GetNumberOfBlocks();
    inputDiagrams.resize(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      inputDiagrams[i] = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(i));
      // if(this->GetMTime() < input[i]->GetMTime()) {
      //  needUpdate_ = true;
      //}
    }
  }

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
  //auto diagramsDistTable = vtkTable::GetData(outputVector);

  auto output_dgm = vtkMultiBlockDataSet::GetData(outputVector , 0);
  auto output_weights = vtkTable::GetData(outputVector , 1 );

  //int numAtom = this->GetAtomNumber();
  output_dgm->SetNumberOfBlocks(numAtom);

  for(int i = 0; i < numAtom; ++i) {
    vtkNew<vtkUnstructuredGrid> vtu;
    vtu->DeepCopy(inputDiagrams[i]);
    output_dgm->SetBlock(i, vtu);
  }

  std::vector<ttk::Diagram> intermediateDiagrams(nDiags);

  double max_dimension_total = 0.0;
  for(int i = 0; i < nDiags; ++i) {
    double max_dimension
      = getPersistenceDiagram(intermediateDiagrams[i], inputDiagrams[i]);
    if(max_dimension < 0.0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
    if(max_dimension_total < max_dimension) {
      max_dimension_total = max_dimension;
    }
  }

  std::vector<ttk::Diagram> dictDiagrams(numAtom);
  double max_dimension_total2 = 0.0;
  for(int i = 0; i < numAtom; ++i) {
    double max_dimension2
      = getPersistenceDiagram(dictDiagrams[i], inputDiagrams[i]);
    if(max_dimension2 < 0.0) {
      this->printErr("Could not read Persistence Diagram");
      return 0;
    }
    if(max_dimension_total2 < max_dimension2) {
      max_dimension_total2 = max_dimension2;
    }
  }

  std::vector<std::vector<double>> vectorWeights(nDiags);
  for (int i = 0; i < vectorWeights.size() ; ++i){
    std::vector<double> weights(numAtom , 1/3);
    vectorWeights[i] = weights;
  }




  const auto diagramsDistMat = this->execute(intermediateDiagrams, nInputs);

  // zero-padd column name to keep Row Data columns ordered
  



  output_weights->SetNumberOfRows(numAtom);

  const auto zeroPad
    = [](std::string &colName, const size_t numberCols, const size_t colIdx) {
        std::string max{std::to_string(numberCols - 1)};
        std::string cur{std::to_string(colIdx)};
        std::string zer(max.size() - cur.size(), '0');
        colName.append(zer).append(cur);
      };
  //output_weights->SetNumberOfTuples(3);
  for(int i=0 ; i < numInputs ; ++i){
    std::string name{"weights"};
    zeroPad(name, i, i);
    // name
    vtkNew<vtkDoubleArray> col{};
    // vtkDoubleArray *col=vtkDoubleArray::New();
    // col->SetNumberOfComponents(1);
    // col->SetNumberOfTuples(3);
    col->SetNumberOfValues(numAtom);
    col->SetName(name.c_str());
    for(int j = 0; j < numAtom; ++j){
      col->SetValue(j , 0.3);
    }
    col->Modified();
    //col->Modified();
    printf("number of values %d" , int(col->GetNumberOfValues()));
    output_weights->AddColumn(col);
  }


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

  diagram.resize(pairingsSize + 1);
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize; ++i) {

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

        diagram[0] = std::make_tuple(
          vertexId1, ttk::CriticalType::Local_minimum, vertexId2,
          ttk::CriticalType::Saddle1, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);
        diagram[pairingsSize] = std::make_tuple(
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
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }

  return max_dimension;
}
