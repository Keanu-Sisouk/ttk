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
//#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
//#include <vtkTable.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPersistenceDiagramDictDecoding);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkPersistenceDiagramDictDecoding::ttkPersistenceDiagramDictDecoding() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

// ttkPersistenceDiagramDictDecoding::~ttkPersistenceDiagramDictDecoding() {
// }

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkPersistenceDiagramDictDecoding::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    /*info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);*/
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkPersistenceDiagramDictDecoding::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    // info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  } else if (port == 1){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  } else {
    return 0;
  }
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
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

  std::vector<ttk::Diagram> dictDiagrams(nDiags);
  double max_dimension_total2 = 0.0;
  for(size_t i = 0; i < nDiags; ++i) {
    ttk::Diagram &atom = dictDiagrams[i];
    // double max_dimension2 = getPersistenceDiagram(
    // atom, vtkUnstructuredGrid::SafeDownCast(inputDiagrams[i]));
    double max_dimension2 = getPersistenceDiagram(atom, inputDiagrams[i]);
    // for(size_t k = 0; k < atom.size(); ++k) {
    //   DiagramTuple &t = atom[k];
    //   std::cout << "Pair atoms: " << std::get<6>(t) << ", " <<
    //   std::get<10>(t)
    //             << std::endl;
    // }
    if(max_dimension2 < 0.0) {
      this->printWrn("Could not read Persistence Diagram");
      // return 0;
    }
    if(max_dimension_total2 < max_dimension2) {
      max_dimension_total2 = max_dimension2;
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
    = [](std::string &colName, const size_t numberCols, const size_t colIdx){
      std::string max{std::to_string(numberCols - 1)};
      std::string cur{std::to_string(colIdx)};
      std::string zer(max.size() -cur.size(), '0');
      colName.append(zer).append(cur);
    };

  std::cout << "PASSED !!!!" << std::endl;
  std::vector<vtkDataArray *> inputWeights;
  int numWeights = weights_vtk->GetNumberOfRows();
  // this->printMsg(std::to_string(numWeights));
  if(weights_vtk != nullptr) {
    // int numWeights = weights_vtk->GetNumberOfColumns();
    // this->printMsg(std::to_string(numWeights));
    inputWeights.resize(nDiags);
    for(size_t i = 0; i < nDiags; ++i) {
      std::string name{"Atom"};
      zeroPad(name, numWeights, i);
      inputWeights[i] = vtkDataArray::SafeDownCast(weights_vtk->GetColumn(i));
      // if(this->GetMTime() < input[i]->GetMTime()) {
      //  needUpdate_ = true;
      //}
    }
  }

  std::cout << "PASSED 2 !!!!!" << std::endl;

  //const int nWeights = ;
  std::vector<std::vector<double>> vectorWeights(numWeights);
  for(int i = 0; i < numWeights; ++i) {
    std::vector<double> &t1 = vectorWeights[i];
    // vtkDoubleArray &t2 = inputWeights[i];
    for(int j = 0; j < nDiags; ++j) {
      //double weight = t1[j];
      std::cout << "ICI???" << std::endl;
      double weight = inputWeights[j]->GetTuple1(i);
      t1.push_back(weight);
    }
  }
  std::cout << "PASSED 3!!!!!!" << std::endl;
  std::vector<ttk::Diagram> Barycenters(numWeights);


  if(!ComputePoints){
    this->execute(dictDiagrams, vectorWeights, Barycenters);
  }
  // this->printMsg("=====ICI?======");
  auto output_dgm = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto output_coordinates = vtkTable::GetData(outputVector, 1);
  output_dgm->SetNumberOfBlocks(numWeights);
  int dim = 2;
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


  outputDiagrams(output_dgm, output_coordinates, Barycenters, dictDiagrams, vectorWeights,Spacing,1);
  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation

  // make a SHALLOW copy of the input
  // outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  // outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}

double ttkPersistenceDiagramDictDecoding::getPersistenceDiagram(
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

void ttkPersistenceDiagramDictDecoding::outputDiagrams(
  vtkMultiBlockDataSet *output,
  vtkTable *output_coordinates,
  const std::vector<ttk::Diagram> &diags,
  const std::vector<ttk::Diagram> &atoms,
  const std::vector<std::vector<double>> &weights,
  const double spacing,
  const double max_persistence) const {

  ttk::SimplexId nDiags = diags.size();
  ttk::SimplexId nAtoms = atoms.size();

  ttk::SimplexId n_existing_blocks = ShowAtoms ? nAtoms : 0;

  output->SetNumberOfBlocks(nDiags+n_existing_blocks);
  std::vector<std::pair<double,double>> coords(nAtoms);
  std::vector<std::pair<double,double>> true_coords(nAtoms);

  if(nAtoms == 3){
    ttk::PersistenceDiagramDistanceMatrix MatrixCalculator;
    std::array<size_t , 2> nInputs{nAtoms, 0};
    MatrixCalculator.setDos(true, true, true);
    MatrixCalculator.setThreadNumber(3);
    std::vector<std::vector<double>> distMatrix = MatrixCalculator.execute(atoms, nInputs);
    coords[0].first = 0.;
    true_coords[0].first = 0.;
    coords[0].second = 0.;
    true_coords[0].second = 0.;
    coords[1].first = spacing  * distMatrix[0][1];
    true_coords[1].first = distMatrix[0][1];
    coords[1].second = 0.;
    true_coords[0].second = 0.;
    double distOpposed = distMatrix[2][1];
    double firstDist = distMatrix[0][1];
    double distAdja = distMatrix[0][2];
    double alpha = std::acos((distOpposed * distOpposed -firstDist * firstDist -distAdja * distAdja)/(-2. * firstDist * distAdja));
    coords[2].first = spacing * distAdja * std::cos(alpha);
    true_coords[2].first = distAdja * std::cos(alpha);
    coords[2].second = spacing * distAdja * std::sin(alpha);
    true_coords[2].second = distAdja * std::sin(alpha);
    
    if(ShowAtoms){
      for(size_t i = 0 ; i < nAtoms ; ++i){
        double X = coords[i].first;
        double Y = coords[i].second;
        vtkNew<vtkUnstructuredGrid> vtu{};
        this->diagramToVTU(vtu, atoms[i], max_persistence);

        vtkNew<vtkTransform> tr{};
        tr->Translate(X,Y,0);

        vtkNew<vtkTransformFilter> trf{};
        trf->SetTransform(tr);
        trf->SetInputData(vtu);
        trf->Update();

        output->SetBlock(i, trf->GetOutputDataObject(0));
      }
    }

  } else {

    for(size_t i = 0; i < nAtoms; ++i) {
      const auto angle = 2.0 * M_PI * static_cast<double>(i)
        / static_cast<double>(nAtoms);
      double X = spacing * max_persistence * std::cos(angle);
      double Y = spacing * max_persistence * std::sin(angle);
      coords[i].first = X;
      coords[i].second = Y;
      true_coords[i].first = max_persistence * std::cos(angle);
      true_coords[i].second = max_persistence * std::sin(angle);

      if(ShowAtoms){
        vtkNew<vtkUnstructuredGrid> vtu{};
        this->diagramToVTU(vtu, atoms[i], max_persistence);

        vtkNew<vtkTransform> tr{};
        tr->Translate(X,Y,0);

        vtkNew<vtkTransformFilter> trf{};
        trf->SetTransform(tr);
        trf->SetInputData(vtu);
        trf->Update();

        output->SetBlock(i, trf->GetOutputDataObject(0));
      }
    }
  }




  for(size_t i = 0; i < diags.size(); ++i) {
    vtkNew<vtkUnstructuredGrid> vtu{};
    if(!ComputePoints){
      this->diagramToVTU(vtu, diags[i], max_persistence);
    }

      double X = 0;
      double Y = 0;
      for(size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        X += weights[i][iAtom]*coords[iAtom].first;
        Y += weights[i][iAtom]*coords[iAtom].second;
      }

      vtkNew<vtkTransform> tr{};
      tr->Translate(X,Y,0);
      vtkNew<vtkTransformFilter> trf{};
      trf->SetTransform(tr);
      trf->SetInputData(vtu);
      trf->Update();

      output->SetBlock(i+n_existing_blocks, trf->GetOutputDataObject(0));
  }


  for(size_t i = 0 ; i < 2 ; ++i){
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfValues(nDiags);
    std::string name;
    std::cout << "NANI THE FUCK" << std::endl;
    if(i == 0){
      name = "X";
    } else {
      name = "Y";
    }
    //name.append(std::to_string(i));
    col->SetName(name.c_str());
    for(size_t j = 0 ; j < nDiags ; ++j){
      double temp = 0;
      for(size_t iAtom = 0; iAtom < nAtoms ; ++iAtom){
        if(i == 0){
          temp += weights[j][iAtom]*true_coords[iAtom].first;
        } else {
          temp += weights[j][iAtom]*true_coords[iAtom].second;
        }
      }
      col->SetValue(j, temp);
    }
    col->Modified();
    output_coordinates->AddColumn(col);
  }


}

void ttkPersistenceDiagramDictDecoding::diagramToVTU(
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
  ttkPersistenceDiagramDictDecoding::getMaxPersistence(ttk::Diagram &diagram) {
  double max_persistence{0};
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &t = diagram[i];
    const double pers = std::get<4>(t);
    max_persistence = std::max(pers, max_persistence);
  }
  return max_persistence;
}
