#include "vtkDataSet.h"
#include "vtkDummyController.h"
#include "vtkMOCReader.h"
#include "vtkNew.h"
#include "vtkPointData.h"

#include <string>

int main(int argc, char* argv[])
{
  if(argc < 2)
    {
    vtkGenericWarningMacro("Must provide the UVCDAT test data directory.");
    return 1;
    }
  vtkNew<vtkDummyController> controller;
  controller->SetGlobalController(controller.GetPointer());
  std::string fileName = argv[1];
  fileName.append("/p94m/in_msf");
  cerr << fileName << " is the filename\n";
  vtkNew<vtkMOCReader> reader;
  reader->SetFileName(fileName.c_str());
  reader->Update();

  vtkDataSet* mocGrid = vtkDataSet::SafeDownCast(reader->GetOutput());
  if(mocGrid == NULL)
    {
    vtkGenericWarningMacro("Cannot find MOC output grid.");
    return 1;
    }
  if(mocGrid->GetNumberOfPoints() != 6440 ||
     mocGrid->GetNumberOfCells() != 6240)
    {
    vtkGenericWarningMacro("MOC grid has wrong number of points or cells");
    return 1;
    }

  // global MOC
  vtkDataArray* mocValues = mocGrid->GetPointData()->GetArray("reader_moc_global");
  double sum = 0;
  double squaredSum = 0;
  for(vtkIdType i=0;i<mocValues->GetNumberOfTuples();i++)
    {
    sum += mocValues->GetTuple1(i);
    squaredSum += mocValues->GetTuple1(i)*mocValues->GetTuple1(i);
    }
  if(sum < -1.42e+36 || sum > -1.4e+36 || sum != sum)
    {
    vtkGenericWarningMacro("Global MOC sum " << sum << " but should be -1.41e+36.");
    return 1;
    }
  if(squaredSum < 1.4e+70 || squaredSum > 1.42e+70 || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Global MOC squared sum " << squaredSum << " but should be 1.41e+70.");
    return 1;
    }

  // atlantic MOC
  mocValues = mocGrid->GetPointData()->GetArray("reader_moc_atl");
  sum = 0;
  squaredSum = 0;
  for(vtkIdType i=0;i<mocValues->GetNumberOfTuples();i++)
    {
    sum += mocValues->GetTuple1(i);
    squaredSum += mocValues->GetTuple1(i)*mocValues->GetTuple1(i);
    }
  if(sum < -1.84e+37 || sum > -1.8e+37 || sum != sum)
    {
    vtkGenericWarningMacro("Atlantic MOC sum " << sum << " but should be -1.811e+37.");
    return 1;
    }
  if(squaredSum < 1.8e+71 || squaredSum > 1.82e+71 || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Atlantic MOC squared sum " << squaredSum << " but should be 1.811e+71.");
    return 1;
    }

  controller->SetGlobalController(0);

  return 0;
}
