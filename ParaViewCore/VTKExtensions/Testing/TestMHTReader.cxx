#include "vtkDataSet.h"
#include "vtkDummyController.h"
#include "vtkMHTReader.h"
#include "vtkNew.h"
#include "vtkPointData.h"

#include <string>

int main(int argc, char* argv[])
{
  if(argc < 3)
    {
    vtkGenericWarningMacro("Must provide the UVCDAT test data directory.");
    return 1;
    }
  vtkNew<vtkDummyController> controller;
  controller->SetGlobalController(controller.GetPointer());
  std::string fileName = argv[2];
  fileName.append("/p94m/in_msf");
  cerr << fileName << " is the filename\n";
  vtkNew<vtkMHTReader> reader;
  reader->SetFileName(fileName.c_str());
  reader->Update();

  vtkDataSet* mhtGrid = vtkDataSet::SafeDownCast(reader->GetOutput());
  if(mhtGrid == NULL)
    {
    vtkGenericWarningMacro("Cannot find MHT output grid.");
    return 1;
    }
  if(mhtGrid->GetNumberOfPoints() != 161 ||
     mhtGrid->GetNumberOfCells() != 160)
    {
    vtkGenericWarningMacro("MHT grid has wrong number of points or cells");
    return 1;
    }

  // global MHT
  vtkDataArray* mhtValues = mhtGrid->GetPointData()->GetArray("reader_mht_global");
  double sum = 0;
  double squaredSum = 0;
  for(vtkIdType i=0;i<mhtValues->GetNumberOfTuples();i++)
    {
    sum += mhtValues->GetTuple1(i);
    squaredSum += mhtValues->GetTuple1(i)*mhtValues->GetTuple1(i);
    }
  if(sum < 1.47e+09 || sum > 1.49e+09 || sum != sum)
    {
    vtkGenericWarningMacro("Global MHT sum " << sum << " but should be 1.47997e+09.");
    return 1;
    }
  if(squaredSum < 1.65e+16 || squaredSum > 1.75e+16 || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Global MHT squared sum " << squaredSum << " but should be 1.69631e+16.");
    return 1;
    }

  // atlantic MHT
  mhtValues = mhtGrid->GetPointData()->GetArray("reader_mht_atl");
  sum = 0;
  squaredSum = 0;
  for(vtkIdType i=0;i<mhtValues->GetNumberOfTuples();i++)
    {
    sum += mhtValues->GetTuple1(i);
    squaredSum += mhtValues->GetTuple1(i)*mhtValues->GetTuple1(i);
    }
  if(sum < 2.5e+08 || sum > 2.6e+08 || sum != sum)
    {
    vtkGenericWarningMacro("Atlantic MHT sum " << sum << " but should be 2.54814e+08.");
    return 1;
    }
  if(squaredSum < 6.43e+14 || squaredSum > 6.44e+14 || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Atlantic MHT squared sum " << squaredSum << " but should be 6.4383e+14.");
    return 1;
    }

  controller->SetGlobalController(0);

  return 0;
}
