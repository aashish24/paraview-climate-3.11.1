#include "vtkDataSet.h"
#include "vtkDummyController.h"
#include "vtkMOCReader.h"
#include "vtkNew.h"
#include "vtkPointData.h"

int main()
{
  vtkNew<vtkDummyController> controller;
  controller->SetGlobalController(controller.GetPointer());
  vtkNew<vtkMOCReader> reader;
  reader->SetFileName("/home/acbauer/DATA/ParaViewData/Data/MHT/mhttest1/no_pbc_in_msf.mht");
  reader->Update();

  vtkDataSet* mhtGrid = vtkDataSet::SafeDownCast(reader->GetOutput(1));
  if(mhtGrid == NULL)
    {
    vtkGenericWarningMacro("Cannot find MHT output grid.");
    return 1;
    }
  if(mhtGrid->GetNumberOfPoints() != 301 ||
     mhtGrid->GetNumberOfCells() != 300)
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
  if(sum < 84.1 || sum > 84.2 || sum != sum)
    {
    vtkGenericWarningMacro("Global MHT sum " << sum << " but should be 84.1643.");
    return 1;
    }
  if(squaredSum < 279. || squaredSum > 280. || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Global MHT squared sum " << squaredSum << " but should be 279.805.");
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
  if(sum < 150. || sum > 151. || sum != sum)
    {
    vtkGenericWarningMacro("Atlantic MHT sum " << sum << " but should be 150.715.");
    return 1;
    }
  if(squaredSum < 121. || squaredSum > 122. || squaredSum != squaredSum)
    {
    vtkGenericWarningMacro("Atlantic MHT squared sum " << squaredSum << " but should be 121.712.");
    return 1;
    }

  controller->SetGlobalController(0);

  return 0;
}
