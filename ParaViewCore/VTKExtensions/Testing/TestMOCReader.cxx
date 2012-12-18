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

  vtkDataArray* mhtValues = mhtGrid->GetPointData()->GetArray("reader_mht_global");
  double sum = 0;
  for(vtkIdType i=0;i<mhtValues->GetNumberOfTuples();i++)
    {
    sum += mhtValues->GetTuple1(i);
    }
  if(sum < 84.1 || sum > 84.2)
    {
    vtkGenericWarningMacro("MHT sum " << sum << " but should be 84.1643.");
    return 1;
    }

  controller->SetGlobalController(0);

  return 0;
}
