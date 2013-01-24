#include "vtkDataArray.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkDummyController.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiBlockTemporalStatistics.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkTimeSourceExample.h"

int main(int argc, char* argv[])
{
  int returnValue = 0;
  vtkNew<vtkDummyController> controller;
  controller->Initialize(&argc, &argv);
  controller->SetGlobalController(controller.GetPointer());

  vtkNew<vtkTimeSourceExample> source;
  vtkNew<vtkMultiBlockTemporalStatistics> statistics;
  statistics->SetInputConnection(source->GetOutputPort());
  statistics->Update();

  vtkMultiBlockDataSet* multiBlock =
    vtkMultiBlockDataSet::SafeDownCast(statistics->GetOutput());
  vtkDataSet* output = vtkDataSet::SafeDownCast(multiBlock->GetBlock(0));

  vtkDataArray* array = output->GetPointData()->GetArray("Point Value_average");
  if(array->GetTuple1(0) < -0.001 || array->GetTuple1(0) > 0.001)
    {
    vtkGenericWarningMacro("Bad average value of " << array->GetTuple1(0) << " which should be 0.");
    returnValue++;
    }

  array = output->GetPointData()->GetArray("Point Value_maximum");
  if(array->GetTuple1(0) < 0.9848 || array->GetTuple1(0) > 0.9849)
    {
    vtkGenericWarningMacro("Bad maximum value of " << array->GetTuple1(0) << " which should be 0.");
    returnValue++;
    }

  array = output->GetPointData()->GetArray("Point Value_minimum");
  if(array->GetTuple1(0) < -0.9849 || array->GetTuple1(0) > -0.9848)
    {
    vtkGenericWarningMacro("Bad minimum value of " << array->GetTuple1(0) << " which should be 0.");
    returnValue++;
    }

  array = output->GetPointData()->GetArray("Point Value_stddev");
  if(array->GetTuple1(0) < 0.67081 || array->GetTuple1(0) > 0.67083)
    {
    vtkGenericWarningMacro("Bad standard deviation value of " << array->GetTuple1(0) << " which should be 0.");
    returnValue++;
    }

  controller->SetGlobalController(0);
  controller->Finalize();

  return returnValue;
}
