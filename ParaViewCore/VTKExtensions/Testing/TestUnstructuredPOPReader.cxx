#include "vtkCellData.h"
#include "vtkCommunicator.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkIntegrateAttributes.h"
#include "vtkMPIController.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredPOPReader.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include "mpi.h"

// returns 0 for success
int TestReader(int stride, const int* voi, vtkIdType numberOfCells,
               double tData, double sData)
{
  int retVal = 0;
  vtkMultiProcessController* controller = vtkMPIController::GetGlobalController();
  vtkNew<vtkUnstructuredPOPReader> reader;
  reader->SetFileName("/home/acbauer/DATA/UVCDAT/TEMP.t.t0.1_42l_oilspill12c.00060101.pop.nc");
  reader->SetStride(stride, stride, stride);
  if(voi)
    {
    reader->SetVOI(voi[0], voi[1], voi[2], voi[3], voi[4], voi[5]);
    }
  vtkNew<vtkIntegrateAttributes> integrate;
  integrate->SetInputConnection(reader->GetOutputPort());
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetNumberOfPieces(controller->GetNumberOfProcesses());
  writer->SetWritePiece(controller->GetLocalProcessId());
  writer->SetInputConnection(integrate->GetOutputPort());
  writer->SetFileName("junk.vtu");
  writer->Update();

  vtkDataSet* grid = reader->GetOutput();
  vtkIdType tmp = grid->GetNumberOfCells();
  vtkIdType numCells = -1;
  controller->Reduce(&tmp, &numCells, 1, vtkCommunicator::SUM_OP, 0);
  grid = integrate->GetOutput();
  if(grid->GetNumberOfPoints() > 0 && grid->GetNumberOfCells() > 0)
    {
    if(numCells != numberOfCells)
      {
      vtkGenericWarningMacro( "I have " << grid->GetNumberOfCells() << " cells but should have "
                              << numberOfCells);
      retVal++;
      }

    if(vtkDataArray* temperature = grid->GetPointData()->GetArray("TEMP"))
      {
      double value = temperature->GetTuple1(0);
      // value is a negative number so the < and > may seem a bit funky
      if(value > tData*0.99999 || value < tData*1.00001)
        {
        vtkGenericWarningMacro("Bad temperature value of " << value << " which should be "
                               << tData);
        retVal++;
        }
      }
    vtkDataArray* size = grid->GetCellData()->GetArray("Area");
    if(!size)
      {
      size = grid->GetCellData()->GetArray("Volume");
      }
    double value = size->GetTuple1(0);
    if(value < sData*0.99999 || value > sData*1.00001)
      {
      vtkGenericWarningMacro("Bad " << size->GetName() << " " << value << " which should be "
                             << sData);
      retVal++;
      }
    }
  return retVal;
}

int main(int argc, char**argv)
{
  //---------------------------------------------------------------------------
  // Initialize MPI.

  // This is here to avoid false leak messages from vtkDebugLeaks when
  // using mpich. It appears that the root process which spawns all the
  // main processes waits in MPI_Init() and calls exit() when
  // the others are done, causing apparent memory leaks for any objects
  // created before MPI_Init().
  MPI_Init(&argc, &argv);

  vtkSmartPointer<vtkMPIController> controller =
    vtkSmartPointer<vtkMPIController>::New();
  controller->Initialize(&argc, &argv, 1);
  controller->SetGlobalController(controller); // this may not be necessary
  // Get information about the group of processes involved.
  int myId = controller->GetLocalProcessId();
  if(myId==0)
    {
    cerr << "=================== FIRST TEST ==========================\n";
    }
  int retVal = 0;
    {
    double tData = -2.44361e+19;
    double sData = 2.83211e+18;
    vtkIdType numberOfCells = 344160;
    int voi[6] = {0, -1, 0, -1, 0, -1};
    retVal += TestReader(10, voi, numberOfCells, tData, sData);
    }
  if(myId==0)
    {
    cerr << "=================== SECOND TEST ==========================\n";
    }
    {
    double tData = -4.60517e+11;
    double sData = 1.02864e+11;
    vtkIdType numberOfCells = 956;
    int voi[6] = {5, 5, 0, -1, 0, -1};
    retVal += TestReader(10, voi, numberOfCells, tData, sData);
    }
  if(myId==0)
    {
    cerr << "=================== THIRD TEST ==========================\n";
    }
    {
    double tData = -1.0303e+12;
    double sData = 5.30292e+10;
    vtkIdType numberOfCells = 1440;
    int voi[6] = {0, -1, 5, 5, 0, -1};
    retVal += TestReader(10, voi, numberOfCells, tData, sData);
    }
  if(myId==0)
    {
    cerr << "=================== FOURTH TEST ==========================\n";
    }
    {
    double tData = -1.83029e+15;
    double sData = 5.0429e+14;
    vtkIdType numberOfCells = 86040;
    int voi[6] = {0, -1, 0, -1, 2, 2};
    retVal += TestReader(10, voi, numberOfCells, tData, sData);
    }

  controller->Finalize();
  return retVal;
}
