/*=========================================================================

  Module:    vtkMHTReader.cxx

  =========================================================================*/

#include "vtkMHTReader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkExtentTranslator.h"
#include "vtkRectilinearGrid.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiProcessStream.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPOPMatrices.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include "vtkMultiProcessController.h"
#include "vtkCommunicator.h"

#include "vtkNew.h"
#include "vtkTimerLog.h"

int LoadData3D(
  std::string& fileName, POPInputInformation& popInfo, int* ext3D1GL,
  int imt1GL, int jmt1GL, int km, Matrix3DFloat& array)
{
  int retVal = vtkAbstractPOPReader::LoadDataBlock3DFloat(&popInfo, fileName, ext3D1GL,
                                                          imt1GL, jmt1GL, km, array);
  if(retVal == 0)
    {
    return 0;
    }
  if(popInfo.byteswap)
    {
    // byteswap all arrays
    array.ByteSwap();
    }
  return 1;
}

class InternalMHTWorkArrays
{
// The default constructor and destructor should be sufficient.
public:
  // Work1 is the local partitioning of the grid with no ghost levels
  Matrix2DFloat Work1;
  Matrix2DFloat WorkY1GL;
  void Clear()
  {
    this->Work1.Clear();
    this->WorkY1GL.Clear();
  }
  void Compute(POPInputInformation& popInfo, int imt1GL, int jmt1GL, int km, int* ext3D1GL,
               Matrix2DFloat& tarea1GL, Matrix2DInt& global_kmt1GL, Matrix1DFloat& dz);

  bool IsComputed()
  {
    return (this->Work1.GetData() != NULL);
  }
};

void InternalMHTWorkArrays::Compute(
  POPInputInformation& popInfo, int imt1GL, int jmt1GL, int km, int* ext3D1GL,
  Matrix2DFloat& tarea1GL, Matrix2DInt& global_kmt1GL, Matrix1DFloat& dz)
{
  Matrix3DFloat uet1GL(imt1GL, jmt1GL, km);
  Matrix3DFloat vnt1GL(imt1GL, jmt1GL, km);
  LoadData3D(popInfo.uet_file, popInfo, ext3D1GL, imt1GL, jmt1GL, km, uet1GL);
  LoadData3D(popInfo.vnt_file, popInfo, ext3D1GL, imt1GL, jmt1GL, km, vnt1GL);

  this->Work1.Allocate(imt1GL-1, jmt1GL-1);
  this->WorkY1GL.Allocate(imt1GL, jmt1GL);

  Matrix2DFloat workX1GL(imt1GL, jmt1GL);

  for(int i=0; i<jmt1GL*imt1GL; i++)
    {
    workX1GL(i) = 0.0;
    this->WorkY1GL(i) = 0.0;
    }

  PBCArrays pbcArrays;
  if(popInfo.use_pbc)
    {
    pbcArrays.Compute(popInfo, imt1GL, jmt1GL, km, ext3D1GL, global_kmt1GL, dz);
    }

  // vertical integration of workX1GL and this->WorkY1GL
  for(int k=0; k < km; k++)
    {
    for (int j=0; j<jmt1GL; j++)
      {
      for (int i=0; i<imt1GL; i++)
        {
        if(popInfo.use_pbc)
          {
          workX1GL(i,j) += uet1GL(i,j,k)*tarea1GL(i,j)*pbcArrays.DZT1GL(i,j,k)*4.186e-15;
          this->WorkY1GL(i,j) += vnt1GL(i,j,k)*tarea1GL(i,j)*pbcArrays.DZT1GL(i,j,k)*4.186e-15;
          }
        else
          {
          workX1GL(i,j) += uet1GL(i,j,k)*tarea1GL(i,j)*dz(k)*4.186e-15;
          this->WorkY1GL(i,j) += vnt1GL(i,j,k)*tarea1GL(i,j)*dz(k)*4.186e-15;
          }
        }
      }
    }
  // find divergence of vertically-integrated heat transport
  for(int i=1; i<imt1GL-1; i++)
    {
    for(int j=1; j<jmt1GL-1; j++)
      {
      int i2 = vtkMHTReader::cshift(i, -1, imt1GL);
      int j2 = vtkMHTReader::cshift(j, -1, jmt1GL);
      this->Work1(i-1,j-1) = workX1GL(i,j) - workX1GL(i2,j) + this->WorkY1GL(i,j) - this->WorkY1GL(i,j2);
      }
    }
}

vtkStandardNewMacro(vtkMHTReader);

//-----------------------------------------------------------------------------
vtkMHTReader::vtkMHTReader()
{
  this->MHTWorkArrays = new InternalMHTWorkArrays;
}

//-----------------------------------------------------------------------------
vtkMHTReader::~vtkMHTReader()
{
  if(this->MHTWorkArrays)
    {
    delete this->MHTWorkArrays;
    this->MHTWorkArrays = NULL;
    }
}

//-----------------------------------------------------------------------------
int vtkMHTReader::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **vtkNotUsed(inputVector),
                              vtkInformationVector *outputVector)
{
  vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
  controller->Barrier();
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  // requests the data
  // return 1 for success, return 0 for failure

  vtkInformation *outInfoMHT = outputVector->GetInformationObject(0);
  int extMHT[6] = {0, -1, 0, -1, 0, -1};
  outInfoMHT->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extMHT);
  cerr << vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() << " has mht extent of "
       << extMHT[0] << " " << extMHT[1] << " "
       << extMHT[2] << " " << extMHT[3] << " "
       << extMHT[4] << " " << extMHT[5] << endl;

  vtkRectilinearGrid* output = vtkRectilinearGrid::GetData(outputVector);
  if(this->CalculateMHT(output, extMHT) == 0)
    {
    vtkErrorMacro("CalculateMHT failed.");
    return 0;
    }

  //controller->Barrier();
  timer->StopTimer();
  cerr << controller->GetLocalProcessId() << " has REQUESTDATA time of " << timer->GetElapsedTime() << endl;

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMHTReader::RequestInformation(vtkInformation *vtkNotUsed(request),
                                     vtkInformationVector **vtkNotUsed(inputVector),
                                     vtkInformationVector *outputVector)
{
  // get information about the data
  // return 1 for success, return 0 for failure
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  POPInputInformation popInfo;
  int retval = this->ParseMetaFile(this->GetFileName(), &popInfo);
  if(retval == 0)
    {
    return 0;
    }

  int ny_mht, z;
  this->GetMOCSize(&popInfo, &ny_mht, &z);
  int ext[6] = {0, ny_mht-1, 0, 0, 0, 0};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext, 6);

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMHTReader::CalculateMHT(vtkRectilinearGrid* output, int* extMHT)
{
  // calculates the MOC.
  // results will be stored as fields in grid.
  // ext is the extents of the final moc array that i'm responsible for.
  // return 1 for success, return 0 for failure

  int rank =
    vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  // printf("proc[%i]: output moc extents [%i %i %i %i %i %i]\n", rank,
  //        ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

  POPInputInformation popInfo;
  int retval = this->ParseMetaFile(this->GetFileName(), &popInfo);
  if(retval == 0)
    {
    return 0;
    }

  // figure out local extents of input data fields
  int ext3D[6];   // local data extents without ghost cells
  int ext3D1GL[6];        // local data extents with 1 layers of ghost cells
  int ext3D2GL[6];       // local data extents with 2 layers of ghost cells
  int numberOfProcesses =
    vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses();
  vtkNew<vtkExtentTranslator> translator;
  translator->SetWholeExtent(0, popInfo.global_imt-1,
                             0, popInfo.global_jmt-1,
                             0, popInfo.global_km-1);
  translator->SetNumberOfPieces(numberOfProcesses);
  translator->SetGhostLevel(0);
  int numberOfSplits = sqrt(numberOfProcesses)+1;
  std::vector<int> splitPath(numberOfSplits);
  for(int i=0;i<numberOfSplits;i++)
    {
    splitPath[i] = i % 2;
    }
  translator->SetSplitPath(numberOfSplits, &splitPath[0]);
  translator->SetPiece(rank);
  translator->PieceToExtentByPoints();
  translator->GetExtent(ext3D);
  translator->GetExtent(ext3D1GL);
  translator->GetExtent(ext3D2GL);
  // store the minimum
  this->GlobalIMT0 = ext3D[0];
  this->GlobalJMT0 = ext3D[2];

  // printf("proc[%i]: data extents [%i %i %i %i %i %i]\n", rank, ext3D1GL[0],
  //        ext3D1GL[1], ext3D1GL[2], ext3D1GL[3], ext3D1GL[4], ext3D1GL[5]);

  // the size of the data block i have
  int imt1GL = ext3D1GL[1] - ext3D1GL[0] + 1;
  int jmt1GL = ext3D1GL[3] - ext3D1GL[2] + 1;
  // we only partition in logical x and y so km should be equal to popInfo.global_km
  int km = ext3D[5] - ext3D[4] + 1;

  // add ghost cells to x and y dimension. if on the border, need to wrap
  // around. use value of -1 to signify the wrap around
  imt1GL += 2;
  jmt1GL += 2;
  ext3D1GL[0]--;
  if(ext3D1GL[1] == popInfo.global_imt-1)
    {
    ext3D1GL[1] = -1;
    }
  else
    {
    ext3D1GL[1]++;
    }
  ext3D1GL[2]--;
  if(ext3D1GL[3] == popInfo.global_jmt-1)
    {
    ext3D1GL[3] = -1;
    }
  else
    {
    ext3D1GL[3]++;
    }
  // second layer of ghost cells
  ext3D2GL[0] = ext3D1GL[0] - 1;
  if(ext3D1GL[1] == popInfo.global_imt-1)
    {
    ext3D2GL[1] = -1;
    }
  else if(ext3D1GL[1] == -1)
    {
    ext3D2GL[1] = -2;
    }
  else
    {
    ext3D2GL[1] = ext3D1GL[1] + 1;
    }
  ext3D2GL[2] = ext3D1GL[2] - 1;
  if(ext3D1GL[3] == popInfo.global_jmt-1)
    {
    ext3D2GL[3] = -1;
    }
  else if(ext3D1GL[3] == -1)
    {
    ext3D2GL[3] = -2;
    }
  else
    {
    ext3D2GL[3] = ext3D1GL[3] + 1;
    }

  // size for two layers of ghost cells
  int imt2GL = imt1GL + 2;
  int jmt2GL = jmt1GL + 2;

  // allocate data that will be read from disk
  // x and y directions are allocated with a layer of ghost cells
  Matrix2DDouble uLat1GL(imt1GL, jmt1GL);
  Matrix2DDouble uLong1GL(imt1GL, jmt1GL);
  Matrix2DDouble htn2GL(imt2GL, jmt2GL);
  Matrix2DDouble hte2GL(imt2GL, jmt2GL);
  Matrix1DFloat dz(km);
  Matrix2DInt global_kmt1GL(imt1GL, jmt1GL);
  Matrix2DInt atl_kmt1GL(imt1GL, jmt1GL);
  Matrix2DInt indopac_kmt1GL(imt1GL, jmt1GL);
  Matrix3DFloat u1GL(imt1GL, jmt1GL, km);
  Matrix3DFloat v1GL(imt1GL, jmt1GL, km);

  // load data from disk
  retval = this->LoadData(&popInfo, ext3D1GL, ext3D2GL, imt1GL, jmt1GL, km, uLat1GL, uLong1GL, htn2GL,
                          hte2GL, dz, global_kmt1GL, atl_kmt1GL, indopac_kmt1GL, u1GL, v1GL,
                          imt2GL, jmt2GL);

  if(retval == 0)
    {
    // some error occurred
    vtkErrorMacro("An error occurred in the function LoadData");
    return 0;
    }

  Matrix2DFloat dxu1GL(imt1GL, jmt1GL);
  Matrix2DFloat dyu1GL(imt1GL, jmt1GL);
  Matrix2DFloat tarea1GL(imt1GL, jmt1GL);

  // calculate other horizontal grid fields
  retval = this->grid_stuff(htn2GL, hte2GL, dxu1GL, dyu1GL, tarea1GL, imt1GL, jmt1GL, imt2GL,jmt2GL);

  // no longer needed
  htn2GL.Clear();
  hte2GL.Clear();

  if(retval == 0)
    {
    // some error occurred
    vtkErrorMacro("An error occurred in the function grid_stuff");
    return 0;
    }

  // convert uLat and uLong to degrees
  for(int i=0; i<imt1GL*jmt1GL; i++)
    {
    uLat1GL(i) = uLat1GL(i) * 180.0 / M_PI;
    uLong1GL(i) = uLong1GL(i) * 180.0 / M_PI;
    }

  Matrix2DFloat tLat1GL(imt1GL, jmt1GL);
  this->sw_4pt(tLat1GL, 0.25, uLat1GL, imt1GL, jmt1GL);

  // not needed anymore
  uLat1GL.Clear();
  uLong1GL.Clear();

  // calculate w from u and v
  Matrix3DFloat w1GL(imt1GL, jmt1GL, km);
  this->wcalc(dz, dxu1GL, dyu1GL, tarea1GL, global_kmt1GL, u1GL, v1GL, w1GL, imt1GL, jmt1GL, km);
  u1GL.Clear();  // not needed anymore

  // set up latitude grid to be used and allocate arrays
  int ny_mht, z;
  this->GetMOCSize(&popInfo, &ny_mht, &z);

  Matrix1DFloat lat_mht(ny_mht);
  for(int j=0; j<ny_mht; j++)
    {
    lat_mht(j) = popInfo.ysouth_mht + j * popInfo.dy_mht;
    }

  Matrix2DFloat psi_temp_old(km, ny_mht);

  // mht[][0] -- global mht
  // mht[][1] -- atlantic mht
  // mht[][2] -- indo-pacific mht
  // mht holds the final mht arrays. note that all processes get a full copy
  // of the mht arrays since they are 1D and then we can keep the same
  // partitioning as is used for MOC. The MHT grid output then uses
  // the appropriate parts of the array for the point data.
  Matrix2DFloat mht(ny_mht, 3);
  // mht_tmp is the output array from meridional_heat() method
  Matrix1DFloat mht_temp(ny_mht);

  // clear out temporary MHT work arrays since we don't know if they're
  // valid.
  this->MHTWorkArrays->Clear();
  this->MHTWorkArrays->Compute(popInfo, imt1GL, jmt1GL, km, ext3D1GL, tarea1GL, global_kmt1GL, dz);

  // dzu is for partial bottom cells
  Matrix3DFloat dzu1GL(imt1GL, jmt1GL, km);

// calculate overturning streamfunctions

// first do global
  if(popInfo.do_global)
    {
    int localJIndexMin1GL = -1; // was local_jj
    bool hasGlobalJIndexMin = false;
    float southern_lat = -1000;
    this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, global_kmt1GL, tLat1GL,
                       &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);

    this->meridional_heat(&popInfo, global_kmt1GL, tLat1GL, lat_mht, imt1GL, jmt1GL, ny_mht,
                          localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, mht_temp);
    for(int y=0; y<ny_mht; y++)
      {
      mht(y,0) = mht_temp(y);
      }
    }

// next do atlantic
  if(popInfo.do_atl)
    {
    int localJIndexMin1GL = -1; // was local_jj
    bool hasGlobalJIndexMin = false;
    float southern_lat = -1000;
    this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, atl_kmt1GL, tLat1GL,
                       &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);
    this->meridional_heat(&popInfo, atl_kmt1GL, tLat1GL, lat_mht, imt1GL, jmt1GL, ny_mht,
                          localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, mht_temp);
    for(int y=0; y<ny_mht; y++)
      {
      mht(y,1) = mht_temp(y);
      }
    }

// now do indo-pacific
  if(popInfo.do_indopac)
    {
    int localJIndexMin1GL = -1; // was local_jj
    bool hasGlobalJIndexMin = false;
    float southern_lat = -1000;
    this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, indopac_kmt1GL, tLat1GL,
                       &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);
    this->meridional_heat(&popInfo, indopac_kmt1GL, tLat1GL, lat_mht, imt1GL, jmt1GL, ny_mht, localJIndexMin1GL, localJIndexMin1GL, southern_lat, mht_temp);
    for(int y=0; y<ny_mht; y++)
      {
      mht(y,2) = mht_temp(y);
      }
    }

// calculate depth
  Matrix1DFloat depth(km);
  depth(0) = 0.0;
  for(int k=1; k<km; k++)
    {
    depth(k) = depth(k-1) + dz(k-1);
    }
  for(int k=1; k<km; k++)
    {
    depth(k) = depth(k) * 0.01;  // convert to m
    }

  // actuall build the mht VTK grid
  output->SetExtent(extMHT);
  vtkNew<vtkFloatArray> xCoords;
  xCoords->SetNumberOfTuples(extMHT[1]-extMHT[0]+1);
  for(int i=extMHT[0]; i<=extMHT[1]; i++)
    {
    xCoords->SetValue(i-extMHT[0], lat_mht(i));
    }
  vtkNew<vtkFloatArray> otherCoords;
  otherCoords->InsertNextValue(0.);
  output->SetXCoordinates(xCoords.GetPointer());
  output->SetYCoordinates(otherCoords.GetPointer()); // a single point
  output->SetZCoordinates(otherCoords.GetPointer()); // a single point
  if(popInfo.do_global)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_mht_global");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    for(int i=extMHT[0]; i<=extMHT[1]; i++)
      {
      data->SetValue(i-extMHT[0], mht(i, 0));
      }
    output->GetPointData()->AddArray(data.GetPointer());
    output->GetPointData()->SetScalars(data.GetPointer());
    }
  if(popInfo.do_atl)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_mht_atl");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    for(int i=extMHT[0]; i<=extMHT[1]; i++)
      {
      data->SetValue(i-extMHT[0], mht(i, 1));
      }
    output->GetPointData()->AddArray(data.GetPointer());
    if(!popInfo.do_global)
      {
      output->GetPointData()->SetScalars(data.GetPointer());
      }
    }
  if(popInfo.do_indopac)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_mht_indopac");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    for(int i=extMHT[0]; i<=extMHT[1]; i++)
      {
      data->SetValue(i-extMHT[0], mht(i, 2));
      }
    output->GetPointData()->AddArray(data.GetPointer());
    if(!popInfo.do_global && !popInfo.do_atl)
      {
      output->GetPointData()->SetScalars(data.GetPointer());
      }
    }

  return 1;
}

//-----------------------------------------------------------------------------
void vtkMHTReader::wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu1GL,
                         Matrix2DFloat& dyu1GL, Matrix2DFloat& tarea1GL,
                         Matrix2DInt& kmtb1GL, Matrix3DFloat& u1GL,
                         Matrix3DFloat& v1GL, Matrix3DFloat& w1GL,
                         int imt1GL, int jmt1GL, int km)
{
  // calculate w ,the vertical velocities, since only u and v are read
  // from file

  int i, j;

  Matrix2DFloat wtk1GL(imt1GL,jmt1GL);
  Matrix2DFloat wtkb1GL(imt1GL,jmt1GL);
  Matrix2DFloat work1GL(imt1GL,jmt1GL);
  Matrix2DFloat fue1GL(imt1GL,jmt1GL);
  Matrix2DFloat fuw1GL(imt1GL,jmt1GL);
  Matrix2DFloat fvn1GL(imt1GL,jmt1GL);
  Matrix2DFloat fvs1GL(imt1GL,jmt1GL);

  for(i=0; i<imt1GL*jmt1GL; i++)
    {
    wtkb1GL(i) = 0.0; // set bottom velocity to zero
    wtk1GL(i) = 0.0;  // initialize
    }

  for(int k=km-1; k>=0; k--)  // integrate from bottom up
    {
    // advection fluxes
    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        work1GL(i,j) = 0.5 * u1GL(i,j,k) * dyu1GL(i,j);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        int j2 = cshift(j, -1, jmt1GL);
        fue1GL(i,j) = work1GL(i,j2);
        fue1GL(i,j) += work1GL(i,j);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        int i2 = cshift(i, -1, imt1GL);
        fuw1GL(i,j) = fue1GL(i2,j);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        work1GL(i,j) = 0.5 * v1GL(i,j,k) * dxu1GL(i,j);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        int i2 = cshift(i, -1, imt1GL);
        fvn1GL(i,j) = work1GL(i2,j);
        fvn1GL(i,j) += work1GL(i,j);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        int j2 = cshift(j, -1, jmt1GL);
        fvs1GL(i,j) = fvn1GL(i,j2);
        }
      }

    // calculate vertical velocity at top of kth level
    // (vertical velocity is zero at bottom of T columns)
    for(i=0; i<imt1GL*jmt1GL; i++)
      {
      work1GL(i) = (fvn1GL(i) - fvs1GL(i) + fue1GL(i) - fuw1GL(i)) / tarea1GL(i);
      }

    for(i=0; i<imt1GL*jmt1GL; i++)
      {
      if(k+1 <= kmtb1GL(i))
        {
        wtk1GL(i) = wtkb1GL(i) - dz(k) * work1GL(i);
        }
      }

    for(j=0; j<jmt1GL; j++)
      {
      for(i=0; i<imt1GL; i++)
        {
        w1GL(i,j,k) = wtk1GL(i,j);
        }
      }

    // top value becomes bottom value for next pass
    for(i=0; i<imt1GL*jmt1GL; i++)
      {
      wtkb1GL(i) = wtk1GL(i);
      }

    }  // for(int k=km-1; k<=0; k--)
}

//-----------------------------------------------------------------------------
void vtkMHTReader::wcalc_pbc(Matrix3DFloat& dzu, Matrix2DFloat& dxu,
                             Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
                             Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w, int imt, int jmt)
{
  // calculate w ,the vertical velocities, since only u and v are read
  // from file.
  // this version is used when using partial bottom cells.

  int i, j, k, im1, jm1;
  double p5 = 0.5;
  double fue, fuw, fvn, fvs;
  double wtkb;

  for(i=0; i<imt*jmt; i++)
    {
    wtkb = 0.0; // set bottom velocity to zero
    }

  for(j=0; j<jmt; j++)
    {
    jm1 = j - 1;
    if(j==0)
      {
      // make the value wrap around, which differs from the fortran code
      //jm1 = jmt-1;
      jm1 = 0;  // fortran code
      }

    for(i=0; i<imt; i++)
      {
      im1 = i - 1;
      if(i==0)
        {
        im1 = imt-1;
        }

      wtkb = 0.0;  // vertical velocity is zero at bottom of cell

      if(kmt(i,j) != 0)
        {
        for(k=kmt(i,j)-1; k>=0; k--)
          {
          // advection fluxes
          fue = p5 * (u(i,j  ,k) * dyu(i,j  ) * dzu(i,j  ,k) +
                      u(i,jm1,k) * dyu(i,jm1) * dzu(i,jm1,k));

          fuw = p5 * (u(im1,j  ,k) * dyu(im1,j  ) * dzu(im1,j  ,k) +
                      u(im1,jm1,k) * dyu(im1,jm1) * dzu(im1,jm1,k));

          fvn = p5 * (v(i  ,j,k) * dxu(i  ,j) * dzu(i  ,j,k) +
                      v(im1,j,k) * dxu(im1,j) * dzu(im1,j,k));

          fvs = p5 * (v(i  ,jm1,k) * dxu(i  ,jm1) * dzu(i  ,jm1,k) +
                      v(im1,jm1,k) * dxu(im1,jm1) * dzu(im1,jm1,k));

          // vertical velocity from top of box from continuity equation
          w(i,j,k) = wtkb - (fvn - fvs + fue - fuw)/tarea(i,j);

          // top value becomes bottom for next pass
          wtkb = w(i,j,k);
          }
        }
      }
    }


  printf("done calculating w\n");
}

//-----------------------------------------------------------------------------
void vtkMHTReader::meridional_heat(
  POPInputInformation* popInfo, Matrix2DInt& kmtb1GL, Matrix2DFloat& tLat1GL,
  Matrix1DFloat& lat_mht, int imt1GL, int jmt1GL, int ny_mht,
  int localJIndexMin1GL, bool hasGlobalJIndexMin, float southern_lat,
  Matrix1DFloat& mht)
{
  for(int i=0; i<ny_mht; i++)
    {
    mht(i) = 0.0;
    }

  printf("southernmost j with 1 ghost layer = %i\n", localJIndexMin1GL);
  printf("southernmost lat = %f\n", southern_lat);

  // zonal (over longitude range) integration to find heat transport
  // across southernmost grid circle in basin
  float global_mht0 = 0.0;
  int j1GL = localJIndexMin1GL-1; // this is jj-1 from original code
  int j2_1GL = cshift(j1GL, 1, popInfo->global_jmt+2);

  //acbauer -- localJIndexMin1GL is out of bounds for process 1 of a 2 proc run

  if(hasGlobalJIndexMin)
    {
    //vtkWarningMacro("!!!!!!!!!!!!!!!!!!!!!!!!! the out of bounds value is " << localJIndexMin1GL);
    for(int i=1; i<imt1GL-1; i++)
      {
      if(kmtb1GL(i,j1GL) == 0 && kmtb1GL(i,j2_1GL) > 0)
        {
        global_mht0 += this->MHTWorkArrays->WorkY1GL(i,j1GL);
        }
      }
    }
  float tmp = global_mht0;
  vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
  controller->AllReduce(&tmp, &global_mht0, 1, vtkCommunicator::SUM_OP);

  printf("mht0: %f for imt_local of %d jmt_local %d jIndexMin1GL %d j2_1gl %d\n",
         global_mht0, popInfo->global_imt,
         popInfo->global_jmt, localJIndexMin1GL, j2_1GL);

  // scan through all points of the layer and find all valid points
  // which are not land. record the latitude of the point as well as its work
  // value
  for(int i=1;i<imt1GL-1;i++)
    {
    for(int j=1;j<jmt1GL-1;j++)
      {
      if(kmtb1GL(i, j) > 0)
        {
        int mhtArrayIndex = this->GetLatitudeIndex(lat_mht, tLat1GL(i, j));
        if(mhtArrayIndex >= 0)
          {
          mht(mhtArrayIndex) += this->MHTWorkArrays->Work1(i-1,j-1);
          }
        }
      }
    }

  std::vector<float> tmpArray(mht.GetSize());
  for(unsigned i=0;i<mht.GetSize();i++)
    {
    tmpArray[i] = mht(i);
    }
  controller->AllReduce(&tmpArray[0], mht.GetData(), tmpArray.size(), vtkCommunicator::SUM_OP);

  for(int j=1;j<ny_mht;j++)
    {
    mht(j) += mht(j-1);
    }

  for(int j=0; j<ny_mht; j++)
    {
    if(lat_mht(j) > southern_lat)
      {
      mht(j) += global_mht0;
      }
    }
}

//-----------------------------------------------------------------------------
void vtkMHTReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
