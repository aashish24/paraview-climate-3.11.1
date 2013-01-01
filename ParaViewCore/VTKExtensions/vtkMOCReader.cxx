/*=========================================================================

  Module:    vtkMOCReader.cxx

  =========================================================================*/

#include "vtkMOCReader.h"

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

namespace
{
  // Description:
  // a point, used for an optimization in computing moc
  struct moc_point_t
  {
    float lat;   // latitude of point
    float work;  // the vertical velocity normalized by some quantities
  };
}

class InternalPBCWorkArrays
{
public:
  Matrix2DFloat DZT1GL; // needed for MHT
  Matrix2DFloat DZU1GL; // needed for MOC

  void Clear()
  {
    this->DZT1GL.Clear();
    this->DZU1GL.Clear();
  }
  void Compute(POPInputInformation& popInfo, int imt1GL, int jmt1GL, int km, Matrix1DFloat& dz);
};

//void InternalPBCWorkArrays::Compute(POPInputInformation& popInfo, int imt1GL, int jmt1GL, int km, Matrix1DFloat& dz)
void InternalPBCWorkArrays::Compute(POPInputInformation& , int, int, int , Matrix1DFloat& )
{
  // // read in pbc file
  // Matrix2DDouble dzbc1GL(imt1GL, jmt1GL);
  // this->DZT1GL->Allocate(imt1GL, jmt1GL, km);
  // this->DZU1GL->Allocate(imt1GL, jmt1GL, km);
  // f = fopen(pbc_file.c_str(), "rb");
  // if(f == NULL)
  //   {
  //   fprintf(stderr,"Error in opening pbc_file: %s\n", pbc_file.c_str());
  //   fprintf(stderr, "Program will now abort...\n");
  //   return 0;
  //   }
  // fread(dzbc1GL.getData(), sizeof(double), imt1GL*jmt1GL, f);
  // fclose(f);
  // printf("done reading in pbc file\n");

  // if(popInfo.byteswap)
  //   {
  //   dzbc1GL.ByteSwap();
  //   }

  // for(int k=0; k<km; k++)
  //   {
  //   for(int j=0; j<jmt1GL; j++)
  //     {
  //     for(int i=0; i<imt1GL; i++)
  //       {
  //       if(kmt(i,j) == k+1)
  //         {
  //         dzt(i,j,k) = dzbc1GL(i,j);
  //         }
  //       else
  //         {
  //         dzt(i,j,k) = dz(k);
  //         }
  //       }
  //     }

  //   // dzu = min of surrounding dzt's
  //   int ip1;
  //   for(int j=0; j<jmt1GL-1; j++)
  //     {
  //     for(int i=0; i<imt1GL; i++)
  //       {
  //       if(i==imt1GL-1)
  //         {
  //         ip1 = 0;
  //       }
  //       else
  //         {
  //         ip1 = i + 1;
  //         }
  //       dzu(i,j,k) = min(min(min(dzt(i,j,k),
  //                                dzt(ip1,j,k)),
  //                            dzt(i,j+1,k)),
  //                        dzt(ip1,j+1,k));
  //       }
  //     }

  //   // assume top row is land
  //   for(int i=0; i<imt1GL; i++)
  //     {
  //     dzu(i,jmt1GL-1,k) = 0;
  //     }
  //   }
}

vtkStandardNewMacro(vtkMOCReader);

//-----------------------------------------------------------------------------
vtkMOCReader::vtkMOCReader()
{
}

//-----------------------------------------------------------------------------
vtkMOCReader::~vtkMOCReader()
{
}

//-----------------------------------------------------------------------------
int vtkMOCReader::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **vtkNotUsed(inputVector),
                              vtkInformationVector *outputVector)
{
  vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
  controller->Barrier();
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  // requests the data
  // return 1 for success, return 0 for failure

  int ext[6] = {0, -1, 0, -1, 0, -1};
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), ext);
  cerr << vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() << " has moc extent of "
       << ext[0] << " " << ext[1] << " "
       << ext[2] << " " << ext[3] << " "
       << ext[4] << " " << ext[5] << endl;

  vtkRectilinearGrid* output = vtkRectilinearGrid::GetData(outputVector);
  if(this->CalculateMOC(output, ext) == 0)
    {
    vtkErrorMacro("CalculateMOC failed.");
    return 0;
    }

  //controller->Barrier();
  timer->StopTimer();
  cerr << controller->GetLocalProcessId() << " has REQUESTDATA time of " << timer->GetElapsedTime() << endl;

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::RequestInformation(vtkInformation *vtkNotUsed(request),
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
  if(popInfo.do_msf)
    {
    int ext[6] = {0, ny_mht-1, 0, z-1, 0, 0};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext, 6);
    }
  else
    {
    int ext[6] = {0, -1, 0, -1, 0, -1};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext, 6);
    }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::CalculateMOC(vtkRectilinearGrid* output, int* ext)
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

  // psi holds the final moc values
  // psi is essentially a 3d array
  // psi[][][0] -- global moc
  // psi[][][1] -- atlantic moc
  // psi[][][2] -- indo-pacific moc
  Matrix3DFloat psi(ny_mht, km, 3);

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

    if(popInfo.do_msf)
      {
      this->moc(&popInfo, v1GL, w1GL, global_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,0) = psi_temp_old(k,y);
          }
        }
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
    if(popInfo.do_msf)
      {
      this->moc(&popInfo, v1GL, w1GL, atl_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,1) = psi_temp_old(k,y);
          }
        }
      } // popInfo.do_msf
    }

// now do indo-pacific
  if(popInfo.do_indopac)
    {
    int localJIndexMin1GL = -1; // was local_jj
    bool hasGlobalJIndexMin = false;
    float southern_lat = -1000;
    this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, indopac_kmt1GL, tLat1GL,
                       &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);
    if(popInfo.do_msf)
      {
      this->moc(&popInfo, v1GL, w1GL, indopac_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,2) = psi_temp_old(k,y);
          }
        }
      } // popInfo.do_msf
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

// these arrays determine the coordinate positions of the grid
  // first do the moc output grid which is a 2D grid
  vtkNew<vtkFloatArray> x_coords;
  x_coords->SetNumberOfTuples(ext[1]-ext[0]+1);
  for(int i=ext[0]; i<=ext[1]; i++)
    {
    x_coords->SetValue(i-ext[0], lat_mht(i));
    }

  vtkNew<vtkFloatArray> y_coords;
  y_coords->SetNumberOfTuples(ext[3]-ext[2]+1);
  for(int i=ext[2]; i<=ext[3]; i++)
    {
    y_coords->SetValue(i-ext[2], depth(i));
    }

  vtkNew<vtkFloatArray> z_coords;
  z_coords->InsertNextValue(0.0);

  output->SetExtent(ext);
  output->SetXCoordinates(x_coords.GetPointer());
  output->SetYCoordinates(y_coords.GetPointer());
  output->SetZCoordinates(z_coords.GetPointer());

// TODO: change the array names
// copy output. remember to only copy the part we want.
  if(popInfo.do_global)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_global");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 0));
        counter++;
        }
      }
    output->GetPointData()->AddArray(data.GetPointer());
    output->GetPointData()->SetScalars(data.GetPointer());
    }
  if(popInfo.do_atl)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_atl");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 1));
        counter++;
        }
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
    data->SetName("reader_moc_indopac");
    data->SetNumberOfTuples(output->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 2));
        counter++;
        }
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
void vtkMOCReader::wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu1GL,
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
void vtkMOCReader::wcalc_pbc(Matrix3DFloat& dzu, Matrix2DFloat& dxu,
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
void vtkMOCReader::moc(POPInputInformation* popInfo, Matrix3DFloat& v1GL, Matrix3DFloat& w1GL,
                       Matrix2DInt& kmtb1GL, Matrix2DFloat& tLat1GL, Matrix2DFloat& dxu1GL,
                       Matrix2DFloat& tarea1GL, Matrix1DFloat& dz, Matrix3DFloat& dzu1GL,
                       Matrix1DFloat& lats, int ny_mht, int localJIndexMin1GL, bool hasGlobalJIndexMin,
                       float southern_lat, int imt1GL, int jmt1GL, Matrix2DFloat& psi)
{
  // calculates the meridional overturning circulation
  Matrix2DFloat work1GL(imt1GL, jmt1GL);
  Matrix1DFloat psi0(popInfo->global_km);
  // initialize psi array
  for(int i=0; i<popInfo->global_km*ny_mht; i++)
    {
    psi(i) = 0.0;
    }
  for(int k=0; k<popInfo->global_km; k++)
    {
    psi0(k) = 0.0;
    }

  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  int rank = controller->GetLocalProcessId();

  // zonal and vertical integration to find streamfunction across southernmost
  // grid circle in basin

  // WORK holds the north-south velocity, modulated by other terms,
  // of one z-slice.
  // psi0(k) holds the sum of all north-south velocities over
  // the latitude jj-1, from depths k to km.
  // everyone computes a local psi0 if necessary, then a reduce is done
  // to get correct psi0. only process 0 will have correct psi0.
  //int j = localJIndexMin1GL-1; not used since localJIndexMin1GL may still be VTK_INT_MAX
  if(hasGlobalJIndexMin)
    {
    for(int i=1; i<imt1GL-1; i++)
      {
      if(popInfo->use_pbc)
        {
        work1GL(i,localJIndexMin1GL-1) = -dzu1GL(i,localJIndexMin1GL,popInfo->global_km-1) * v1GL(i,localJIndexMin1GL-1,popInfo->global_km-1) * dxu1GL(i,localJIndexMin1GL-1);
        }
      else
        {
        work1GL(i,localJIndexMin1GL-1) = -dz(popInfo->global_km-1) * v1GL(i,localJIndexMin1GL-1,popInfo->global_km-1) * dxu1GL(i,localJIndexMin1GL-1);
        }
      }

    psi0(popInfo->global_km-1) = 0.0;
    for(int i=1; i<imt1GL-1; i++)
      {
      int j2 = cshift(localJIndexMin1GL-1, 1, jmt1GL);
      if(kmtb1GL(i,localJIndexMin1GL-1) == 0 && kmtb1GL(i, j2) > 0)
        {
        psi0(popInfo->global_km-1) += work1GL(i,localJIndexMin1GL-1);
        }
      }
    }

  for(int k=popInfo->global_km-1; k>=1; k--)
    {
    if(hasGlobalJIndexMin)
      {
      for(int i=1; i<imt1GL-1; i++)
        {
        if(popInfo->use_pbc)
          {
          work1GL(i,localJIndexMin1GL-1) = -dzu1GL(i,localJIndexMin1GL,k-1) * v1GL(i,localJIndexMin1GL-1,k-1) * dxu1GL(i,localJIndexMin1GL-1);
          }
        else
          {
          work1GL(i,localJIndexMin1GL-1) = -dz(k-1) * v1GL(i,localJIndexMin1GL-1,k-1) * dxu1GL(i,localJIndexMin1GL-1);
          }
        }

      for(int i=1; i<imt1GL-1; i++)
        {
        int j2 = cshift(localJIndexMin1GL-1, 1, jmt1GL);
        if(kmtb1GL(i,localJIndexMin1GL-1) == 0 && kmtb1GL(i, j2) > 0)
          {
          psi0(k-1) += work1GL(i,localJIndexMin1GL-1);
          }
        }
      }
    psi0(k-1) += psi0(k);
    }
  Matrix1DFloat tempArray(popInfo->global_km);
  controller->Reduce(psi0.GetData(), tempArray.GetData(), popInfo->global_km,
                     vtkCommunicator::SUM_OP, 0);

  for(int k=0; k<popInfo->global_km; k++)
    {
    psi0(k) = tempArray(k);
    }
  tempArray.Clear();

  // psi(k,j) holds the moc value of depth k and lat j
  // compute my local moc
  std::vector<moc_point_t> points(imt1GL*jmt1GL);
  for(int k=0; k<popInfo->global_km; k++)
    {
    // first scan through all points of the layer and find all valid points
    // which are not land. record the latitude of the point as well as its
    // work value.
    int npoints = 0;     // number of points at this layer
    int k2 = k+1;  // the actual depth value that should be used
    for(int j=1; j<jmt1GL-1; j++)
      {
      for(int i=1; i<imt1GL-1; i++)
        {
        if(k2 <= kmtb1GL(i,j))
          {
          points[npoints].work = w1GL(i,j,k) * tarea1GL(i,j);
          points[npoints].lat = tLat1GL(i, j);
          npoints++;
          }
        }
      }

    // sort all valid points by latitude
    qsort(&points[0], npoints, sizeof(moc_point_t), vtkMOCReader::compare_latitude);

    // step through latitudes from the bottom up, accumulating all points which
    // are less than or equal to the current latitude. keep track of where in
    // the points array the accumulation has gone to. start each search through
    // the points at the point where we left off before.
    int index = 0;  // the current place of the points array
                    // indicates the first element that has not been read yet.
    for(int j=0; j<ny_mht; j++)
      {
      psi(k,j) = 0.0;
      if(j > 0)
        {
        psi(k,j) = psi(k,j-1);
        }
      float lats_j = lats(j);
      while(index < npoints && points[index].lat < lats_j)
        {
        psi(k,j) += points[index].work;
        index++;
        }
      }
    }

  // the actual moc is the elementwise-sum of all local moc's.
  // use a tree structure to send local mocs to parent node.
  // composite mocs from children and then send to parent.
  // final moc is gathered at the root, process 0.
  Matrix2DFloat psi_temp(popInfo->global_km, ny_mht);
  controller->Reduce(psi.GetData(), psi_temp.GetData(), ny_mht*popInfo->global_km, vtkCommunicator::SUM_OP, 0);

  // at this point process 0 should have the accumulated moc,
  // perform some more processing to get final moc.
  if(rank == 0)
    {
    for(int j=0; j<ny_mht; j++)
      {
      for(int k=0; k<popInfo->global_km; k++)
        {
        psi(k,j) = psi_temp(k,j);
        }
      }
    // add in the baseline of the southernmost latitude
    for(int j=0; j<ny_mht; j++)
      {
      if(lats(j) >= southern_lat)
        {
        for(int k=0; k<popInfo->global_km; k++)
          {
          psi(k,j) += psi0(k);
          }
        }
      }

    // smooth grid-point noise over y
    for(int j=1; j<ny_mht-1; j++)
      {
      for(int k=0; k<popInfo->global_km; k++)
        {
        psi(k,j) = 0.25*(psi(k,j-1) + psi(k,j+1)) + 0.5*psi(k,j);
        }
      }

    // normalize to Sv
    for(int i=0; i<popInfo->global_km*ny_mht; i++)
      {
      psi(i) *= 1e-12;
      }

    float special_value = -1e34;
    // replace any zeroes with the special value
    for(int i=0; i<popInfo->global_km*ny_mht; i++)
      {
      if(psi(i) == 0.0)
        {
        psi(i) = special_value;
        }
      }
    }

  // process 0 broadcasts results to everyone
  controller->Broadcast(psi.GetData(), ny_mht*popInfo->global_km, 0);
}

//-----------------------------------------------------------------------------
void vtkMOCReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
