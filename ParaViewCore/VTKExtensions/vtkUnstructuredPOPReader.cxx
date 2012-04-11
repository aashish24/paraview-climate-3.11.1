/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkUnstructuredPOPReader.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*
  Ghost cells
  Read in field data from multiple files
  Adjust vector quantities to this coordinate system
*/


#include "vtkUnstructuredPOPReader.h"
#include "vtkCallbackCommand.h"
#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkDataArraySelection.h"
#include "vtkExtentTranslator.h"
#include "vtkFloatArray.h"
#include "vtkGradientFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"

#include "vtk_netcdf.h"
#include "vtk_netcdfcpp.h"

#include "vtkMath.h"
#include <string>
#include <vector>

#include "vtkCleanUnstructuredGrid.h"

#include "vtkIntArray.h"
#include "vtksys/SystemTools.hxx"
#include "vtkMultiProcessController.h"

namespace
{
//-----------------------------------------------------------------------------
  template<class DataType>
  void FixArray(DataType* array, size_t* start,
                size_t* count, ptrdiff_t* rStride, int numberOfDimensions,
                int fileId, int varidp)
  {
    // first get the new values
    DataType newValues[count[0]*count[1]];
    size_t newCount[3] = {count[0], count[1], 1};
    size_t newStart[3] = {start[0], start[1], 0};
    ptrdiff_t newStride[3] = { rStride[0], rStride[1], 1 }; // last newStride shouldn't matter since we're only reading 1 entry in that dimension anyways
    int offset = numberOfDimensions == 2 ? 1 : 0;
    // nc_get_vars_float(fileId, varidp, newStart+offset, newCount+offset,
    //                   newStride+offset, newValues);
    nc_get_vars_float(fileId, varidp, newStart+offset, newCount+offset,
                      newStride, newValues);

    // shift the values and replace
    DataType value;
    if(numberOfDimensions == 3)
      {
      vtkIdType shiftAmount = static_cast<vtkIdType>(count[0]*count[1])-1; // the amount each tuple has to shift
      for(vtkIdType iz=static_cast<vtkIdType>(count[0])-1;iz>=0;iz--)
        {
        for(vtkIdType iy=static_cast<vtkIdType>(count[1])-1;iy>=0;iy--)
          {
          if(shiftAmount != 0)
            {
            for(vtkIdType ix=static_cast<vtkIdType>(count[2])-1;ix>=0;ix--)
              {
              vtkIdType index = ix+iy*count[2]+iz*count[1]*count[2];
              std::copy(array+index, array+index+1, array+index+shiftAmount);
              if(shiftAmount < 0)
                {
                cerr << "i don't even know how to shift!!!!\n";
                }
              }
            }
          // this lat/depth line is shifted so I can now add in the value
          array[count[2]+iy*(count[2]+1)+iz*(count[2]+1)*count[1]] = newValues[iy+count[1]*iz];
          shiftAmount--;
          }
        }
      cerr << "the final shiftAmount is " << shiftAmount << " " << vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() << endl;
      }
    else if(numberOfDimensions == 2)
      {
      vtkIdType shiftAmount = static_cast<vtkIdType>(count[1])-1; // the amount each tuple has to shift
      for(vtkIdType iy=static_cast<vtkIdType>(count[1])-1;iy>=0;iy--)
        {
        if(shiftAmount != 0)
          {
          for(vtkIdType ix=static_cast<vtkIdType>(count[2])-1;ix>=0;ix--)
            {
            vtkIdType index = ix+iy*count[2];
            std::copy(array+index, array+index+1, array+index+shiftAmount);
            if(shiftAmount < 0)
              {
              cerr << "i don't even know how to shift!!!!\n";
              }
            }
          }
        // this lat/depth line is shifted so I can now add in the value
        array[count[2]+iy*(count[2]+1)] = newValues[iy];
        shiftAmount--;
        }
      cerr << "the final shiftAmount is " << shiftAmount << " " << vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() << endl;
      }
  }

//-----------------------------------------------------------------------------
  void ConvertScalarsToVector(vtkPointData* pointData,
                              std::vector<std::string> & scalarArrayNames)
  {
    vtkFloatArray* vectorArray = vtkFloatArray::New();
    vectorArray->SetName("velocity");
    vectorArray->SetNumberOfComponents(3);
    vtkIdType numberOfTuples =
      pointData->GetArray(scalarArrayNames[0].c_str())->GetNumberOfTuples();
    vectorArray->SetNumberOfTuples(numberOfTuples);
    pointData->AddArray(vectorArray);
    vectorArray->Delete();
    std::vector<vtkFloatArray*> scalarArrays;
    for(std::vector<std::string>::iterator it=scalarArrayNames.begin();
        it!=scalarArrayNames.end();it++)
      {
      scalarArrays.push_back(vtkFloatArray::SafeDownCast(pointData->GetArray(it->c_str())));
      }
    float values[3] = {0, 0, 0};
    for(vtkIdType i=0;i<numberOfTuples;i++)
      {
      for(size_t j=0;j<scalarArrays.size();j++)
        {
        scalarArrays[j]->GetTupleValue(i, values+j);
        }
      if(values[0] < -1e+31)
        {
        values[0] = values[1] = values[2] = 0;
        }
      vectorArray->SetTupleValue(i, values);
      }
    // remove the old arrays
    for(std::vector<std::string>::iterator it=scalarArrayNames.begin();
        it!=scalarArrayNames.end();it++)
      {
      pointData->RemoveArray(it->c_str());
      }
  }

//-----------------------------------------------------------------------------
  size_t GetPOPIndexFromGridIndices(
    int numDimensions, size_t* count, size_t* start, size_t* stride,
    int i, int j, int k)
  {
    if(numDimensions == 3)
      {
      size_t iIndex = start[2]+i*stride[2];
      size_t jIndex = start[1]+j*stride[1];
      size_t kIndex = start[0]+k*stride[0];
      return iIndex + jIndex*count[2] + kIndex*count[2]*count[1];
      }
    size_t iIndex = start[1]+i*stride[1];
    size_t jIndex = start[0]+j*stride[0];
    return iIndex + jIndex*count[1];
  }

} // end anonymous namespace

vtkStandardNewMacro(vtkUnstructuredPOPReader);

//============================================================================
#define CALL_NETCDF(call) \
{ \
  int errorcode = call; \
  if (errorcode != NC_NOERR) \
  { \
    vtkErrorMacro(<< "netCDF Error: " << nc_strerror(errorcode)); \
    return 0; \
  } \
}
//============================================================================

class vtkUnstructuredPOPReaderInternal
{
public:
  vtkSmartPointer<vtkDataArraySelection> VariableArraySelection;
  // a mapping from the list of all variables to the list of available
  // point-based variables
  std::vector<int> VariableMap;
  vtkUnstructuredPOPReaderInternal()
    {
      this->VariableArraySelection =
        vtkSmartPointer<vtkDataArraySelection>::New();
    }
};

//----------------------------------------------------------------------------
//set default values
vtkUnstructuredPOPReader::vtkUnstructuredPOPReader()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->FileName = NULL;
  this->Stride[0] = this->Stride[1] = this->Stride[2] = 1;
  this->NCDFFD = 0;
  this->OpenedFileName = NULL;
  this->SelectionObserver = vtkCallbackCommand::New();
  this->SelectionObserver->SetCallback
    (&vtkUnstructuredPOPReader::SelectionModifiedCallback);
  this->SelectionObserver->SetClientData(this);
  this->Internals = new vtkUnstructuredPOPReaderInternal;
  this->Internals->VariableArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this->SelectionObserver);
  // radius of the earth in meters assuming the earth is a perfect
  // sphere.  it isn't.  the range going from 6,353 km to 6,384 km.
  // we can revisit this later.
  this->Radius = 6371000;
  for(int i=0;i<3;i++)
    {
    this->VOI[i*2] = 0;
    this->VOI[i*2+1] = -1;
    }
  this->SubsettingXMin = this->SubsettingXMax = false;
  this->VectorGrid = 0;
  this->VerticalVelocity = false;
}

//----------------------------------------------------------------------------
//delete filename and netcdf file descriptor
vtkUnstructuredPOPReader::~vtkUnstructuredPOPReader()
{
  this->SetFileName(0);
  if(this->OpenedFileName)
    {
    nc_close(this->NCDFFD);
    this->SetOpenedFileName(NULL);
    }
  if(this->SelectionObserver)
    {
    this->SelectionObserver->Delete();
    this->SelectionObserver = NULL;
    }
  if(this->Internals)
    {
    delete this->Internals;
    this->Internals = NULL;
    }
}

//----------------------------------------------------------------------------
void vtkUnstructuredPOPReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;
  os << indent << "OpenedFileName: "
     << (this->OpenedFileName ? this->OpenedFileName : "(NULL)") << endl;
  os << indent << "Stride: {" << this->Stride[0] << ", "
     << this->Stride[1] << ", " << this->Stride[2] << ", "
     << "}" << endl;
  os << indent << "NCDFFD: " << this->NCDFFD << endl;

  this->Internals->VariableArraySelection->PrintSelf(os, indent.GetNextIndent());
}

// NC_MAX_VAR_DIMS comes from the nc library
bool vtkUnstructuredPOPReader::ReadMetaData(int wholeExtent[6])
{
  if(this->FileName == NULL)
    {
    vtkErrorMacro("FileName not set.");
    return false;
    }

  if(this->OpenedFileName == NULL || strcmp(this->OpenedFileName, this->FileName) != 0)
    {
    if(this->OpenedFileName)
      {
      nc_close(this->NCDFFD);
      this->SetOpenedFileName(NULL);
      }
    int retval = nc_open(this->FileName, NC_NOWRITE, &this->NCDFFD);//read file
    if (retval != NC_NOERR)//checks if read file error
      {
      // we don't need to close the file if there was an error opening the file
      vtkErrorMacro(<< "Can't read file " << nc_strerror(retval));
      return false;
      }
    this->SetOpenedFileName(this->FileName);
    }
  // get number of variables from file
  int numberOfVariables;
  nc_inq_nvars(this->NCDFFD, &numberOfVariables);
  int dimidsp[NC_MAX_VAR_DIMS];
  int dataDimension;
  size_t dimensions[4]; //dimension value
  this->Internals->VariableMap.resize(numberOfVariables);
  char variableName[NC_MAX_NAME+1];
  int actualVariableCounter = 0;
  // For every variable in the file
  for(int i=0;i<numberOfVariables;i++)
    {
    this->Internals->VariableMap[i] = -1;
    //get number of dimensions
    CALL_NETCDF(nc_inq_varndims(this->NCDFFD, i, &dataDimension));
    //Variable Dimension ID's containing x,y,z coords for the rectilinear
    //grid spacing
    CALL_NETCDF(nc_inq_vardimid(this->NCDFFD, i, dimidsp));
    if(dataDimension == 3)
      {
      this->Internals->VariableMap[i] = actualVariableCounter++;
      //get variable name
      CALL_NETCDF(nc_inq_varname(this->NCDFFD, i, variableName));
      this->Internals->VariableArraySelection->AddArray(variableName);
      for(int m=0;m<dataDimension;m++)
        {
        CALL_NETCDF(nc_inq_dimlen(this->NCDFFD, dimidsp[m], dimensions+m));
        //acquire variable dimensions
        }
      wholeExtent[0] = wholeExtent[2] = wholeExtent[4] =0; //set extent
      wholeExtent[1] = static_cast<int>((dimensions[2] -1) / this->Stride[0]);
      wholeExtent[3] = static_cast<int>((dimensions[1] -1) / this->Stride[1]);
      wholeExtent[5] = static_cast<int>((dimensions[0] -1) / this->Stride[2]);
      }
    }
  this->SubsettingXMin = this->SubsettingXMax = false;
  if(this->VOI[1] != -1 || this->VOI[3] != -1 || this->VOI[5] != -1)
    {
    for(int i=0;i<3;i++)
      {
      if(wholeExtent[2*i] < this->VOI[2*i] &&
         this->VOI[2*i] <= wholeExtent[2*i+1])
        {
        wholeExtent[2*i] = this->VOI[2*i];
        if(i == 0)
          {
          this->SubsettingXMin = true;
          }
        }
      if(wholeExtent[2*i+1] > this->VOI[2*i+1] &&
         this->VOI[2*i+1] >= wholeExtent[2*i])
        {
        wholeExtent[2*i+1] = this->VOI[2*i+1];
        if(i == 0)
          {
          this->SubsettingXMax = true;
          }
        }
      }
    }

  return true;
}

//-----------------------------------------------------------------------------
int vtkUnstructuredPOPReader::RequestInformation(vtkInformation *request,
                                                 vtkInformationVector **inputVector,
                                                 vtkInformationVector *outputVector)
{
  // Let the superclass do the heavy lifting.
  if (!this->Superclass::RequestInformation(request, inputVector, outputVector))
    {
    return 0;
    }

  // Superclass understands structured data, but we have to handle unstructured
  // "extents" (pieces).
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  return 1;
}

//----------------------------------------------------------------------------
// Setting extents of the rectilinear grid
int vtkUnstructuredPOPReader::RequestData(vtkInformation* request,
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector  )
{
  this->UpdateProgress(0);
  // the default implementation is to do what the old pipeline did find what
  // output is requesting the data, and pass that into ExecuteData
  // which output port did the request come from
  int outputPort = request->Get(vtkDemandDrivenPipeline::FROM_OUTPUT_PORT());
  // if output port is negative then that means this filter is calling the
  // update directly, in that case just assume port 0
  if (outputPort == -1)
    {
    outputPort = 0;
    }

  vtkInformation *outInfo = outputVector->GetInformationObject(outputPort);
  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numberOfPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int numberOfGhostLevels = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  // we use the follwing to get the output grid instead of this->GetOutput() since the
  // vtkInformation object passed in here may be different than the vtkInformation
  // object used in GetOutput() to get the grid pointer.
  vtkUnstructuredGrid* outputGrid = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  int retVal = 0;
  bool merge = true;
  if(merge)
    {
    cerr << "WE ARE MERGING\n";
    vtkNew<vtkUnstructuredGrid> tempGrid;
    retVal = this->ProcessGrid(tempGrid.GetPointer(), piece, numberOfPieces, numberOfGhostLevels);
    vtkNew<vtkCleanUnstructuredGrid> cleanToGrid;
    cleanToGrid->SetInput(tempGrid.GetPointer());
    cleanToGrid->Update();
    outputGrid->ShallowCopy(cleanToGrid->GetOutput());
    }
  else
    {
    cerr << "WE ARE NOT MERGING\n";
    retVal = this->ProcessGrid(outputGrid, piece, numberOfPieces, numberOfGhostLevels);
    }

  if(numberOfGhostLevels > 0)
    {
    cerr << "the number of GHOSTS is " << numberOfGhostLevels << endl;
    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), numberOfGhostLevels);
    }
  else
    {
    cerr << "NO GHOSTS\n";
    }

  return retVal;
}

//----------------------------------------------------------------------------
int vtkUnstructuredPOPReader::ProcessGrid(
  vtkUnstructuredGrid* grid, int piece, int numberOfPieces, int & numberOfGhostLevels)
{
  // Determine if the point data is going to be scalar or vector fields
  std::string fileName = this->FileName;
  size_t position = fileName.find("UVEL");
  this->VectorGrid = 1;
  int netCDFFD = -1; // netcdf file descriptor for VVEL, also used at the end of this method
  if(position != std::string::npos)
    {
    fileName.replace(position, 4, "VVEL");
    cerr << "trying to open " << fileName << endl;
    int retVal = nc_open(fileName.c_str(), NC_NOWRITE, &netCDFFD);
    if (NC_NOERR != retVal )
      {
      // we don't need to close the file if there was an error opening the file
      vtkErrorMacro(<< "Can't read file " << nc_strerror(retVal));
      netCDFFD = -1; // this may not be needed
      }
    else
      {
      int varidp = -1;
      char variableName[] = "VVEL";
      if(NC_NOERR == nc_inq_varid(netCDFFD, variableName, &varidp))
        {
        this->VectorGrid = 2;
        }
      else
        {
        vtkWarningMacro("Could not find VVEL variable");
        nc_close(netCDFFD);
        netCDFFD = -1;
        }
      }
    }
  if(this->VectorGrid == 2 && this->VerticalVelocity == true &&
     numberOfGhostLevels == 0)
    {
    // we need to compute gradients in order to compute the vertical
    // velocity component and to make sure that the we avoid partitioning
    // effects we add a level of ghosts.
    numberOfGhostLevels = 1;
    }
  int subExtent[6], wholeExtent[6];
  this->GetExtentInformation(piece, numberOfPieces, numberOfGhostLevels,
                             wholeExtent, subExtent);

  cerr << "x whole extent " << wholeExtent[0] << " to " << wholeExtent[1] << endl
       << "y whole extent " << wholeExtent[2] << " to " << wholeExtent[3] << endl
       << "z whole extent " << wholeExtent[4] << " to " << wholeExtent[5] << endl;
  //setup extents for netcdf library to read the netcdf data file
  size_t start[]= {subExtent[4]*this->Stride[2], subExtent[2]*this->Stride[1],
                   subExtent[0]*this->Stride[0]};

  // the points do not wrap.  e.g. for the lowest latitude the longitude
  // starts at -179.9 and finishes at 180.0
  // wrapped keeps track of whether or not this proc wraps.  this really
  // keeps track of whether or not we add on the extra cell in the x-dir
  // or not without adding extending the points in the x-direction.
  int wrapped = static_cast<int>(!this->SubsettingXMin && !this->SubsettingXMax &&
                                 subExtent[0] == wholeExtent[0] &&
                                 subExtent[1] == wholeExtent[1]);

  size_t count[]= {subExtent[5]-subExtent[4]+1, subExtent[3]-subExtent[2]+1,
                   subExtent[1]-subExtent[0]+1};

  ptrdiff_t rStride[3] = { (ptrdiff_t)this->Stride[2], (ptrdiff_t)this->Stride[1],
                           (ptrdiff_t)this->Stride[0] };

  //initialize memory (raw data space, x y z axis space) and rectilinear grid
  bool firstPass = true;
  for(size_t i=0;i<this->Internals->VariableMap.size();i++)
    {
    if(this->Internals->VariableMap[i] != -1 &&
       this->Internals->VariableArraySelection->GetArraySetting(
         this->Internals->VariableMap[i]) != 0)
      {
      // varidp is probably i in which case nc_inq_varid isn't needed
      int varidp;
      nc_inq_varid(this->NCDFFD,
                   this->Internals->VariableArraySelection->GetArrayName(
                     this->Internals->VariableMap[i]), &varidp);

      if(firstPass == true)
        {
        int dimidsp[3];
        nc_inq_vardimid(this->NCDFFD, varidp, dimidsp);
        firstPass = false;
        std::vector<float> z(count[0]);
        std::vector<float> y(count[1]);
        std::vector<float> x(count[2]);
        //gets data from x,y,z axis (spacing)
        nc_get_vars_float(this->NCDFFD, dimidsp[0],
                          start, count, rStride, &(z[0]));
        nc_get_vars_float(this->NCDFFD, dimidsp[1],
                          start+1, count+1, rStride+1, &(y[0]));
        nc_get_vars_float(this->NCDFFD, dimidsp[2],
                          start+2, count+2, rStride+2, &(x[0]));

        vtkNew<vtkPoints> points;
        grid->SetPoints(points.GetPointer());
        points->SetDataTypeToFloat();
        points->SetNumberOfPoints(count[0]*count[1]*count[2]);
        float point[3];
        vtkIntArray* indexArray = vtkIntArray::New();
        indexArray->SetNumberOfComponents(2);
        indexArray->SetNumberOfTuples(points->GetNumberOfPoints());
        indexArray->SetName("indices");
        int ind[2];
        for(size_t iz=0;iz<count[0];iz++)
          {
          point[2] = -z[iz];
          for(size_t iy=0;iy<count[1];iy++)
            {
            point[1] = y[iy];
            ind[1] = iy;
            for(size_t ix=0;ix<count[2];ix++)
              {
              point[0] = x[ix];
              ind[0] = ix;
              vtkIdType id = ix + iy*count[2] + iz*count[2]*count[1];
              indexArray->SetTupleValue(id, ind);
              points->SetPoint(id, point);
              }
            }
          }
        grid->GetPointData()->AddArray(indexArray);
        indexArray->Delete();

        // need to create the cells
        vtkIdType pointIds[8] = {0,0,0,0,0,0,0,0};
        // we make sure that value is at least 1 so that we can do both quads and hexes
        size_t count2Plus = std::max(count[2],static_cast<size_t>(1));
        size_t count1Plus = std::max(count[1],static_cast<size_t>(1));
        grid->Allocate(std::max(count[0]-1,static_cast<size_t>(1)) *
                       std::max(count[1]-1,static_cast<size_t>(1)) *
                       std::max(count[2]-1+wrapped,static_cast<size_t>(1)));
        size_t iz=0;
        do // make sure we loop through iz once
          {
          size_t iy=0;
          do // make sure we loop through iy once
            {
            size_t ix=0;
            do // make sure we loop through ix once
              {
              pointIds[0] = ix+iy*count2Plus+(1+iz)*count2Plus*count1Plus;
              pointIds[1] = 1+ix+iy*count2Plus+(1+iz)*count2Plus*count1Plus;
              pointIds[2] = 1+ix+(1+iy)*count2Plus+(1+iz)*count2Plus*count1Plus;
              pointIds[3] = ix+(1+iy)*count2Plus+(1+iz)*count2Plus*count1Plus;
              pointIds[4] = ix+iy*count2Plus+iz*count2Plus*count1Plus;
              pointIds[5] = 1+ix+iy*count2Plus+iz*count2Plus*count1Plus;
              pointIds[6] = 1+ix+(1+iy)*count2Plus+iz*count2Plus*count1Plus;
              pointIds[7] = ix+(1+iy)*count2Plus+iz*count2Plus*count1Plus;

              if(wrapped && ix == count[2]-1)
                {
                pointIds[1] = iy*count2Plus+(1+iz)*count2Plus*count1Plus;
                pointIds[2] = (1+iy)*count2Plus+(1+iz)*count2Plus*count1Plus;
                pointIds[5] = iy*count2Plus+iz*count2Plus*count1Plus;
                pointIds[6] = (1+iy)*count2Plus+iz*count2Plus*count1Plus;
                }

              if(count[0] < 2)
                { // constant depth/logical z
                for(int jj=4;jj<8;jj++)
                  {
                  if(pointIds[jj] >= grid->GetNumberOfPoints() || pointIds[jj] < 0)
                    {
                    vtkErrorMacro("this is the problem with id " << pointIds[jj] << " index " << jj);
                    }
                  }
                grid->InsertNextCell(VTK_QUAD, 4, pointIds+4);
                }
              else if(count[1] < 2)
                { // constant latitude/logical y
                pointIds[6] = pointIds[1];
                pointIds[7] = pointIds[0];
                for(int jj=4;jj<8;jj++)
                  {
                  if(pointIds[jj] >= grid->GetNumberOfPoints() || pointIds[jj] < 0)
                    {
                    vtkErrorMacro("this is the problem with id " << pointIds[jj] << " index " << jj);
                    }
                  }
                grid->InsertNextCell(VTK_QUAD, 4, pointIds+4);
                }
              else if(count[2] < 2)
                { // constant longitude/logical x
                pointIds[6] = pointIds[0];
                pointIds[7] = pointIds[1];
                for(int jj=4;jj<8;jj++)
                  {
                  if(pointIds[jj] >= grid->GetNumberOfPoints() || pointIds[jj] < 0)
                    {
                    vtkErrorMacro("this is the problem with id " << pointIds[jj] << " index " << jj);
                    }
                  }
                grid->InsertNextCell(VTK_QUAD, 4, pointIds+4);
                }
              else
                {
                for(int jj=0;jj<8;jj++)
                  {
                  if(pointIds[jj] >= grid->GetNumberOfPoints() || pointIds[jj] < 0)
                    {
                    vtkErrorMacro("this is the problem with id " << pointIds[jj] << " index " << jj);
                    }
                  }
                grid->InsertNextCell(VTK_HEXAHEDRON, 8, pointIds);
                }
              ix++;
              } while (ix<count[2]-1+wrapped);
            iy++;
            } while (iy<count[1]-1);
          iz++;
          } while (iz<count[0]-1);
        } // done creating the points and cells
      //create vtkFloatArray and get the scalars into it
      this->LoadPointData(grid, this->NCDFFD, varidp, start, count, rStride,
                          this->Internals->VariableArraySelection->GetArrayName(
                            this->Internals->VariableMap[i]));
      }
    this->UpdateProgress((i+1.0)/this->Internals->VariableMap.size());
    }

  if(netCDFFD != -1)
    {
    int varidp = -1;
    char variableName[] = "VVEL";
    if(NC_NOERR == nc_inq_varid(netCDFFD, variableName, &varidp))
      {
      this->LoadPointData(grid, netCDFFD, varidp, start, count, rStride, variableName);
      std::vector<std::string> scalarArrayNames;
      scalarArrayNames.push_back("UVEL");
      scalarArrayNames.push_back("VVEL");
      ConvertScalarsToVector(grid->GetPointData(), scalarArrayNames);
      this->VectorGrid = 2;
      }
    nc_close(netCDFFD);
    }

  // transfrom from logical tripolar coordinates to a sphere.  also transforms
  // any vector quantities
  this->Transform(grid, start, count);

  int retVal = 1;
  if(numberOfGhostLevels > 0)
    {
    retVal = this->BuildGhostInformation(
      grid, numberOfGhostLevels, wholeExtent, subExtent, wrapped);
    if(this->VectorGrid && this->VerticalVelocity)
      {
      this->ComputeVerticalVelocity(grid, wholeExtent, subExtent);
      }
    }
  return 1;
}

//----------------------------------------------------------------------------
//following 5 functions are used for paraview user interface
void vtkUnstructuredPOPReader::SelectionModifiedCallback(vtkObject*, unsigned long,
                                                   void* clientdata, void*)
{
  static_cast<vtkUnstructuredPOPReader*>(clientdata)->Modified();
}

//-----------------------------------------------------------------------------
int vtkUnstructuredPOPReader::GetNumberOfVariableArrays()
{
  return this->Internals->VariableArraySelection->GetNumberOfArrays();
}

//-----------------------------------------------------------------------------
const char* vtkUnstructuredPOPReader::GetVariableArrayName(int index)
{
  if(index < 0 || index >= this->GetNumberOfVariableArrays())
    {
    return NULL;
    }
  return this->Internals->VariableArraySelection->GetArrayName(index);
}

//-----------------------------------------------------------------------------
int vtkUnstructuredPOPReader::GetVariableArrayStatus(const char* name)
{
  return this->Internals->VariableArraySelection->ArrayIsEnabled(name);
}

//-----------------------------------------------------------------------------
void vtkUnstructuredPOPReader::SetVariableArrayStatus(const char* name, int status)
{
  vtkDebugMacro("Set cell array \"" << name << "\" status to: " << status);
  if(this->Internals->VariableArraySelection->ArrayExists(name) == 0)
    {
    vtkErrorMacro(<< name << " is not available in the file.");
    return;
    }
  int enabled = this->Internals->VariableArraySelection->ArrayIsEnabled(name);
  if(status != 0 && enabled == 0)
    {
    this->Internals->VariableArraySelection->EnableArray(name);
    this->Modified();
    }
  else if(status == 0 && enabled != 0)
    {
    this->Internals->VariableArraySelection->DisableArray(name);
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
bool vtkUnstructuredPOPReader::Transform(
  vtkUnstructuredGrid* grid, size_t* start, size_t* count)
{
  if(this->VectorGrid != 1 && this->VectorGrid != 2)
    {
    vtkErrorMacro("Don't know if this should be a scalar or vector field grid.");
    return 0;
    }

  int latlonFileId = 0;
  std::string gridFileName =
    vtksys::SystemTools::GetFilenamePath(this->FileName) + "/GRID.nc";
  //cerr << "opening up GRID file " << gridFileName << endl;
  int retval = nc_open(gridFileName.c_str(), NC_NOWRITE, &latlonFileId);
  if (retval != NC_NOERR)//checks if read file error
    {
    // we don't need to close the file if there was an error opening the file
    vtkErrorMacro(<< "Can't read file " << nc_strerror(retval));
    return 0;
    }

  int varidp;
  if(this->VectorGrid == 2)
    {
    nc_inq_varid(latlonFileId, "U_LON_2D", &varidp);
    }
  else
    {
    nc_inq_varid(latlonFileId, "T_LON_2D", &varidp);
    }
  int dimensionIds[3];
  nc_inq_vardimid(latlonFileId, varidp, dimensionIds);
  size_t zeros[2] = {0, 0};
  size_t dimensions[3];
  nc_inq_dimlen(latlonFileId, dimensionIds[0], dimensions);
  nc_inq_dimlen(latlonFileId, dimensionIds[1], dimensions+1);
  size_t latlonCount[2] = {dimensions[0], dimensions[1]};

  std::vector<float> realLongitude(dimensions[0]*dimensions[1]);
  nc_get_vara_float(latlonFileId, varidp,
                    zeros, latlonCount, &(realLongitude[0]));

  if(this->VectorGrid == 2)
    {
    nc_inq_varid(latlonFileId, "U_LAT_2D", &varidp);
    }
  else
    {
    nc_inq_varid(latlonFileId, "T_LAT_2D", &varidp);
    }

  std::vector<float> realLatitude(dimensions[0]*dimensions[1]);
  nc_get_vara_float(latlonFileId, varidp,
                    zeros, latlonCount, &(realLatitude[0]));

  nc_inq_varid(latlonFileId, "depth_t", &varidp);
  nc_inq_vardimid(latlonFileId, varidp, dimensionIds+2);
  nc_inq_dimlen(latlonFileId, dimensionIds[2], dimensions+2);

  std::vector<float> realHeight(dimensions[2]);
  nc_get_vara_float(latlonFileId, varidp,
                    zeros, &(dimensions[2]), &(realHeight[0]));

  cerr << "DIMENSIONS: (" << dimensions[0] << ", " << dimensions[1] << ", " << dimensions[2] << ")" << endl;


  // the vector arrays that need to be manipulated
  std::vector<vtkFloatArray*> vectorArrays;
  for(int i=0;i<grid->GetPointData()->GetNumberOfArrays();i++)
    {
    if(vtkFloatArray* array = vtkFloatArray::SafeDownCast(
         grid->GetPointData()->GetArray(i)))
      {
      if(array->GetNumberOfComponents() == 3)
        {
        vectorArrays.push_back(array);
        }
      }
    }

  size_t rStride[2] = { (size_t)this->Stride[1],
                        (size_t)this->Stride[0] };

  vtkPoints* points = grid->GetPoints();

  for(size_t j=0;j<count[1];j++) // y index
    {
    for(size_t k=0;k<count[0];k++) // z index
      {
      for(size_t i=0;i<count[2];i++) // x index
        {
        vtkIdType index =i + j*count[2] + k*count[2]*count[1];
        if(index >= points->GetNumberOfPoints())
          {
          vtkErrorMacro("doooh");
          }
        size_t latlonIndex = GetPOPIndexFromGridIndices(2, dimensions, start+1, rStride, i, j, k);
        if(latlonIndex < 0 || latlonIndex >= dimensions[0]*dimensions[1])
          {
          vtkErrorMacro("lat lon index doooooh");
          }
        double point[3];
        points->GetPoint(index, point);
        point[0] = realLongitude[latlonIndex];
        point[1] = realLatitude[latlonIndex];

        // convert to spherical
        double radius = this->Radius - realHeight[k];
        double lonRadians = vtkMath::RadiansFromDegrees(point[0]);
        double latRadians = vtkMath::RadiansFromDegrees(point[1]);
        bool sphere = true;
        if(sphere)
          {
          point[0] = radius * cos(latRadians) * cos(lonRadians);
          point[1] = radius * cos(latRadians) * sin(lonRadians);
          point[2] = radius * sin(latRadians);
          }
        points->SetPoint(index, point);

        for(std::vector<vtkFloatArray*>::iterator vit=vectorArrays.begin();
            vit!=vectorArrays.end();vit++)
          {
          float values[3];
          (*vit)->GetTupleValue(index, values);

          size_t startIndex = latlonIndex;
          size_t endIndex = latlonIndex+1;
          if(start[2]+i*rStride[1] >= dimensions[1]-2)
            {
            startIndex = latlonIndex-1;
            endIndex = latlonIndex;
            }
          float vals[3];

          double startLon = vtkMath::RadiansFromDegrees(realLongitude[startIndex]);
          double startLat = vtkMath::RadiansFromDegrees(realLatitude[startIndex]);

          double startPos[3] = {radius * cos(startLat) * cos(startLon),
                             radius * cos(startLat) * sin(startLon),
                             radius * sin(startLat)};

          double endLon = vtkMath::RadiansFromDegrees(realLongitude[endIndex]);
          double endLat = vtkMath::RadiansFromDegrees(realLatitude[endIndex]);

          double endPos[3] = {radius * cos(endLat) * cos(endLon),
                             radius * cos(endLat) * sin(endLon),
                             radius * sin(endLat)};

          double norm = sqrt((endPos[0] - startPos[0]) * (endPos[0] - startPos[0]) + (endPos[1] - startPos[1]) * (endPos[1] - startPos[1]) + (endPos[2] - startPos[2]) * (endPos[2] - startPos[2]));

          vals[0] = values[0] * (endPos[0] - startPos[0]) / norm;
          vals[1] = values[0] * (endPos[1] - startPos[1]) / norm;
          vals[2] = values[0] * (endPos[2] - startPos[2]) / norm;



          startIndex = latlonIndex;
          endIndex = latlonIndex+dimensions[1];
          if(start[1]+j*rStride[0] >= dimensions[0]-2)
            {
            startIndex = latlonIndex-dimensions[1];
            endIndex = latlonIndex;
            }
          vals[3];

          startLon = vtkMath::RadiansFromDegrees(realLongitude[startIndex]);
          startLat = vtkMath::RadiansFromDegrees(realLatitude[startIndex]);

          startPos[0] = radius * cos(startLat) * cos(startLon);
          startPos[1] = radius * cos(startLat) * sin(startLon);
          startPos[2] = radius * sin(startLat);

          endLon = vtkMath::RadiansFromDegrees(realLongitude[endIndex]);
          endLat = vtkMath::RadiansFromDegrees(realLatitude[endIndex]);

          endPos[0] = radius * cos(endLat) * cos(endLon);
          endPos[1] = radius * cos(endLat) * sin(endLon);
          endPos[2] = radius * sin(endLat);

          norm = sqrt((endPos[0] - startPos[0]) * (endPos[0] - startPos[0]) + (endPos[1] - startPos[1]) * (endPos[1] - startPos[1]) + (endPos[2] - startPos[2]) * (endPos[2] - startPos[2]));

          vals[0] += values[1] * (endPos[0] - startPos[0]) / norm;
          vals[1] += values[1] * (endPos[1] - startPos[1]) / norm;
          vals[2] += values[1] * (endPos[2] - startPos[2]) / norm;

          (*vit)->SetTupleValue(index, vals);
          }

        }
      }
    }

  //cerr << "counts are " << count[2] << " " << count[1] << " " << count[0] << endl;

  nc_close(latlonFileId);
  return true;
}

namespace
{
  // iterate over depth columns and then over the points in a depth column
  class PointIterator
  {
  public:
    PointIterator(int *subExtent)
    {
      this->CurrentColumn = 0;
      this->CurrentColumnIndex = 0;
      for(int i=0;i<6;i++)
        {
        this->SubExtent[i] = subExtent[i];
        }
    }
    vtkIdType BeginColumn()
    {
      this->CurrentColumn = 0;
      this->CurrentColumnIndex = 0;
      return 0;
    }
    vtkIdType EndColumn()
    {
      return (this->SubExtent[1]-this->SubExtent[0])*(this->SubExtent[3]-this->SubExtent[2]);
    }
    vtkIdType NextColumn()
    {
      this->CurrentColumn++;
      return this->CurrentColumn;
    }
    vtkIdType BeginColumnIndex()
    {
      this->CurrentColumnIndex = 0;
      return this->CurrentColumn;
    }
    vtkIdType EndColumnIndex()
    {
      return this->SubExtent[5]-this->SubExtent[4];
    }
    vtkIdType NextColumnIndex()
    {
      this->CurrentColumnIndex++;
      return this->CurrentColumnIndex;
    }
    vtkIdType GetCurrentId()
    {
      return this->CurrentColumn+this->EndColumn()*this->CurrentColumnIndex;
    }

  private:
    int CurrentColumn;
    int CurrentColumnIndex;
    int SubExtent[6];
  };
}

//-----------------------------------------------------------------------------
void vtkUnstructuredPOPReader::LoadPointData(
  vtkUnstructuredGrid* grid, int netCDFFD, int varidp,
  size_t* start, size_t* count, ptrdiff_t* rStride, const char* arrayName)
{
  vtkFloatArray *scalars = vtkFloatArray::New();
  vtkIdType numberOfTuples = grid->GetNumberOfPoints();
  float* data = new float[numberOfTuples];
  nc_get_vars_float(netCDFFD, varidp, start, count, rStride, data);
  scalars->SetArray(data, numberOfTuples, 0, 1);
  //set list of variables to display data on grid
  scalars->SetName(arrayName);
  grid->GetPointData()->AddArray(scalars);
  scalars->Delete();
}

void ComputeDirectionCosines(const double point[3], double directionCosines[3][3])
{
  memcpy(directionCosines[0], point, sizeof(double)*3);
  // first direction is aligned with the radial direction, i.e. the one we want
  vtkMath::Normalize(directionCosines[0]);
  // the other 2 directions are arbitrarily set. i don't really care
  // what those directions are as long as they're perpendicular
  vtkMath::Perpendiculars(directionCosines[0], directionCosines[1], directionCosines[2], 0);
}

//-----------------------------------------------------------------------------
void vtkUnstructuredPOPReader::ComputeVerticalVelocity(
  vtkUnstructuredGrid* grid, int* wholeExtent, int* subExtent)
{
  vtkNew<vtkUnstructuredGrid> tempGrid;
  tempGrid->ShallowCopy(grid);
  vtkNew<vtkGradientFilter> gradientFilter;
  gradientFilter->SetInput(tempGrid.GetPointer());
  gradientFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, "velocity");
  gradientFilter->Update();
  grid->ShallowCopy(gradientFilter->GetOutput());

  std::vector<double> dwdr(grid->GetNumberOfPoints()); // change in velocity in the radial direction
  PointIterator pointIterator(subExtent);
  for(vtkIdType column=pointIterator.BeginColumn();column!=pointIterator.EndColumn();
      column=pointIterator.NextColumn())
    {
    // assuming a sphere, the direction cosines for each column is the same
    // so no need to recompute it if we iterate through this way.
    vtkIdType pointId = pointIterator.GetCurrentId();
    double point[3];
    grid->GetPoint(pointId, point);
    double directionCosines[3][3];
    ComputeDirectionCosines(point, directionCosines);
    for(vtkIdType index=pointIterator.BeginColumnIndex();
        index<pointIterator.EndColumnIndex();index=pointIterator.NextColumnIndex())
      {
      pointId = pointIterator.GetCurrentId();
      if(pointId < 0 || pointId >= grid->GetNumberOfPoints())
        {
        vtkErrorMacro("SCREWED this up!!!");  
        continue;
        }
      double gradient[3][3];
      grid->GetPointData()->GetArray("Gradients")->GetTuple(pointId, &gradient[0][0]);
      double tmp1[3][3], tmp2[3][3];
      vtkMath::Multiply3x3(directionCosines, gradient, tmp1);
      vtkMath::Multiply3x3(directionCosines, tmp1, tmp2);
      // assuming incompressible flow du/dx+dv/dy+dw/dz = 0.
      // remember that I already set the first direction to the one I want
      dwdr[pointId] = -tmp2[1][1] - tmp2[2][2];
      }
    }

  // I now have the change in velocity in the direction I want.  It's time
  // to integrate locally.  If this is in parallel, I'll go back later and
  // add in the values from integrations on other procs

  std::vector<double> w(grid->GetNumberOfPoints()); // vertical velocity
  for(vtkIdType column=pointIterator.BeginColumn();column!=pointIterator.EndColumn();
      column=pointIterator.NextColumn())
    {
    double integratedVelocity = 0; // assume no slip/zero velocity at the "bottom"
    double lastPoint[3] = {0,0,0};
    vtkIdType index=pointIterator.BeginColumnIndex();
    vtkIdType pointId = pointIterator.GetCurrentId();
    if(pointId < 0 || pointId >= grid->GetNumberOfPoints())
      {
      vtkErrorMacro("SCREWED this up!!!");  
      continue;
      }
    grid->GetPoint(pointId, lastPoint);
    w[pointId] = integratedVelocity;
    double lastdwdr = dwdr[pointId];
    for(index=pointIterator.NextColumnIndex();
        index<pointIterator.EndColumnIndex();index=pointIterator.NextColumnIndex())
      {
      if(pointId < 0 || pointId >= grid->GetNumberOfPoints())
        {
        vtkErrorMacro("SCREWED this up!!!");  
        continue;
        }
      double currentPoint[3];
      grid->GetPoint(pointId, currentPoint);
      double currentdwdr = dwdr[pointId];
      double average = (currentdwdr+lastdwdr)*.5;
      double length = sqrt(vtkMath::Distance2BetweenPoints(currentPoint, lastPoint));
      integratedVelocity += average*length;
      w[pointId] = integratedVelocity;
      memcpy(lastPoint, currentPoint, sizeof(double)*3);
      lastdwdr = currentdwdr;
      }
    }
  if(subExtent[4] != wholeExtent[4])
    { // need to receive and adjust values

    }
  if(subExtent[5] != wholeExtent[5])
    { // need to send to the next set of processes so they can adjust values
    }
  // now w[] should have the proper values and we need to add them back in
  // to the point data array
  vtkDataArray* velocityArray = grid->GetPointData()->GetArray("velocity");
  for(vtkIdType column=pointIterator.BeginColumn();column!=pointIterator.EndColumn();
      column=pointIterator.NextColumn())
    {
    vtkIdType pointId = pointIterator.GetCurrentId();
    double direction[3];
    grid->GetPoint(pointId, direction);
    vtkMath::Normalize(direction);
    for(vtkIdType index=pointIterator.BeginColumnIndex();
        index<pointIterator.EndColumnIndex();index=pointIterator.NextColumnIndex())
      {
      double velocity[3];
      vtkIdType pointId = pointIterator.GetCurrentId();
      velocityArray->GetTuple(pointId, velocity);
      for(int i=0;i<3;i++)
        {
        velocity[i] += w[pointId]*direction[i];
        }
      velocityArray->SetTuple(pointId, velocity);
      }
    }
  grid->GetPointData()->RemoveArray("Gradients");
}

//-----------------------------------------------------------------------------
bool vtkUnstructuredPOPReader::GetExtentInformation(
  int piece, int numberOfPieces, int & numberOfGhostLevels,
  int* wholeExtent, int* subExtent)
{
  if(this->ReadMetaData(wholeExtent) == false)
    {
    vtkErrorMacro("Error getting meta data.");
    return false;
    }

  vtkNew<vtkExtentTranslator> extentTranslator;
  extentTranslator->SetPiece(piece);
  extentTranslator->SetNumberOfPieces(numberOfPieces);
  // we only split in the y and z topological directions to make
  // wrapping easier for topological y max.
  std::vector<int> splitPath;
  int numberOfProcesses = vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses();
  int count = 0;
  int extents[2] = {wholeExtent[3]-wholeExtent[2], wholeExtent[5]-wholeExtent[4]};
  do
    {
    if(extents[0] >= extents[1])
      {
      extents[0] /= 2;
      splitPath.push_back(1);
      }
    else
      {
      extents[1] /= 2;
      splitPath.push_back(2);
      }
    count = (count == 0 ? 1 : count*2);
    } while(count < numberOfProcesses);
  extentTranslator->SetSplitPath(static_cast<int>(splitPath.size()), &(splitPath[0]));

  extentTranslator->SetWholeExtent(wholeExtent);

  extentTranslator->SetGhostLevel(numberOfGhostLevels);
  extentTranslator->PieceToExtent();
  extentTranslator->GetExtent(subExtent);

//  vtkWarningMacro(<< getpid() << " " << piece << " " << numberOfPieces << " PPP||subextents are " << subExtent[0] << " " << subExtent[1] << " " << subExtent[2] << " " << subExtent[3] << " " << subExtent[4] << " " << subExtent[5]);
  return true;
}

//-----------------------------------------------------------------------------
bool vtkUnstructuredPOPReader::BuildGhostInformation(
  vtkUnstructuredGrid* grid, int numberOfGhostLevels, int* wholeExtent,
  int* subExtent, int wrapped)
{
  if(numberOfGhostLevels == 0)
    {
    return true;
    }
  vtkUnsignedCharArray* cellGhostLevels = vtkUnsignedCharArray::New();
  cellGhostLevels->SetName("vtkGhostLevels");
  cellGhostLevels->SetNumberOfTuples(grid->GetNumberOfCells());
  grid->GetCellData()->AddArray(cellGhostLevels);
  cellGhostLevels->Delete();

  unsigned char *ia = cellGhostLevels->GetPointer(0);
  for(vtkIdType i=0;i<grid->GetNumberOfCells();i++)
    {
    ia[i] = (unsigned char)0;
    }
  int actualXDimension = subExtent[1]-subExtent[0]+wrapped;
  int k = 0;
  do  // k loop with at least 1 pass
    {
    int kGhostLevel = 0;
    if(k < numberOfGhostLevels && subExtent[4] != wholeExtent[4])
      {
      kGhostLevel = numberOfGhostLevels-k;
      }
    else if(k>=subExtent[5]-subExtent[4]-numberOfGhostLevels &&
            subExtent[5] != wholeExtent[5])
      {
      kGhostLevel = k-subExtent[5]+subExtent[4]+numberOfGhostLevels+1;
      }
    int j=0;
    do // j loop with at least 1 pass
      {
      int jGhostLevel = 0;
      if(j < numberOfGhostLevels && subExtent[2] != wholeExtent[2])
        {
        jGhostLevel = numberOfGhostLevels-j;
        }
      else if(j>=subExtent[3]-subExtent[2]-numberOfGhostLevels &&
              subExtent[3] != wholeExtent[3])
        {
        jGhostLevel = j-subExtent[3]+subExtent[2]+numberOfGhostLevels+1;
        }
      unsigned char ghostLevel =
        static_cast<unsigned char>(std::max(kGhostLevel, jGhostLevel));
      // we never split in the logical x-direction so we'll never
      // need ghost cells in that direction
      if( ! ((k>=numberOfGhostLevels || subExtent[4] == wholeExtent[4])
             && (k<(subExtent[5]-subExtent[4]-numberOfGhostLevels) || subExtent[5] == wholeExtent[5]) &&
             (j>=numberOfGhostLevels || subExtent[2] == wholeExtent[2]) &&
             (j<(subExtent[3]-subExtent[2]-numberOfGhostLevels) || subExtent[3] == wholeExtent[3]) ) )
        {
        int i=0;
        do // i loop with at least 1 pass
          {
          int index = i+j*std::max(actualXDimension,1) +
            k*std::max(actualXDimension,1)*std::max(subExtent[3]-subExtent[2],1);
          if(index < 0 || index >= grid->GetNumberOfCells())
            {
            cerr << index << " CELL ghostlevel ERROR " << grid->GetNumberOfCells() << endl;
            }
          else
            {
            ia[index] = ghostLevel;
            }
          i++;
          }
        while(i<actualXDimension);
        }
      j++;
      }
    while (j<subExtent[3]-subExtent[2]);
    k++;
    }
  while (k<subExtent[5]-subExtent[4]);


  vtkUnsignedCharArray* pointGhostLevels = vtkUnsignedCharArray::New();
  pointGhostLevels->SetName("vtkGhostLevels");
  pointGhostLevels->SetNumberOfTuples(grid->GetNumberOfPoints());
  grid->GetPointData()->AddArray(pointGhostLevels);
  pointGhostLevels->Delete();

  ia = pointGhostLevels->GetPointer(0);
  for(vtkIdType i=0;i<grid->GetNumberOfPoints();i++)
    {
    ia[i] = (unsigned char)0;
    }
  actualXDimension = subExtent[1]-subExtent[0]+1;
  k = 0;
  do // k loop with at least 1 pass
    {
    int kGhostLevel = 0;
    if(k < numberOfGhostLevels && subExtent[4] != wholeExtent[4])
      {
      kGhostLevel = numberOfGhostLevels-k;
      }
    else if(k>subExtent[5]-subExtent[4]-numberOfGhostLevels &&
            subExtent[5] != wholeExtent[5])
      {
      kGhostLevel = k-subExtent[5]+subExtent[4]+numberOfGhostLevels;
      }
    int j = 0;
    do // j loop with at least 1 pass
      {
      int jGhostLevel = 0;
      if(j < numberOfGhostLevels && subExtent[2] != wholeExtent[2])
        {
        jGhostLevel = numberOfGhostLevels-j;
        }
      else if(j>subExtent[3]-subExtent[2]-numberOfGhostLevels &&
              subExtent[3] != wholeExtent[3])
        {
        jGhostLevel = j-subExtent[3]+subExtent[2]+numberOfGhostLevels;
        }
      unsigned char ghostLevel =
        static_cast<unsigned char>(std::max(kGhostLevel, jGhostLevel));
      // we never split in the logical x-direction so we'll never
      // need ghost cells in that direction
      if( ! ((k>=numberOfGhostLevels || subExtent[4] == wholeExtent[4])
             && (k<(subExtent[5]-subExtent[4]-numberOfGhostLevels) || subExtent[5] == wholeExtent[5]) &&
             (j>=numberOfGhostLevels || subExtent[2] == wholeExtent[2]) &&
             (j<(subExtent[3]-subExtent[2]-numberOfGhostLevels) || subExtent[3] == wholeExtent[3]) ) )
        {
        int i=0;
        do // i loop with at least 1 pass
          {
          int index = i+j*actualXDimension +
            k*actualXDimension*(subExtent[3]-subExtent[2]+1);
          if(index < 0 || index >= grid->GetNumberOfPoints())
            {
            cerr << "POINT ghostlevel ERROR\n";
            }
          else
            {
            ia[index] = ghostLevel;
            }
          i++;
          } while (i<actualXDimension);
        }
      j++;
      } while (j<subExtent[3]-subExtent[2]+1);
    k++;
    } while (k<subExtent[5]-subExtent[4]+1);

  return true;
}
