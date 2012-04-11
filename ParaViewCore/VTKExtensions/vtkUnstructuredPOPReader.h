/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkUnstructuredPOPReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkUnstructuredPOPReader - read NetCDF files
// .Author Joshua Wu 09.15.2009
// .SECTION Description
// vtkUnstructuredPOPReader is a source object that reads NetCDF files.
// It should be able to read most any NetCDF file that wants to output a
// rectilinear grid.  The ordering of the variables is changed such that
// the NetCDF x, y, z directions correspond to the vtkRectilinearGrid
// z, y, x directions, respectively.  The striding is done with
// respect to the vtkRectilinearGrid ordering.  Additionally, the
// z coordinates of the vtkRectilinearGrid are negated so that the
// first slice/plane has the highest z-value and the last slice/plane
// has the lowest z-value.

#ifndef __vtkUnstructuredPOPReader_h
#define __vtkUnstructuredPOPReader_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkDataArraySelection;
class vtkCallbackCommand;
class vtkUnstructuredPOPReaderInternal;

#include <iostream>

class VTK_IO_EXPORT vtkUnstructuredPOPReader : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(vtkUnstructuredPOPReader,vtkUnstructuredGridAlgorithm);
  static vtkUnstructuredPOPReader *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //The file to open
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  //Description:
  //Enable subsampling in i,j and k dimensions in the vtkRectilinearGrid
  vtkSetVector3Macro(Stride, int);
  vtkGetVector3Macro(Stride, int);

  //Description:
  //Enable subsampling in i,j and k dimensions in the vtkRectilinearGrid
  vtkSetVector6Macro(VOI, int);
  vtkGetVector6Macro(VOI, int);

  // Description:
  // Variable array selection.
  virtual int GetNumberOfVariableArrays();
  virtual const char *GetVariableArrayName(int idx);
  virtual int GetVariableArrayStatus(const char *name);
  virtual void SetVariableArrayStatus(const char *name, int status);

  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);

  vtkGetMacro(VectorGrid, int);

  vtkSetMacro(VerticalVelocity, bool);
  vtkGetMacro(VerticalVelocity, bool);

protected:
  vtkUnstructuredPOPReader();
  ~vtkUnstructuredPOPReader();

  int RequestData(vtkInformation*,vtkInformationVector**,
                  vtkInformationVector*);

  static void SelectionModifiedCallback(vtkObject *caller, unsigned long eid,
                                        void *clientdata, void *calldata);

  static void EventCallback(vtkObject* caller, unsigned long eid,
                            void* clientdata, void* calldata);

  vtkCallbackCommand* SelectionObserver;

  char *FileName;

  int Stride[3];

  // Description:
  // The NetCDF file descriptor.
  int NCDFFD;

  // Description:
  // The file name of the opened file.
  char* OpenedFileName;

  vtkSetStringMacro(OpenedFileName);

  // Description:
  // The radius of the sphere to be outputted.
  double Radius;

  int VOI[6];
  bool SubsettingXMin;
  bool SubsettingXMax;

  // Description:
  // Specify whether the grid points are at the vector field (U_LAT and U_LON) locations
  // or the scalar field (T_LAT and T_LON) locations or unset.
  // 0 means unset, 2 means vector field, and 1 means scalar field.
  int VectorGrid;

  // Description:
  // If it is a vector grid (i.e. VectorGrid=2) then specify whether or not
  // to compute the vertical velocity.  This can be a costly computation
  // so if the vertical velocity is not needed then keep this off
  // (the default).
  bool VerticalVelocity;

  bool Transform(vtkUnstructuredGrid* grid,
                 size_t* start, size_t* count);

  int ProcessGrid(vtkUnstructuredGrid* grid, int piece, int numberOfPieces, int &numberOfGhostLevels);

  // Returns true for success.  Also fills out wholeExtent.
  bool ReadMetaData(int wholeExtent[6]);

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector);

  void LoadPointData(vtkUnstructuredGrid* grid, int netCDFFD, int varidp,
                     size_t* start, size_t* count, ptrdiff_t* rStride, const char* arrayName);

  void ComputeVerticalVelocity(vtkUnstructuredGrid* grid, int* wholeExtent, int* subExtent);

  // Returns true for success.  If numberOfGhostLevels is 0 and we need to
  // compute gradients, it will increase the number of ghost levels to 1.
  bool GetExtentInformation(int piece, int numberOfPieces, int &numberOfGhostLevels,
                            int* wholeExtent, int *subExtent);


  // returns true for success.
  bool BuildGhostInformation(vtkUnstructuredGrid* grid, int ghostLevels,
                             int* wholeExtent, int* subExtent, int wrapped);

private:
  vtkUnstructuredPOPReader(const vtkUnstructuredPOPReader&);  // Not implemented.
  void operator=(const vtkUnstructuredPOPReader&);  // Not implemented.

  vtkUnstructuredPOPReaderInternal* Internals;
};
#endif
