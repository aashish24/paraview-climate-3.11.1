/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPVMultiBlockTemporalStatisticsExecutive.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPVMultiBlockTemporalStatisticsExecutive - Executive supporting composite datasets.
// .SECTION Description
// vtkPVMultiBlockTemporalStatisticsExecutive is an executive that supports the processing of
// composite dataset. It supports algorithms that are aware of composite
// dataset as well as those that are not. Type checking is performed at run
// time. Algorithms that are not composite dataset-aware have to support
// all dataset types contained in the composite dataset. The pipeline
// execution can be summarized as follows:
//
// * REQUEST_INFORMATION: The producers have to provide information about
// the contents of the composite dataset in this pass.
// Sources that can produce more than one piece (note that a piece is
// different than a block; each piece consistes of 0 or more blocks) should
// set MAXIMUM_NUMBER_OF_PIECES to -1.
//
// * REQUEST_DATA: This is where the algorithms execute. If the
// vtkPVMultiBlockTemporalStatisticsExecutive is assigned to a simple filter,
// it will invoke the  vtkStreamingDemandDrivenPipeline passes in a loop,
// passing a different block each time and will collect the results in a
// composite dataset.
// .SECTION See also
//  vtkCompositeDataSet

#ifndef __vtkPVMultiBlockTemporalStatisticsExecutive_h
#define __vtkPVMultiBlockTemporalStatisticsExecutive_h

#include "vtkPVCompositeDataPipeline.h"

class VTK_EXPORT vtkPVMultiBlockTemporalStatisticsExecutive : public vtkPVCompositeDataPipeline
{
public:
  static vtkPVMultiBlockTemporalStatisticsExecutive* New();
  vtkTypeMacro(vtkPVMultiBlockTemporalStatisticsExecutive,vtkPVCompositeDataPipeline);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual int ProcessRequest(vtkInformation *request, vtkInformationVector **inInfo,
                             vtkInformationVector *outInfo);

protected:
  vtkPVMultiBlockTemporalStatisticsExecutive();
  ~vtkPVMultiBlockTemporalStatisticsExecutive();

private:
  vtkPVMultiBlockTemporalStatisticsExecutive(const vtkPVMultiBlockTemporalStatisticsExecutive&);  // Not implemented.
  void operator=(const vtkPVMultiBlockTemporalStatisticsExecutive&);  // Not implemented.
};

#endif
