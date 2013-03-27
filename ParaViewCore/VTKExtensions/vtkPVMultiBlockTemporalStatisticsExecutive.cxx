/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkPVMultiBlockTemporalStatisticsExecutive.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPVMultiBlockTemporalStatisticsExecutive.h"

#include "vtkMultiBlockTemporalStatistics.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkPVMultiBlockTemporalStatisticsExecutive);

//----------------------------------------------------------------------------
vtkPVMultiBlockTemporalStatisticsExecutive::vtkPVMultiBlockTemporalStatisticsExecutive()
{
}

//----------------------------------------------------------------------------
vtkPVMultiBlockTemporalStatisticsExecutive::~vtkPVMultiBlockTemporalStatisticsExecutive()
{
}

//----------------------------------------------------------------------------
void vtkPVMultiBlockTemporalStatisticsExecutive::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkPVMultiBlockTemporalStatisticsExecutive::ProcessRequest(
  vtkInformation *request, vtkInformationVector **inInfo,
  vtkInformationVector *outInfo)
{
  cerr << "=================vtkpvmbtse::processrequest begin\n";
  vtkMultiBlockTemporalStatistics* filter =
    vtkMultiBlockTemporalStatistics::SafeDownCast(this->GetAlgorithm());
  vtkMultiProcessController* controller = filter->GetTimeCompartmentController();
  controller->SetGlobalController(controller);
  int retVal = this->Superclass::ProcessRequest(request, inInfo, outInfo);
  controller->SetGlobalController(filter->GetRealGlobalController());
  cerr << "=================vtkpvmbtse::processrequest end\n";
  return retVal;
}

