/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiBlockTemporalStatistics.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

/*
 * Copyright 2008 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "vtkMultiBlockTemporalStatistics.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkExtentTranslator.h"
#include "vtkInformation.h"
#include "vtkInformationObjectBaseKey.h"
#include "vtkInformationVector.h"
#ifdef PARAVIEW_USE_MPI
#include "vtkMPIController.h"
#else
#include "vtkMultiProcessController.h"
#endif
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStdString.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <sstream>

vtkInformationKeyMacro(vtkMultiBlockTemporalStatistics,MPI_SUBCOMMUNICATOR,ObjectBase);

namespace
{
  int MonthLengths[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int LeapYearMonthLengths[] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const char* MonthSuffixes[] = {"_JAN", "_FEB", "_MAR", "_APR", "_MAY", "_JUN", "_JUL",
                                 "_AUG", "_SEP", "_OCT", "_NOV", "_DEC"};
  const char* SeasonSuffixes[] = {"_DJF", "_MAM", "_JJA", "_SON"};
}

//=============================================================================
vtkStandardNewMacro(vtkMultiBlockTemporalStatistics);

class vtkMultiBlockTemporalStatisticsInternal
{
public:
#ifdef PARAVIEW_USE_MPI
  vtkSmartPointer<vtkMPIController> SubController;
#endif
  vtkSmartPointer<vtkMultiProcessController> GlobalController;
};
//=============================================================================

namespace
{
  const char * const AVERAGE_SUFFIX = "average";
  const char * const MINIMUM_SUFFIX = "minimum";
  const char * const MAXIMUM_SUFFIX = "maximum";
  const char * const STANDARD_DEVIATION_SUFFIX = "stddev";

  inline vtkStdString vtkMultiBlockTemporalStatisticsMangleName(
    const char *originalName, const char *statisticsSuffix,
    const char *climatologicalSuffix)
  {
    if (!originalName)
      {
      return vtkStdString(statisticsSuffix) + vtkStdString(climatologicalSuffix);
      }
    return vtkStdString(originalName) + "_" + vtkStdString(statisticsSuffix) +
      vtkStdString(climatologicalSuffix);
  }

  template<class T>
  inline void InitializeValues(T *array, vtkIdType arraySize, double initialValue)
  {
    for(vtkIdType i=0;i<arraySize;i++)
      {
      array[i] = initialValue;
      }
  }

  void BuildNewArray(vtkDataArray* array, vtkFieldData* outFd,
                     double initialValue, const char* statisticsSuffix,
                     const char* climatologicalSuffix)
  {
    vtkSmartPointer<vtkDataArray> newArray;
    newArray.TakeReference(vtkDataArray::SafeDownCast(
                             vtkAbstractArray::CreateArray(array->GetDataType())));
    newArray->SetNumberOfComponents(array->GetNumberOfComponents());
    newArray->CopyComponentNames( array );

    newArray->SetNumberOfTuples(array->GetNumberOfTuples()+1);
    newArray->SetName(vtkMultiBlockTemporalStatisticsMangleName(array->GetName(),
                                                                statisticsSuffix, climatologicalSuffix));
    if (outFd->HasArray(newArray->GetName()))
      {
      vtkGenericWarningMacro(<< "Input has two arrays named " << array->GetName()
                             << ".  Output statistics will probably be wrong.");
      return;
      }
    switch (array->GetDataType())
      {
      vtkTemplateMacro(
        InitializeValues(static_cast<VTK_TT*>(newArray->GetVoidPointer(0)),
                         newArray->GetNumberOfComponents()*newArray->GetNumberOfTuples(),
                         initialValue));
      }
    outFd->AddArray(newArray);
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsAccumulateAverage(
    const T *inArray, T *outArray, vtkIdType inArraySize)
  { // this just adds up the values.  at the end we will divide appropriately
    //vtkGenericWarningMacro("average " << inArray[0] << " " << outArray[0]);
    for (vtkIdType i = 0; i < inArraySize; i++)
      {
      outArray[i] += inArray[i];
      }
    // keep track of the number of times we add to this array. the outarray
    // has one more tuple than the in array
    outArray[inArraySize] += 1;
  }

  template<class T>
  inline void vtkMultiBlockTemporalStatisticsAccumulateMinimum(
    const T *inArray, T *outArray, vtkIdType inArraySize)
  {
    for (vtkIdType i = 0; i < inArraySize; i++)
      {
      if (outArray[i] > inArray[i])
        {
        outArray[i] = inArray[i];
        }
      }
  }

  template<class T>
  inline void vtkMultiBlockTemporalStatisticsAccumulateMaximum(
    const T *inArray, T *outArray, vtkIdType inArraySize)
  {
    for (vtkIdType i = 0; i < inArraySize; i++)
      {
      if (outArray[i] < inArray[i])
        {
        outArray[i] = inArray[i];
        }
      }
  }

  // standard deviation one-pass algorithm from
  // http://www.cs.berkeley.edu/~mhoemmen/cs194/Tutorials/variance.pdf
  // this is numerically stable!
  template<class T>
  inline void vtkTemporalStatisticsAccumulateStdDev(
    const T *inArray, T *outArray, const T *previousAverage,
    vtkIdType inArraySize)
  {
    double pass = static_cast<double>(outArray[inArraySize]);
    if(pass != 0)
      {
      for (vtkIdType i = 0; i < inArraySize; i++)
        {
        double temp = inArray[i]-previousAverage[i]/static_cast<double>(pass);
        outArray[i] = outArray[i] + static_cast<T>(
          pass*temp*temp/static_cast<double>(pass+1) );
        }
      }
    // keep track of the number of times we add to this array. the outarray
    // has one more tuple than the in array
    outArray[inArraySize] += 1;
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsFinishAverage(T *outArray, vtkIdType inArraySize)
  {
    double sumSize = static_cast<double>(outArray[inArraySize]);
    if(sumSize > 0.)
      {
      for (vtkIdType i = 0; i < inArraySize; i++)
        {
        outArray[i] /= sumSize;
        }
      }
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsSubtractBFromA(
    T *a, T *b, vtkIdType arraySize)
  {
    for (vtkIdType i = 0; i < arraySize; i++)
      {
      a[i] -= b[i];
      }
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsUpdateGlobalAverage(
    T *globalAverage, T *diffAverage, vtkIdType arraySize,
    int remoteSampleSize, int totalSampleSize)
  {
    double ratio = static_cast<double>(remoteSampleSize) /
      static_cast<double>(totalSampleSize);
    for (vtkIdType i = 0; i < arraySize; i++)
      {
      globalAverage[i] += diffAverage[i]*ratio;
      }
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsUpdateGlobalStdDev(
    T *globalStdDev, T *remoteStdDev, T * diffAverage, vtkIdType arraySize,
    int remoteSampleSize, int currentSampleSize)
  {
    double ratio = static_cast<double>(remoteSampleSize)*static_cast<double>(currentSampleSize) /
      static_cast<double>(remoteSampleSize+currentSampleSize);
    for (vtkIdType i = 0; i < arraySize; i++)
      {
      globalStdDev[i] += remoteStdDev[i]+diffAverage[i]*diffAverage[i]*ratio;
      }
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void vtkMultiBlockTemporalStatisticsFinishGlobalStdDev(
    T *globalStdDev, vtkIdType arraySize,int sampleSize)
  {
    double denominator = 1./sampleSize;
    for (vtkIdType i = 0; i < arraySize; i++)
      {
      globalStdDev[i] = sqrt(globalStdDev[i]*denominator);
      }
  }
}  // end anonymous namespace

//=============================================================================
vtkMultiBlockTemporalStatistics::vtkMultiBlockTemporalStatistics()
{
  this->ComputeAverage = 1;
  this->ComputeMinimum = 1;
  this->ComputeMaximum = 1;
  this->ComputeStandardDeviation = 1;
  this->TimeCompartmentSize = 0;
  this->TimeSpan = 0;
  this->SamplingMethod = 0;
  this->StartDate = 0;
  this->StartYear = 0;
  this->TimeStepLength = 1;
  this->TimeStepType = 0;
  this->Internal = new vtkMultiBlockTemporalStatisticsInternal;

  this->Internal->GlobalController = vtkMultiProcessController::GetGlobalController();
  int numberOfProcesses = this->Internal->GlobalController->GetNumberOfProcesses();

  // below set CurrentTimeIndex as well
  this->SetTimeCompartmentSize(numberOfProcesses);
}

//-----------------------------------------------------------------------------
vtkMultiBlockTemporalStatistics::~vtkMultiBlockTemporalStatistics()
{
  if(this->Internal)
    {
    delete this->Internal;
    this->Internal = NULL;
    }
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "ComputeAverage: " << this->ComputeAverage << endl;
  os << indent << "ComputeMinimum: " << this->ComputeMinimum << endl;
  os << indent << "ComputeMaximum: " << this->ComputeMaximum << endl;
  os << indent << "ComputeStandardDeviation: " << this->ComputeStandardDeviation << endl;
  os << indent << "TimeCompartmentSize: " << this->TimeCompartmentSize<<endl;
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::SetTimeCompartmentSize(int size)
{
  if(size == this->TimeCompartmentSize)
    {
    return;
    }
  int numberOfProcesses = 1;
  if(this->Internal->GlobalController)
    {
    numberOfProcesses = this->Internal->GlobalController->GetNumberOfProcesses();
    }
  if(size < 0 || size > numberOfProcesses || numberOfProcesses % size != 0)
    {
    vtkWarningMacro("Bad TimeCompartmentSize of " << size << " NumberOfProcesses "
                    << numberOfProcesses << " must be an exact multiple. Value not changed.");
    return;
    }
  this->TimeCompartmentSize = size;
#ifdef PARAVIEW_USE_MPI
  if(this->TimeCompartmentSize == numberOfProcesses)
    { // just use the default controller
    this->Internal->SubController = vtkMPIController::SafeDownCast(this->Internal->GlobalController);
    this->Modified();
    this->CurrentTimeIndex = this->GetTimeCompartmentIndex();
    return;
    }
  // set up the time compartment groups
  vtkMPIController* controller = vtkMPIController::SafeDownCast(
    this->Internal->GlobalController);
  if(controller == NULL)
    {
    vtkErrorMacro("I need a global mpi controller and don't have one!");
    return;
    }
  this->Internal->SubController.TakeReference(controller->PartitionController(
                                                controller->GetLocalProcessId()/this->TimeCompartmentSize, 0));

  this->CurrentTimeIndex = this->GetTimeCompartmentIndex();
  this->Modified();
#endif
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::GetTimeCompartmentIndex()
{
  return this->Internal->GlobalController->GetLocalProcessId()/this->TimeCompartmentSize;
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::GetNumberOfTimeCompartments()
{
  int size = this->Internal->GlobalController->GetNumberOfProcesses() /
    this->TimeCompartmentSize;
  if(size == 0)
    {
    size = 1;
    }
  return size;
}

#ifdef PARAVIEW_USE_MPI
//-----------------------------------------------------------------------------
vtkMPIController* vtkMultiBlockTemporalStatistics::GetTimeCompartmentController()
{
  return this->Internal->SubController;
}
#endif

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::GetTimeCompartmentControllerLocalProcessId()
{
#ifdef PARAVIEW_USE_MPI
  return this->Internal->SubController->GetLocalProcessId();
#else
  return 0;
#endif
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation *info)
{
  info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // The output data of this filter has no time assoicated with it.  It is the
  // result of computations that happen over all time.
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

#ifdef PARAVIEW_USE_MPI
  if(vtkMPIController* subController = this->GetTimeCompartmentController())
    {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    inInfo->Set(this->MPI_SUBCOMMUNICATOR(), subController);
    outInfo->Set(this->MPI_SUBCOMMUNICATOR(), subController);
    }
#endif

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector,
  vtkInformationVector *vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  // The RequestData method will tell the pipeline executive to iterate the
  // upstream pipeline to get each time step in order.  The executive in turn
  // will call this method to get the extent request for each iteration (in this
  // case the time step).
  double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  if (inTimes)
    {
    if(this->CurrentTimeIndex >=
       inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
      {
      vtkErrorMacro("Bad time step.  Adding in garbage time step.");
      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS(),
                  &inTimes[0], 1);
      }
    else
      {
      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS(),
                  &inTimes[this->CurrentTimeIndex], 1);
      }
    }
  int piece = 0;
#ifdef PARAVIEW_USE_MPI
  if(vtkMPIController* subController = this->GetTimeCompartmentController())
    {
    piece = subController->GetLocalProcessId();
    }
#endif
  // the parts below are probably not correct for every situation.  i probably
  // should check for a 3d whole extent for structured grids and may need
  // some dustup for unstructured grids.
  if(vtkExtentTranslator* et = vtkExtentTranslator::SafeDownCast(
       inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())))
    {
    int extent[6];
    et->SetNumberOfPieces(this->TimeCompartmentSize);
    et->SetWholeExtent(inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
    et->SetPiece(piece);
    et->SetGhostLevel(0);
    et->PieceToExtent();
    et->GetExtent(extent);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);
    }
  else
    {
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), this->TimeCompartmentSize);
    }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMultiBlockTemporalStatistics::RequestData(
  vtkInformation *request, vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkDataSet *input = vtkDataSet::GetData(inInfo);
  vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::GetData(outInfo);
  int numberOfProcesses = this->Internal->GlobalController->GetNumberOfProcesses();
  int numberOfTimeSteps =
    inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  if (this->CurrentTimeIndex == this->GetTimeCompartmentIndex())
    {
    // First execution initialize output grids
    if(input->IsA("vtkUnstructuredGrid") || input->IsA("vtkPolyData"))
      {
      output->SetNumberOfBlocks(0); // get rid of any existing data sets
      output->SetNumberOfBlocks(numberOfProcesses);
      this->Grid.TakeReference(input->NewInstance());
      if(this->GetTimeCompartmentIndex() == 0)
        {
        output->SetBlock(this->GetTimeCompartmentControllerLocalProcessId(), this->Grid);
        }
      else
        {
        // this may not be the most efficient but should work
        vtkDataSet* temp = this->Grid->NewInstance();
        output->SetBlock(this->Internal->GlobalController->GetLocalProcessId(), temp);
        temp->Delete();
        }
      }
    else // a structured data set input with a multipiece output
      {
      vtkNew<vtkMultiPieceDataSet> multiPiece;
      output->SetNumberOfBlocks(1);
      output->SetBlock(0, multiPiece.GetPointer());
      multiPiece->SetNumberOfPieces(numberOfProcesses);
      this->Grid.TakeReference(input->NewInstance());
      if(this->GetTimeCompartmentIndex() == 0)
        {
        multiPiece->SetPiece(this->GetTimeCompartmentControllerLocalProcessId(), this->Grid);
        }
      else
        { // controller should not be null here
        vtkDataSet* temp = this->Grid->NewInstance();
        multiPiece->SetPiece(this->Internal->GlobalController->GetLocalProcessId(), temp);
        temp->Delete();
        }
      }
    this->InitializeStatistics(numberOfTimeSteps, input, this->Grid);
    }
  if(this->CurrentTimeIndex <
     inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
    // Accumulate new data if this is a valid time step.
    this->AccumulateStatistics(input, this->Grid);
    }

  this->CurrentTimeIndex += this->GetNumberOfTimeCompartments();

  if (  this->CurrentTimeIndex < numberOfTimeSteps)
    {
    // There is still more to do.
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    }
  else
    {
    // We are done.  Finish up including doing the global sum.
    int numberOfLocalTimeSteps =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) /
      this->GetNumberOfTimeCompartments();

    int difference = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) -
      numberOfLocalTimeSteps*this->GetNumberOfTimeCompartments();
    if(difference && this->GetTimeCompartmentIndex() < difference)
      {
      numberOfLocalTimeSteps++;
      }

    this->PostExecute(numberOfTimeSteps, input, this->Grid);
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    // reset back to this processes first time step
    this->CurrentTimeIndex = this->GetTimeCompartmentIndex();
    this->Grid = NULL; // get rid of the reference to this grid
    }

  return 1;
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::AccumulateStatistics(
  vtkDataSet *input, vtkDataSet *output)
{
  this->AccumulateArrays(input->GetFieldData(), output->GetFieldData());
  this->AccumulateArrays(input->GetPointData(), output->GetPointData());
  this->AccumulateArrays(input->GetCellData(), output->GetCellData());
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::AccumulateArrays(
  vtkFieldData *inFd, vtkFieldData *outFd)
{
  int numArrays = inFd->GetNumberOfArrays();
  for (int i = 0; i < numArrays; i++)
    {
    vtkDataArray *inArray = inFd->GetArray(i);
    if (!inArray)
      {
      continue;
      }

    vtkDataArray* outArray = this->GetArray(outFd, inArray, AVERAGE_SUFFIX);
    if (outArray)
      {
      vtkDataArray* stdevOutArray =
        this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX);
      if (stdevOutArray)
        {
        switch (inArray->GetDataType())
          {
          // standard deviation must be called before average since the one-pass
          // algorithm uses the average up to the previous time step
          vtkTemplateMacro(vtkTemporalStatisticsAccumulateStdDev(
                             static_cast<const VTK_TT*>(inArray->GetVoidPointer(0)),
                             static_cast<VTK_TT*>(stdevOutArray->GetVoidPointer(0)),
                             static_cast<const VTK_TT*>(outArray->GetVoidPointer(0)),
                             inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
          }
        // Alert change in data.
        stdevOutArray->DataChanged();
        }

      switch (inArray->GetDataType())
        {
        vtkTemplateMacro(vtkMultiBlockTemporalStatisticsAccumulateAverage(
                           static_cast<const VTK_TT*>(inArray->GetVoidPointer(0)),
                           static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                           inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
        }
      // Alert change in data.
      outArray->DataChanged();
      }

    outArray = this->GetArray(outFd, inArray, MINIMUM_SUFFIX);
    if (outArray)
      {
      switch (inArray->GetDataType())
        {
        vtkTemplateMacro(vtkMultiBlockTemporalStatisticsAccumulateMinimum(
                           static_cast<const VTK_TT*>(inArray->GetVoidPointer(0)),
                           static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                           inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
        }
      // Alert change in data.
      outArray->DataChanged();
      }

    outArray = this->GetArray(outFd, inArray, MAXIMUM_SUFFIX);
    if (outArray)
      {
      switch (inArray->GetDataType())
        {
        vtkTemplateMacro(vtkMultiBlockTemporalStatisticsAccumulateMaximum(
                           static_cast<const VTK_TT*>(inArray->GetVoidPointer(0)),
                           static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                           inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
        }
      // Alert change in data.
      outArray->DataChanged();
      }
    }
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::PostExecute(
  int numberOfTimeSteps, vtkDataSet *input, vtkDataSet *output)
{
  this->FinishArrays(numberOfTimeSteps, input->GetFieldData(), output->GetFieldData());
  this->FinishArrays(numberOfTimeSteps, input->GetPointData(), output->GetPointData());
  this->FinishArrays(numberOfTimeSteps, input->GetCellData(), output->GetCellData());
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::FinishArrays(
  int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd)
{
  int numArrays = inFd->GetNumberOfArrays();
  if(numArrays == 0)
    {
    return;
    }
  // create new communicators so that we can efficiently do a reduce instead
  // of point to point communications. this subcontroller is orthogonal to
  // this->Internal->SubController and only exists if it is different
  // than the global controller.
  int numberOfTimeCompartments = this->GetNumberOfTimeCompartments();

  if(numberOfTimeCompartments == 1)
    {
    this->FinishArraysSerial(numberOfTimeSteps, inFd, outFd);
    }
#ifdef PARAVIEW_USE_MPI
  else
    {
    this->FinishArraysParallel(numberOfTimeSteps, inFd, outFd);
    }
#endif

  // now we go through and set the number of tuples for the output arrays to the proper length
  std::set<std::string> climatologicalSuffixes;
  this->GetAllClimatologicalSuffixes(numberOfTimeSteps, climatologicalSuffixes);
  for(std::set<std::string>::iterator it=climatologicalSuffixes.begin();
      it!=climatologicalSuffixes.end();it++)
    {
    for (int i = 0; i < numArrays; i++)
      {
      vtkDataArray *inArray = inFd->GetArray(i);
      if (!inArray)
        {
        continue;
        }
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, AVERAGE_SUFFIX,
                                                 it->c_str()))
        {
        outArray->SetNumberOfTuples(inArray->GetNumberOfTuples());
        }
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, MINIMUM_SUFFIX,
                                                 it->c_str()))
        {
        outArray->SetNumberOfTuples(inArray->GetNumberOfTuples());
        }
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, MAXIMUM_SUFFIX,
                                                 it->c_str()))
        {
        outArray->SetNumberOfTuples(inArray->GetNumberOfTuples());
        }
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX,
                                                 it->c_str()))
        {
        outArray->SetNumberOfTuples(inArray->GetNumberOfTuples());
        }

      }
    }
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::FinishArraysSerial(
  int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd)
{
  int numArrays = inFd->GetNumberOfArrays();
  std::set<std::string> climatologicalSuffixes;
  this->GetAllClimatologicalSuffixes(numberOfTimeSteps, climatologicalSuffixes);
  for(std::set<std::string>::iterator it=climatologicalSuffixes.begin();
      it!=climatologicalSuffixes.end();it++)
    {
    for (int i = 0; i < numArrays; i++)
      {
      vtkDataArray *inArray = inFd->GetArray(i);
      if (!inArray)
        {
        continue;
        }

      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, AVERAGE_SUFFIX, it->c_str()))
        {
        // first compute this process's average.  every process already
        // knows enough to finish its local std dev computation
        switch (inArray->GetDataType())
          {
          vtkTemplateMacro(vtkMultiBlockTemporalStatisticsFinishAverage(
                             static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                             inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
          }
        if(vtkDataArray* tempOutArray = this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX, it->c_str()))
          {
          double numStepsDouble = tempOutArray->GetComponent(tempOutArray->GetNumberOfTuples()-1, 0);
          switch (inArray->GetDataType())
            {
            vtkTemplateMacro(vtkMultiBlockTemporalStatisticsFinishGlobalStdDev(
                               static_cast<VTK_TT*>(tempOutArray->GetVoidPointer(0)),
                               inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples(),
                               numStepsDouble));
            }
          }
        }
      } // iterating over input field names
    } // iterating over climatological suffixes
}

#ifdef PARAVIEW_USE_MPI
//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::FinishArraysParallel(
  int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd)
{
  int numArrays = inFd->GetNumberOfArrays();
  if(numArrays == 0)
    {
    return;
    }
  // create new communicators so that we can efficiently do a reduce instead
  // of point to point communications. this subcontroller is orthogonal to
  // this->Internal->SubController and only exists if it is different
  // than the global controller.
  int numberOfTimeCompartments = this->GetNumberOfTimeCompartments();
  vtkSmartPointer<vtkMPIController> mySubController;

  if(this->Internal->GlobalController)
    {
    if(numberOfTimeCompartments > 1 &&
       numberOfTimeCompartments == this->Internal->GlobalController->GetNumberOfProcesses())
      {
      mySubController = vtkMPIController::SafeDownCast(this->Internal->GlobalController);
      }
    else if(numberOfTimeCompartments > 1)
      {
      vtkMPIController* controller = vtkMPIController::SafeDownCast(
        this->Internal->GlobalController);
      mySubController.TakeReference(controller->PartitionController(
                                      this->Internal->GlobalController->GetLocalProcessId() %
                                      this->TimeCompartmentSize, 0));
      }
    }

  std::set<std::string> climatologicalSuffixes;
  this->GetAllClimatologicalSuffixes(numberOfTimeSteps, climatologicalSuffixes);
  for(std::set<std::string>::iterator it=climatologicalSuffixes.begin();
      it!=climatologicalSuffixes.end();it++)
    {
    for (int i = 0; i < numArrays; i++)
      {
      vtkDataArray *inArray = inFd->GetArray(i);
      if (!inArray)
        {
        continue;
        }

      // minimum.
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, MINIMUM_SUFFIX, it->c_str()))
        {
        vtkDataArray* tempArray = outArray->NewInstance();
        tempArray->DeepCopy(outArray);
        mySubController->Reduce(tempArray, outArray, vtkCommunicator::MIN_OP, 0);
        tempArray->Delete();
        }
      // maximum.
      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, MAXIMUM_SUFFIX, it->c_str()))
        {
        vtkDataArray* tempArray = outArray->NewInstance();
        tempArray->DeepCopy(outArray);
        mySubController->Reduce(tempArray, outArray, vtkCommunicator::MAX_OP, 0);
        tempArray->Delete();
        }

      if(vtkDataArray* outArray = this->GetArray(outFd, inArray, AVERAGE_SUFFIX, it->c_str()))
        {
        // first compute this process's average.  every process already
        // knows enough to finish its local std dev computation
        switch (inArray->GetDataType())
          {
          vtkTemplateMacro(vtkMultiBlockTemporalStatisticsFinishAverage(
                             static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                             inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
          }
        if(mySubController->GetLocalProcessId() == 0)
          {
          // keep track of the size of our statistical sample set stored in
          // outArray
          vtkDataArray* diffAverageArray = outArray->NewInstance();
          diffAverageArray->SetNumberOfComponents(outArray->GetNumberOfComponents());
          diffAverageArray->SetNumberOfTuples(outArray->GetNumberOfTuples());
          vtkSmartPointer<vtkDataArray> remoteStdDevArray;
          if(this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX), it->c_str())
            {
            remoteStdDevArray.TakeReference(outArray->NewInstance());
            remoteStdDevArray->SetNumberOfComponents(outArray->GetNumberOfComponents());
            remoteStdDevArray->SetNumberOfTuples(outArray->GetNumberOfTuples());
            }

          double currentSampleSize = outArray->GetComponent(outArray->GetNumberOfTuples()-1, 0);
          for(int j=1;j<mySubController->GetNumberOfProcesses();j++)
            {
            double remoteSampleSize = 0;
            mySubController->Receive(&remoteSampleSize, 1, j, j+234);
            if(remoteSampleSize > 0)
              {
              mySubController->Receive(diffAverageArray, j, j+4234);
              // get the difference between the remote average and current global average
              switch (inArray->GetDataType())
                {
                vtkTemplateMacro(vtkMultiBlockTemporalStatisticsSubtractBFromA(
                                   static_cast<VTK_TT*>(diffAverageArray->GetVoidPointer(0)),
                                   static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                                   inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples()));
                }

              // update the current global average
              switch (inArray->GetDataType())
                {
                vtkTemplateMacro(vtkMultiBlockTemporalStatisticsUpdateGlobalAverage(
                                   static_cast<VTK_TT*>(outArray->GetVoidPointer(0)),
                                   static_cast<VTK_TT*>(diffAverageArray->GetVoidPointer(0)),
                                   inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples(),
                                   remoteSampleSize, remoteSampleSize+currentSampleSize));
                }
              // update the current global standard deviation
              if(vtkDataArray* tempArray = this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX, it->c_str()))
                {
                mySubController->Receive(remoteStdDevArray, j, j+8234);
                switch (inArray->GetDataType())
                  {
                  vtkTemplateMacro(vtkMultiBlockTemporalStatisticsUpdateGlobalStdDev(
                                     static_cast<VTK_TT*>(tempArray->GetVoidPointer(0)),
                                     static_cast<VTK_TT*>(remoteStdDevArray->GetVoidPointer(0)),
                                     static_cast<VTK_TT*>(diffAverageArray->GetVoidPointer(0)),
                                     inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples(),
                                     remoteSampleSize, currentSampleSize));
                  }
                }
              }

            currentSampleSize += remoteSampleSize;
            } // iterating over process that need to send to their master proc
          diffAverageArray->Delete();
          // finish the standard deviation computation
          if(vtkDataArray* stdOutArray = this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX, it->c_str()))
            {
            switch (inArray->GetDataType())
              {
              vtkTemplateMacro(vtkMultiBlockTemporalStatisticsFinishGlobalStdDev(
                                 static_cast<VTK_TT*>(stdOutArray->GetVoidPointer(0)),
                                 inArray->GetNumberOfComponents()*inArray->GetNumberOfTuples(),
                                 currentSampleSize));
              }
            }
          }
        else // send array to local process 0
          {
          double temp = outArray->GetComponent(outArray->GetNumberOfTuples()-1, 0);
          mySubController->Send(&temp, 1, 0, mySubController->GetLocalProcessId()+234);
          if(temp > 0)
            {
            mySubController->Send(outArray, 0, mySubController->GetLocalProcessId()+4234);
            if(vtkDataArray* sendOutArray = this->GetArray(outFd, inArray, STANDARD_DEVIATION_SUFFIX, it->c_str()))
              {
              mySubController->Send(sendOutArray, 0, mySubController->GetLocalProcessId()+8234);
              }
            }
          }
        } // if the average array exists
      } // iterating over input field names
    } // iterating over climatological suffixes
}
#endif

//-----------------------------------------------------------------------------
vtkDataArray *vtkMultiBlockTemporalStatistics::GetArray(
  vtkFieldData* fieldData, vtkDataArray *inArray, const char *statisticsSuffix)
{
  vtkStdString outArrayName = vtkMultiBlockTemporalStatisticsMangleName(
    inArray->GetName(), statisticsSuffix, this->GetClimatologicalSuffix().c_str());
  vtkDataArray *outArray = fieldData->GetArray(outArrayName.c_str());
  if (!outArray)
    {
    return NULL;
    }

  // output array has one more tuple than input array
  if (   (inArray->GetNumberOfComponents() != outArray->GetNumberOfComponents())
         || (inArray->GetNumberOfTuples()+1 != outArray->GetNumberOfTuples()) )
    {
    vtkWarningMacro(<< "Size of array " << outArray->GetName()
                    << " has changed.  Does the source change the topology "
                    << " over time? " << inArray->GetNumberOfTuples() << " is the in and the out is "
                    << outArray->GetNumberOfTuples());
    fieldData->RemoveArray(outArray->GetName());
    return NULL;
    }

  return outArray;
}

//-----------------------------------------------------------------------------
vtkDataArray *vtkMultiBlockTemporalStatistics::GetArray(
  vtkFieldData* fieldData, vtkDataArray *inArray, const char *statisticsSuffix,
  const char* climatologicalSuffix)
{
  vtkStdString outArrayName = vtkMultiBlockTemporalStatisticsMangleName(
    inArray->GetName(), statisticsSuffix, climatologicalSuffix);
  vtkDataArray *outArray = fieldData->GetArray(outArrayName.c_str());
  if (!outArray)
    {
    return NULL;
    }

  // output array has one more tuple than input array
  if (   (inArray->GetNumberOfComponents() != outArray->GetNumberOfComponents())
         || (inArray->GetNumberOfTuples()+1 != outArray->GetNumberOfTuples()) )
    {
    vtkWarningMacro(<< "Size of array " << outArray->GetName()
                    << " has changed.  Does the source change the topology "
                    << " over time?");
    fieldData->RemoveArray(outArray->GetName());
    return NULL;
    }

  return outArray;
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::GetDateFromTimeIndex(
  int timeIndex, int& month, int& year, bool& isLeapYear)
{
  if(this->TimeStepType == 0) // daily time steps
    {
    int numberOfDays = timeIndex*TimeStepLength;
    year = this->StartYear;
    while(1)
      {
      if(year % 4 != 0 || (year % 100 == 0 && year % 400 != 0))
        {
        if(numberOfDays < 365)
          {
          isLeapYear = false;
          break;
          }
        year++;
        numberOfDays -= 365;
        }
      else // a leap year
        {
        if(numberOfDays < 366)
          {
          isLeapYear = true;
          break;
          }
        year++;
        numberOfDays -= 366;
        }
      }

    int* months = (isLeapYear == true ? LeapYearMonthLengths : MonthLengths);
    int sum = months[0];
    month = 0;
    while(numberOfDays >= sum)
      {
      month++;
      sum += months[month];
      }
    }
  else if(this->TimeStepType == 1) // monthly time steps
    {
    int numberOfMonths = timeIndex*TimeStepLength;
    year = this->StartYear;
    while(1)
      {
      if(numberOfMonths >= 12)
        {
        numberOfMonths += 12;
        year++;
        }
      else
        {
        isLeapYear = !(year % 4 != 0 || (year % 100 == 0 && year % 400 != 0));
        month = numberOfMonths;
        break;
        }
      }
    }

  return;
}

//-----------------------------------------------------------------------------
std::string vtkMultiBlockTemporalStatistics::GetClimatologicalSuffix()
{
  return this->GetClimatologicalSuffix(this->CurrentTimeIndex);
}

//-----------------------------------------------------------------------------
std::string vtkMultiBlockTemporalStatistics::GetClimatologicalSuffix(
  int timeIndex)
{
  if(this->TimeSpan == AllTimeSteps)
    {
    return "";
    }
  if(this->SamplingMethod == OddEven)
    {
    if(timeIndex % 2 == 0)
      {
      return "_EVEN";
      }
    return "_ODD";
    }
  int month, year;
  bool isLeapYear;
  this->GetDateFromTimeIndex(timeIndex, month, year, isLeapYear);
  if(this->SamplingMethod == Climatology)
    {
    if(this->TimeSpan == Month)
      {
      return MonthSuffixes[month];
      }
    if(this->TimeSpan == Season)
      {
      return SeasonSuffixes[month/4];
      }
    vtkWarningMacro("TimeSpan must be either Month or Season for Climatology Sampling Method.");
    }
  else if(this->SamplingMethod == Consecutive)
    {
    if(this->TimeSpan == Month)
      {
      std::stringstream suffix;
      suffix << MonthSuffixes[month] << month+year*12;
      return suffix.str();
      }
    else if(this->TimeSpan == Year)
      {
      std::stringstream suffix;
      suffix << "Year" << year;
      return suffix.str();
      }
    else if(this->TimeSpan == Decade)
      {
      std::stringstream suffix;
      suffix << "Decade" << year/10;
      return suffix.str();
      }
    }

  vtkWarningMacro("Not able to handle the current settings");
  return "";
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::GetAllClimatologicalSuffixes(
  int numberOfTimeSteps, std::set<std::string>& climatologicalSuffixes)
{
  climatologicalSuffixes.clear();
  if(numberOfTimeSteps == 0 || this->TimeSpan == AllTimeSteps)
    {
    climatologicalSuffixes.insert("");
    }
  for(int i=0;i<numberOfTimeSteps;i++)
    {
    climatologicalSuffixes.insert(this->GetClimatologicalSuffix(i));
    }
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::InitializeStatistics(
  int numberOfTimeSteps, vtkDataSet *input, vtkDataSet *output)
{
  output->CopyStructure(input);
  this->InitializeArrays(numberOfTimeSteps, input->GetFieldData(),
                         output->GetFieldData());
  this->InitializeArrays(numberOfTimeSteps, input->GetPointData(),
                         output->GetPointData());
  this->InitializeArrays(numberOfTimeSteps, input->GetCellData(),
                         output->GetCellData());
}

//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::InitializeArrays(
  int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd)
{
  // Because we need to do mathematical operations, we require all arrays we
  // process to be numeric data (i.e. a vtkDataArray).  We also handle global
  // ids and pedigree ids special (we just pass them).  Ideally would just let
  // vtkFieldData or vtkDataSetAttributes handle this for us, but no such method
  // that fits our needs here.  Thus, we pass data a bit differently then other
  // filters.  If I miss something important, it should be added here.

  outFd->Initialize();

  vtkDataSetAttributes *inDsa = vtkDataSetAttributes::SafeDownCast(inFd);
  vtkDataSetAttributes *outDsa = vtkDataSetAttributes::SafeDownCast(outFd);
  if (inDsa)
    {
    vtkDataArray *globalIds = inDsa->GetGlobalIds();
    vtkAbstractArray *pedigreeIds = inDsa->GetPedigreeIds();
    if (globalIds)
      {
      outDsa->SetGlobalIds(globalIds);
      }
    if (pedigreeIds)
      {
      outDsa->SetPedigreeIds(pedigreeIds);
      }
    }

  int numArrays = inFd->GetNumberOfArrays();
  for (int i = 0; i < numArrays; i++)
    {
    vtkDataArray *array = inFd->GetArray(i);
    if (!array)
      {
      continue;                               // Array not numeric.
      }
    if (outFd->HasArray(array->GetName()))
      {
      continue;    // Must be Ids.
      }

    this->InitializeArray(numberOfTimeSteps, array, outFd);
    }
}


//-----------------------------------------------------------------------------
void vtkMultiBlockTemporalStatistics::InitializeArray(
  int numberOfTimeSteps, vtkDataArray *array, vtkFieldData *outFd)
{
  std::set<std::string> climatologicalSuffixes;
  this->GetAllClimatologicalSuffixes(numberOfTimeSteps, climatologicalSuffixes);

  for(std::set<std::string>::iterator it=climatologicalSuffixes.begin();
      it!=climatologicalSuffixes.end();it++)
    {
    if (this->ComputeAverage)
      {
      BuildNewArray(array, outFd, 0, AVERAGE_SUFFIX, it->c_str());
      }

    if (this->ComputeMinimum)
      {
      BuildNewArray(array, outFd, VTK_FLOAT_MAX, MINIMUM_SUFFIX, it->c_str());
      }

    if (this->ComputeMaximum)
      {
      BuildNewArray(array, outFd, VTK_FLOAT_MIN, MAXIMUM_SUFFIX, it->c_str());
      }
    if (this->ComputeStandardDeviation)
      {
      BuildNewArray(array, outFd, 0, STANDARD_DEVIATION_SUFFIX, it->c_str());
      }
    }
}
