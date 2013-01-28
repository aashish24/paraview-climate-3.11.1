/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiBlockTemporalStatistics.h

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

// .NAME vtkMultiBlockTemporalStatistics - Compute statistics of point or cell data as it changes over time using Time Compartments
//
// .SECTION Description
//
// Given an input that changes over time, vtkMultiBlockTemporalStatistics looks at the
// data for each time step and computes some statistical information of how a
// point or cell variable changes over time.  For example, vtkMultiBlockTemporalStatistics
// can compute the average value of "pressure" over time of each point.
//
// Note that this filter will require the upstream filter to be run on every
// time step that it reports that it can compute.  This may be a time consuming
// operation.
//
// vtkMultiBlockTemporalStatistics ignores the temporal spacing.  Each timestep will be
// weighted the same regardless of how long of an interval it is to the next
// timestep.  Thus, the average statistic may be quite different from an
// integration of the variable if the time spacing varies.
//
// The standard deviation is computed with the following algorithm:
// www.janinebennett.org/index_files/ParallelStatisticsAlgorithms.pdf
// in equation II.4.
//
// Note that the output is a multiblock because not every process will
// have points or cells to output since the input has a different
// partitioning than the output because of the way that we play
// around with the global controller.
//
// .SECTION Thanks
// This class was originally written by Kenneth Moreland (kmorel@sandia.gov)
// from Sandia National Laboratories.
//


#ifndef _vtkMultiBlockTemporalStatistics_h
#define _vtkMultiBlockTemporalStatistics_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkSmartPointer.h" // for vtkDataSet
#include <set>
#include <string>

class vtkDataSet;
class vtkFieldData;
class vtkInformationObjectBaseKey;
class vtkMultiProcessController;
class vtkMultiBlockTemporalStatisticsInternal;

class VTK_EXPORT vtkMultiBlockTemporalStatistics : public vtkMultiBlockDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkMultiBlockTemporalStatistics, vtkMultiBlockDataSetAlgorithm);
  static vtkMultiBlockTemporalStatistics *New();
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // The averages that we do may be over subsets of times. For example, we
  // may want averages for each month or each season instead of over the
  // entire time. Also, we have time and time step information but don't
  // really know what that corresponds to.

  // Description:
  // What time span to compute statistics over. The default is to compute over the entire time
  // span (Total). Month specifies to compute information over each month
  // of the year. For example, if the time range was 5 years, this would
  // compute statistics for all 5 Januaries into a single value.
  // Season specifies to compute information over the 4 seasons of DJF, MAM,
  // JJA, SON. The outputted arrays for Month and Season are appended
  // to indicate the time that it is computed over (e.g. "_DJF" would be
  // appended if Season was specified for the Winter season or _").
  enum TimeSpanType
  {
    AllTimeSteps = 0,
    Month,
    Season,
    Year,
    Decade
  };
  vtkGetMacro(TimeSpan, int);
  vtkSetClampMacro(TimeSpan, int, 0, 4);

  // Description:
  // The number of days in a year. See
  // http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.0/cf-conventions.html#calendar
  // for details.
  /* enum CalendarType */
  /* { */
  /*   Gregorian = 0,         //!< also called standard, 365 day years except every 4th year unless divided by 100 and not by 400 */
  /*   NoLeap,                //!< all 365 day years */
  /*   AllLeap,               //!< all 366 day years */
  /*   ThreeHundredSixty,     //!< all 366 day years with each month being 30 days long */
  /*   Julian                 //!< 365 day years except every 4th year being a leap year with 366 days */
  /* }; */
  /* vtkGetMacro(Calendar, int); */
  /* vtkSetClampMacro(Calendar, int, 0, 4); */

  // Description:
  // Whether to do average over subsets of the time steps.
  // By default his is false.
  enum SamplingMethodType
  {
    Climatology = 0,
    Consecutive,
    OddEven //!< used for testing only -- needs to be removed
  };
  vtkGetMacro(SamplingMethod, int);
  vtkSetClampMacro(SamplingMethod, int, 0, 2);

  vtkGetMacro(TimeStepLength, int);
  vtkSetClampMacro(TimeStepLength, int, 1, VTK_INT_MAX);

  vtkGetMacro(TimeStepType, int);
  vtkSetClampMacro(TimeStepType, int, 0, 1);

  // Description:
  // Turn on/off the computation of the average values over time.  On by
  // default.  The resulting array names have "_average" appended to them.
  vtkGetMacro(ComputeAverage, int);
  vtkSetMacro(ComputeAverage, int);
  vtkBooleanMacro(ComputeAverage, int);

  // Description:
  // Turn on/off the computation of the minimum values over time.  On by
  // default.  The resulting array names have "_minimum" appended to them.
  vtkGetMacro(ComputeMinimum, int);
  vtkSetMacro(ComputeMinimum, int);
  vtkBooleanMacro(ComputeMinimum, int);

  // Description:
  // Turn on/off the computation of the maximum values over time.  On by
  // default.  The resulting array names have "_maximum" appended to them.
  vtkGetMacro(ComputeMaximum, int);
  vtkSetMacro(ComputeMaximum, int);
  vtkBooleanMacro(ComputeMaximum, int);

  // Description:
  // Turn on/off the computation of the standard deviation over time.  On by
  // default.  The resulting array names have "_stddev" appended to them.
  vtkGetMacro(ComputeStandardDeviation, int);
  vtkSetMacro(ComputeStandardDeviation, int);
  vtkBooleanMacro(ComputeStandardDeviation, int);

  // Definition:
  // Set/get the time compartment size.  This is the number of processes
  // working together to process a single time step.
  vtkGetMacro(TimeCompartmentSize, int);
  void SetTimeCompartmentSize(int size);

  // Description:
  // Get the index of the time compartment that this process belongs to.  It will
  // be between 0 and NumberOfTimeCompartments-1.  The first TimeCompartmentSize
  // processes have index 0, the next TimeCompartmentSize processes have index 1, etc.
  int GetTimeCompartmentIndex();

  // Description:
  // Compute dynamically the number of time compartments and return that value.
  int GetNumberOfTimeCompartments();

  // Description:
  // Each time compartment group has their own MPI subcontroller.
  vtkMultiProcessController* GetTimeCompartmentController();

  static vtkInformationObjectBaseKey* MPI_SUBCOMMUNICATOR();

protected:
  vtkMultiBlockTemporalStatistics();
  ~vtkMultiBlockTemporalStatistics();

  int ComputeAverage;
  int ComputeMaximum;
  int ComputeMinimum;
  int ComputeStandardDeviation;

  // Used when iterating the pipeline to keep track of which timestep we are on.
  int CurrentTimeIndex;

  // The number of processes working together on a single data set.
  int TimeCompartmentSize;

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector);
  virtual int RequestUpdateExtent(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector);
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  virtual void AccumulateStatistics(vtkDataSet *input, vtkDataSet *output);
  virtual void AccumulateArrays(vtkFieldData *inFd, vtkFieldData *outFd);

  virtual void PostExecute(int numberOfTimeSteps, vtkDataSet *input, vtkDataSet *output);

  virtual void FinishArrays(int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd);
  virtual void FinishArraysSerial(int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd);
#ifdef PARAVIEW_USE_MPI
  virtual void FinishArraysParallel(int numberOfTimeSteps, vtkFieldData *inFd, vtkFieldData *outFd);
#endif

  void InitializeStatistics(int numberOfTimeSteps, vtkDataSet *input,
                            vtkDataSet *output);
  void InitializeArrays(int numberOfTimeSteps, vtkFieldData *inFd,
                        vtkFieldData *outFd);
  void InitializeArray(int numberOfTimeSteps, vtkDataArray *array,
                       vtkFieldData *outFd);

  // Description:
  // This gets an array with a given statistics suffix. It uses CurrentTimeStep
  // to figure out the proper climatological suffix, if requested.
  virtual vtkDataArray *GetArray(vtkFieldData* outFd, vtkDataArray *inArray,
                                 const char *statisticsSuffix);

  // Description:
  // This gets an array with a given statistics suffix and climatologicalSuffix.
  virtual vtkDataArray *GetArray(vtkFieldData* outFd, vtkDataArray *inArray,
                                 const char *statisticsSuffix, const char* climatologicalSuffix);

  // Description:
  // If a climatology is being computed (i.e. computing statistics
  // over a repeating time frame such as each month), this returns
  // the suffice to add to the output array to specify which time
  // frame is being computed (e.g. FEB for February).
  virtual std::string GetClimatologicalSuffix();
  virtual std::string GetClimatologicalSuffix(int timeIndex);
  virtual void GetAllClimatologicalSuffixes(int numberOfTimeSteps,
                                            std::set<std::string>& climatologicalSuffixes);

  // Description:
  // Given a time index, return the month (0-11), year, and whether
  // or not it is a leap year.
  void GetDateFromTimeIndex(int timeIndex, int& month, int& year, bool& isLeapYear);

private:
  vtkMultiBlockTemporalStatistics(const vtkMultiBlockTemporalStatistics &); // Not implemented.
  void operator=(const vtkMultiBlockTemporalStatistics &);        // Not implemented.

  vtkMultiBlockTemporalStatisticsInternal* Internal;

  // Description:
  // Used to temporarily store the statistical data.
  vtkSmartPointer<vtkDataSet> Grid;

  int TimeSpan;
  //int Calendar;
  int SamplingMethod;

  // Description:
  // The day of the year the first time step corresponds to.
  // Default value is 0 and possible values are between 0 and
  // 365, inclusive.
  int StartDate;

  int StartYear;

  // Description:
  // The number of days, months, etc. between time steps.
  int TimeStepLength;

  // Description:
  // Default is day (0) but can be month (1) also.
  int TimeStepType;
};

#endif //_vtkMultiBlockTemporalStatistics_h
