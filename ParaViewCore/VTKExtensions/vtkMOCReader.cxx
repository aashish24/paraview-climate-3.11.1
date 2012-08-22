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
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkFloatArray.h"
#include "vtkMultiProcessController.h"
#include "vtkCommunicator.h"


#include "vtkNew.h"
#include "vtkTimerLog.h"

namespace
{
  inline void Swap4(char *n)
  {
    char *n1;
    char c;

    n1 = n + 3;
    c = *n;
    *n = *n1;
    *n1 = c;

    n++;
    n1--;
    c = *n;
    *n = *n1;
    *n1 = c;
  }
  inline void Swap8(char *n)
  {
    char *n1;
    char c;

    n1 = n + 7;
    c = *n;
    *n = *n1;
    *n1 = c;

    n++;
    n1--;
    c = *n;
    *n = *n1;
    *n1 = c;

    n++;
    n1--;
    c = *n;
    *n = *n1;
    *n1 = c;

    n++;
    n1--;
    c = *n;
    *n = *n1;
    *n1 = c;
  }

  // Description:
  // a point, used for an optimization in computing moc
  struct moc_point_t
  {
    float lat;   // latitude of point
    float work;  // the vertical velocity normalized by some quantities
  };
}

// matrix classes used in the MOC code

//-----------------------------------------------------
//  1D float array
//-----------------------------------------------------

class Matrix1DFloat {
public:
  Matrix1DFloat(unsigned int xDim);
  float& operator() (unsigned int i);
  float  operator() (unsigned int i) const;
  float* GetData();
  unsigned GetSize()
  {
    return this->XDim;
  }
  void Clear();
  void ByteSwap();

  ~Matrix1DFloat();                          // Destructor
  //Matrix1DFloat(const Matrix1DFloat& m);    // Copy constructor

private:
  unsigned XDim;
  float* Data;
};

inline
Matrix1DFloat::Matrix1DFloat(unsigned int xDim)
{
  this->XDim = xDim;
  this->Data = new float[xDim];
}

inline
Matrix1DFloat::~Matrix1DFloat()
{
  this->Clear();
}

inline
float* Matrix1DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

inline
void Matrix1DFloat::Clear()
{
  if(this->Data != NULL)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
}

inline
void Matrix1DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

inline
float& Matrix1DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
float Matrix1DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D int array
//-----------------------------------------------------

class Matrix2DInt {
public:
  Matrix2DInt(unsigned int xDim, unsigned int yDim);
  int& operator() (unsigned int i, unsigned j);
  int  operator() (unsigned int i, unsigned int j) const;
  int& operator() (unsigned int i);
  int  operator() (unsigned int i) const;
  int* GetData();
  unsigned GetSize()
  {
    return this->XDim*this->YDim;
  }
  void Clear();
  void ByteSwap();

  ~Matrix2DInt();                          // Destructor
  //Matrix2DInt(const Matrix2DInt& m);    // Copy constructor

private:
  unsigned XDim, YDim;
  int* Data;
};

inline
Matrix2DInt::Matrix2DInt(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new int[xDim * yDim];
}

inline
Matrix2DInt::~Matrix2DInt()
{
  this->Clear();
}

inline
int* Matrix2DInt::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

inline
void Matrix2DInt::Clear()
{
  // deletes the data
  if(this->Data)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
  this->YDim = 0;
}

inline
void Matrix2DInt::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

inline
int& Matrix2DInt::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

inline
int Matrix2DInt::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

inline
int& Matrix2DInt::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
int Matrix2DInt::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D float array
//-----------------------------------------------------

class Matrix2DFloat {
public:
  Matrix2DFloat(unsigned int xDim, unsigned int yDim);
  float& operator() (unsigned int xDim, unsigned yDim);
  float  operator() (unsigned int xDim, unsigned int yDim) const;
  float& operator() (unsigned int i);
  float  operator() (unsigned int i) const;
  float* GetData();
  unsigned GetSize()
  {
    return this->XDim*this->YDim;
  }
  void Clear();
  void ByteSwap();

  ~Matrix2DFloat();                          // Destructor
  //Matrix2DFloat(const Matrix2DFloat& m);    // Copy constructor

private:
  unsigned XDim, YDim;
  float* Data;
};

inline
Matrix2DFloat::Matrix2DFloat(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new float[xDim * yDim];
}

inline
Matrix2DFloat::~Matrix2DFloat()
{
  this->Clear();
}

inline
float* Matrix2DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

inline
void Matrix2DFloat::Clear()
{
  // deletes the data
  if(this->Data)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
  this->YDim = 0;
}

inline
void Matrix2DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

inline
float& Matrix2DFloat::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

inline
float Matrix2DFloat::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

inline
float& Matrix2DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
float Matrix2DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D double array
//-----------------------------------------------------

class Matrix2DDouble {
public:
  Matrix2DDouble(unsigned int xDim, unsigned int yDim);
  double& operator() (unsigned int i, unsigned j);
  double  operator() (unsigned int i, unsigned int j) const;
  double& operator() (unsigned int i);
  double  operator() (unsigned int i) const;
  double* GetData();
  unsigned GetSize()
  {
    return this->XDim*this->YDim;
  }
  void Clear();
  void ByteSwap();

  ~Matrix2DDouble();                          // Destructor
  //Matrix2DDouble(const Matrix2DDouble& m);    // Copy constructor

private:
  unsigned XDim, YDim;
  double* Data;
};

inline
Matrix2DDouble::Matrix2DDouble(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new double[xDim * yDim];
}

inline
Matrix2DDouble::~Matrix2DDouble()
{
  this->Clear();
}

inline
double* Matrix2DDouble::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

inline
void Matrix2DDouble::Clear()
{
  // deletes the data
  if(this->Data)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
  this->YDim = 0;
}

inline
void Matrix2DDouble::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap8((char*)&this->Data[i]);
    }
}

inline
double& Matrix2DDouble::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

inline
double Matrix2DDouble::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

inline
double& Matrix2DDouble::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
double Matrix2DDouble::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  3D float array
//-----------------------------------------------------

class Matrix3DFloat {
public:
  Matrix3DFloat(unsigned int xDim, unsigned int yDim, unsigned int zDim);
  float& operator() (unsigned int i, unsigned j, unsigned int k);
  float  operator() (unsigned int i, unsigned int j, unsigned int k) const;
  float& operator() (unsigned int i);
  float  operator() (unsigned int i) const;
  float* GetData();
  unsigned GetSize()
  {
    return this->XDim*this->YDim*this->ZDim;
  }
  void Clear();
  void ByteSwap();

  ~Matrix3DFloat();                          // Destructor
  //Matrix3DFloat(const Matrix3DFloat& m);    // Copy constructor

private:
  unsigned XDim, YDim, ZDim;
  float* Data;
};

inline
Matrix3DFloat::Matrix3DFloat(unsigned int xDim, unsigned int yDim, unsigned int zDim)
{
  this->XDim = xDim;
  this->YDim = yDim;
  this->ZDim = zDim;

  this->Data = new float[xDim * yDim * zDim];
}

inline
Matrix3DFloat::~Matrix3DFloat()
{
  this->Clear();
}

inline
float* Matrix3DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

inline
void Matrix3DFloat::Clear()
{
  // deletes the data
  if(this->Data)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
  this->YDim = 0;
  this->ZDim = 0;
}

inline
void Matrix3DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim*this->ZDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

inline
float& Matrix3DFloat::operator() (unsigned int i, unsigned int j, unsigned int k)
{
  return this->Data[this->XDim*this->YDim*k + this->XDim*j + i];
}

inline
float Matrix3DFloat::operator() (unsigned int i, unsigned int j, unsigned int k) const
{
  return this->Data[this->XDim*this->YDim*k + this->XDim*j + i];
}

inline
float& Matrix3DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
float Matrix3DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}

// Description:
// holds all configuration information for an MOC calculation
class MOCInfo {

public:
  MOCInfo()
  {
    this->imt = 320;
    this->jmt = 384;
    this->km = 40;
    this->ysouth_mht = -80.0;
    this->ynorth_mht =  80.0;
    this->dy_mht = 1.0;
    this->kmt_global_file = "";
    this->kmt_atl_file = "";
    this->kmt_indopac_file = "";
    this->in_depths = "";
    this->grid_file = "";
    this->u_file = "";
    this->v_file = "";
    this->do_global = true;
    this->do_atl = true;
    this->do_indopac = true;
    this->byteswap = false;
  }

  // store the internal variable values in data to be streamed
  // to other processes.
  void Serialize(vtkMultiProcessStream& data)
  {
    data.Reset();
    data << this->imt << this->jmt << this->km;
    data << this->ysouth_mht << this->ynorth_mht << this->dy_mht;
    data << this->kmt_global_file << this->kmt_atl_file << this->kmt_indopac_file;
    data << this->in_depths << this->grid_file << this->u_file << this->v_file;
    data << static_cast<int>(this->do_global) << static_cast<int>(this->do_atl);
    data << static_cast<int>(this->do_indopac) << static_cast<int>(this->byteswap);
  }

  // set the values from data in the internal variables.
  void Deserialize(vtkMultiProcessStream& data)
  {
    data >> this->imt >> this->jmt >> this->km;
    data >> this->ysouth_mht >> this->ynorth_mht >> this->dy_mht;
    data >> this->kmt_global_file >> this->kmt_atl_file >> this->kmt_indopac_file;
    data >> this->in_depths >> this->grid_file >> this->u_file >> this->v_file;
    int tmp;
    data >> tmp;
    this->do_global = static_cast<bool>(tmp);
    data >> tmp;
    this->do_atl = static_cast<bool>(tmp);
    data >> tmp;
    this->do_indopac = static_cast<bool>(tmp);
    data >> tmp;
    this->byteswap = static_cast<bool>(tmp);
  }

  int imt;                     // size of grid in x dimension
  int jmt;                     // size of grid in y dimension
  int km;                      // size of grid in z dimension
  float ysouth_mht;            // min latitude to compute moc over
  float ynorth_mht;            // max latitude to compute moc over
  float dy_mht;                // step size of each latitude row in moc
  std::string kmt_global_file;             // global kmt
  std::string kmt_atl_file;         // atlantic
  std::string kmt_indopac_file;     // indian and pacific
  std::string in_depths;            // index of max depth
  std::string grid_file;            // grid information
  std::string u_file;               // u-velocities
  std::string v_file;               // v-velocities
  bool do_global;              // compute global moc or not
  bool do_atl;                 // compute atlantic moc or not
  bool do_indopac;             // compute indian-pacific moc or not
  bool byteswap;               // byteswap the binary input files
};

// TODO: change to vtk coding style


vtkStandardNewMacro(vtkMOCReader);

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
vtkMOCReader::vtkMOCReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//-----------------------------------------------------------------------------
vtkMOCReader::~vtkMOCReader()
{
  this->SetFileName(0);
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

  vtkRectilinearGrid* grid = vtkRectilinearGrid::GetData(outputVector);
  if(this->CalculateMOC(grid, ext) == 0)
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

  MOCInfo mocInfo;
  int retval = this->ParseHeader(this->FileName, &mocInfo);
  if(retval == 0)
    {
    return 0;
    }

  int y, z;
  this->GetMOCSize(&mocInfo, &y, &z);
  int ext[6] = {0, y-1, 0, z-1, 0, 0};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext, 6);

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::CanReadFile(const char* fname )
{
  // return 1 for success, return 0 for failure
  vtkWarningMacro("testing read file!");
  MOCInfo mocInfo;
  return this->ParseHeader(fname, &mocInfo);
}

//-----------------------------------------------------------------------------
int vtkMOCReader::ParseHeader(const char* fileName, MOCInfo* mocInfo)
{
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  int retVal = 1;
  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  if(controller->GetLocalProcessId() == 0)
    {
    retVal = this->SingleProcessParseHeader(fileName, mocInfo);
    }
  if(controller->GetNumberOfProcesses() > 1)
    {
    vtkMultiProcessStream data;
    if(controller->GetLocalProcessId() == 0)
      {
      mocInfo->Serialize(data);
      }
    controller->Broadcast(data, 0);
    if(controller->GetLocalProcessId() > 0)
      {
      mocInfo->Deserialize(data);
      }
    }
  timer->StopTimer();
  cerr << controller->GetLocalProcessId() << " has PARSEHEADER time of " << timer->GetElapsedTime() << endl;
  return retVal;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::SingleProcessParseHeader(
  const char* fileName, MOCInfo* mocInfo)
{
  // parses a namelist header file.
  // places all information in mocInfo.
  // return 1 for success, return 0 for failure
  //
  // TODO: support having file paths relative to where the file is?

  ifstream file(fileName);
  std::string line;

  if(file.fail())
    {
    vtkErrorMacro("Error: opening the namelist file, please ensure path is correct\n");
    return 0;
    }

  // first try to find the first namelist with "msf_stuff"
  std::string namelist_name = "msf_stuff";

  bool found = false;
  line = getNextLine(file);
  while(line.compare("") != 0 && !found)
    {
    if(line[0] == '&')
      {
      std::string name = line.substr(1);
      if(name.compare(namelist_name) == 0)
        {
        found = true;
        break;
        }
      }
    line = getNextLine(file);
    }

  if(!found)
    {
    vtkErrorMacro("Error: did not find namelist with name " << namelist_name.c_str());
    return 0;
    }

  // found the right namelist, parse lines until end of namelist block
  std::stringstream line2;
  std::string name;
  int retval;
  while(line.compare("") != 0 && line.compare("/") != 0)
    {
    line = getNextLine(file);
    line2.clear();
    line2.str(line);
    line2 >> name;
    if(name.compare("imt") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->imt;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("jmt") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->jmt;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("km") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->km;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("ysouth_mht") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->ysouth_mht;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("ynorth_mht") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->ynorth_mht;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("dy_mht") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->dy_mht;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("kmt_global_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->kmt_global_file;
      mocInfo->kmt_global_file = mocInfo->kmt_global_file.substr(1, mocInfo->kmt_global_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("kmt_atl_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->kmt_atl_file;
      mocInfo->kmt_atl_file = mocInfo->kmt_atl_file.substr(1, mocInfo->kmt_atl_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("kmt_indopac_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->kmt_indopac_file;
      mocInfo->kmt_indopac_file = mocInfo->kmt_indopac_file.substr(1,mocInfo->kmt_indopac_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("in_depths") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->in_depths;
      mocInfo->in_depths = mocInfo->in_depths.substr(1, mocInfo->in_depths.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("grid_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->grid_file;
      mocInfo->grid_file = mocInfo->grid_file.substr(1, mocInfo->grid_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("u_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->u_file;
      mocInfo->u_file = mocInfo->u_file.substr(1, mocInfo->u_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("v_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->v_file;
      mocInfo->v_file = mocInfo->v_file.substr(1, mocInfo->v_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("do_global") == 0)
      {
      std::string equal, dostr;
      line2 >> equal;
      line2 >> dostr;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      dostr = dostr.substr(1, dostr.length()-2);
      if(dostr[0] == 'f' || dostr[0] == 'F')
        {
        mocInfo->do_global = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        mocInfo->do_global = true;
        }
      }
    if(name.compare("do_atl") == 0)
      {
      std::string equal, dostr;
      line2 >> equal;
      line2 >> dostr;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      dostr = dostr.substr(1, dostr.length()-2);
      if(dostr[0] == 'f' || dostr[0] == 'F')
        {
        mocInfo->do_atl = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        mocInfo->do_atl = true;
        }
      }
    if(name.compare("do_indopac") == 0)
      {
      std::string equal, dostr;
      line2 >> equal;
      line2 >> dostr;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      dostr = dostr.substr(1, dostr.length()-2);
      if(dostr[0] == 'f' || dostr[0] == 'F')
        {
        mocInfo->do_indopac = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        mocInfo->do_indopac = true;
        }
      }
    if(name.compare("byteswap") == 0)
      {
      std::string equal, bytestr;
      line2 >> equal;
      line2 >> bytestr;
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      bytestr = bytestr.substr(1, bytestr.length()-2);
      if(bytestr[0] == 'f' || bytestr[0] == 'F')
        {
        mocInfo->byteswap = false;
        }
      if(bytestr[0] == 't' || bytestr[0] == 'T')
        {
        mocInfo->byteswap = true;
        }
      }
    }

  file.close();
  return 1;
}

//-----------------------------------------------------------------------------
std::string vtkMOCReader::getNextLine(ifstream& file)
{
  // returns the next non-comment line in the file
  while(true)
    {
    std::string line;
    getline(file, line);
    line = trim(line);
    if(line.size() > 0 and line[0] != '!')
      {
      return line;
      }
    if(file.eof())
      {
      // reached end of file
      return std::string("");
      }
    }
}

//-----------------------------------------------------------------------------
const std::string vtkMOCReader::trim(const std::string& pString,
                                     const std::string& pWhitespace)
{
  const size_t beginStr = pString.find_first_not_of(pWhitespace);
  if (beginStr == std::string::npos)
    {
    // no content
    return "";
    }

  const size_t endStr = pString.find_last_not_of(pWhitespace);
  const size_t range = endStr - beginStr + 1;

  return pString.substr(beginStr, range);
}

//-----------------------------------------------------------------------------
int vtkMOCReader::checkParse(std::string& line, std::ios_base::iostate state)
{
  // given a state, see if any error flags have been set.
  // return 1 for good, 0 for any errors found
  if((state & (std::stringstream::failbit | std::stringstream::badbit)) != 0)
    {
    vtkErrorMacro("Error parsing line: " << line.c_str());
    return 0;
    }
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::CalculateMOC(vtkRectilinearGrid* grid, int* ext)
{
  // calculates the MOC.
  // results will be stored as fields in grid.
  // ext is the extents of the final moc array that i'm responsible for.
  // return 1 for success, return 0 for failure

  int rank =
    vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  // printf("proc[%i]: output moc extents [%i %i %i %i %i %i]\n", rank,
  //        ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

  MOCInfo mocInfo;
  int retval = this->ParseHeader(this->FileName, &mocInfo);
  if(retval == 0)
    {
    return 0;
    }

  // figure out local extents of input data fields
  int real_ext3D[6];   // local data extents without ghost cells
  int ext3D[6];        // local data extents with 1 layers of ghost cells
  int ext3D2[6];       // local data extents with 2 layers of ghost cells
  int numberOfProcesses =
    vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses();
  vtkNew<vtkExtentTranslator> translator;
  translator->SetWholeExtent(0, mocInfo.imt-1,
                             0, mocInfo.jmt-1,
                             0, mocInfo.km-1);
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
  translator->GetExtent(real_ext3D);
  translator->GetExtent(ext3D);
  translator->GetExtent(ext3D2);

  // printf("proc[%i]: data extents [%i %i %i %i %i %i]\n", rank, ext3D[0],
  //        ext3D[1], ext3D[2], ext3D[3], ext3D[4], ext3D[5]);

  // the size of the data block i have
  int imt = ext3D[1] - ext3D[0] + 1;
  int jmt = ext3D[3] - ext3D[2] + 1;
  int km = ext3D[5] - ext3D[4] + 1;

  // add ghost cells to x and y dimension. if on the border, need to wrap
  // around. use value of -1 to signify the wrap around
  imt += 2;
  jmt += 2;
  ext3D[0]--;
  if(ext3D[1] == mocInfo.imt-1)
    {
    ext3D[1] = -1;
    }
  else
    {
    ext3D[1]++;
    }
  ext3D[2]--;
  if(ext3D[3] == mocInfo.jmt-1)
    {
    ext3D[3] = -1;
    }
  else
    {
    ext3D[3]++;
    }
  // second layer of ghost cells
  ext3D2[0] = ext3D[0] - 1;
  if(ext3D[1] == mocInfo.imt-1)
    {
    ext3D2[1] = -1;
    }
  else if(ext3D[1] == -1)
    {
    ext3D2[1] = -2;
    }
  else
    {
    ext3D2[1] = ext3D[1] + 1;
    }
  ext3D2[2] = ext3D[2] - 1;
  if(ext3D[3] == mocInfo.jmt-1)
    {
    ext3D2[3] = -1;
    }
  else if(ext3D[3] == -1)
    {
    ext3D2[3] = -2;
    }
  else
    {
    ext3D2[3] = ext3D[3] + 1;
    }

  // size for two layers of ghost cells
  int imt2 = imt + 2;
  int jmt2 = jmt + 2;

  // allocate data that will be read from disk
  Matrix2DDouble uLat(imt, jmt);
  Matrix2DDouble uLong(imt, jmt);
  Matrix2DDouble htn(imt2, jmt2);
  Matrix2DDouble hte(imt2, jmt2);
  Matrix1DFloat dz(km);
  Matrix2DInt global_kmt(imt, jmt);
  Matrix2DInt atl_kmt(imt, jmt);
  Matrix2DInt indopac_kmt(imt, jmt);
  Matrix3DFloat u(imt, jmt, km);
  Matrix3DFloat v(imt, jmt, km);

  // load data from disk
  retval = this->LoadData(&mocInfo, ext3D, ext3D2, imt, jmt, km, uLat, uLong, htn,
                          hte, dz, global_kmt, atl_kmt, indopac_kmt, u, v,
                          imt2, jmt2);

  if(retval == 0)
    {
    // some error occurred
    vtkErrorMacro("An error occurred in the function LoadData");
    return 0;
    }

  Matrix2DFloat dxu(imt, jmt);
  Matrix2DFloat dyu(imt, jmt);
  Matrix2DFloat tarea(imt, jmt);

  // calculate other horizontal grid fields
  retval = this->grid_stuff(htn, hte, dxu, dyu, tarea, imt, jmt, imt2,jmt2);

  // no longer needed
  htn.Clear();
  hte.Clear();

  if(retval == 0)
    {
    // some error occurred
    vtkErrorMacro("An error occurred in the function grid_stuff");
    return 0;
    }

  // convert uLat and uLong to degrees
  for(int i=0; i<imt*jmt; i++)
    {
    uLat(i) = uLat(i) * 180.0 / M_PI;
    uLong(i) = uLong(i) * 180.0 / M_PI;
    }

  Matrix2DFloat tLat(imt, jmt);
  this->sw_4pt(tLat, 0.25, uLat, imt, jmt);

  // not needed anymore
  uLat.Clear();
  uLong.Clear();

  // calculate w from u and v
  Matrix3DFloat w(imt, jmt, km);
  this->wcalc(dz, dxu, dyu, tarea, global_kmt, u, v, w, imt, jmt, km);
  u.Clear();  // not needed anymore

  // set up latitude grid to be used and allocate arrays
  int ny_mht, z;
  this->GetMOCSize(&mocInfo, &ny_mht, &z);

  Matrix1DFloat lat_mht(ny_mht);
  for(int j=0; j<ny_mht; j++)
    {
    lat_mht(j) = mocInfo.ysouth_mht + j * mocInfo.dy_mht;
    }

  Matrix2DFloat psi_temp_old(km, ny_mht);

  // psi holds the final moc values
  // psi is essentially a 3d array
  // psi[][][0] -- global moc
  // psi[][][1] -- atlantic moc
  // psi[][][2] -- indo-pacific moc
  Matrix3DFloat psi(ny_mht, km, 3);

  // calculate overturning streamfunctions

  // first do global
  if(mocInfo.do_global)
    {
    this->moc(&mocInfo, v, w, global_kmt, tLat, dxu, tarea, dz, lat_mht, ny_mht,
              psi_temp_old, ext3D, real_ext3D, imt, jmt);

    // copy values to correct array
    for(int k=0; k<km; k++)
      {
      for(int y=0; y<ny_mht; y++)
        {
        psi(y,k,0) = psi_temp_old(k,y);
        }
      }
    }

  // next do atlantic
  if(mocInfo.do_atl)
    {
    this->moc(&mocInfo, v, w, atl_kmt, tLat, dxu, tarea, dz, lat_mht, ny_mht,
              psi_temp_old, ext3D, real_ext3D, imt, jmt);

    // copy values to correct array
    for(int k=0; k<km; k++)
      {
      for(int y=0; y<ny_mht; y++)
        {
        psi(y,k,1) = psi_temp_old(k,y);
        }
      }
    }

  // now do indo-pacific
  if(mocInfo.do_indopac)
    {
    this->moc(&mocInfo, v, w, indopac_kmt, tLat, dxu, tarea, dz, lat_mht, ny_mht,
              psi_temp_old, ext3D, real_ext3D, imt, jmt);

    // copy values to correct array
    for(int k=0; k<km; k++)
      {
      for(int y=0; y<ny_mht; y++)
        {
        psi(y,k,2) = psi_temp_old(k,y);
        }
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

  // these arrays determine the coordinate positions of the grid
  vtkNew<vtkFloatArray> x_coords;
  for(int i=ext[0]; i<=ext[1]; i++)
    {
    x_coords->InsertNextValue(lat_mht(i));
    }

  vtkNew<vtkFloatArray> y_coords;
  for(int i=ext[2]; i<=ext[3]; i++)
    {
    y_coords->InsertNextValue(depth(i));
    }

  vtkNew<vtkFloatArray> z_coords;
  z_coords->InsertNextValue(0.0);

  grid->SetExtent(ext);
  grid->SetXCoordinates(x_coords.GetPointer());
  grid->SetYCoordinates(y_coords.GetPointer());
  grid->SetZCoordinates(z_coords.GetPointer());

  // TODO: change the array names
  // copy output. remember to only copy the part we want.
  int numPoints = (ext[1]+1) * (ext[3]+1) * (ext[5]+1);
  if(mocInfo.do_global)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_global");
    data->Allocate(numPoints);
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->InsertNextValue(psi(i, j, 0));
        }
      }
    grid->GetPointData()->AddArray(data.GetPointer());
    grid->GetPointData()->SetScalars(data.GetPointer());
    }
  if(mocInfo.do_atl)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_atl");
    data->Allocate(numPoints);
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->InsertNextValue(psi(i, j, 1));
        }
      }
    grid->GetPointData()->AddArray(data.GetPointer());
    if(!mocInfo.do_global)
      {
      grid->GetPointData()->SetScalars(data.GetPointer());
      }
    }
  if(mocInfo.do_indopac)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_indopac");
    data->Allocate(numPoints);
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->InsertNextValue(psi(i, j, 2));
        }
      }
    grid->GetPointData()->AddArray(data.GetPointer());
    if(!mocInfo.do_global && !mocInfo.do_atl)
      {
      grid->GetPointData()->SetScalars(data.GetPointer());
      }
    }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::grid_stuff(Matrix2DDouble& htn,
                             Matrix2DDouble& hte, Matrix2DFloat& dxu,
                             Matrix2DFloat& dyu, Matrix2DFloat& tarea,
                             int imt, int jmt, int imt2, int jmt2)
{
  // create grid-related quantities needed for future calculations.
  // dxu, dyu, and tarea are calculated.

  // allocate memory (local variables)
  Matrix2DFloat dxt(imt,jmt);
  Matrix2DFloat dyt(imt,jmt);

  int i, j;

  // dxu, dyu, dxt, dyt will have a ghost layer of 1
  for(j=1; j<jmt2-1; j++)
    {
    for(i=1; i<imt2-1; i++)
      {
      dxu(i-1,j-1) = 0.5 * (htn(i,j) + htn(i+1,j));
      dyu(i-1,j-1) = 0.5 * (hte(i,j) + hte(i,j+1));
      dxt(i-1,j-1) = 0.5 * (htn(i,j) + htn(i,j-1));
      dyt(i-1,j-1) = 0.5 * (hte(i,j) + hte(i-1,j));
      }
    }

  // tarea has ghost cell layer of 1
  for(i=0; i<imt*jmt; i++)
    {
    tarea(i) = dxt(i) * dyt(i);
    }
  return 1;
}

//-----------------------------------------------------------------------------
void vtkMOCReader::sw_4pt(Matrix2DFloat& xout, float factor, Matrix2DDouble& x,
                          int imt, int jmt)
{
  // perform averaging over a neighborhood
  int i, j;

  for(j=1; j<jmt; j++)
    {
    for(i=1; i<imt; i++)
      {
      xout(i,j) = factor * x(i,j)   +
        factor * x(i,j-1) +
        factor * x(i-1,j) +
        factor * x(i-1,j-1);
      }
    }

  j = 0;
  for(i=1; i<imt; i++)
    {
    xout(i,j) = xout(i,j+1);
    }

  i = 0;
  for(j=1; j<jmt; j++)
    {
    xout(i,j) = factor * x(i,j) +
      factor * x(i,j-1) +
      factor * x(imt-1,j) +
      factor * x(imt-1,j-1);
    }

  xout(0, 0) = xout(1,1);
}

//-----------------------------------------------------------------------------
void vtkMOCReader::wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu,
                         Matrix2DFloat& dyu, Matrix2DFloat& tarea,
                         Matrix2DInt& kmt, Matrix3DFloat& u,
                         Matrix3DFloat& v, Matrix3DFloat& w,
                         int imt, int jmt, int km)
{
  // calculate w ,the vertical velocities, since only u and v are read
  // from file

  int i, j;

  Matrix2DFloat wtk(imt,jmt);
  Matrix2DFloat wtkb(imt,jmt);
  Matrix2DFloat work(imt,jmt);
  Matrix2DFloat fue(imt,jmt);
  Matrix2DFloat fuw(imt,jmt);
  Matrix2DFloat fvn(imt,jmt);
  Matrix2DFloat fvs(imt,jmt);

  for(i=0; i<imt*jmt; i++)
    {
    wtkb(i) = 0.0; // set bottom velocity to zero
    wtk(i) = 0.0;  // initialize
    }

  for(int k=km-1; k>=0; k--)  // integrate from bottom up
    {
    // advection fluxes
    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        work(i,j) = 0.5 * u(i,j,k) * dyu(i,j);
        }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        int j2 = cshift(j, -1, jmt);
        fue(i,j) = work(i,j2);
        fue(i,j) += work(i,j);
        }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        int i2 = cshift(i, -1, imt);
        fuw(i,j) = fue(i2,j);
        }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        work(i,j) = 0.5 * v(i,j,k) * dxu(i,j);
        }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        int i2 = cshift(i, -1, imt);
        fvn(i,j) = work(i2,j);
        fvn(i,j) += work(i,j);
        }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        int j2 = cshift(j, -1, jmt);
        fvs(i,j) = fvn(i,j2);
        }
      }

    // calculate vertical velocity at top of kth level
    // (vertical velocity is zero at bottom of T columns)
    for(i=0; i<imt*jmt; i++)
      {
      work(i) = (fvn(i) - fvs(i) + fue(i) - fuw(i)) / tarea(i);
      }

    for(i=0; i<imt*jmt; i++)
      {
      if(k+1 <= kmt(i))   {
      wtk(i) = wtkb(i) - dz(k) * work(i);
      }
      }

    for(j=0; j<jmt; j++)
      {
      for(i=0; i<imt; i++)
        {
        w(i,j,k) = wtk(i,j);
        }
      }

    // top value becomes bottom value for next pass
    for(i=0; i<imt*jmt; i++)
      {
      wtkb(i) = wtk(i);
      }

    }  // for(int k=km-1; k<=0; k--)
}

//-----------------------------------------------------------------------------
void vtkMOCReader::moc(MOCInfo* mocInfo, Matrix3DFloat& v, Matrix3DFloat& w,
                       Matrix2DInt& kmtb, Matrix2DFloat& tLat, Matrix2DFloat& dxu,
                       Matrix2DFloat& tarea, Matrix1DFloat& dz,
                       Matrix1DFloat& lats, int ny, Matrix2DFloat& psi,
                       int* ext3D, int* real_ext3D, int imt, int jmt)
{
  int i, j, k;
  int ii, jj, local_jj, true_jj;

  // calculates the meridional overturning circulation
  Matrix2DFloat work(imt, jmt);
  Matrix1DFloat psi0(mocInfo->km);

  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  int numprocs = controller->GetNumberOfProcesses();
  int rank = controller->GetLocalProcessId();

  // initialize psi array
  for(i=0; i<mocInfo->km*ny; i++)
    {
    psi(i) = 0.0;
    }
  for(k=0; k<mocInfo->km; k++)
    {
    psi0(k) = 0.0;
    }

  // find j index of southernmost ocean point in basin
  bool found = false;
  for(j=1; j<jmt-1; j++)
    {
    for(i=1; i<imt-1; i++)
      {
      if(kmtb(i,j) != 0)
        {
        jj = j + ext3D[2];
        local_jj = j;
        found = true;
        break;
        }
      }
    if(found)
      {
      break;
      }
    }

  // the minimum j-index is the true j-index
  controller->AllReduce(&jj, &true_jj, 1, vtkCommunicator::MIN_OP);

  // find southern_lat, and collect it at process 0
  float my_southern_lat;
  if(jj == true_jj && real_ext3D[0] == 0)
    {
    // ydeg(jj-1)
    my_southern_lat = 0.5*(tLat(0, local_jj) + tLat(0, local_jj-1));
    }
  else
    {
    my_southern_lat = std::numeric_limits<float>::max();
    }
  float southern_lat;
  controller->Reduce(&my_southern_lat, &southern_lat, 1,
                     vtkCommunicator::MIN_OP, 0);

  // zonal and vertical integration to find streamfunction across southernmost
  // grid circle in basin

  // WORK holds the north-south velocity, modulated by other terms,
  // of one z-slice.
  // psi0(k) holds the sum of all north-south velocities over
  // the latitude jj-1, from depths k to km.
  // everyone computes a local psi0 if necessary, then a reduce is done
  // to get correct psi0. only process 0 will have correct psi0.
  j = local_jj-1;
  if(jj == true_jj)
    {
    for(i=1; i<imt-1; i++)
      {
      work(i,j) = -dz(mocInfo->km-1) * v(i,j,mocInfo->km-1) * dxu(i,j);
      }

    psi0(mocInfo->km-1) = 0.0;
    for(i=1; i<imt-1; i++)
      {
      int j2 = cshift(j, 1, jmt);
      if(kmtb(i,j) == 0 && kmtb(i, j2) > 0)
        {
        psi0(mocInfo->km-1) += work(i,j);
        }
      }
    }

  for(k=mocInfo->km-1; k>=1; k--)
    {
    if(jj == true_jj)
      {
      for(i=1; i<imt-1; i++)
        {
        work(i,j) = -dz(k-1) * v(i,j,k-1) * dxu(i,j);
        }

      for(i=1; i<imt-1; i++)
        {
        int j2 = cshift(j, 1, jmt);
        if(kmtb(i,j) == 0 && kmtb(i, j2) > 0)
          {
          psi0(k-1) += work(i,j);
          }
        }
      }
    psi0(k-1) += psi0(k);
    }
  Matrix1DFloat tempArray(mocInfo->km);
  controller->Reduce(psi0.GetData(), tempArray.GetData(), mocInfo->km,
                     vtkCommunicator::SUM_OP, 0);
  for(k=0; k<mocInfo->km; k++)
    {
    psi0(k) = tempArray(k);
    }
  tempArray.Clear();

  // psi(k,j) holds the moc value of depth k and lat j
  // compute my local moc
  std::vector<moc_point_t> points(imt*jmt);
  for(k=0; k<mocInfo->km; k++)
    {
    // first scan through all points of the layer and find all valid points
    // which are not land. record the latitude of the point as well as its
    // work value.
    int npoints = 0;     // number of points at this layer
    int k2 = k+1;  // the actual depth value that should be used
    for(jj=1; jj<jmt-1; jj++)
      {
      for(ii=1; ii<imt-1; ii++)
        {
        if(k2 <= kmtb(ii,jj))
          {
          points[npoints].work = w(ii,jj,k) * tarea(ii,jj);
          points[npoints].lat = tLat(ii, jj);
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
    for(j=0; j<ny; j++)
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

  // determine parent and children in binary tree.
  // used for compositing mocs together
  int parent = (rank-1) / 2;
  int child[2];
  int numChildren = 0;
  if((rank+1) * 2 - 1 < numprocs)
    {
    child[0] = (rank+1) * 2 - 1;
    numChildren++;
    }
  if((rank+1) * 2 < numprocs)
    {
    child[1] = (rank+1) * 2;
    numChildren++;
    }

  // receive moc from children, and sum up with my moc
  if(numChildren > 0)
    {
    Matrix2DFloat psi_child(mocInfo->km, ny);
    for(i=0; i<numChildren; i++)
      {
      controller->Receive(psi_child.GetData(), ny*mocInfo->km, child[i], 0);
      for(j=0; j<ny*mocInfo->km; j++)
        {
        psi(j) += psi_child(j);
        }
      }
    }

  // send my moc to parent
  if(rank > 0)
    {
    controller->Send(psi.GetData(), ny*mocInfo->km, parent, 0);
    }

  // at this point process 0 should have the accumulated moc,
  // perform some more processing to get final moc.
  if(rank == 0)
    {
    // add in the baseline of the southernmost latitude
    for(j=0; j<ny; j++)
      {
      if(lats(j) >= southern_lat)
        {
        for(k=0; k<mocInfo->km; k++)
          {
          psi(k,j) += psi0(k);
          }
        }
      }

    // smooth grid-point noise over y
    for(j=1; j<ny-1; j++)
      {
      for(k=0; k<mocInfo->km; k++)
        {
        psi(k,j) = 0.25*(psi(k,j-1) + psi(k,j+1)) + 0.5*psi(k,j);
        }
      }

    // normalize to Sv
    for(i=0; i<mocInfo->km*ny; i++)
      {
      psi(i) *= 1e-12;
      }

    float special_value = -1e34;
    // replace any zeroes with the special value
    for(i=0; i<mocInfo->km*ny; i++)
      {
      if(psi(i) == 0.0)
        {
        psi(i) = special_value;
        }
      }
    }

  // process 0 broadcasts results to everyone
  controller->Broadcast(psi.GetData(), ny*mocInfo->km, 0);
}

//-----------------------------------------------------------------------------
int vtkMOCReader::cshift(int i, int offset, int size)
{
  // find the correct index when a circular shift is performed
  // i is current index
  // offset is how many times to shift
  // size is the full size of this array

  int n = (i + offset) % size;
  if(n < 0)
    {
    n += size;
    }
  return n;
}

//-----------------------------------------------------------------------------
void vtkMOCReader::GetMOCSize(MOCInfo* mocInfo, int* y, int* z)
{
  // calculates the final size of the moc array.
  // note that when moc is usually
  // shown, the y (latitude) is shown along the horizontal axis, and the z
  // (depth) is usually shown along the vertical axis

  float temp = ((mocInfo->ynorth_mht - mocInfo->ysouth_mht) / mocInfo->dy_mht);
  int ny_mht = floor(temp + 0.5) + 1;
  *y = ny_mht;
  *z = mocInfo->km;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadData(MOCInfo* mocInfo, int* ext3D, int* ext3D2, int imt,
                           int jmt, int km, Matrix2DDouble& uLat,
                           Matrix2DDouble& uLong,
                           Matrix2DDouble& htn, Matrix2DDouble& hte,
                           Matrix1DFloat& dz, Matrix2DInt& global_kmt,
                           Matrix2DInt& atl_kmt, Matrix2DInt& indopac_kmt,
                           Matrix3DFloat& u, Matrix3DFloat& v, int imt2, int jmt2)
{
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();

  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  // load data from disk

  // uLat
  int offset = mocInfo->imt * mocInfo->jmt * 0;
  this->LoadDataBlock2DDouble(mocInfo, mocInfo->grid_file, ext3D, offset,
                              imt, jmt, uLat);

  // cerr << "imt " << imt << " jmt " << jmt << " filename " << mocInfo->grid_file << " ulatsize "
  //      << uLat.GetSize() << " extents " << ext3D[0] << " " << ext3D[1] << " " << ext3D[2]
  //      << " " << ext3D[3] << endl;

  // uLong
  offset = mocInfo->imt * mocInfo->jmt * 1;
  this->LoadDataBlock2DDouble(mocInfo, mocInfo->grid_file, ext3D, offset,
                              imt, jmt, uLong);

  // htn
  offset = mocInfo->imt * mocInfo->jmt * 2;
  this->LoadDataBlock2DDouble2(mocInfo, mocInfo->grid_file, ext3D2, offset,
                               imt2, jmt2, htn);

  // hte
  offset = mocInfo->imt * mocInfo->jmt * 3;
  this->LoadDataBlock2DDouble2(mocInfo, mocInfo->grid_file, ext3D2, offset,
                               imt2, jmt2, hte);

  // read depths file. read the first number in each row.
  // process 0 reads it and broadcasts the values.
  if(controller->GetLocalProcessId() == 0)
    {
    FILE* f = fopen(mocInfo->in_depths.c_str(), "r");
    if(f == NULL)
      {
      vtkErrorMacro("Error in opening in_depths file: " << mocInfo->in_depths.c_str());
      vtkMultiProcessStream data;
      controller->Broadcast(data, 0); // sending an error message
      return 0;
      }
    char buffer[1024];
    for(int k=0; k<km; k++)
      {
      fgets(buffer, 1024, f);
      char* tok = strtok(buffer, " ");
      dz(k) = atof(tok);
      }
    fclose(f);
    vtkMultiProcessStream data;
    data.Push(dz.GetData(), dz.GetSize());
    controller->Broadcast(data, 0);
    }
  else
    {
    vtkMultiProcessStream data;
    controller->Broadcast(data, 0);
    if(data.Empty() == true)
      { // empty stream indicates error on proc 0 reading data
      return 0;
      }
    float * tmp = NULL;
    unsigned int size = -1;
    data.Pop(tmp, size);
    if(size != dz.GetSize())
      {
      vtkErrorMacro("Bad array size for dz.");
      }
    memcpy(dz.GetData(), tmp, sizeof(float)*size);
    delete []tmp;
    }

  // read kmt files
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_global_file, ext3D, imt, jmt, global_kmt);
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_atl_file, ext3D, imt, jmt, atl_kmt);
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_indopac_file, ext3D, imt, jmt,
                           indopac_kmt);

  // read velocity files
  this->LoadDataBlock3DFloat(mocInfo, mocInfo->u_file, ext3D, imt, jmt, km, u);
  this->LoadDataBlock3DFloat(mocInfo, mocInfo->v_file, ext3D, imt, jmt, km, v);

  if(mocInfo->byteswap)
    {
    // byteswap all arrays
    uLat.ByteSwap();
    uLong.ByteSwap();
    htn.ByteSwap();
    hte.ByteSwap();
    global_kmt.ByteSwap();
    atl_kmt.ByteSwap();
    indopac_kmt.ByteSwap();
    u.ByteSwap();
    v.ByteSwap();
    }

  timer->StopTimer();
  cerr << controller->GetLocalProcessId() << " has LOADDATA time of " << timer->GetElapsedTime() << endl;

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadDataBlock2DDouble(MOCInfo* mocInfo, std::string filename,
                                        int* ext3D, int offset, int imt, int jmt,
                                        Matrix2DDouble& data)
{
  // read a 2D block from a file of doubles

  // local array indices
  int x_start;
  int x_end;
  int y_start;
  int y_end;

  // indices relative to file
  int x_start_file;
  int y_start_file;

  if(ext3D[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3D[0];
    }
  if(ext3D[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3D[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3D[2];
    }
  if(ext3D[3] == -1)
    {
    y_end = jmt - 2;
    }
  else
    {
    y_end = jmt - 1;
    }

  int row_length = x_end - x_start + 1;
  int col_length = y_end - y_start + 1;

  FILE* f = fopen(filename.c_str(), "rb");
  if(f == NULL)
    {
    vtkErrorMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = (offset +
                        (y_start_file + y) * mocInfo->imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(double), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3D[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (offset +
                          (y_start_file + y) * mocInfo->imt +
                          mocInfo->imt-1) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start+y), sizeof(double), 1, f);
      }
    }
  if(ext3D[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (offset +
                          (y_start_file + y) * mocInfo->imt +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, y_start+y), sizeof(double), 1, f);
      }
    }
  if(ext3D[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    int total_offset = (offset +
                        (mocInfo->jmt-1)*mocInfo->imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, 0), sizeof(double), row_length, f);
    }
  if(ext3D[3] == -1)
    {
    // top ghost cells wrap around
    // copy values of the row over
    int total_offset = (offset +
                        0 +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, jmt-1), sizeof(double), row_length, f);
    }

  // take into account corner cases
  if(ext3D[0] == -1 && ext3D[2] == -1)
    {
    // lower left corner needs upper right value
    int total_offset = (offset +
                        (mocInfo->jmt-1)*mocInfo->imt +
                        mocInfo->imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(double), 1, f);
    }
  if(ext3D[0] == -1 && ext3D[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (offset +
                        0 +
                        mocInfo->imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(double), 1, f);
    }
  if(ext3D[1] == -1 && ext3D[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = (offset +
                        (mocInfo->jmt-1)*mocInfo->imt +
                        0) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, 0), sizeof(double), 1, f);
    }
  if(ext3D[1] == -1 && ext3D[3] == -1)
    {
    // upper right corner needs lower left value
    int total_offset = (offset +
                        0 +
                        0) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, jmt-1), sizeof(double), 1, f);
    }
  fclose(f);
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadDataBlock2DDouble2(MOCInfo* mocInfo, std::string filename,
                                         int* ext3D, int offset, int imt, int jmt,
                                         Matrix2DDouble& data)
{
  // read a 2D block from a file of doubles.
  // this one loads 2 layers of ghost cells.

  // local array indices
  int x_start;
  int x_end;
  int y_start;
  int y_end;

  // indices relative to file
  int x_start_file;
  int y_start_file;

  if(ext3D[0] < 0)
    {
    x_start = abs(ext3D[0]);
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3D[0];
    }
  if(ext3D[1] < 0)
    {
    x_end = imt - 1 + ext3D[1];
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3D[2] < 0)
    {
    y_start = abs(ext3D[2]);
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3D[2];
    }
  if(ext3D[3] < 0)
    {
    y_end = jmt - 1 + ext3D[3];
    }
  else
    {
    y_end = jmt - 1;
    }

  int row_length = x_end - x_start + 1;
  int col_length = y_end - y_start + 1;

  FILE* f = fopen(filename.c_str(), "rb");
  if(f == NULL)
    {
    vtkErrorMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = (offset +
                        (y_start_file + y) * mocInfo->imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(double), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3D[0] < 0)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int x=0; x<abs(ext3D[0]); x++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (offset +
                            (y_start_file + y) * mocInfo->imt +
                            mocInfo->imt-1-x) * sizeof(double);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(x_start-1-x, y_start+y), sizeof(double), 1, f);
        }
      }
    }
  if(ext3D[1] < 0)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int x=0; x<abs(ext3D[1]); x++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (offset +
                            (y_start_file + y) * mocInfo->imt +
                            x) * sizeof(double);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(x_end+1+x, y_start+y), sizeof(double), 1, f);
        }
      }
    }
  if(ext3D[2] < 0)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    for(int y=0; y<abs(ext3D[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->jmt-1-y)*mocInfo->imt +
                          x_start_file) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_start-1-y), sizeof(double), row_length, f);
      }
    }
  if(ext3D[3] < 0)
    {
    // top ghost cells wrap around
    // copy values of the row over
    for(int y=0; y<abs(ext3D[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          x_start_file) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_end+1+y), sizeof(double), row_length, f);
      }
    }

  // take into account corner cases
  if(ext3D[0] < 0 && ext3D[2] < 0)
    {
    // lower left corner needs upper right value
    for(int y=0; y<abs(ext3D[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->jmt-1-y)*mocInfo->imt +
                          mocInfo->imt+ext3D[0]) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start-1-y), sizeof(double), abs(ext3D[0]), f);
      }
    }
  if(ext3D[0] < 0 && ext3D[3] < 0)
    {
    // upper left corner needs lower right value
    for(int y=0; y<abs(ext3D[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          mocInfo->imt+ext3D[0]) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_end+1+y), sizeof(double), abs(ext3D[0]), f);
      }
    }
  if(ext3D[1] < 0 && ext3D[2] < 0)
    {
    // lower right corner needs upper left value
    for(int y=0; y<abs(ext3D[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->jmt-1-y)*mocInfo->imt +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt+ext3D[1], y_start-1-y), sizeof(double), abs(ext3D[1]), f);
      }
    }
  if(ext3D[1] < 0 && ext3D[3] < 0)
    {
    // upper right corner needs lower left value
    for(int y=0; y<abs(ext3D[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt+ext3D[1], y_end+1+y),sizeof(double),abs(ext3D[1]),f);
      }
    }
  fclose(f);
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadDataBlock2DInt(MOCInfo* mocInfo, std::string filename, int* ext3D,
                                     int imt, int jmt, Matrix2DInt& data)
{
  // load a 2D block of ints from a file

  // local array indices
  int x_start;
  int x_end;
  int y_start;
  int y_end;

  // indices relative to file
  int x_start_file;
  int y_start_file;

  if(ext3D[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3D[0];
    }
  if(ext3D[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3D[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3D[2];
    }
  if(ext3D[3] == -1)
    {
    y_end = jmt - 2;
    }
  else
    {
    y_end = jmt - 1;
    }

  int row_length = x_end - x_start + 1;
  int col_length = y_end - y_start + 1;

  FILE* f = fopen(filename.c_str(), "rb");
  if(f == NULL)
    {
    vtkErrorMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = ((y_start_file + y) * mocInfo->imt +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(int), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3D[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = ((y_start_file + y) * mocInfo->imt +
                          mocInfo->imt-1) * sizeof(int);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start+y), sizeof(int), 1, f);
      }
    }
  if(ext3D[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = ((y_start_file + y) * mocInfo->imt +
                          0) * sizeof(int);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, y_start+y), sizeof(int), 1, f);
      }
    }
  if(ext3D[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    int total_offset = ((mocInfo->jmt-1)*mocInfo->imt +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, 0), sizeof(int), row_length, f);
    }
  if(ext3D[3] == -1)
    {
    // top ghost cells wrap around
    // copy values of the row over
    int total_offset = (0 +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, jmt-1), sizeof(int), row_length, f);
    }

  // take into account corner cases
  if(ext3D[0] == -1 && ext3D[2] == -1)
    {
    // lower left corner needs upper right value
    int total_offset = ((mocInfo->jmt-1)*mocInfo->imt +
                        mocInfo->imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(int), 1, f);
    }
  if(ext3D[0] == -1 && ext3D[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (0 +
                        mocInfo->imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(int), 1, f);
    }
  if(ext3D[1] == -1 && ext3D[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = ((mocInfo->jmt-1)*mocInfo->imt +
                        0) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, 0), sizeof(int), 1, f);
    }
  if(ext3D[1] == -1 && ext3D[3] == -1)
    {
    // upper right corner needs lower left value
    int total_offset = (0 +
                        0) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, jmt-1), sizeof(int), 1, f);
    }
  fclose(f);
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadDataBlock3DFloat(MOCInfo* mocInfo, std::string filename,
                                       int* ext3D, int imt, int jmt, int km,
                                       Matrix3DFloat& data)
{
  // load a 3D block of floats from a file

  // local array indices
  int x_start;
  int x_end;
  int y_start;
  int y_end;

  // indices relative to file
  int x_start_file;
  int y_start_file;

  if(ext3D[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3D[0];
    }
  if(ext3D[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3D[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3D[2];
    }
  if(ext3D[3] == -1)
    {
    y_end = jmt - 2;
    }
  else
    {
    y_end = jmt - 1;
    }

  int row_length = x_end - x_start + 1;
  int col_length = y_end - y_start + 1;

  FILE* f = fopen(filename.c_str(), "rb");
  if(f == NULL)
    {
    vtkErrorMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int z=0; z<km; z++)
    {
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          (y_start_file + y) * mocInfo->imt +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_start+y, z), sizeof(float), row_length, f);
      }
    }

  // fill in wraparound ghost cells
  if(ext3D[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int z=0; z<km; z++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                            (y_start_file + y) * mocInfo->imt +
                            mocInfo->imt-1) * sizeof(float);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(0, y_start+y, z), sizeof(float), 1, f);
        }
      }
    }
  if(ext3D[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int z=0; z<km; z++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                            (y_start_file + y) * mocInfo->imt +
                            0) * sizeof(float);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(imt-1, y_start+y, z), sizeof(float), 1, f);
        }
      }
    }
  if(ext3D[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          (mocInfo->jmt-1)*mocInfo->imt +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, 0, z), sizeof(float), row_length, f);
      }
    }
  if(ext3D[3] == -1)
    {
    // top ghost cells wrap around
    // copy values of the row over
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          0 +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, jmt-1, z), sizeof(float), row_length, f);
      }
    }

  // take into account corner cases
  if(ext3D[0] == -1 && ext3D[2] == -1)
    {
    // lower left corner needs upper right value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          (mocInfo->jmt-1)*mocInfo->imt +
                          mocInfo->imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, 0, z), sizeof(float), 1, f);
      }
    }
  if(ext3D[0] == -1 && ext3D[3] == -1)
    {
    // upper left corner needs lower right value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          0 +
                          mocInfo->imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, jmt-1, z), sizeof(float), 1, f);
      }
    }
  if(ext3D[1] == -1 && ext3D[2] == -1)
    {
    // lower right corner needs upper left value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          (mocInfo->jmt-1)*mocInfo->imt +
                          0) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, 0, z), sizeof(float), 1, f);
      }
    }
  if(ext3D[1] == -1 && ext3D[3] == -1)
    {
    // upper right corner needs lower left value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->imt * mocInfo->jmt +
                          0 +
                          0) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, jmt-1, z), sizeof(float), 1, f);
      }
    }
  fclose(f);
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::compare_latitude(const void* x, const void* y)
{
  float lat_x = ((moc_point_t*)x)->lat;
  float lat_y = ((moc_point_t*)y)->lat;

  if(lat_x < lat_y)
    {
    return -1;
    }
  if(lat_x > lat_y)
    {
    return 1;
    }
  return 0;
}

//-----------------------------------------------------------------------------
void vtkMOCReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;

}
