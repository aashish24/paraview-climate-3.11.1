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
// a point, used for an optimization in computing moc and mht
  struct moc_mht_point_t
  {
    float lat;   // latitude of point
    float work;  // quantity being considered
    moc_mht_point_t(float l, float w) : lat(l), work(w)
    {}
  };

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

  int seekFile(FILE* f, int first_record, int imt, int jmt)
  {
    // move the file pointer to the correct position for the given record.
    // have to be careful since the amount to move may be larger than
    // the long int used in fseek.
    // after the function, the file pointer should be in the correct position.
    // returns 0 if successful, otherwise returns a non-zero value.

    int ret_val;
    long long int offset;
    offset = (long long int)imt*jmt*(first_record-1)*sizeof(float);

    ret_val = fseek(f, 0, SEEK_SET);
    if(ret_val != 0)
      {
      return ret_val;
      }

    while(offset > 0)
      {
      long int chunk;
      if(offset > std::numeric_limits<long int>::max())
        {
        chunk = std::numeric_limits<long int>::max();
        }
      else
        {
        chunk = offset;
        }
      offset = offset - chunk;
      ret_val = fseek(f, chunk, SEEK_CUR);
      if(ret_val != 0)
        {
        return ret_val;
        }
      }

    return 0;
  }
}

// matrix classes used in the MOC code

//-----------------------------------------------------
//  1D float array
//-----------------------------------------------------

class Matrix1DFloat {
public:
  Matrix1DFloat();
  Matrix1DFloat(unsigned int xDim);
  void Allocate(unsigned int xDim);
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
Matrix1DFloat::Matrix1DFloat()
{
  this->XDim = 0;
  this->Data = NULL;
}

inline
Matrix1DFloat::Matrix1DFloat(unsigned int xDim)
{
  this->XDim = xDim;
  this->Data = new float[xDim];
}

inline
void Matrix1DFloat::Allocate(unsigned int xDim)
{
  if(this->XDim != xDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->Data = new float[xDim];
    }
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
  Matrix2DInt();
  Matrix2DInt(unsigned int xDim, unsigned int yDim);
  void Allocate(unsigned int xDim, unsigned int yDim);
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
Matrix2DInt::Matrix2DInt()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

inline
Matrix2DInt::Matrix2DInt(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new int[xDim * yDim];
}

inline
void Matrix2DInt::Allocate(unsigned int xDim, unsigned int yDim)
{
  if(this->XDim != xDim || this->YDim != yDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->YDim = yDim;
    this->Data = new int[xDim * yDim];
    }
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
  Matrix2DFloat();
  Matrix2DFloat(unsigned int xDim, unsigned int yDim);
  void Allocate(unsigned int xDim, unsigned int yDim);
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
Matrix2DFloat::Matrix2DFloat()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

inline
Matrix2DFloat::Matrix2DFloat(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new float[xDim * yDim];
}

inline
void Matrix2DFloat::Allocate(unsigned int xDim, unsigned int yDim)
{
  if(xDim != this->XDim || yDim != this->YDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->YDim = yDim;
    this->Data = new float[xDim * yDim];
    }
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
  Matrix2DDouble();
  Matrix2DDouble(unsigned int xDim, unsigned int yDim);
  void Allocate(unsigned int xDim, unsigned int yDim);
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
Matrix2DDouble::Matrix2DDouble()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

inline
Matrix2DDouble::Matrix2DDouble(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new double[xDim * yDim];
}

inline
void Matrix2DDouble::Allocate(unsigned int xDim, unsigned int yDim)
{
  if(this->XDim != xDim || this->YDim != yDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->YDim = yDim;
    this->Data = new double[xDim * yDim];
    }
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
  Matrix3DFloat();
  Matrix3DFloat(unsigned int xDim, unsigned int yDim, unsigned int zDim);
  void Allocate(unsigned int xDim, unsigned int yDim, unsigned int zDim);
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
Matrix3DFloat::Matrix3DFloat()
{
  this->XDim = 0;
  this->YDim = 0;
  this->ZDim = 0;

  this->Data = NULL;
}

inline
Matrix3DFloat::Matrix3DFloat(unsigned int xDim, unsigned int yDim, unsigned int zDim)
{
  this->XDim = xDim;
  this->YDim = yDim;
  this->ZDim = zDim;

  this->Data = new float[xDim * yDim * zDim];
}

inline
void Matrix3DFloat::Allocate(unsigned int xDim, unsigned int yDim, unsigned int zDim)
{
  if(this->XDim != xDim || this->YDim != yDim || this->ZDim != zDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->YDim = yDim;
    this->ZDim = zDim;
    this->Data = new float[xDim * yDim * zDim];
    }
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
    this->global_imt = 320;
    this->global_jmt = 384;
    this->global_km = 40;
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
    this->uet_file = "";
    this->vnt_file = "";
    this->u_first_record = 0;
    this->v_first_record = 0;
    this->uet_first_record = 0;
    this->vnt_first_record = 0;
    this->do_global = false;
    this->do_atl = false;
    this->do_indopac = false;
    this->do_msf = false;
    this->do_mht = false;
    this->use_pbc = false;
    this->byteswap = false;
  }

  // store the internal variable values in data to be streamed
  // to other processes.
  void Serialize(vtkMultiProcessStream& data)
  {
    data.Reset();
    data << this->global_imt << this->global_jmt << this->global_km;
    data << this->ysouth_mht << this->ynorth_mht << this->dy_mht;
    data << this->kmt_global_file << this->kmt_atl_file << this->kmt_indopac_file;
    data << this->in_depths << this->grid_file << this->u_file << this->v_file;
    data << this->uet_file << this->vnt_file;
    data << this->u_first_record << this->v_first_record;
    data << this->uet_first_record << this->vnt_first_record;
    data << static_cast<int>(this->do_global) << static_cast<int>(this->do_atl);
    data << static_cast<int>(this->do_indopac) << static_cast<int>(this->do_msf);
    data << static_cast<int>(this->do_mht) << static_cast<int>(this->use_pbc);
    data << static_cast<int>(this->byteswap);
  }

  // set the values from data in the internal variables.
  void Deserialize(vtkMultiProcessStream& data)
  {
    data >> this->global_imt >> this->global_jmt >> this->global_km;
    data >> this->ysouth_mht >> this->ynorth_mht >> this->dy_mht;
    data >> this->kmt_global_file >> this->kmt_atl_file >> this->kmt_indopac_file;
    data >> this->in_depths >> this->grid_file >> this->u_file >> this->v_file;
    data >> this->uet_file >> this->vnt_file;
    data >> this->u_first_record >> this->v_first_record;
    data >> this->uet_first_record >> this->vnt_first_record;
    int tmp;
    data >> tmp;
    this->do_global = static_cast<bool>(tmp);
    data >> tmp;
    this->do_atl = static_cast<bool>(tmp);
    data >> tmp;
    this->do_indopac = static_cast<bool>(tmp);
    data >> tmp;
    this->do_msf = static_cast<bool>(tmp);
    data >> tmp;
    this->do_mht = static_cast<bool>(tmp);
    data >> tmp;
    this->use_pbc = static_cast<bool>(tmp);
    data >> tmp;
    this->byteswap = static_cast<bool>(tmp);
  }

  int global_imt;                     // size of grid in x dimension
  int global_jmt;                     // size of grid in y dimension
  int global_km;                      // size of grid in z dimension
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
  std::string uet_file;             // for MHT computation
  std::string vnt_file;             // for MHT computation
  int u_first_record;           // where the u-velocities start in the file
  int v_first_record;           // where the v-velocities start in the file
  int uet_first_record;         // where the u-velocities start in the file
  int vnt_first_record;         // where the v-velocities start in the file
  bool do_global;              // compute global quantities
  bool do_atl;                 // compute atlantic quantities
  bool do_indopac;             // compute indian-pacific quantities
  bool do_msf;                 // compute meridional overturning circulation
  bool do_mht;                 // compute meridional heat transport
  bool use_pbc;                // use partial bottom cells
  bool byteswap;               // byteswap the binary input files
};


double getvalacb(Matrix1DFloat& m, int i)
{
  cerr << m(i) << endl;
  return m(i);
}

double getvalacb(Matrix2DFloat& m, int i, int j)
{
  cerr << m(i,j) << endl;
  return m(i,j);
}

double getvalacb(Matrix2DDouble& m, int i, int j)
{
  cerr << m(i,j) << endl;
  return m(i,j);
}

double getvalacb(Matrix3DFloat& m, int i, int j, int k)
{
  cerr << m(i,j, k) << endl;
  return m(i,j, k);
}

class InternalMHTWorkArrays
{
// The default constructor and destructor should be sufficient.
public:
  // currenlty Work1 and WorkY are in global indices with no ghost layers
  Matrix2DFloat Work1;
  Matrix2DFloat WorkY;
  void Clear()
  {
    this->Work1.Clear();
    this->WorkY.Clear();
  }
  void Compute(int globalIimt, int globalJmt, int km, bool use_pbc, Matrix3DFloat& uet,
               Matrix3DFloat& vnt, Matrix2DFloat& tarea1GL, Matrix3DFloat& dzt,
               Matrix1DFloat& dz);
  bool IsComputed()
  {
    return (this->Work1.GetData() != NULL);
  }
};

void InternalMHTWorkArrays::Compute(
  int globalImt, int globalJmt, int km, bool use_pbc, Matrix3DFloat& uet,
  Matrix3DFloat& vnt, Matrix2DFloat& tarea1GL, Matrix3DFloat& dzt,
  Matrix1DFloat& dz)
{
  this->Work1.Allocate(globalImt, globalJmt);
  this->WorkY.Allocate(globalImt, globalJmt);

  Matrix2DFloat workX(globalImt, globalJmt);

  for(int i=0; i<globalJmt*globalImt; i++)
    {
    workX(i) = 0.0;
    this->WorkY(i) = 0.0;
    }

  getvalacb(uet, 1, 1, 1);
  getvalacb(tarea1GL, 1, 1);
  // vertical integration of workX and this->WorkY
  for(int k=0; k < km; k++)
    {
    for (int j=0; j<globalJmt; j++)
      {
      for (int i=0; i<globalImt; i++)
        {
        if(use_pbc)
          {
          workX(i,j) += uet(i,j,k)*tarea1GL(i+1,j+1)*dzt(i,j,k)*4.186e-15;
          this->WorkY(i,j) += vnt(i,j,k)*tarea1GL(i+1,j+1)*dzt(i,j,k)*4.186e-15;
          }
        else
          {
          workX(i,j) += uet(i,j,k)*tarea1GL(i+1,j+1)*dz(k)*4.186e-15;
          this->WorkY(i,j) += vnt(i,j,k)*tarea1GL(i+1,j+1)*dz(k)*4.186e-15;
          }
        }
      }
    }

  // find divergence of vertically-integrated heat transport
  for(int i=0; i<globalImt; i++)
    {
    for(int j=0; j<globalJmt; j++)
      {
      int i2 = vtkMOCReader::cshift(i, -1, globalImt);
      int j2 = vtkMOCReader::cshift(j, -1, globalJmt);
      this->Work1(i,j) = workX(i,j) - workX(i2,j) + this->WorkY(i,j) - this->WorkY(i,j2);
      }
    }
}

vtkStandardNewMacro(vtkMOCReader);

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
vtkMOCReader::vtkMOCReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(2);
  this->MHTWorkArrays = new InternalMHTWorkArrays;
  this->GlobalIMT0 = VTK_INT_MIN;
  this->GlobalJMT0 = VTK_INT_MIN;
}

//-----------------------------------------------------------------------------
vtkMOCReader::~vtkMOCReader()
{
  this->SetFileName(0);
  if(this->MHTWorkArrays)
    {
    delete this->MHTWorkArrays;
    this->MHTWorkArrays = NULL;
    }
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

  vtkRectilinearGrid* mocGrid = vtkRectilinearGrid::GetData(outputVector);
  vtkRectilinearGrid* mhtGrid = vtkRectilinearGrid::GetData(outputVector, 1);
  if(this->CalculateMOC(mocGrid, mhtGrid, ext) == 0)
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
  int retval = this->ParseMetaFile(this->FileName, &mocInfo);
  if(retval == 0)
    {
    return 0;
    }

  int ny_mht, z;
  this->GetMOCSize(&mocInfo, &ny_mht, &z);
  int ext[6] = {0, ny_mht-1, 0, z-1, 0, 0};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), ext, 6);

  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::CanReadFile(const char* fname )
{
  // return 1 for success, return 0 for failure
  MOCInfo mocInfo;
  return this->ParseMetaFile(fname, &mocInfo);
}

//-----------------------------------------------------------------------------
int vtkMOCReader::ParseMetaFile(const char* fileName, MOCInfo* mocInfo)
{
  int retVal = 1;
  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  if(controller->GetLocalProcessId() == 0)
    {
    retVal = this->SingleProcessParseMetaFile(fileName, mocInfo);
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
  return retVal;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::SingleProcessParseMetaFile(
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
      line2 >> mocInfo->global_imt;
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
      line2 >> mocInfo->global_jmt;
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
      line2 >> mocInfo->global_km;
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
    if(name.compare("do_msf") == 0)
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
        mocInfo->do_msf = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        mocInfo->do_msf = true;
        }
      }
    if(name.compare("do_mht") == 0)
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
        mocInfo->do_mht = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        mocInfo->do_mht = true;
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
    if(name.compare("uet_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->uet_file;
      mocInfo->uet_file = mocInfo->uet_file.substr(1, mocInfo->uet_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("vnt_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> mocInfo->vnt_file;
      mocInfo->vnt_file = mocInfo->vnt_file.substr(1, mocInfo->vnt_file.length()-2);
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
int vtkMOCReader::CalculateMOC(vtkRectilinearGrid* mocGrid, vtkRectilinearGrid* mhtGrid, int* ext)
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
  int retval = this->ParseMetaFile(this->FileName, &mocInfo);
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
  translator->SetWholeExtent(0, mocInfo.global_imt-1,
                             0, mocInfo.global_jmt-1,
                             0, mocInfo.global_km-1);
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
  // we only partition in logical x and y so km should be equal to mocInfo.global_km
  int km = ext3D[5] - ext3D[4] + 1;

  // add ghost cells to x and y dimension. if on the border, need to wrap
  // around. use value of -1 to signify the wrap around
  imt1GL += 2;
  jmt1GL += 2;
  ext3D1GL[0]--;
  if(ext3D1GL[1] == mocInfo.global_imt-1)
    {
    ext3D1GL[1] = -1;
    }
  else
    {
    ext3D1GL[1]++;
    }
  ext3D1GL[2]--;
  if(ext3D1GL[3] == mocInfo.global_jmt-1)
    {
    ext3D1GL[3] = -1;
    }
  else
    {
    ext3D1GL[3]++;
    }
  // second layer of ghost cells
  ext3D2GL[0] = ext3D1GL[0] - 1;
  if(ext3D1GL[1] == mocInfo.global_imt-1)
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
  if(ext3D1GL[3] == mocInfo.global_jmt-1)
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
  retval = this->LoadData(&mocInfo, ext3D1GL, ext3D2GL, imt1GL, jmt1GL, km, uLat1GL, uLong1GL, htn2GL,
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
  // mht holds the final mht arrays
  Matrix2DFloat mht(ny_mht, 3);






  int localJIndexMin1GL = -1; // was local_jj
  bool hasGlobalJIndexMin = false;
  float southern_lat = -1000;
  this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, global_kmt1GL, tLat1GL,
                     &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);

  // clear out temporary MHT work arrays since we don't know if they're
  // valid.
  this->MHTWorkArrays->Clear();
  if(mocInfo.do_mht && (mocInfo.do_global || mocInfo.do_atl || mocInfo.do_indopac) )
    {
    Matrix3DFloat uet(ext3D[1]-ext3D[0]+1, ext3D[3]-ext3D[2]+1, km);
    Matrix3DFloat vnt(ext3D[1]-ext3D[0]+1, ext3D[3]-ext3D[2]+1, km);
    FILE* f = fopen(mocInfo.uet_file.c_str(), "rb");
    if(f == NULL)
      {
      vtkErrorMacro("Error in opening uet_file: " << mocInfo.uet_file);
      return 0;
      }
    retval = seekFile(f, mocInfo.uet_first_record, ext3D[1]-ext3D[0]+1, ext3D[3]-ext3D[2]+1);
    if(retval != 0)
      {
      vtkErrorMacro("Error during file seek for uet_file: " << mocInfo.uet_file);
      return 0;
      }
    fread(uet.GetData(), sizeof(float), (ext3D[1]-ext3D[0]+1)*(ext3D[3]-ext3D[2]+1)*km, f);
    fclose(f);

    // read in vnt
    f = fopen(mocInfo.vnt_file.c_str(), "rb");
    if(f == NULL)
      {
      vtkErrorMacro("Error in opening vnt_file: " << mocInfo.vnt_file);
      return 0;
      }
    retval = seekFile(f, mocInfo.vnt_first_record, ext3D[1]-ext3D[0]+1, ext3D[3]-ext3D[2]+1);
    if(retval != 0)
      {
      vtkErrorMacro("Error during file seek for vnt_file: " << mocInfo.vnt_file);
      return 0;
      }
    fread(vnt.GetData(), sizeof(float), (ext3D[1]-ext3D[0]+1)*(ext3D[3]-ext3D[2]+1)*km, f);
    fclose(f);

    if(mocInfo.byteswap)
      {
      uet.ByteSwap();
      vnt.ByteSwap();
      }

    Matrix3DFloat dzt(mocInfo.global_imt, mocInfo.global_jmt, mocInfo.global_km); // acbauer still needs to read in dzt which is used with partial bottom cells only
    this->MHTWorkArrays->Compute(mocInfo.global_imt, mocInfo.global_jmt, mocInfo.global_km,
                                 mocInfo.use_pbc, uet, vnt, tarea1GL, dzt, dz);
    }

  // mht_tmp is the output array from meridional_heat() method
  Matrix1DFloat mht_temp(ny_mht);


  // dzu is for partial bottom cells
  Matrix3DFloat dzu1GL(imt1GL, jmt1GL, km);



// calculate overturning streamfunctions

// first do global
  if(mocInfo.do_global)
    {
    if(mocInfo.do_msf)
      {
      this->moc(&mocInfo, v1GL, w1GL, global_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,0) = psi_temp_old(k,y);
          }
        }
      } // mocInfo.do_msf
    if(mocInfo.do_mht)
      {
      if(mocInfo.do_msf == false)
        {
        this->FindSouthern(imt1GL, jmt1GL, ext3D1GL, ext3D, global_kmt1GL, tLat1GL,
                           &localJIndexMin1GL, &hasGlobalJIndexMin, &southern_lat);
        }
      // acbauer -- need to still figure out proper jIndexMin1GL value
      int jIndexMin1GL = localJIndexMin1GL;  // j index of southernmost ocean point in basin

      this->meridional_heat(&mocInfo, global_kmt1GL, tLat1GL, lat_mht, ny_mht, jIndexMin1GL, southern_lat, mht_temp);
      for(int y=0; y<ny_mht; y++)
        {
        mht(y,0) = mht_temp(y);
        }
      }
    }

// next do atlantic
  if(mocInfo.do_atl)
    {
    if(mocInfo.do_msf)
      {
      this->moc(&mocInfo, v1GL, w1GL, atl_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,1) = psi_temp_old(k,y);
          }
        }
      } // mocInfo.do_msf
    if(mocInfo.do_mht)
      {
      vtkWarningMacro("still need to do atl mht");
      }
    }

// now do indo-pacific
  if(mocInfo.do_indopac)
    {
    if(mocInfo.do_msf)
      {
      this->moc(&mocInfo, v1GL, w1GL, indopac_kmt1GL, tLat1GL, dxu1GL, tarea1GL, dz, dzu1GL, lat_mht,
                ny_mht, localJIndexMin1GL, hasGlobalJIndexMin, southern_lat, imt1GL, jmt1GL, psi_temp_old);

      // copy values to correct array
      for(int k=0; k<km; k++)
        {
        for(int y=0; y<ny_mht; y++)
          {
          psi(y,k,2) = psi_temp_old(k,y);
          }
        }
      } // mocInfo.do_msf
    if(mocInfo.do_mht)
      {
      vtkWarningMacro("still need to do indopac mht");
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

  mocGrid->SetExtent(ext);
  mocGrid->SetXCoordinates(x_coords.GetPointer());
  mocGrid->SetYCoordinates(y_coords.GetPointer());
  mocGrid->SetZCoordinates(z_coords.GetPointer());

// TODO: change the array names
// copy output. remember to only copy the part we want.
  if(mocInfo.do_global)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_global");
    data->SetNumberOfTuples(mocGrid->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 0));
        counter++;
        //cerr << "moc["<<i<<","<<j<<"] = " << psi(i, j, 0) << endl;
        }
      }
    mocGrid->GetPointData()->AddArray(data.GetPointer());
    mocGrid->GetPointData()->SetScalars(data.GetPointer());
    }
  if(mocInfo.do_atl)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_atl");
    data->SetNumberOfTuples(mocGrid->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 1));
        counter++;
        }
      }
    mocGrid->GetPointData()->AddArray(data.GetPointer());
    if(!mocInfo.do_global)
      {
      mocGrid->GetPointData()->SetScalars(data.GetPointer());
      }
    }
  if(mocInfo.do_indopac)
    {
    vtkNew<vtkFloatArray> data;
    data->SetName("reader_moc_indopac");
    data->SetNumberOfTuples(mocGrid->GetNumberOfPoints());
    vtkIdType counter = 0;
    for(int j=ext[2]; j<=ext[3]; j++)
      {
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(counter, psi(i, j, 2));
        counter++;
        }
      }
    mocGrid->GetPointData()->AddArray(data.GetPointer());
    if(!mocInfo.do_global && !mocInfo.do_atl)
      {
      mocGrid->GetPointData()->SetScalars(data.GetPointer());
      }
    }

  // next do the mht grid which is a 1D grid
  mhtGrid->SetExtent(ext[0], ext[1], 0, 0, 0, 0);
  mhtGrid->SetXCoordinates(x_coords.GetPointer());
  mhtGrid->SetYCoordinates(z_coords.GetPointer()); // a single point
  mhtGrid->SetZCoordinates(z_coords.GetPointer()); // a single point
  if(mocInfo.do_mht)
    {
    if(mocInfo.do_global)
      {
      vtkNew<vtkFloatArray> data;
      data->SetName("reader_mht_global");
      data->SetNumberOfTuples(mhtGrid->GetNumberOfPoints());
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(i-ext[0], mht(i, 0));
        cerr << "mht["<<i<<"] = " << mht(i, 0) << " at " << x_coords->GetValue(i) << endl;
        }
      mhtGrid->GetPointData()->AddArray(data.GetPointer());
      mhtGrid->GetPointData()->SetScalars(data.GetPointer());
      }
    if(mocInfo.do_atl)
      {
      vtkNew<vtkFloatArray> data;
      data->SetName("reader_mht_atl");
      data->SetNumberOfTuples(mhtGrid->GetNumberOfPoints());
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(i-ext[0], mht(i, 1));
        }
      mhtGrid->GetPointData()->AddArray(data.GetPointer());
      if(!mocInfo.do_global)
        {
        mhtGrid->GetPointData()->SetScalars(data.GetPointer());
        }
      }
    if(mocInfo.do_indopac)
      {
      vtkNew<vtkFloatArray> data;
      data->SetName("reader_mht_indopac");
      data->SetNumberOfTuples(mhtGrid->GetNumberOfPoints());
      for(int i=ext[0]; i<=ext[1]; i++)
        {
        data->SetValue(i-ext[0], mht(i, 2));
        }
      mhtGrid->GetPointData()->AddArray(data.GetPointer());
      if(!mocInfo.do_global && !mocInfo.do_atl)
        {
        mhtGrid->GetPointData()->SetScalars(data.GetPointer());
        }
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

  // for(i=1;i<imt1GL-1;i++)
  //   for(j=1;j<jmt1GL-1;j++)
  //     {
  //     cerr << "noGL w[" << i-1 << "," << j-1 << ",0] = " << w1GL(i,j, 0) << endl;
  //     }
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
void vtkMOCReader::moc(MOCInfo* mocInfo, Matrix3DFloat& v1GL, Matrix3DFloat& w1GL,
                       Matrix2DInt& kmtb1GL, Matrix2DFloat& tLat1GL, Matrix2DFloat& dxu1GL,
                       Matrix2DFloat& tarea1GL, Matrix1DFloat& dz, Matrix3DFloat& dzu1GL,
                       Matrix1DFloat& lats, int ny_mht, int localJIndexMin1GL, bool hasGlobalJIndexMin,
                       float southern_lat, int imt1GL, int jmt1GL, Matrix2DFloat& psi)
{
  // calculates the meridional overturning circulation
  Matrix2DFloat work1GL(imt1GL, jmt1GL);
  Matrix1DFloat psi0(mocInfo->global_km);
  // initialize psi array
  for(int i=0; i<mocInfo->global_km*ny_mht; i++)
    {
    psi(i) = 0.0;
    }
  for(int k=0; k<mocInfo->global_km; k++)
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
      if(mocInfo->use_pbc)
        {
        work1GL(i,localJIndexMin1GL-1) = -dzu1GL(i,localJIndexMin1GL,mocInfo->global_km-1) * v1GL(i,localJIndexMin1GL-1,mocInfo->global_km-1) * dxu1GL(i,localJIndexMin1GL-1);
        }
      else
        {
        work1GL(i,localJIndexMin1GL-1) = -dz(mocInfo->global_km-1) * v1GL(i,localJIndexMin1GL-1,mocInfo->global_km-1) * dxu1GL(i,localJIndexMin1GL-1);
        }
      }

    psi0(mocInfo->global_km-1) = 0.0;
    for(int i=1; i<imt1GL-1; i++)
      {
      int j2 = cshift(localJIndexMin1GL-1, 1, jmt1GL);
      if(kmtb1GL(i,localJIndexMin1GL-1) == 0 && kmtb1GL(i, j2) > 0)
        {
        psi0(mocInfo->global_km-1) += work1GL(i,localJIndexMin1GL-1);
        }
      }
    }

  for(int k=mocInfo->global_km-1; k>=1; k--)
    {
    if(hasGlobalJIndexMin)
      {
      for(int i=1; i<imt1GL-1; i++)
        {
        if(mocInfo->use_pbc)
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
  Matrix1DFloat tempArray(mocInfo->global_km);
  controller->Reduce(psi0.GetData(), tempArray.GetData(), mocInfo->global_km,
                     vtkCommunicator::SUM_OP, 0);

  for(int k=0; k<mocInfo->global_km; k++)
    {
    psi0(k) = tempArray(k);
    }
  tempArray.Clear();

  // psi(k,j) holds the moc value of depth k and lat j
  // compute my local moc
  std::vector<moc_point_t> points(imt1GL*jmt1GL);
  for(int k=0; k<mocInfo->global_km; k++)
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
  Matrix2DFloat psi_temp(mocInfo->global_km, ny_mht);
  controller->Reduce(psi.GetData(), psi_temp.GetData(), ny_mht*mocInfo->global_km, vtkCommunicator::SUM_OP, 0);

  // at this point process 0 should have the accumulated moc,
  // perform some more processing to get final moc.
  if(rank == 0)
    {
    for(int j=0; j<ny_mht; j++)
      {
      for(int k=0; k<mocInfo->global_km; k++)
        {
        psi(k,j) = psi_temp(k,j);
        }
      }
    // add in the baseline of the southernmost latitude
    for(int j=0; j<ny_mht; j++)
      {
      if(lats(j) >= southern_lat)
        {
        for(int k=0; k<mocInfo->global_km; k++)
          {
          psi(k,j) += psi0(k);
          }
        }
      }

    // smooth grid-point noise over y
    for(int j=1; j<ny_mht-1; j++)
      {
      for(int k=0; k<mocInfo->global_km; k++)
        {
        psi(k,j) = 0.25*(psi(k,j-1) + psi(k,j+1)) + 0.5*psi(k,j);
        }
      }

    // normalize to Sv
    for(int i=0; i<mocInfo->global_km*ny_mht; i++)
      {
      psi(i) *= 1e-12;
      }

    float special_value = -1e34;
    // replace any zeroes with the special value
    for(int i=0; i<mocInfo->global_km*ny_mht; i++)
      {
      if(psi(i) == 0.0)
        {
        psi(i) = special_value;
        }
      }
    }

  // process 0 broadcasts results to everyone
  controller->Broadcast(psi.GetData(), ny_mht*mocInfo->global_km, 0);
}

//-----------------------------------------------------------------------------
void vtkMOCReader::meridional_heat(
  MOCInfo* mocInfo,
  Matrix2DInt& kmtb1GL, Matrix2DFloat& tLat1GL,
  Matrix1DFloat& lat_mht, int ny_mht, int jIndexMin1GL, float southern_lat,
  Matrix1DFloat& mht)
{
  for(int i=0; i<ny_mht; i++)
    {
    mht(i) = 0.0;
    }

  // for(int i=0;i<3600;i++)
  //   for(int j=0;j<2400;j++)
  //     cerr << "workY[" << i << "," << j << "] = " << this->MHTWorkArrays->WorkY(i,j) << endl;

  printf("southernmost j with 1 ghost layer = %i\n", jIndexMin1GL);
  printf("southernmost lat = %f\n", southern_lat);

  // zonal (over longitude range) integration to find heat transport
  // across southernmost grid circle in basin
  float global_mht0 = 0.0;
  int j1GL = jIndexMin1GL-1; // this is jj-1 from original code
  int j2_1GL = cshift(j1GL, 1, mocInfo->global_jmt+2);
  for(int i=0; i<mocInfo->global_imt; i++)
    {
    if(kmtb1GL(i+1,j1GL) == 0 && kmtb1GL(i+1,j2_1GL) > 0)
      {
      global_mht0 += this->MHTWorkArrays->WorkY(i,j1GL-1);
      cerr << "WorkY[" << i << "," << j1GL-1 << "] = " << this->MHTWorkArrays->WorkY(i,j1GL-1) << endl;
      }
    }

  printf("mht0: %f for imt_local of %d jmt_local %d jIndexMin1GL %d j2_1gl %d\n", global_mht0, mocInfo->global_imt,
         mocInfo->global_jmt, jIndexMin1GL, j2_1GL);

  // optimized way to find mht
  // scan through all points of the layer and find all valid points
  // which are not land. record the latitude of the point as well as its work
  // value
  moc_point_t* points = new moc_point_t[mocInfo->global_imt*mocInfo->global_jmt];
  int npoints = 0;
  for(int i=0;i<mocInfo->global_imt;i++)
    {
    for(int j=0;j<mocInfo->global_jmt;j++)
      {
      if(kmtb1GL(i+1, j+1) > 0)
        {
        points[npoints].work = this->MHTWorkArrays->Work1(i,j);
        points[npoints].lat = tLat1GL(i+1,j+1);
        npoints++;
        }
      }
    }

  // sort all points by their latitude
  qsort(points, npoints, sizeof(moc_mht_point_t), compare_latitude);

  // step through latitude from the bottom up, accumulating all points with are
  // less than or equal to the current latitude. keep track of where in the
  // points array the accumulation has gone to. start each search through the
  // points at the point where we left off before.
  int index = 0; // the current place of the points array
  for(int j=0; j<ny_mht; j++)
    {
    if(j > 0)
      {
      mht(j) = mht(j-1);
      }
    float lats_j = lat_mht(j);
    while(index < npoints && points[index].lat < lats_j)
      {
      mht(j) += points[index].work;
      index++;
      }
    }

  for(int j=0; j<ny_mht; j++)
    {
    if(lat_mht(j) > southern_lat)
      {
      mht(j) += global_mht0;
      }
    }

  delete[] points;
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
void vtkMOCReader::GetMOCSize(MOCInfo* mocInfo, int* ny_mht, int* z)
{
  // calculates the final size of the moc array.
  // note that when moc is usually
  // shown, the y (latitude) is shown along the horizontal axis, and the z
  // (depth) is usually shown along the vertical axis

  float temp = ((mocInfo->ynorth_mht - mocInfo->ysouth_mht) / mocInfo->dy_mht);
  *ny_mht = floor(temp + 0.5) + 1;
  *z = mocInfo->global_km;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadData(MOCInfo* mocInfo, int* ext3D1GL, int* ext3D2GL, int imt,
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
  int offset = mocInfo->global_imt * mocInfo->global_jmt * 0;
  this->LoadDataBlock2DDouble(mocInfo, mocInfo->grid_file, ext3D1GL, offset,
                              imt, jmt, uLat);

  // uLong
  offset = mocInfo->global_imt * mocInfo->global_jmt * 1;
  this->LoadDataBlock2DDouble(mocInfo, mocInfo->grid_file, ext3D1GL, offset,
                              imt, jmt, uLong);

  // htn
  offset = mocInfo->global_imt * mocInfo->global_jmt * 2;
  this->LoadDataBlock2DDouble2(mocInfo, mocInfo->grid_file, ext3D2GL, offset,
                               imt2, jmt2, htn);

  // hte
  offset = mocInfo->global_imt * mocInfo->global_jmt * 3;
  this->LoadDataBlock2DDouble2(mocInfo, mocInfo->grid_file, ext3D2GL, offset,
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
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_global_file, ext3D1GL, imt, jmt, global_kmt);
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_atl_file, ext3D1GL, imt, jmt, atl_kmt);
  this->LoadDataBlock2DInt(mocInfo, mocInfo->kmt_indopac_file, ext3D1GL, imt, jmt,
                           indopac_kmt);

  // read velocity files
  this->LoadDataBlock3DFloat(mocInfo, mocInfo->u_file, ext3D1GL, imt, jmt, km, u);
  this->LoadDataBlock3DFloat(mocInfo, mocInfo->v_file, ext3D1GL, imt, jmt, km, v);

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
                                        int* ext3DArbitrary, int offset, int imt, int jmt,
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

  if(ext3DArbitrary[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3DArbitrary[0];
    }
  if(ext3DArbitrary[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3DArbitrary[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3DArbitrary[2];
    }
  if(ext3DArbitrary[3] == -1)
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
                        (y_start_file + y) * mocInfo->global_imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(double), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3DArbitrary[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (offset +
                          (y_start_file + y) * mocInfo->global_imt +
                          mocInfo->global_imt-1) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start+y), sizeof(double), 1, f);
      }
    }
  if(ext3DArbitrary[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (offset +
                          (y_start_file + y) * mocInfo->global_imt +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, y_start+y), sizeof(double), 1, f);
      }
    }
  if(ext3DArbitrary[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    int total_offset = (offset +
                        (mocInfo->global_jmt-1)*mocInfo->global_imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, 0), sizeof(double), row_length, f);
    }
  if(ext3DArbitrary[3] == -1)
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
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower left corner needs upper right value
    int total_offset = (offset +
                        (mocInfo->global_jmt-1)*mocInfo->global_imt +
                        mocInfo->global_imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(double), 1, f);
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (offset +
                        0 +
                        mocInfo->global_imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(double), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = (offset +
                        (mocInfo->global_jmt-1)*mocInfo->global_imt +
                        0) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, 0), sizeof(double), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[3] == -1)
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
                                         int* ext3DArbitrary, int offset, int imt, int jmt,
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

  if(ext3DArbitrary[0] < 0)
    {
    x_start = abs(ext3DArbitrary[0]);
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3DArbitrary[0];
    }
  if(ext3DArbitrary[1] < 0)
    {
    x_end = imt - 1 + ext3DArbitrary[1];
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3DArbitrary[2] < 0)
    {
    y_start = abs(ext3DArbitrary[2]);
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3DArbitrary[2];
    }
  if(ext3DArbitrary[3] < 0)
    {
    y_end = jmt - 1 + ext3DArbitrary[3];
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
                        (y_start_file + y) * mocInfo->global_imt +
                        x_start_file) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(double), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3DArbitrary[0] < 0)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int x=0; x<abs(ext3DArbitrary[0]); x++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (offset +
                            (y_start_file + y) * mocInfo->global_imt +
                            mocInfo->global_imt-1-x) * sizeof(double);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(x_start-1-x, y_start+y), sizeof(double), 1, f);
        }
      }
    }
  if(ext3DArbitrary[1] < 0)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int x=0; x<abs(ext3DArbitrary[1]); x++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (offset +
                            (y_start_file + y) * mocInfo->global_imt +
                            x) * sizeof(double);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(x_end+1+x, y_start+y), sizeof(double), 1, f);
        }
      }
    }
  if(ext3DArbitrary[2] < 0)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    for(int y=0; y<abs(ext3DArbitrary[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->global_jmt-1-y)*mocInfo->global_imt +
                          x_start_file) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_start-1-y), sizeof(double), row_length, f);
      }
    }
  if(ext3DArbitrary[3] < 0)
    {
    // top ghost cells wrap around
    // copy values of the row over
    for(int y=0; y<abs(ext3DArbitrary[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          x_start_file) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_end+1+y), sizeof(double), row_length, f);
      }
    }

  // take into account corner cases
  if(ext3DArbitrary[0] < 0 && ext3DArbitrary[2] < 0)
    {
    // lower left corner needs upper right value
    for(int y=0; y<abs(ext3DArbitrary[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->global_jmt-1-y)*mocInfo->global_imt +
                          mocInfo->global_imt+ext3DArbitrary[0]) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start-1-y), sizeof(double), abs(ext3DArbitrary[0]), f);
      }
    }
  if(ext3DArbitrary[0] < 0 && ext3DArbitrary[3] < 0)
    {
    // upper left corner needs lower right value
    for(int y=0; y<abs(ext3DArbitrary[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          mocInfo->global_imt+ext3DArbitrary[0]) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_end+1+y), sizeof(double), abs(ext3DArbitrary[0]), f);
      }
    }
  if(ext3DArbitrary[1] < 0 && ext3DArbitrary[2] < 0)
    {
    // lower right corner needs upper left value
    for(int y=0; y<abs(ext3DArbitrary[2]); y++)
      {
      int total_offset = (offset +
                          (mocInfo->global_jmt-1-y)*mocInfo->global_imt +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt+ext3DArbitrary[1], y_start-1-y), sizeof(double), abs(ext3DArbitrary[1]), f);
      }
    }
  if(ext3DArbitrary[1] < 0 && ext3DArbitrary[3] < 0)
    {
    // upper right corner needs lower left value
    for(int y=0; y<abs(ext3DArbitrary[3]); y++)
      {
      int total_offset = (offset +
                          y +
                          0) * sizeof(double);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt+ext3DArbitrary[1], y_end+1+y),sizeof(double),abs(ext3DArbitrary[1]),f);
      }
    }
  fclose(f);
  return 1;
}

//-----------------------------------------------------------------------------
int vtkMOCReader::LoadDataBlock2DInt(MOCInfo* mocInfo, std::string filename, int* ext3DArbitrary,
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

  if(ext3DArbitrary[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3DArbitrary[0];
    }
  if(ext3DArbitrary[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3DArbitrary[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3DArbitrary[2];
    }
  if(ext3DArbitrary[3] == -1)
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
    int total_offset = ((y_start_file + y) * mocInfo->global_imt +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, y_start+y), sizeof(int), row_length, f);
    }

  // fill in wraparound ghost cells
  if(ext3DArbitrary[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = ((y_start_file + y) * mocInfo->global_imt +
                          mocInfo->global_imt-1) * sizeof(int);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, y_start+y), sizeof(int), 1, f);
      }
    }
  if(ext3DArbitrary[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = ((y_start_file + y) * mocInfo->global_imt +
                          0) * sizeof(int);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, y_start+y), sizeof(int), 1, f);
      }
    }
  if(ext3DArbitrary[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    int total_offset = ((mocInfo->global_jmt-1)*mocInfo->global_imt +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, 0), sizeof(int), row_length, f);
    }
  if(ext3DArbitrary[3] == -1)
    {
    // top ghost cells wrap around
    // copy values of the row over
    int total_offset = (0 +
                        x_start_file) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(x_start, jmt-1), sizeof(int), row_length, f);
    }

  // take into account corner cases
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower left corner needs upper right value
    int total_offset = ((mocInfo->global_jmt-1)*mocInfo->global_imt +
                        mocInfo->global_imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(int), 1, f);
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (0 +
                        mocInfo->global_imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(int), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = ((mocInfo->global_jmt-1)*mocInfo->global_imt +
                        0) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(imt-1, 0), sizeof(int), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[3] == -1)
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
                                       int* ext3DArbitrary, int imt, int jmt, int km,
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

  if(ext3DArbitrary[0] == -1)
    {
    x_start = 1;
    x_start_file = 0;
    }
  else
    {
    x_start = 0;
    x_start_file = ext3DArbitrary[0];
    }
  if(ext3DArbitrary[1] == -1)
    {
    x_end = imt - 2;
    }
  else
    {
    x_end = imt - 1;
    }
  if(ext3DArbitrary[2] == -1)
    {
    y_start = 1;
    y_start_file = 0;
    }
  else
    {
    y_start = 0;
    y_start_file = ext3DArbitrary[2];
    }
  if(ext3DArbitrary[3] == -1)
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
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          (y_start_file + y) * mocInfo->global_imt +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, y_start+y, z), sizeof(float), row_length, f);
      }
    }

  // fill in wraparound ghost cells
  if(ext3DArbitrary[0] == -1)
    {
    // left ghost cells wrap around
    // copy values of the column over
    for(int z=0; z<km; z++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                            (y_start_file + y) * mocInfo->global_imt +
                            mocInfo->global_imt-1) * sizeof(float);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(0, y_start+y, z), sizeof(float), 1, f);
        }
      }
    }
  if(ext3DArbitrary[1] == -1)
    {
    // right ghost cells wrap around
    // copy values of the column over
    for(int z=0; z<km; z++)
      {
      for(int y=0; y<col_length; y++)
        {
        // get to correct place in file
        int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                            (y_start_file + y) * mocInfo->global_imt +
                            0) * sizeof(float);
        fseek(f, total_offset, SEEK_SET);
        fread(&data(imt-1, y_start+y, z), sizeof(float), 1, f);
        }
      }
    }
  if(ext3DArbitrary[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          (mocInfo->global_jmt-1)*mocInfo->global_imt +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, 0, z), sizeof(float), row_length, f);
      }
    }
  if(ext3DArbitrary[3] == -1)
    {
    // top ghost cells wrap around
    // copy values of the row over
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          0 +
                          x_start_file) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(x_start, jmt-1, z), sizeof(float), row_length, f);
      }
    }

  // take into account corner cases
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower left corner needs upper right value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          (mocInfo->global_jmt-1)*mocInfo->global_imt +
                          mocInfo->global_imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, 0, z), sizeof(float), 1, f);
      }
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          0 +
                          mocInfo->global_imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, jmt-1, z), sizeof(float), 1, f);
      }
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
                          (mocInfo->global_jmt-1)*mocInfo->global_imt +
                          0) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, 0, z), sizeof(float), 1, f);
      }
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper right corner needs lower left value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * mocInfo->global_imt * mocInfo->global_jmt +
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
void vtkMOCReader::FindSouthern(int imt1GL, int jmt1GL, int* ext3D1GL, int* ext3D,
                                Matrix2DInt& kmtb1GL, Matrix2DFloat& tLat1GL, int* localJIndexMin1GL,
                                bool* hasGlobalJIndexMin, float* southern_lat)
{
  // need to initialize these variables since on some processes
  // they may not get set before they're used
  *localJIndexMin1GL = VTK_INT_MAX;

  // find j index of southernmost ocean point in basin
  bool found = false;
  for(int j=1; j<jmt1GL-1; j++)
    {
    for(int i=1; i<imt1GL-1; i++)
      {
      if(kmtb1GL(i,j) != 0)
        {
        *localJIndexMin1GL = j;
        found = true;
        break;
        }
      }
    if(found)
      {
      break;
      }
    }

  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();

  // the minimum j-index is the true j-index
  int tempJ = *localJIndexMin1GL + ext3D1GL[2];
  int globalJIndexMin1GL = VTK_INT_MAX;
  controller->AllReduce(&tempJ, &globalJIndexMin1GL, 1, vtkCommunicator::MIN_OP);
  *hasGlobalJIndexMin = (tempJ == globalJIndexMin1GL);

  // find southern_lat, and collect it at process 0
  float my_southern_lat;
  if(*hasGlobalJIndexMin)
    {
    my_southern_lat = 0.5*(tLat1GL(0, *localJIndexMin1GL+ext3D[2]) + tLat1GL(0, *localJIndexMin1GL-1+ext3D[2]));
    }
  else
    {
    my_southern_lat = std::numeric_limits<float>::max();
    }
  controller->Reduce(&my_southern_lat, southern_lat, 1,
                     vtkCommunicator::MIN_OP, 0);
}

//-----------------------------------------------------------------------------
void vtkMOCReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;
}
