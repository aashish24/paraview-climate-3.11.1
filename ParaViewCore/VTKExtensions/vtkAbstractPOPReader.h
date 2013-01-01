/*=========================================================================

  Module:    vtkAbstractPOPReader.h

=========================================================================*/
// .NAME vtkAbstractPOPReader - abstract reader for POP raw data
// .SECTION Description
// .SECTION Caveats
// Files may need to be byte swapped.

#ifndef __vtkAbstractPOPReader_h
#define __vtkAbstractPOPReader_h

#include <string>
#include <iostream>

#include "vtkRectilinearGridAlgorithm.h"
#include "vtkMultiProcessStream.h"

class vtkRectilinearGrid;
class Matrix1DFloat;
class Matrix1DInt;
class Matrix2DDouble;
class Matrix2DFloat;
class Matrix2DInt;
class Matrix3DFloat;
class InternalMHTWorkArrays;

// Description:
// holds all configuration information for an MOC calculation
class POPInputInformation {

public:
  POPInputInformation()
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
    this->pbc_file = "";
    this->grid_file = "";
    this->u_file = "";
    this->v_file = "";
    this->uet_file = "";
    this->vnt_file = "";
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
    data << this->in_depths << this->pbc_file << this->grid_file << this->u_file << this->v_file;
    data << this->uet_file << this->vnt_file;
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
    data >> this->in_depths >> this->pbc_file >> this->grid_file >> this->u_file >> this->v_file;
    data >> this->uet_file >> this->vnt_file;
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
  std::string pbc_file;             // partial bottom cell file
  std::string grid_file;            // grid information
  std::string u_file;               // u-velocities
  std::string v_file;               // v-velocities
  std::string uet_file;             // for MHT computation
  std::string vnt_file;             // for MHT computation
  bool do_global;              // compute global quantities
  bool do_atl;                 // compute atlantic quantities
  bool do_indopac;             // compute indian-pacific quantities
  bool do_msf;                 // compute meridional overturning circulation
  bool do_mht;                 // compute meridional heat transport
  bool use_pbc;                // use partial bottom cells
  bool byteswap;               // byteswap the binary input files
};



class VTK_EXPORT vtkAbstractPOPReader : public vtkRectilinearGridAlgorithm
{
public:
  vtkTypeMacro(vtkAbstractPOPReader, vtkRectilinearGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // set and get the namelist header filename
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // check file for suitability to this reader
  int CanReadFile(const char* fname);

  static int cshift(int i, int offset, int size);

  // Description:
  static int LoadDataBlock2DDouble(POPInputInformation* popinfo, std::string filename, int* ext3D,
                                   int offset, int imt, int jmt, Matrix2DDouble& data);
  static int LoadDataBlock2DDouble2(POPInputInformation* popinfo, std::string filename, int* ext3D,
                                    int offset, int imt, int jmt,Matrix2DDouble& data);
  static int LoadDataBlock2DInt(POPInputInformation* popinfo, std::string filename, int* ext3D,
                                int imt, int jmt, Matrix2DInt& data);
  static int LoadDataBlock3DFloat(POPInputInformation* popinfo, std::string filename, int* ext3D,
                                  int imt, int jmt, int km, Matrix3DFloat& data);
protected:
  vtkAbstractPOPReader();
  ~vtkAbstractPOPReader();

  int ParseMetaFile(const char* fileName, POPInputInformation* popinfo);
  int SingleProcessParseMetaFile(const char* fileName, POPInputInformation* popinfo);
  std::string getNextLine(ifstream& file);
  int checkParse(std::string& line, std::ios_base::iostate state);
  const std::string trim(const std::string& pstring, const std::string& pWhitespace = " \t");

  void GetMOCSize(POPInputInformation*, int* y, int* z);

  int LoadData(POPInputInformation* popinfo, int* ext3D, int* ext3D2, int imt, int jmt,
               int km, Matrix2DDouble& uLat,
               Matrix2DDouble& uLong, Matrix2DDouble& htn, Matrix2DDouble& hte,
               Matrix1DFloat& dz, Matrix2DInt& global_kmt,
               Matrix2DInt& atl_kmt, Matrix2DInt& indopac_kmt,
               Matrix3DFloat& u, Matrix3DFloat& v, int imt2, int jmt2);

  static int compare_latitude(const void* x, const void* y);

  void FindSouthern(int imt, int jmt, int* ext3D, int* real_ext3D, Matrix2DInt& kmtb,
                    Matrix2DFloat& tLat, int* local_jj, bool* has_global_jj, float* southern_lat);

  // Description:
  // Computes tarea with one layer of ghost cells.
  int grid_stuff(Matrix2DDouble& htn1GL,
                 Matrix2DDouble& hte1GL, Matrix2DFloat& dxu1GL,
                 Matrix2DFloat& dyu1GL, Matrix2DFloat& tarea1GL, int imt1GL, int jmt1GL,
                 int imt2GL, int jmt2GL);

  void sw_4pt(Matrix2DFloat& xout, float factor, Matrix2DDouble& x,
              int imt, int jmt);

  // Description:
  // Conversions from global to local indexing. For serial the values will be the same.
  int GetLocalIMT(int global_imt)
  {
    return global_imt - this->GlobalIMT0;
  }
  int GetLocalJMT(int global_jmt)
  {
    return global_jmt - this->GlobalJMT0;
  }

  // Description:
  // Conversions from local to global indexing. For serial the values will be the same.
  int GetGlobalIMT(int local_imt)
  {
    return local_imt+this->GlobalIMT0;
  }
  int GetGlobalJMT(int local_jmt)
  {
    return local_jmt+this->GlobalJMT0;
  }

  // Description:
  // Given a latitude and an array of latitude values, return the
  // index of the smallest value of lat_mht(index) such that
  // latitude < lat_mht(index). Returns -1 if the index couldn't
  // be found.
  int GetLatitudeIndex(Matrix1DFloat& lat_mht, float latitude);

  // Description:
  // Values that correspond to imt=0 and jmt=0.
  int GlobalIMT0;
  int GlobalJMT0;

private:
  vtkAbstractPOPReader(const vtkAbstractPOPReader&);          // not implemented
  void operator = (const vtkAbstractPOPReader&);   // not implemented

  // Description:
  // the namelist file containing all input variables
  char* FileName;
};

#endif
