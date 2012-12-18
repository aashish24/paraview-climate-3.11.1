/*=========================================================================

  Module:    vtkMOCReader.h

=========================================================================*/
// .NAME vtkMOCReader - read in MOC data file
// .SECTION Description
// vtkMOCReader is a source object that calculates the meridional overturning
// circulation. The initial file is a text file listing paths to all the
// required binary data files. The output of this reader is a 2D Imagedata
// object
// .SECTION Caveats
// Files may need to be byte swapped.
// .SECTION See Also
// vtkMOCReader

#ifndef __vtkMOCReader_h
#define __vtkMOCReader_h

#include <string>
#include <iostream>

#include "vtkRectilinearGridAlgorithm.h"

class vtkRectilinearGrid;
class Matrix1DFloat;
class Matrix1DInt;
class Matrix2DDouble;
class Matrix2DFloat;
class Matrix2DInt;
class Matrix3DFloat;
class MOCInfo;
class InternalMHTWorkArrays;

class VTK_EXPORT vtkMOCReader : public vtkRectilinearGridAlgorithm
{
public:
  static vtkMOCReader *New();
  vtkTypeMacro(vtkMOCReader, vtkRectilinearGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // set and get the namelist header filename
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // check file for suitability to this reader
  int CanReadFile(const char* fname);

  static int cshift(int i, int offset, int size);

protected:
  vtkMOCReader();
  ~vtkMOCReader();
  int RequestInformation(vtkInformation *,
                         vtkInformationVector **,
                         vtkInformationVector *);
  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *);

  int ParseMetaFile(const char* fileName, MOCInfo* mocinfo);
  int SingleProcessParseMetaFile(const char* fileName, MOCInfo* mocinfo);
  std::string getNextLine(ifstream& file);
  int checkParse(std::string& line, std::ios_base::iostate state);
  const std::string trim(const std::string& pstring, const std::string& pWhitespace = " \t");

  int CalculateMOC(vtkRectilinearGrid* mocGrid, vtkRectilinearGrid* mhtGrid, int* ext);
  int grid_stuff(Matrix2DDouble& htn,
                 Matrix2DDouble& hte, Matrix2DFloat& dxu,
                 Matrix2DFloat& dyu, Matrix2DFloat& tarea, int imt, int jmt,
                 int imt2, int jmt2);
  void sw_4pt(Matrix2DFloat& xout, float factor, Matrix2DDouble& x,
              int imt, int jmt);
  void wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu,
             Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
             Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w,
             int imt, int jmt, int km);
  void wcalc_pbc(Matrix3DFloat& dzu, Matrix2DFloat& dxu,
                 Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
                 Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w, int imt, int jmt);

  void moc(MOCInfo* mocInfo, Matrix3DFloat& v, Matrix3DFloat& w,
           Matrix2DInt& kmtb, Matrix2DFloat& tLat, Matrix2DFloat& dxu,
           Matrix2DFloat& tarea, Matrix1DFloat& dz, Matrix3DFloat& dzu, Matrix1DFloat& lats,
           int ny, int local_jj, bool has_global_jj, float southern_lat,
           int imt, int jmt, Matrix2DFloat& psi);

  void meridional_heat(MOCInfo* mocInfo, Matrix2DInt& kmtb, Matrix2DFloat& tLat,
                       Matrix1DFloat& lats, int ny, int jj, float southern_lat,
                       Matrix1DFloat& mht);

  void GetMOCSize(MOCInfo*, int* y, int* z);

  int LoadData(MOCInfo* mocinfo, int* ext3D, int* ext3D2, int imt, int jmt,
               int km, Matrix2DDouble& uLat,
               Matrix2DDouble& uLong, Matrix2DDouble& htn, Matrix2DDouble& hte,
               Matrix1DFloat& dz, Matrix2DInt& global_kmt,
               Matrix2DInt& atl_kmt, Matrix2DInt& indopac_kmt,
               Matrix3DFloat& u, Matrix3DFloat& v, int imt2, int jmt2);
  int LoadDataBlock2DDouble(MOCInfo* mocinfo, std::string filename, int* ext3D,
                            int offset, int imt, int jmt, Matrix2DDouble& data);
  int LoadDataBlock2DDouble2(MOCInfo* mocinfo, std::string filename, int* ext3D,
                             int offset, int imt, int jmt,Matrix2DDouble& data);
  int LoadDataBlock2DInt(MOCInfo* mocinfo, std::string filename, int* ext3D,
                         int imt, int jmt, Matrix2DInt& data);
  int LoadDataBlock3DFloat(MOCInfo* mocinfo, std::string filename, int* ext3D,
                           int imt, int jmt, int km, Matrix3DFloat& data);

  static int compare_latitude(const void* x, const void* y);

  void FindSouthern(int imt, int jmt, int* ext3D, int* real_ext3D, Matrix2DInt& kmtb,
                    Matrix2DFloat& tLat, int* local_jj, bool* has_global_jj, float* southern_lat);

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

private:
  vtkMOCReader(const vtkMOCReader&);          // not implemented
  void operator = (const vtkMOCReader&);   // not implemented

  // Description:
  // the namelist file containing all input variables
  char* FileName;
  InternalMHTWorkArrays* MHTWorkArrays;

  // Description:
  // Values that correspond to imt=0 and jmt=0.
  int GlobalIMT0;
  int GlobalJMT0;
};



#endif
