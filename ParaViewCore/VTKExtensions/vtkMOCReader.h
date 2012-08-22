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


protected:
  vtkMOCReader();
  ~vtkMOCReader();
  int RequestInformation(vtkInformation *,
                         vtkInformationVector **,
                         vtkInformationVector *);
  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *);

  int ParseHeader(const char* fileName, MOCInfo* mocinfo);
  int SingleProcessParseHeader(const char* fileName, MOCInfo* mocinfo);
  std::string getNextLine(ifstream& file);
  int checkParse(std::string& line, std::ios_base::iostate state);
  const std::string trim(const std::string& pstring, const std::string& pWhitespace = " \t");

  int CalculateMOC(vtkRectilinearGrid* grid, int* ext);
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
  void moc(MOCInfo* mocinfo, Matrix3DFloat& v, Matrix3DFloat& w,
           Matrix2DInt& kmtb, Matrix2DFloat& tLat, Matrix2DFloat& dxu,
           Matrix2DFloat& tarea, Matrix1DFloat& dz, Matrix1DFloat& lats,
           int ny, Matrix2DFloat& psi, int* ext3D, int* real_ext3D,
           int imt, int jmt);
  void GetMOCSize(MOCInfo*, int* y, int* z);
  int cshift(int i, int offset, int size);

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

  // Description:
  // the namelist file containing all input variables
  char* FileName;

private:
  vtkMOCReader(const vtkMOCReader&);          // not implemented
  void operator = (const vtkMOCReader&);   // not implemented
};



#endif
