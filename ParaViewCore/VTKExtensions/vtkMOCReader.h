/*=========================================================================

  Module:    vtkMOCReader.h

=========================================================================*/
// .NAME vtkMOCReader - read in POP data files and output MOC
// .SECTION Description
// vtkMOCReader is a source object that calculates the meridional overturning
// circulation. The initial file is a text file listing paths to all the
// required binary data files. The output of this reader is a 2D rectilinear grid
// object
// .SECTION Caveats
// Files may need to be byte swapped.
// .SECTION See Also
// vtkMHTReader

#ifndef __vtkMOCReader_h
#define __vtkMOCReader_h

#include <string>
#include <iostream>

#include "vtkAbstractPOPReader.h"

class vtkRectilinearGrid;
class Matrix1DFloat;
class Matrix1DInt;
class Matrix2DDouble;
class Matrix2DFloat;
class Matrix2DInt;
class Matrix3DFloat;
class InternalMHTWorkArrays;

class VTK_EXPORT vtkMOCReader : public vtkAbstractPOPReader
{
public:
  static vtkMOCReader *New();
  vtkTypeMacro(vtkMOCReader, vtkAbstractPOPReader);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkMOCReader();
  ~vtkMOCReader();
  int RequestInformation(vtkInformation *,
                         vtkInformationVector **,
                         vtkInformationVector *);
  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *);

  int CalculateMOC(vtkRectilinearGrid* output, int* ext);
  void wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu,
             Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
             Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w,
             int imt, int jmt, int km);
  void wcalc_pbc(Matrix3DFloat& dzu, Matrix2DFloat& dxu,
                 Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
                 Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w, int imt, int jmt);

  void moc(POPInputInformation* popInfo, Matrix3DFloat& v, Matrix3DFloat& w,
           Matrix2DInt& kmtb, Matrix2DFloat& tLat, Matrix2DFloat& dxu,
           Matrix2DFloat& tarea, Matrix1DFloat& dz, Matrix3DFloat& dzu, Matrix1DFloat& lats,
           int ny, int local_jj, bool has_global_jj, float southern_lat,
           int imt, int jmt, Matrix2DFloat& psi);

private:
  vtkMOCReader(const vtkMOCReader&);          // not implemented
  void operator = (const vtkMOCReader&);   // not implemented
};
#endif
