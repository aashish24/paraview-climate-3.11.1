/*=========================================================================

  Module:    vtkMHTReader.h

=========================================================================*/
// .NAME vtkMHTReader - read in POP data files and output MHT
// .SECTION Description
// vtkMHTReader is a source object that calculates the meridional heat
// transport. The initial file is a text file listing paths to all the
// required binary data files. The output of this reader is a 1D rectilinear grid
// object
// .SECTION Caveats
// Files may need to be byte swapped.
// .SECTION See Also
// vtkMOCReader

#ifndef __vtkMHTReader_h
#define __vtkMHTReader_h

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
class POPInputInformation;
class InternalMHTWorkArrays;

class VTK_EXPORT vtkMHTReader : public vtkAbstractPOPReader
{
public:
  static vtkMHTReader *New();
  vtkTypeMacro(vtkMHTReader, vtkAbstractPOPReader);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkMHTReader();
  ~vtkMHTReader();
  int RequestInformation(vtkInformation *,
                         vtkInformationVector **,
                         vtkInformationVector *);
  int RequestData(vtkInformation *,
                  vtkInformationVector **,
                  vtkInformationVector *);

  int CalculateMHT(vtkRectilinearGrid* output, int* extMHT);
  void wcalc(Matrix1DFloat& dz, Matrix2DFloat& dxu,
             Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
             Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w,
             int imt, int jmt, int km);
  void wcalc_pbc(Matrix3DFloat& dzu, Matrix2DFloat& dxu,
                 Matrix2DFloat& dyu, Matrix2DFloat& tarea, Matrix2DInt& kmt,
                 Matrix3DFloat& u, Matrix3DFloat& v, Matrix3DFloat& w, int imt, int jmt);

  void meridional_heat(POPInputInformation* popInfo, Matrix2DInt& kmtb, Matrix2DFloat& tLat,
                       Matrix1DFloat& lats, int imt1GL, int jmt1GL, int ny, int jj, bool hasGlobalJIndexMin, float southern_lat,
                       Matrix1DFloat& mht);

private:
  vtkMHTReader(const vtkMHTReader&);          // not implemented
  void operator = (const vtkMHTReader&);   // not implemented
  InternalMHTWorkArrays* MHTWorkArrays;
};
#endif
