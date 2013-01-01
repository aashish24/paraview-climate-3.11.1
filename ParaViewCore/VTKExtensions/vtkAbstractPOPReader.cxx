/*=========================================================================

  Module:    vtkAbstractPOPReader.cxx

  =========================================================================*/

#include "vtkAbstractPOPReader.h"

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
#include "vtkPOPMatrices.h"
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

//-----------------------------------------------------------------------------
vtkAbstractPOPReader::vtkAbstractPOPReader()
{
  this->SetNumberOfInputPorts(0);
  this->FileName = 0;
  this->GlobalIMT0 = VTK_INT_MIN;
  this->GlobalJMT0 = VTK_INT_MIN;
}

//-----------------------------------------------------------------------------
vtkAbstractPOPReader::~vtkAbstractPOPReader()
{
  this->SetFileName(0);
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::CanReadFile(const char* fname )
{
  // return 1 for success, return 0 for failure
  POPInputInformation popinfo;
  return this->ParseMetaFile(fname, &popinfo);
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::ParseMetaFile(const char* fileName, POPInputInformation* popinfo)
{
  int retVal = 1;
  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  if(controller->GetLocalProcessId() == 0)
    {
    retVal = this->SingleProcessParseMetaFile(fileName, popinfo);
    }
  if(controller->GetNumberOfProcesses() > 1)
    {
    vtkMultiProcessStream data;
    if(controller->GetLocalProcessId() == 0)
      {
      popinfo->Serialize(data);
      }
    controller->Broadcast(data, 0);
    if(controller->GetLocalProcessId() > 0)
      {
      popinfo->Deserialize(data);
      }
    }
  return retVal;
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::SingleProcessParseMetaFile(
  const char* fileName, POPInputInformation* popinfo)
{
  // parses a namelist header file.
  // places all information in popinfo.
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
      line2 >> popinfo->global_imt;
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
      line2 >> popinfo->global_jmt;
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
      line2 >> popinfo->global_km;
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
      line2 >> popinfo->ysouth_mht;
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
      line2 >> popinfo->ynorth_mht;
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
      line2 >> popinfo->dy_mht;
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
        popinfo->do_msf = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->do_msf = true;
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
        popinfo->do_mht = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->do_mht = true;
        }
      }
    if(name.compare("kmt_global_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> popinfo->kmt_global_file;
      popinfo->kmt_global_file = popinfo->kmt_global_file.substr(1, popinfo->kmt_global_file.length()-2);
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
      line2 >> popinfo->kmt_atl_file;
      popinfo->kmt_atl_file = popinfo->kmt_atl_file.substr(1, popinfo->kmt_atl_file.length()-2);
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
      line2 >> popinfo->kmt_indopac_file;
      popinfo->kmt_indopac_file = popinfo->kmt_indopac_file.substr(1,popinfo->kmt_indopac_file.length()-2);
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
      line2 >> popinfo->in_depths;
      popinfo->in_depths = popinfo->in_depths.substr(1, popinfo->in_depths.length()-2);
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
      line2 >> popinfo->grid_file;
      popinfo->grid_file = popinfo->grid_file.substr(1, popinfo->grid_file.length()-2);
      retval = checkParse(line, line2.rdstate());
      if(retval == 0)
        {
        return 0;
        }
      }
    if(name.compare("pbc_file") == 0)
      {
      std::string equal;
      line2 >> equal;
      line2 >> popinfo->pbc_file;
      popinfo->pbc_file = popinfo->pbc_file.substr(1, popinfo->pbc_file.length()-2);
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
      line2 >> popinfo->u_file;
      popinfo->u_file = popinfo->u_file.substr(1, popinfo->u_file.length()-2);
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
      line2 >> popinfo->v_file;
      popinfo->v_file = popinfo->v_file.substr(1, popinfo->v_file.length()-2);
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
      line2 >> popinfo->uet_file;
      popinfo->uet_file = popinfo->uet_file.substr(1, popinfo->uet_file.length()-2);
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
      line2 >> popinfo->vnt_file;
      popinfo->vnt_file = popinfo->vnt_file.substr(1, popinfo->vnt_file.length()-2);
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
        popinfo->do_global = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->do_global = true;
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
        popinfo->do_atl = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->do_atl = true;
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
        popinfo->do_indopac = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->do_indopac = true;
        }
      }
    if(name.compare("use_pbc") == 0)
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
        popinfo->use_pbc = false;
        }
      if(dostr[0] == 't' || dostr[0] == 'T')
        {
        popinfo->use_pbc = true;
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
        popinfo->byteswap = false;
        }
      if(bytestr[0] == 't' || bytestr[0] == 'T')
        {
        popinfo->byteswap = true;
        }
      }
    }

  file.close();
  return 1;
}

//-----------------------------------------------------------------------------
std::string vtkAbstractPOPReader::getNextLine(ifstream& file)
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
const std::string vtkAbstractPOPReader::trim(const std::string& pString,
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
int vtkAbstractPOPReader::checkParse(std::string& line, std::ios_base::iostate state)
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
int vtkAbstractPOPReader::grid_stuff(
  Matrix2DDouble& htn2GL, Matrix2DDouble& hte2GL, Matrix2DFloat& dxu1GL,
  Matrix2DFloat& dyu1GL, Matrix2DFloat& tarea1GL, int imt1GL, int jmt1GL,
  int imt2GL, int jmt2GL)
{
  // create grid-related quantities needed for future calculations.
  // dxu1GL, dyu, and tarea are calculated.

  // allocate memory (local variables)
  Matrix2DFloat dxt1GL(imt1GL,jmt1GL);
  Matrix2DFloat dyt1GL(imt1GL,jmt1GL);

  // dxu1GL, dyu, dxt, dyt will have a ghost layer of 1
  for(int j=1; j<jmt2GL-1; j++)
    {
    for(int i=1; i<imt2GL-1; i++)
      {
      dxu1GL(i-1,j-1) = 0.5 * (htn2GL(i,j) + htn2GL(i+1,j));
      dyu1GL(i-1,j-1) = 0.5 * (hte2GL(i,j) + hte2GL(i,j+1));
      dxt1GL(i-1,j-1) = 0.5 * (htn2GL(i,j) + htn2GL(i,j-1));
      dyt1GL(i-1,j-1) = 0.5 * (hte2GL(i,j) + hte2GL(i-1,j));
      }
    }

  // tarea has ghost cell layer of 1
  for(int i=0; i<imt1GL*jmt1GL; i++)
    {
    tarea1GL(i) = dxt1GL(i) * dyt1GL(i);
    }
  // htn2GL is correct. dxt1GL is wrong every 2400 counts but that shouldn't matter because that's always the lowest latitude which doesn't get used anyways

  return 1;
}

//-----------------------------------------------------------------------------
void vtkAbstractPOPReader::sw_4pt(Matrix2DFloat& xout1GL, float factor, Matrix2DDouble& x1GL,
                                  int imt1GL, int jmt1GL)
{
  // perform averaging over a neighborhood
  int i, j;

  for(j=1; j<jmt1GL; j++)
    {
    for(i=1; i<imt1GL; i++)
      {
      xout1GL(i,j) = factor * x1GL(i,j)   +
        factor * x1GL(i,j-1) +
        factor * x1GL(i-1,j) +
        factor * x1GL(i-1,j-1);
      }
    }

  j = 0;
  for(i=1; i<imt1GL; i++)
    {
    xout1GL(i,j) = xout1GL(i,j+1);
    }

  i = 0;
  for(j=1; j<jmt1GL; j++)
    {
    xout1GL(i,j) = factor * x1GL(i,j) +
      factor * x1GL(i,j-1) +
      factor * x1GL(imt1GL-1,j) +
      factor * x1GL(imt1GL-1,j-1);
    }

  xout1GL(0, 0) = xout1GL(1,1);
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::GetLatitudeIndex(Matrix1DFloat& lat_mht, float latitude)
{
  for(unsigned j=0; j<lat_mht.GetSize();j++)
    {
    if(latitude < lat_mht(j))
      {
      return static_cast<int>(j);
      }
    }
  // couldn't find one
  return -1;
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::cshift(int i, int offset, int size)
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
void vtkAbstractPOPReader::GetMOCSize(POPInputInformation* popinfo, int* ny_mht, int* z)
{
  // calculates the final size of the moc array.
  // note that when moc is usually
  // shown, the y (latitude) is shown along the horizontal axis, and the z
  // (depth) is usually shown along the vertical axis

  float temp = ((popinfo->ynorth_mht - popinfo->ysouth_mht) / popinfo->dy_mht);
  *ny_mht = floor(temp + 0.5) + 1;
  *z = popinfo->global_km;
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::LoadData(
  POPInputInformation* popinfo, int* ext3D1GL, int* ext3D2GL, int imt1GL,
  int jmt1GL, int km, Matrix2DDouble& uLat1GL, Matrix2DDouble& uLong1GL,
  Matrix2DDouble& htn2GL, Matrix2DDouble& hte2GL,
  Matrix1DFloat& dz, Matrix2DInt& global_kmt1GL,
  Matrix2DInt& atl_kmt1GL, Matrix2DInt& indopac_kmt1GL,
  Matrix3DFloat& u1GL, Matrix3DFloat& v1GL, int imt2GL, int jmt2GL)
{
  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();

  vtkMultiProcessController* controller =
    vtkMultiProcessController::GetGlobalController();
  // load data from disk

  // uLat
  int offset = popinfo->global_imt * popinfo->global_jmt * 0;
  this->LoadDataBlock2DDouble(popinfo, popinfo->grid_file, ext3D1GL, offset,
                              imt1GL, jmt1GL, uLat1GL);

  // uLong
  offset = popinfo->global_imt * popinfo->global_jmt * 1;
  this->LoadDataBlock2DDouble(popinfo, popinfo->grid_file, ext3D1GL, offset,
                              imt1GL, jmt1GL, uLong1GL);

  // htn
  offset = popinfo->global_imt * popinfo->global_jmt * 2;
  this->LoadDataBlock2DDouble2(popinfo, popinfo->grid_file, ext3D2GL, offset,
                               imt2GL, jmt2GL, htn2GL);

  // hte
  offset = popinfo->global_imt * popinfo->global_jmt * 3;
  this->LoadDataBlock2DDouble2(popinfo, popinfo->grid_file, ext3D2GL, offset,
                               imt2GL, jmt2GL, hte2GL);

  // read depths file. read the first number in each row.
  // process 0 reads it and broadcasts the values.
  if(controller->GetLocalProcessId() == 0)
    {
    FILE* f = fopen(popinfo->in_depths.c_str(), "r");
    if(f == NULL)
      {
      vtkErrorMacro("Error in opening in_depths file: " << popinfo->in_depths.c_str());
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
  this->LoadDataBlock2DInt(popinfo, popinfo->kmt_global_file, ext3D1GL, imt1GL, jmt1GL, global_kmt1GL);
  this->LoadDataBlock2DInt(popinfo, popinfo->kmt_atl_file, ext3D1GL, imt1GL, jmt1GL, atl_kmt1GL);
  this->LoadDataBlock2DInt(popinfo, popinfo->kmt_indopac_file, ext3D1GL, imt1GL, jmt1GL,
                           indopac_kmt1GL);

  // read velocity files
  this->LoadDataBlock3DFloat(popinfo, popinfo->u_file, ext3D1GL, imt1GL, jmt1GL, km, u1GL);
  this->LoadDataBlock3DFloat(popinfo, popinfo->v_file, ext3D1GL, imt1GL, jmt1GL, km, v1GL);

  if(popinfo->byteswap)
    {
    // byteswap all arrays
    uLat1GL.ByteSwap();
    uLong1GL.ByteSwap();
    htn2GL.ByteSwap();
    hte2GL.ByteSwap();
    global_kmt1GL.ByteSwap();
    atl_kmt1GL.ByteSwap();
    indopac_kmt1GL.ByteSwap();
    u1GL.ByteSwap();
    v1GL.ByteSwap();
    }

  timer->StopTimer();
  cerr << controller->GetLocalProcessId() << " has LOADDATA time of " << timer->GetElapsedTime() << endl;

  return 1;
}

//-----------------------------------------------------------------------------
int vtkAbstractPOPReader::LoadDataBlock2DDouble(POPInputInformation* popinfo, std::string filename,
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
    vtkGenericWarningMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = (offset +
                        (y_start_file + y) * popinfo->global_imt +
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
                          (y_start_file + y) * popinfo->global_imt +
                          popinfo->global_imt-1) * sizeof(double);
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
                          (y_start_file + y) * popinfo->global_imt +
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
                        (popinfo->global_jmt-1)*popinfo->global_imt +
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
                        (popinfo->global_jmt-1)*popinfo->global_imt +
                        popinfo->global_imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(double), 1, f);
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (offset +
                        0 +
                        popinfo->global_imt-1) * sizeof(double);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(double), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = (offset +
                        (popinfo->global_jmt-1)*popinfo->global_imt +
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
int vtkAbstractPOPReader::LoadDataBlock2DDouble2(
  POPInputInformation* popinfo, std::string filename, int* ext3DArbitrary,
  int offset, int imt, int jmt, Matrix2DDouble& data)
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
    vtkGenericWarningMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = (offset +
                        (y_start_file + y) * popinfo->global_imt +
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
                            (y_start_file + y) * popinfo->global_imt +
                            popinfo->global_imt-1-x) * sizeof(double);
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
                            (y_start_file + y) * popinfo->global_imt +
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
                          (popinfo->global_jmt-1-y)*popinfo->global_imt +
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
                          (popinfo->global_jmt-1-y)*popinfo->global_imt +
                          popinfo->global_imt+ext3DArbitrary[0]) * sizeof(double);
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
                          popinfo->global_imt+ext3DArbitrary[0]) * sizeof(double);
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
                          (popinfo->global_jmt-1-y)*popinfo->global_imt +
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
int vtkAbstractPOPReader::LoadDataBlock2DInt(POPInputInformation* popinfo, std::string filename, int* ext3DArbitrary,
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
    vtkGenericWarningMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int y=0; y<col_length; y++)
    {
    // get to correct place in file
    int total_offset = ((y_start_file + y) * popinfo->global_imt +
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
      int total_offset = ((y_start_file + y) * popinfo->global_imt +
                          popinfo->global_imt-1) * sizeof(int);
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
      int total_offset = ((y_start_file + y) * popinfo->global_imt +
                          0) * sizeof(int);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(imt-1, y_start+y), sizeof(int), 1, f);
      }
    }
  if(ext3DArbitrary[2] == -1)
    {
    // bottom ghost cells wrap around
    // copy values of the row over
    int total_offset = ((popinfo->global_jmt-1)*popinfo->global_imt +
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
    int total_offset = ((popinfo->global_jmt-1)*popinfo->global_imt +
                        popinfo->global_imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, 0), sizeof(int), 1, f);
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    int total_offset = (0 +
                        popinfo->global_imt-1) * sizeof(int);
    fseek(f, total_offset, SEEK_SET);
    fread(&data(0, jmt-1), sizeof(int), 1, f);
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    int total_offset = ((popinfo->global_jmt-1)*popinfo->global_imt +
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
int vtkAbstractPOPReader::LoadDataBlock3DFloat(POPInputInformation* popinfo, std::string filename,
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
    vtkGenericWarningMacro("Error in opening file: " << filename.c_str());
    return 0;
    }

  // read data one row at a time
  for(int z=0; z<km; z++)
    {
    for(int y=0; y<col_length; y++)
      {
      // get to correct place in file
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                          (y_start_file + y) * popinfo->global_imt +
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
        int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                            (y_start_file + y) * popinfo->global_imt +
                            popinfo->global_imt-1) * sizeof(float);
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
        int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                            (y_start_file + y) * popinfo->global_imt +
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
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                          (popinfo->global_jmt-1)*popinfo->global_imt +
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
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
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
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                          (popinfo->global_jmt-1)*popinfo->global_imt +
                          popinfo->global_imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, 0, z), sizeof(float), 1, f);
      }
    }
  if(ext3DArbitrary[0] == -1 && ext3DArbitrary[3] == -1)
    {
    // upper left corner needs lower right value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                          0 +
                          popinfo->global_imt-1) * sizeof(float);
      fseek(f, total_offset, SEEK_SET);
      fread(&data(0, jmt-1, z), sizeof(float), 1, f);
      }
    }
  if(ext3DArbitrary[1] == -1 && ext3DArbitrary[2] == -1)
    {
    // lower right corner needs upper left value
    for(int z=0; z<km; z++)
      {
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
                          (popinfo->global_jmt-1)*popinfo->global_imt +
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
      int total_offset = (z * popinfo->global_imt * popinfo->global_jmt +
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
int vtkAbstractPOPReader::compare_latitude(const void* x, const void* y)
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
void vtkAbstractPOPReader::FindSouthern(int imt1GL, int jmt1GL, int* ext3D1GL, int* ext3D,
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
  int tempJ = (found == true ? *localJIndexMin1GL + ext3D1GL[2] : VTK_INT_MAX);
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
  controller->AllReduce(&my_southern_lat, southern_lat, 1,
                        vtkCommunicator::MIN_OP);
}

//-----------------------------------------------------------------------------
void vtkAbstractPOPReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;
}
