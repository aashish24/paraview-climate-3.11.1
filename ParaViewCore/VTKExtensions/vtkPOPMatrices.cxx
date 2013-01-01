#include "vtkPOPMatrices.h"
#include <iostream>

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
} // end anonymous namespace

Matrix1DFloat::Matrix1DFloat()
{
  this->XDim = 0;
  this->Data = NULL;
}

Matrix1DFloat::Matrix1DFloat(unsigned int xDim)
{
  this->XDim = xDim;
  this->Data = new float[xDim];
}

void Matrix1DFloat::Allocate(unsigned int xDim)
{
  if(this->XDim != xDim)
    {
    this->Clear();
    this->XDim = xDim;
    this->Data = new float[xDim];
    }
}

Matrix1DFloat::~Matrix1DFloat()
{
  this->Clear();
}

float* Matrix1DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

void Matrix1DFloat::Clear()
{
  if(this->Data != NULL)
    {
    delete[] this->Data;
    this->Data = NULL;
    }
  this->XDim = 0;
}

void Matrix1DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

float& Matrix1DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

float Matrix1DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D int array methods
//-----------------------------------------------------
Matrix2DInt::Matrix2DInt()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

Matrix2DInt::Matrix2DInt(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new int[xDim * yDim];
}

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

Matrix2DInt::~Matrix2DInt()
{
  this->Clear();
}

int* Matrix2DInt::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

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

void Matrix2DInt::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

int& Matrix2DInt::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

int Matrix2DInt::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

int& Matrix2DInt::operator() (unsigned int i)
{
  return this->Data[i];
}

int Matrix2DInt::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D float array methods
//-----------------------------------------------------
Matrix2DFloat::Matrix2DFloat()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

Matrix2DFloat::Matrix2DFloat(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new float[xDim * yDim];
}

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

Matrix2DFloat::~Matrix2DFloat()
{
  this->Clear();
}

float* Matrix2DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

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

void Matrix2DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

float& Matrix2DFloat::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

float Matrix2DFloat::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

float& Matrix2DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

float Matrix2DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  2D double array methods
//-----------------------------------------------------
Matrix2DDouble::Matrix2DDouble()
{
  this->XDim = 0;
  this->YDim = 0;

  this->Data = NULL;
}

Matrix2DDouble::Matrix2DDouble(unsigned int xDim, unsigned int yDim)
{
  this->XDim = xDim;
  this->YDim = yDim;

  this->Data = new double[xDim * yDim];
}

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

Matrix2DDouble::~Matrix2DDouble()
{
  this->Clear();
}

double* Matrix2DDouble::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

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

void Matrix2DDouble::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim; i++)
    {
    Swap8((char*)&this->Data[i]);
    }
}

double& Matrix2DDouble::operator() (unsigned int i, unsigned int j)
{
  return this->Data[this->XDim*j + i];
}

double Matrix2DDouble::operator() (unsigned int i, unsigned int j) const
{
  return this->Data[this->XDim*j + i];
}

double& Matrix2DDouble::operator() (unsigned int i)
{
  return this->Data[i];
}

double Matrix2DDouble::operator() (unsigned int i) const
{
  return this->Data[i];
}

//-----------------------------------------------------
//  3D float array methods
//-----------------------------------------------------
Matrix3DFloat::Matrix3DFloat()
{
  this->XDim = 0;
  this->YDim = 0;
  this->ZDim = 0;

  this->Data = NULL;
}

Matrix3DFloat::Matrix3DFloat(unsigned int xDim, unsigned int yDim, unsigned int zDim)
{
  this->XDim = xDim;
  this->YDim = yDim;
  this->ZDim = zDim;

  this->Data = new float[xDim * yDim * zDim];
}

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

Matrix3DFloat::~Matrix3DFloat()
{
  this->Clear();
}

float* Matrix3DFloat::GetData()
{
  // return the address of the actual data stored
  return this->Data;
}

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

void Matrix3DFloat::ByteSwap()
{
  // ByteSwap data in place
  for(unsigned i=0; i<this->XDim*this->YDim*this->ZDim; i++)
    {
    Swap4((char*)&this->Data[i]);
    }
}

float& Matrix3DFloat::operator() (unsigned int i, unsigned int j, unsigned int k)
{
  return this->Data[this->XDim*this->YDim*k + this->XDim*j + i];
}

float Matrix3DFloat::operator() (unsigned int i, unsigned int j, unsigned int k) const
{
  return this->Data[this->XDim*this->YDim*k + this->XDim*j + i];
}

float& Matrix3DFloat::operator() (unsigned int i)
{
  return this->Data[i];
}

inline
float Matrix3DFloat::operator() (unsigned int i) const
{
  return this->Data[i];
}
