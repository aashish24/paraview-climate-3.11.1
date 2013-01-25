#ifndef __MATRIX_HELPER__
#define __MATRIX_HELPER__

// matrix classes used in the POP readers

//-----------------------------------------------------
//  1D float array
//-----------------------------------------------------

class Matrix1DFloat
{
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

 private:
  unsigned XDim;
  float* Data;
};

//-----------------------------------------------------
//  2D int array
//-----------------------------------------------------

class Matrix2DInt
{
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

//-----------------------------------------------------
//  2D float array
//-----------------------------------------------------

class Matrix2DFloat
{
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

//-----------------------------------------------------
//  2D double array
//-----------------------------------------------------

class Matrix2DDouble
{
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

//-----------------------------------------------------
//  3D float array
//-----------------------------------------------------

class Matrix3DFloat
{
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

#endif
