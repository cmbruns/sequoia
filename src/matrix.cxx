#include "matrix.hxx"
#include "vector.cxx"

// $Id: matrix.cxx,v 1.3 2001/12/12 03:52:06 cdputnam Exp $
// $Header: /usr/data/cvs/sequoia1/matrix.cxx,v 1.3 2001/12/12 03:52:06 cdputnam Exp $
// $Log: matrix.cxx,v $
// Revision 1.3  2001/12/12 03:52:06  cdputnam
// Comment out friend functions: operator<<, operator>>, Matrix<Type> operator*
// Reimplement only operator<< using public fn Matrix<Type>::print(ostream&)
// Now compiles with gcc 2.95 and 2.96.
//
// Revision 1.2  2001/11/28 23:08:17  bruns
// Added cvs header tags
// Removed ^M characters
//

// template <class Type>
// Matrix<Type>::Matrix() : Vector< Vector<Type> >()
// {
// }

// template <class Type>
// Matrix<Type>::Matrix(uint r, uint c, Type val)
// {
//   init_mat(r,c,val);
// }

// copy constructor
// template <class Type>
// Matrix<Type>::Matrix(const Matrix<Type> & M) : Vector< Vector<Type> >(M)
// {
// }

template <class Type>
Matrix<Type>::Matrix(int r, int c, const Type* M)
{
  init_mat(r,c);
  Matrix<Type> & t = *this;
  uint i, j;
  for (i=0; i < m(); ++i)
    for (j=0; j < n(); ++j)
      t[i][j] = M[i*c + j];
}

// template <class Type>
// Matrix<Type>::~Matrix()
// {
// }

template<class Type>
void Matrix<Type>::init_mat(uint r, uint c)
{
  Vector<Type> bit (c);
  init_array(r, bit);
}

template<class Type>
void Matrix<Type>::init_mat(uint r, uint c, Type val)
{
  Vector<Type> bit (c, val);
  init_array(r, bit);
}


template<class Type>
Matrix<Type> Matrix<Type>::T() const // transpose
{
  Matrix<Type> answer(n(), m());
  const Matrix<Type> & M = *this;
  uint i, j;
  for (i=0; i < m(); ++i)
    for (j=0; j < n(); ++j)
      answer[j][i] = M[i][j];
  return answer;
}

template<class Type>
const Matrix<Type> Matrix<Type>::csub(uint i1, uint i2, uint j1, uint j2) const
{
  const Matrix<Type> & a = *this;
  Matrix<Type> answer (i2-i1+1, j2-j1+1);
  uint i, j;
  for (i=0; i < answer.m(); ++i)
    for (j=0; j < answer.n(); ++j)
      answer[i][j] = a[i+i1][j+j1];
  return answer;
}

template<class Type>
subMatrix<Type> Matrix<Type>::sub(uint i1, uint i2, uint j1, uint j2)
{
  subMatrix<Type> answer;
  answer.i_offset = i1;
  answer.j_offset = j1;
  answer.i_dim = i2 - i1 + 1;
  answer.j_dim = j2 - j1 + 1;
  answer.parent = this;
  return answer;
}

// template<class Type>
// matrix_column<Type> Matrix<Type>::col(uint i)
// {
//   matrix_column<Type> answer;
//   answer.column = i;
//   answer.parent = this;
//   return answer;
// }

template<class Type>
const Vector<Type> Matrix<Type>::ccol(uint j) const
{
  const Matrix<Type> & M = *this;
  Vector<Type> answer (m());
  uint i;
  for (i=0; i < m(); ++i)
    answer[i] = M[i][j];
  return answer;
}


template<class Type>
Matrix<Type> Matrix<Type>::operator*(const Matrix<Type> & M2) const
{
  const Matrix<Type> & M1 = *this;
  Matrix<Type> answer (M1.m(), M2.n());
  uint i, j, k;
  uint r = MIN(M1.n(), M2.m());
  for (i=0; i < answer.m(); ++i)
    for (j=0; j < answer.n(); ++j)
      {
	answer[i][j] = 0;
	for (k=0; k < r; ++k)
	  answer[i][j] += M1[i][k] * M2[k][j];
      }
  return answer;
}


template<class Type>
Vector<Type> Matrix<Type>::operator*(const Vector<Type> & v) const
{
  const Matrix<Type> & M = *this;
  Vector<Type> answer (M.m());
  uint r = MIN(M.n(), v.dim());
  uint i, j;
  for (i=0; i < M.m(); ++i)
    {
      answer[i] = 0;
      for (j=0; j < r; ++j)
	answer[i] += M[i][j] * v[j];
    }
  return answer;
}


// outer product
template<class Type>
Matrix<Type> operator*(const Vector<Type> & v1, const row_vector<Type> & v2)
{
  Matrix<Type> answer(v1.dim(), v2.dim());
  uint i, j;
  for (i=0; i < answer.m(); ++i)
    for (j=0; j < answer.n(); ++j)
      answer[i][j] = v1[i] * v2[j];
  return answer;
}


// friend function
//template<class Type>
//Matrix<Type> operator*(double r, const Matrix<Type> & M)
//{
//  return r * ((Vector< Vector<Type> >)M);
//}

template <class Type>
Type Matrix<Type>::mean() const
{
  uint i,j;
  Type sum, answer;
  const Matrix<Type> & M = *this;
  
  sum = 0;
  for (i=0;i < m();++i) 
    for (j=0;j < n();++j)
      sum += M[i][j];
  answer = sum/ ( (Type) (m()*n()) );
  return answer;
}

template <class Type>
void Matrix<Type>::print(ostream & os) const
{
  const Matrix<Type> & M = *this;
  uint i, j;
  for (i=0; i < M.m(); ++i)
    {
      for (j=0; j < M.n(); ++j)
	{
	  os.precision(3);
	  os.width(7);
	  os << M[i][j];
	}
      os << "\n";
    }
}

//template <class Type>
//istream & operator>>(istream & is, Matrix<Type> & M)
//{
//  return is;
//}

template <class Type>
Type Matrix<Type>::std_dev() const
{
  uint i,j;
  Type diff, answer;
  const Matrix & s = *this;
  Type sum = 0;
  Type mean = this->mean();

  for (i=0;i < m();++i) for (j=0;j < n();++j)
  {
    diff = mean - s[i][j];
    sum += (diff * diff);
  }

  answer = (Type) sqrt(sum / ( (Type) (m()*n()) ));
  return answer;
}


template <class Type>
const Vector<Type> Matrix<Type>::diag() const
{
  Vector<Type> answer (MIN(m(),n()));
  uint i;
  for (i=0; i < answer.dim(); ++i)
    answer[i] = this->operator[](i)[i];
  return answer;
}

template <class Type>
subMatrix<Type> & subMatrix<Type>::operator=(const Matrix<Type> & M)
{
  Matrix<Type> & M1 = *parent;
  uint i, j;
  for (i=0; i < i_dim; ++i)
    for (j=0; j < j_dim; ++j)
      M1[i + i_offset][j + j_offset] = M[i][j];
  return *this;
}


template <class Type>
matrix_column<Type> & matrix_column<Type>::operator=(const Vector<Type> & v)
{
  Matrix<Type> & M1 = *parent;
  uint i;
  for (i=0; i < parent->m(); ++i)
    M1[i][column] = v[i];
  return *this;
}

template <class Type>
const Vector<Type> matrix_column<Type>::sub(uint i, uint j) const
{
  Vector<Type> answer(j - i);
  uint k;
  for (k = i; k <= j; ++k)
    answer[k-i] = this->operator[](k);
  return answer;
}

// template class Matrix<int>;
//
// #include "array.cxx"
// #include "vector.cxx"
//
// template class Abstract_array< Vector<int> >;

template<class Type>
ostream& operator<< (ostream& os, const Matrix<Type> & M )
{
  M.print(os);
 return os;
}
