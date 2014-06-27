/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

#include "MatrixTemplate.h"
//#include "MatrixPrinter.h"
#include "vector.h"

namespace Math {
  typedef MatrixTemplate<Real> Matrix;
}

#if 0
#include "vector.h"

namespace Math {

class Matrix
{
public:
	Matrix();
	Matrix(const Matrix&);
	Matrix(int m, int n);
	Matrix(int m, int n, Real initval);
	Matrix(int m, int n, const Real* vals);
	Matrix(int m, int n, const Real** vals);
	Matrix(int m, int n, const Vector* rows);
	~Matrix();

	void clear();
	void resize(int m, int n);
	void resize(int m, int n, Real initval);
	inline int numRows() const;
	inline int numCols() const;

	const Matrix& operator = (const Matrix&);
	bool operator == (const Matrix&) const;
	inline bool operator != (const Matrix&) const;
	inline operator Real** ();
	inline operator Real* ();
	inline const Real& operator() (int,int) const;
	inline Real& operator() (int,int);
	void operator += (const Matrix&);
	void operator -= (const Matrix&);
	void operator *= (const Matrix&);
	void operator *= (Real);
	void operator /= (Real);

	void add(const Matrix&, const Matrix&);
	void sub(const Matrix&, const Matrix&);
	void mul(const Matrix&, const Matrix&);
	void mulTransposeA(const Matrix& a, const Matrix& b);
	void mulTransposeB(const Matrix& a, const Matrix& b);
	void mul(const Vector&, Vector&) const;
	void mulTranspose(const Vector&, Vector&) const;
	void mul(const Matrix&, Real);
	void div(const Matrix&, Real);
	void madd(const Matrix&, Real);
	void madd(const Vector&, Vector&) const;
	void maddTranspose(const Vector&, Vector&) const;

	void set(const Matrix&);
	void set(Real);
	void set(const Real* vals);
	void set(const Real** vals);
	void setRows(const Vector* rows);
	void setCols(const Vector* cols);
	void setZero();
	void setIdentity();
	void setNegative(const Matrix&);
	void setTranspose(const Matrix&);
	void setInverse(const Matrix&);			//uses LU decomposition
	void setSubMatrix(int i, int j, const Matrix&);
	void setDiagonal(const Vector&);
	void setDiagonal(Real c);

	void inplaceNegative();
	void inplaceScale(Real);
	void inplaceTranspose();
	void inplaceInverse();

	void getSubMatrix(int i, int j, Matrix&) const;
	void getDiagonal(Vector&) const;

	inline bool isEmpty() const;
	inline bool hasDims(int m, int n) const;
	inline bool isValidRow(int) const;
	inline bool isValidCol(int) const;
	inline bool isSquare() const;
	bool isZero() const;
	bool isIdentity() const;
	bool isDiagonal() const;
	bool isSymmetric() const;
	bool isLowerTriangular() const;
	bool isUpperTriangular() const;
	bool isDiagonallyDominant() const;
	inline bool isOrthogonal() const;
	inline bool isInvertible() const;

	Real trace() const;
	Real determinant() const;
	Real diagonalProduct() const;
	Real minElement(int*i=NULL,int*j=NULL) const;
	Real maxElement(int*i=NULL,int*j=NULL) const;
	Real minAbsElement(int*i=NULL,int*j=NULL) const;
	Real maxAbsElement(int*i=NULL,int*j=NULL) const;

	void setRow(int i, Real c);
	void setCol(int j, Real c);
	void setRow(int i, const Vector&);
	void setCol(int j, const Vector&);
	void setRow(int i, const Real*);
	void setCol(int j, const Real*);
	void getRow(int i, Vector&) const;
	void getCol(int j, Vector&) const;
	void addRow(int i,const Vector& v) const;
	void addCol(int j,const Vector& v) const;
	void addRow(int i,const Matrix& m,int im) const;
	void addCol(int j,const Matrix& m,int jm) const;
	void scaleRow(int i,Real c) const;
	void scaleCol(int j,Real c) const;
	void maddRow(int i,const Vector& v,Real c) const;
	void maddCol(int j,const Vector& v,Real c) const;
	void maddRow(int i,const Matrix& m,int im,Real c) const;
	void maddCol(int j,const Matrix& m,int jm,Real c) const;
	Real dotRow(int i,const Vector& v) const;
	Real dotCol(int j,const Vector& v) const;
	Real dotRow(int i,const Matrix& m,int im) const;
	Real dotCol(int j,const Matrix& m,int jm) const;

        void print(std::ostream& out=std::cout,char delim=' ',char bracket='[') const;
	bool load(File&);
	bool save(File&);

	int m,n;
	Real** vals;
};

std::ostream& operator << (std::ostream&, const Matrix&);
std::istream& operator >> (std::istream&, Matrix&);

extern const char* MatrixError_IncompatibleDimensions;
extern const char* MatrixError_ArgIncompatibleDimensions;
extern const char* MatrixError_DestIncompatibleDimensions;
extern const char* MatrixError_SizeZero;
extern const char* MatrixError_NotSquare;
extern const char* MatrixError_NotSymmetric;
extern const char* MatrixError_InvalidRow;
extern const char* MatrixError_InvalidCol;


inline int Matrix::numRows() const { return m; }
inline int Matrix::numCols() const { return n; }
inline Matrix::operator Real** () {	return vals; }
inline Matrix::operator Real* () { return vals[0]; }

inline bool Matrix::operator != (const Matrix& m) const
{
  return !operator == (m);
}

inline const Real& Matrix::operator() (int i,int j) const
{
#ifdef _DEBUG
	if(!isValidRow(i))
		Error(MatrixError_InvalidRow);
	if(!isValidCol(j))
		Error(MatrixError_InvalidCol);
#endif
	return vals[i][j];
}

inline Real& Matrix::operator() (int i,int j)
{
#ifdef _DEBUG
	if(!isValidRow(i))
		Error(MatrixError_InvalidRow);
	if(!isValidCol(j))
		Error(MatrixError_InvalidCol);
#endif
	return vals[i][j];
}

inline bool Matrix::isEmpty() const
{
	return m == 0 && n == 0;
}

inline bool Matrix::hasDims(int _m, int _n) const
{
	return m == _m && n == _n;
}

inline bool Matrix::isValidRow(int i) const
{
	return i >= 0 && i < m;
}

inline bool Matrix::isValidCol(int j) const
{
	return j >= 0 && j < n;
}

inline bool Matrix::isSquare() const
{
	return m == n;
}


inline bool Matrix::isInvertible() const
{
	if(isEmpty())
	{
		Error(MatrixError_SizeZero);
	}
	if(!isSquare())
		return false;

	return !FuzzyZero(determinant());
}

} //namespace Math

namespace std
{
  template<> inline void swap(Math::Matrix& a, Math::Matrix& b)
  {
    Math::Real** vals_tmp = a.vals;
    int m_tmp = a.m;
    int n_tmp = a.n;
    a.vals = b.vals;
    a.m = b.m;
    a.n = b.n;
    b.vals = vals_tmp;
    b.m = m_tmp;
    b.n = n_tmp;
  }
} //namespace std

#endif

#endif
