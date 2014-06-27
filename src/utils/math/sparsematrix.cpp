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

#if 0

#include "sparsematrix.h"
#include "fastarray.h"
#include <stdincludes.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "errors.h"

namespace Math {


SparseMatrix::SparseMatrix()
  :row_offsets(NULL), col_indices(NULL), val_array(NULL),
   m(0), n(0), num_entries(0)
{
}

SparseMatrix::~SparseMatrix()
{
	cleanup();
}

void SparseMatrix::init(int _m, int _n, int _num_entries)
{
  Assert(_m > 0 && _n > 0);
  SafeArrayDelete(row_offsets);
  row_offsets = new int[_m+1];
  row_offsets[m] = _num_entries;

  SafeArrayDelete(col_indices);
  SafeArrayDelete(val_array);
  col_indices = new int[_num_entries];
  val_array = new Real [_num_entries];

  m = _m;
  n = _n;
  num_entries = _num_entries;
}

void SparseMatrix::resize(int _m, int _n, int _num_entries)
{
  if(_m != m || num_entries != _num_entries)
    init(_m,_n,_num_entries);
}

void SparseMatrix::cleanup()
{
  SafeArrayDelete(row_offsets);
  SafeArrayDelete(col_indices);
  SafeArrayDelete(val_array);
  m = n = num_entries = 0;
}

void SparseMatrix::operator = (const SparseMatrix& m)
{
  copy(m);
}

Real SparseMatrix::operator () (int i, int j) const
{
  int* ind;
  Real* f;
  int* min, *max;
  f = rowValues(i);
  min = rowIndices(i);
  max = rowIndices(i+1);
  //binary search
  ind = std::lower_bound(min,max,j);
  if((ind!=max) && (*ind==j)) return f[ind-min];
  return Zero;
}




void SparseMatrix::setZero()
{
  array_zero(val_array, num_entries);
}

void SparseMatrix::setMatrix(const Matrix& mat,Real zeroTol)
{
  int i,j;
  cleanup();
  m=mat.m;
  n=mat.n;
  row_offsets=new int[m+1];
  num_entries=0;
  for(i=0;i<mat.m;i++) {
	  row_offsets[i]=num_entries;
	  for(j=0;j<mat.n;j++)
		  if(!FuzzyZero(mat(i,j),zeroTol)) num_entries++;
  }
  row_offsets[i]=num_entries;
  col_indices = new int[num_entries];
  val_array = new Real [num_entries];

  //everything's set up now, copy the entries
  num_entries=0;
  for(i=0;i<mat.m;i++) {
	  for(j=0;j<mat.n;j++)
		  if(!FuzzyZero(mat(i,j),zeroTol)) {
			  col_indices[num_entries]=j;
			  val_array[num_entries]=mat(i,j);
			  num_entries++;
		  }
  }
}

void SparseMatrix::getMatrix(Matrix& mat) const
{
  int* j;
  int* min, *max;
  Real* f;

  mat.resize(m, n, Zero);
  for(int i=0; i<m; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    for(j=min; j<max; j++, f++) {
      mat(i,*j) = *f;
    }
  }
}

void SparseMatrix::makeSimilar(const SparseMatrix& mat)
{
  int i;
  init(mat.m, mat.n, mat.num_entries);
  for(i=0; i<=m; i++) {
    row_offsets[i] = mat.row_offsets[i];
  }
  for(i=0; i<num_entries; i++) {
    col_indices[i] = mat.col_indices[i];
  }
}

void SparseMatrix::copy(const SparseMatrix& m)
{
  makeSimilar(m);
  array_equal(val_array, m.val_array, num_entries);
}

void SparseMatrix::scale(Real s)
{
  array_scale(val_array, s, num_entries);
}

void SparseMatrix::scaleCols(const Vector& s)
{
  if(s.n != n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int i;
  int* j;
  int* min, *max;
  Real* f;
  for(i=0; i<m; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    for(j=min; j<max; j++, f++) {
      *f *= s[*j];
    }
  }
}



void SparseMatrix::scaleRows(const Vector& s)
{
  if(s.n != m)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int i;
  //int* j;
  //int* min, max;
  Real* f;
  Real* fmax;
  for(i=0; i<m; i++) {
    //min = col_indices + row_offsets[i];
    //max = col_indices + row_offsets[i+1];
    f = rowValues(i);
    fmax = rowValues(i+1);
    //for(j=min; j<max; j++, f++) {
    for(; f<fmax; f++) {
      *f *= s(i);
    }
  }
}

void SparseMatrix::mul(const SparseMatrix& m, Real s)
{
  Assert(num_entries == m.num_entries);
  array_scale(val_array, m.val_array, s, num_entries);
}

Real SparseMatrix::dotCol(int j_row, const Vector& v) const
{
  if(!isValidCol(j_row))
	Error(MatrixError_InvalidCol);
  if(n != v.n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int* j;
  int* min, *max;
  Real* f;
  Real sum;

  sum = 0.0;
  for(int i=0; i<m; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    j = std::lower_bound(min,max,j_row);
    if((j!=max) && (*j==j_row))
      sum += v(i)*(f[j-min]);
  }
  return sum;
}

Real SparseMatrix::dotRow(int i, const Vector& v) const
{
  if(!isValidRow(i))
	Error(MatrixError_InvalidRow);
  if(n != v.n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int* j;
  int* min, *max;
  Real* f;
  Real sum;
  Assert(i >= 0 && i < m);
  min = rowIndices(i);
  max = rowIndices(i+1);
  f = rowValues(i);
  sum = 0.0;
  for(j=min; j<max; j++, f++) {
    sum += (*f) * v(*j);
  }
  return sum;
}

Real SparseMatrix::dotSymmL(int ind, const Vector& v) const
{
  if(!isSquare())
  {
	  Error(MatrixError_IncompatibleDimensions);
  }
  if(m != v.n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }
  if(!isValidRow(ind))
	  Error(MatrixError_InvalidRow);

  int* j;
  int* min, *max;
  Real* f;
  Real sum;
  min = rowIndices(ind);
  max = rowIndices(ind+1);
  f = rowValues(ind);
  //the row part (lower)
  sum = 0.0;
  for(j=min; j<max; j++, f++) {
    Assert(*j > ind);
    sum += (*f) * v(*j);
  }
  //the column part (upper)
  for(int i=0; i<ind; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    j=std::lower_bound(min,max,ind);
    if((j!=max) && (*j==ind))
      sum += v(i)*(f[j-min]);
  }
  return sum;
}

void SparseMatrix::mul(const Vector& w, Vector& v) const
{
  if(w.n != n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }
  if(v.n == 0)
  {
	  v.resize(m);
  }
  else
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int* j;
  int* min, *max;
  Real* f;

  v.set(Zero);
  for(int i=0; i<m; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    for(j=min; j<max; j++, f++) {
      v(i) += *f * w(*j);
    }
  }
} 



void SparseMatrix::mul(const Matrix& w, Matrix& v) const
{
	/*
	QUESTION: this is a lot faster with the matrices being lists of column vectors -- do we want to do it this way?
	//v = M*w
	//so v.m=M.m, v.n=w.n, M.n=w.m
  if(w.n != v.n)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }
  if(w.m != n || v.m != m)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int* j;
  int k;
  int* min, *max;
  Real* f;

  v.setZero();
  for(int i=0; i<m; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    for(k=0; k<w.m; k++) {
      f=rowValues(i);
      for(j=min; j<max; j++, f++) {
        v(k,i) += *f * w(k,*j);
      }
    }
  }
  */
	printf("sparse mul vector not done yet\n");
	Abort();
}




void SparseMatrix::mulTranspose(const Vector& v, Vector& w) const
{
  if(v.n != m)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }
  if(w.n == 0)
  {
	  w.resize(m);
  }
  else if(w.n != m)
  {
	  Error(MatrixError_ArgIncompatibleDimensions);
  }

  int* j;
  int* min, *max;
  Real* f;

  w.set(Zero);
  for(int i=0; i<n; i++) {
    min = rowIndices(i);
    max = rowIndices(i+1);
    f = rowValues(i);
    for(j=min; j<max; j++, f++) {
      w[*j] += *f * v[i];
    }
  }
} 

bool SparseMatrix::isValid() const
{
  bool result = true;
  if(row_offsets[0] != 0) {
    printf("row offsets don't start at 0\n");
    result = false;
  }
  int last = 0;
  int i;
  for(i=1; i<m; i++) {
    if(row_offsets[i] < last) {
      printf("row offset index %d is not in ascending order\n", i);
      result = false;
    }
    last = row_offsets[i];
  }
  if(row_offsets[m] != num_entries) {
    printf("columns don't end at num_entries\n");
    result = false;
  }
  for(i=0; i<num_entries; i++) {
    if(!isValidCol(col_indices[i])) {
      printf("invalid column entry %d in index position %d\n", col_indices[i], i);
      result = false;
    }
    //if(fabs(val_array[i]) < 2e-6) {
    if(fabs(val_array[i]) == 0.0) {
      printf("zero-value %f in index position %d\n", val_array[i], i);
      result = false;
    }
  }
  return result;
}

void SparseMatrix::self_test()
{
  self_test(2,2,4);
  self_test(10,10, 50);
  self_test(5,10, 30);
}

int rand(int i) {
  return ::rand()%i;
}

float frand(float min, float max) {
  float u = float(::rand())/RAND_MAX;
  return min + u*(max-min);
}

#define EPSILON 2e-3
#include <math.h>

void SparseMatrix::self_test(int m, int n, int nnz) {
  Matrix matdense (m,n,Zero);
  Matrix matdense2 (m,n);
  Vector x (n);
  Vector y (m), y2 (m);

  SparseMatrixBuilder matbuild;
  matbuild.initialize(m,n);

  int i,p,q;
  float v;
  for(i=0; i<nnz; i++) {
    p = rand(m);
    q = rand(n);
    v = frand(-100, 100);
    matdense(p,q) += v;
    matbuild.addEntry(p,q,v);
  }
  for(i=0; i<n; i++)
    x[i] = frand(-100, 100);

  SparseMatrix mat;
  matbuild.buildMatrix(mat);
  mat.getMatrix(matdense2);
  for(i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      if(fabs(matdense(i,j) - matdense2(i,j)) > EPSILON) {
	printf("Error testing sparse matrix %d %d\n", m,n);
	printf("Matrices differ in %d,%d\n", i,j);
	std::cout << matdense << std::endl;
	std::cout << matdense2 << std::endl;
	exit(1);
      }
    }
  }

  matdense.mul(x,y);
  mat.mul(x, y2);
  for(i=0; i<m; i++) {
    if(fabs(y2[i]- y[i]) > EPSILON) {
      printf("Error testing sparse matrix %d %d\n", m,n);
      printf("Matrix multiply differs, %f vs %f!\n",y[i],y2[i]);
      exit(1);
    }
  }
}






SparseMatrixBuilder::SparseMatrixBuilder()
  :rows(NULL), m(0), n(0)
{}

SparseMatrixBuilder::~SparseMatrixBuilder()
{
  SafeArrayDelete(rows);
}

void SparseMatrixBuilder::initialize(int _m, int _n)
{
  SafeArrayDelete(rows);
  m = _m;
  n = _n;
  rows = new RowT[m];
}

void SparseMatrixBuilder::resize(int _m, int _n)
{
  n = _n;
  if(_m != m)
    initialize(_m,_n);
}

void SparseMatrixBuilder::setZero()
{
  for(int i=0; i<m; i++) {
    rows[i].clear();
  }
}

void SparseMatrixBuilder::set(const SparseMatrix& mat)
{
  resize(mat.m, mat.n);
  setZero();
  for(int i=0; i<m; i++) {
    int* min = mat.rowIndices(i);
    int* max = mat.rowIndices(i+1);
    Real* f = mat.rowValues(i);
    for(int* j=min; j<max; j++, f++)
      insertEntry(i,*j,*f);
  }
}

void SparseMatrixBuilder::zeroEntry(int i, int j)
{
  for(RowIterator e=rows[i].begin();e!=rows[i].end();e++) {
    if(e->col == j) {
      rows[i].erase(e);
      return;
    }
  }
  printf("Error- entry is already zero\n");
  AssertNotReached();
}

void SparseMatrixBuilder::insertEntry(int i, int j, Real val)
{
  if(val != 0.0) {
    matrix_entry_t e;
    e.col=j;
    e.val=val;
    rows[i].push_back(e);
  }
}


void SparseMatrixBuilder::addEntry(int i, int j, Real val)
{
  if(val == 0.0) {
    return;
  }

  Real* f = operator() (i,j);
  if(!f)
    insertEntry(i,j,val);
  else {
    *f += val;
    if(*f == 0.0)
      zeroEntry(i,j);
  }
}

void SparseMatrixBuilder::mul(const SparseMatrix& a, const SparseMatrix& b)
{
  Assert(a.n == b.m);
  resize(a.m, b.n);
  /*
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      Real sum = 0.0;
      for(int k=0; k<a.n; k++) {
	sum += a(i,k)*b(k,j);
      }
      insertEntry(i,j,sum);
    }
  }
  */

  setZero();
  for(int i=0; i<m; i++) {
    int *amin, *amax;
    Real* e;
    amin = a.rowIndices(i);
    amax = a.rowIndices(i+1);
    e = a.rowValues(i);
    for(int* kp=amin; kp<amax; kp++, e++) {
      int k = *kp;
      int* bmin, *bmax;
      Real* f;
      bmin = b.rowIndices(k);
      bmax = b.rowIndices(k+1);
      f = b.rowValues(i);
      for(int* jp=bmin; jp<bmax; jp++, f++) {
        //entry goes into (i,j)
        int j = *jp;
        addEntry(i,j,(*e)*(*f));
      }
    }
  }
}

void SparseMatrixBuilder::mul(const SparseMatrixBuilder& a, const SparseMatrixBuilder& b)
{
  Assert(a.n == b.m);
  resize(a.m, b.n);
  /*
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      Real sum = 0.0;
      for(int k=0; k<a.n; k++) {
	sum += a(i,k)*b(k,j);
      }
      insertEntry(i,j,sum);
    }
  }
  */

  setZero();
  for(int i=0; i<m; i++) {
    for(RowIterator e=a.rows[i].begin();e!=a.rows[i].end();e++) {
      int k = e->col;
      for(RowIterator f=b.rows[k].begin();f!=b.rows[k].end();f++) {
        //entry goes into (i,j)
        int j = f->col;
        addEntry(i,j,e->val*f->val);
      }
    }
  }
}

Real* SparseMatrixBuilder::operator () (int i, int j) {
  for(RowIterator e=rows[i].begin();e!=rows[i].end();e++) {
    if(e->col == j) {
      return &e->val;
    }
  }
  return NULL;
}

void SparseMatrixBuilder::buildMatrix(SparseMatrix& mat) const
{
  int i,j;
  int num_entries = 0;
  RowIterator e;

  //test for duplicates
  int* cols = new int [n];
  for(i=0; i<m; i++) {
    for(j=0; j<n; j++)
      cols[j] = 0;
    for(e=rows[i].begin();e!=rows[i].end();e++) {
      if(cols[e->col]) {
        printf("Duplicated row %d %d\n", i, e->col);
	Abort();
      }
      cols[e->col] = 1;
    }
  }
  delete [] cols;

  //count num nonzeros
  for(i=0; i<m; i++) {
    num_entries += rows[i].size();
  }
  mat.init(m,n,num_entries);
  j = 0;
  for(i=0; i<m; i++) {
    mat.row_offsets[i] = j;
    for(e=rows[i].begin();e!=rows[i].end();e++) {
      mat.col_indices[j] = e->col;
      mat.val_array[j] = e->val;
      j++;
    }
  }
  mat.row_offsets[i] = j;
  Assert(j == num_entries);
}




std::istream& operator >> (std::istream& in, SparseMatrix& mat)
{
  int m, n, num_entries;
  in >> m >> n >> num_entries;
  if(!in) return in;
  mat.init(m,n,num_entries);

  int i;
  for(i=0; i<m; i++)
    in >> mat.row_offsets[i];
  mat.row_offsets[n] = num_entries;
  for(i=0; i<num_entries; i++)
    in >> mat.col_indices[i];
  for(i=0; i<num_entries; i++)
    in >> mat.val_array[i];
  return in;
}

std::ostream& operator << (std::ostream& out, const SparseMatrix& mat)
{
  out << mat.m << " " << mat.n << " " << mat.num_entries << std::endl;
  int i;
  for(i=0; i<mat.m; i++)
    out << mat.row_offsets[i] << " ";
  out << std::endl;
  for(i=0; i<mat.num_entries; i++)
    out << mat.col_indices[i] << " ";
  out << std::endl;
  for(i=0; i<mat.num_entries; i++)
    out << mat.val_array[i] << " ";
  out << std::endl;
  return out;
}


void SparseMatrix::print(std::ostream& out) const
{
  for(int i=0; i<m; i++) {
    out << "Row " << i << ": ";
    int j = row_offsets[i];
    int max = row_offsets[i+1];
    for(;j<max;j++) 
      out << "(" << col_indices[j] << "," << val_array[j] << ") ";
    out << std::endl;
  }
}

void SparseMatrixBuilder::print(std::ostream& out) const
{
  for(int i=0; i<m; i++) {
    out << "Row " << i << ": ";
    for(RowIterator e=rows[i].begin();e!=rows[i].end();e++) {
      out << "(" << e->col << "," << e->val << ") ";
    }
    out << std::endl;
  }
}

} //namespace Math

#endif // 0
