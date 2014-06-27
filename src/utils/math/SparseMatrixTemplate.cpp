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

#include "SparseMatrixTemplate.h"
#include "complex.h"
using namespace Math;
using namespace std;

template <class T>
SparseMatrixTemplate_RM<T>::SparseMatrixTemplate_RM()
  :m(0),n(0)
{}

template <class T>
void SparseMatrixTemplate_RM<T>::initialize(int _m, int _n)
{
  clear();
  resize(_m,_n);
}

template <class T>
void SparseMatrixTemplate_RM<T>::resize(int _m, int _n)
{
  if(m != _m || n != _n) {
    m = _m;
    n = _n;
    rows.resize(m);
    for(size_t i=0;i<rows.size();i++)
      rows[i].n = n;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::clear()
{
  m=n=0;
  rows.clear();
}

template <class T>
void SparseMatrixTemplate_RM<T>::setZero()
{
  for(size_t i=0;i<rows.size();i++)
    rows[i].clear();
}

template <class T>
T& SparseMatrixTemplate_RM<T>::operator () (int i, int j)
{
  Assert(isValidRow(i));
  Assert(isValidCol(j));
  RowIterator entry=rows[i].find(j);
  if(entry != rows[i].end()) return entry->second;
  else {
    entry = rows[i].insert(j,(T)0);
    return entry->second;
  }
}

template <class T>
T* SparseMatrixTemplate_RM<T>::getEntry(int i, int j)
{
  Assert(isValidRow(i));
  Assert(isValidCol(j));
  RowIterator entry=rows[i].find(j);
  if(entry != rows[i].end()) {
    Assert(entry->first == j);
    return &(entry->second);
  }
  else return NULL;
}

template <class T>
const T* SparseMatrixTemplate_RM<T>::getEntry(int i,int j) const
{
  Assert(isValidRow(i));
  Assert(isValidCol(j));
  ConstRowIterator entry=rows[i].find(j);
  if(entry != rows[i].end()) {
    Assert(entry->first == j);
    return &(entry->second);
  }
  else return NULL;
}

template <class T>
void SparseMatrixTemplate_RM<T>::insertEntry(int i, int j, const T& val)
{
  Assert(isValidRow(i));
  Assert(isValidCol(j));
  rows[i].insert(j,val);
}

template <class T>
void SparseMatrixTemplate_RM<T>::eraseEntry(int i, int j)
{
  Assert(isValidRow(i));
  Assert(isValidCol(j));
  bool res = rows[i].erase(j);
  if(!res) {
    cerr<<"Warning, entry "<<i<<","<<j<<" doesn't exist"<<endl;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::copy(const MyT& A)
{
  m=A.m;
  n=A.n;
  rows=A.rows;
}

template <class T>
template <class T2>
void SparseMatrixTemplate_RM<T>::copy(const SparseMatrixTemplate_RM<T2>& A)
{
  initialize(A.m,A.n);
  for(int i=0;i<m;i++) {
    typename SparseMatrixTemplate_RM<T2>::ConstRowIterator it;
    for(it=A.rows[i].begin();it!=A.rows[i].end();it++)
      insertEntry(i,it->first,(T)it->second);
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::set(const MatrixT& A,T zeroTol)
{
  resize(A.m,A.n);
  setZero();
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      if(!FuzzyZero(A(i,j),zeroTol)) rows[i].insert(j,A(i,j));
}

template <class T>
void SparseMatrixTemplate_RM<T>::setTranspose(const MyT& A)
{
  resize(A.n,A.m);
  setZero();
  for(int i=0;i<A.m;i++) {
    for(ConstRowIterator it=A.rows[i].begin();it!=A.rows[i].end();it++)
      insertEntry(it->first,i,it->second);
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::get(MatrixT& A) const
{
  A.resize(m,n,Zero);
  for(int i=0;i<m;i++) {
    for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
      A(i,it->first) = it->second;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::mul(const MyT& a, T s)
{
  copy(a);
  inplaceMul(s);
}

template <class T>
void SparseMatrixTemplate_RM<T>::mul(const VectorT& a,VectorT& x) const
{
  if(x.n == 0) x.resize(m);
  if(x.n != m) {
    FatalError("Destination vector has incorrect dimensions");
  }
  if(a.n != n) {
    FatalError("Source vector has incorrect dimensions");
  }
  for(int i=0;i<m;i++) {
    T sum=0;
    for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
      sum += it->second*a(it->first);
    x(i) = sum;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::madd(const VectorT& a,VectorT& x) const
{
  if(x.n != m) {
    FatalError("Destination vector has incorrect dimensions");
  }
  if(a.n != n) {
    FatalError("Source vector has incorrect dimensions");
  }
  for(int i=0;i<m;i++) {
    T sum=0;
    for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
      sum += it->second*a(it->first);
    x(i) += sum;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::mulTranspose(const VectorT& a,VectorT& x) const
{
  if(x.n == 0) x.resize(n);
  if(x.n != n) {
    FatalError("Destination vector has incorrect dimensions");
  }
  if(a.n != m) {
    FatalError("Source vector has incorrect dimensions");
  }
  x.setZero();
  for(int i=0;i<m;i++) {
    for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
      x(it->first) += it->second*a(i);
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::maddTranspose(const VectorT& a,VectorT& x) const
{
  if(x.n != n) {
    FatalError("Destination vector has incorrect dimensions");
  }
  if(a.n != m) {
    FatalError("Source vector has incorrect dimensions");
  }
  for(int i=0;i<m;i++) {
    for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
      x(it->first) += it->second*a(i);
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::mul(const MatrixT& a,MatrixT& x) const
{
  if(a.m != m) {
    FatalError("A matrix has incorrect # of rows");
  }
  if(x.isEmpty()) x.resize(m,a.n);
  if(m != x.m) {
    FatalError("X matrix has incorrect # of rows");
  }
  if(a.n != x.n) {
    FatalError("X matrix has incorrect # of columns");
  }
  for(int i=0;i<a.n;i++) {
    VectorT ai,xi;
    a.getColRef(i,ai);
    x.getColRef(i,xi);
    mul(ai,xi);
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::mulTranspose(const MatrixT& a,MatrixT& x) const
{
  if(a.m != n) {
    FatalError("A matrix has incorrect # of rows");
  }
  if(x.isEmpty()) x.resize(n,a.n);
  if(n != x.m) {
    FatalError("X matrix has incorrect # of rows");
  }
  if(a.n != x.n) {
    FatalError("X matrix has incorrect # of columns");
  }
  for(int i=0;i<a.n;i++) {
    VectorT ai,xi;
    a.getColRef(i,ai);
    x.getColRef(i,xi);
    mulTranspose(ai,xi);
  }
}

template <class T>
T SparseMatrixTemplate_RM<T>::dotRow(int i, const VectorT& v) const
{
  Assert(isValidRow(i));
  Assert(v.n == n);
  T sum=0;
  for(ConstRowIterator it=rows[i].begin();it!=rows[i].end();it++)
    sum += v(it->first)*it->second;
  return sum;
}
template <class T>
T SparseMatrixTemplate_RM<T>::dotCol(int j, const VectorT& v) const
{
  Assert(isValidCol(j));
  Assert(v.n == m);
  T sum=0;
  for(int i=0;i<m;i++) {
    ConstRowIterator it=rows[i].find(j);
    if(it != rows[i].end()) sum += v(i)*it->second;
  }
  return sum;
}

template <class T>
void SparseMatrixTemplate_RM<T>::inplaceMul(T c)
{
  for(int i=0;i<m;i++) {
    for(RowIterator it=rows[i].begin();it!=rows[i].end();it++)
      it->second *= c;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::inplaceDiv(T c)
{
  for(int i=0;i<m;i++) {
    for(RowIterator it=rows[i].begin();it!=rows[i].end();it++)
      it->second /= c;
  }
}

template <class T>
void SparseMatrixTemplate_RM<T>::inplaceMulRow(int i,T c)
{
  Assert(isValidRow(i));
  for(RowIterator it=rows[i].begin();it!=rows[i].end();it++)
    it->second *= c;
}

template <class T>
void SparseMatrixTemplate_RM<T>::inplaceMulCol(int j,T c)
{
  Assert(isValidCol(j));
  for(int i=0;i<m;i++) {
    RowIterator it=rows[i].find(j);
    if(it != rows[i].end()) it->second *= c;
  }
}

template <class T>
bool SparseMatrixTemplate_RM<T>::isValid() const
{
  if(m != (int)rows.size()) return false;
  for(size_t i=0;i<rows.size();i++) {
    if((int)rows[i].n != n) return false;
    if(!rows[i].isValid()) return false;
  }
  return true;
}

template <class T>
size_t SparseMatrixTemplate_RM<T>::numNonZeros() const
{
  size_t nnz=0;
  for(size_t i=0;i<rows.size();i++)
    nnz += rows[i].size();
  return nnz;
}

template <class T>
std::ostream& operator << (std::ostream& out, const SparseMatrixTemplate_RM<T>& A)
{
  out<<A.m<<" "<<A.n<<" "<<A.numNonZeros()<<endl;
  for(size_t i=0;i<A.rows.size();i++) {
    typename SparseMatrixTemplate_RM<T>::ConstRowIterator it;
    for(it=A.rows[i].begin();it!=A.rows[i].end();it++)
      out<<i<<" "<<it->first<<"   "<<it->second<<endl;
  }
  return out;
}

template <class T>
std::istream& operator >> (std::istream& in, SparseMatrixTemplate_RM<T>& A)
{
  int m,n,nnz;
  in >> m >> n >> nnz;
  if(in.bad()) return in;
  A.resize(m,n);
  for(int i=0;i<nnz;i++) {
    int row,col; T val;
    in >> row >> col >> val;
    if(in.bad()) return in;
    A(row,col) = val;
  }
  return in;
}


template class SparseMatrixTemplate_RM<float>;
template class SparseMatrixTemplate_RM<double>;
template class SparseMatrixTemplate_RM<Complex>;
template ostream& operator << (ostream& out, const SparseMatrixTemplate_RM<float>& v);
template ostream& operator << (ostream& out, const SparseMatrixTemplate_RM<double>& v);
template ostream& operator << (ostream& out, const SparseMatrixTemplate_RM<Complex>& v);
template istream& operator >> (istream& in, SparseMatrixTemplate_RM<float>& v);
template istream& operator >> (istream& in, SparseMatrixTemplate_RM<double>& v);
template istream& operator >> (istream& in, SparseMatrixTemplate_RM<Complex>& v);

template void SparseMatrixTemplate_RM<float>::copy(const SparseMatrixTemplate_RM<double>& a);
template void SparseMatrixTemplate_RM<double>::copy(const SparseMatrixTemplate_RM<float>& a);
template void SparseMatrixTemplate_RM<Complex>::copy(const SparseMatrixTemplate_RM<float>& a);
template void SparseMatrixTemplate_RM<Complex>::copy(const SparseMatrixTemplate_RM<double>& a);
