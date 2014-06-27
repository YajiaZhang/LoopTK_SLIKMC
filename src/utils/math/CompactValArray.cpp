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

#include "CompactValArray.h"
#include "fastarray.h"
#include "complex.h"
using namespace std;

namespace Math {

#define CHECKEMPTY() { Assert(!hasSize(0)); }

template <class T>
CompactValArray<T>::CompactValArray()
:vals(NULL),size(0)
{}

template <class T>
CompactValArray<T>::CompactValArray(const MyT& v)
:vals(NULL),size(0)
{
  copy(v);
}

template <class T>
CompactValArray<T>::CompactValArray(int n)
:vals(NULL),size(0)
{
  resize(n);
}

template <class T>
CompactValArray<T>::CompactValArray(int n,T initval)
:vals(NULL),size(0)
{
  resize(n,initval);
}

template <class T>
CompactValArray<T>::CompactValArray(int n,const T* vals)
:vals(NULL),size(0)
{
  resize(n);
  copy(vals);
}

template <class T>
CompactValArray<T>::~CompactValArray()
{
  clear();
}

template <class T>
void CompactValArray<T>::resize(int _n)
{
  if(_n > size) {
    clear();
    vals = new T[_n];
  }
  size = _n;
}

template <class T>
void CompactValArray<T>::resize(int _n,T initval)
{
  resize(_n);
  set(initval);
}

template <class T>
void CompactValArray<T>::clear()
{
  SafeArrayDelete(vals);
  size = 0;
}

template <class T>
bool CompactValArray<T>::operator == (const MyT& v) const
{
  if(this == &v) return true;
  if(!hasSize(v.size)) return false;
  T* va=vals;
  T* vb=v.vals;
  for(int i=0;i<size;i++,va++,vb++)
    if(*va != *vb) return false;
  return true;
}

template <class T>
void CompactValArray<T>::swap(MyT& a)
{
  std::swap(vals,a.vals);
  std::swap(size,a.size);
}

template <class T>
void CompactValArray<T>::copy(const MyT& a)
{
  if(this == &a) return;
  Assert(hasSize(a.size));
	array_equal(vals, a.vals, size);
}

template <class T>
void CompactValArray<T>::copy(const T* _vals)
{
	CHECKEMPTY();
  array_equal(vals, _vals, size);
}

template <class T>
void CompactValArray<T>::copySubArray(int i,const MyT& a)
{
  Assert(this != &a);
  Assert(isValidIndex(i));
  Assert(isValidIndex(i+a.size-1));
	array_equal(vals+i,a.vals, a.size);
}

template <class T>
void CompactValArray<T>::add(const MyT& a, const MyT& b)
{
  Assert(a.size == b.size);
  Assert(size == a.size);
	array_add(vals, a.vals, b.vals, size);
}

template <class T>
void CompactValArray<T>::sub(const MyT& a, const MyT& b)
{
  Assert(a.size == b.size);
  Assert(size == a.size);
	array_sub(vals, a.vals, b.vals, size);
}

template <class T>
void CompactValArray<T>::mul(const MyT& a, T c)
{
  Assert(size == a.size);
	array_multiply(vals,a.vals, c, size);
}

template <class T>
void CompactValArray<T>::div(const MyT& a, T c)
{
	mul(a,Inv(c));
}

template<class T>
void CompactValArray<T>::axpby(T a,const MyT& x,T b,const MyT& y)
{
  Assert(x.size == y.size);
  Assert(size == x.size);
  array_axpby(vals, a, x.vals, b, y.vals,size);
}

template <class T>
void CompactValArray<T>::inc(const T& c)
{
	CHECKEMPTY();
  gen_array_acc(vals,1, &c,0, size);
}

template <class T>
void CompactValArray<T>::inc(const MyT& v)
{
	Assert(hasSize(v.size));
  array_acc(vals, v.vals, size);
}

template <class T>
void CompactValArray<T>::dec(const MyT& v)
{
	Assert(hasSize(v.size));
  array_dec(vals, v.vals, size);
}

template <class T>
void CompactValArray<T>::madd(const MyT& v, T c)
{
	Assert(hasSize(v.size));
  array_madd(vals, v.vals, c,size);
}






template <class T>
void CompactValArray<T>::set(T c)
{
	CHECKEMPTY();
  array_fill(vals,c,size);
}

template <class T>
void CompactValArray<T>::setZero()
{
	set((T)0);
}

template <class T>
void CompactValArray<T>::setNegative(const MyT& a)
{
	Assert(hasSize(a.size));
  array_negate(vals, a.vals, size);
}




template <class T>
void CompactValArray<T>::inplaceNegative()
{
  CHECKEMPTY();
  array_negate(vals, vals, size);
}

template <class T>
void CompactValArray<T>::inplaceMul(T c)
{
	CHECKEMPTY();
  array_scale(vals,c,size);
}

template <class T>
void CompactValArray<T>::inplaceDiv(T c)
{
	CHECKEMPTY();
  array_scale(vals,Inv(c),size);
}

template <class T>
void CompactValArray<T>::getCopy(MyT& v) const
{
  v.copy(*this);
}

template <class T>
void CompactValArray<T>::getCopy(T* _vals) const
{
  CHECKEMPTY();
  array_equal(_vals,vals,size);
}

template <class T>
void CompactValArray<T>::getSubArrayCopy(int i,MyT& a) const
{
  Assert(&a != this);
  Assert(isValidIndex(i));
  Assert(isValidIndex(i+a.size-1));
	array_equal(a.vals, vals+i, a.size);
}



template <class T>
bool CompactValArray<T>::isZero(T eps) const
{
  CHECKEMPTY();
  T* v=vals;
  for(int i=0;i<size;i++,v++)
    if(!FuzzyZero(*v,eps)) return false;
  return true;
}

template <class T>
bool CompactValArray<T>::isEqual(const MyT& a,T eps) const
{
	CHECKEMPTY();
	Assert(hasSize(a.size));
  T* v=vals;
  T* va=a.vals;
  for(int i=0;i<size;i++,v++,va++)
    if(!FuzzyEquals(*v,*va,eps)) return false;
  return true;
}


template <class T>
T CompactValArray<T>::dot(const MyT& a) const
{
  CHECKEMPTY();
  Assert(hasSize(a.size));
	return array_dot(vals, a.vals, size);
}

template <class T>
T CompactValArray<T>::dotSelf() const 
{
  CHECKEMPTY();
  return array_norm_squared(vals,size);
}

template <class T>
T CompactValArray<T>::minElement(int* index) const
{
	CHECKEMPTY();
  T* v=vals;
	T b=*v;
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(*v < b) {
	    b=*v;
	    if(index) *index=i;
	  }
	return b;
}

template <class T>
T CompactValArray<T>::maxElement(int *index) const
{
	CHECKEMPTY();
  T* v=vals;
	T b=*v;
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(*v > b) {
	    b=*v;
	    if(index) *index=i;
	  }
	return b;
}

template <class T>
T CompactValArray<T>::minAbsElement(int* index) const
{
	CHECKEMPTY();
  T* v=vals;
	T b=Abs(*v);
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(Abs(*v) < b) {
	    b=Abs(*v);
	    if(index) *index=i;
	  }
	return b;
}

template <class T>
T CompactValArray<T>::maxAbsElement(int *index) const
{
	CHECKEMPTY();
  T* v=vals;
	T b=Abs(*v);
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(Abs(*v) > b) {
	    b=Abs(*v);
	    if(index) *index=i;
	  }
	return b;
}


template <class T>
bool CompactValArray<T>::load(File& f)
{
  int _n;
  if(!ReadFile(f,_n)) return false;
  resize(_n);
  if(!ReadArrayFile(f,vals,_n)) return false;
  return true;
}

template <class T>
bool CompactValArray<T>::save(File& f) const
{
  if(!WriteFile(f,size)) return false;
  if(!WriteArrayFile(f,vals,size)) return false;
  return true;
}

template <class T>
ostream& operator << (ostream& out, const CompactValArray<T>& v)
{
  cout<<v.n()<<"\t";
	for(int i=0; i<v.n(); i++)
		out << v[i] << " ";
	return out;
}

template <class T>
istream& operator >> (istream& in, CompactValArray<T>& v)
{
  int n;
  in >> n;
  v.resize(n);
	for(int i=0; i<n; i++)
		in >> v[i];
	return in;
}


//template instantiation for Complex
template<> bool CompactValArray<Complex>::isZero(Complex eps) const
{
  CHECKEMPTY();
  Complex* v=vals;
  for(int i=0;i<size;i++,v++)
    if(!FuzzyZero(*v,eps.x)) return false;
  return true;
}

template<> bool CompactValArray<Complex>::isEqual(const MyT& a,Complex eps) const
{
	CHECKEMPTY();
	Assert(hasSize(a.size));
  Complex* v=vals;
  Complex* va=a.vals;
  for(int i=0;i<size;i++,v++,va++)
    if(!FuzzyEquals(*v,*va,eps.x)) return false;
  return true;
}

template<> Complex CompactValArray<Complex>::minElement(int* index) const
{
  AssertNotReached();
  return Zero;
}

template<> Complex CompactValArray<Complex>::maxElement(int* index) const
{
  AssertNotReached();
  return Zero;
}

template<> Complex CompactValArray<Complex>::minAbsElement(int* index) const
{
	CHECKEMPTY();
  Complex* v=vals;
	Real b=Abs(*v);
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(Abs(*v) < b) {
	    b=Abs(*v);
	    if(index) *index=i;
	  }
	return b;
}

template<> Complex CompactValArray<Complex>::maxAbsElement(int* index) const
{
	CHECKEMPTY();
  Complex* v=vals;
	Real b=Abs(*v);
	if(index) *index=0;
	for(int i=1;i<size;i++,v++)
	  if(Abs(*v) > b) {
	    b=Abs(*v);
	    if(index) *index=i;
	  }
	return b;
}

template class CompactValArray<float>;
template class CompactValArray<double>;
template class CompactValArray<Complex>;


} //namespace Math
