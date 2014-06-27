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

#include "sparsevector.h"
#include "fastarray.h"
#include <algorithm>
#include <assert.h>
#include "errors.h"
using namespace Math;
using namespace std;

template <class T>
void SparseVector<T>::print(std::ostream& out) const
{
  for(const_iterator i=this->begin();i!=this->end();i++)
    out<<i->first<<":"<<i->second<<" ";
  out<<endl;
}

template <class T>
T SparseVector<T>::operator() (int i) const
{
  const_iterator it=this->entries.find(i);
  if(it!=this->end()) return it->second;
  return 0;
}

template <class T>
void SparseVector<T>::set(const VectorT& v,T zeroTol)
{
  BaseT::resize(v.n);
  this->entries.clear();
  for(int i=0;i<v.n;i++) {
    if(!FuzzyZero(v(i),zeroTol)) push_back(i,v(i));
  }
}

template <class T>
void SparseVector<T>::set(const T* v,int n,T zeroTol)
{
  this->resize(n);
  this->entries.clear();
  for(int i=0;i<n;i++) {
    if(!FuzzyZero(v[i],zeroTol)) push_back(i,v[i]);
  }
}

template <class T>
void SparseVector<T>::get(T* v) const
{
  int k=0;
  for(const_iterator i=this->begin();i!=this->end();i++) {
    while(k < i->first) { v[k]=0; k++; }
    v[k] = i->second;
    k=i->first+1;
  }
  while(k < (int)this->n) { v[k]=0; k++; }
}

template <class T>
void SparseVector<T>::inplaceNegative()
{
  for(iterator i=this->begin();i!=this->end();i++)
    i->second = -i->second;
}

template <class T>
void SparseVector<T>::inplaceMul(T s)
{
  for(iterator i=this->begin();i!=this->end();i++)
    i->second *= s;
}

template <class T>
void SparseVector<T>::inplaceDiv(T s)
{
  for(iterator i=this->begin();i!=this->end();i++)
    i->second /= s;
}

/*
template <class T>
void SparseVector<T>::add(const MyT& a, const MyT& b);
template <class T>
void SparseVector<T>::sub(const MyT&, const MyT&);
*/

template <class T>
void SparseVector<T>::mul(const MyT& a, T s)
{
  operator = (a);
  inplaceMul(s);
}

template <class T>
void SparseVector<T>::div(const MyT& b, T s)
{
  operator = (b);
  inplaceDiv(s);
}

template <class T>
T SparseVector<T>::dot(const VectorT& v) const
{
  T sum=0;
  Assert(v.n == (int)this->n);
  for(const_iterator i=this->begin();i!=this->end();i++)
    sum += i->second*v(i->first);
  return sum;
}

template <class T>
T SparseVector<T>::dot(const MyT& b) const
{
  Assert(this->n==b.n);
  const_iterator k=this->begin(),bk=b.begin();
  int i,bi;
  T sum=Zero;
  while(k!=this->end() || bk!=b.end()) {
    i = (k!=this->end()? k->first : this->n);
    bi = (bk!=b.end()? bk->first : this->n);
    if(i < bi) k++;
    else if (bi < i) bk++;
    else { 
      sum+=k->second*bk->second;
      k++; bk++; 
    }
  }
  return sum;
}

template <class T>
T SparseVector<T>::norm() const { return Sqrt(normSquared()); }

template <class T>
T SparseVector<T>::normSquared() const
{
  T sum=0;
  for(const_iterator i=this->begin();i!=this->end();i++)
    sum += Sqr(i->second);
  return sum;
}

template <class T>
T SparseVector<T>::distance(const MyT& b) const { return Sqrt(distanceSquared(b)); }

template <class T>
T SparseVector<T>::distanceSquared(const MyT& b) const
{
  return normSquared() + b.normSquared() - Two*dot(b);
}

template <class T>
T SparseVector<T>::minElement(int* index) const
{
  T vmin=Inf;
  int zeroIndex=-1;
  for(const_iterator i=this->begin();i!=this->end();i++) {
    if(i->second < vmin) {
      vmin = i->second;
      if(index) *index=i->first;
    }
    if(zeroIndex == -1 && i!=this->begin()) {
      const_iterator prev=i; --prev;
      if(prev->first < i->first-1) //a zero element between i-1 and i
	zeroIndex = i->first-1;
    }
  }
  if(vmin < 0) return vmin;
  if(index) *index = zeroIndex;
  return 0;
}

template <class T>
T SparseVector<T>::maxElement(int* index) const
{
  T vmax=-Inf;
  int zeroIndex=-1;
  for(const_iterator i=this->begin();i!=this->end();i++) {
    if(i->second > vmax) {
      vmax = i->second;
      if(index) *index=i->first;
    }
    if(zeroIndex == -1 && i!=this->begin()) {
      const_iterator prev=i; --prev;
      if(prev->first < i->first-1) //a zero element between i-1 and i
	zeroIndex = i->first-1;
    }
  }
  if(vmax > 0) return vmax;
  if(index) *index = zeroIndex;
  return 0;
}

template <class T>
T SparseVector<T>::minAbsElement(int* index) const
{
  T vmin=Inf;
  int zeroIndex=-1;
  for(const_iterator i=this->begin();i!=this->end();i++) {
    if(Abs(i->second) < vmin) {
      vmin = Abs(i->second);
      if(index) *index=i->first;
    }
    if(zeroIndex == -1 && i!=this->begin()) {
      const_iterator prev=i; --prev;
      if(prev->first < i->first-1) //a zero element between i-1 and i
	zeroIndex = i->first-1;
    }
  }
  if(zeroIndex == -1) return vmin;
  if(index) *index = zeroIndex;
  return Zero;
}

template <class T>
T SparseVector<T>::maxAbsElement(int* index) const
{
  T vmax=-Inf;
  for(const_iterator i=this->begin();i!=this->end();i++) {
    if(Abs(i->second) > vmax) {
      vmax = Abs(i->second);
      if(index) *index=i->first;
    }
  }
  return vmax;
}

template<class T>
SparseVectorCompressed<T>::SparseVectorCompressed()
  :indices(NULL),vals(NULL),num_entries(0),n(0)
{}

template<class T>
SparseVectorCompressed<T>::SparseVectorCompressed(int _n, int _num_entries)
  :indices(NULL),vals(NULL),num_entries(0),n(0)
{
  init(_n,_num_entries);
}

template<class T>
SparseVectorCompressed<T>::SparseVectorCompressed(const MyT& v)
{
  set(v);
}

template<class T>
SparseVectorCompressed<T>::~SparseVectorCompressed()
{
  cleanup();
}

template<class T>
void SparseVectorCompressed<T>::init(int _n, int _num_entries)
{
  SafeArrayDelete(indices);
  SafeArrayDelete(vals);
  n=_n;
  num_entries=_num_entries;
  indices=new int[num_entries];
  vals=new T[num_entries];
}

template<class T>
void SparseVectorCompressed<T>::resize(int _n, int _num_entries)
{
  if(num_entries != _num_entries)
    init(_n,_num_entries);
  n=_n;
}

template<class T>
void SparseVectorCompressed<T>::makeSimilar(const MyT& v)
{
  resize(v.n,v.num_entries);
  for(int i=0;i<num_entries;i++) indices[i]=v.indices[i];
}

template<class T>
void SparseVectorCompressed<T>::cleanup()
{
  SafeArrayDelete(indices);
  SafeArrayDelete(vals);
  n=num_entries=0;
}

template<class T>
const SparseVectorCompressed<T>& SparseVectorCompressed<T>::operator =(const MyT& v)
{
  set(v);
  return *this;
}

template<class T>
T SparseVectorCompressed<T>::operator() (int i) const
{
  int* ind=std::lower_bound(indices,indices+num_entries,i);
  if((ind!=indices+num_entries)&&(*ind==i)) return vals[ind-indices];
  return Zero;
}

template<class T>
void SparseVectorCompressed<T>::setZero()
{
  array_zero(vals,num_entries);
}

template<class T>
void SparseVectorCompressed<T>::set(const MyT& v)
{
  makeSimilar(v);
  array_equal(vals,v.vals,num_entries);
}

template<class T>
void SparseVectorCompressed<T>::set(const VectorT& v,T zeroTol)
{
  int nnz=0;
  for(int i=0;i<this->n;i++)
    if(!FuzzyZero(v[i],zeroTol)) nnz++;
  resize(n,nnz);
  nnz=0;
  for(int i=0;i<this->n;i++) {
    if(!FuzzyZero(v[i],zeroTol)) {
      vals[nnz] = v[i];
      indices[nnz]=i;
      nnz++;
    }
  }
}

template<class T>
void SparseVectorCompressed<T>::set(const T* v,int n,T zeroTol)
{
  int nnz=0;
  for(int i=0;i<this->n;i++)
    if(!FuzzyZero(v[i],zeroTol)) nnz++;
  resize(n,nnz);
  nnz=0;
  for(int i=0;i<this->n;i++) {
    if(!FuzzyZero(v[i],zeroTol)) {
      vals[nnz] = v[i];
      indices[nnz]=i;
      nnz++;
    }
  }
}

template<class T>
void SparseVectorCompressed<T>::get(T* v) const
{
  int i=0;
  for(int k=0;k<num_entries;k++) {
    for(;i<indices[k];i++)
      v[i]=Zero;
    v[i]=vals[k];
  }
  for(;i<this->n;i++)
    v[i]=Zero;
}

template<class T>
void SparseVectorCompressed<T>::inplaceNegative()
{
  for(int i=0;i<num_entries;i++)
    vals[i]=-vals[i];
}

template<class T>
void SparseVectorCompressed<T>::inplaceMul(T s)
{
  array_scale(vals,s,num_entries);
}

template<class T>
void SparseVectorCompressed<T>::inplaceDiv(T s)
{
  array_scale(vals,Inv(s),num_entries);
}

template<class T>
void SparseVectorCompressed<T>::add(const MyT& a, const MyT& b)
{
  Assert(this != &a && this != &b);
  Assert(a.n == b.n);
  Assert(a.isValid() && b.isValid());
  int ak=0,bk=0;
  int ai,bi;
  int nnz=0;

  while(ak<a.num_entries || bk<b.num_entries) {
    ai = (ak<a.num_entries? a.indices[ak] : n);
    bi = (bk<b.num_entries? b.indices[bk] : n);
    if(ai < bi) ak++;
    else if (bi < ai) bk++;
    else { ak++; bk++; }
    nnz++;
  }
  resize(a.n,nnz);

  nnz=0; ak=0;bk=0;
  while(ak<a.num_entries || bk<b.num_entries) {
    ai = (ak<a.num_entries? a.indices[ak] : n);
    bi = (bk<b.num_entries? b.indices[bk] : n);
    if(ai < bi) {
      indices[nnz] = ai;
      vals[nnz] = a.vals[ak];
      ak++;
    }
    else if(bi < ai) {
      indices[nnz] = bi;
      vals[nnz] = b.vals[bk];
      bk++;
    }
    else {
      indices[nnz] = ai;
      vals[nnz] = a.vals[ak]+b.vals[bk];
      ak++; bk++;
    }
    nnz++;
  }
}

template<class T>
void SparseVectorCompressed<T>::sub(const MyT& a, const MyT& b)
{
  Assert(this != &a && this != &b);
  Assert(a.n == b.n);
  Assert(a.isValid() && b.isValid());
  int ak=0,bk=0;
  int ai,bi;
  int nnz=0;

  while(ak<a.num_entries || bk<b.num_entries) {
    ai = (ak<a.num_entries? a.indices[ak] : n);
    bi = (bk<b.num_entries? b.indices[bk] : n);
    if(ai < bi) ak++;
    else if (bi < ai) bk++;
    else { ak++; bk++; }
    nnz++;
  }
  resize(a.n,nnz);

  nnz=0; ak=0;bk=0;
  while(ak<a.num_entries || bk<b.num_entries) {
    ai = (ak<a.num_entries? a.indices[ak] : n);
    bi = (bk<b.num_entries? b.indices[bk] : n);
    if(ai < bi) {
      indices[nnz] = ai;
      vals[nnz] = a.vals[ak];
      ak++;
    }
    else if(bi < ai) {
      indices[nnz] = bi;
      vals[nnz] = b.vals[bk];
      bk++;
    }
    else {
      indices[nnz] = ai;
      vals[nnz] = a.vals[ak]-b.vals[bk];
      ak++; bk++;
    }
    nnz++;
  }
}

template<class T>
void SparseVectorCompressed<T>::mul(const MyT& v, T s)
{
  makeSimilar(v);
  array_multiply(vals,v.vals,s,num_entries);
}

template<class T>
void SparseVectorCompressed<T>::div(const MyT& v, T s)
{
  makeSimilar(v);
  array_multiply(vals,v.vals,Inv(s),num_entries);
}

template<class T>
T SparseVectorCompressed<T>::dot(const VectorT& v) const
{
  Assert(n==v.n);
  T sum=Zero;
  for(int i=0;i<num_entries;i++) sum+=v(indices[i])*vals[i];
  return sum;
}

template<class T>
T SparseVectorCompressed<T>::dot(const MyT& b) const
{
  Assert(n==b.n);
  int k=0,bk=0;
  int i,bi;
  T sum=Zero;

  while(k<num_entries || bk<b.num_entries) {
    i = (k<num_entries? indices[k] : n);
    bi = (bk<b.num_entries? b.indices[bk] : n);
    if(i < bi) k++;
    else if (bi < i) bk++;
    else { 
      sum+=vals[k]*b.vals[bk];
      k++; bk++; 
    }
  }
  return sum;
}

template<class T>
T SparseVectorCompressed<T>::normSquared() const
{
  return array_dot(vals,vals,num_entries);
}

template<class T>
bool SparseVectorCompressed<T>::isValid() const
{
  for(int i=0;i<num_entries;i++) {
    if(indices[i] < 0 || indices[i] >= n) return false;
    if(i!=0 && (indices[i-1] >= indices[i])) return false;
  }
  return true;
}

template<class T>
T SparseVectorCompressed<T>::minElement(int* index) const
{
  T vmin=Inf;
  int zeroIndex=-1;
  for(int i=0;i<num_entries;i++) {
    if(vals[i] < vmin) {
      vmin = vals[i];
      if(index) *index=indices[i];
    }
    if(zeroIndex == -1 && i!=0) {
      if(indices[i-1] < indices[i]-1) //a zero element between i-1 and i
	zeroIndex = indices[i]-1;
    }
  }
  if(vmin < Zero) return vmin;
  if(index) *index = zeroIndex;
  return Zero;
}

template<class T>
T SparseVectorCompressed<T>::maxElement(int* index) const
{
  T vmax=-Inf;
  int zeroIndex=-1;
  for(int i=0;i<num_entries;i++) {
    if(vals[i] < vmax) {
      vmax = vals[i];
      if(index) *index=indices[i];
    }
    if(zeroIndex == -1 && i!=0) {
      if(indices[i-1] < indices[i]-1) //a zero element between i-1 and i
	zeroIndex = indices[i]-1;
    }
  }
  if(vmax > Zero) return vmax;
  if(index) *index = zeroIndex;
  return Zero;
}

template<class T>
T SparseVectorCompressed<T>::minAbsElement(int* index) const
{
  T vmin=Inf;
  int zeroIndex=-1;
  for(int i=0;i<num_entries;i++) {
    if(Abs(vals[i]) < vmin) {
      vmin = Abs(vals[i]);
      if(index) *index=indices[i];
    }
    if(zeroIndex == -1 && i!=0) {
      if(indices[i-1] < indices[i]-1) //a zero element between i-1 and i
	zeroIndex = indices[i]-1;
    }
  }
  if(zeroIndex == -1) return vmin;
  if(index) *index = zeroIndex;
  return Zero;
}

template<class T>
T SparseVectorCompressed<T>::maxAbsElement(int* index) const
{
  T vmax=-Inf;
  for(int i=0;i<num_entries;i++) {
    if(Abs(vals[i]) > vmax) {
      vmax = Abs(vals[i]);
      if(index) *index=indices[i];
    }
  }
  return vmax;
}

template<class T>
void SparseVectorCompressed<T>::print(std::ostream& out) const
{
  for(int i=0;i<num_entries;i++) {
    out<<"("<<indices[i]<<" , "<<vals[i]<<") ";
  }
  out<<endl;
}

template class SparseVector<float>;
template class SparseVector<double>;
template class SparseVectorCompressed<float>;
template class SparseVectorCompressed<double>;
