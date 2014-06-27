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

#ifndef MATH_SPARSE_VECTOR_H
#define MATH_SPARSE_VECTOR_H

#include "VectorTemplate.h"
#include <structs/SparseArray.h>

namespace Math {

template <class T>
class SparseVector : public SparseArray<T>
{
public:
  typedef SparseVector<T> MyT;
  typedef SparseArray<T> BaseT;
  typedef typename BaseT::iterator iterator;
  typedef typename BaseT::const_iterator const_iterator;
  typedef VectorTemplate<T> VectorT;

  SparseVector() {}
  SparseVector(size_t n):BaseT(n) {}
  SparseVector(const BaseT& v):BaseT(v) {}

  void print(std::ostream&) const;

  T operator() (int i) const;
  inline void set(int i,const T& t) { BaseT::insert(i,t); }

  inline void setZero() { BaseT::entries.clear(); }
  inline void set(const BaseT& v) { BaseT::operator = (v); }
  void set(const VectorT&,T zeroTol=Zero);
  void set(const T*,int n,T zeroTol=Zero);
  void get(T*) const;
  void inplaceNegative();
  void inplaceMul(T s);
  void inplaceDiv(T s);

  void add(const MyT&, const MyT&);
  void sub(const MyT&, const MyT&);
  void mul(const MyT&, T s);
  void div(const MyT&, T s);

  T dot(const VectorT&) const;
  T dot(const MyT&) const;
  T norm() const;
  T normSquared() const;
  T distance(const MyT&) const;
  T distanceSquared(const MyT&) const;

  inline bool isEmpty() const { return BaseT::empty(); }
  inline bool hasDims(size_t _n) const { return BaseT::n==_n; }

  T minElement(int* index=NULL) const;
  T maxElement(int* index=NULL) const;
  T minAbsElement(int* index=NULL) const;
  T maxAbsElement(int* index=NULL) const;
};

template <class T>
class SparseVectorCompressed 
{
public:
  typedef SparseVectorCompressed<T> MyT;
  typedef VectorTemplate<T> VectorT;

  SparseVectorCompressed();
  SparseVectorCompressed(int n, int num_entries);
  SparseVectorCompressed(const MyT&);
  ~SparseVectorCompressed();
  void init(int n, int num_entries);
  void resize(int n, int num_entries);
  void makeSimilar(const MyT&);
  void cleanup();

  bool Read(File&);
  bool Write(File&) const;
  void print(std::ostream&) const;

  const MyT& operator =(const MyT&);
  T operator() (int i) const;

  void setZero();
  void set(const MyT&);
  void set(const VectorT&,T zeroTol=Zero);
  void set(const T*,int n,T zeroTol=Zero);
  void get(T*) const;
  void inplaceNegative();
  void inplaceMul(T s);
  void inplaceDiv(T s);

  void add(const MyT&, const MyT&);
  void sub(const MyT&, const MyT&);
  void mul(const MyT&, T s);
  void div(const MyT&, T s);

  T dot(const VectorT&) const;
  T dot(const MyT&) const;
  T norm() const;
  T normSquared() const;
  T distance(const MyT&) const;
  T distanceSquared(const MyT&) const;

  bool isValid() const;
  inline bool isEmpty() const { return n == 0; }
  inline bool hasDims(int _n) const { return n==_n; }
  inline bool isValidIndex(int i) const { return 0<=i&&i<n; }

  T minElement(int* index=NULL) const;
  T maxElement(int* index=NULL) const;
  T minAbsElement(int* index=NULL) const;
  T maxAbsElement(int* index=NULL) const;

  /***********************************************************
   * Compressed vector format:
   * There are num_entries nonzero entries .
   * The index of entry i is indices[i] (0<=indices[i]<indices[i+1]<n)
   * The value of entry i is vals[i]
   **********************************************************/
  int* indices;
  T* vals;
  int num_entries;
  int n;
};

} //namespace Math

#endif