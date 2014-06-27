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

#ifndef MATH_VAL_ARRAY_H
#define MATH_VAL_ARRAY_H

#include "math.h"
#include <myfile.h>
#include <iostream>

namespace Math {

/* CompactValArray
 * A very basic container for pointer data.  Intended for use with
 *   other classes.
 * Most operators (besides allocate/deallocate, shallow) assume 
 *   correctly sized arguments and destination.
 * Allocate must only be used on an uninitialized object.  It only
 *   sets the capacity.  Size must be adjusted manually.
 * Deallocate uninitializes an allocated object.
 * Clear uninitializes an object without deallocation (using it on an 
 *   allocated object will result in a memory leak!)
 */
template <class T>
class CompactValArray
{
public:
  typedef CompactValArray<T> MyT;

  CompactValArray();
  CompactValArray(const CompactValArray&);
  CompactValArray(int size);
  CompactValArray(int size,T);
  CompactValArray(int size,const T*);
  ~CompactValArray();

  void resize(int size);
  void resize(int size,T initVal);
  void clear();

  bool operator == (const MyT&) const;
  inline bool operator != (const MyT& a) const { return !operator==(a); }
  inline T& operator[](int i) { return vals[i]; }
  inline const T& operator[](int i) const { return vals[i]; }
  inline int n() const { return size; }
  inline T* ptr() const { return vals; }

  void copy(const MyT&);
  void copy(const T* vals);
  void copySubArray(int i,const MyT&);
  void swap(MyT&);
  void add(const MyT&, const MyT&);
  void sub(const MyT&, const MyT&);
  void mul(const MyT&, T);
  void div(const MyT&, T);
  void axpby(T a,const MyT& x,T b,const MyT& y);
  //the following are increment-type methods
  void inc(const T&);
  void inc(const MyT&);
  void dec(const MyT&);
  void madd(const MyT&, T);

  void set(T);
  void setZero();
  void setNegative(const MyT&);

  void inplaceNegative();
  void inplaceMul(T);
  void inplaceDiv(T);

  void getCopy(MyT&) const;
  void getCopy(T* vals) const;
  void getSubArrayCopy(int i,MyT&) const;

  inline bool isEmpty() const { return vals==NULL; }
  inline bool hasSize(int s) const { return s==size; }
  inline bool isValidIndex(int i) const { return i>=0 && i<size; }
  bool isZero(T eps=0) const;
  bool isEqual(const MyT&,T eps=0) const;

  T dot(const MyT&) const;
  T dotSelf() const;
  T minElement(int* index=NULL) const;
  T maxElement(int* index=NULL) const;
  T minAbsElement(int* index=NULL) const;
  T maxAbsElement(int* index=NULL) const;

  bool load(File&);
  bool save(File&) const;

protected:
  T* vals;
  int size;
};

class Complex;
typedef class CompactValArray<float> fCompactValArray;
typedef class CompactValArray<double> dCompactValArray;
typedef class CompactValArray<Complex> cCompactValArray;

template <class T>
std::ostream& operator << (std::ostream&, const CompactValArray<T>&);
template <class T>
std::istream& operator >> (std::istream&, CompactValArray<T>&);




template <class T>
inline bool FuzzyEquals(const CompactValArray<T>& a, const CompactValArray<T>& b,T eps)
{
  return a.isEqual(b,eps);
}

template <class T>
inline bool FuzzyZero(const CompactValArray<T>& a,T eps)
{
  return a.isZero(eps);
}

template <class T>
inline T dot(const CompactValArray<T>& a, const CompactValArray<T>& b)
{
	return a.dot(b);
}

} //namespace Math

namespace std
{
  template<class T> inline void swap(Math::CompactValArray<T>& a, Math::CompactValArray<T>& b)
  {
    a.swap(b);
  }
} //namespace std

template <class T>
inline bool ReadFile(File&f,Math::CompactValArray<T>& v) { return v.load(f); }
template <class T>
inline bool WriteFile(File&f,const Math::CompactValArray<T>& v) { return v.save(f); }

#endif
