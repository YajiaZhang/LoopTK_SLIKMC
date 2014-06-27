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

#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include "VectorTemplate.h"
#include "VectorPrinter.h"

namespace Math {
  typedef VectorTemplate<Real> Vector;
}

#if 0
#include "math.h"
#include <myfile.h>
#include <stdincludes.h>

namespace Math {

class Vector
{
public:
	Vector();
	Vector(const Vector&);
	Vector(int n);
	Vector(int n, Real initval);
	Vector(int n, const Real* vals);
	~Vector();

	void clear();
	void resize(int n);
	void resize(int n, Real initval);
	inline int size() const { return n; }

	const Vector& operator = (const Vector&);
	bool operator == (const Vector&) const;
	inline bool operator != (const Vector& v) const;
	inline operator Real* ();
	inline operator const Real* () const;
	inline const Real& operator() (int) const;
	inline Real& operator() (int);
	inline const Real& operator[] (int) const;
	inline Real& operator[] (int);
	void operator += (const Vector&);
	void operator -= (const Vector&);
	void operator *= (Real);
	void operator /= (Real);

	void add(const Vector&, const Vector&);
	void sub(const Vector&, const Vector&);
	void mul(const Vector&, Real);
	void div(const Vector&, Real);
	void madd(const Vector&, Real);
	void axpby(Real a,const Vector& x,Real b,const Vector& y);

	void set(const Vector&);
	void set(Real);
	void set(const Real* vals);
	void setZero();
	void setNegative(const Vector&);
	void setNormalized(const Vector&);

	void inplaceNegative();
	void inplaceScale(Real);
	void inplaceNormalize();

	inline void get(Vector&) const;
	void get(Real* vals) const;
	inline void getNegative(Vector&) const;
	inline void getNormalized(Vector&) const;

	inline bool isEmpty() const;
	inline bool hasDims(int n) const;
	inline bool isValidIndex(int i) const;
	bool isZero(Real eps=Zero) const;
	bool isEqual(const Vector&,Real eps=Zero) const;

	Real dot(const Vector&) const;
	inline Real norm() const;
	inline Real normSquared() const;
	Real minElement(int* index=NULL) const;
	Real maxElement(int* index=NULL) const;
	Real minAbsElement(int* index=NULL) const;
	Real maxAbsElement(int* index=NULL) const;

	void print(std::ostream& out=std::cout,char delim=' ',char bracket='[') const;
	bool load(File&);
	bool save(File&);

	int n;
	Real* vals;
};

inline bool FuzzyEquals(const Vector&,const Vector&,Real eps=Epsilon);
inline bool FuzzyZero(const Vector&, Real eps=Epsilon);
inline Real dot(const Vector&, const Vector&);
inline Real norm(const Vector&);
//the following are very inefficient, but included for completeness
inline Vector operator + (const Vector&, const Vector&);
inline Vector operator - (const Vector&, const Vector&);
inline Vector operator * (const Vector&, Real);
inline Vector operator * (Real, const Vector&);
inline Vector operator / (const Vector&, Real);
inline Vector operator / (Real, const Vector&);

std::ostream& operator << (std::ostream&, const Vector&);
std::istream& operator >> (std::istream&, Vector&);



//Vector inlined methods

inline bool Vector::operator != (const Vector& v) const
{
  return !operator==(v);
}


inline Vector::operator Real* ()
{
	return vals;
}

inline Vector::operator const Real* () const
{
	return vals;
}

inline const Real& Vector::operator() (int i) const
{
	return operator [] (i);
}

inline Real& Vector::operator() (int i)
{
	return operator [] (i);
}

inline const Real& Vector::operator[] (int i) const
{
#ifdef _DEBUG
	if(!isValidIndex(i))
		Error("Vector index out of range");
#endif
	return vals[i];

}

inline Real& Vector::operator[] (int i)
{
#ifdef _DEBUG
	if(!isValidIndex(i))
		Error("Vector index out of range");
#endif
	return vals[i];
}

inline void Vector::get(Vector& v) const
{
	v.set(*this);
}

inline void Vector::getNegative(Vector& v) const
{
	v.setNegative(*this);
}

inline void Vector::getNormalized(Vector& v) const
{
	v.setNormalized(*this);
}

inline bool Vector::isEmpty() const
{
	return n == 0;
}

inline bool Vector::hasDims(int _n) const
{
	return n == _n;
}

inline bool Vector::isValidIndex(int i) const
{
	return i >= 0 && i < n;
}

inline Real Vector::norm() const
{
	return Sqrt(normSquared());
}

inline Real norm(const Vector& a)
{
	return a.norm();
}

inline Real Vector::normSquared() const
{
	return dot(*this);
}

inline bool FuzzyEquals(const Vector& a,const Vector& b,Real eps)
{
  return a.isEqual(b,eps);
}

inline bool FuzzyZero(const Vector& a, Real eps)
{
  return a.isZero(eps);
}

inline Real dot(const Vector& a, const Vector& b)
{
	return a.dot(b);
}


inline Vector operator + (const Vector& a, const Vector& b)
{
	Vector v;
	v.add(a,b);
	return v;
}

inline Vector operator - (const Vector& a, const Vector& b)
{
	Vector v;
	v.sub(a,b);
	return v;
}

inline Vector operator * (const Vector& a, Real c)
{
	Vector v;
	v.mul(a,c);
	return v;
}

inline Vector operator * (Real c, const Vector& a)
{
	Vector v;
	v.mul(a,c);
	return v;
}

inline Vector operator / (const Vector& a, Real c)
{
	Vector v;
	v.div(a,c);
	return v;
}

inline Vector operator / (Real c, const Vector& a)
{
	Vector v;
	v.div(a,c);
	return v;
}

} //namespace Math

namespace std
{
  template<> inline void swap(Math::Vector& a, Math::Vector& b)
  {
    Math::Real* vals_tmp = a.vals;
    int n_tmp = a.n;
    a.vals = b.vals;
    a.n = b.n;
    b.vals = vals_tmp;
    b.n = n_tmp;
  }
} //namespace std

#endif

#endif
