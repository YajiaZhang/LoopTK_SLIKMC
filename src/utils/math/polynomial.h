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

#ifndef MATH_POLYNOMIAL_H
#define MATH_POLYNOMIAL_H

#include "misc.h"

namespace Math {

namespace Polynomial {

//alias for x^i
template <class T>
inline T BasisTerm(const T x, int i) { return IntegerPower(x,i); }

//kth derivative of x^i
template <class T>
T BasisTerm_Dk(const T x, int i,int k)
{
	if(k>i) return 0;
	return (T)FactorialTruncated(i,k)*BasisTerm(x,i-k);
}

//fills in p[i] with x^i, p must have size degree+1
template <class T>
void Basis(const T x, int degree, T p[])
{
	p[0]=(T)1;
	for(int i=1;i<=degree;i++) p[i]=x*p[i-1];
}

//fills in p[i] with kth derivative of x^i, p must have size degree+1
template <class T>
void Basis_Dx(const T x, int degree, T p[])
{
	Basis(x,degree,p);
	for(int i=degree;i>=1;i--) p[i]=Real(i)*x*p[i-1];
	p[0]=0;
}

template <class T>
void Basis_DDx(const T x, int degree, T p[])
{
	Basis(x,degree,p);
	for(int i=degree;i>=2;i--) p[i]=Real(i*(i-1))*x*p[i-2];
	p[0]=0; if(degree<1) return;
	p[1]=0;
}

//The kth derivative of the basis
template <class T>
void Basis_Dkx(const T x, int degree, T p[], int k)
{
	Basis(x,degree,p);
	for(int i=degree;i>=k;i--) p[i]=Real(FactorialTruncated(i,k))*x*p[i-k];
	int end=(k<degree?k:degree);
	for(int i=0;i<end;i++) p[i]=0;
}


} //namespace Polynomial
} //namespace Math
#endif
