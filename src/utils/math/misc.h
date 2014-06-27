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

#ifndef MATH_MISC_H
#define MATH_MISC_H

#include "math.h"

namespace Math {

  /** @addtogroup Math */
  /*@{*/

int quadratic(double a, double b, double c, double& x1, double& x2);
int cubic(double a, double b, double c, double d, double x[3]);
///returns c where a^2+b^2=c^2
double pythag(double a, double b);
///returns b where a^2+b^2=c^2
double pythag_leg(double a,double c);
double Sinc(double x);
double Sinc_Dx(double x);


///returns -1 if x<eps, 0 if |x|<=eps, 1 if x>eps
inline int FuzzySign(double x,double eps=dEpsilon) { return (x>eps?1:(x<-eps?-1:0)); }
///Returns true of a and be have different signs, treating 0 as 
///a positive number
inline bool OpposingSigns_pos(double a,double b) { return (a<dZero) != (b<dZero); }
///Same as above, but treats 0 as a negative number
inline bool OpposingSigns_neg(double a,double b) { return (a<=dZero) != (b<=dZero); }


int quadratic(float a, float b, float c, float& x1, float& x2);
int cubic(float a, float b, float c, float d, float x[3]);
float pythag(float a, float b);
float pythag_leg(float a,float c);
float Sinc(float x);
float Sinc_Dx(float x);

inline int FuzzySign(float x,float eps=fEpsilon) { return (x>eps?1:(x<-eps?-1:0)); }
inline bool OpposingSigns_pos(float a,float b) { return (a<fZero) != (b<fZero); }
inline bool OpposingSigns_neg(float a,float b) { return (a<=fZero) != (b<=fZero); }



///calculates x^i where i is an integer (in O(log i) time)
template <class T>
T IntegerPower(const T x, int i)
{
  if(i<=0) return 1;
  if(i&1)	//odd
    { T y=IntegerPower(x,(i-1)>>1); return x*y*y; }
  else //even
    { T y=IntegerPower(x,i>>1); return y*y; }
}

///returns n! in O(n) time
inline int Factorial(int n)
{
  //assert(n <= 12); //otherwise, it'll cause overflow with 32-bit arith
  int x=1;
  for(int i=1;i<=n;i++) x*=i;
  return x;
}
 
///returns n!/(n-k)! in O(k) time
inline int FactorialTruncated(int n,int k)
{
  int fact=1;
  for(int j=0;j<k;j++) fact*=(n-j);
  return fact;
}

///returns `n choose k' = n!/k!(n-k)! in O(k) time
inline int Choose(int n,int k)
{
  if(2*k > n) return FactorialTruncated(n,n-k)/Factorial(n-k);
  else return FactorialTruncated(n,k)/Factorial(k);
}

  /*@}*/
} //namespace Math

#endif

