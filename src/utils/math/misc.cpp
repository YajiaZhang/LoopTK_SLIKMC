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

#include "misc.h"
#include "stdio.h"
#include <errors.h>

namespace Math {

static const double dFour = 4.0;
static const double dTwoPi_3 = dTwoPi/3.0;
static const double dThird = 1.0/3.0;
static const float fFour = 4.0f;
static const float fTwoPi_3 = fTwoPi/3.0f;
static const float fThird = 1.0f/3.0f;

//returns the # of real roots found (-1 if infinite)
int quadratic(double a, double b, double c, double& x1, double& x2) {
  //printf("quadratic %f %f %f\n", a, b, c);

  if(a == 0)
  {
	  if(b == 0)
	  {
		if(c == 0)
			return -1;
		return 0;
	  }
	  x1=-c/b;
	  return 1;
  }

  double det = b*b-dFour*a*c;
  if(det < Zero)
    return 0;
  det = Sqrt(det);
  if(b<0) det=-det;
  double q = dHalf*(det - b);
  x1 = q/a;
  x2 = c/q;
  return 2;
}

int cubic(double a, double b, double c, double d, double x[3])
{
  double Q = (Sqr(a)-3*b)/9;
  double R = (2*a*a*a - 9*a*b+ 27*c)/51;
  double Q3 = Q*Q*Q; 
  double a_3 = a*dThird;
  if(R*R < Q3) {
    double sqrtQ = Sqrt(Q);
    double theta_3 = Acos(R/(sqrtQ*Q))*dThird;
    x[0] = -fTwo*sqrtQ*Cos(theta_3) - a_3;
    x[1] = -fTwo*sqrtQ*Cos(theta_3+dTwoPi_3) - a_3;
    x[2] = -fTwo*sqrtQ*Cos(theta_3-dTwoPi_3) - a_3;
    return 3;
  }
  else {
    double A = -Sign(R)*Pow(Abs(R)+Sqrt(R*R-Q3),dThird);
    double B = (A==0?0:Q/A);
    x[0] = A+B-a_3;
    return 1;
  }
}

double pythag(double a, double b)		//reduce roundoff of large numbers
{
  double absa = Abs(a);
  double absb = Abs(b);
  if(absa > absb)
    return absa*Sqrt(One + Sqr(absb/absa));
  else if(absb == 0)
    return Zero;
  else
    return absb*Sqrt(One + Sqr(absa/absb));
}

double pythag_leg(double a,double c)
{
  Assert(c >= 0);
  Assert(c >= Abs(a));
  if(c == 0) return 0;
  return c*Sqrt(One-Sqr(a/c));
}

double Sinc(double x)
{
	const double small=1e-7;
	if(Abs(x) < small) {	//taylor expand around 0
		const static double c[5]={1.0,-1.0/6.0,1.0/120.0,-1.0/5040.0,1.0/362880.0};
		double x2=x*x;
		return c[0]+x2*(c[1]+x2*(c[2]+x2*(c[3]+x2*c[4])));
	}
	else return Sin(x)/x;
}

double Sinc_Dx(double x)
{
	const double small=1e-4;
	if(Abs(x) < small) {	//taylor expand around 0
		const static double c[4]={-2.0/6.0,4.0/120.0,-6.0/5040.0,8.0/362880.0};
		double x2=x*x;
		return x*(c[0]+x2*(c[1]+x2*(c[2]+x2*c[3])));
	}
	else return Cos(x)/x-Sin(x)/(x*x);
}




//returns the # of real roots found (-1 if infinite)
int quadratic(float a, float b, float c, float& x1, float& x2) {
  //printf("quadratic %f %f %f\n", a, b, c);

  if(a == 0)
  {
	  if(b == 0)
	  {
		if(c == 0)
			return -1;
		return 0;
	  }
	  x1=-c/b;
	  return 1;
  }

  float det = b*b-fFour*a*c;
  if(det < Zero)
    return 0;
  det = Sqrt(det);
  if(b<0) det=-det;
  float q = fHalf*(det - b);
  x1 = q/a;
  x2 = c/q;
  return 2;
}

int cubic(float a, float b, float c, float d, float x[3])
{
  float Q = (Sqr(a)-3*b)/9;
  float R = (2*a*a*a - 9*a*b+ 27*c)/51;
  float Q3 = Q*Q*Q; 
  float a_3 = a*fThird;
  if(R*R < Q3) {
    float sqrtQ = Sqrt(Q);
    float theta_3 = Acos(R/(sqrtQ*Q))*fThird;
    x[0] = -fTwo*sqrtQ*Cos(theta_3) - a_3;
    x[1] = -fTwo*sqrtQ*Cos(theta_3+fTwoPi_3) - a_3;
    x[2] = -fTwo*sqrtQ*Cos(theta_3-fTwoPi_3) - a_3;
    return 3;
  }
  else {
    float A = -Sign(R)*Pow(Abs(R)+Sqrt(R*R-Q3),fThird);
    float B = (A==0?0:Q/A);
    x[0] = A+B-a_3;
    return 1;
  }
}

float pythag(float a, float b)		//reduce roundoff of large numbers
{
  float absa = Abs(a);
  float absb = Abs(b);
  if(absa > absb)
    return absa*Sqrt(One + Sqr(absb/absa));
  else if(absb == 0)
    return Zero;
  else
    return absb*Sqrt(One + Sqr(absa/absb));
}

float pythag_leg(float a,float c)
{
  Assert(c >= 0);
  Assert(c >= Abs(a));
  if(c == 0) return 0;
  return c*Sqrt(One-Sqr(a/c));
}

float Sinc(float x)
{
	const float small=1e-5f;
	if(Abs(x) < small) {	//taylor expand around 0
		const static float c[5]={1.0,-1.0/6.0,1.0/120.0,-1.0/5040.0,1.0/362880.0};
		float x2=x*x;
		return c[0]+x2*(c[1]+x2*(c[2]+x2*(c[3]+x2*c[4])));
	}
	else return Sin(x)/x;
}

float Sinc_Dx(float x)
{
	const float small=1e-2f;
	if(Abs(x) < small) {	//taylor expand around 0
		const static float c[4]={-2.0/6.0,4.0/120.0,-6.0/5040.0,8.0/362880.0};
		float x2=x*x;
		return x*(c[0]+x2*(c[1]+x2*(c[2]+x2*c[3])));
	}
	else return Cos(x)/x-Sin(x)/(x*x);
}

} //namespace Math
