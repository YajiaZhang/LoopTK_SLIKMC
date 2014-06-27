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

#ifndef MATH_REAL_FUNCTION_H
#define MATH_REAL_FUNCTION_H

#include "function.h"

namespace Math { 

//y=a*t+b
class LinearFunction : public RealFunction
{
public:
  LinearFunction(Real _a=One,Real _b=Zero) : a(_a),b(_b) {}
  virtual std::string Label() const { return "<a*t+b>"; }
  virtual Real Eval(Real t) { return a*t+b; }
  virtual Real Deriv(Real t) { return a; }
  virtual Real Deriv2(Real t) { return 0; }  

  Real a,b;
};

//y=1/t
class InverseFunction : public RealFunction
{
public:
  virtual std::string Label() const { return "<1/t>"; }
  virtual Real Eval(Real t) { return One/t; }
  virtual Real Deriv(Real t) { return -One/(t*t); }
  virtual Real Deriv2(Real t) { return Two/(t*t*t); }
};

//y=f(g(t))
class ComposeFunction : public RealFunction
{
public:
  ComposeFunction(RealFunction* _f,RealFunction* _g) : f(_f), g(_g) {}
  virtual std::string Label() const;
  virtual void PreEval(Real t);
  virtual Real Eval(Real t);
  virtual Real Deriv(Real t);
  virtual Real Deriv2(Real t);
  
  RealFunction *f,*g;
  Real gt;
};

//y=f(t)+g(t)
class AddFunction : public RealFunction
{
public:
  AddFunction(RealFunction* _f,RealFunction* _g) : f(_f), g(_g) {}
  virtual std::string Label() const;
  virtual void PreEval(Real t);
  virtual Real Eval(Real t);
  virtual Real Deriv(Real t);
  virtual Real Deriv2(Real t);
  
  RealFunction *f,*g;
};

//y=f(t)*g(t)
class MulFunction : public RealFunction
{
public:
  MulFunction(RealFunction* _f,RealFunction* _g) : f(_f), g(_g) {}
  virtual std::string Label() const;
  virtual void PreEval(Real t);
  virtual Real Eval(Real t);
  virtual Real Deriv(Real t);
  virtual Real Deriv2(Real t);
  
  RealFunction *f,*g;
  Real ft,gt;
};

} //namespace Math

#endif
