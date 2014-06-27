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

#include "BlockVector.h"
#include <utils/stringutils.h>
#include <errors.h>
using namespace Math;
using namespace std;

BlockVector::BlockVector()
{}

BlockVector::BlockVector(int numBlocks)
{
  resize(numBlocks);
}

BlockVector::BlockVector(int numBlocks,int vecSize)
{
  resize(numBlocks,vecSize);
}

BlockVector::BlockVector(int numBlocks,int vecSize,Real initVal)
{
  resize(numBlocks,vecSize,initVal);
}

BlockVector::BlockVector(int numBlocks,const Vector& initVal)
  :BaseT(numBlocks,initVal)
{}

void BlockVector::resize(int numBlocks)
{
  BaseT::resize(numBlocks);
}

void BlockVector::resize(int numBlocks,int vecSize)
{
  BaseT::resize(numBlocks);
  for(int i=0;i<numBlocks;i++)
    operator[](i).resize(vecSize);
}

void BlockVector::resize(int numBlocks,int vecSize,Real initVal)
{
  BaseT::resize(numBlocks);
  for(int i=0;i<numBlocks;i++)
    operator[](i).resize(vecSize,initVal);
}

void BlockVector::resizeSimilar(const BlockVector& v)
{
  BaseT::resize(v.size());
  for(size_t i=0;i<size();i++)
    operator[](i).resize(v[i].n);
}

void BlockVector::copy(const BlockVector& a)
{
  BaseT::operator=(a);
}

void BlockVector::swap(BlockVector& a)
{
  std::swap((BaseT&)*this,(BaseT&)a);
}

void BlockVector::add(const BlockVector& a,const BlockVector& b)
{
  Assert(a.size()==b.size());
  resize(a.size());
  for(size_t i=0;i<a.size();i++)
    operator[](i).add(a[i],b[i]);
}

void BlockVector::sub(const BlockVector& a,const BlockVector& b)
{
  Assert(a.size()==b.size());
  resize(a.size());
  for(size_t i=0;i<a.size();i++)
    operator[](i).sub(a[i],b[i]);
}

void BlockVector::mul(const BlockVector& a,Real c)
{
  resize(a.size());
  for(size_t i=0;i<a.size();i++)
    operator[](i).mul(a[i],c);
}

void BlockVector::div(const BlockVector& a,Real c)
{
  resize(a.size());
  for(size_t i=0;i<a.size();i++)
    operator[](i).div(a[i],c);
}

void BlockVector::inc(const BlockVector& a)
{
  Assert(hasDims(a));
  for(size_t i=0;i<a.size();i++)
    operator[](i) += a[i];
}

void BlockVector::dec(const BlockVector& a)
{
  Assert(hasDims(a));
  for(size_t i=0;i<a.size();i++)
    operator[](i) -= a[i];
}

void BlockVector::madd(const BlockVector& a,Real c)
{
  Assert(hasDims(a));
  for(size_t i=0;i<a.size();i++)
    operator[](i).madd(a[i],c);
}

void BlockVector::setZero()
{
  set(Zero);
}

void BlockVector::set(Real c)
{
  for(size_t i=0;i<size();i++)
    operator[](i).set(c);
}

void BlockVector::set(const BlockVector& v)
{
  BaseT::operator =(v);
}

void BlockVector::setNegative(const BlockVector& v)
{
  resize(v.size());
  for(size_t i=0;i<size();i++)
    operator[](i).setNegative(v[i]);
}

void BlockVector::inplaceNegative()
{
  for(size_t i=0;i<size();i++)
    operator[](i).inplaceNegative();
}

void BlockVector::inplaceMul(Real c)
{
  for(size_t i=0;i<size();i++)
    operator[](i) *= c;
}

void BlockVector::inplaceDiv(Real c)
{
  for(size_t i=0;i<size();i++)
    operator[](i) /= c;
}

Real BlockVector::dot(const BlockVector& v) const
{
  Assert(hasDims(v));
  Real sum=Zero;
  for(size_t i=0;i<size();i++)
    sum += operator[](i).dot(v[i]);
  return sum;
}

Real BlockVector::normSquared() const
{
  Real sum=Zero;
  for(size_t i=0;i<size();i++)
    sum += operator[](i).normSquared();
  return sum;
}

bool BlockVector::hasDims(const BlockVector& v) const
{
  if(size()!=v.size()) return false;
  for(size_t i=0;i<size();i++)
    if(operator[](i).n != v[i].n) return false;
  return true;
}

bool BlockVector::hasDims(int numBlocks,int vecSize) const
{
  if((int)size()!=numBlocks) return false;
  for(size_t i=0;i<size();i++)
    if(operator[](i).n != vecSize) return false;
  return true;
}



