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

#ifndef ITERATOR_UTILS_H
#define ITERATOR_UTILS_H

/** @file utils/iteratorutils.h
 * @ingroup Utils
 * @brief Provides basic utilities to increment/decrement an 
 * iterator by a certain amount.
 */

template <class T>
inline T increment(const T& a,int n=1)
{
  T x=a;
  for(int i=0;i<n;i++) ++x;
  return x;
}

template <class T>
inline T decrement(const T& a,int n=1)
{
  T x=a;
  for(int i=0;i<n;i++) --x;
  return x;
}

//returns a-b, the # of times you'd need to increment b to get to a
template <class T>
inline int iterator_diff(const T& a,const T& b)
{
  T x=b;
  int n=0;
  while(x!=a) { ++x; ++n; }
  return n;
}

#endif
