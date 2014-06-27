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

#ifndef UTILS_PERMUTATION_H
#define UTILS_PERMUTATION_H

#include <vector>

/** @file utils/permutation.h
 * @brief Various permutation utilities.
 *
 * A permutation is defined by an array v (or a vector) of length n.
 * y[i] = x[v[i]] is the action of the permutation.
 */

inline void IdentityPermutation(int v[],int n)
{
  for(int i=0;i<n;i++) v[i]=i;
}

inline void RandomlyPermute(int v[],int n)
{
  for(int i=0;i<n;i++) {
    int k=i+rand()%(n-i);
    std::swap(v[i],v[k]);
  }
}

inline void RandomPermutation(int v[],int n)
{
  IdentityPermutation(v,n);
  RandomlyPermute(v,n);
}

inline void IdentityPermutation(std::vector<int>& v)
{
  int n=(int)v.size();
  for(int i=0;i<n;i++) v[i]=i;
}

inline void RandomlyPermute(std::vector<int>& v)
{
  int n=(int)v.size();
  for(int i=0;i<n;i++) {
    int k=i+rand()%(n-i);
    std::swap(v[i],v[k]);
  }
}

inline void RandomPermutation(std::vector<int>& v)
{
  IdentityPermutation(v);
  RandomlyPermute(v);
}



///Gets the next lexicographically ordered combination.
///Returns 1 if it reaches the end.
inline int NextCombination(int v[],int n)
{
  for(int i=n-1;i>=0;i--) {
    v[i]++;
    if(v[i] >= n) v[i]=0;
    else return 0;
  }
  return 1;
}

///Gets the next lexicographically ordered combination.
///Returns 1 if it reaches the end
inline int NextCombination(std::vector<int>& v)
{
  int n=(int)v.size();
  for(int i=n-1;i>=0;i--) {
    v[i]++;
    if(v[i] >= n) v[i]=0;
    else return 0;
  }
  return 1;
}


#endif
