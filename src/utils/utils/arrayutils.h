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

#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

#include <algorithm>
#include <functional>
#if defined (__GNUC__) && (__GNUC__ > 2)
#include <ext/algorithm>
namespace std {
  using __gnu_cxx::copy_n;
}
#endif
#include <assert.h>

namespace ArrayUtils {

template <typename type>
inline void copy(const type& a, type* out, int n)
{
  std::fill(out,out+n,a);
}

template <typename type>
inline void copy(const type* a, type* out, int n)
{
  std::copy_n(a,n,out);
}

template <typename type,typename ftype>
inline void foreach(type* a, ftype f,int n)
{
  std::for_each(a,a+n,f);
}

template <typename type>
inline void reverse (type *a, int n)
{
  std::reverse(a,a+n);
}


template <typename t1,typename t2,typename ftype>
inline void transform(const t1* a, t2* out, ftype f,int n)
{
  std::transform(a,a+n,out,f);
}

template <typename t1,typename t2,typename t3,typename ftype>
inline void binary_transform(const t1* a, const t2* b, t3* out, ftype f,int n)
{
  std::transform(a,a+n,b,out,f);
}

template <typename type>
inline void add(const type* a, const type* b, type* out, int n)
{
  binary_transform(a,b,out,std::plus<type>(),n);
}

template <typename type>
inline void sub(const type* a, const type* b, type* out, int n)
{
  binary_transform(a,b,out,std::minus<type>(),n);
}

template <typename type>
inline void mul(const type* a, const type* b, type* out, int n)
{
  binary_transform(a,b,out,std::multiplies<type>(),n);
}

template <typename type>
inline void div(const type* a, const type* b, type* out, int n)
{
  binary_transform(a,b,out,std::divides<type>(),n);
}

///returns the nth largest element in the array a
template <typename type>
inline type nth_element (const std::vector<type>& S, size_t n)
{
  assert(n < S.size());
  size_t i=rand()%S.size();
  const type& m=S[i];
  std::vector<type>S1,S2;
  S1.reserve(n);
  S2.reserve(n);
  for(i=0;i<S.size();i++) {
    if(S[i] < m) S1.push_back(S[i]);
    else if(m < S[i]) S2.push_back(S[i]);
  }
  if(S1.size() > n) return nth_element(S1,n,std::less<type>());
  else if(S.size()-S2.size()>=n) return m;
  else return nth_element(S2,n-(S.size()-S2.size()),std::less<type>());
}

///returns the nth largest element in the array a
template <typename type,typename fless>
inline type nth_element (const std::vector<type>& S, size_t n, fless less)
{
  assert(n < S.size());
  size_t i=rand()%S.size();
  const type& m=S[i];
  std::vector<type>S1,S2;
  S1.reserve(n);
  S2.reserve(n);
  for(i=0;i<S.size();i++) {
    if(less(S[i],m)) S1.push_back(S[i]);
    else if(less(m,S[i])) S2.push_back(S[i]);
  }
  if(S1.size() > n) return nth_element(S1,n,less);
  else if(S.size()-S2.size()>=n) return m;
  else return nth_element(S2,n-(S.size()-S2.size()),less);
}

template <typename type>
inline bool is_sorted(type* a, int n)
{
  for(int i=1;i<n;i++)
    if(a[i]<a[i-1]) return false;
  return true;
  //return std::is_sorted(a,a+n);
}

template <typename type,typename fless>
inline bool is_sorted(type* a, int n, fless f)
{
  for(int i=1;i<n;i++)
    if(f(a[i],a[i-1])) return false;
  return true;
  //return std::is_sorted(a,a+n,f);
}

template <typename type>
inline void quicksort(type* a, int p, int r)
{
  if(p < r) {
    type x = a[p];
    type temp;
    int i = p;
    int j = p + 1;
    while (j <= r) {
      if (a[j] < x) {
        i++;
        temp=a[j];a[j]=a[i];a[i]=temp;
      }
      j++;
    }
    temp=a[p];a[p]=a[i];a[i]=temp;

    quicksort(a,p,i-1);
    quicksort(a,i+1,r);
  }
}

template <typename type,typename fless>
inline void quicksort(type* a, int p, int r,fless f)
{
  if(p < r) {
    type x = a[p];
    type temp;
    int i = p;
    int j = p + 1;
    while (j <= r) {
      if (f(a[j],x)) {
        i++;
        temp=a[j];a[j]=a[i];a[i]=temp;
      }
      j++;
    }
    temp=a[p];a[p]=a[i];a[i]=temp;

    quicksort(a,p,i-1,f);
    quicksort(a,i+1,r,f);
  }
}

template <typename type>
inline void sort(type* a, int n)
{
  quicksort(a,0,n-1);
}

template <typename type,typename fless>
inline void sort(type* a, int n, fless f)
{
  quicksort(a,0,n-1,f);
}

} //namespace ArrayUtils

#endif
