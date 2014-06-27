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

#ifndef BASIC_UTILS_H
#define BASIC_UTILS_H

#include <stdlib.h>
#include <stdio.h>

/** @defgroup Standard 
 * @brief Commonly used classes/routines.
 * @{
 */

/** @file utils.h
 * @brief Utilities commonly used throughout a program.  
 *
 * Many of these are here just to provide a cross-platform interface
 * to non-standard build environments (e.g. windows)
 */

/** @def SafeDelete
 * @brief Delete a non-NULL pointer and set it to NULL.
 */

/** @def SafeArrayDelete
 * @brief Delete a non-NULL pointer to an array and set it to NULL.
 */

/** @def SafeDeleteProc
 * @brief Delete a non-NULL pointer, using the provided function,
 *  and set it to NULL.
 */

/** @def ARRAYSIZE
 * @brief Returns the size of the array x[].
 *
 * Returns n for constant-sized arrays declared as type x[n].  Does not
 * work for dynamic arrays allocated on the heap.
 */

#ifdef  _MSC_VER
template<class type>
inline type Max(const type& a,const type& b) { return (a>b?a:b); }
template<class type>
inline type Min(const type& a,const type& b) { return (a<b?a:b); }
template<class type>
inline void Swap(type& a,type& b) { type temp=a; a=b; b=temp; }
#else
#include <algorithm>
template<class type>
inline type Max(const type& a,const type& b) { return std::max(a,b); }
template<class type>
inline type Min(const type& a,const type& b) { return std::min(a,b); }
template<class type>
inline void Swap(type& a,type& b) { std::swap(a,b); }
#endif //_MSC_VER

template <class type>
inline type Max(const type& a,const type& b,const type& c)
{
  return Max(Max(a,b),c);
}

template <class type>
inline type Max(const type& a,const type& b,const type& c,const type& d)
{
  return Max(Max(a,b),Max(c,d));
}

template <class type>
inline type Min(const type& a,const type& b,const type& c)
{
  return Min(Min(a,b),c);
}

template <class type>
inline type Min(const type& a,const type& b,const type& c,const type& d)
{
  return Min(Min(a,b),Min(c,d));
}



inline void toggle(bool& bit)
{
  bit=!bit;
}

#define SafeDelete(x) { if (x) delete x; x=NULL; }
#define SafeArrayDelete(x) { if (x) delete [] x; x=NULL; }
#define SafeDeleteProc(x,proc) { if (x) proc(x); x=NULL; }

#define ARRAYSIZE(x) (sizeof(x)/sizeof(x[0]))

/*@}*/

#endif

