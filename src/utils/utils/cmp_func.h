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

#ifndef UTILS_CMP_FUNC_H
#define UTILS_CMP_FUNC_H

#include <functional>

namespace std {

/** @ingroup Utils
 * @brief A templated function class like strcmp: 
 * returns -1 if a<b, 0 if a==b, 1 if a>b.
 */
template <class T,class Less=less<T> >
struct cmp_func
{
  int operator()(const T& a,const T& b) {
    if(lessFunc(a,b)) return -1;
    else if(lessFunc(b,a)) return 1;
    return 0;
  }
  Less lessFunc;
};

} //namespace std

#endif
