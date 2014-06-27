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

#ifndef ARRAY_MAPPING_H
#define ARRAY_MAPPING_H

#include <vector>

/* ArrayMapping
 * Maps array indices i in [0,max) to either:
 *  - (i+offset) 
 *  - (mapping[i]) (if the mapping is defined)
 * Also, Map : x->mapx is defined so mapx[Map(i)] <- x[i]
 */
struct ArrayMapping
{
  inline ArrayMapping();

  inline bool IsOffset() const { return mapping.empty(); }
  inline int Size() const { return (IsOffset() ?  max : (int)mapping.size()); } 
  inline void SetOffset(int offset,int max);

  inline int Map(int i) const;
  inline int InvMap(int imap) const;
  template <class A>
  inline void Map(const A& x,A& mapx) const;
  template <class A>
  inline void InvMap(const A& mapx,A& x) const;  

  std::vector<int> mapping;    
  int max;
  int offset;
};


inline ArrayMapping::ArrayMapping()
  :max(0),offset(0)
{}

inline void ArrayMapping::SetOffset(int _offset,int _max) 
{
  mapping.clear(); 
  offset=_offset; 
  max=_max; 
}

inline int ArrayMapping::Map(int i) const 
{
  return (IsOffset() ? i+offset : mapping[i]);
}

inline int ArrayMapping::InvMap(int imap) const 
{
  if(IsOffset()) 
    return imap-offset;
  else {
    for(size_t i=0;i<mapping.size();i++)
      if(mapping[i] == imap) return (int)i;
    abort();
    return -1;
  }
}

template <class A>
inline void ArrayMapping::Map(const A& x,A& mapx) const 
{
  if(mapping.empty()) {
    for(int i=0;i<max;i++)
      mapx[i+offset] = x[i];
  }
  else {
    for(int i=0;i<(int)mapping.size();i++)
      mapx[mapping[i]] = x[i];
  }
}

template <class A>
void ArrayMapping::InvMap(const A& mapx,A& x) const 
{
  if(IsOffset()) {
    for(int i=0;i<max;i++)
      x[i] = mapx[i+offset];
  }
  else {
    for(int i=0;i<(int)mapping.size();i++)
      x[i] = mapx[mapping[i]];
  }
}
  


#endif
