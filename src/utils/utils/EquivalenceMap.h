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

#ifndef UTILS_EQUIVALENCE_MAP_H
#define UTILS_EQUIVALENCE_MAP_H

#include "unionfind.h"
#include <map>

//for a mapping f:i->fwd[i], this gets the inverse multi-mapping
//g:j->{i|fwd[i]=j}
inline void InverseMapping(const std::vector<int>& fwd,std::vector<vector<int> >& bwd) {
  bwd.clear();
  bwd.reserve(fwd.size());
  int N=*max_element(fwd.begin(),fwd.end());
  bwd.resize(N+1);
  for(size_t i=0;i<fwd.size();i++) {
    assert(0<=fwd[i]&&fwd[i]<=N);
    bwd[fwd[i]].push_back(i);
  }
  size_t n=0;
  for(size_t i=0;i<bwd.size();i++) {
    n+=bwd[i].size();
  }
  if(n != fwd.size()) {
    cout<<"InverseMapping error!"<<endl;
    cout<<"inserted only "<<n<<" of "<<fwd.size()<<" items!"<<endl;
    abort();
  }
}

//each entry in eq is a vector containing all
//indices such that x[i]=x[j] for some i,j
//n^2 algorithm
template <class T,class EqFn>
void EquivalenceMap(const std::vector<T>& x,std::vector<vector<int> >& eq,EqFn& Eq)
{
  int n=(int)x.size();
  UnionFind uf(n);
  for(int i=0;i<n;i++) {
    int seti=i;
    for(int j=0;j!=i;j++) {
      int setj = uf.FindSet(j);
      if(seti != setj) {
	if(Eq(x[i],x[j])) {
	  //merge i->j
	  uf.Union(j,i);
	  assert(uf.FindSet(j) == setj);
	  assert(uf.FindSet(i) == setj);
	  seti = setj;
	}
      }
    }
  }
  vector<int> sets;
  uf.EnumerateSets(sets);
  InverseMapping(sets,eq);
  //remove empty sets
  vector<vector<int> >::iterator it,prev;
  int i=0;
  int num=0;
  for(it=eq.begin();it!=eq.end();it++,i++) {
    num += (int)it->size();
    if(it->empty()) {
      prev = it; prev--;
      eq.erase(it);
      it = prev;
    }
  }
  assert(num==n);
}


#endif
