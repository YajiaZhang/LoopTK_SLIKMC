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

#include "unionfind.h"

UnionFind::UnionFind(int entries)
{
    parents.resize(entries,-1);
}

void UnionFind::Initialize(const int entries)
{
	parents.clear();
	parents.resize(entries,-1);
}

int UnionFind::AddEntry()
{
	parents.push_back(-1);
	return parents.size()-1;
}

int UnionFind::FindSet(const int i)
{
	int root=FindRoot(i);
	PathCompress(i,root);
	return root;
}

//unions the sets to which i and j belong
int UnionFind::Union(const int i,const int j)
{
	int root_i=FindSet(i),root_j=FindRoot(j);
	PathCompress(j,root_i);
	if(root_i!=root_j)
		parents[root_j]=root_i;
	return root_i;
}

//since on a merge, we don't go through all the children of a set
//to change their roots, we have to do some path compression
void UnionFind::CompressAll()
{
  for(int i=0;i<(int)parents.size();i++) PathCompress(i,FindRoot(i));
}

//enumerates the sets to which the items belong
void UnionFind::EnumerateSets(std::vector<int>& sets)
{
  CompressAll();
  sets.resize(parents.size());
  for(int i=0;i<(int)parents.size();i++) sets[i]=FindRoot(i);
}

int UnionFind::FindRoot(const int i) const
{
  int j=i;
  while(!IsRoot(j)) { j=parents[j]; }
  return j;
}

void UnionFind::PathCompress(const int i,const int root)
{
  int j=i,k;
  while(!IsRoot(j)) { k=parents[j]; parents[j]=root; j=k;}
}
