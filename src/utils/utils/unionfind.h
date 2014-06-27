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

#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <vector>

/** @ingroup Utils
 * @brief From an indexed set of elments, allows fast unioning 
 * of the elements into sets.
 *
 * From a set X={x1,x2,...} (such that |X| = entries), the sets
 * are initialized to Si = {xi}.  Larger sets are created by
 * unioning the set of xi with that of xj using the Union() method.
 *
 * The sets are referenced using unique integer identifiers
 * from 0 to entries-1.
 */
class UnionFind
{
public:
  UnionFind(int entries);
  ///Resize X to size entries, sets all sets Si={xi}
  void Initialize(const int entries);
  ///Increment the size of X by 1, sets Sn={xn}
  int AddEntry();
  ///Returns the id of the set to which xi belongs 
  int FindSet(const int i);  
  ///Unions the sets to which xi and xj belong
  int Union(const int i,const int j);
  ///Returns the sets to which the items belong, such that sets[i] = set(xi)
  void EnumerateSets(std::vector<int>& sets);

private:
  std::vector<int> parents;

  inline bool IsRoot(const int i) const { return parents[i]==-1; }
  int FindRoot(const int i) const;
  void PathCompress(const int i,const int root);
  void CompressAll();
};
#endif

