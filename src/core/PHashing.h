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

#ifndef PHASHING_H
#define PHASHING_H

#include <string>
#include <vector>
#include <ext/hash_set>
#include <ext/hash_map>
#include <math3d/primitives.h>

using namespace std;
using namespace __gnu_cxx;

typedef pair<string,string> StringPair;

class PAtom;
class PBond;

struct vectorHash  {
  public:
    size_t operator()(const Vector3 &v) const
    {
      hash<int> intHash;
      int x = int(v.x) + 101, y = int(v.y) + 201, z = int(v.z) + 301;
      return intHash(x * y * z);
    }
};
                                                                                                                                                             
struct vectorEq {
  public:
    bool operator()(const Vector3 &v1, const Vector3 &v2) const
    {
      return(v1 == v2); // use predefined Vector3.operator ==
    }
};

struct vectorPointerHash {
  public:
    size_t operator()(const Vector3 *v) const
    {
      hash<size_t> h;
      size_t vectorPtr = (size_t) v;
      return h(vectorPtr);
    }
};

struct vectorPointerEq {
  public:
  bool operator()(const Vector3 *v1, const Vector3 *v2) const
  {
    return (v1 == v2);
  }
};

class stringHasher {
 public:
  size_t operator()(string const &str) const
    {
      hash<char const *> h;
      return (h(str.c_str()));
    }
};

struct eqString
{
  bool operator()(string s1, string s2) const
  {
    return(s1 == s2);
  }
};


class stringPairHasher {
 public:
  size_t operator()(StringPair const &s) const
    {
      stringHasher h;
      return (h(s.first+s.second));
    }
};

struct eqStringPairExclusive 
{
  bool operator()(const StringPair &s1, const StringPair &s2) const
  {
    return s1.first == s2.first && s1.second == s2.second;
  }

};

struct eqStringPairEither
{
  bool operator()(const StringPair &s1, const StringPair &s2) const
  {
    if (s1.first == s2.first) return s1.second == s2.second;
    else if (s1.first == s2.second) return s1.second==s2.first;
    else return false;
  }

};

/*
 * Hashing functors.
 */

struct pairAtomHash {
  public:
    size_t operator()(const pair<PAtom *,PAtom *> a) const {
      return (size_t) a.first + (size_t) a.second;
    }
};

struct pairAtomEq {
  public:
    bool operator()(const pair<PAtom *,PAtom *> a1, const pair<PAtom *,PAtom *> a2) const {
      return (a1.first == a2.first) && a1.second == a2.second ||
              a1.first==a2.second && a1.second == a2.first;
    }
};

struct atomHash {
  public:
    size_t operator()(const PAtom *a) const {
      hash<size_t> h;
      size_t atomPtr = (size_t) a;
      return h(atomPtr);
    }
};

struct atomEq {
  public:
    bool operator()(const PAtom *a1, const PAtom *a2) const {
      return (a1 == a2);
    }
};

/*
 * Data structures to store atom information,
 * e.g. sets of atoms, positions of atoms, etc.
 */

typedef hash_set<PAtom*, atomHash, atomEq> AtomSet;
typedef hash_set<const PAtom *, atomHash, atomEq> ConstAtomSet;
typedef hash_set<const Vector3 *, vectorPointerHash, vectorPointerEq> ConstVectorSet;

typedef hash_set<pair<PAtom *, PAtom *>, pairAtomHash, pairAtomEq> AtomCollisions;
typedef hash_map<pair<PAtom *, PAtom *>, int, pairAtomHash, pairAtomEq> AtomSeparations;
typedef hash_map<Vector3, AtomSet, vectorHash, vectorEq> CollisionMap;
typedef hash_map<Vector3, bool, vectorHash, vectorEq> OccupancyMap;


//MACROS for making string and string pair hashmaps

#define HASH_MAP_STR(t) hash_map<string, t, stringHasher, eqString>
#define HASH_MAP_STRPAIR_EX(t) hash_map<StringPair, t, stringPairHasher, eqStringPairExclusive>
#define HASH_MAP_STRPAIR_OR(t) hash_map<StringPair, t, stringPairHasher, eqStringPairEither>


typedef HASH_MAP_STR(vector<PBond *>) DOF_Cache;
typedef HASH_MAP_STR(vector<PAtom *>) Atom_Cache;
typedef HASH_MAP_STR(Vector3) AtomPositions;
typedef HASH_MAP_STR(string) StringMap;
typedef HASH_MAP_STRPAIR_OR(PBond *) ID_To_DOF_Map;


#endif
