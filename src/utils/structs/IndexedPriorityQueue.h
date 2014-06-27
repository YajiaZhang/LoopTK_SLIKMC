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

#ifndef INDEXED_PRIORITY_QUEUE_H
#define INDEXED_PRIORITY_QUEUE_H

#include <set>
#include <map>

namespace std {

/** @brief Contains a priority queue with an associated index, allowing
 * updates of the priority value.
 *
 * An index is provided with each priority queue entry, allowing a priority
 * value
 * to be located from an index in logarithmic time.  
 * 
 * IT is the type of the index.  It is assumed that the < operator has been defined.
 * PT is the type of the priority value.  It is assumed that the < operator has
 * been defined.
 *
 * Most members operate like STL's priority_queue<pair<PT,IT> >.  Be sure
 * to remember that the returned pair has PT as its first type, and IT
 * as its second!
 *
 * The new members are as follows:
 * - iterator insert(const IT& index,const PT& p): inserts a priority queue
 *   entry p, indexed by index.  Returns the inserted iterator.
 * - iterator find(const IT& index): returns the iterator corresponding
 *   to index.
 * - void refresh(const IT& index,const PT& p): sets the priority value
 *   of index to p.
 * - bool is_valid() const: Returns true if the internal structures are valid.
 *
 * One index value can only map to one priority value (that is, 
 * insert(2,0.5) and insert(2,0.7) is an invalid operation).
 *
 * Undefined behavior will result if any iterator is modified.
 */
template <class IT,class PT>
class IndexedPriorityQueue
{
 public:
  typedef pair<PT,IT> value_type;
  typedef typename set<pair<PT,IT> >::iterator iterator;
  typedef typename set<pair<PT,IT> >::const_iterator const_iterator;

  iterator begin() { return q.begin(); }
  const_iterator begin() const { return q.begin(); }
  iterator end() { return q.end(); }
  const iterator end() const { return q.end(); }
  const value_type& front() const { return *q.begin(); }
  const value_type& back() const { const_iterator i=q.end(); --i; return *i; }
  bool empty() const { return q.empty(); }
  void clear() { q.clear(); indices.clear(); }
  size_t size() const { return q.size(); }
  const value_type& top() const { return *q.begin(); }
  void erase(iterator i) {
    size_t n=indices.erase(i->second);
    Assert(n==1);
    q.erase(i);
  }
  void pop() { 
    erase(q.begin());
  }
  iterator insert(const value_type& item) {
    iterator i=q.insert(item).first;
    Assert(indices.count(item.second)==0);
    indices[item.second] = i;
    return i;
  }
  template <class InputIterator>
  void insert(InputIterator i,InputIterator j) {
    for(InputIterator k=i;k!=j;++k)
      insert(*k);
  }
  iterator insert(const IT& index,const PT& p) {
    iterator i=q.insert(value_type(p,index)).first;
    Assert(indices.count(index)==0);
    indices[index] = i;
    return i;
  }
  iterator find(const IT& index) {
    typename map<IT,iterator>::iterator i=indices.find(index);
    if(i==indices.end()) return end();
    return i->second;
  }
  iterator refresh(const IT& index,const PT& p) {
    iterator i=find(index);
    erase(i);
    insert(p,index);
  }
  bool is_valid() const {
    for(typename map<IT,iterator>::const_iterator i=indices.begin();i!=indices.end();++i) {
      if(i->second->second != i->first) return false;
    }
    return true;
  }

 private:
  map<IT,iterator> indices;
  set<pair<PT,IT> > q;
};

} //namespace std

#endif
