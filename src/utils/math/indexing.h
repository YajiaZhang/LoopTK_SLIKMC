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

#ifndef MATH_INDEXING_H
#define MATH_INDEXING_H

#include "VectorTemplate.h"
#include "MatrixTemplate.h"
#include <vector>

/** @ingroup Math
 * @file math/indexing.h
 * @brief Utilities for matrix/vector access / manipulation with
 * index vectors.
 *
 * Be careful, little error checking is performed.  Each index
 * must be valid for the matrix/vector operated on.
 */

namespace Math {

/** @addtogroup Math */
/*@{*/

///Sets A[indices] = fill
template <class T>
inline void SetElements(VectorTemplate<T>& A,const std::vector<int>& indices,T fill)
{
  for(size_t i=0;i<indices.size();i++)
    A(indices[i]) = fill;
}

/** @brief Sets A[indices] = fill
 *
 * The fill vector is of size indices.size().  For copying from
 * a vector to another using the same set of indices, use
 * CopyElements().
 */
template <class T>
inline void SetElements(VectorTemplate<T>& A,const std::vector<int>& indices,const VectorTemplate<T>& fill)
{
  for(size_t i=0;i<indices.size();i++)
    A(indices[i]) = fill(i);
}

/// Sets B = A[indices]
template <class T>
inline void GetElements(const VectorTemplate<T>& A,const std::vector<int>& indices,VectorTemplate<T>& B)
{
  for(size_t i=0;i<indices.size();i++)
    B(i) = A(indices[i]);
}

/// Sets A[aindices] = B[bindices]
template <class T>
inline void CopyElements(VectorTemplate<T>& A,const std::vector<int>& aindices,
			 const VectorTemplate<T>& B,const std::vector<int>& bindices)
{
  for(size_t i=0;i<aindices.size();i++) 
    A(aindices[i]) = B(bindices[i]);
}

///Removes the indexed elements of A (shrinking A)
template <class T>
void RemoveElements(VectorTemplate<T>& A,const std::vector<int>& indices)
{
  for(size_t i=0;i<indices.size();i++) {
    int start=indices[i]+1;
    int end=(i+1==indices.size() ? A.n : indices[i+1]);
    assert(end >= start);
    assert(start > 0 && end <= A.n);
    for(int j=start;j<end;j++) {
      A(j-i-1) = A(j);
    }
  }
  A.n = A.n - (int)indices.size();
}

/** @brief The inverse of RemoveElements.
 *
 * Expands A such that the elements indexed by indices are set to fill.
 * A must be able to expand to size A.n + indices.size()
 */
template <class T>
void AddElements(VectorTemplate<T>& A,const std::vector<int>& indices,T fill)
{
  assert(A.getCapacity() >= A.base + A.stride*(A.n+(int)indices.size()));
  A.n = A.n + (int)indices.size();
  assert(A.isValid());
  for(int i=(int)indices.size()-1;i>=0;i--) {
    int start=indices[i]+1;
    int end=(i+1==(int)indices.size() ? A.n : indices[i+1]);
    assert(end >= start);
    assert(start > 0 && end <= A.n);
    for(int j=end-1;j>=start;j--) {
      A(j) = A(j-i-1);
    }
    A(indices[i]) = fill;
  }
}

///Sets the indexed rows of A to fill 
template <class T>
inline void SetRows(MatrixTemplate<T>& A,const std::vector<int>& indices,T fill)
{
  for(size_t i=0;i<indices.size();i++)
    A.setRow(indices[i],fill);
}

///Sets the indexed columns of A to fill 
template <class T>
inline void SetColumns(MatrixTemplate<T>& A,const std::vector<int>& indices,T fill)
{
  for(size_t i=0;i<indices.size();i++)
    A.setCol(indices[i],fill);
}

///Sets the indexed rows of A to the vector fill 
template <class T>
inline void SetRows(MatrixTemplate<T>& A,const std::vector<int>& indices,const VectorTemplate<T>& fill)
{
  for(size_t i=0;i<indices.size();i++)
    A.copyRow(indices[i],fill);
}

///Sets the indexed columns of A to the vector fill 
template <class T>
inline void SetColumns(MatrixTemplate<T>& A,const std::vector<int>& indices,const VectorTemplate<T>& fill)
{
  for(size_t i=0;i<indices.size();i++)
    A.copyColumn(indices[i],fill);
}

/** @brief Sets the indexed rows of A to the columns of the matrix fill
 *
 * The fill matrix is indices.size() by A.n.  For copying from
 * a matrix to another using the same set of indices, use
 * CopyRows().
 */
template <class T>
inline void SetRows(MatrixTemplate<T>& A,const std::vector<int>& indices,const MatrixTemplate<T>& fill)
{
  VectorTemplate<T> filli;
  for(size_t i=0;i<indices.size();i++) {
    fill.getRowRef(i,filli);
    A.copyRow(indices[i],filli);
  }
}

/** @brief Sets the indexed columns of A to the columns of the matrix fill
 *
 * The fill matrix is A.m by indices.size().  For copying from
 * a matrix to another using the same set of indices, use
 * CopyColumns().
 */
template <class T>
inline void SetColumns(MatrixTemplate<T>& A,const std::vector<int>& indices,const MatrixTemplate<T>& fill)
{
  VectorTemplate<T> filli;
  for(size_t i=0;i<indices.size();i++) {
    fill.getColRef(i,filli);
    A.copyColumn(indices[i],filli);
  }
}

/** @brief Copies the indexed rows in A into B
 *
 * B has dimensions indices.size() by A.n.
 */
template <class T>
inline void GetRows(MatrixTemplate<T>& A,const std::vector<int>& indices,const MatrixTemplate<T>& B)
{
  VectorTemplate<T> Bi;
  for(size_t i=0;i<indices.size();i++) {
    B.getRowRef(i,Bi);
    A.getRowCopy(indices[i],Bi);
  }
}

/** @brief Copies the indexed columns in A into B
 *
 * B has dimensions A.m by indices.size()
 */
template <class T>
inline void GetColumns(const MatrixTemplate<T>& A,const std::vector<int>& indices,VectorTemplate<T>& B)
{
  VectorTemplate<T> Bi;
  for(size_t i=0;i<indices.size();i++) {
    B.getColRef(i,Bi);
    A.getColCopy(indices[i],Bi);
  }
}

///Removes the rows of A indexed by the <em>sorted</em> list indices
template <class T>
void RemoveRows(MatrixTemplate<T>& A,const std::vector<int>& indices)
{
  VectorTemplate<T> src,dest;
  for(size_t i=0;i<indices.size();i++) {
    int start=indices[i]+1;
    int end=(i+1==indices.size() ? A.m : indices[i+1]);
    assert(end >= start);
    assert(start > 0 && end <= A.m);
    for(int j=start;j!=end;j++) {
      A.getRowRef(j,src);
      A.getRowRef(j-i-1,dest);
      dest.copy(src);
    }
  }
  A.m = A.m - (int)indices.size();
}

///Removes the columns of A indexed by the <em>sorted</em> indices
template <class T>
void RemoveColumns(MatrixTemplate<T>& A,const std::vector<int>& indices)
{
  VectorTemplate<T> src,dest;
  for(size_t i=0;i<indices.size();i++) {
    int start=indices[i]+1;
    int end=(i+1==indices.size() ? A.n : indices[i+1]);
    assert(end >= start);
    assert(start > 0 && end <= A.n);
    for(int j=start;j<end;j++) {
      A.getColRef(j,src);
      A.getColRef(j-i-1,dest);
      dest.copy(src);
    }
  }
  A.n = A.n - (int)indices.size();
}

/*@}*/

} //namespace Math

#endif




