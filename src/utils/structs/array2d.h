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

#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <myfile.h>
#include <iostream>
#include <assert.h>

template <class T>
class Array2D
{
public:
	Array2D();
	Array2D(int m,int n);
	Array2D(int m,int n,const T& initVal);
	Array2D(const Array2D<T>& rhs);
	~Array2D();

	inline T& operator()(int i,int j) { assert(i>=0&&i<m); assert(j>=0&&j<n); return items[i*n+j]; }
	inline const T& operator()(int i,int j) const { assert(i>=0&&i<m); assert(j>=0&&j<n); return items[i*n+j]; }
	const Array2D<T>& operator =(const Array2D<T>& rhs);

	bool Read(File& f);
	bool Write(File& f) const;

	inline int numRows() const { return m; }
	inline int numCols() const { return n; }
	inline bool empty() const { return m==0&&n==0; }
	void initialize(int m,int n);
	void initialize(int m,int n,const T& initVal);
	void resize(int m,int n);
	void resize(int m,int n,const T& initVal);
	void clear();

	bool find(const T& item,int& i,int& j) const;
	inline bool contains(const T& item) const { int i,j; return find(item,i,j); }
	void set(const T& item);
	void set(const Array2D<T>&);
	void swap(Array2D<T>&);
	inline T* getData() const { return items; }
	inline T* getRowData(int i) const { return &items[i*n]; }

	//READ ONLY
	int m,n;
protected:
	T* items;
};

template <class T>
Array2D<T>::Array2D()
:m(0),n(0),items(0)
{}

template <class T>
Array2D<T>::Array2D(int _m,int _n)
:m(0),n(0),items(0)
{
	initialize(_m,_n);
}

template <class T>
Array2D<T>::Array2D(int _m,int _n,const T& initVal)
:m(0),n(0),items(0)
{
	initialize(_m,_n,initVal);
}

template <class T>
Array2D<T>::Array2D(const Array2D<T>& rhs)
:m(0),n(0),items(0)
{
	set(rhs);
}

template <class T>
Array2D<T>::~Array2D()
{
	clear();
}

template <class T>
const Array2D<T>& Array2D<T>::operator =(const Array2D<T>& rhs)
{
	set(rhs);
	return *this;
}

template <class T>
bool Array2D<T>::Read(File& f)
{
	if(!ReadFile(f,m)) return false;
	if(!ReadFile(f,n)) return false;
	initialize(m,n);
	if(!ReadArrayFile(f,items,m*n)) return false;
	return true;
}

template <class T>
bool Array2D<T>::Write(File& f) const
{
	if(!WriteFile(f,m)) return false;
	if(!WriteFile(f,n)) return false;
	if(!WriteArrayFile(f,items,m*n)) return false;
	return true;
}

template <class T>
void Array2D<T>::initialize(int _m,int _n)
{
	clear();
	m=_m;
	n=_n;
	items=new T[m*n];
}

template <class T>
void Array2D<T>::initialize(int _m,int _n,const T& initVal)
{
	clear();
	m=_m;
	n=_n;
	items=new T[m*n];
	set(initVal);
}

template <class T>
void Array2D<T>::resize(int _m,int _n)
{
	if(_m*_n != m*n)
		initialize(_m,_n);
	m=_m;
	n=_n;
}

template <class T>
void Array2D<T>::resize(int _m,int _n,const T& initVal)
{
	if(_m*_n != m*n)
		initialize(_m,_n);
	m=_m;
	n=_n;
	set(initVal);
}

template <class T>
void Array2D<T>::clear()
{
  if(items) {
    delete [] items;
    items = NULL;
  }
  m=n=0;
}

template <class T>
bool Array2D<T>::find(const T& item,int& i,int& j) const
{
	for(int p=0;p<m;p++)
		for(int q=0;q<n;q++)
			if(operator()(p,q)==item) {
				i=p; j=q;
				return true;
			}
	return false;
}

template <class T>
void Array2D<T>::set(const T& item)
{
	for(int i=0;i<m*n;i++) items[i]=item;
}

template <class T>
void Array2D<T>::set(const Array2D<T>& rhs)
{
	resize(rhs.m,rhs.n);
	for(int i=0;i<m*n;i++) items[i]=rhs.items[i];
}

template <class T>
void Array2D<T>::swap(Array2D<T>& b)
{
  int tempm=m,tempn=n;
  T* tempitems=items;
  m=b.m; n=b.n;
  items=b.items;
  b.m=tempm; b.n=tempn;
  b.items=tempitems;
}


template <class T>
std::ostream& operator << (std::ostream& out,const Array2D<T>& a)
{
  out<<a.m<<" "<<a.n<<std::endl;
  for(int i=0;i<a.m;i++) {
    for(int j=0;j<a.n;j++) {
      out<<a(i,j);
      if(j+1 < a.n) out<<" ";
    }
    if(i+1 < a.m) out<<std::endl;
  }
  return out;
}

template <class T>
std::istream& operator >> (std::istream& in,Array2D<T>& a)
{
  int m,n;
  in>>m>>n;
  if(in.bad()) return in;
  assert(m>=0 && n>=0);
  a.resize(m,n);
  for(int i=0;i<m;i++) {
    for(int j=0;j<n;j++) {
      in>>a(i,j);
    }
  }
  return in;
}

namespace std {

  template <class T>
  void swap(Array2D<T>& a,Array2D<T>&b)
  {
    a.swap(b);
  }

} //namespace std

#endif
