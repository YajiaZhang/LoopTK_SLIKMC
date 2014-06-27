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

#ifndef ARRAY3D_H
#define ARRAY3D_H

#include "myfile.h"
#include <assert.h>

template <class type>
class Array3D
{
public:
	Array3D();
	Array3D(int m,int n,int p);
	Array3D(int m,int n,int p,const type& initVal);
	Array3D(const Array3D<type>& rhs);
	~Array3D();

	inline type& operator()(int i,int j,int k) { assert(i>=0&&i<m); assert(j>=0&&j<n); assert(k>=0&&k<p); return items[i*n*p+j*p+k]; }
	inline const type& operator()(int i,int j,int k) const { assert(i>=0&&i<m); assert(j>=0&&j<n); return items[i*n*p+j*p+k]; }
	const Array3D<type>& operator =(const Array3D<type>& rhs);

	bool Read(File& f);
	bool Write(File& f) const;

	inline int dim1() const { return m; }
	inline int dim2() const { return n; }
	inline int dim3() const { return p; }
	inline bool empty() const { return m==0&&n==0&&p==0; }
	void initialize(int m,int n,int p);
	void initialize(int m,int n,int p,const type& initVal);
	void resize(int m,int n,int p);
	void resize(int m,int n,int p,const type& initVal);
	void clear();

	bool find(const type& item,int& i,int& j,int& k) const;
	inline bool contains(const type& item) const { int i,j,k; return find(item,i,j); }
	void set(const type& item);
	void set(const Array3D<type>&);
	void swap(Array3D<type>&);
	inline type* getData() const { return items; }

	//READ ONLY
	int m,n,p;
protected:
	type* items;
};

template <class type>
Array3D<type>::Array3D()
:m(0),n(0),p(0),items(0)
{}

template <class type>
Array3D<type>::Array3D(int _m,int _n,int _p)
:m(0),n(0),p(0),items(0)
{
	initialize(_m,_n,_p);
}

template <class type>
Array3D<type>::Array3D(int _m,int _n,int _p,const type& initVal)
:m(0),n(0),p(0),items(0)
{
	initialize(m,n,p,initVal);
}

template <class type>
Array3D<type>::Array3D(const Array3D<type>& rhs)
:m(0),n(0),p(0),items(0)
{
	set(rhs);
}

template <class type>
Array3D<type>::~Array3D()
{
	clear();
}

template <class type>
const Array3D<type>& Array3D<type>::operator =(const Array3D<type>& rhs)
{
	set(rhs);
	return *this;
}

template <class type>
bool Array3D<type>::Read(File& f)
{
	if(!ReadFile(f,m)) return false;
	if(!ReadFile(f,n)) return false;
	if(!ReadFile(f,p)) return false;
	initialize(m,n,p);
	if(!ReadArrayFile(f,items,m*n*p)) return false;
	return true;
}

template <class type>
bool Array3D<type>::Write(File& f) const
{
	if(!WriteFile(f,m)) return false;
	if(!WriteFile(f,n)) return false;
	if(!WriteFile(f,p)) return false;
	if(!WriteArrayFile(f,items,m*n*p)) return false;
	return true;
}

template <class type>
void Array3D<type>::initialize(int _m,int _n,int _p)
{
	clear();
	m=_m;
	n=_n;
	p=_p;
	items=new type[m*n*p];
}

template <class type>
void Array3D<type>::initialize(int _m,int _n,int _p,const type& initVal)
{
	clear();
	m=_m;
	n=_n;
	p=_p;
	items=new type[m*n*p];
	set(initVal);
}

template <class type>
void Array3D<type>::resize(int _m,int _n,int _p)
{
	if(_m*_n*_p != m*n*p)
		initialize(_m,_n,_p);
	m=_m;
	n=_n;
	p=_p;
}

template <class type>
void Array3D<type>::resize(int _m,int _n,int _p,const type& initVal)
{
	if(_m*_n*_p != m*n*p)
		initialize(_m,_n,_p);
	m=_m;
	n=_n;
	p=_p;
	set(initVal);
}

template <class type>
void Array3D<type>::clear()
{
	SafeArrayDelete(items);
	m=n=p=0;
}

template <class type>
bool Array3D<type>::find(const type& item,int& i,int& j,int& k) const
{
  for(int s=0;s<m;s++)
    for(int t=0;t<n;t++)
      for(int u=0;u<p;u++)
	if(operator()(s,t,u)==item) {
	  i=s; j=t; k=u;
	  return true;
	}
  return false;
}

template <class type>
void Array3D<type>::set(const type& item)
{
	for(int i=0;i<m*n*p;i++) items[i]=item;
}

template <class type>
void Array3D<type>::set(const Array3D<type>& rhs)
{
	resize(rhs.m,rhs.n,rhs.p);
	for(int i=0;i<m*n*p;i++) items[i]=rhs.items[i];
}

template <class type>
void Array3D<type>::swap(Array3D<type>& b)
{
  int tempm=m,tempn=n,tempp=p;
  type* tempitems=items;
  m=b.m; n=b.n; p=b.p;
  items=b.items;
  b.m=tempm; b.n=tempn; b.p=tempp;
  b.items=tempitems;
}

namespace std {

  template <class type>
  void swap(Array3D<type>& a,Array3D<type>&b)
  {
    a.swap(b);
  }

} //namespace std

#endif
