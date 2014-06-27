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

#include "gaussian.h"
#include "CholeskyDecomposition.h"
#include "backsubstitute.h"
#include "random.h"
using namespace std;

namespace Math {

template <class T>
Gaussian<T>::Gaussian()
{
}

template <class T>
Gaussian<T>::Gaussian(int d)
{
  resize(d);
}

template <class T>
Gaussian<T>::Gaussian(const MatrixT& sigma, const VectorT& _mu)
{
  setCovariance(sigma);
  mu = _mu;
}

template <class T>
Gaussian<T>::Gaussian(const Gaussian& g)
  :L(g.L),mu(g.mu)
{
}

template <class T>
void Gaussian<T>::resize(int d)
{
  L.resize(d,d);
  mu.resize(d,Zero);
}

template <class T>
int Gaussian<T>::numDims() const
{
  return mu.n;
}

template <class T>
bool Gaussian<T>::setCovariance(const MatrixT& sigma)
{
  if(sigma.m != mu.n) {
    cerr<<"Invalid dimensions on covariance matrix"<<endl;
    return false;
  }
  CholeskyDecomposition<T> chol;
  if(!chol.set(sigma)) {
    cerr<<"Unable to set cholesky decomposition of covariance matrix"<<endl;
    return false;
  }
  L = chol.L;
  return true;
}

template <class T>
void Gaussian<T>::getCovariance(MatrixT& sigma) const
{
  sigma.mulTransposeB(L,L);
}

template <class T>
void Gaussian<T>::setMean(const VectorT& _mu)
{
  if(_mu.n != L.m) {
    cerr<<"Invalid dimensions on mean vector"<<endl;
    abort();
  }
  mu = _mu;
}

template <class T>
T Gaussian<T>::probability(const VectorT& x) const
{
  int d = numDims();
  T det=One;
  for(int i=0;i<d;i++) det*=L(i,i);
  T invc = Pow(2.0*Pi,0.5*T(d))*det;
  VectorT x_mu,y;
  x_mu.sub(x,mu);
  y.resize(L.n);
  LBackSubstitute(L,x_mu,y);
  return Exp(-Half*y.normSquared())/invc;
}

template <class T>
void Gaussian<T>::generate(VectorT& x) const
{
  //generate a normalized gaussian distributed variable y
  VectorT y(numDims());
  for(int i=0;i<y.n;i++) y(i) = RandGaussian(0,1);
  L.mul(y,x);
  x += mu;
}

template class Gaussian<float>;
template class Gaussian<double>;

} //namespace Math

