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

#include "rotation.h"
#include <math/misc.h>
using namespace std;

namespace Math3D {

const static Real angleEps=(Real)1e-3;

EulerAngleRotation::EulerAngleRotation()
{}

EulerAngleRotation::EulerAngleRotation(const EulerAngleRotation& r)
  :Vector3(r)
{}

EulerAngleRotation::EulerAngleRotation(const Vector3& v)
  :Vector3(v)
{}

EulerAngleRotation::EulerAngleRotation(Real a, Real b, Real c)
  :Vector3(a,b,c)
{}


void EulerAngleRotation::setMatrix(int u,int v,int w,const Matrix3& m)
{
  if(u==0&&v==1&&w==2) {
    setMatrixXYZ(m);
  }
  else if(u==2&&v==1&&w==0) {
    setMatrixZYX(m);
  }
  else if(u==2&&v==0&&w==1) {
    setMatrixZXY(m);
  }
  else if(u==1&&v==0&&w==2) {
    setMatrixYXZ(m);
  }
  else {
    cerr<<"Not done with general euler angle rotation setMatrix"<<endl;
    abort();
  }
}

void EulerAngleRotation::setMatrixXYZ(const Matrix3& m)
{
  Real a,b,c;
  b=Asin(m(0,2));  //m(0,2)=sb
  Real cb = Cos(b);
  if (Abs(cb) > Epsilon) {
    Real ca, cc;
    ca = m(2,2)/cb;   //m(2,2)=ca*cb
    ca = Clamp(ca,-One,One);
    if (Sign(m(1,2)) == Sign(cb)) //m(1,2)=-sa*cb
      a = TwoPi - Acos(ca);
    else
 	    a = Acos(ca);

    cc = m(0,0) / cb;  //m(0,0)=cb*cc
    cc = Clamp(cc,-One,One);
    if (Sign(m(0,1)) != Sign(cb)) //m(0,1)=cb*sc
  		c = Acos(cc);
    else
      c = TwoPi - Acos(cc);
  }
  else {
    // b is close to 90 degrees, i.e. cb=0
    // this reduces the degrees of freedom, so we can set c=0
    c = 0;
    a = Acos(m(1,1)); //m(1,1)=ca
    if(Sign(Sin(a)) != Sign(m(2,1)))  //m(2,1)=sa
      a = TwoPi - a;
  }

  x=(Real)a;
  y=(Real)b;
  z=(Real)c;
}

void EulerAngleRotation::setMatrixZYX(const Matrix3& m)
{
  Real a,b,c;
  b=-Asin(m(2,0));  //m(2,0)=-sb
  Real cb = Cos(b);
  if (Abs(cb) > Epsilon) {
    Real ca, cc;
    ca = m(0,0)/cb;   //m(0,0)=ca*cb
    ca = Clamp(ca,-One,One);
    if (Sign(m(1,0)) == Sign(cb)) //m(1,0)=sa*cb
      a = Acos(ca);
    else
      a = TwoPi - Acos(ca);

    cc = m(2,2) / cb;  //m(2,2)=cb*cc
    cc = Clamp(cc,-One,One);
    if (Sign(m(2,1)) == Sign(cb)) //m(2,1)=cb*sc
  		c = Acos(cc);
    else
      c = TwoPi - Acos(cc);
  }
  else {
    // b is close to 90 degrees, i.e. cb=0
    // this reduces the degrees of freedom, so we can set c=0
    c = 0;
    //m(0,1)=-sa
    a = -Asin(m(0,1));
    if(Sign(Cos(a)) != Sign(m(1,1))) //m(1,1)=ca
      a = Pi - a;
  }

  x=(Real)a;
  y=(Real)b;
  z=(Real)c;
}

void EulerAngleRotation::setMatrixZXY(const Matrix3& m)
{
  Real a,b,c;
  b=Asin(m(2,1));  //m(2,1)=sb
  Real cb = Cos(b);
  if (Abs(cb) > Epsilon) {
    Real ca, cc;
    ca = m(1,1)/cb;   //m(1,1)=ca*cb
    ca = Clamp(ca,-One,One);
    if (Sign(m(0,1)) != Sign(cb)) //m(0,1)=-sa*cb
      a = Acos(ca);
    else
      a = TwoPi - Acos(ca);

    cc = m(2,2) / cb;  //m(2,2)=ca*cc
    cc = Clamp(cc,-One,One);
    if (Sign(m(2,0)) != Sign(cb)) //m(2,0)=-cb*sc
      c = Acos(cc);
    else
      c = TwoPi - Acos(cc);
  }
  else {
    // b is close to 90 degrees, i.e. cb=0
    // this reduces the degrees of freedom, so we can set c=0
    c = 0;
    //m(1,0)=sa
    a = Asin(m(1,0));
    if(Sign(Cos(a)) != Sign(m(0,0))) //m(0,0)=ca
      a = Pi - a;
  }

  x=(Real)a;
  y=(Real)b;
  z=(Real)c;
}

void EulerAngleRotation::setMatrixYXZ(const Matrix3& m)
{
  Real a,b,c;
  b=-Asin(m(1,2));  //m(1,2)=-sb
  Real cb = Cos(b);
  if (Abs(cb) > Epsilon) {
    Real ca, cc;
    ca = m(2,2)/cb;   //m(2,2)=ca*cb
    ca = Clamp(ca,-One,One);
    if (Sign(m(0,2)) == Sign(cb)) //m(0,2)=sa*cb
      a = Acos(ca);
    else
      a = TwoPi - Acos(ca);

    cc = m(1,1) / cb;  //m(1,1)=cb*cc
    cc = Clamp(cc,-One,One);
    if (Sign(m(1,0)) == Sign(cb)) //m(1,0)=cb*sc
      c = Acos(cc);
    else
      c = TwoPi - Acos(cc);
  }
  else {
    // b is close to 90 degrees, i.e. cb=0
    // this reduces the degrees of freedom, so we can set c=0
    c = 0;
    //m(2,0)=-sa
    a = -Asin(m(2,0));
    if(Sign(Cos(a)) != Sign(m(0,0))) //m(0,0)=ca
      a = Pi - a;
  }

  x=(Real)a;
  y=(Real)b;
  z=(Real)c;
}

void EulerAngleRotation::getMatrix(int u,int v,int w,Matrix3& m) const
{
  Matrix3 Ru,Rv,Rw;
  switch(u) {
  case 0:  Ru.setRotateX(x);  break;
  case 1:  Ru.setRotateY(x);  break;
  case 2:  Ru.setRotateZ(x);  break;
  default: cerr<<"Invalid axis "<<u<<endl; break;
  }
  switch(v) {
  case 0:  Rv.setRotateX(y);  break;
  case 1:  Rv.setRotateY(y);  break;
  case 2:  Rv.setRotateZ(y);  break;
  default: cerr<<"Invalid axis "<<v<<endl; break;
  }
  switch(w) {
  case 0:  Rw.setRotateX(z);  break;
  case 1:  Rw.setRotateY(z);  break;
  case 2:  Rw.setRotateZ(z);  break;
  default: cerr<<"Invalid axis "<<w<<endl; break;
  }
  m = Ru*Rv*Rw;
}


void EulerAngleRotation::getMatrixXYZ(Matrix3& m) const
{
  Real ca = Cos(x);
  Real cb = Cos(y);
  Real cc = Cos(z);
  Real sa = Sin(x);
  Real sb = Sin(y);
  Real sc = Sin(z);

  //remember data is column major
  m.data[0][0]=cb*cc;
  m.data[0][1]=sa*sb*cc+ca*sc;
  m.data[0][2]=sa*sc-ca*sb*cc; 
  m.data[1][0]=-cb*sc;
  m.data[1][1]=ca*cc-sa*sb*sc; 
  m.data[1][2]=ca*sb*sc+sa*cc;
  m.data[2][0]=sb; 
  m.data[2][1]=-sa*cb; 
  m.data[2][2]=ca*cb;
}

void EulerAngleRotation::getMatrixZYX(Matrix3& m) const
{
  Real ca = Cos(x);
  Real cb = Cos(y);
  Real cc = Cos(z);
  Real sa = Sin(x);
  Real sb = Sin(y);
  Real sc = Sin(z);

  //remember data is column major
  m.data[0][0]=ca*cb;
  m.data[0][1]=sa*cb;
  m.data[0][2]=-sb; 
  m.data[1][0]=ca*sb*sc-sa*cc;
  m.data[1][1]=sa*sb*sc+ca*cc; 
  m.data[1][2]=cb*sc;
  m.data[2][0]=ca*sb*cc+sa*sc; 
  m.data[2][1]=sa*sb*cc-ca*sc; 
  m.data[2][2]=cb*cc;
}

void EulerAngleRotation::getMatrixZXY(Matrix3& m) const
{
  Real ca = Cos(x);
  Real cb = Cos(y);
  Real cc = Cos(z);
  Real sa = Sin(x);
  Real sb = Sin(y);
  Real sc = Sin(z);

  //remember data is column major
  m.data[0][0]=-sa*sb*sc+ca*cc;
  m.data[0][1]=ca*sb*sc+sa*cc;
  m.data[0][2]=-cb*sc; 
  m.data[1][0]=-sa*cb;
  m.data[1][1]=ca*cb; 
  m.data[1][2]=sb;
  m.data[2][0]=sa*sb*cc+ca*sc; 
  m.data[2][1]=-ca*sb*cc+sa*sc; 
  m.data[2][2]=cb*cc;
}

void EulerAngleRotation::getMatrixYXZ(Matrix3& m) const
{
  Matrix3 ma,mb,mc;
  ma.setRotateY(x);
  mb.setRotateX(y);
  mc.setRotateZ(z);
  m = ma*mb*mc;
}


AngleAxisRotation::AngleAxisRotation()
{}

AngleAxisRotation::AngleAxisRotation(const AngleAxisRotation& r)
  :angle(r.angle), axis(r.axis)
{}

AngleAxisRotation::AngleAxisRotation(Real _angle, const Vector3& _axis)
{
  set(_angle, _axis);
}

AngleAxisRotation::AngleAxisRotation(const MomentRotation& r) 
{
  setMoment(r);
}

void AngleAxisRotation::set(const AngleAxisRotation& r)
{
  angle = r.angle;
  axis = r.axis;
}

void AngleAxisRotation::set(Real _angle, const Vector3& _axis)
{
  angle = _angle;
  setAxis(_axis);
}

void AngleAxisRotation::setAxis(const Vector3& _axis)
{
  Real l2 = _axis.normSquared();
  if(FuzzyEquals(l2,One))
    axis = _axis;
  else
    axis.div(_axis, Inv(Sqrt(l2)));
}

void AngleAxisRotation::setIdentity()
{
  angle = Zero;
  axis.set(One,Zero,Zero);
}

void AngleAxisRotation::transformPoint(const Vector3& in,Vector3& out) const
{
  Real cm = Cos(angle);
  Real sm = Sin(angle);

  //m = s[r]-c[r][r]+rrt = s[r]-c(rrt-I)+rrt = cI + rrt(1-c) + s[r]
  //=> mv = cv + (1-c)*(r.v)r + s*(r x v)
  out.setCross(axis,in);
  out *= sm;
  out.madd(axis,(One-cm)*dot(axis,in));
  out.madd(in,cm);
}

void AngleAxisRotation::setMatrix(const Matrix3& r)
{
  MomentRotation m;
  m.setMatrix(r);
  setMoment(m);
}

void AngleAxisRotation::getMatrix(Matrix3& m) const
{
  if(angle == Zero) {
    m.setIdentity();
    return;
  }

  Real cm = Cos(angle);
  Real sm = Sin(angle);

  //m = s[r]-c[r][r]+rrt = s[r]-c(rrt-I)+rrt = cI + rrt(1-c) + s[r]
  m.setCrossProduct(axis);
  m.inplaceScale(sm);
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++)
      m(i,j) += axis[i]*axis[j]*(One-cm);
    m(i,i) += cm;
  }
  /*
  Matrix3 rmat;
  rmat.setCrossProduct(axis);

  m.mul(rmat, rmat);
  m.inplaceScale(-cm);

  rmat.inplaceScale(sm);
  m += rmat;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      m.data[i][j] += axis[i]*axis[j];
  */
}

void AngleAxisRotation::setMoment(const MomentRotation& r)
{
  Real l = r.norm();
  angle = l;
  if(l == Zero)
    axis.setZero();
  else
    axis.div(r, l);
}

void AngleAxisRotation::getMoment(MomentRotation& r) const
{
  r.mul(axis, angle);
}

MomentRotation::MomentRotation()
{}

MomentRotation::MomentRotation(const MomentRotation& r)
  :Vector3(r)
{}

MomentRotation::MomentRotation(const Vector3& v)
  :Vector3(v)
{}

MomentRotation::MomentRotation(const AngleAxisRotation& a)
{
  setAngleAxis(a);
}

void MomentRotation::transformPoint(const Vector3& in,Vector3& out) const
{
  AngleAxisRotation a;
  getAngleAxis(a);
  a.transformPoint(in,out);
}

void MomentRotation::setMatrix(const Matrix3& r)
{
  //calculate theta
  Real theta;
  Real s=(r.trace() - One)*Half;
  if(s >= One) {
    if(s > One+Epsilon) {
      cerr<<"MomentRotation::setMatrix(): Warning- trace of matrix is greater than 3"<<endl;
    }
    theta = Zero;
  }
  else if(s <= -One) {
    if(s < -One-Epsilon) {
      cerr<<"MomentRotation::setMatrix(): Warning- trace of matrix is less than -1"<<endl;
    }
    theta = Pi;
  }
  else 
    theta = Acos(s);

  if(FuzzyEquals(theta,Pi,angleEps)) {
    //can't do normal version because the scale factor reaches a singularity
    x = Pi*Sqrt((r(0,0)+One)*Half);
    y = Pi*Sqrt((r(1,1)+One)*Half);
    z = Pi*Sqrt((r(2,2)+One)*Half);
    //determined up to a scale factor, we know r12=2xy,r13=2xz,r23=2yz
    Real xy=r(0,1),xz=r(0,2),yz=r(1,2);
    if(xy > Zero) {
      //assert(Sign(xz) == Sign(yz));
      if(xz < Zero) z=-z;
    }
    else {
      //assert(Sign(xz) != Sign(yz));
      if(xz > Zero) //x and z have same sign
	y=-y;
      else //y and z have same sign
	x=-x;
    }
    return;
  }

  //Matrix3 r_cross;
  //r_cross.setTranspose(r);
  //r_cross.sub(r, r_cross);
  //r_cross.inplaceScale(scale);

  //Real scale = Half*theta/Sin(theta);
  Real scale = Half/Sinc(theta);  //avoids the singularity at 0

  x = (r.data[1][2]-r.data[2][1]) * scale;
  y = (r.data[2][0]-r.data[0][2]) * scale;
  z = (r.data[0][1]-r.data[1][0]) * scale;
}

void MomentRotation::getMatrix(Matrix3& m) const
{
  AngleAxisRotation a;
  a.setMoment(*this);
  a.getMatrix(m);
}

void MomentRotation::setAngleAxis(const AngleAxisRotation& a)
{
  mul(a.axis, a.angle);
}

void MomentRotation::getAngleAxis(AngleAxisRotation& a) const
{
  a.setMoment(*this);
}

QuaternionRotation::QuaternionRotation()
{}

QuaternionRotation::QuaternionRotation(Real w, Real x, Real y, Real z)
:Quaternion(w,x,y,z)
{}

QuaternionRotation::QuaternionRotation(const QuaternionRotation& q)
:Quaternion(q)
{}

QuaternionRotation::QuaternionRotation(const Quaternion& q)
:Quaternion(q)
{}


void QuaternionRotation::setAngleAxis(const AngleAxisRotation& r)
{
  Real imscale = Sin(r.angle*Half);
  w = Cos(r.angle*Half);
  x = imscale*r.axis.x;
  y = imscale*r.axis.y;
  z = imscale*r.axis.z;
}


void QuaternionRotation::setMoment(const MomentRotation& m)
{
  setAngleAxis(AngleAxisRotation(m));
}

void QuaternionRotation::setMatrix(const Matrix3& m)
{
  Real tr, s;
  tr = m.trace() + One;

  // check the diagonal
  if (tr > 0.00001f) {
    s = Sqrt (tr);
    w = s * Half;
    s = Half / s;
    x = (m.data[1][2] - m.data[2][1]) * s;
    y = (m.data[2][0] - m.data[0][2]) * s;
    z = (m.data[0][1] - m.data[1][0]) * s;
  }
  else {             
    const static int nxt[3] = {1, 2, 0};
    int i, j, k;
    // diagonal is negative
    i = 0;
    if (m.data[1][1] > m.data[0][0]) i = 1;
    if (m.data[2][2] > m.data[i][i]) i = 2;
    j = nxt[i];
    k = nxt[j];
    
    Real q[4];
    s = Sqrt ((m.data[i][i] - (m.data[j][j] + m.data[k][k])) + One);
    q[i] = s * Half;
    
    if (s != Zero) s = Half / s;
    
    q[3] = (m.data[j][j] - m.data[k][j]) * s;
    q[j] = (m.data[j][i] + m.data[i][j]) * s;
    q[k] = (m.data[k][i] + m.data[i][i]) * s;

    x = q[0];
    y = q[1];
    z = q[2];
    w = q[3];
  }
}

void QuaternionRotation::getAngleAxis(AngleAxisRotation& r) const
{
  r.angle = Two*Acos(w);
  Real imscale = One/imNorm();
  r.axis.set(x*imscale, y*imscale, z*imscale);
}

void QuaternionRotation::getMoment(MomentRotation& r) const
{
  if(w>=One||w<=-One) r.setZero();
  else {
    Real angle = Two*Acos(w);
    Real den=pythag_leg(w,One);
    if(FuzzyEquals(den,Zero)) {
      r.setZero();
    }
    else {
      Real imscale = angle/den;
      r.set(x*imscale, y*imscale, z*imscale);
    }
  }
}

void QuaternionRotation::getMatrix(Matrix3& m) const {
  Real wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;

  // calculate coefficients
  x2 = x + x; y2 = y + y; z2 = z + z;
  xx = x * x2;   xy = x * y2;   xz = x * z2;
  yy = y * y2;   yz = y * z2;   zz = z * z2;
  wx = w * x2;   wy = w * y2;   wz = w * z2;

  m.data[0][0] = One - (yy + zz);    m.data[1][0] = xy - wz; 			m.data[2][0] = xz + wy;
  m.data[0][1] = xy + wz;            m.data[1][1] = One - (xx + zz);	m.data[2][1] = yz - wx;
  m.data[0][2] = xz - wy;            m.data[1][2] = yz + wx;			m.data[2][2] = One - (xx + yy);
}

void QuaternionRotation::transform(const Vector3& v, Vector3& out) const
{
  Quaternion qinv, x(0,v), tmp;
  qinv.setConjugate(*this);
  tmp.mul(x,*this);
  x.mul(qinv,tmp);
  out.x=x.x;
  out.y=x.y;
  out.z=x.z;
}

void QuaternionRotation::slerp(const Quaternion& a, const Quaternion& b, Real t)
{
  //a + b unit quaternions?
  /* angle theta is defined as
     cos(theta)*|a||b| = <a,b>
  */
  Real d = dot(a,b);
  Real bscale=One;
  if(d < 0) { bscale=-One; d=-d; }
  if(FuzzyEquals(d,One,angleEps)) {	//axes are the same axis
    set(b);
    return;
  }

  Real theta = Acos(d);
  Real sininv = Sin(theta);
  sininv = One/sininv;

  //this = (Sin((One-t)*theta)*sininv) * a +  (Sin(t*theta)*sininv) * b;
  Real a_coeff = Sin((One-t)*theta)*sininv;
  Real b_coeff = Sin(t*theta)*sininv;
  mul(a, a_coeff);
  madd(b, bscale*b_coeff);
}

void QuaternionRotation::mag(const Quaternion& a, Real t)
{
  //just like slerp(identity transform, a)

  if(a.w == One) {	//a is identity
    set(a);
    return;
  }
  else if(a.w == -One) { //a is a rotation of 180 degrees
	cerr<<"QuaternionRotation.mag(): Quaternion is a rotation of 180 degrees"<<endl;
    return;
  }


  Real theta = Acos(a.w);
  Real sininv = Sin(theta);
  sininv = One/sininv;

  Real a_coeff = Sin((One-t)*theta)*sininv;
  Real b_coeff = (Sin(t*theta)*sininv);
  mul(a, b_coeff);
  w += a_coeff;
}


void SLerp(const Quaternion& a,
	   const Quaternion& b,
	   Quaternion& out,
	   Real t)
{
  //a + b unit quaternions?
  /* angle theta is defined as
     cos(theta)*|a||b| = <a,b>
  */
  Real d = dot(a,b);
  if(d == One) {	//axes are the same axis
    out.set(b);
    return;
  }
  else if(d == -One) {	//axes are opposing axis
	cerr<<"SLerp(): Quaternions on opposing sides of unit sphere"<<endl;
    return;
  }

  Real theta = Acos(d);
  Real sininv = Sin(theta);
  sininv = One/sininv;

  //out = (Sin((One-t)*theta)*sininv) * a +  (Sin(t*theta)*sininv) * b;
  Real a_coeff = Sin((One-t)*theta)*sininv;
  Real b_coeff = Sin(t*theta)*sininv;
  out.mul(a, a_coeff);
  out.madd(b, b_coeff);
}

void SCerp(const Quaternion& q_1,
	   const Quaternion& q0,
	   const Quaternion& q1,
	   const Quaternion& q2,
	   Quaternion& out,
	   Real t)
{
  //calc control pts
  Quaternion c0,c1;
  Quaternion p,q;

  static const Real Third = (Real)0.3333333333333333333333333333;
  SLerp(q0, q_1, p, -1.0);
  SLerp(p, q1, q, 0.5);
  SLerp(q0, q, c0, Third);

  SLerp(q1, q2, p, -1.0);
  SLerp(p, q0, q, 0.5);
  SLerp(q1, q, c1, Third);

  SBezier(q0, c0, c1, q1, out, t);
}

void SBezier(const Quaternion& q0,
	     const Quaternion& c0,
	     const Quaternion& c1,
	     const Quaternion& q1,
	     Quaternion& out,
	     Real t)
{
  Quaternion a,b,c,d,e;
  SLerp(q0,c0, a,t);
  SLerp(c0,c1, b,t);
  SLerp(c1,q1, c,t);

  SLerp(a,b, d,t);
  SLerp(b,c, e,t);

  SLerp(d,e, out,t);
}



void SetMatrixRotationZYX(Matrix3& m, const Vector3& v)
{
	EulerAngleRotation e(v);
	e.getMatrixZYX(m);
}

void SetMatrixRotationZYX(Matrix4& m, const Vector3& v)
{
	Matrix3 m3;
	SetMatrixRotationZYX(m3,v);
	m.set(m3);
}

void SetMatrixRotationVector(Matrix3& m, const Vector3& v)
{
	MomentRotation r(v);
	r.getMatrix(m);
}

void SetMatrixRotationVector(Matrix4& m, const Vector3& v)
{
	Matrix3 m3;
	SetMatrixRotationVector(m3,v);
	m.set(m3);
}

void SetMatrixRotationVector(Matrix3& m, const Vector3& v, Real angle)
{
	AngleAxisRotation r(angle,v);
	r.getMatrix(m);
}

void SetMatrixRotationVector(Matrix4& m, const Vector3& v, Real angle)
{
	Matrix3 m3;
	SetMatrixRotationVector(m3,v,angle);
	m.set(m3);
}

void SetMatrixRotationQuaternion(Matrix3& m, const Quaternion& q)
{
	QuaternionRotation r(q);
	r.getMatrix(m);
}

void SetMatrixRotationQuaternion(Matrix4& m, const Quaternion& q)
{
	Matrix3 m3;
	SetMatrixRotationQuaternion(m3,q);
	m.set(m3);
}

} //namespace Math3D
