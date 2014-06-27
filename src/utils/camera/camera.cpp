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

#include "camera.h"
#include <math3d/rotation.h>
#include <iostream>
using namespace std;

void Camera::GetOrientationMatrix(Orientation o,Matrix3& mat)
{
  switch(o) {
  case XYZ:
    mat.setIdentity();
    break;
  case XYnZ:
    mat.setIdentity();
    mat(2,2) = -One;
    break;
  case XZY:
    mat.setZero();
    mat(0,0) = mat(1,2) = mat(2,1) = One;
    break;
  case XZnY:
    mat.setZero();
    mat(0,0) = mat(1,2) = One;
    mat(2,1) = -One;
    break;
  default:
    cerr<<"Unknown orientation to Camera::GetOrientationMatrix"<<endl;
    abort();
    break;
  }
}

void Camera::Orient(Orientation o,Matrix3& mat)
{
  Vector3 ry,rz;
  switch(o) {
  case XYZ:
    break;
  case XYnZ:
    mat.getRow3(rz);
    rz.inplaceNegative();
    mat.setRow3(rz);
    break;
  case XZY:
    mat.getRow2(ry);
    mat.getRow3(rz);
    mat.setRow3(ry);
    mat.setRow2(rz);
    break;
  case XZnY:
    mat.getRow2(ry);
    mat.getRow3(rz);
    rz.inplaceNegative();
    mat.setRow3(ry);
    mat.setRow2(rz);
    break;
  default:
    cerr<<"Unknown orientation to Camera::Orient"<<endl;
    abort();
    break;
  }
}

void Camera::Unorient(Orientation o,Matrix3& mat)
{
  Vector3 ry,rz;
  switch(o) {
  case XYZ:
    break;
  case XYnZ:
    mat.getRow3(rz);
    rz.inplaceNegative();
    mat.setRow3(rz);
    break;
  case XZY:
    mat.getRow2(ry);
    mat.getRow3(rz);
    mat.setRow3(ry);
    mat.setRow2(rz);
    break;
  case XZnY:
    mat.getRow2(ry);
    mat.getRow3(rz);
    ry.inplaceNegative();
    mat.setRow3(ry);
    mat.setRow2(rz);
    break;
  default:
    cerr<<"Unknown orientation to Camera::Unorient"<<endl;
    abort();
    break;
  }
}

void Camera::setFree(const Vector3& pos, const Vector3& rot, Orientation o)
{
  //rot is (pitch yaw roll) = (x y z)
  //world->camera transform is Rz(roll)Rx(pitch)Ry(yaw)
  //so inverse is or Ry(yaw)Rx(pitch)Rz(roll)
  xform.t = pos;
  EulerAngleRotation r(rot.y,rot.x,rot.z);
  r.getMatrixYXZ(xform.R);
  Orient(o,xform.R);
}

void Camera::setOrbit(const Vector3& rot, const Vector3& target, Real dist, Orientation o)
{
  EulerAngleRotation r(rot.y,rot.x,rot.z);
  r.getMatrixYXZ(xform.R);
  Orient(o,xform.R);

  Vector3 zb;
  xform.R.getCol3(zb);
  xform.t = target + dist*zb;
}

void Camera::setTarget(const Vector3& pos, const Vector3& tgt, const Vector3& up)
{
  xform.t = pos;
  Vector3 z = tgt-pos;
  normalize(z);
  Vector3 x;
  x.setCross(z,up);
  normalize(x);
  Vector3 y;
  y.setCross(z,x);
  xform.R.set(x,y,z);
}

void Camera::setCameraMatrix(const Matrix4& m)
{
  Matrix4 minv;
  if(!minv.setInverse(m)) {
    cerr<<"Camera modelview matrix not invertible"<<endl;
    return;
  }
  xform.set(minv);
}

void Camera::getCameraMatrix(Matrix4& mat) const
{
  RigidTransform inv;
  inv.setInverse(xform);
  inv.get(mat);
}



void Camera::getFree(Vector3& pos, Vector3& rot, Orientation o) const
{
  Matrix3 R=xform.R;
  Unorient(o,R);
  EulerAngleRotation r;
  r.setMatrixYXZ(R);
  rot.set(r.y,r.x,r.z);
  pos = xform.t;
}

void Camera::getOrbit(Vector3& rot, Vector3& target, Real tgtdist, Orientation o) const
{
  Matrix3 R=xform.R;
  Unorient(o,R);
  EulerAngleRotation r;
  r.setMatrixYXZ(R);
  rot.set(r.y,r.x,r.z);

  Vector3 z;
  xform.R.getCol3(z);
  target = xform.t - z*tgtdist;
}

void Camera::getTarget(Vector3& pos, Vector3& tgt, Vector3& up, Real tgtDist) const
{
  //y dir
  xform.R.getCol2(up);

  //position
  pos = xform.t;

  //z dir
  Vector3 z;
  xform.R.getCol3(z);

  tgt = pos + z*tgtDist;
}


void CameraController_Free::toCamera(Camera& c) const
{
  c.setFree(pos,rot,ori);
}

void CameraController_Free::fromCamera(const Camera& c)
{
  c.getFree(pos,rot,ori);
}

void CameraController_Target::toCamera(Camera& c) const
{
  c.setTarget(pos,tgt,up);
}

void CameraController_Target::fromCamera(const Camera& c, Real tgtdist)
{
  c.getTarget(pos,tgt,up,tgtdist);
}

void CameraController_Orbit::toCamera(Camera& c) const
{
  c.setOrbit(rot,tgt,dist,ori);
}

void CameraController_Orbit::fromCamera(const Camera& c, Real tgtdist)
{
  c.getOrbit(rot,tgt,tgtdist,ori);
  dist = tgtdist;
}
