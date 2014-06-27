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

#include "PExtension.h"
#include "PLibraries.h"

Real PMath::getAngle(const Vector3 &v1, const Vector3 &v2)
{
  Real dotProduct = dot(v1, v2);

  dotProduct /= v1.norm();
  dotProduct /= v2.norm();
  return acos(dotProduct);
}

bool PMath::sphereIntersectsCube(const Vector3 &sphereCenter, Real sphereRadius,
         const Vector3 &cubeVertex1, const Vector3 &cubeVertex2)
{
  Sphere3D sphere;
  AABB3D cube(cubeVertex1, cubeVertex2);

  sphere.center = sphereCenter;
  sphere.radius = sphereRadius;

  return(sphere.intersects(cube));
}

Vector3 PMath::ComputePos(const Vector3 &prev, const Vector3 &cur, Real angle, Real bondLength)
{
  Matrix3 rotMat;
  Vector3 a = cur - prev, b, c, pt;
  Real aLen = a.length(), newALen = bondLength * cos(DtoR(180 - angle)) + aLen;

  a.inplaceNormalize();
  a.inplaceScale(newALen);

  Real perpLength = bondLength * sin(DtoR(180-angle));

  b = ArbitraryNonParallel(a);
  c = cross(a, b);

  c.inplaceNormalize();
  c.inplaceScale(perpLength);

  pt = a + c + prev;
  return pt;
}

Matrix3 PMath::FindRotationMatrix(const Vector3 &axis, Real angle)
{
  Matrix3 rotMat;
  Vector3 scaledAxis(axis);
  Real c = cos(angle), s = sin(angle), omc = 1 - c, len = scaledAxis.length();

  if (len != 1) scaledAxis /= len;

  float x = scaledAxis[0];
  float y = scaledAxis[1];
  float z = scaledAxis[2];
  float xs = x * s;
  float ys = y * s;
  float zs = z * s;
  float xyomc = x * y * omc;
  float xzomc = x * z * omc;
  float yzomc = y * z * omc;

  rotMat(0,0) = x*x*omc + c;
  rotMat(0,1) = xyomc + zs;
  rotMat(0,2) = xzomc - ys;

  rotMat(1,0) = xyomc - zs;
  rotMat(1,1) = y*y*omc + c;
  rotMat(1,2) = yzomc + xs;

  rotMat(2,0) = xzomc + ys;
  rotMat(2,1) = yzomc - xs;
  rotMat(2,2) = z*z*omc + c;

  return rotMat;
}

Vector3 PMath::ArbitraryNonParallel(const Vector3 &v)
{
  Vector3 output = v;

  if (v.x == 0) output.x = 1;
  else if (v.y == 0) output.y = 1;
  else output.z = v.z + 1;

  return output;
}

//negative is clockwise and positive is counterclockwise.
Real PMath::AngleBetweenPlanes(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &p4)
{
  Vector3 u = p1 - p2, v = p4 - p3, w = p3 - p2;
//  return -(180-TorsionAngle(u, v, w));
	return TorsionAngle(u,v,w);
}

Vector3 PMath::LocalToGlobalCoord(const vector<Vector3> definedRelativeCoordinates, const vector<Vector3> definedGlobalCoordinates){
  //setup local coordinate frames
  Vector3 s,t,u,v;
  s.sub(definedRelativeCoordinates[1],definedRelativeCoordinates[0]);
  t.sub(definedRelativeCoordinates[2],definedRelativeCoordinates[0]);
  u.setCross(s,t);
  v.mul(definedRelativeCoordinates[0],-1);
  Real sR,tR,uR; 
  sR = (s.dot(v)) / (s.norm() * v.norm());
  tR = (t.dot(v)) / (t.norm() * v.norm());
  if (u.norm()==0){ uR = 0;}
  else{uR = (u.dot(v)) / (u.norm() * v.norm());}
  
  //setup global coordinate frame.
  Vector3 l,m,n;
  l.sub(definedGlobalCoordinates[1],definedGlobalCoordinates[0]);
  m.sub(definedGlobalCoordinates[2],definedGlobalCoordinates[0]);
  n.setCross(l,m);
  Matrix3 definedGlobalCoordinatesM;
  definedGlobalCoordinatesM.setRow(0,l);
  definedGlobalCoordinatesM.setRow(1,m);
  definedGlobalCoordinatesM.setRow(2,n);
  
  //singularity check.
  definedGlobalCoordinatesM.inplaceInverse();

  //covert local to global.
  Vector3 globalCoordinatesConstants;
  globalCoordinatesConstants.x = sR*l.norm()*v.norm();
  globalCoordinatesConstants.y = tR*m.norm()*v.norm();
  globalCoordinatesConstants.z = uR*n.norm()*v.norm();
  Vector3 toDefineGlobalCoordinates;
  definedGlobalCoordinatesM.mul(globalCoordinatesConstants,toDefineGlobalCoordinates);
  Vector3 finalGlobalCoordinates;
  finalGlobalCoordinates.add(toDefineGlobalCoordinates,definedGlobalCoordinates[0]);  

  return finalGlobalCoordinates;
}

Matrix4 PMath::Find2PtTransform(Vector3 s1, Vector3 s2, Vector3 e1, Vector3 e2) {
  const double eps = .001;
  Matrix4 ret;
  Matrix4 trans;
  Vector3 s = s2-s1;
  Vector3 e = e2-e1;
  trans.setTranslate(e1-s1);
  ret = trans;
  trans.setTranslate(s1);
  ret = ret*trans;
  double angle = acos(s.dot(e)/(s.length()*e.length()));
  Vector3 axis;
  axis.setCross(s,e);
  if (angle<eps&&angle>-eps||angle>Pi-eps&&angle<Pi+eps) {
    axis = ArbitraryNonParallel(s);
  }
  Matrix4 rot(FindRotationMatrix(axis,angle));
  ret = ret*rot;
  trans.setTranslate(-s1);
  ret = ret*trans;
  return ret;
}

void PMath::ConvertToQuaternion(double** input, Vector4* output,vector<int> sign){
  if (sign.size()==0) for(int i=1;i<=4;i++) sign.push_back(1);
  
  double e1,e2,e3,e4;
  e1 = 0.5*sqrt(1+input[1][1]+input[2][2]+input[3][3]);
  e2 = 0.5*sqrt(1+input[1][1]-input[2][2]-input[3][3]);
  e3 = 0.5*sqrt(1-input[1][1]+input[2][2]-input[3][3]);
  e4 = 0.5*sqrt(1-input[1][1]-input[2][2]+input[3][3]);
  
  if ((e1>=e2)&&(e1>=e3)&&(e1>=e4)){
    e1 = sign[0]*e1;
    e2 = (input[3][2]-input[2][3])/(4.0*e1);
    e3 = (input[1][3]-input[3][1])/(4.0*e1);
    e4 = (input[2][1]-input[1][2])/(4.0*e1);
  }
  else if ((e2>=e3)&&(e2>=e4)&&(e2>=e1)){
    e2 = sign[1]*e2;
    e3=(input[2][1]+input[1][2])/(4.0*e2);    
    e4=(input[1][3]+input[3][1])/(4.0*e2);
    e1=(input[3][2]-input[2][3])/(4.0*e2);
  }
  else if ((e3>=e4)&&(e3>=e1)&&(e3>=e2)){
    e3 = sign[2]*e3;
    e4=(input[3][2]+input[2][3])/(4.0*e3);    
    e1=(input[1][3]-input[3][1])/(4.0*e3);
    e2=(input[2][1]+input[1][2])/(4.0*e3);
  }
  else if ((e4>=e1)&&(e4>=e2)&&(e4>=e3)){
    e4 = sign[3]*e4;
    e1=(input[2][1]-input[1][2])/(4.0*e4);    
    e2=(input[1][3]+input[3][1])/(4.0*e4);
    e3=(input[3][2]+input[2][3])/(4.0*e4);
  }
  
  output->set(e1,e2,e3,e4);
}


void PMath::InterpolateQuaternion(Vector4 qi, Vector4 qf, Vector4* qu, double u){
  double cost = qi.x*qf.x+qi.y*qf.y+qi.z*qf.z+qi.w*qf.w;
  double theta = acos(cost);
  //cout<<"sin(theta):"<<sin(theta)<<endl;
  double c1 = sin(theta*(1.0-u))/sin(theta);
  double c2 = sin(theta*u)/sin(theta);
  qu->set(c1*qi.x+c2*qf.x,c1*qi.y+c2*qf.y,c1*qi.z+c2*qf.z,c1*qi.w+c2*qf.w);
}

void PMath::DeltaQuaterToDeltaVel(Vector4 qi, Vector4 qf,Vector3* deltaV){
  double a = 2.0*(-qi.y*(qf-qi).x + qi.x*(qf-qi).y - qi.w*(qf-qi).z + qi.z*(qf-qi).w);
  double b = 2.0*(-qi.z*(qf-qi).x + qi.w*(qf-qi).y + qi.x*(qf-qi).z - qi.y*(qf-qi).w);
  double c = 2.0*(-qi.w*(qf-qi).x - qi.z*(qf-qi).y + qi.y*(qf-qi).z + qi.x*(qf-qi).w);
  deltaV->set(a,b,c);
  
}

int PMath::signum(double a){
  if (a>=0.0) return 1;
  if (a<0.0) return -1;
}

void PMath::normalize(vector<double> &v) {
  double sum = 0;

  // Compute sum
  for (int i = 0; i < v.size(); ++i) {
    sum += v[i];
  }

  // Divide all elements by sum
  for (int i = 0; i < v.size(); ++i) {
    v[i] /= sum;
  }
}

double PMath::TorsionAngle(Vector3 u,Vector3 v,Vector3 w){
  int sign = 1;
  Vector3 destu,destv,destf;
  destu.setCross(w,u);
  destv.setCross(w,v);
  destf.setCross(destu,destv);
  double sinv = destf.norm()/(destu.norm()*destv.norm());
  double cosv = destu.dot(destv)/(destu.norm()*destv.norm());
  destu.setCross(u,v);
  if (destu.dot(w)<0) sign = -1;
  return sign*atan2(sinv,cosv)*180.0f/PI;
}

double PMath::abs(double d) {
	double value=0;
	if (d>=0)
		value = d;
	else
		value = -d;
	return value;
}
