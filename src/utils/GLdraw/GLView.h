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

#ifndef GL_VIEW_H
#define GL_VIEW_H

#include <math3d/primitives.h>
#include <camera/viewport.h>
using namespace Math3D;

//pscreen = V*homogeneous(P*M*pworld)
//pcamera=M*pworld
//pobject=P*M*pworld
class GLView
{
public:
	GLView();
	void setCurrentGL();
	void getCurrentGL();
	void updateInverses();  //must call this after each change of modelview/projection matrices to update the inverse transforms

	void setViewport(const Viewport& v);
	bool getViewport(Viewport& v) const;

	static bool getFrustumMatrix(Real l,Real r,Real b,Real t,Real n,Real f,Matrix4& m);
	static bool getOrthoMatrix(Real l,Real r,Real b,Real t,Real n,Real f,Matrix4& m);

	void worldToScreen(const Vector4& p, Vector3& out);
	void worldToObject(const Vector4& p, Vector4& out);
	void worldToCamera(const Vector4& p, Vector4& out);
	void screenToWorld(const Vector3& p, Vector4& out);
	void objectToWorld(const Vector4& p, Vector4& out);
	void cameraToWorld(const Vector4& p, Vector4& out);

	void pointWorldToObject(const Vector3& p, Vector3& out);
	void pointWorldToCamera(const Vector3& p, Vector3& out);
	void pointObjectToWorld(const Vector3& p, Vector3& out);
	void pointCameraToWorld(const Vector3& p, Vector3& out);
	void vectorWorldToObject(const Vector3& v, Vector3& out);
	void vectorWorldToCamera(const Vector3& p, Vector3& out);
	void vectorObjectToWorld(const Vector3& v, Vector3& out);
	void vectorCameraToWorld(const Vector3& p, Vector3& out);
	void normalWorldToObject(const Vector3& v, Vector3& out);
	void normalWorldToCamera(const Vector3& p, Vector3& out);
	void normalObjectToWorld(const Vector3& v, Vector3& out);
	void normalCameraToWorld(const Vector3& p, Vector3& out);

	Real x,y,w,h;   //viewport
	Matrix4 modelview;
	Matrix4 projection;

	Matrix4 modelviewInverse;
	Matrix4 projectionInverse;
};

#endif
