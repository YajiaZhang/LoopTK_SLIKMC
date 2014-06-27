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

#ifndef GL_DRAWEXTRA_H
#define GL_DRAWEXTRA_H


#include "GL.h"
#include <math3d/primitives.h>
using namespace Math3D;

void drawPoint(const Vector3& pt);
void drawLineSegment(const Vector3& a,const Vector3& b);
void drawTriangle(const Vector3& a,const Vector3& b,const Vector3& c);
void drawQuad(const Vector3& a,const Vector3& b,const Vector3& c,const Vector3& d);

void drawWireCircle2D(const Vector2& center,float r,int numIncrements=16);
void drawWireCircle(const Vector3& axis,float r, int numIncrements=16);
void drawWireSphere(float r, int numIncrements=16);
void drawWireArc(float r, const Vector3& axis, const Vector3& dir, float min, float max);
void drawOrientedWireSphere(float r, const Matrix4& basis);
void drawWirePyramid(const Vector3& c, const Matrix4& basis, float l, float w, float h);
void drawWireBox(float l, float w, float h);
void drawWireBoxCorner(float l, float w, float h);
void drawOrientedWireBox(float l, float w, float h, const Matrix4& basis);
void drawWireBoundingBox(const Vector3& bmin,const Vector3& bmax);
void drawWireOctahedron(float l, float w, float h);
void drawWireCone(const Vector3& h, float r, int numSteps=16);
void drawWireConeFlipped(const Vector3& h, float r, int numSteps=16);

void drawCircle2D(const Vector2& center,float r,int numIncrements=16);
void drawCircle(const Vector3& axis,float r, int numIncrements=16);
void drawSphere(float r, int numSlices, int numStacks);
void drawPyramid(const Vector3& c, const Matrix4& basis, float l, float w, float h);
void drawBox(float l, float w, float h);
void drawBoxCorner(float l, float w, float h);
void drawOrientedBox(float l, float w, float h, const Matrix4& basis);
void drawBoundingBox(const Vector3& bmin,const Vector3& bmax);
void drawOctahedron(float l, float w, float h);
void drawCone(const Vector3& h, float r, int numSteps=16);
void drawConeFlipped(const Vector3& h, float r, int numSteps=16);

void drawCoords(float len);
void drawXYGrid(int n, float spacing);
void drawXZGrid(int n, float spacing);
void drawXYCheckerboard(int n, float spacing, float col1[4],float col2[4]);


//overloaded aliases for float/double objects

inline void glNormal3v(const GLfloat* v)
{
	glNormal3fv(v);
}

inline void glNormal3v(const GLdouble* v)
{
	glNormal3dv(v);
}

inline void glVertex3v(const GLfloat* v)
{
	glVertex3fv(v);
}

inline void glVertex3v(const GLdouble* v)
{
	glVertex3dv(v);
}

inline void glLoadMatrix(const GLfloat* m)
{
	glLoadMatrixf(m);
}

inline void glLoadMatrix(const GLdouble* m)
{
	glLoadMatrixd(m);
}

inline void glMultMatrix(const GLfloat* m)
{
	glMultMatrixf(m);
}

inline void glMultMatrix(const GLdouble* m)
{
	glMultMatrixd(m);
}

inline void glTranslate(const GLfloat* t)
{
	glTranslatef(t[0],t[1],t[2]);
}

inline void glTranslate(const GLdouble* t)
{
	glTranslated(t[0],t[1],t[2]);
}

#endif
