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

#ifndef CAMERA_VIEWPORT_H
#define CAMERA_VIEWPORT_H

#include "camera.h"

class Viewport : public Camera
{
public:
	Viewport();
	void setPerspective(bool);
	void setLensAngle(float rads);

	bool clicked(int mx, int my) const;
	void zoom(float s);
	void scroll(float x, float y, float z=0);

	void getViewVector(Vector3&) const;
	void getMovementVector(float x, float y, Vector3&) const;
	void getMovementVectorAtDistance(float x, float y, float dist, Vector3&) const;
	
	void getClickVector(float mx, float my, Vector3&) const;
	void getClickSource(float mx, float my, Vector3&) const;

	bool perspective;
	float scale;
 
	int x,y,w,h;
	float n, f;
};

std::ostream& operator << (std::ostream& out, const Viewport&);
std::istream& operator >> (std::istream& in, Viewport&);

#endif
