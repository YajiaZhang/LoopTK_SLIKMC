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

#ifndef GL_LIGHT_H
#define GL_LIGHT_H

#include "GLColor.h"
#include <math3d/primitives.h>
using namespace Math3D;

class GLLight
{
public:
	GLLight();
	GLLight(const Vector3& direction);
	GLLight(const Vector3& position,const Vector3& direction);
	GLLight(const GLLight& light);

	void setColor(const GLColor& col);
	void setPointLight(const Vector3& pos);
	void setDirectionalLight(const Vector3& dir);
	void setSpotLight(const Vector3& pos,const Vector3& dir,float exponent=0,float cutoff=180);

	void setCurrentGL(const int id=0);

	Vector4 position;
	float att2,att1,att0;	//quadratic,linear,constant attenuation
	GLColor diffuse, specular;
	Vector3 spot_direction;
	float spot_exponent,spot_cutoff;
};

#endif
