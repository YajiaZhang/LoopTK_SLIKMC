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

#include "GLLight.h"
#include "GL.h"

inline void glLightv(GLenum light, GLenum pname,const GLfloat* v)
{
  glLightfv(light,pname,v);
}

inline void glLightv(GLenum light, GLenum pname,const GLdouble* v)
{
  GLfloat vf[4];
  vf[0]=(GLfloat)v[0];
  vf[1]=(GLfloat)v[1];
  vf[2]=(GLfloat)v[2];
  vf[3]=(GLfloat)v[3];
  glLightfv(light,pname,vf);
}

GLLight::GLLight()
:position(Zero),att2(0),att1(0),att0(1),
spot_direction(Zero),spot_exponent(0),spot_cutoff(0)
{}

GLLight::GLLight(const Vector3& direction)
:att2(0),att1(0),att0(1)
{
	setDirectionalLight(direction);
}

GLLight::GLLight(const Vector3& position,const Vector3& direction)
:att2(0),att1(0),att0(1)
{
	setSpotLight(position,direction);
}

GLLight::GLLight(const GLLight& light)
:position(light.position),att2(light.att2),att1(light.att2),att0(light.att0),
diffuse(light.diffuse), specular(light.diffuse),
spot_direction(light.spot_direction),spot_exponent(light.spot_exponent),spot_cutoff(light.spot_cutoff)
{
}

void GLLight::setColor(const GLColor& col) { diffuse=specular=col; }

void GLLight::setPointLight(const Vector3& pos)
{
	position.set(pos.x,pos.y,pos.z,1);
	spot_direction.setZero();
	spot_exponent=0;
	spot_cutoff=180;
}

void GLLight::setDirectionalLight(const Vector3& dir)
{
	position.set(dir.x,dir.y,dir.z,0);
	spot_direction.setZero();
	spot_exponent=0;
	spot_cutoff=180;
}

void GLLight::setSpotLight(const Vector3& pos,const Vector3& dir,float exponent,float cutoff)
{
	position.set(pos.x,pos.y,pos.z,1);
	spot_direction=dir;
	assert(exponent>=0 && exponent<=128); //max GL range
	assert(cutoff>=0 && (cutoff<=90 || cutoff==180)); //max GL range
	spot_exponent=exponent;
	spot_cutoff=cutoff;
}

void GLLight::setCurrentGL(const int id)
{
	GLenum gl_light_id=GL_LIGHT0+id;
	glLightv(gl_light_id,GL_DIFFUSE,diffuse.rgba);
	glLightv(gl_light_id,GL_SPECULAR,specular.rgba);
	glLightv(gl_light_id,GL_POSITION,position);
	glLightf(gl_light_id,GL_QUADRATIC_ATTENUATION,att2);
	glLightf(gl_light_id,GL_LINEAR_ATTENUATION,att2);
	glLightf(gl_light_id,GL_CONSTANT_ATTENUATION,att0);
	glLightv(gl_light_id,GL_SPOT_DIRECTION,spot_direction);
	glLightf(gl_light_id,GL_SPOT_EXPONENT,spot_exponent);
	glLightf(gl_light_id,GL_SPOT_CUTOFF,spot_cutoff);
	glEnable(gl_light_id);
}
