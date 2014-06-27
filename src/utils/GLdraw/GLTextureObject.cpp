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

#include "GLTextureObject.h"
#include "GL.h"
#include <iostream>
using namespace std;

GLTextureObject::GLTextureObject()
:glName(0)
{}

GLTextureObject::GLTextureObject(GLTextureObject& obj)
{
	stealObject(obj);
}

GLTextureObject::~GLTextureObject()
{
	cleanup();
}

bool GLTextureObject::isNull() const
{
	return glName == 0;
}

void GLTextureObject::stealObject(GLTextureObject& obj)
{
	cleanup();
	glName = obj.glName;
	obj.glName = 0;
}

void GLTextureObject::generate()
{
  if(glName == 0)
    glGenTextures(1, &glName);
  else
    cout<<"Warning, GLTextureObject.generate() called on a non-null object"<<endl;
}

void GLTextureObject::cleanup()
{
	if(glName != 0)
	{
		glDeleteTextures(1, &glName);
		glName = 0;
	}
}

void GLTextureObject::bind(unsigned int target) const
{
	glBindTexture(target,glName);
}

void GLTextureObject::unbind(unsigned int target) const
{
	glBindTexture(target,0);
}
