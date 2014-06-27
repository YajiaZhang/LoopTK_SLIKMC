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

#include "GLMaterial.h"

GLMaterial::GLMaterial()
:specularExponent(0)
{}

void GLMaterial::setCurrentGL(GLenum tgt) const
{
	glMaterialfv(tgt, GL_AMBIENT, ambient.rgba);
	glMaterialfv(tgt, GL_DIFFUSE, diffuse.rgba);
	glMaterialfv(tgt, GL_SPECULAR, specular.rgba);
	glMaterialf(tgt, GL_SHININESS, specularExponent);
	glMaterialfv(tgt, GL_EMISSION, emission.rgba);
}


void GLMaterial::setCurrentGL_Front() const
{
	setCurrentGL(GL_FRONT);
}


void GLMaterial::setCurrentGL_Back() const
{
	setCurrentGL(GL_BACK);
}

void GLMaterial::setCurrentGL_FrontAndBack() const
{
	setCurrentGL(GL_FRONT_AND_BACK);
}
