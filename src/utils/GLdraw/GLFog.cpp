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

#include "GLFog.h"
#include "GL.h"

GLFog::GLFog()
:type(None),start(0),end(1),density(1)
{
}

void GLFog::setCurrent()
{
	glFogfv(GL_FOG_COLOR, color.rgba);
	switch(type) {
	case None:
		glDisable(GL_FOG);
		break;
	case Linear:
		glEnable(GL_FOG);
		glFogi(GL_FOG_MODE, GL_LINEAR);
		glFogf(GL_FOG_START, start);
		glFogf(GL_FOG_END, end);
		break;
	case Exponential:
		glEnable(GL_FOG);
		glFogi(GL_FOG_MODE, GL_EXP);
		glFogf(GL_FOG_DENSITY, density);
		break;
	case ExponentialSquared:
		glEnable(GL_FOG);
		glFogi(GL_FOG_MODE, GL_EXP2);
		glFogf(GL_FOG_DENSITY, density);
		break;
	}
}

void GLFog::unsetCurrent()
{
	glDisable(GL_FOG);
}

void GLFog::setLinear(Real s, Real e)
{
	type=Linear;
	start=s;
	end=e;
}

void GLFog::setExponential(Real d)
{
	type=Exponential;
	density=d;
}

void GLFog::setExponentialSquared(Real d)
{
	type=ExponentialSquared;
	density=d;
}

