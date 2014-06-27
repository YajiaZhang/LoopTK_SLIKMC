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

#ifndef GLDRAW_GL_ERROR_H
#define GLDRAW_GL_ERROR_H

#include "GL.h"

//Returns the string corresponding to the GL error code err
const char* GLErrorString(GLenum err);

//Checks for GL errors, prints them with the header string name,
//pauses if pause=true.  Returns false if no error occurred.
bool CheckGLErrors(const char* name="GL error",bool pause=true);

#endif
