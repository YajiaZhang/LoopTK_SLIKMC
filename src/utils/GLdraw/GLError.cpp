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

#include "GLError.h"
#include <iostream>
using namespace std;

const char* GLErrorString(GLenum err)
{
  switch(err) {
  case GL_NO_ERROR: return "GL_NO_ERROR";
  case GL_INVALID_ENUM: return "GL_INVALID_ENUM";
  case GL_INVALID_VALUE: return "GL_INVALID_VALUE"; 
  case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
  case GL_STACK_OVERFLOW: return "GL_STACK_OVERFLOW";
  case GL_STACK_UNDERFLOW: return "GL_STACK_UNDERFLOW";
  case GL_OUT_OF_MEMORY: return "GL_OUT_OF_MEMORY";
  case GL_TABLE_TOO_LARGE: return "GL_TABLE_TOO_LARGE";
  default: return "GLErrorString(): invalid error code";
  }
}

bool CheckGLErrors(const char* name,bool pause)
{
  bool res=false;
  GLenum err;
  while((err=glGetError()) != GL_NO_ERROR) {
    cout<<name<<" "<<GLErrorString(err)<<endl;
    res=true;
  }
  if(res&&pause)
    getchar();
  return res;
}

