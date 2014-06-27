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

#ifndef PLIBRARIES_H
#define PLIBRARIES_H

/*
 * Basic libraries from the STL.
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <set>
using namespace std;

/*
 * Advanced libraries from the STL.
 */

#include <ext/hash_set>
#include <ext/hash_map>
using namespace __gnu_cxx;

/*
 * File/system libraries.
 */

#include <sys/types.h> 
#include <sys/stat.h>
#include <sys/param.h> 
#include <dirent.h>
#include <unistd.h>

/*
 * Math and OpenGL libraries.
 */

#include <math3d/primitives.h>
#include <math3d/rotation.h>
#include <math.h>
#include <GLdraw/GLUTNavigationProgram.h>
#include <GLdraw/GLUINavigationProgram.h>
#include <GLdraw/GLLight.h>
#include <GLdraw/GLMaterial.h>
#include <GLdraw/drawextra.h>
#include <GLdraw/GL.h>
#include <GL/gl.h>
#include <glui.h>
using namespace Math3D;

/*
 * Other program libraries.
 */

#include "PHashing.h"
#include "PUtilities.h"

#endif
