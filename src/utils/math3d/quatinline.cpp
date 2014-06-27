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

#include "quatinline.h"
#include <math/math.h>
#include <stdio.h>

namespace Math3D {

void quat_slerp(quat_t out, const quat_t a, const quat_t b, Real t)
{
	//a + b unit quaternions?
/* angle theta is defined as
cos(theta)*|a||b| = <a,b>
*/
	Real dot = quat_dot(a,b);
	if(FuzzyEquals(dot,One))		//axes are the same axis
	{
		quat_equal(out,b);
		return;
	}
	else if(FuzzyEquals(dot,-One))	//axes are opposing axis
	{
		fprintf(stderr, "Quaternions on opposing sides of unit sphere\n");
		return;
	}
	Real theta = Acos(dot);
	Real sininv = Sin(theta);
	sininv = Inv(sininv);

	quat_multiply(out, a, Sin((One-t)*theta)*sininv);
	quat_t temp;
	quat_multiply(temp, b, Sin(t*theta)*sininv);
	quat_add(out, out, temp);
}

} //namespace Math3D
