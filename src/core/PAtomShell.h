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

#ifndef __P_ATOM_SHELL_H
#define __P_ATOM_SHELL_H

#include <string>
using std::string;

#include <math3d/primitives.h>
#include <GLdraw/GLMaterial.h>
using namespace Math3D;
//@package Main Infrastructure
/**
 *
 *
 * Encapsulates data common to all atoms
 * of the same element.
 */

class PAtomShell {
 public:

  /**
   * Constructs a new PAtomShell.  Caller must specify both
   * the covalent and van der Waals radius, in Angstroms; the
   * <code>name</code> of the atom; and the <code>color</code>
   * to be used when drawing the atom.
   */

  PAtomShell(Real covalentRadius, Real vanDerWaalsRadius, const string &name, GLColor color);

  /**
   * Returns the covalent radius of this atom.
   */

  Real getCovalentRadius() const { return m_covalentRadius; }

  /**
   * Returns the van der Waals radius of this atom.
   */

  Real getVanDerWaalsRadius() const { return m_vanDerWaalsRadius; }

  /**
   * Returns the name of this atom.
   */

  string getName() const { return m_name; }

  /**
   * Returns the color to be drawn when this atom
   * is displayed graphically.
   */

  GLColor getColor() const { return m_color; }

 private:
  Real m_covalentRadius, m_vanDerWaalsRadius;
  string m_name;
  GLColor m_color;
};

#endif  // __P_ATOM_SHELL_H
