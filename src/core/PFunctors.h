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

#ifndef __P_FUNCTORS_H
#define __P_FUNCTORS_H

class PAtom;
class PBond;
 //@package Functors
/**
 *
 * A functor to be called on an atom and its
 * source bond (<i>not</i> destination bond) during a
 * graph traversal of atoms and bonds. Create
 * a subclass of <code>AtomFunctor</code> to 
 * specify custom behavior.
 */

class AtomFunctor {
  public:
    virtual void operator()(PAtom *atom, PBond *bondFrom) = 0;
};
 //@package Functors
/**
 *
 * A functor to be called on each bond and its
 * source and destination atoms during a
 * graph traversal of atoms and bonds. Create 
 * a subclass of <code>BondFunctor</code> to
 * specify custom behavior.
 */

class BondFunctor {
  public:
    virtual void operator()(PBond *bond, PAtom *atomFrom, PAtom *atomTo) = 0;
};
 //@package Functors
/**
 *
 * A functor which contains the objective function
 * to be minimized over certain DOFs of loop. Create 
 * a subclass of <code>FunctFunctor</code> to
 * specify custom behavior.
 */


class FunctFunctor {
  public:
    virtual double operator()(double p[]) = 0;
    virtual double operator()(double p[], int i) = 0;
};
 //@package Functors
/**
 *
 * A functor which contains the derivative of the objective function
 * to be minimized over certain DOFs of loop. Create 
 * a subclass of <code>FunctFunctor</code> to
 * specify custom behavior.
 */
class DerivFunctor {
  public:
    virtual void operator()(double p[], double deriv[]) = 0;
};

#endif  // __P_FUNCTORS_H
