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

#include "PBasic.h"
#include "PConstants.h"
#include "PIKAlgorithms.h"
#include "PTools.h"

PProteinCCDSolver::PProteinCCDSolver(PProtein *loop, ProteinSide side)
{
  Vector3 goalPrior, goalEnd;
  PTools::GetGoals(loop, side, goalPrior, goalEnd);

  ConstructorHelper(loop, side, goalPrior, goalEnd);
}

PProteinCCDSolver::PProteinCCDSolver(PProtein *loop)
{
  Vector3 goalPrior, goalEnd;
  ProteinSide side = PTools::GetGoals(loop, goalPrior, goalEnd);

  ConstructorHelper(loop, side, goalPrior, goalEnd);
}

void PProteinCCDSolver::ConstructorHelper(PProtein *loop, ProteinSide side, Vector3 goalPrior, Vector3 goalEnd)
{
  PAtom *effectorEnd, *effectorPrior;
  loop->GetEndEffectors(side, effectorEnd, effectorPrior);

  m_solver = new PCCDSolver(loop, effectorPrior, effectorEnd, goalPrior, goalEnd);
}

struct Transformer: AtomFunctor {
public:
  Transformer(Matrix4 toApply) {
    m_toApply = toApply;
    cerr<<toApply<<endl;
  }

  void operator()(PAtom *atom, PBond *bondFrom) {
    atom->ApplyTransform(m_toApply);
  }

private:
  Matrix4 m_toApply;
};

