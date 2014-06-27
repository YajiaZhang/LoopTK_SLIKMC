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
#include "PExtension.h"
#include "PTools.h"

/*
 * Initializes some of the member variables needed before
 * CCD computations can be performed. Notes:
 * 
 *  - len is length between goalEnd and anchorPrior
 *  - angle1 is angle between goalPrior-goalEnd-anchorPrior
 *  - angle2 is angle between goalEnd-anchorPrior-anchorEnd
 */

PCCDSolver::PCCDSolver(PChain *loop, PAtom *effectorPrior, PAtom *effectorEnd, Vector3 goalPrior, Vector3 goalEnd)
{
  m_chain = loop;
  m_DOFs = m_chain->GetDOFs(PID::BACKBONE);
  m_ccdIter = 0;

  m_effectorPrior = effectorPrior;
  m_effectorEnd = effectorEnd;

  m_goalPrior = goalPrior;
  m_goalEnd = goalEnd;
}

/*
 * Returns the goals for this PCCDSolver
 * in the specified references.
 */

void PCCDSolver::GetGoals(Vector3 &goalPrior, Vector3 &goalEnd) {
  goalPrior = m_goalPrior;
  goalEnd = m_goalEnd;
}

void PCCDSolver::UpdateSides(const Vector3 &axisOfRotation, const Vector3 &basePrior, const Vector3 &atom,
        const Vector3 &goal, Real &k1, Real &k2)
{
  Vector3 vecToCur = atom - basePrior, toCircleCenter = project(vecToCur, axisOfRotation),
    circleCenter = basePrior + toCircleCenter, r = atom - circleCenter, g = goal - basePrior, s;
  Real rLength = r.length();
  
  r.inplaceNormalize();
  s = cross(r, axisOfRotation);
  s.inplaceNormalize();

  k1 += rLength * dot(r, g);
  k2 += rLength * dot(s, g);
}

/*
 * Performs one iteration of CCD.
 */

Real PCCDSolver::DoDescent()
{
  if (m_ccdIter >= m_DOFs.size()) {
    m_ccdIter = 0;
  }

  Real k1 = 0, k2 = 0, tanVal, theta, angle3;
  PBond *curBond = m_DOFs[m_ccdIter];
  PAtom *forwardAtom = curBond->getAtom(forward), *backAtom = curBond->getAtom(backward);

  Vector3 axisOfRotation = forwardAtom->getPos() - backAtom->getPos(),
    basePrior = backAtom->getPos();

  UpdateSides(axisOfRotation, basePrior, m_effectorEnd->getPos(), m_goalEnd, k1, k2);
  UpdateSides(axisOfRotation, basePrior, m_effectorPrior->getPos(), m_goalPrior, k1, k2);

  tanVal = k2 / k1;
  theta = atan(tanVal);
  angle3 = RtoD(theta);

  if (k1 * cos(theta) + k2 * sin(theta) < 0) {
    angle3 += 180;
  }

  if (angle3 > 180) {
    angle3 -= 360;
  }

  // check to avoid "nan"
  if (!isNAN(angle3)) {
    m_chain->RotateChain(PID::BACKBONE, m_ccdIter, forward, angle3);
  }

  m_ccdIter++;

  Real d1 = m_effectorEnd->getPos().distance(m_goalEnd), d2 = m_effectorPrior->getPos().distance(m_goalPrior);
  return(d1*d1 + d2*d2);
}
