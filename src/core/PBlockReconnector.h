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

#ifndef __P_BLOCK_RECONNECTOR_H
#define __P_BLOCK_RECONNECTOR_H

/**
 * <code>PBlockReconnector</code> is used to generate a 
 * homogeneous transformation matrix needed to reconnect a detached PBlock.
 */

class PBlockReconnector {
  public:

    /**
     * Constructs a new <code>PBlockReconnector</code>. Required
     * are: the PBlock <code>anchor</code> and the PBlock <code>toDetach</code> 
     */
    PBlockReconnector(PBlock *anchor, PBlock *toDetach);
    
    /**
     * Reconnects the detached PBlock <code>toDetach</code> to the 
     * PBlock <code>anchor</code> (in the new position).
     */
    void ReconnectBlocks();

    /**
     * Returns the homogeneous transformation matrix that realigns the atoms.
     */
    Matrix4 GetReconnectTransform();

    /**
     * Returns the anchor <code>PBlock</code>.
     */
    PBlock *getAnchor() {return m_anchor; }

    /**
     * Returns the detached <code>PBlock</code>.
     */
    PBlock *getDetached() {return m_detached; }

  private:
    PBlock *m_anchor;
    PBlock *m_detached;

    Matrix4 m_transformWA,m_transformAB,m_transformAW;

    Vector3 m_originalCA;
    Vector3 m_xOrig,m_yOrig,m_zOrig;
    Vector3 m_xPr,m_yPr,m_zPr;

};

#endif  // __P_BLOCK_RECONNECTOR_H
