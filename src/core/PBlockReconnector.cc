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
#include "PExtension.h"

PBlockReconnector::PBlockReconnector(PBlock *anchor, PBlock *toDetach) {
  //store anchor PBlock and detached PBlock
  m_anchor = anchor;
  m_detached = toDetach;
  
  //define local coordinate using anchor atoms.
  vector<Vector3> anchorPos;
  for(HASH_MAP_STR(PAtom *)::iterator it = m_anchor->m_atoms.begin(); it!=m_anchor->m_atoms.end(); ++it) {
    anchorPos.push_back(it->second->getPos());
  }

  m_originalCA = anchorPos[0];

  m_xOrig =  anchorPos[1] - anchorPos[0];
  m_xOrig.inplaceNormalize();
  Vector3 temp =  anchorPos[2] - anchorPos[0];
  temp.inplaceNormalize();  
  m_yOrig.setCross(m_xOrig,temp);
  m_yOrig.inplaceNormalize();
  m_zOrig.setCross(m_xOrig,m_yOrig);
  m_zOrig.inplaceNormalize();
  
  //define local coord to global coord transform; define its inverse
  m_transformAW.setCol(0,Vector4(m_xOrig));
  m_transformAW.setCol(1,Vector4(m_yOrig));
  m_transformAW.setCol(2,Vector4(m_zOrig));
  m_transformAW.setCol(3,Vector4(anchorPos[0]));
  m_transformAW.setRow(3,Vector4(0,0,0,1));
  m_transformAW.getInverse(m_transformWA);

}

void PBlockReconnector::ReconnectBlocks() {
  Matrix4 transform = GetReconnectTransform();
  m_detached->ApplyTransform(transform);
}

Matrix4 PBlockReconnector::GetReconnectTransform() {

  //define new local coordinate using anchor atoms.
  vector<Vector3> anchorPos;
  for(HASH_MAP_STR(PAtom *)::iterator it = m_anchor->m_atoms.begin(); it!=m_anchor->m_atoms.end(); ++it) {
    anchorPos.push_back(it->second->getPos());
  }

  m_xPr = anchorPos[1] - anchorPos[0];
  m_xPr.inplaceNormalize();
  Vector3 temp =  anchorPos[2] - anchorPos[0];
  temp.inplaceNormalize();  
  m_yPr.setCross(m_xPr,temp);
  m_yPr.inplaceNormalize();
  m_zPr.setCross(m_xPr,m_yPr);
  m_zPr.inplaceNormalize();

  Vector3 p = anchorPos[0] - m_originalCA;

  //define original local coord to new local coord transform
  m_transformAB.setRow(0,Vector4(m_xPr.dot(m_xOrig),m_yPr.dot(m_xOrig),m_zPr.dot(m_xOrig),p.dot(m_xOrig)));
  m_transformAB.setRow(1,Vector4(m_xPr.dot(m_yOrig),m_yPr.dot(m_yOrig),m_zPr.dot(m_yOrig),p.dot(m_yOrig)));
  m_transformAB.setRow(2,Vector4(m_xPr.dot(m_zOrig),m_yPr.dot(m_zOrig),m_zPr.dot(m_zOrig),p.dot(m_zOrig)));
  m_transformAB.setRow(3,Vector4(0,0,0,1));

  //multiply transforms to get transform of anchors in global coordinate.
  Matrix4 transformWB,transformW;
  transformWB.mul(m_transformAB,m_transformWA);
  transformW.mul(m_transformAW,transformWB);

  return transformW;

}
