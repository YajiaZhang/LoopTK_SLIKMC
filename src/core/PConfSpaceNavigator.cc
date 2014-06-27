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

#include "PConfSpaceNavigator.h"
#include "PConstants.h"

PConfSpaceNavigator::PConfSpaceNavigator(PConformationSpace *pcs) {
  m_space = pcs;
  m_onlyBackbone = false;
  m_onlyLoops = false;
  m_currPos = m_space->begin();
  m_dispOneOnly = false;
  m_sectionPos = 0;
  m_sectionDisp = false;
  m_quality = 4;
  TOTAL_BACKBONE_ATOMS = m_space->front()->size()*3;
}

PConfSpaceNavigator::~PConfSpaceNavigator() {
  delete m_space;
}

bool PConfSpaceNavigator::Initialize() {
  glClearColor(0,0,0,1);
  return GLUTNavigationProgram::Initialize();
}

void PConfSpaceNavigator::DrawAtom(PAtom *a, bool drawOrange) {
  if (a->getParentBlock()->getType()==PID::BACKBONE||!m_onlyBackbone) {
    glPushMatrix();
    Vector3 pos = a->getPos();
    Real x,y,z;
    pos.get(x,y,z);
    glTranslated(x,y,z);
    GLMaterial mat;
    mat.emission.setBlack();
    mat.specular.setWhite();
    mat.specularExponent = 10;
    if (drawOrange) {
      mat.ambient.set(1,0.5,0);    
      mat.diffuse.set(1,0.5,0);
    } else {
      mat.ambient = a->getColor();
      mat.diffuse = a->getColor();
    }
    mat.setCurrentGL();
    drawSphere(a->getCovalentRadius(),m_quality,m_quality);
    glPopMatrix();
  }
}

void PConfSpaceNavigator::DrawLightChain(PLightChain *chain, bool drawOrange) {
  if (!m_sectionDisp||drawOrange) {
    for(int i=0;i<chain->size();i++) {
      PResidue *res = chain->getResidue(i);
      vector<PAtom *> *atoms = res->getAtoms();
      for(int i=0;i<atoms->size();i++) {
  PAtom *a = (*atoms)[i];
  DrawAtom(a,drawOrange);
      }
    }
  } else {
    int resNum = m_sectionPos/3;
    PResidue *res = chain->getResidue(resNum);
    int atomNum = m_sectionPos%3;
    PAtom *a;
    switch(atomNum) {
    case 0:
      a = res->getAtom(PID::N);
      break;
    case 1:
      a = res->getAtom(PID::C_ALPHA);
      break;
    case 2:
      a = res->getAtom(PID::C);
      break;
    }
    DrawAtom(a,false);
  }
}

void PConfSpaceNavigator::RenderWorld() {
  GLLight light;
  light.setDirectionalLight(Vector3(0,0,1));
  light.setCurrentGL();
  glEnable(GL_LIGHTING);
  if (!m_onlyLoops) {
    DrawLightChain(m_space->getStaticPortions().first,true);
    DrawLightChain(m_space->getStaticPortions().second,true);
  }
  if (!m_dispOneOnly) {
    for(PConformationSpace::iterator it = m_space->begin();it!=m_space->end();++it) {
      DrawLightChain(*it,false);
    }
  } else {
    DrawLightChain(*m_currPos,false);
  }

}

void PConfSpaceNavigator::Handle_Keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'b':
    m_onlyBackbone = !m_onlyBackbone;
    break;
  case 'l':
    m_onlyLoops = !m_onlyLoops;
    break;
  case '-':
    if (m_currPos==m_space->begin())
      m_currPos = m_space->end();
    m_currPos--;
    break;
  case '=':
    m_currPos++;
    if (m_currPos==m_space->end())
      m_currPos = m_space->begin();
    break;
  case 'd':
    m_dispOneOnly = !m_dispOneOnly;
    break;
  case 's':
    m_sectionDisp = !m_sectionDisp;
    break;
  case '[':
    m_sectionPos = (m_sectionPos+TOTAL_BACKBONE_ATOMS-1)%TOTAL_BACKBONE_ATOMS;
    break;
  case ']':
    m_sectionPos = (m_sectionPos+1)%TOTAL_BACKBONE_ATOMS;
    break;
  case '+':
    m_quality++;
    break;
  case '_':
    m_quality--;
    if (m_quality==1) m_quality = 2;
    break;
  }
  Refresh();
}
