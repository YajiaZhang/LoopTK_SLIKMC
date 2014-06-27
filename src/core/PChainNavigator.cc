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

#include <stdio.h>
#include <stdlib.h>
#include "PChainNavigator.h"
#include "PExtension.h"
#include "PConstants.h"
#include "PTools.h"

/* Do we want to show the occupancy map at all (0),
 * only the occupied grid cells (1), or all cells (2)? */
#define    NO_OCCUPANCY_MAP  0
#define    ONLY_OCCUPIED_CELLS  1
#define    ALL_OCCUPANCY_MAP  2

/* Tracks the last pair of atoms in collision. */
pair<PAtom *, PAtom *> lastColl;

/* Color constants for drawing atoms and grid cells. */
const GLColor occupiedCell = GLColor(1, 0, 0, 0.75);
const GLColor unoccupiedCell = GLColor(0, 0, 1, 0.75);
const GLColor sphereColor = GLColor(1, 0, 0, 0.55);
const GLColor outerColor = GLColor(0.8, 0.1, 0, 0.45);
const GLColor innerColor = GLColor(0, 0.1, 0.8, 0.45);
const GLColor pathColor = GLColor(0, 0.5, 0.5, 0.75);
const GLColor nearestPathColor = GLColor(0, 0.5, 0.5, 0.85);
//const double NEAR_DISTANCE = 10.0;
const int MAX_PATH_LENGTH = 10000;
const int MAX_EXTEND_LENGTH = 5;
const int MIN_MOVE_PATH_LENGTH = 6;

void PChainNavigator::ConstructorHelper(PProtein *chain, PProtein *manipulate)
{
  m_chain = chain;
  m_manipulate = manipulate;

  if (m_chain == NULL) {
    PUtilities::AbortProgram("Error: cannot visualize a null chain!");
  }

  if (m_manipulate == NULL) {
    PUtilities::AbortProgram("Error: cannot visualize a null subchain!");
  }

  if (m_chain != m_manipulate) {
    m_ccd = new PProteinCCDSolver(m_manipulate,end);
  } else {
    m_ccd = NULL;
  }

  m_noLoop = false;
  m_onlyBackbone = false;
  m_residueHighlight = false;
  m_index = 0;
  m_side_index = 0;
  m_residue_index = 0;

  m_showProtein = true;
  m_showOccupancy = NO_OCCUPANCY_MAP;

  m_autoTweakTimes = 0;
  m_animateIndex = 0;

  FindBoundingBox();
  FindCenteringTranslate();
}

PChainNavigator::PChainNavigator(PProtein *p) {
  ConstructorHelper(p, p);
}

PChainNavigator::PChainNavigator(PProtein *protein, PProtein *subchain) {
  ConstructorHelper(protein, subchain);
}

struct CenteringGetter: AtomFunctor {
public:
  CenteringGetter() {
    m_vec = Vector3(0,0,0);
    n = 0;
  }

  void operator()(PAtom *atom, PBond *bondFrom) {
    m_vec = m_vec + atom->getPos();
    n++;
  }

  Vector3 GetCenteringTranslate() { return -1*m_vec/n; }

private:
  Vector3 m_vec;
  int n;
};

void PChainNavigator::FindBoundingBox()
{
  m_boundingBox = m_chain->getSpaceManager()->getGlobalBoundingBox();
  m_curX = int(m_boundingBox.second.x) + 1;
  m_curY = int(m_boundingBox.second.y) + 1;
  m_curZ = int(m_boundingBox.second.z) + 1;
}

void PChainNavigator::FindCenteringTranslate() {
  CenteringGetter cg;
  m_manipulate->traverseFromStart(&cg);
  m_centeringTranslate = cg.GetCenteringTranslate();
}

bool PChainNavigator::Initialize() {
  if (!GLUINavigationProgram::Initialize()) return false;

  glClearColor(0,0,0,1);
  lastColl = make_pair<PAtom *, PAtom *>(NULL,NULL);

  glui = GLUI_Master.create_glui_subwindow(main_window, GLUI_SUBWINDOW_BOTTOM);
  glui->set_main_gfx_window(main_window);
}

void PChainNavigator::Handle_Keypress(unsigned char key, int x, int y) {
  int iter, numIter;
  Vector3 MovementVector(0.0,0.1,0.0);
  vector<PProtein*> loops;
  vector<int> Atom;
  vector<Vector3> MovementDirn;
  vector<vector<CDof> > Dofs;
  std::string file_str;

  switch(key) {
    case 'r':
      m_manipulate->RotateChain("backbone",m_index,forward,5);
      break;
    case 'e':
      m_manipulate->RotateChain("backbone",m_index,forward,-5);
      break;
    case 'a':
      if (m_ccd!=NULL)
      cout << "Sum of square distances is " << m_ccd->DoDescent() << endl;
      break;
    case 'c':
      Real curSum;
      numIter = 0;
      if (m_ccd!=NULL) {
        do {
          numIter++;
          curSum = m_ccd->DoDescent();
          cout << "Iteration " << numIter << ": sum of squares is " << curSum << endl;
        } while(curSum >= 0.001);
      }
      break;
    case '>':
      m_index++;
      break;
    case '<':
      m_index--;
      break;
    case 'b':
      m_onlyBackbone = !m_onlyBackbone;
      break;
    case '[':
      m_index=0;
      break;
    case ']':
      m_index = m_manipulate->NumDOF("backbone")-1;
      break;
    case 's':
      PDBIO::writeToFile(m_chain, "pdbfiles/output.pdb");
      break;
    case 'v':
      m_manipulate->RandomizeDOFs(forward);
      break;
    case 'D':
      m_manipulate->DisableSidechains();
      break;
    case 'E':
      m_manipulate->EnableSidechains();
      break;
    case 'l':
      m_noLoop = !m_noLoop;
      break;

  /*rotate side chain; applpy rotamer options*/
    case '1':
      m_side_index=0;
      break;
    case '2':
      m_side_index++;
      break;
    case '3':
      m_side_index--;
      break;
    case '4':
      m_manipulate->RotateChain(PID::SIDECHAIN,m_side_index,forward,5);
      break;
    case '5':
      m_manipulate->RotateChain(PID::SIDECHAIN,m_side_index,forward,-5);
      break;
    case '6':
      m_residue_index=0;
      break;
    case '7':
      m_residue_index++;
      break;
    case '8':
      m_residue_index--;
      break;
    case '9':
      m_manipulate->getResidue(m_residue_index)->ApplyRotamer();
      break;
    case '0':
      m_manipulate->getResidue(m_residue_index)->ResetSideChain();
      break;
    case '/':
      m_residueHighlight=!m_residueHighlight;
      break;
  }

  //boundary for index.
  if (m_index<0) m_index =0;
  if (m_residue_index<0) m_residue_index =0;
  if (m_side_index<0) m_side_index =0;

  if (m_index>=m_manipulate->NumDOF("backbone")) m_index = m_manipulate->NumDOF("backbone")-1;
  if (m_side_index>=m_manipulate->NumDOF(PID::SIDECHAIN)) m_side_index = m_manipulate->NumDOF(PID::SIDECHAIN)-1;
  if (m_residue_index>=m_manipulate->size()) m_residue_index = m_manipulate->size()-1;

  Refresh();
}

void PChainNavigator::RenderWorld()
{
  glPushMatrix();

  CenterWorld();
  RenderLight();

  if (m_showProtein) {
    RenderProtein();
  }

  if (m_showOccupancy != NO_OCCUPANCY_MAP) {
    RenderOccupancy();
  }

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glPopMatrix();
}

void PChainNavigator::CenterWorld()
{
  glTranslatef(m_centeringTranslate[0], m_centeringTranslate[1], m_centeringTranslate[2]);
}

void PChainNavigator::RenderLight()
{
  GLLight light;

  light.setDirectionalLight(Vector3(0,0,1));
  light.setCurrentGL();
  glEnable(GL_LIGHTING);
}

void PChainNavigator::RenderProtein()
{
  pair<PAtom *, PAtom *> coll = m_chain->FindAnyCollision();
  if (coll.first != NULL) lastColl = coll;


  for(int j = 0; j < m_chain->size(); j++) {
    PResidue *res = m_chain->getResidue(j);
    vector<PAtom *> *atoms = res->getAtoms();
    for(int i = 0; i < atoms->size(); i++) {
      PAtom *a = (*atoms)[i];
      if ((m_noLoop==false || m_manipulate==NULL || a->getChain()->IsSubChainOf(m_manipulate))&&(a->getParentBlock()->getType()=="backbone"||!m_onlyBackbone)) {
  glPushMatrix();

  /* Set the current sphere's position. */
  Vector3 pos = a->getPos();
  Real x, y, z;
  pos.get(x, y, z);
  glTranslated(x, y, z);

  /* Set the current sphere's color. */
  GLMaterial mat;
  mat.emission.setBlack();
  mat.specular.setWhite();
  mat.specularExponent = 10;

  if (a != coll.first && a != coll.second) {
    if (m_manipulate==NULL || a->getChain()->IsSubChainOf(m_manipulate)) {
      mat.ambient = a->getColor();  /* Use the atom's color to draw. */
      mat.diffuse = a->getColor();
    } else {
      mat.ambient.set(1,.5,0);
      mat.diffuse.set(1,.5,0);
    }
  } else {
    mat.ambient.setGreen();    /* Use the collision color (green). */
    mat.diffuse.setGreen();
  }
  mat.setCurrentGL();
  drawSphere(a->getCovalentRadius(), 5, 5);

  /* Highlight current residue */
  if (j==m_residue_index&&m_residueHighlight){
    //mat.ambient.set(1,.5,0);
    mat.diffuse.set(1,.5,0,0.5);
    mat.setCurrentGL();
    drawSphere(a->getCovalentRadius(),10 , 10);
  }

  glPopMatrix();
      }
    }
  }
}

void PChainNavigator::RenderOccupancy()
{
  OccupancyMap occupancy = m_chain->getSpaceManager()->getOccupancy();
  pair<Vector3, Vector3> curBox;

  for(OccupancyMap::const_iterator it = occupancy.begin(); it != occupancy.end(); ++it) {
    /* Do we only care about visualizing occupied cells? */
    if (m_showOccupancy == ONLY_OCCUPIED_CELLS && !it->second) continue;

    /* Do we only care about visualizing one cross-section? */
    if (m_curX != int(m_boundingBox.second.x) + 1 && m_curX != int(it->first.x)) continue;
    if (m_curY != int(m_boundingBox.second.y) + 1 && m_curY != int(it->first.y)) continue;
    if (m_curZ != int(m_boundingBox.second.z) + 1 && m_curZ != int(it->first.z)) continue;

    glPushMatrix();

    /* Get the coordinates of the bounding box vertices. */
    curBox = m_chain->getSpaceManager()->getCellBoundingBox(it->first);

    /* Set up the material to be drawn. */
    GLMaterial mat;
    mat.emission.setBlack();
    mat.specular.setWhite();
    mat.specularExponent = 20;

    if (it->second) {
      mat.ambient = mat.diffuse = occupiedCell;
    } else {
      mat.ambient = mat.diffuse = unoccupiedCell;
    }

    /* Install the material and draw the bounding box. */
    mat.setCurrentGL();
    drawBoundingBox(curBox.first, curBox.second);
    drawWireBoundingBox(curBox.first, curBox.second);

    glPopMatrix();
  }
}
