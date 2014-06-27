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

#ifndef PCONF_SPACE_NAV
#define PCONF_SPACE_NAV

#include "PBasic.h"
#include "PExtension.h"




class PConfSpaceNavigator: public GLUTNavigationProgram {
 public:

   /**
   *  Constructs a PConfSpaceNavigator given a <code>PConformationSpace</code>.
   */
  PConfSpaceNavigator(PConformationSpace *pcs);

   /**
   *  Deletes <code>PConformationSpace</code> and 
   * destroys <code>PConfSpaceNavigator</code>.
   */
  ~PConfSpaceNavigator();

   /**
   * Prepares for rendering.
   */
  virtual bool Initialize();

   /**
   * Renders the conformation space..
   */
  void RenderWorld();
 
  /**
   * Handles the keyboard input.
   */
  void Handle_Keypress(unsigned char key, int x, int y);

 private:
  void DrawLightChain(PLightChain *chain, bool drawOrange);
  void DrawAtom(PAtom *a, bool drawOrange);
  PConformationSpace *m_space;
  bool m_onlyBackbone;
  bool m_onlyLoops;
  PConformationSpace::iterator m_currPos;
  bool m_dispOneOnly;
  int m_sectionPos;
  bool m_sectionDisp;
  int m_quality;

  int TOTAL_BACKBONE_ATOMS;
  

};








#endif
