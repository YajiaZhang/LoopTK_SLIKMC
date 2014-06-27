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

#ifndef PCHAIN_NAV
#define PCHAIN_NAV

#include "PBasic.h"
#include "PResources.h"
#include "PIKAlgorithms.h"
#include "PLibraries.h"
#include "PExtension.h"

class PChainNavigator: public GLUINavigationProgram {
  public:

   /**
   *  Constructs a PChainNavigator for a given protein p.
   */
    PChainNavigator(PProtein *p);
 
   /**
   *  Constructs a PChainNavigator for a given protein and a subchain.
   */
   PChainNavigator(PProtein *protein, PProtein *subchain);
    
   /**
   * Destroys the PChainNavigator, <code>PProtein</code>, 
   * and frees <code>PResources</code>.
   */
    ~PChainNavigator() {
      delete m_chain;
      PResources::FreeResources();
    }

   /**
   * Prepares for rendering.
   */
    virtual bool Initialize();

   /**
   * Renders the protein.
   */
    void RenderWorld();

   /**
   * Handles the keyboard input.
   */
    void Handle_Keypress(unsigned char key, int x, int y);

  protected:
    void ConstructorHelper(PProtein *chain, PProtein *manipulate);
    void FindCenteringTranslate();
    void FindBoundingBox();
    
    /*
     * This method will move current atom path inward, so please press
     * 'I', 'O', 'R' first to see the effect.
     */
    void MovePathInward();
    
    /* Drawing helper methods. */
    void CenterWorld();
    void RenderLight();
    void RenderProtein();
    void RenderOccupancy();

    /* Private member variables. */
    GLUI *glui;

    PProtein *m_chain;
    PProtein *m_manipulate;

    pair<Vector3, Vector3> m_boundingBox;
    int m_curX, m_curY, m_curZ;

    bool m_onlyBackbone;
    bool m_noLoop;
    bool m_showProtein;
    int m_showOccupancy;

    bool m_residueHighlight;
    
    PProteinCCDSolver *m_ccd;
    Vector3 m_centeringTranslate;
    //backbone
    int m_index;
    //sidechain
    int m_side_index;
    //residue
    int m_residue_index;
    int m_autoTweakTimes;  
    int m_animateIndex;    
  
};

#endif
