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

#ifndef GL_GLUT_NAVIGATION_PROGRAM_H
#define GL_GLUT_NAVIGATION_PROGRAM_H


#include <math3d/geometry3d.h>
#include <camera/viewport.h>
#include <Timer.h>
#include "GLUTProgram.h"

class GLUTNavigationProgram : public GLUTProgramBase
{
public:
  GLUTNavigationProgram();
  virtual bool Initialize();
  virtual void Handle_Display();
  virtual void Handle_Reshape(int w,int h);
  virtual void Handle_Click(int button,int state,int x,int y);
  virtual void Handle_Drag(int x,int y);
  virtual void Handle_Keypress(unsigned char key,int x,int y);
  virtual void Handle_Idle();
  
  //overrideable
  virtual void SetWorldLights() {}
  virtual void RenderWorld() {}
  virtual void RenderScreen() {}
  
  virtual void BeginDrag(int x,int y,int button,int modifiers);
  virtual void DoDrag(int dx,int dy,int button,int modifiers);
  virtual void EndDrag(int x,int y,int button,int modifiers);
  virtual void DoFreeDrag(int dx,int dy,int button);
  virtual void DoCtrlDrag(int dx,int dy,int button);
  virtual void DoAltDrag(int dx,int dy,int button);
  virtual void DoShiftDrag(int dx,int dy,int button);

  void DragPan(int dx,int dy);
  void DragRotate(int dx,int dy);
  void DragZoom(int dx,int dy);
  void DragTruck(int dx,int dy);

  void Set2DMode(bool mode=true);
  void DisplayCameraTarget();
  void CenterCameraOn(const AABB3D& bbox);

  void WriteDisplaySettings(std::ostream& out) const;
  void ReadDisplaySettings(std::istream& in);
  
  Viewport viewport;
  CameraController_Orbit camera;
  int oldmousex,oldmousey;
  int clickButton, clickModifiers;
  bool stereo_mode;
  float stereo_offset;
  Timer timer;
  int show_view_target; float t_hide_view_target;
  float frames_per_second; bool show_frames_per_second;
  int frames_rendered;
  bool mode_2d;
};

#endif
