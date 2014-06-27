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

#ifndef GL_GLUI_PROGRAM_H
#define GL_GLUI_PROGRAM_H

class GLUIProgramBase
{
public:
  GLUIProgramBase(int width=800,int height=600);
  ///if displayMode is non-zero, initializes glut with that display mode
  int Run(const char *window_title="OpenGL Viewer",unsigned int displayMode=0);

  ///overrideable
  virtual bool Initialize();
  virtual void Handle_Display() {}
  virtual void Handle_Reshape(int w,int h) { width=w; height=h; }
  virtual void Handle_Keypress(unsigned char key,int x,int y){}
  virtual void Handle_Special(int key,int x,int y) {}
  virtual void Handle_Click(int button,int state,int x,int y){}
  virtual void Handle_Drag(int x,int y){}
  virtual void Handle_Motion(int x,int y){}
  virtual void Handle_Idle(){}

  ///override this to handle GLUI control callbacks
  virtual void Handle_Control(int id) {}

  ///Refreshes the screen (equivalent to glutPostRedisplay())
  void Refresh();
  ///Turns on fullscreen mode
  void SetFullscreen(bool fullscreen_on);

  int main_window;
  int width,height; // window size
  bool fullscreen_mode;
  int saved_width,saved_height;

  ///pass this as a callback for GLUI controls
  static void ControlFunc (int);

private:
  static GLUIProgramBase* current_program;
  static void DisplayFunc();
  static void ReshapeFunc(int w,int h);
  static void KeyboardFunc(unsigned char key,int x,int y);
  static void SpecialFunc(int key,int x,int y);
  static void MouseFunc(int button,int state,int x,int y);
  static void MotionFunc(int x,int y);
  static void PassiveMotionFunc(int x,int y);
  static void IdleFunc();
};

#endif
