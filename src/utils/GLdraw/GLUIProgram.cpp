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

#include "GLUIProgram.h"
#include "GL.h"
#include "GLError.h"
#include "glui.h"
#include <assert.h>
#include <stdio.h>

#define DEBUG(x) { x; }
/*
#define DEBUG(x) { \
  fprintf(stderr,"Debug begin %s\n",#x); \
  x; \
  fprintf(stderr,"Debug end %s\n",#x); \
}
*/

GLUIProgramBase::GLUIProgramBase(int w,int h)
:main_window(0),width(w),height(h),fullscreen_mode(false)
{}

GLUIProgramBase* GLUIProgramBase::current_program=NULL;
void GLUIProgramBase::DisplayFunc() {
  DEBUG(current_program->Handle_Display());
  CheckGLErrors("DisplayFunc");
}
void GLUIProgramBase::ReshapeFunc(int w,int h) { DEBUG(current_program->Handle_Reshape(w,h)); }
void GLUIProgramBase::KeyboardFunc(unsigned char c,int x,int h) { DEBUG(current_program->Handle_Keypress(c,x,h)); }
void GLUIProgramBase::SpecialFunc(int key,int x,int h) { DEBUG(current_program->Handle_Special(key,x,h)); }
void GLUIProgramBase::MouseFunc(int button,int state,int x,int y) { DEBUG(current_program->Handle_Click(button,state,x,y)); }
void GLUIProgramBase::MotionFunc(int x,int y) { DEBUG(current_program->Handle_Drag(x,y)); }
void GLUIProgramBase::PassiveMotionFunc(int x,int y) { DEBUG(current_program->Handle_Motion(x,y)); }
void GLUIProgramBase::IdleFunc() 
{ 
  glutSetWindow(current_program->main_window);
  DEBUG(current_program->Handle_Idle());
}
void GLUIProgramBase::ControlFunc (int id) { current_program->Handle_Control(id); }

int GLUIProgramBase::Run(const char *window_title,unsigned int mode)
{
	current_program=this;
	int argc=1;char *(argv[1]);argv[0]="Program";
	glutInit(&argc,(char**)argv);
	if(mode == 0) mode=GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH;
	glutInitDisplayMode(mode);
	glutInitWindowSize(width,height);
	main_window = glutCreateWindow(window_title);
	GLUI_Master.set_glutDisplayFunc(DisplayFunc);
	GLUI_Master.set_glutReshapeFunc(ReshapeFunc);
	GLUI_Master.set_glutKeyboardFunc(KeyboardFunc);
	GLUI_Master.set_glutSpecialFunc(SpecialFunc);
	GLUI_Master.set_glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
	glutPassiveMotionFunc(PassiveMotionFunc);
	GLUI_Master.set_glutIdleFunc(IdleFunc);

	if(!Initialize()) return -1;

	glutMainLoop();
	return 0;
}

bool GLUIProgramBase::Initialize()
{
	glEnable (GL_DEPTH_TEST);
	glEnable (GL_CULL_FACE);
	return true;
}

void GLUIProgramBase::Refresh()
{
	glutPostRedisplay();
}

void GLUIProgramBase::SetFullscreen(bool fullscreen_on)
{
	if(fullscreen_mode != fullscreen_on) {
		fullscreen_mode=fullscreen_on;
		if(fullscreen_mode) { glutFullScreen(); saved_width=width; saved_height=height; }
		else glutReshapeWindow(saved_width,saved_height);
	}
}



