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

#include "GLUTProgram.h"
#include "GL.h"
#include <GL/glut.h>
#include <assert.h>
#include <stdio.h>

GLUTProgramBase::GLUTProgramBase(int w,int h)
:main_window(0),width(w),height(h),fullscreen_mode(false)
{}

GLUTProgramBase* GLUTProgramBase::current_program=NULL;
void GLUTProgramBase::DisplayFunc() { current_program->Handle_Display(); }
void GLUTProgramBase::ReshapeFunc(int w,int h) { current_program->Handle_Reshape(w,h); }
void GLUTProgramBase::KeyboardFunc(unsigned char c,int x,int h) { current_program->Handle_Keypress(c,x,h); }
void GLUTProgramBase::SpecialFunc(int key,int x,int h) { current_program->Handle_Special(key,x,h); }
void GLUTProgramBase::MouseFunc(int button,int state,int x,int y) { current_program->Handle_Click(button,state,x,y); }
void GLUTProgramBase::MotionFunc(int x,int y) { current_program->Handle_Drag(x,y); }
void GLUTProgramBase::PassiveMotionFunc(int x,int y) { current_program->Handle_Motion(x,y); }
void GLUTProgramBase::IdleFunc() { current_program->Handle_Idle(); }

int GLUTProgramBase::Run(const char *window_title,unsigned int mode)
{
	current_program=this;
	int argc=1;char *(argv[1]);argv[0]="Program";
	glutInit(&argc,(char**)argv);
	if(mode == 0) mode=GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH;
	glutInitDisplayMode(mode);
	glutInitWindowSize(width,height);
	main_window = glutCreateWindow(window_title);
	glutDisplayFunc(DisplayFunc);
	glutReshapeFunc(ReshapeFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutSpecialFunc(SpecialFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
	glutPassiveMotionFunc(PassiveMotionFunc);
	glutIdleFunc(IdleFunc);

	if(!Initialize()) return -1;

	glutMainLoop();
	return 0;
}

bool GLUTProgramBase::Initialize()
{
	glEnable (GL_DEPTH_TEST);
	glEnable (GL_CULL_FACE);

	return true;
}

void GLUTProgramBase::Refresh()
{
	glutPostRedisplay();
}

void GLUTProgramBase::SetFullscreen(bool fullscreen_on)
{
	if(fullscreen_mode != fullscreen_on) {
		fullscreen_mode=fullscreen_on;
		if(fullscreen_mode) { glutFullScreen(); saved_width=width; saved_height=height; }
		else glutReshapeWindow(saved_width,saved_height);
	}
}



