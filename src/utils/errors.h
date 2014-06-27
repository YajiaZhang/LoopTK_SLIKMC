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

#ifndef ERRORS_H
#define ERRORS_H

#include <stdio.h>
#include <stdarg.h>

/* Aborts and asserts don't give a proper stack trace on cygwin. */
#ifdef CYGWIN

#ifdef NDEBUG           /* required by ANSI standard */
  #define Assert(p)  	((void)0)
#else

  #ifdef __STDC__
  #define Assert(e)       ((e) ? (void)0 : __Assert(__FILE__, __LINE__, #e))
  #else   /* PCC */
  #define Assert(e)       ((e) ? (void)0 : __Assert(__FILE__, __LINE__, "e"))
  #endif

#endif // NDEBUG

#define Abort() { fprintf(stderr,"Aborting with segfault...\n"); *((int*)0)=1; }

inline void __Assert(const char * file, int line, const char * e)
{
  fprintf(stderr,"Assertion \"%s\" failed: file \"%s\", line %d\n",e,file,line);
  Abort();
}


#else

  #include "assert.h"
  #include "stdlib.h"
  #define Assert assert
  #define Abort abort
  #define __Assert __assert

#endif //CYGWIN

inline void RaiseErrorFmt(const char* func, const char* file, int line, const char* fmt, ...)
{
  fprintf(stderr,"Error in %s (%s:%d): ", func,file,line); 
  va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
  fprintf(stderr,"\n");
  Abort();
}

inline void RaiseErrorFmt(const char* fmt,...)
{
  fprintf(stderr,"Error (unknown function): ");
  va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
  fprintf(stderr,"\n");
  Abort();
}


inline void RaiseError(const char* func, const char* file, int line, const char* text)
{
  fprintf(stderr,"Error in %s (%s:%d): %s\n", func,file,line,text); 
  Abort();
}



//the following is bending over backwards to support MS's lack of 
//variable argument macro support

#ifdef HAVE_PRETTY_FUNCTION
#define WHERE_AM_I __PRETTY_FUNCTION__, __FILE__, __LINE__
#else
#define WHERE_AM_I __FUNCTION__, __FILE__, __LINE__
#endif

//Error1 is guaranteed to print line numbers with a single-argument error
#define FatalError1(text) RaiseError(WHERE_AM_I,text)

#if HAVE_VARARGS_MACROS
#define FatalError(fmt,...) RaiseErrorFmt(WHERE_AM_I,fmt,__VA_ARGS__)
#else
//if no variable arguments, can't get any line info 
#define FatalError RaiseErrorFmt
#endif

#define PrintLocation(file)  fprintf(file,"%s (%s:%d): ", WHERE_AM_I)
#define AssertNotReached() RaiseError(WHERE_AM_I,"Code should not be reached")

#endif

