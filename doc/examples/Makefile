##
## Sample makefile for a project that links against LoopTK.
##

# Edit these to point where GLUT(Mesa), GLUI, and GSL are located
ROOT = /home/Yajia
SRCDIR = $(ROOT)/Library
GLUT_INCLUDE = $(SRCDIR)/Mesa-8.0.2/include
GLUI_INCLUDE = $(SRCDIR)/glui-2.35/src/include/GL/
GSL_INCLUDE = /home/Yajia/Library/gsl-1.15/gsl

# Edit these to point where the libraries are located
LIBROOT = $(ROOT)/Library
#GLUILIBDIR = $(LIBROOT)/glui-2.35
GLUILIBDIR = /lib
#GLUTLIBDIR = $(LIBROOT)/glut-3.7.1
GLUTLIBDIR = /lib
#XLIBDIR = /usr/X11R6/lib /usr/X11R6/lib/modules/extensions /usr/X11R6/lib64
XLIBDIR = /lib
SYSLIBDIR = /lib
GSLLIBDIR = /usr/local/lib
#LIBDIRS = $(TK) $(GLUILIBDIR) $(GLUTLIBDIR) $(XLIBDIR) $(SYSLIBDIR) $(GSLLIBDIR)
LIBDIRS = $(TK) $(GLUILIBDIR) $(GLUTLIBDIR) $(XLIBDIR) $(SYSLIBDIR) $(GSLLIBDIR)

GLUTLIBS = glui glut GLU GL
GSLLIBS = gsl gslcblas
SYSLIBS = Xm Xi Xext Xmu X11 SM ICE pthread
LIBS = $(GLUTLIBS) $(GSLLIBS) $(SYSLIBS)

# Where additional libraries (special thanks to Kris Hauser) is located
UTILS_INCLUDE = ../../src/utils
LOOPTK = ../../src/core
MY_INCLUDE = $(UTILS_INCLUDE) $(LOOPTK) $(GLUT_INCLUDE) $(GLUI_INCLUDE) $(GSL_INCLUDE)

CXX = g++
LDFLAGS = ../../lib/looptk.a $(addprefix -L, $(LIBDIRS)) $(addprefix -l, $(LIBS)) $(shell xml2-config --libs)
CPPFLAGS = -O3 $(addprefix -I,$(MY_INCLUDE)) $(shell xml2-config --cflags) -Wno-deprecated
OBJS = $(patsubst %.cc,%.o,$(wildcard *.cc))

# What the program will be named
EXECUTABLE = example

default : $(OBJS)
	$(CXX) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)

clean : 
	/bin/rm -f *.o a.out $(EXECUTABLE) $(EXECUTABLE).purify core Makefile.dependencies

immaculate: clean
	rm -fr *~
