# Edit these to point where GLUT (Mesa) and GLUI are located.
#ROOT = /afs/cs/group/latombe
#SRCDIR = $(ROOT)/src
ROOT = /home/Yajia
SRCDIR = $(ROOT)/Library
GLUT_INCLUDE = $(SRCDIR)/Mesa-8.0.2/include
GLUI_INCLUDE = $(SRCDIR)/glui-2.36/src/include/GL/

# Where additional libraries (special thanks to Kris Hauser) is located
UTILS_INCLUDE = ./src/utils
MY_INCLUDE = $(UTILS_INCLUDE) $(GLUT_INCLUDE) $(GLUI_INCLUDE)

DOCDIR = doc

# Directory in which to create the LoopTK library.
OUTLIBDIR = ./lib
 
# The primary flags for the compiler and linker.
AR = ar
ARFLAGS = -rcs

CPPFLAGS = -O3 $(addprefix -I,$(MY_INCLUDE)) $(shell xml2-config --cflags) -Wno-deprecated
CXX = g++

# Doxygen configuration file
DOXYGEN_CONF = doxygen.conf

LOOPTK_OBJECTS = $(patsubst %.cc,%.o,$(wildcard src/core/*.cc)) 

OBJECTS = $(LOOPTK_OBJECTS)
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard src/utils/*.cpp)) 
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard src/utils/camera/*.cpp))
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard src/utils/math3d/*.cpp))
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard src/utils/math/*.cpp))
OBJECTS += $(patsubst %.cpp,%.o,$(wildcard src/utils/GLdraw/*.cpp))

# Build targets.

lib : $(OBJECTS)
	mkdir -p $(OUTLIBDIR)
	$(AR) $(ARFLAGS) $(OUTLIBDIR)/looptk.a $(OBJECTS)

docs:
	doxygen $(DOXYGEN_CONF)

default: lib
 
# Cleaning targets.

coreclean :
	rm -f src/core/*.o lib/*

clean : coreclean
	rm -f src/utils/camera/*.o src/utils/*.o src/utils/math3d/*.o src/utils/math/*.o src/utils/GLdraw/*.o

realclean : clean
	rm -rf $(DOCDIR)/html $(DOCDIR)/latex

main:
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)
