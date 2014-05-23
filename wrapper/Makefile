#!/bin/csh

#
# $Id: Makefile,v 1.1 2011/03/12 09:17:36 jlang Exp $
#
# $Log: Makefile,v $
#

CXX = g++
CC = gcc

ifndef OSTYPE 
OSTYPE=linux-gnu
endif

# CXXFLAGS=-fpermissive

COMMON=geodesic
INCDIR = -I$(COMMON)


# OPENMESH_ROOT=/usr/local/include/OpenMesh
OPENMESH_ROOT=/home/langj/projects/defoCap/tracking/OpenMesh

ifeq ($(OSTYPE),linux-gnu)
INCDIR += -I$(OPENMESH_ROOT)/include
LDLIBDIR += -L$(OPENMESH_ROOT)/lib/OpenMesh
LDLIBS += -lm -lOpenMeshCore -lOpenMeshTools
endif



.PHONY: all clean


SRC =  tester.cpp g_Part.cpp g_Element.cpp g_Node.cpp g_PEdge.cpp g_Plane.cpp g_Vector.cpp TriangleMesh.cpp


# Default build : 
all: dependencies tester 


dependencies: $(SRC:%.cpp=%.d) $(SRC:%.c=%.d)


tester: $(SRC:%.cpp=%.o) $(SRC_C:%.c=%.o)
	$(CXX) $(INCDIR) $(LDLIBDIR) $(DEFINES) -o $@ $^ $(LDLIBS)

# sed line based on John Graham-Cumming, Dependency Management,
# Dr.Dobb's Portal, Apr 01, 2006.
%.d:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC:.cpp=.d)

%.d:%.c
	$(CC) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC_C:.c=.d)


%.o:%.c
	$(CC) $(INCDIR) -c $*.c -o $@ 

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -c $*.cpp -o $@ 

clean:
	rm -f *~ *.o *.d




