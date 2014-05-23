#!/bin/csh

#
# $Id: Makefile 3862 2013-04-25 15:25:51Z jlang $
#
# $Log: Makefile,v $
#

CXX = g++
CC = gcc

ifndef OSTYPE 
OSTYPE=linux-gnu
endif

CXXFLAGS= -DUSE_WRAPPER -g
# -O3
# -fpermissive

COMMON= wrapper
INCDIR = -I$(COMMON)

GEODESIC=wrapper/geodesic
INCDIR += -I$(GEODESIC)

# ITK_ROOT=/usr/local/include/InsightToolkit
ITK_ROOT=/usr/include/InsightToolkit
# OPENMESH_ROOT=/usr/local/include/OpenMesh
OPENMESH_ROOT=/home/langj/projects/defoCap/tracking/OpenMesh
ANN_ROOT=/home/langj/projects/defoCap/tracking/ann_1.1.2


ifeq ($(OSTYPE),linux-gnu)
DEFINES = -DNO_F2C
DEFINES += -DUSE_MOELLER
# LDLIBDIR += -LCLAPACK-3.2.1
# INCDIR += -ICLAPACK-3.2.1/INCLUDE 
# INCDIR += -I$(ITK_ROOT)/Algorithms
# INCDIR += -I$(ITK_ROOT)/BasicFilters
# INCDIR += -I$(ITK_ROOT)/Common
# INCDIR += -I$(ITK_ROOT)/Numerics
# INCDIR += -I$(ITK_ROOT)/IO
# INCDIR += -I$(ITK_ROOT)/Numerics/FEM
# INCDIR += -I$(ITK_ROOT)/Numerics/NeuralNetworks
# INCDIR += -I$(ITK_ROOT)/Numerics/Statistics
# INCDIR += -I$(ITK_ROOT)/SpatialObject
# INCDIR += -I$(ITK_ROOT)/Utilities/MetaIO
# INCDIR += -I$(ITK_ROOT)/Utilities/NrrdIO
# INCDIR += -I$(ITK_ROOT)/Utilities/DICOMParser
# INCDIR += -I$(ITK_ROOT)/Utilities/expat
# INCDIR += -I$(ITK_ROOT)/Utilities/itkExtHdrs
# INCDIR += -I$(ITK_ROOT)/Utilities
# INCDIR += -I$(ITK_ROOT)/Utilities/vxl/v3p/netlib
INCDIR += -I$(ITK_ROOT)/Utilities/vxl/vcl
INCDIR += -I$(ITK_ROOT)/Utilities/vxl/core
INCDIR += -Itetgen1.4.3 
INCDIR += -I$(OPENMESH_ROOT)/include
INCDIR += -I$(ANN_ROOT)/include
LDLIBDIR += -Ltetgen1.4.3
LDLIBDIR += -Lann_1.1.2/lib
LDLIBDIR += -L$(OPENMESH_ROOT)/lib/OpenMesh
LDLIBS += -llapack -lblas -lm -lANN -litkvcl -litkvnl -litkvnl_algo -ltet -lOpenMeshCore -lOpenMeshTools
endif

.PHONY: all clean

SRC_CORE =  asa047.cpp util_wrapper.cpp FineFittingTrackingLocalFrame.cpp LinearFemDeformation.cpp main.cpp MdLibrary.cpp MeshSimplification.cpp CorrespondSurfaceAndVolumetricTemplate.cpp tritri.cpp ElasticParameterOptimization.cpp RBFDeform.cpp

SRC_WRAPPER = g_Part.cpp g_Element.cpp g_Node.cpp g_PEdge.cpp g_Vector.cpp g_Plane.cpp TriangleMesh.cpp

SRC = $(SRC_CORE) $(SRC_WRAPPER:%=$(COMMON)/%)

# Default build : 
all: dependencies tracker 

dependencies: $(SRC:%.cpp=%.d) $(SRC:%.c=%.d)
	

# sed line based on John Graham-Cumming, Dependency Management,
# Dr.Dobb's Portal, Apr 01, 2006.
%.d:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC:.cpp=.d)

%.d:%.c
	$(CC) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC_C:.c=.d)


tracker: $(SRC:%.cpp=%.o) $(SRC_C:%.c=%.o)
	$(CXX) $(INCDIR) $(LDLIBDIR) $(DEFINES) -o $@ $^ $(LDLIBS)

%.o:%.c
	$(CC) $(INCDIR) -c $*.c -o $@ 

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -c $*.cpp -o $@ 

clean:
	rm -f *~ *.o *.d




