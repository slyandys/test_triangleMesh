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

COMMON=wrapper
INCDIR = -I$(COMMON)


ifeq ($(OSTYPE),linux-gnu)
# LDLIBDIR += -Lann_1.1.2/lib
LDLIBS += -lm 
endif




.PHONY: all clean


SRC = example1.cpp  


# Default build : 
all: dependencies example


dependencies: $(SRC:%.cpp=%.d) $(SRC:%.c=%.d)


# sed line based on John Graham-Cumming, Dependency Management,
# Dr.Dobb's Portal, Apr 01, 2006.
%.d:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC:.cpp=.d)

%.d:%.c
	$(CC) $(INCDIR) $(DEFINES) -MG -MM $< | sed 's,\($*\.o\)[ :]*\(.*\), $@ : $$\(wildcard \2\)\n\1 : \2,g' > $@
-include $(SRC_C:.c=.d)


example: $(SRC:%.cpp=%.o) $(SRC_C:%.c=%.o)
	$(CXX) $(INCDIR) $(LDLIBDIR) $(DEFINES) -o $@ $^ $(LDLIBS)

%.o:%.c
	$(CC) $(INCDIR) -c $*.c -o $@ 

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) $(DEFINES) -c $*.cpp -o $@ 

clean:
	rm -f *~ *.o *.d




