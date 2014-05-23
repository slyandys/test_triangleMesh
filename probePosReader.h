#ifndef _PROBE_POS_READER_H_
#define _PROBE_POS_READER_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "g_Vector.h"

using std::endl;
using std::cerr;

class ProbePosReader {
 protected:
  std::vector<g_Vector> d_pos;
  std::vector<g_Vector> d_orient;
  std::vector<double> d_depth;
  std::vector<g_Vector> d_force;
  double d_radius;

 public:
  inline ProbePosReader( std::string fileName );

  inline void getPos(std::vector<g_Vector>& pos); 
  inline void getOrient(std::vector<g_Vector>& orient);
  inline void getDepth(std::vector<double>& depth);
  inline void getForce(std::vector<g_Vector>& force);
  inline double getRadius();
};


ProbePosReader::ProbePosReader( std::string fileName ) : d_radius(-1.0) {
  std::ifstream in(fileName.c_str());
  if ( !in ) {
    cerr << "Could not open file for input: " << fileName << endl;
    throw fileName;
  }
  std::string line;
  while ( in ) {
    std::getline(in,line);
    // Skip comment lines
    if ( line[0] == '#') continue;
    // Skip empty lines
    if ( line[0] == 'R'|| line[0] == 'r') {
        std::istringstream streamRadius(line.substr(1));
        streamRadius >> d_radius;
    	continue;
    }
    if ( line.size() == 0 ) continue;
    std::istringstream streamLine(line);
    double x,y,z;
    g_Vector vec;
    // position
    streamLine >> x >> y >> z;
    vec.Set( x, y, z );
    d_pos.push_back(vec);
    // orientation
    streamLine >> x >> y >> z;
    vec.Set( x, y, z );
    d_orient.push_back(vec);
    // depth
    streamLine >> x;
    d_depth.push_back(x);
    // force
    streamLine >> x >> y >> z;
    vec.Set( x, y, z );
    d_force.push_back(vec);
    if ( streamLine.fail() ) {
      cerr << "Could not parse input: " << line << endl;
      in.close();
      throw fileName;
    }
  }
}


void ProbePosReader::getPos(std::vector<g_Vector>& pos) {
  pos = d_pos;
  return;
}

void ProbePosReader::getOrient(std::vector<g_Vector>& orient) {
  orient = d_orient;
  return;
}

void ProbePosReader::getDepth(std::vector<double>& depth) {
  depth = d_depth;
  return;
}

void ProbePosReader::getForce(std::vector<g_Vector>& force) {
  force = d_force;
  return;
}

double ProbePosReader::getRadius() {
  return d_radius;
}

  
#endif





