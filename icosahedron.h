#ifndef _ICOSAHEDRON_H_
#define _ICOSAHEDRON_H_

#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>

#include "g_Vector.h"

#define ICOSAHEDRON_EPSILON 1e-5


class Icosahedron {
 protected: 
  int d_q; // frequency of tesselation
  // offset for coordinates
  g_Vector d_offset;
  double d_radius;
  int d_maxLevels, d_maxVert;
  g_Vector d_rot[3];
  bool d_rotate;


 public:
  inline Icosahedron( int q, g_Vector offset=g_Vector(), double radius=1.0 );
  inline void setRotation( g_Vector orientation );
  inline int generateTesselation(std::vector<g_Vector>& vertices, std::vector<g_Vector>& normals, double z_min=-1e9 );
};


Icosahedron::Icosahedron( int q, g_Vector offset, double radius ) :
d_q(q), d_offset(offset), d_radius(radius), d_rotate(false) {
  d_maxLevels = 1 << q;
  // maximum number of points per horizontal circle
  d_maxVert = 5 * d_maxLevels;
  // Correct maxLevels by one
  d_maxLevels++;
  // number of points
  int n = 10 * ( 1 << q ) * ( 1 << q ) + 2;
}

/*
 * Expects the direction, i.e., orientation of the probe
 */
void Icosahedron::setRotation( g_Vector orientation ) {
  // new z-axis
  orientation.Normalize();
  // calculate rotation axis
  g_Vector r(-orientation.y(), orientation.x(), 0.0 );
  if ( r.Length() < ICOSAHEDRON_EPSILON ) {
    if ( orientation.z() <= -1.0 + ICOSAHEDRON_EPSILON ) {
      // flip the z and x-axis
      d_rot[0] = g_Vector( -1.0, 0.0, 0.0 );
      d_rot[1] = g_Vector( 0.0, 1.0, 0.0 );
      d_rot[2] = g_Vector( 0.0, 0.0, -1.0 );
      d_rotate = true;
      return;
    } else  {
      d_rotate = false;
      return;
    }
  }
  r.Normalize();
  g_Vector ry = orientation.Cross(r);
  d_rot[0].Set( r.x(), ry.x(), orientation.x() );
  d_rot[1].Set( r.y(), ry.y(), orientation.y() );
  d_rot[2].Set( r.z(), ry.z(), orientation.z() );
  d_rotate = true;
}


int Icosahedron::generateTesselation(std::vector<g_Vector>& vertices, std::vector<g_Vector>& normals, double depth ) {
  int totalLevels = d_maxLevels + 2 * (d_maxVert/5);
  int noVert = 1; 
  double elevationStep = M_PI / ( totalLevels - 1 );
  double elevation = M_PI/2;
  double orientOffset = 0; 

  for ( int level = 0; level < totalLevels; level++ ) {
    double z = d_radius * sin(elevation);
    if ( d_radius - depth > z ) return  vertices.size();
    double cosE = d_radius * cos(elevation);
    // Generate vertices for the current level
    for ( int i=0; i<noVert; ++i ) {
      // elevation
      double azimuth = 2.0*M_PI *((double)i+orientOffset)/((double)noVert);
      double y = cosE * sin(azimuth);
      double x = cosE * cos(azimuth);
      if ( d_rotate ) {
    	  g_Vector pos(x,y,z);
    	  vertices.push_back(g_Vector(d_offset.x()+pos.Dot( d_rot[0] ),d_offset.y()+pos.Dot( d_rot[1] ),d_offset.z()+pos.Dot( d_rot[2] )));
    	  g_Vector norm(-x,-y, -z);
    	  norm.Normalize();
    	  normals.push_back(g_Vector(norm.Dot( d_rot[0] ), norm.Dot( d_rot[1] ), norm.Dot( d_rot[2] )));
      } else {
    	  vertices.push_back(g_Vector( x+d_offset.x(),y+d_offset.y(),z+d_offset.z()));
    	  normals.push_back(g_Vector(-x, -y, -z ));
      }
    }
    // Calculate elevation angle for next level
    elevation -= elevationStep;
    // Adjust orientation offset
    if ( orientOffset > 0 + ICOSAHEDRON_EPSILON ) {
      orientOffset = 0; 
    } else {
      orientOffset = 0.5; 
    }
    // Calculate the number of vertices on the next level
    if  ( level+1 < ( d_maxVert/5 )) {
      // increasing
      if ( noVert == 1 ) noVert = 0;
      noVert += 5;
    } else {
      if  ( level+1 < d_maxLevels + ( d_maxVert/5 )) {
	// maximum 
	noVert = d_maxVert;
      } else {
	// decreasing
	noVert -= 5;
	if ( noVert == 0 ) { noVert = 1; orientOffset = 0;}
      }
    }
  }
  return vertices.size();
}
  
#endif





