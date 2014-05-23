// ==========================================================================
// $Id: g_Plane.cpp 3497 2012-11-30 17:21:59Z jlang $
// Wrapper code to interface g_plane
// ==========================================================================
// (C)opyright:
//
//   Jochen Lang
//   SITE, University of Ottawa
//   800 King Edward Ave.
//   Ottawa, On., K1N 6N5
//   Canada.
//   http://www.site.uottawa.ca
//
// Creator: Jochen Lang
// Email:   jlang@site.uottawa.ca
// ==========================================================================
// $Rev: 3497 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-11-30 12:21:59 -0500 (Fri, 30 Nov 2012) $
// ==========================================================================
#ifndef WRAPPER_G_PLANE_CPP
#define WRAPPER_G_PLANE_CPP

#include "cmath"

#include "g_Plane.h"

// Constructor from three points
g_Plane::g_Plane(const g_Vector& _p0, const g_Vector& _p1, const g_Vector& _p2) {
  // n = (_p1 - _p0) x (_p2 - _p0)
  g_Vector dir0(_p1), dir1(_p2);
  dir0 -= _p0;
  dir1 -= _p0;
  d_normal = dir0.Cross(dir1);
  d_normal.Normalize();
  d_dist = - d_normal.Dot(_p0);
}

g_Plane::g_Plane(const g_Vector& _normal, const g_Vector& _p ) {
  d_normal = _normal;
  d_normal.Normalize();
  d_dist = - d_normal.Dot(_p);
}

// Ouch by ptr
g_Vector* g_Plane::normal() {
  return &d_normal;
}

g_Vector g_Plane::project( const g_Vector& _p ) const {
  g_Vector res = -distanceToPoint(_p) * d_normal;
  res += _p;
  return res;
}

double g_Plane::distanceToPoint(const g_Node& _n ) const {
  return distanceToPoint( _n.coordinate() );
}

double g_Plane::distanceToPoint(const g_Vector& _p ) const {
  return  (d_normal.Dot( _p ) + d_dist); 
}

bool g_Plane::line_plane_intersection(const g_Vector& _start, const g_Vector& _end, 
				      g_Vector& _is, double& tVal ) const {
  // n . ((e-s)*t + s) + d = 0
  // t = (n . s + d)/(n . (s-e))
 
  double num = (d_normal.Dot(_start) + d_dist);
  double div = d_normal.Dot(_start-_end);
  if ( std::abs(div)>1e-18 ) {
    // Avoid problem when the start or end point are in the plane
    if (std::abs(num) < 1e-9 * std::abs(div)) { 
      tVal=0.0;
    } else {
      if (num - div < 1e-18 ) { 
	tVal=1.0;
      } else { 
	tVal = num/div;
      }
    }
    _is = (1-tVal) * _start + tVal * _end; 
  } else {
    // Line is parallel -- return start point
    _is = _start;
    tVal = 0.0;
    return false;
  }
  if ( tVal >= 0.0 && tVal <= 1.0 )  return true;
  else return false;
}
  

#endif
