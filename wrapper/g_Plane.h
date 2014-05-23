// ==========================================================================
// $Id: g_Plane.h 3389 2012-10-19 08:23:25Z jlang $
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
// $Rev: 3389 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-10-19 04:23:25 -0400 (Fri, 19 Oct 2012) $
// ==========================================================================
#ifndef WRAPPER_G_PLANE_H
#define WRAPPER_G_PLANE_H

#include "g_Vector.h"
#include "g_Element.h"
#include "g_Node.h"

// Should use a namespace

class g_Plane {
 protected:
  g_Vector d_normal;
  double d_dist;
 public:
  g_Plane(const g_Vector&,const g_Vector&,const g_Vector&);
  g_Plane(const g_Vector&,const g_Vector&);
  g_Vector* normal();
  g_Vector project( const g_Vector& ) const;
  double distanceToPoint(const g_Node&) const;
  double distanceToPoint(const g_Vector&) const;
  bool line_plane_intersection(const g_Vector&, const g_Vector&, 
			       g_Vector&, double& tVal ) const;
};
  

#endif
