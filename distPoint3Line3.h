// Geometric Tools, LLC
// Copyright (c) 1998-2012
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef DIST_POINT3_LINE3_H
#define DIST_POINT3_LINE3_H

#include "g_Vector.h"

namespace Wm5
{
//----------------------------------------------------------------------------
inline double distPoint3Line3Squared( const g_Vector& point, 
				      const g_Vector& lineOrigin, 
				      const g_Vector& lineDirection, 
				      g_Vector& closestPoint )
{
  g_Vector diff = point - lineOrigin;
  double lineParameter = lineDirection.Dot(diff);
  closestPoint = lineOrigin + lineParameter*lineDirection;
  diff = closestPoint - point;
  return diff.SquaredLength();
}



inline double distPoint3Segment3Squared( const g_Vector& point, 
					 const g_Vector& segOrigin, 
					 const g_Vector& segDirection,
					 double segLength,
					 g_Vector& closestPoint )
{
    g_Vector diff = point - segOrigin;
    double segmentParameter = segDirection.Dot(diff);

    if ( segmentParameter < segLength/2.0 ) {
        if (segmentParameter > 0.0 ) {
	  closestPoint = segOrigin + segmentParameter*segDirection;
        } else {
	  closestPoint = segOrigin;
        } 
    } else {
      if ( segmentParameter > segLength ) {
        closestPoint = segOrigin + segLength*segDirection;
      } else {
	closestPoint = segOrigin + segmentParameter*segDirection;
      }
    }

    diff = closestPoint - point;
    return diff.SquaredLength();
}



} // end of namespace
#endif
