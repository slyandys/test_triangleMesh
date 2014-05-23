// Geometric Tools, LLC
// Copyright (c) 1998-2012
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef DIST_POINT3_TRIANGLE3_H
#define DIST_POINT3_TRIANGLE3_H


#include "g_Element.h"
#include "g_Vector.h"
#include "g_Vector.h"

namespace Wm5
{
//----------------------------------------------------------------------------
double computeClosestPoint( g_Node& point,  g_Element& triangle, g_Vector& closestPoint )
{
  g_Vector diff = triangle.nodes()[0]->coordinate() - point.coordinate();
  g_Vector edge0 = triangle.nodes()[1]->coordinate() - triangle.nodes()[0]->coordinate();
  g_Vector edge1 = triangle.nodes()[2]->coordinate() - triangle.nodes()[0]->coordinate();
  double a00 = edge0.SquaredLength();
  double a01 = edge0.Dot(edge1);
  double a11 = edge1.SquaredLength();
  double b0 = diff.Dot(edge0);
  double b1 = diff.Dot(edge1);
  double c = diff.SquaredLength();
  double det = fabs(a00*a11 - a01*a01);
  double s = a01*b1 - a11*b0;
  double t = a01*b0 - a00*b1;
  double sqrDistance;
  
  if (s + t <= det)
    {
        if (s < 0.0)
        {
            if (t < 0.0)  // region 4
            {
                if (b0 < 0.0)
                {
                    t = 0.0;
                    if (-b0 >= a00)
                    {
                        s = 1.0;
                        sqrDistance = a00 + 2.0*b0 + c;
                    }
                    else
                    {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
                else
                {
                    s = 0.0;
                    if (b1 >= 0.0)
                    {
                        t = 0.0;
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11)
                    {
                        t = 1.0;
                        sqrDistance = a11 + 2.0*b1 + c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else  // region 3
            {
                s = 0.0;
                if (b1 >= 0.0)
                {
                    t = 0.0;
                    sqrDistance = c;
                }
                else if (-b1 >= a11)
                {
                    t = 1.0;
                    sqrDistance = a11 + 2.0*b1 + c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < 0.0)  // region 5
        {
            t = 0.0;
            if (b0 >= 0.0)
            {
                s = 0.0;
                sqrDistance = c;
            }
            else if (-b0 >= a00)
            {
                s = 1.0;
                sqrDistance = a00 + 2.0*b0 + c;
            }
            else
            {
                s = -b0/a00;
                sqrDistance = b0*s + c;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double invDet = (1.0)/det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                t*(a01*s + a11*t + 2.0*b1) + c;
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if (s < 0.0)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0*a01 + a11;
                if (numer >= denom)
                {
                    s = 1.0;
                    t = 0.0;
                    sqrDistance = a00 + 2.0*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = 1.0 - s;
                    sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                        t*(a01*s + a11*t + 2.0*b1) + c;
                }
            }
            else
            {
                s = 0.0;
                if (tmp1 <= 0.0)
                {
                    t = 1.0;
                    sqrDistance = a11 + 2.0*b1 + c;
                }
                else if (b1 >= 0.0)
                {
                    t = 0.0;
                    sqrDistance = c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < 0.0)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0*a01 + a11;
                if (numer >= denom)
                {
                    t = 1.0;
                    s = 0.0;
                    sqrDistance = a11 + 2.0*b1 + c;
                }
                else
                {
                    t = numer/denom;
                    s = 1.0 - t;
                    sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                        t*(a01*s + a11*t + 2.0*b1) + c;
                }
            }
            else
            {
                t = 0.0;
                if (tmp1 <= 0.0)
                {
                    s = 1.0;
                    sqrDistance = a00 + 2.0*b0 + c;
                }
                else if (b0 >= 0.0)
                {
                    s = 0.0;
                    sqrDistance = c;
                }
                else
                {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= 0.0)
            {
                s = 0.0;
                t = 1.0;
                sqrDistance = a11 + 2.0*b1 + c;
            }
            else
            {
                denom = a00 - 2.0*a01 + a11;
                if (numer >= denom)
                {
                    s = 1.0;
                    t = 0.0;
                    sqrDistance = a00 + 2.0*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = 1.0 - s;
                    sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
                        t*(a01*s + a11*t + 2.0*b1) + c;
                }
            }
        }
    }

    // Account for numerical round-off error.
    if (sqrDistance < 0.0)
    {
        sqrDistance = 0.0;
    }
    closestPoint =  triangle.nodes()[0]->coordinate() + s*edge0 + t*edge1;
    return sqrDistance;
}

} // end namespace

#endif
