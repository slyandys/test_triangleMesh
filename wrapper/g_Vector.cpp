// ==========================================================================
// $Id: g_Vector.cpp 3389 2012-10-19 08:23:25Z jlang $
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
#ifndef WRAPPER_G_VECTOR_CPP
#define WRAPPER_G_VECTOR_CPP

// Should use a namespace

#include "g_Vector.h"

g_Vector operator*(double _s,const g_Vector& _oVec) {
  g_Vector res(_oVec);
  res *= _s;
  return res;
}

#endif
