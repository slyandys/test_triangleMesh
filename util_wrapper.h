// ==========================================================================
// $Id: util_wrapper.h 3505 2012-12-05 10:38:15Z jlang $
// Utility to keep wrapper and tracking code separate
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
// $Rev: 3505 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2012-12-05 05:38:15 -0500 (Wed, 05 Dec 2012) $
// ==========================================================================
#ifndef UTIL_WRAPPER
#define UTIL_WRAPPER
 
#include <string>
#include <fstream>

#include "TriangleMesh.h"

bool exportMeshWrapper( const std::string, TriangleMesh* );
bool importMeshWrapper(  const std::string, TriangleMesh* );

#endif
