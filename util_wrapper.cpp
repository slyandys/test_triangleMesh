// ==========================================================================
// $Id: util_wrapper.cpp 3505 2012-12-05 10:38:15Z jlang $
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
#include "util_wrapper.h"

using std::string;
using std::ofstream;
using std::ios_base;

bool exportMeshWrapper( const string _buffer, TriangleMesh* _mesh ) {
  ofstream ofs;
#ifdef USE_WRAPPER
  string fName(_buffer);
  int pos = fName.rfind( ".wrl" );
  if (pos!=string::npos)
    fName.replace (pos,4,".ply");
  ofs.open(fName.c_str(), ios_base::out|ios_base::binary);
#else
  ofs.open(_buffer.c_str());
#endif
  if ( ofs.fail() ) {
    std::cerr << "File open failed: " << _buffer << std::endl;
    return false;
  }
  _mesh->exportVRML(ofs);
  ofs.close();
  return true;
}


bool importMeshWrapper(  const string _buffer, TriangleMesh* _mesh ) {
  ifstream ifs;
#ifdef USE_WRAPPER
  string fName(_buffer);
  int pos = fName.rfind( ".wrl" );
  if (pos!=string::npos)
    fName.replace (pos,4,".ply");
  ifs.open(fName.c_str(), ios_base::in|ios_base::binary);
#else
  ifs.open(_buffer.c_str());
#endif
  if ( ifs.fail() ) {
    std::cerr << "File open failed: " << _buffer << std::endl;
    return false;
  }
  _mesh->init(ifs);
  ifs.close();
  return true;
}
  
