// ==========================================================================
// $Id: g_Vector.h 3737 2013-03-08 13:56:27Z jlang $
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
// $Rev: 3737 $
// $LastChangedBy: jlang $
// $LastChangedDate: 2013-03-08 08:56:27 -0500 (Fri, 08 Mar 2013) $
// ==========================================================================
#ifndef WRAPPER_G_VECTOR_H
#define WRAPPER_G_VECTOR_H

// Should use a namespace

#include <OpenMesh/Core/Geometry/VectorT.hh>
using OpenMesh::Vec3d;

// #include "g_Node.h"

using std::cerr;
using std::endl;

class g_Vector : public Vec3d {
 public:
  inline double x() {
    return (*this)[0];
  }
  inline double y() {
    return (*this)[1];
  }
  inline double z() {
    return (*this)[2];
  }
  // cctor
  g_Vector( const g_Vector& _oVec ) : Vec3d( _oVec ) {}
  g_Vector( const Vec3d& _oVec ) : Vec3d( _oVec ) {}
  g_Vector() {}
  // try cast to Vec3d

  inline g_Vector operator=(const g_Vector& _oVec ) {
    if (&_oVec != this ) {
      Vec3d::operator=(_oVec);
    }
    return *this;
  }

  // g_Vector( float&, float&, float& );
  // g_Vector( double&, double&, double& );
  inline g_Vector( double _x, double _y, double _z) : Vec3d( _x, _y, _z ) {}  
  // void Set( const float&, const float&, const float& );
  // void Set( const double&, const double&, const double& );
  inline void Set( double _x, double _y, double _z) {
    data()[0] = _x;
    data()[1] = _y;
    data()[2] = _z;
  }
  inline void Normalize() {
    normalize();
  }
  // scaling the vector
  inline g_Vector operator*(double _s) {
    Vec3d::operator*=( _s );
    return *this;
  }
  inline g_Vector operator/(double _s) {
    Vec3d::operator/=( _s );
    return *this;
  }
  inline g_Vector operator-(const g_Vector& _oVec) const {
    g_Vector res( *this );
    res -= _oVec;
    return res;
  }
  inline g_Vector operator+(const g_Vector& _oVec) const {
    g_Vector res( *this );
    res += _oVec;
    return res;
  }
  inline double SquaredLength() const {
    return sqrnorm();
  }
  // double DistanceTo( g_Node );
  // double DistanceTo(const g_Vector&);
  inline double DistanceTo( const Vec3d& _oVec ) const {
    // cerr << static_cast<Vec3d>(*this) << " to " << static_cast<Vec3d>(_oVec) << endl;
    Vec3d res(*this);
    res -= _oVec;
    return res.norm();
  }   
  inline double Dot(const Vec3d& _oVec) const {
    return operator|( _oVec );
  }
  inline double Dot(const g_Vector* _oVec) const {
    return operator|( *_oVec );
  }
  inline g_Vector Cross(const g_Vector& _oVec ) const {
    return g_Vector( operator%( _oVec ));
  }
  inline double AngleBetween(const g_Vector& _oVec ) const {
    return acos((operator|(_oVec))/norm()/_oVec.norm());
  }
  inline double Length() const {
    return norm();
  }
};
  
g_Vector operator*(double _s,const g_Vector& _oVec);

#endif
