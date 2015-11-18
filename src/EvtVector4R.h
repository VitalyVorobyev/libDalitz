#ifndef EVTVECTOR4R_H
#define EVTVECTOR4R_H

#include <iostream>
#include <cmath>
#include "EvtVector3R.h"

class EvtVector4R {
  friend EvtVector4R rotateEuler(const EvtVector4R& rs,const double& alpha,const double& beta,const double& gamma);
  friend EvtVector4R boostTo(const EvtVector4R& rs,const EvtVector4R& p4);
  friend EvtVector4R boostTo(const EvtVector4R& rs,const EvtVector3R& boost);

  inline friend EvtVector4R operator*(double d,const EvtVector4R& v2);
  inline friend EvtVector4R operator*(const EvtVector4R& v2,double d);
  inline friend EvtVector4R operator/(const EvtVector4R& v2,double d);
  inline friend double operator*(const EvtVector4R& v1,const EvtVector4R& v2);
  inline friend EvtVector4R operator+(const EvtVector4R& v1,const EvtVector4R& v2);
  inline friend EvtVector4R operator-(const EvtVector4R& v1,const EvtVector4R& v2);

public:
  EvtVector4R(){}
  EvtVector4R(const double& e,const double& p1,const double& p2,const double& p3);
  inline void set(int i,double d);
  inline void set(double e,double px,double py ,double pz);
  inline EvtVector4R& operator*=(double c);
  inline EvtVector4R& operator/=(double c);
  inline EvtVector4R& operator=(const EvtVector4R& v2);
  inline EvtVector4R& operator+=(const EvtVector4R& v2);
  inline EvtVector4R& operator-=(const EvtVector4R& v2);
  inline double get(int i) const;
  inline double cont(const EvtVector4R& v4) const;
  friend std::ostream& operator<<(std::ostream& s, const EvtVector4R& v);
  double mass2() const;
  double mass() const;
  void applyRotateEuler(const double &alpha, const double &beta, const double &gamma);
  void applyBoostTo(const EvtVector4R& p4);
  void applyBoostTo(const EvtVector3R& boost);
  double px(void) const {return get(1);}
  double py(void) const {return get(2);}
  double pz(void) const {return get(3);}
  double e(void) const {return get(0);}

  void px(const double& x) {return set(1,x);}
  void py(const double& x) {return set(2,x);}
  void pz(const double& x) {return set(3,x);}
  void e(const double& x)  {return set(0,x);}

  EvtVector4R cross(const EvtVector4R& v2) const;
  double dot(const EvtVector4R& v2) const;
  double d3mag() const;
private:
  double v[4];
};

inline EvtVector4R& EvtVector4R::operator=(const EvtVector4R& v2){
  v[0]=v2.v[0];
  v[1]=v2.v[1];
  v[2]=v2.v[2];
  v[3]=v2.v[3];
  return *this;
}

inline EvtVector4R& EvtVector4R::operator+=(const EvtVector4R& v2){
  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  v[3]+=v2.v[3];
  return *this;
}

inline EvtVector4R& EvtVector4R::operator-=(const EvtVector4R& v2){
  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  v[3]-=v2.v[3];
  return *this;
}

inline double EvtVector4R::mass2() const{
  return v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3];
}

inline EvtVector4R operator*(double c,const EvtVector4R& v2){
  return EvtVector4R(v2)*=c;
}

inline EvtVector4R operator*(const EvtVector4R& v2,double c){
  return EvtVector4R(v2)*=c;
}

inline EvtVector4R operator/(const EvtVector4R& v2,double c){
  return EvtVector4R(v2)/=c;
}

inline EvtVector4R& EvtVector4R::operator*=(double c){
  v[0]*=c;
  v[1]*=c;
  v[2]*=c;
  v[3]*=c;
  return *this;
}

inline EvtVector4R& EvtVector4R::operator/=(double c){
  double cinv=1.0/c;
  v[0]*=cinv;
  v[1]*=cinv;
  v[2]*=cinv;
  v[3]*=cinv;
  return *this;
}

inline double operator*(const EvtVector4R& v1,const EvtVector4R& v2){
  return v1.v[0]*v2.v[0]-v1.v[1]*v2.v[1]-
         v1.v[2]*v2.v[2]-v1.v[3]*v2.v[3];
}

inline double EvtVector4R::cont(const EvtVector4R& v4) const {
  return v[0]*v4.v[0]-v[1]*v4.v[1]-
         v[2]*v4.v[2]-v[3]*v4.v[3];
}

inline EvtVector4R operator-(const EvtVector4R& v1,const EvtVector4R& v2){
  return EvtVector4R(v1)-=v2;
}

inline EvtVector4R operator+(const EvtVector4R& v1,const EvtVector4R& v2){
  return EvtVector4R(v1)+=v2;
}

inline double EvtVector4R::get(int i) const {
  return v[i];
}

inline void EvtVector4R::set(int i,double d){
  v[i]=d;
}

inline void EvtVector4R::set(double e,double p1,double p2, double p3){
  v[0]=e;
  v[1]=p1;
  v[2]=p2;
  v[3]=p3;
}

#endif // EVTVECTOR4R_H
