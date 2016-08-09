#ifndef ABSPROPAGATOR_H
#define ABSPROPAGATOR_H

#include "consts.h"

class ResPropType{
public:
  static const int NR       = 0;
  static const int RBW      = 1;
  static const int GS       = 2;
  static const int RhoOmega = 3;
  static const int Bugg     = 4;
  static const int VDst     = 5;
  static const int Flatte   = 6;
  static const int VDst2    = 7;
};

class AbsPropagator{
public:
//  AbsPropagator(const double &m, const double &p0, AbsVarWidth* width);
  AbsPropagator(const double &m, const double &p0);
  virtual ~AbsPropagator() {}
  virtual compld operator()(const double& s, const double& p) const = 0;

  double m(void)  const {return m_m;}
  double p0(void) const {return m_p0;}
//  AbsVarWidth* width(void) const {return m_width;}

private:
  double m_m;
  double m_p0;
//  AbsVarWidth* m_width;
};

#endif // ABSPROPAGATOR_H
