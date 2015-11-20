#include "virtualresff.h"
#include <math.h>

VirtualResFF::VirtualResFF(const double& _r, const double& _p0):
  FormFactor(_r,_p0)
{
}

double VirtualResFF::operator()(const double& p) const {
  const double _p0 = sqrt(p0sq());
  const double _p = sqrt(p);
  return exp(-r()*(_p-_p0));
}
