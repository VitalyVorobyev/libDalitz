#include "virtualresff.h"
#include <math.h>

VirtualResFF::VirtualResFF(const double& _r, const double& _p0):
  FormFactor(_r,_p0)
{
}

double VirtualResFF::operator()(const double& p) const {
  return exp(-r()*(p-p0()));
}
