#include "virtualresff.h"
#include <math.h>

VirtualResFF::VirtualResFF(const double& _r, const double& _p0):
  FormFactor(_r,_p0)//,bwff(new BlattWeisskopf(1,_p0,FFType::FFResonance))
{
}

double VirtualResFF::operator()(const double& p) const {
  return exp(-1.*fabs(p-p0()));//*bwff->operator()(p);
}
