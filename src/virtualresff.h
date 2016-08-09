#ifndef VIRTUALRESFF_H
#define VIRTUALRESFF_H

#include "formfactor.h"
//#include "blattweisskopf.h"

///
/// \brief The VirtualResFF class
///
class VirtualResFF : public FormFactor{
public:
  VirtualResFF(const double& _r, const double& _p0);
  double operator()(const double& p) const;

private:
//  BlattWeisskopf* bwff;
};

#endif // VIRTUALRESFF_H
