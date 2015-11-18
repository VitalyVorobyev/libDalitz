#ifndef BWWIDTH_H
#define BWWIDTH_H

#include "absvarwidth.h"
#include "blattweisskopf.h"

// Eq.(10) in Phys. Rev D92, 032002 (2015)

class BWWidth : public AbsVarWidth{
public:
  BWWidth(const double& G0, const double& m, const double& p0, const int mom);

  double operator()(const double& s, const double& p) const;

private:
  int m_mom;
  BlattWeisskopf* m_ff;
  double m_precalc;
};

#endif // BWWIDTH_H
