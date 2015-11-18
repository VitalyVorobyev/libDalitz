#include "buggpropagator.h"
#include "buggwidth.h"

BuggPropagator::BuggPropagator(void):
  AbsPropagator(0,0,new BuggWidth())
{
}

EvtComplex BuggPropagator::operator()(const double& s, const double& p) const{
  BuggWidth* wdth = (BuggWidth*)width();
  double G1,GTot;
  wdth->GetWidths(s,G1,GTot);
  const double sA   = wdth->sA();
  const double mrsq = wdth->mrsq();
  const double z    = wdth->z();
  const double g1sq = wdth->g1sq();
  const EvtComplex ione = EvtComplex(0,1);

  return G1/(mrsq-s-g1sq*(s-sA)*z/(mrsq-sA)-ione*sqrt(mrsq)*GTot);
}
