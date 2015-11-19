#include "buggpropagator.h"
#include "buggwidth.h"

BuggPropagator::BuggPropagator(void):
  AbsPropagator(0,0),m_width(new BuggWidth())
{
}

EvtComplex BuggPropagator::operator()(const double& s, const double& p) const{
  double G1,GTot;
  m_width->GetWidths(s,G1,GTot);
  const double sA   = m_width->sA();
  const double mrsq = m_width->mrsq();
  const double z    = m_width->z();
  const double g1sq = m_width->g1sq();
  const EvtComplex ione = EvtComplex(0,1);

  return G1/(mrsq-s-g1sq*(s-sA)*z/(mrsq-sA)-ione*sqrt(mrsq)*GTot);
}

BuggPropagator::~BuggPropagator(void){
  delete m_width;
  return;
}
