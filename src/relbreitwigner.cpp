#include "relbreitwigner.h"

RelBreitWigner::RelBreitWigner(const double &G0, const double &m, const double &p0, const int mom, const bool constwidth):
  AbsPropagator(m,p0,constwidth ? (AbsVarWidth*) new ConstWidth(G0) : new BWWidth(G0,m,p0,mom)),
  m_const_width(constwidth)
{
}

EvtComplex RelBreitWigner::operator()(const double& s, const double& p) const{
  const EvtComplex ione(0,1);
  const BWWidth& G = *((BWWidth*)width());
  const double& mass = m();
  return 1./(mass*mass-s-ione*mass*G(s,p));
}

RelBreitWigner::~RelBreitWigner(){
  if(m_const_width) delete (ConstWidth*)width();
  else              delete (BWWidth*)width();
}
