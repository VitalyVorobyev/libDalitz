#include "dalitzvetodst2010.h"
#include "math.h"

DalitzVetoDst2010::DalitzVetoDst2010(const double cut):
  m_mdst(2.0103),m_cut(cut)
{
}

bool DalitzVetoDst2010::operator()(const double& AB,const double& BC) const {
  return fabs(sqrt(AB)-m_mdst)<m_cut;
}
