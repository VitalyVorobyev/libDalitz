#include "relbreitwigner.h"
#include "bwwidth.h"
#include "constwidth.h"
#include "flattewidth.h"

//const int VarWType::Const  = 0;
//const int VarWType::BW     = 1;
//const int VarWType::GS     = 2;
//const int VarWType::Flatte = 3;
//const int VarWType::Bugg   = 4;

RelBreitWigner::RelBreitWigner(const double &G0, const double &m, const double &p0, const int mom, const int wtype):
  AbsPropagator(m,p0), m_wtype(wtype)
{
  switch(m_wtype){
  case VarWType::Const:
    m_width = new ConstWidth(G0);
//    std::cout << "ConstWisdth is defined " << m_wtype << " " << wtype << std::endl;
    break;
  case VarWType::BW:
    m_width = new BWWidth(G0,m,p0,mom);
//    std::cout << "BWWidth is defined " << m_wtype << std::endl;
    break;
  case VarWType::Flatte:
    m_width = new FlatteWidth(m);
    break;
  default:
    std::cout << "RelBreitWigner: wrong VarWType " << m_wtype << std::endl;
    break;
  }
}

EvtComplex RelBreitWigner::operator()(const double& s, const double& p) const{
  const EvtComplex ione(0,1);
  const double& mass = m();
//  std::cout << "RBW: m: " << mass << " s: " << s << " width: " << (*m_width)(s,p) << std::endl;
  return 1./(mass*mass-s-ione*mass*(*m_width)(s,p));
}

RelBreitWigner::~RelBreitWigner(){
  delete m_width;
  return;
}
