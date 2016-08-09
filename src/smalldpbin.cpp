#include "smalldpbin.h"

using namespace std;

SmallDPBin::SmallDPBin(const double &mAB, const double &mAC):
  m_val(0), m_valb(0),
  m_wght(0), m_dph(0), m_p(0), m_pb(0),
  m_h(0), m_gsize(0), m_norm(1),
  m_dpsize(0), m_dpmin(0),
  m_mAB(mAB), m_mAC(mAC),
  m_ABbin(GetBin(m_mAB)),
  m_ACbin(GetBin(m_mAC))
{
}

SmallDPBin::SmallDPBin(DalitzModel* model, const double& mAB, const double& mAC):
  SmallDPBin(mAB,mAC)
{
  SetModel(model);
  Calc();
}

SmallDPBin::SmallDPBin(const SmallDPBin& x){
  *this = x;
}

//SmallDPBin& SmallDPBin::operator=(const SmallDPBin& x){
//  m_model  = x.m_model;
//  m_mAB    = x.m_mAB;
//  m_mAC    = x.m_mAC;
//  m_val    = x.m_val;
//  m_h      = x.m_h;
//  m_gsize  = x.m_gsize;
//  m_norm   = x.m_norm;
//  m_dpsize = x.m_dpsize;
//  return *this;
//}

int SmallDPBin::SetModel(DalitzModel* model){
  m_model  = model;
  m_dpmin  = m_model->mABsq_min();
  m_dpsize = m_model->mABsq_max() - m_dpmin;
  return 0;
}

int SmallDPBin::SetGridSize(const int gsize){
  m_gsize = gsize;
  m_h = m_dpsize/m_gsize;
  return 0;
}

int SmallDPBin::SetGridStep(const double& gstep){
  m_h = gstep;
  m_gsize = (int)(m_dpsize/m_h);
  return 0;
}

int SmallDPBin::SetPoint(const double& mAB, const double& mAC){
  if(!m_model->IsInPlot(mAB,mAC)) return -1;
  m_mAB = mAB; m_mAC = mAC;
  m_ABbin = GetBin(m_mAB);
  m_ACbin = GetBin(m_mAC);
  Calc();
  return 0;
}

int SmallDPBin::SetBin(const int i, const int j){
  m_ABbin = i; m_ACbin = j;
  m_mAB = GetVal(m_ABbin);
  m_mAC = GetVal(m_ACbin);
  if(!m_model->IsInPlot(m_mAB,m_mAC)) return -1;
  Calc();
  return 0;
}

int SmallDPBin::SetNorm(const double& norm){
  m_norm = norm;
  return 0;
}

bool SmallDPBin::IsInBin(const double& mAB, const double& mAC) const{
  if(fabs(mAB-m_mAB)>m_h) return false;
  if(fabs(mAC-m_mAC)>m_h) return false;
  return true;
}

double SmallDPBin::Phase(void)  const{ return arg(m_val);}
double SmallDPBin::Amp(void)    const{ return sqrt(m_p);}
double SmallDPBin::Ampb(void)   const{ return sqrt(m_pb);}
double SmallDPBin::P(void)      const{ return m_p;}
double SmallDPBin::Pb(void)     const{ return m_pb;}
compld SmallDPBin::CAmp(void)   const{ return m_val;}
compld SmallDPBin::CAmpb(void)  const{ return m_valb;}
int SmallDPBin::GSize(void)     const{ return m_gsize;}
double SmallDPBin::Weight(void) const{ return m_wght;}

int SmallDPBin::GetCurrentPoint(double& mAB, double& mAC) const{
  mAB = m_mAB; mAC = m_mAC;
  return 0;
}

int SmallDPBin::GetBin(const double& x) const{
  if(x<= 0 || x>(m_dpmin+m_dpsize) || x<m_dpmin) return -1;
  return (x-m_dpmin)/m_h;
}

double SmallDPBin::GetVal(const int bin) const{
  if(bin<0 || bin>=m_gsize) return 0;
  return m_dpmin + m_h*((double)bin+0.5);
}

int SmallDPBin::Calc(void){
  if(!m_model->IsInPlot(m_mAB,m_mAC)) return -1;
  m_val  = m_model->Amp(m_mAB,m_mAC);
  m_valb = m_model->Amp(m_mAC,m_mAB);
  m_dph  = arg(m_valb/m_val);
  m_p    = norm(m_val);
  m_pb   = norm(m_valb);
//  m_wght = log(1.+2.*sqrt(m_p*m_pb)*(m_p+m_pb));
  m_wght = sqrt(m_p*m_pb)*(m_p+m_pb);
  return 0;
}

