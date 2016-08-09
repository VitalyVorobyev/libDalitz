#include "eqphasebin.h"

#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

EqPhaseBin::EqPhaseBin(const AbsDalitzModel* model,const int NBins):
  m_model(model), m_nbins(NBins),
  m_del_min(lEdge(1)), m_del_max(rEdge(8))
{
}

int EqPhaseBin::Bin(cdouble& mABsq, cdouble& mACsq){
  if(!m_model->IsInPlot(mABsq,mACsq)) return 0;
  cdouble delt = delta(mABsq,mACsq);
  for(int i=1; i<=m_nbins; i++){
    if(lEdge(i) < delt && delt < rEdge(i)){
      return (mABsq>mACsq ? i : -i);
    }
  }
  return 0;
}

double EqPhaseBin::lEdge(const int i){ return 2.*M_PI*(i-1.5)/m_nbins;}
double EqPhaseBin::rEdge(const int i){ return 2.*M_PI*(i-0.5)/m_nbins;}

double EqPhaseBin::delta(cdouble& mABsq, cdouble& mACsq){
  double del = -(m_model->Arg(mACsq,mABsq) - m_model->Arg(mABsq,mACsq));
  if(std::isnan(del)){ cout << "delta: " << del << endl; return 0; }
  if(del<m_del_min) return del + 2.*M_PI;
  if(del>m_del_max) return del - 2.*M_PI;
  return del;
}
