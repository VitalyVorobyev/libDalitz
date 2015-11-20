#include "symdalitzmodel.h"

SymDalitzModel::SymDalitzModel(const double& mmo,const double& mcha,const double& mchb,const double& delmin,const double& delmax):
DalitzModel(mmo,mcha,mchb,mchb),del_min(delmin),del_max(delmax),nbins(8)
{
}

int SymDalitzModel::GetBin(const double& mp, const double& mm){
  if(!IsInPlot(mp,mm)) return 0;
  const double delt = mp>mm ? delta(mp,mm) : delta(mm,mp);
  for(int i=1; i<=nbins; i++){
    if(2.*M_PI*(i-1.5)/nbins < delt && delt < 2.*M_PI*(i-0.5)/nbins){
      return mp>mm ? i : -i;
    }
  }
  std::cout << "GetBin: delta = " << delt << ", mp: " << mp << ", mm: " << mm << std::endl;
  return 0;
}

void SymDalitzModel::PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta){
  EvtComplex A  = Amp(mp,mm);
  EvtComplex Ab = Amp(mm,mp);
  delta = arg(A) - arg(Ab);
  if(delta<del_min) delta += 2.*M_PI;
  if(delta>del_max) delta -= 2.*M_PI;
  P = abs2(A); Pbar = abs2(Ab);
  return;
}

double SymDalitzModel::delta(const double& mp,const double& mm){
  double del = - (arg(Amp(mm,mp)) - arg(Amp(mp,mm)));
  if(std::isnan(del)){
    std::cout << "delta: " << del << ", mp:" << mp << ", mm: " << mm << std::endl;
  }
  if(del<del_min) return del + 2.*M_PI;
  if(del>del_max) return del - 2.*M_PI;
  return del;
}
