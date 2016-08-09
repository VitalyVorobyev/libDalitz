#include "symdalitzmodel.h"
#include "consts.h"

#include <iostream>

using namespace std;

const int sign = -1;

SymDalitzModel::SymDalitzModel(const double& mmo,const double& mcha,const double& mchb,const double& delmin,const double& delmax):
DalitzModel(mmo,mcha,mchb,mchb),del_min(delmin),del_max(delmax),nbins(8)
{
}

int SymDalitzModel::GetBin(const double& mp, const double& mm){
  if(!IsInPlot(mp,mm)) return 0;
  double delt = delta(mp,mm);
  for(int i=1; i<=nbins; i++){
    if(2.*M_PI*(i-1.5)/nbins < delt && delt < 2.*M_PI*(i-0.5)/nbins){
      return mm>mp ? i : -i;
    }
  }
  cout << "GetBin: delta = " << delt << ", mp: " << mp << ", mm: " << mm << endl;
  return 0;
}

void SymDalitzModel::PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta){
  compld A  = Amp(mp,mm);
  compld Ab = Amp(mm,mp);
  delta = sign*(arg(Ab) - arg(A));
  if(delta<del_min) delta += 2.*M_PI;
  if(delta>del_max) delta -= 2.*M_PI;
  P = norm(A); Pbar = norm(Ab);
  return;
}

double SymDalitzModel::delta(const double& mp,const double& mm){
  double del = sign*(arg(Amp(mm,mp)) - arg(Amp(mp,mm)));
  if(std::isnan(del)){
    cout << "delta: " << del << ", mp:" << mp << ", mm: " << mm << endl;
    return 0;
  }
  if(del<del_min) return del + 2.*M_PI;
  if(del>del_max) return del - 2.*M_PI;
  return del;
}
