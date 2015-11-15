#include "symdalitzmodel.h"
using namespace std;
//const int BINS[8] = {8,7,6,5,4,3,2,1};
//const int BINS[8] = {8,1,2,3,4,5,6,7};
//const int BINS[8] = {1,2,3,4,5,6,7,8};

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
  cout << "GetBin: delta = " << delt << ", mp: " << mp << ", mm: " << mm << endl;
  return 0;
}

void SymDalitzModel::PPbarDelta(const double& mp, const double& mm, double& P, double& Pbar, double& delta){
  EvtVector4R pM;
  EvtVector4R pA;
  EvtVector4R pB;
  EvtVector4R pC;
  GetLVs(mp,mm,pM,pA,pB,pC);
  EvtComplex A  = Amp(pM,pA,pB,pC);
  EvtComplex Ab = Amp(pM,pA,pC,pB);
  delta = arg(A) - arg(Ab);
  if(delta<del_min) delta += 2.*M_PI;
  if(delta>del_max) delta -= 2.*M_PI;
  P = abs2(A); Pbar = abs2(Ab);
  return;
}

double SymDalitzModel::delta(const double& mp,const double& mm){
  EvtVector4R pM;
  EvtVector4R pA;
  EvtVector4R pB;
  EvtVector4R pC;
  GetLVs(mp,mm,pM,pA,pB,pC);
  double del = - (arg(Amp(pM,pA,pC,pB)) - arg(Amp(pM,pA,pB,pC)));
  if(std::isnan(del)){
    cout << "delta: " << del << ", mp:" << mp << ", mm: " << mm << endl;
    cout << " pM:" << pM << endl;
    cout << " pA:" << pA << endl;
    cout << " pB:" << pB << endl;
    cout << " pC:" << pC << endl;
  }
  if(del<del_min) return del + 2.*M_PI;
  if(del>del_max) return del - 2.*M_PI;
  return del;
}
