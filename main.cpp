#include "libDalitz/libdalitz.h"
#include "drawbdparams.h"

int main(int argc, char** argv){
  KspipiModel* kspipimodel = new KspipiModel();
  B0toD0pipiModel* d0pipimodel = new B0toD0pipiModel();
  DrawBDParams drawpar;
  ModelIntegral kspipiintegral(kspipimodel);
  kspipiintegral.SetGridSize(1000);

  vector<double> C,S,K,Kb;
  const string fname("kspipi_binning.txt");
  kspipiintegral.Calculate(fname,C,S,K,Kb);
  const string dpname("kspipi_binning");
  drawpar.DrawBinsmABmAC(fname,dpname);
  drawpar.DrawBinsmABmBC(fname,dpname);
//  const string csname("kspipi_cs_model");
//  drawpar.DrawCS(C,S,csname);
//  const string kname("kspipi_k_model");
//  drawpar.DrawK(K,Kb,kname);

  ModelIntegral d0pipiintegral(d0pipimodel);
  d0pipiintegral.SetGridSize(1000);
  vector<double> Cd,Sd,Kd,Kbd;
  const string dfname("d0pipi_binning.txt");
  d0pipiintegral.Calculate(dfname,Cd,Sd,Kd,Kbd);
  const string ddpname("d0pipi_binning");
  drawpar.DrawBinsmABmAC(dfname,ddpname);
  drawpar.DrawBinsmABmBC(dfname,ddpname);

  return 0;
}
