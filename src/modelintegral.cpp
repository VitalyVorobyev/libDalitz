#include "modelintegral.h"
#include "eqphasebin.h"

#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

ModelIntegral::ModelIntegral(const AbsDalitzModel *model,const int NBins, const int gridsize):
  m_model(model), m_nbins(NBins), m_gsize(gridsize)
{
}

int ModelIntegral::Calculate(cstr& label,vectd& C,vectd& S,vectd& Kp,vectd& Kn){
   C.resize(m_nbins,0);
   S.resize(m_nbins,0);
  Kp.resize(m_nbins,0);
  Kn.resize(m_nbins,0);
  EqPhaseBin eqbin(m_model,m_nbins);

  double mmMax;
  double mmMin;
  cdouble M_min = m_model->mABsq_min();
  cdouble M_max = m_model->mABsq_max();
  cdouble dm = (M_max - M_min)/m_gsize;
  const str fname = str("params/") + label + str("_binning.txt");
  ofstream ofile(fname.c_str());
  ofile << "GridSize " << m_gsize << endl;
  ofile << m_model->mM() << " -> " << m_model->mA();
  ofile << " " << m_model->mB();
  ofile << " " << m_model->mC() << endl;
  ofile << "mAB: " << m_model->mABsq_min() << " " << m_model->mABsq_max() << endl;
  ofile << "mAC: " << m_model->mACsq_min() << " " << m_model->mACsq_max() << endl;
  ofile << "mBC: " << m_model->mBCsq_min() << " " << m_model->mBCsq_max() << endl;
  ofile << m_model->ABaxis() << endl;
  ofile << m_model->ACaxis() << endl;
  ofile << m_model->BCaxis() << endl;
  double mp = M_min;
  for(int i=0; i<m_gsize; i++){
    mp += dm;
    m_model->mABsqRange_AC(mp,mmMin,mmMax);
    if(mp>mmMax) break;
    double mm = mp;
    for(int j=0; mm<mmMax; j++){
      mm += dm;
      if(!m_model->IsInPlot(mp,mm)) continue;
      int bin = abs(eqbin.Bin(mp,mm));
      ofile << mp << " " << mm << " " << m_model->GetmBCsq(mp,mm) << " " << bin << endl;
      double delta, P, Pbar;
      if(PPbarDelta(mp,mm,P,Pbar,delta)) continue;
      bin -= 1;

      C[bin]  += sqrt(P*Pbar)*cos(delta);
      S[bin]  += sqrt(P*Pbar)*sin(delta);
      Kp[bin] += P;
      Kn[bin] += Pbar;
    }
  }
  ofile.close();

  double norm = 0;
  for(int i=0; i<8; i++){
    C[i] /= sqrt(Kp[i]*Kn[i]);
    S[i] /= sqrt(Kp[i]*Kn[i]);
    norm += Kp[i] + Kn[i];
  }

  cout << "Hello, God! " << label << endl;
  cstr fname1("params/test.txt");//= str("params/") + label + str("_eq_phase_params.txt");
  cout << fname1 << endl;
  ofstream ofile1(fname1.c_str());
  cout << "Norm     = " << norm << endl;
  for(int i=0; i<m_nbins; i++){
    Kp[i] /= norm;
    Kn[i] /= norm;

    cout << i+1;
    cout << ": C = "  << C[i];
    cout << ", S = "  << S[i];
    cout << ", Kp = " << Kp[i];
    cout << ", Kn = " << Kn[i];
    cout << ", Q = "  << C[i]*C[i]+S[i]*S[i] << endl;

    ofile1 << i+1;
    ofile1 << ": C = "  << C[i];
    ofile1 << ", S = "  << S[i];
    ofile1 << ", Kp = " << Kp[i];
    ofile1 << ", Kn = " << Kn[i];
    ofile1 << ", Q = "  << C[i]*C[i]+S[i]*S[i] << endl;
  }
  ofile1.close();
  return 0;
}

int ModelIntegral::PPbarDelta(cdouble& mABsq, cdouble& mACsq, double& P, double& Pbar, double& delta){
  compld A  = m_model->Amp(mABsq,mACsq);
  compld Ab = m_model->Amp(mACsq,mABsq);
  delta = -(std::arg(Ab) - std::arg(A));
  if(std::isnan(delta)){ cout << "delta: " << delta << endl; return -1;}
  P = norm(A); Pbar = norm(Ab);
  return 0;
}
