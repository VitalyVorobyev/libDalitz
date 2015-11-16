#include "modelintegral.h"

ModelIntegral::ModelIntegral(SymDalitzModel *model){
  m_model = model;
  m_nbins = model->GetNBins();
  m_gsize = 500;
}

double ModelIntegral::Calculate(const string& label,vector<double>& C,vector<double>& S,vector<double>& K,vector<double>& Kb){
  C.clear();
  S.clear();
  K.clear();
  Kb.clear();
  for(int i=0; i<m_nbins; i++){
    C.push_back(0);
    S.push_back(0);
    K.push_back(0);
    Kb.push_back(0);
  }
  double mmMax;
  double mmMin;
  const double M_min = m_model->mAB_min();
  const double M_max = m_model->mAB_max();
  const double dm = (M_max - M_min)/m_gsize;
  double mp = M_min;
  m_majorant = 0;
  const string fname = label + string("_binning.txt");
  ofstream ofile(fname.c_str());
  ofile << "GridSize " << m_gsize << endl;
  ofile << m_model->mM() << " -> " << m_model->mA();
  ofile << " " << m_model->mB();
  ofile << " " << m_model->mC() << endl;
  ofile << "mAB: " << m_model->mAB_min() << " " << m_model->mAB_max() << endl;
  ofile << "mAC: " << m_model->mAC_min() << " " << m_model->mAC_max() << endl;
  ofile << "mBC: " << m_model->mBC_min() << " " << m_model->mBC_max() << endl;
  ofile << m_model->ABaxis() << endl;
  ofile << m_model->ACaxis() << endl;
  ofile << m_model->BCaxis() << endl;
  for(int i=0; i<m_gsize; i++){
    mp += dm;
    m_model->mAB_range(mp,mmMin,mmMax);
    if(mp>mmMax) break;
//    mmMin = mmMin > mp ? mmMin : mp;
    double mm = mmMin;
    for(int j=0; mm<mmMax; j++){
      mm += dm;
      if(!m_model->IsInPlot(mp,mm)) continue;
      int bin = abs(m_model->GetBin(mp,mm));
      ofile << mp << " " << mm << " " << m_model->mBC(mp,mm) << " " << bin << endl;
      double delta, P, Pbar;
      m_model->PPbarDelta(mp,mm,P,Pbar,delta);
      if(std::isnan(delta)) continue;
      if(P+Pbar > m_majorant) m_majorant = P+Pbar;
      bin -= 1;
      C[bin]  += sqrt(P*Pbar)*cos(delta);
      S[bin]  += sqrt(P*Pbar)*sin(delta);
      K[bin]  += P;
      Kb[bin] += Pbar;
    }
  }
  ofile.close();
  double norm = 0;
  for(int i=0; i<8; i++){
    C[i] /= sqrt(K[i]*Kb[i]);
    S[i] /= sqrt(K[i]*Kb[i]);
    norm += K[i] + Kb[i];
  }

  const string fname1 = label + string("_eq_phase_params.txt");
  ofstream ofile1(fname1.c_str());
  cout << "Norm     = " << norm << endl;
  cout << "Majorant = " << m_majorant << endl;
  for(int i=0; i<m_nbins; i++){
    K[i]  /= norm;
    Kb[i] /= norm;
    cout << i+1 << ": C = " << C[i] << ", S = " << S[i] << ", K = " << K[i] << ", Kb = " << Kb[i] << ", Q = " << C[i]*C[i]+S[i]*S[i] << endl;
    ofile1 << i+1 << ": C = " << C[i] << ", S = " << S[i] << ", K = " << K[i] << ", Kb = " << Kb[i] << ", Q = " << C[i]*C[i]+S[i]*S[i] << endl;
  }
  ofile.close();
  return m_majorant;
}
