#include "dalitzgenerator.h"

#include <iostream>

using namespace std;
typedef DalitzGenerator DGen;

DalitzGenerator::DalitzGenerator(const AbsDalitzModel *model):
  RandomDalitzPoint(model),
  m_ntries(1e8),
  m_maj_counts(1e4),
  m_maj(-1),
  m_model(model),
  ren(new rndmeng())
{
  SetMajorant(CalcMajorant());
}

void DGen::SetMajorant(cdouble& p){
  m_maj = p;
  cout << "DalitzGenerator: majorant is set to " << m_maj << endl;
}

double DGen::CalcMajorant(void) const{
  cout << "CalcMajorant..." << endl;
  double cur_val;
  double maj = 0;
  double mABsq,mACsq;
  for(int i=0; i<m_maj_counts; i++){
    GetPoint(mABsq,mACsq);
    cur_val = m_model->P(mABsq,mACsq);
    if(cur_val>maj) maj = cur_val;
  }
  maj *= 1.05;
//  cout << "Majorant is set to " << maj << endl;
  return maj;
}

int DGen::Generate(double &mABsq, double &mACsq) const{
  U64 tries = 0;
  double xi = 0;
  unifdist unifMaj(0.,m_maj);
  while(tries < m_ntries){
    tries++;
    GetPoint(mABsq,mACsq);
    xi = unifMaj(*ren);
    if(xi < m_model->P(mABsq,mACsq)) break;
  }
  if(tries == m_ntries){
    cout << "DalitzGenerator::Generate: tries limit exceeded!" << endl;
    return 0;
  }
  return tries;
}

int DGen::Generate(const int NEv,vectd& mABv,vectd& mACv,const bool silent) const{
  mABv.clear(); mACv.clear();

  U64 tries = 0;
  int Evtn = 0;
  double mAB,mAC;
  cout << "Generating " << NEv << " events..." << endl;
  while(Evtn < NEv){
    const int ntries = Generate(mAB,mAC);
    if(!ntries) break;
    tries += ntries;
    mABv.push_back(mAB);
    mACv.push_back(mAC);
    Evtn++;
    if(!(Evtn%1000)) cout << Evtn << " events" << endl;
  }
  if(!silent){
    cdouble Eff = (double)Evtn/tries;
    if(tries == m_ntries){
      cout << "DalitzGenerator: tries limit exceeded!";
      cout << "  efficiency: " << Eff << endl;
      cout << "  events generated: " << Evtn << endl;
      cout << "  tries limit: " << m_ntries << endl;
      return -1;
    }
    cout << "Done!. Efficiency " << Eff << endl;
  }
  return 0;
}
