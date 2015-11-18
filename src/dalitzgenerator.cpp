#include "dalitzgenerator.h"

using namespace std;

DalitzGenerator::DalitzGenerator(const DalitzModel &_dm):
  RandomDalitzPoint(_dm),m_ntries(1e8),m_maj_counts(1e4),m_maj(-1)
{
  dm = const_cast<DalitzModel*>(&_dm);
}

double DalitzGenerator::CalcMajorant(void){
  double cur_val;
  double mAB,mAC;
  for(int i=0; i<m_maj_counts; i++){
    GetPoint(mAB,mAC);
    cur_val = dm->P(mAB,mAC);
    if(cur_val>m_maj) m_maj = cur_val;
  }
  m_maj *= 1.05;
  cout << "Majorant is set equals " << m_maj << endl;
  return m_maj;
}

int DalitzGenerator::Generate(const int NEv,vector<double>& mABv,vector<double>& mACv){
  if(m_maj<0) CalcMajorant();
  mABv.clear(); mACv.clear();
  uniform_real_distribution<double> unifMaj(0.,m_maj);
  default_random_engine ren;

  long tries = 0;
  int Evtn = 0;
  double mAB,mAC,xi;
  cout << "Generating " << NEv << " events..." << endl;
  while(Evtn < NEv && tries < m_ntries){
    GetPoint(mAB,mAC);
    xi = unifMaj(ren);
    if(xi<dm->P(mAB,mAC)){
      mABv.push_back(mAB);
      mACv.push_back(mAC);
      Evtn++;
      if(!(Evtn%1000)) cout << Evtn << " events" << endl;
    }
    tries++;
  }
  const double Eff = (double)Evtn/tries;
  if(tries == m_ntries){
    cout << "DalitzGenerator: tries limit exceed!";
    cout << "  efficiency: " << Eff << endl;
    cout << "  events generated: " << Evtn << endl;
    cout << "  tries limit: " << m_ntries << endl;
    return -1;
  }
  cout << "Done!. Efficiency " << Eff << endl;
  return 0;
}
