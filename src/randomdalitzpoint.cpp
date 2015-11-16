#include "randomdalitzpoint.h"

RandomDalitzPoint::RandomDalitzPoint(const DalitzPhaseSpace &phsp):
 DalitzPhaseSpace(phsp)
{
  init();
}

RandomDalitzPoint::RandomDalitzPoint(const double& mmo, const double& mca, const double& mcb, const double& mcc):
  DalitzPhaseSpace(mmo,mca,mcb,mcc)
{
  init();
}

RandomDalitzPoint::~RandomDalitzPoint(){
  delete unifAB;
  delete unifAC;
}

void RandomDalitzPoint::GetPoint(double& mAB, double& mAC){
  do{
    mAB = (*unifAB)(re); mAC = (*unifAC)(re);
  } while(!IsInPlot(mAB,mAC));
  return;
}

void RandomDalitzPoint::SetSeed(const int seed){
  m_seed = seed;
  re.seed(m_seed);
}

void RandomDalitzPoint::init(void){
  unifAB = new std::uniform_real_distribution<double>(mAB_min(),mAB_max());
  unifAC = new std::uniform_real_distribution<double>(mAC_min(),mAC_max());
  m_seed = chrono::system_clock::now().time_since_epoch().count();
  SetSeed(m_seed);
}
