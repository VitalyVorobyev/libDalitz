#include "dalitzmcintegral.h"
#include <math.h>
#include <iomanip>

DalitzMCIntegral::DalitzMCIntegral(const DalitzModel &model):
  RandomDalitzPoint(model),m_ncounts(1e6),m_int(0),m_err(0)
{
  m_model = const_cast<DalitzModel*>(&model);
}

int DalitzMCIntegral::CalcIntegral(double &val, double &err, const long& nc){
  if(nc) m_ncounts = nc;
  val = 0; err = 0;
  double mAB,mAC;
  std::cout << "DalitzMCIntegral: integrating with " << m_ncounts << " points..." << std::endl;
  for(int i=0; i<m_ncounts; i++){
    if(!(i%100000)) std::cout << i << " counts val = " << val << " val2 = " << err << std::endl;
    GetPoint(mAB,mAC);
    const double res = m_model->P(mAB,mAC);
    val += res/m_ncounts;
    err += res*res/m_ncounts;
  }
  err = sqrt((err - val*val)/m_ncounts);

  m_int = val; m_err = err;
  std::cout << "Done! Int = " << m_int << " +- " << m_err << " (" << 100.*m_err/m_int << "%)" << std::endl;
  return 0;
}

double DalitzMCIntegral::CalcBranchings(std::vector<double>& brvec, std::vector<double> &brerr, const long& nc){
  if(nc) m_ncounts = nc;
  const int NRes = m_model->ResNum();
  brvec.clear(); brerr.clear();
  std::vector<EvtComplex> amps;
  for(int i=0; i<NRes; i++){ brvec.push_back(0); brerr.push_back(0);}

  double mAB,mAC;
  double totint = 0;
  double toterr = 0;

  std::cout << "DalitzMCIntegral: calculating branchings with " << m_ncounts << " points..." << std::endl;
  for(int i=0; i<m_ncounts; i++){
    if(!(i%100000)) std::cout << i << " counts" << std::endl;
    GetPoint(mAB,mAC);
    const double res = abs2(m_model->GetAmplitudes(amps,mAB,mAC));
    totint += res/m_ncounts;
    toterr += res*res/m_ncounts;
    for(int j=0; j<NRes; j++){
      const double res2 = abs2(amps[j]);
      brvec[j] += res2/m_ncounts;
      brerr[j] += res2*res2/m_ncounts;
    }
  }
  toterr = sqrt((toterr - totint*totint)/m_ncounts);
  std::cout << "Full integral = " << totint << " +- " << toterr << " (" << 100.*toterr/totint << "%)" << std::endl;
  std::cout << std::setprecision(2) << std::scientific;
  for(int j=0; j<NRes; j++){
    brerr[j] = sqrt((brerr[j]-brvec[j]*brvec[j])/m_ncounts);
    std::cout << " int = " << brvec[j] << " +- " << brerr[j] << " (" << 100.*brerr[j]/brvec[j] << "%) for " << m_model->Res(j)->Name() << std::endl;
  }

  double BrSum = 0;
  std::cout << std::endl;
  for(int j=0; j<NRes; j++){
    brvec[j] /= totint;
    const double err1 = brerr[j]/totint;
    const double err2 = brvec[j]*toterr/totint;
    brvec[j] *= 100;
    brerr[j] = sqrt(err1*err1+err2*err2)*100.;
    BrSum += brvec[j];
    std::cout << "Br = (" << brvec[j] << " +- " << brerr[j] << ")%, rel err " << 100.*brerr[j]/brvec[j] << "% for " << m_model->Res(j)->Name() << std::endl;
  }
  std::cout << "Branchings sum: " << BrSum << std::endl;
  return totint;
}
