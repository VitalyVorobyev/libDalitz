#include "dalitzmcintegral.h"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> // max

using namespace std;

DalitzMCIntegral::DalitzMCIntegral(const AbsDalitzModel *model):
  RandomDalitzPoint(model),
  m_model(model),
  m_ncounts(1e6),
  m_prefix("")
{
}

double DalitzMCIntegral::GetDalitzPlotArea(const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  cdouble mABsqMin = m_model->mABsq_min();
  cdouble mABsqMax = m_model->mABsq_max();

  cdouble mACsqMin = m_model->mACsq_min();
  cdouble mACsqMax = m_model->mACsq_max();

  U64 nhits = 0;
  double mABsq, mACsq;
  for(U64 count=0; count<Counts; count++){
    GetUnconstrainedPoint(mABsq,mACsq);
    if(IsInPlot(mABsq,mACsq)) nhits++;
  }
  cdouble square = (mABsqMax - mABsqMin)*(mACsqMax - mACsqMin);
  cdouble p = ((double)nhits)/Counts;
  cdouble area  = square*p;
  cdouble darea = square*sqrt(p*(1.-p)/Counts);
  cout << "Dalitz Plot area = " << area << " +- " << darea << " (" << Counts << ") counts" << endl;
  return area;
}

int DalitzMCIntegral::CalcIntegral(double &val, double &err, const U64& nc, const bool norm_flag) const{
  const U64 Counts = nc ? nc : m_ncounts;
  cdouble DPArea  = m_model->Area() == 0 ? GetDalitzPlotArea(1e7) : m_model->Area();
  val = 0; err = 0;
  double mABsq,mACsq;
  double res;
  cout << "DalitzMCIntegral: integrating with " << Counts << " points..." << endl;
  for(U64 i=0; i<Counts; i++){
    if(!(i%100000)) cout << i << " counts val = " << val << " val2 = " << err << endl;
    GetPoint(mABsq,mACsq);
    res = m_model->P(mABsq,mACsq);
    if(norm_flag) res /= m_model->Norm();
    val += res/Counts;
    err += res*res/Counts;
  }
  err = sqrt((err - val*val)/Counts);
  val *= DPArea;
  err *= DPArea;
  cout << "Done! Int = " << val << " +- " << err << " (" << 100.*err/val << "%)" << endl;
  return 0;
}

void DalitzMCIntegral::CheckNormalization(const U64 &nc) const{
  cout << "Checking normalization... " << m_model->Norm() << endl;
  double val, err;
  CalcIntegral(val,err,nc,true);
}

vector<vectcd> DalitzMCIntegral::CalcNormMap(str& ofile,const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  const int ncoef = m_model->AmpNum();
  cdouble DPArea  = m_model->Area() == 0 ? GetDalitzPlotArea(1e7) : m_model->Area();
  // DPArea = 174.2 for B0 -> D0 pi+ pi-

  vector<vectcd> int_matrix;
  vector<vectcd> err_matrix;
  vector<vectd>  err_matrix_real;
  vector<vectd>  err_matrix_imag;
  int_matrix.resize(ncoef);
  err_matrix.resize(ncoef);
  err_matrix_real.resize(ncoef);
  err_matrix_imag.resize(ncoef);
  for(unsigned i=0; i<int_matrix.size(); i++){
    err_matrix[i].resize(ncoef,0);
    int_matrix[i].resize(ncoef,0);
    err_matrix_real[i].resize(ncoef,0);
    err_matrix_imag[i].resize(ncoef,0);
  }
  vectcd res_amp_vec;
  double mABsq;
  double mACsq;
  for(U64 count=0; count<Counts; count++){
    if(count && !(count%100000)) cout << count << " counts" << endl;
    GetPoint(mABsq,mACsq);
    m_model->GetAmpVals(res_amp_vec,mABsq,mACsq);
    for(int i=0; i<ncoef; i++){
      for(int j=i; j<ncoef; j++){
        const compld cvalue = res_amp_vec[i]*conj(res_amp_vec[j]);
        int_matrix[i][j] += cvalue/compld(Counts);
        err_matrix_real[i][j] += pow(real(cvalue),2)/Counts;
        err_matrix_imag[i][j] += pow(imag(cvalue),2)/Counts;
      }
    }
  }
  for(int i=0; i<ncoef; i++){
    for(int j=ncoef-1; j>=i; j--){
      cdouble reerr = err_matrix_real[i][j];
      cdouble imerr = err_matrix_imag[i][j];
      const compld val      = int_matrix[i][j];
      int_matrix[i][j]     *= DPArea;
      err_matrix_real[i][j] = DPArea*sqrt((reerr - pow(real(val),2))/Counts);
      err_matrix_imag[i][j] = DPArea*sqrt((imerr - pow(imag(val),2))/Counts);
      err_matrix[i][j] = compld(err_matrix_real[i][j],err_matrix_imag[i][j]);
      if(i != j){
        int_matrix[j][i]      = int_matrix[i][j];
        err_matrix[j][i]      = err_matrix[i][j];
        err_matrix_real[j][i] = err_matrix_real[i][j];
        err_matrix_imag[j][i] = err_matrix_imag[i][j];
      }
    }
  }
  ShowIntMatrix(int_matrix,err_matrix);
  stringstream out;
  out.str(""); out << "params/NaiveNormMap" << Counts;
  if(m_prefix != str("")) out << "_" << m_prefix;
  out << ".txt";
  WriteIntMatrix(out.str(),int_matrix,err_matrix);
  ofile = out.str();
  return int_matrix;
}


double DalitzMCIntegral::CalcBranchings(vectd& brvec,vectd& brerr, const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  const int NRes = m_model->ResNum();
  brvec.clear(); brerr.clear();
  vectcd amps;
  for(int i=0; i<NRes; i++){ brvec.push_back(0); brerr.push_back(0);}

  double mAB,mAC;
  double totint = 0;
  double toterr = 0;

  const double scale = 1.e5;
  cout << "DalitzMCIntegral: calculating branchings with " << Counts << " points..." << endl;
  for(U64 i=0; i<Counts; i++){
    if(!(i%100000)) cout << i << " counts" << endl;
    GetPoint(mAB,mAC);
    const double res = scale*norm((compld)m_model->GetAmplitudes(amps,mAB,mAC));
    totint += res/Counts;
    toterr += res*res/Counts;
    for(int j=0; j<NRes; j++){
      const double res2 = scale*norm(amps[j]);
      brvec[j] += res2/Counts;
      brerr[j] += res2*res2/Counts;
    }
  }
  toterr = sqrt((toterr - totint*totint)/Counts);
  cout << "Full integral = " << totint << " +- " << toterr << " (" << 100.*toterr/totint << "%)" << endl;
  cout << setprecision(2) << scientific;
  for(int j=0; j<NRes; j++){
    brerr[j] = sqrt((brerr[j]-brvec[j]*brvec[j])/Counts);
    cout << " int = " << brvec[j] << " +- " << brerr[j];
    cout << " (" << (brvec[j]>0 ? 100.*brerr[j]/brvec[j] : 0) << "%) for " << m_model->ResName(j) << endl;
  }

  double BrSum = 0;
  cout << endl;
  for(int j=0; j<NRes; j++){
    brvec[j] /= totint;
    const double err1 = brerr[j]/totint;
    const double err2 = brvec[j]*toterr/totint;
    brvec[j] *= 100;
    brerr[j] = sqrt(err1*err1+err2*err2)*100.;
    BrSum += brvec[j];
    cout << "Br = (" << brvec[j] << " +- " << brerr[j];
    cout << ")%, rel err " << (brvec[j]>0 ? 100.*brerr[j]/brvec[j] : 0);
    cout << "% for " << m_model->ResName(j) << endl;
  }
  cout << "Branchings sum: " << BrSum << endl;
  return totint;
}

vector<vectcd> DalitzMCIntegral::SmartNormMap(const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  const int ncoef  = m_model->AmpNum();

  vector<vectcd> int_matrix;
  int_matrix.resize(ncoef);
  for(unsigned i=0; i<int_matrix.size(); i++) int_matrix[i].resize(ncoef,0);

  for(int i=0; i<ncoef; i++){
    cout << "res1: " << i+1 << endl;
    int_matrix[i][i] = GetGaussIntegral(i,Counts);
    for(int j=i+1; j<ncoef; j++){
      cout << "res2: " << j+1 << endl;
      int_matrix[i][j] = GetGaussIntegral(i,j,Counts);
      int_matrix[j][i] = int_matrix[i][j];
    }
  }

  ShowIntMatrix(int_matrix);
  stringstream out;
  out.str(""); out << "params/SmartNormMap" << Counts << ".txt";
  WriteIntMatrix(out.str(),int_matrix);
  return int_matrix;
}

double normal_pdf(cdouble x, cdouble m, cdouble s){
  cdouble inv_sqrt_2pi = 0.3989422804014327;
  cdouble a = (x-m)/s;
  return inv_sqrt_2pi/s*exp(-0.5*a*a);
}

double DalitzMCIntegral::GetGaussIntegral(const int resnum, const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  DStrip* strip = m_model->GetResAreas()[resnum];
  vectd mABsqV;
  vectd mACsqV;
  vectd mBCsqV;
  if(GetGaussPoints(mABsqV,mACsqV,mBCsqV,Counts,strip)) return 0;
  const vectd& xV = strip->type == 1 ? mABsqV : (strip->type == 2 ? mACsqV : mBCsqV);
  double res = 0;
  cdouble mean = strip->mean();
  cdouble sigm = strip->sigma();
  for(U64 count=0; count<Counts; count++){
    if(count && !(count%100000)) cout << count << " counts" << endl;
    res += norm(m_model->GetResVal(mABsqV[count],mACsqV[count],resnum))/normal_pdf(xV[count],mean,sigm);
  }
  return res/Counts;
}

compld DalitzMCIntegral::GetGaussIntegral(const int resnum1,const int resnum2, const U64& nc) const{
  const U64 Counts = nc ? nc : m_ncounts;
  DStrip* strip1 = m_model->GetResAreas()[resnum1];
  DStrip* strip2 = m_model->GetResAreas()[resnum2];
  vectd mABsqV;
  vectd mACsqV;
  vectd mBCsqV;
  compld res = 0;

  if(strip1->type == strip2->type){
    const int rtype = strip1->type;
    cdouble mean = 0.5*(strip1->mean() + strip2->mean());
    cdouble sigm = max(0.5*fabs(strip1->mean()-strip2->mean()),max(strip1->sigma(),strip2->sigma()));
    DStrip strip(mean-sigm,mean+sigm,rtype);
    if(GetGaussPoints(mABsqV,mACsqV,mBCsqV,Counts,&strip)) return 0;
    const vectd& xV = rtype == 1 ? mABsqV : (rtype == 2 ? mACsqV : mBCsqV);
    for(U64 count=0; count<Counts; count++){
      if(count && !(count%100000)) cout << count << " counts" << endl;
      cdouble gdenc = normal_pdf(xV[count],mean,sigm);
      const compld res1 = m_model->GetResVal(mABsqV[count],mACsqV[count],resnum1);
      const compld res2 = m_model->GetResVal(mABsqV[count],mACsqV[count],resnum2);
      res += res1*conj(res2)/gdenc;
    }
  } else{
    if(GetGaussPoints(mABsqV,mACsqV,mBCsqV,Counts,strip1,strip2)) return 0;
    const vectd& xV = strip1->type == 1 ? mABsqV : (strip1->type == 2 ? mACsqV : mBCsqV);
    const vectd& yV = strip2->type == 1 ? mABsqV : (strip2->type == 2 ? mACsqV : mBCsqV);
    cdouble meanX = strip1->mean();
    cdouble sigmX = strip1->sigma();
    cdouble meanY = strip2->mean();
    cdouble sigmY = strip2->sigma();
    for(U64 count=0; count<Counts; count++){
      if(!(count%100000)) cout << count << " counts" << endl;
      cdouble gdenc = normal_pdf(xV[count],meanX,sigmX) * normal_pdf(yV[count],meanY,sigmY);
      const compld res1 = m_model->GetResVal(mABsqV[count],mACsqV[count],resnum1);
      const compld res2 = m_model->GetResVal(mABsqV[count],mACsqV[count],resnum2);
      res += res1*conj(res2)/gdenc;
    }
  }
  return res/(compld)Counts;
}

int DalitzMCIntegral::WriteIntMatrix(const str& fname,const vector<vectcd>& vals) const {
  ofstream ofile(fname.c_str(),ofstream::out);
  if(!ofile.is_open()){
    cout << "Can't create file " << fname << endl;
    return -1;
  }
  const unsigned size = vals.size();
  ofile << "Size: " << size << endl;
  for(unsigned i=0; i<size; i++){
    for(unsigned j=0; j<size; j++){
      ofile << i << " " << j << " " << vals[i][j] << endl;
    }
  }
  ofile.close();
  return 0;
}

int DalitzMCIntegral::WriteIntMatrix(const str& fname,const vector<vectcd>& vals,const vector<vectcd>& errs) const {
  ofstream ofile(fname.c_str(),ofstream::out);
  if(!ofile.is_open()){
    cout << "Can't create file " << fname << endl;
    return -1;
  }
  const unsigned size = vals.size();
  ofile << "Size: " << size << endl;
  for(unsigned i=0; i<size; i++){
    for(unsigned j=0; j<size; j++){
      ofile << i << " " << j << " " << vals[i][j] << " +- " << errs[i][j] << endl;
    }
  }
  ofile.close();
  return 0;
}

void DalitzMCIntegral::ShowIntMatrix(const vector<vectcd>& vals) const{
  const unsigned size = vals.size();
  cout << "Size: " << size << endl;
  cout << setprecision(2) << scientific;
  for(unsigned i=0; i<size; i++){
    for(unsigned j=0; j<size; j++){
      cout << vals[i][j] << " ";
    }
    cout << endl;
  }
}

void DalitzMCIntegral::ShowIntMatrix(const vector<vectcd>& vals,const vector<vectcd>& errs) const{
  const unsigned size = vals.size();
  cout << "Size: " << size << endl;
  cout << setprecision(2) << scientific;
  for(unsigned i=0; i<size; i++){
    for(unsigned j=0; j<size; j++){
      cout << i << " " << j << " " << vals[i][j] << " +- " << errs[i][j] << endl;
    }
  }
}
