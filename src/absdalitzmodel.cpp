#include "absdalitzmodel.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
typedef AbsDalitzModel AbsDM;

AbsDalitzModel::AbsDalitzModel(const DalitzPhaseSpace& phsp):
  DalitzPhaseSpace(phsp)
{
}

AbsDalitzModel::AbsDalitzModel(cdouble& mM, cdouble& mA, cdouble& mB, cdouble& mC):
  DalitzPhaseSpace(mM,mA,mB,mC)
{
}

compld AbsDM::Amp(cdouble& mABsq, cdouble& mACsq) const{
  vectcd resvals;
  compld amp = 0;
  GetAmpVals(resvals,mABsq,mACsq);
  for(unsigned i=0; i<m_ampl.size(); i++)
    amp += m_ampl[i]*resvals[i];
  return amp;
}

vector<MnPar> AbsDM::MnPars(void) const{
  cdouble phmin = -2.0*M_PI;
  cdouble phmax =  2.0*M_PI;
  cdouble pherr =  0.1*M_PI;
  vector<MnPar> parv;
  stringstream out;
  for(unsigned i=0; i<m_amp_names.size(); i++){
    out.str(""); out << "amp_" << m_amp_names[i];
    cdouble amp = abs(m_ampl[i]);
    parv.push_back(MnPar(out.str(),amp,0.1*amp,amp*0.1,amp*10.));
    out.str(""); out << "arg_" << m_amp_names[i];
    cdouble pha = arg(m_ampl[i]);
    parv.push_back(MnPar(out.str(),pha,pherr,phmin,phmax));
  }
  return parv;
}

double AbsDM::P(  cdouble& mABsq, cdouble& mACsq) const{ return norm(Amp(mABsq,mACsq));}
double AbsDM::Arg(cdouble& mABsq, cdouble& mACsq) const{ return  arg(Amp(mABsq,mACsq));}

compld AbsDM::GetAmplitudes(vectcd& resv, cdouble& mABsq, cdouble& mACsq) const{
  resv.clear();
  vectcd resvals;
  GetResVals(resvals,mABsq,mACsq);
  if(m_ampl.size() == resvals.size()){
    for(unsigned i=0; i<m_amp_signature.size(); i++) resv.push_back(m_ampl[i]*resvals[i]);
  } else{
    int res = 0;
    for(unsigned i=0; i<m_amp_signature.size(); i++){
       for(int j=0; j<m_amp_signature[i]; j++){
         resv.push_back(m_ampl[i]*resvals[res++]);
       }
    }
  }
  compld sum(0);
  for(compld amp : resv) sum += amp;
  return sum;
}

void AbsDM::GetAmpVals(vectcd& resv, cdouble& mABsq,cdouble& mACsq) const{
  resv.clear();
  vectcd resvals;
  GetResVals(resvals,mABsq,mACsq);
  if(m_ampl.size() == resvals.size()){
    resv = resvals;
    return;
  } else{
    int resnum = 0;
    for(unsigned i=0; i<m_amp_signature.size(); i++){
       resv.push_back(0);
       for(int j=0; j<m_amp_signature[i]; j++){
         resv[i] += resvals[resnum++];
       }
    }
  }
}

int AbsDM::OpenCachedIntegrals(const std::string& fname,const bool silent){
  ifstream ifile(fname.c_str(),ifstream::in);
  if(!ifile.is_open()){
    cout << "Can't open file " << fname << endl;
    return -1;
  }
  if(!silent) cout << "Reading cashed normalization from file " << fname << endl;
  string line;
  int size = 0;
  getline(ifile,line);
  sscanf(line.c_str(),"Size: %d",&size);
  if(!silent) cout << size << " amplitudes" << endl;
  m_res_int.resize(size);
  for(int i=0; i<size; i++) m_res_int[i].resize(size,0);
  int num1, num2;
  double val_re, val_im, err_re, err_im;
  while(getline(ifile,line)){
    sscanf(line.c_str(),"%d %d (%lf,%lf) +- (%lf,%lf)",&num1,&num2,&val_re,&val_im,&err_re,&err_im);
    m_res_int[num1][num2] = compld(val_re,val_im);
  }
  if(!silent){
    cout << "Success!" << endl;
    for(int i=0; i<size; i++){
      for(int j=0; j<size; j++){
        cout << i << " " << j << " " << m_res_int[i][j] << endl;
      }
    }
  }
  return 0;
}

str AbsDM::GetAmpStr(void) const{
  stringstream out;
  out.str("");
  out << setprecision(2) << fixed;
  for(compld amp : m_ampl) out << abs(amp) << "*exp(i*" << arg(amp) << ") ";
  return out.str();
}

void AbsDM::ShowAmpls(void) const{
  cout << GetAmpStr() << endl;
}

double AbsDalitzModel::NormWithCache(void) const{
  if(!m_ampl.size() || m_ampl.size() != m_res_int.size()){
    cout << "Can't calculate normalizations with cached integrals: " << m_ampl.size() << ", " << m_res_int.size() << endl;
    return 1;
  }
  double res = 0;
  for(unsigned i=0; i<m_ampl.size(); i++){
    // I += |a_i|^2 * I_i
    res += norm(m_ampl[i])*real(m_res_int[i][i]);
    for(unsigned j=i+1; j<m_ampl.size(); j++){
      // I += 2*Re( a_i * a_j^* * I_{ij} )
      res += 2.*real(m_ampl[i]*conj(m_ampl[j])*m_res_int[i][j]);
    }
  }
  return res;
}

void AbsDalitzModel::GetCoefficients(vectcd& coefv) const{coefv = m_ampl;}

void AbsDalitzModel::SetResAreas(const vectd& ledge,const vectd& redge, const vecti& types){
  m_res_areas.clear();
  for(unsigned i=0; i<types.size(); i++) m_res_areas.push_back(new DStrip(ledge[i],redge[i],types[i]));
}
