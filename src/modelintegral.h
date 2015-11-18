#ifndef MODELINTEGRAL_H
#define MODELINTEGRAL_H

#include "symdalitzmodel.h"
#include <fstream>

class ModelIntegral{
public:
  ModelIntegral(SymDalitzModel* model);

  void SetGridSize(const int gsize) {m_gsize = gsize; return;}
  void SetNBins(const int nbins) {m_model->SetNBins(nbins); m_nbins = nbins; return;}
  double Calculate(const std::string& label, std::vector<double>& C, std::vector<double>& S, std::vector<double>& K, std::vector<double>& Kb);

private:
  void SetMajorant(void);
  SymDalitzModel* m_model;
  int m_gsize;
  int m_nbins;
  double m_majorant;
};

#endif // MODELINTEGRAL_H
