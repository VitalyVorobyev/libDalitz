#ifndef DALITZGENERATOR_H
#define DALITZGENERATOR_H

#include <vector>

#include "absdalitzmodel.h"
#include "randomdalitzpoint.h"

typedef std::uniform_real_distribution<double> unifdist;
typedef std::default_random_engine rndmeng;
typedef unsigned long long U64;

///
/// \brief The DalitzGenerator class for generation of Dalitz plot distribution.
/// Can generate events for any amplitude defined with AbsDalitzModel descendant class.
///
class DalitzGenerator : public RandomDalitzPoint{
public:
  ///
  /// \brief DalitzGenerator
  /// \param model
  ///
  DalitzGenerator(const AbsDalitzModel* model);
  ///
  /// \brief Generate method generates Dalitz distribution via Neumann method
  /// \param NEv --- number of events to be generated
  /// \param mABv --- vector of m^2(AB) to be filled
  /// \param mACv --- vector of m^2(AC) to be filled
  /// \param silent --- flag for suppression of terminal output. Suppression is
  /// done if flag is true (false by default).
  /// \return Returns 0 if worked properly
  ///
  int Generate(const int NEv, vectd& mABv, vectd& mACv, const bool silent = false) const;
  ///
  /// \brief Generate method generates a single Dalitz plot point according to the
  /// PDF specified by model
  /// \param mABsq m^2(AB) value to be written in this variable
  /// \param mACsq m^2(AC) value to be written in this variable
  /// \return Returns 0 if worked properly
  ///
  int Generate(double& mABsq, double& mACsq) const;
  ///
  /// \brief SetMaxTries
  /// \param p
  ///
  void SetMaxTries(const U64& p) {m_ntries = p; return;}
  ///
  /// \brief GetMaxTries
  /// \return
  ///
  U64 GetMaxTries(void) const {return m_ntries;}
  ///
  /// \brief SetNMajCounts
  /// \param p
  ///
  void SetNMajCounts(const int p){m_maj_counts = p; return;}
  ///
  /// \brief SetMajorant
  /// \param p
  ///
  void SetMajorant(cdouble& p);

private:
  ///
  /// \brief CalcMajorant
  /// \return
  ///
  double CalcMajorant(void) const;
  ///
  /// \brief m_ntries
  ///
  U64    m_ntries;
  int    m_maj_counts;
  double m_maj;
  const AbsDalitzModel* m_model;

//  unifdist* unifMaj;
  rndmeng*  ren;
};

#endif // DALITZGENERATOR_H
