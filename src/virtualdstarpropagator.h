#ifndef VIRTUALDSTARPROPAGATOR_H
#define VIRTUALDSTARPROPAGATOR_H

#include "consts.h"
#include "abspropagator.h"

///
/// \brief The VirtualDstarPropagator class
///
class VirtualDstarPropagator : public AbsPropagator{
public:
  VirtualDstarPropagator(const double& beta1, const double& beta2);
//  ~VirtualDstarPropagator(void);

  compld operator()(const double& s,const double& p = 0) const;
private:
  double m_b1;
  double m_b2;
};

#endif // VIRTUALDSTARPROPAGATOR_H
