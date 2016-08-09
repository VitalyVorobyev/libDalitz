#ifndef DALITZVETO_H
#define DALITZVETO_H

///
/// \brief The DalitzVeto class. Abstrack class for implementation of vetoes on Dalitz phase space
///
class DalitzVeto{
public:
  DalitzVeto();
  virtual bool operator()(const double& AB, const double& BC) const = 0;
};

#endif // DALITZVETO_H
