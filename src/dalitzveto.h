/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzveto.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DALITZVETO_H_
#define SRC_DALITZVETO_H_

///
/// \brief The DalitzVeto class. Abstrack class for implementation
/// of vetoes on Dalitz phase space
///
class DalitzVeto {
 public:
    DalitzVeto();
    virtual bool operator()(const double& AB, const double& BC) const = 0;
};

#endif  // SRC_DALITZVETO_H_
