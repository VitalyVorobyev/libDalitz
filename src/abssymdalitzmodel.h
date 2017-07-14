/** Copyright 2017 Vitaly Vorobyev
 ** @file abssymdalitzmodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_ABSSYMDALITZMODEL_H_
#define SRC_ABSSYMDALITZMODEL_H_

#include <string>
#include <cmath>

#include "./absdalitzmodel.h"

class AbsSymDalitzModel : virtual public AbsDalitzModel {
 public:
    /**
     * @brief AbsSymDalitzModel. Abstract class for decay amplitude models of
     * three-body decays M -> ABC, where m(B) = m(C).
     * Some interface for binned analysis provided.
     * @param mmo
     * @param mcha
     * @param mchb
     * @param delmin
     * @param delmax
     */
    AbsSymDalitzModel(const double& mm, const double& ma, const double& mb,
                      const double& delmin = 0,
                      const double& delmax = 2.*M_PI);

    double delta(const double& mp, const double& mm) const;
    void ppdelt(const double& mp, const double& mm,
                double* const pp, double* const pn, double* const del) const;
    /**
     * @brief bin. Dalitz plot bin for a given Dalitz variables
     * @param mp
     * @param mm
     * @return
     */
    int bin(const double& mp, const double& mm) const;
    /**
     * @brief bin. Signed Dalitz plot bin for a given phase difference
     * @param mp
     * @param mm
     * @param dphi
     * @return
     */
    int bin(const double& mp, const double& mm, const double& dphi) const;
    /**
     * @brief bin. Dalitz plot bin for a given phase difference
     * @param dphi
     * @return
     */
    int bin(const double& dphi, const bool sign) const;

    void nbins(const unsigned nb) {m_nbins = nb;}
    unsigned nbins(void) const {return m_nbins;}

    /**
     * @brief Tabulate. Calculate amplitude on a mAB x mAC half grid
     * and save the result in text file
     * @param fname. File name
     * @param grid_size. Grid size
     */
    void TabulateSymABAC(const std::string& fname,
                         const unsigned grid_size = 1000) const;

    /**
     * @brief Tabulate. Calculate amplitude on a mAB x mBC grid
     * and save the result in text file
     * @param fname. File name
     * @param grid_size. Grid size
     */
    void TabulateSymABBC(const std::string& fname,
                         const unsigned grid_size = 1000) const;

 private:
    double del_min, del_max;
    unsigned m_nbins;

    static const double m_2pi;
    static int sign(const double& x);
    double delta(const std::complex<double>& ap,
                 const std::complex<double>& an) const;
};

#endif  // SRC_ABSSYMDALITZMODEL_H_
