/** Copyright 2017 Vitaly Vorobyev
 ** @file blattweisskopf.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_BLATTWEISSKOPF_H_
#define SRC_BLATTWEISSKOPF_H_

#include "./formfactor.h"

/**
 * @brief The BlattWeisskopf class. Blatt-Weisskopf penetration form factor
 * for a resonance R->AB. Taken from CLEO preprint 00-23 (hep-ex/0011065)
 * See original paper J. Blatt and V. Weisskopf, "Theoretical Nuclear Physics"
 * (Wiley, New York, 1952).
 */
class BlattWeisskopf : public FormFactor {
 public:
    BlattWeisskopf(const int spin, const double& p0, const int type);

    double operator()(const double& p) const;
    /**
     * @brief m_r_meson radius of a meson state
     */
    static double m_r_meson;
    /**
     * @brief m_r_resonance. radius of a resonance
     */
    static double m_r_resonance;
    static const int FFMeson;
    static const int FFResonance;

 private:
    /**
     * @brief m_spin Angular momentum of resonance
     */
    int m_spin;
    /**
     * @brief m_type FF type (meson or resonance)
     */
    int m_type;
    /**
     * @brief m_F0 formula evaluated at p0
     */
    double m_F0;

    double compute(const double& psq) const;
};

#endif  // SRC_BLATTWEISSKOPF_H_
