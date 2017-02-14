/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzresonance.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DALITZRESONANCE_H_
#define SRC_DALITZRESONANCE_H_

#include <string>
#include <complex>

#include "./dalitzplotobject.h"
#include "./resdecayangulardistribution.h"
#include "./formfactor.h"
#include "./abspropagator.h"
#include "./dalitzphasespace.h"

///
/// \brief The DalitzResonance class. Class describes complex amplitude of a
/// three-body decay through an intermediate resonance (M -> RC, R -> AB)
///
class DalitzResonance : public DalitzPlotObject {
 public:
    /// Constuctor for RBW and GS resonances
    DalitzResonance(const std::string& name, const int PropType,
                    const int WidthType, const DalitzPhaseSpace* phsp,
                    const int respath, const double& mres, const double& wres,
                    const int spin, const std::complex<double>& camp);
    /// Constuctor for RBW and GS resonances
    DalitzResonance(const std::string& name, const int PropType,
                    const int WidthType, const DalitzPhaseSpace* phsp,
                    const int respath, const double& mres, const double& wres,
                    const int spin, const double& amp, const double& phi);
    /// Constuctor for RBW and GS resonances
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& mres, const double& wres, const int spin,
                    const std::complex<double>& camp);
    /// Constuctor for RBW and GS resonances
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& mres, const double& wres, const int spin,
                    const double& amp, const double& phi);
    /// Constuctor for NR and Flatte
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& alpha, const std::complex<double>& camp);
    /// Constuctor for NR and Flatte
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& alpha, const double& amp, const double& phi);
    /// Constuctor for Bugg f0(500) and virtual D*0(2010) by Belle
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const std::complex<double>& camp);
    /// Constuctor for Bugg f0(500) and virtual D*0(2010) by Belle
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& amp, const double& phi);
    /// Constructor for virtual D*(2010) by LHCb and rho-omega interference
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& beta1, const double& beta2,
                    const std::complex<double>& amp);
    /// Constructor for virtual D*(2010) by LHCb and rho-omega interference
    DalitzResonance(const std::string& name, const int PropType,
                    const DalitzPhaseSpace* phsp, const int respath,
                    const double& beta1, const double& beta2,
                    const double& amp, const double& phi);

    std::complex<double> evaluate(const double& mACsq,
                                  const double& mBCsq) const;

    int Path(void) const {return m_path;}

    /// Set if the resonance amplitude was changed
    void SetAmpUpd(const bool x) {m_amp_upd = x;}
    /// Set if the resonance parameters were changed
    void SetParUpd(const bool x) {m_par_upd = x;}
    /// Check if the resonance amplitude was changed
    bool IsAmpUpd(void) const {return m_amp_upd;}
    /// Check if the resonance parameters were changed
    bool IsParUpd(void) const {return m_par_upd;}

    static const int AB = 0;
    static const int AC = 1;
    static const int BC = 2;

 private:
    /**
     * @brief SetFFAngAmp Initializes form factors and angular distribution
     * @param phsp
     * @param mres
     * @param vdst
     * @return
     */
    double SetFFAngAmp(const DalitzPhaseSpace* phsp, const double &mres,
                       const bool vdst = false);

    ResDecayAngularDistribution* m_ang_amp;
    FormFactor* m_mff;
    FormFactor* m_rff;
    AbsPropagator* m_prop;

    int m_ptype;
    int m_wtype;
    int m_spin;
    int m_path;

    bool m_amp_upd;
    bool m_par_upd;
};

#endif  // SRC_DALITZRESONANCE_H_
