/** Copyright 2017 Vitaly Vorobyev
 ** @file absdalitzmodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_ABSDALITZMODEL_H_
#define SRC_ABSDALITZMODEL_H_

#include <complex>
#include <string>
#include <vector>

#include "./dalitzphasespace.h"
#include "./dstrip.h"
#include "./mnpar.h"

/**
 * @brief The AbsDalitzModel class
 */
class AbsDalitzModel : public DalitzPhaseSpace {
 public:
    explicit AbsDalitzModel(const DalitzPhaseSpace& phsp);
    AbsDalitzModel(const double& mM, const double& mA, const double& mB,
                   const double& mC);

    /// Get DalitzModel complex amplitude
    std::complex<double> Amp(const double& mABsq, const double& mACsq) const;
    /// Get modulo squared DalitzModel amplitude
    double P(const double& mABsq, const double& mACsq) const;
    /// Get complex phase (radians) of DalitzModel amplitude
    double Arg(const double& mABsq, const double& mACsq) const;

    /// Set caption for mAB^2 axis
    void SetABaxis(const std::string& s) {mABaxis = s;}
    /// Set caption for mAC^2 axis
    void SetACaxis(const std::string& s) {mACaxis = s;}
    /// Set caption for mBC^2 axis
    void SetBCaxis(const std::string& s) {mBCaxis = s;}
    /**
     * @brief SetModelTitle. Set model Title
     * @param s. The title string
     */
    void SetModelTitle(const std::string& s) {m_title = s;}
    /**
     * @brief SetParState
     * @param v
     */
    void SetParState(const std::vector<bool>& v) {m_parstate = v;}

    /// Get caption for mAB^2 axis
    std::string ABaxis(void) const { return mABaxis;}
    /// Get caption for mAC^2 axis
    std::string ACaxis(void) const { return mACaxis;}
    /// Get caption for mBC^2 axis
    std::string BCaxis(void) const { return mBCaxis;}

    /// Get name of the i'th DalitzPlotObject
    std::string ResName(const int resn) const {return m_res_names[resn];}
    /// Get name of the i'th DalitzPlotObject
    std::string AmpName(const int resn) const {return m_amp_names[resn];}

    /// To fill std::vector<std::complex<double>> resv with amplitudes of
    /// each DalitzPlotObject of a DalitzModel. Returns size full amplitude
    std::complex<double> GetAmplitudes(std::vector<std::complex<double>>* resv,
                                       const double& mABsq,
                                       const double& mACsq) const;
    /// Dimention of the normalization matrix
    unsigned AmpNum(void) const { return m_ampl.size();}
    ///
    const std::vector<bool>& PState(void) const {return m_parstate;}

    ///
    int OpenCachedIntegrals(const std::string& fname,
                            const bool silent = true);

    /**
     * @brief NormWithCache. Speed up the computation of normalization using
     * the relation
     * I = \sum_i |a_i|^2 I_i + 2Re(\sum_{i>j} a_i a_j^* I_{ij}), where
     * I_i = \int |A_i|^2 dm_+^2 dm_-^2 and
     * I_{ij} = \int A_i A_j^* dm_+^2 dm_-^2,
     * where A_i is a complex amplitude of i^th resonance
     * @return Value of normalization integral
     */
    double NormWithCache(void) const;

    /**
     * @brief Tabulate. Calculate amplitude on a mAB x mAC grid
     * and save the result in text file
     * @param fname. File name
     * @param grid_size. Grid size
     */
    void TabulateABAC(const std::string& fname,
                      const unsigned grid_size = 1000) const;

    /**
     * @brief Tabulate. Calculate amplitude on a mAB x mBC grid
     * and save the result in text file
     * @param fname. File name
     * @param grid_size. Grid size
     */
    void TabulateABBC(const std::string& fname,
                      const unsigned grid_size = 1000) const;

    /**
     * @brief GetCoefficients. Get vector of coefs
     * @param coefv. Vector to assign
     */
    void GetCoefficients(std::vector<std::complex<double>>* coefv) const;

    ///
    std::vector<DStrip*> GetResAreas(void) const {return m_res_areas;}

    ///
    const DStrip* GetResArea(const int i) const {return m_res_areas[i];}

    /// Get number of DalitzPlotObject in DalitzModel
    unsigned ResNum(void) const {return m_res_names.size();}

    /// Get vector of complex amplitudes for precalculation of
    /// the normalization integrals
    void GetAmpVals(std::vector<std::complex<double>>* resv,
                    const double& mABsq,
                    const double& mACsq) const;

    /**
     * @brief GetAmpStr. Generates text with values of complex amplitudes
     * @return std::string
     */
    std::string GetAmpStr(void) const;

    /**
     * @brief AsText. String with amplitude info.
     * @return std::string
     */
    std::string AsText(void) const;

    /// Get vector of complex amplitudes for all resonances
    virtual void GetResVals(std::vector<std::complex<double>>* resv,
                            const double& mABsq,
                            const double& mACsq) const = 0;
    /// Get
    virtual std::complex<double> GetResVal(const double& mABsq,
                                           const double& mACsq,
                                           const int resnum) const = 0;
    /// Set parameters
    virtual void SetParams(const std::vector<double>& pars) = 0;
    /// Get normalization
    virtual double Norm(void) const = 0;
    /// Vector of parameter for minuit2
    std::vector<MnPar> MnPars(void) const;

 protected:
    std::string mABaxis;
    std::string mACaxis;
    std::string mBCaxis;

    void SetResNames(const std::vector<std::string>& vnames) {
        m_res_names = vnames;
    }
    void SetAmpNames(const std::vector<std::string>& vnames) {
        m_amp_names = vnames;
    }
    void SetCoefficients(const std::vector<std::complex<double>>& coefv) {
        m_ampl = coefv;
    }
    void SetAmpSignature(const std::vector<unsigned>& vals) {
        m_amp_signature = vals;
    }
    void SetResAreas(const std::vector<DStrip*>& vec) {m_res_areas = vec;}
    void SetResAreas(const std::vector<double> &ledge,
                     const std::vector<double> &redge,
                     const std::vector<int>& types);

 private:
    /**
     * @brief m_title. Decay model title
     */
    std::string m_title;
    /// Matrix of cached normalization integrals
    std::vector<std::vector<std::complex<double>>> m_res_int;
    /// List of complex coefficients (Amplitudes) for normalization units
    std::vector<std::complex<double>> m_ampl;
    /// List of bools. 1 - parameter is fixed, 0 - parameter is released
    std::vector<bool> m_parstate;
    /// Correspondence between Amplitudes and Resonances
    std::vector<unsigned> m_amp_signature;
    /// List of resonance names
    std::vector<std::string> m_res_names;
    /// List of normalization unit names
    std::vector<std::string> m_amp_names;
    /// List of Dalitz std::stringips for normalization units
    std::vector<DStrip*> m_res_areas;
};

#endif  // SRC_ABSDALITZMODEL_H_
