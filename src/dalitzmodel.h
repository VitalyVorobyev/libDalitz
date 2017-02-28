/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzmodel.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_DALITZMODEL_H_
#define SRC_DALITZMODEL_H_

#include <vector>
#include <string>
#include <complex>
#include <cstdint>

#include "./absdalitzmodel.h"
#include "./dalitzplotobject.h"
#include "./dalitzresonance.h"
#include "./dalitzveto.h"

/**
 * \brief Class for full description of a three-body decay model
 * A DalitzModel class object contains std::vector of objects inheritanced
 * from DalitzPlotObject class. Method Amp(const double& mAB, const
 * double& mAC) returns coherent sum of complex amplitudes of
 * DalitzPlotObject's.
 */
class DalitzModel : virtual public AbsDalitzModel {
 public:
    DalitzModel(const double& mm, const double& ma,
                const double& mb, const double& mc);
    /// Get DalitzModel complex amplitude
    std::complex<double> Amp(const double& mAB, const double& mAC) const;
    /// Get resonance complex amplitude
    std::complex<double> ResAmp(const unsigned n, const double& mAB,
                                const double& mAC) const;
    /// Check if the Dalitz plot point is vetoed
    bool IsVetoed(const double& mAB, const double& mAC) const;
    /// To fill std::vector<std::complex<double>> vec with amplitudes of each
    /// DalitzPlotObject of a DalitzModel. Returns size full amplitude
    std::complex<double> GetAmplitudes(std::vector<std::complex<double>>* vec,
                                       const double& mAB,
                                       const double& mAC) const;
    /// To add DalitzResonance (or DalitzPlotObject) to a DalitzModel
    void AddRes(DalitzPlotObject* res) {m_res_v.push_back(res);}

    /// Get pointer to i'th DalitzPlotObject of DalitzModel
    const DalitzPlotObject* Res(const unsigned resn) {
        return resn < m_res_v.size() ? m_res_v[resn] : m_res_v[0];
    }
    ///
    double Norm(void) const {return 1;}
    ///
    void SetParams(const std::vector<double>& pars) {}
    /// Get vector of complex amplitudes for precalculation of the
    /// normalization integrals
    void GetResVals(std::vector<std::complex<double>>* resv,
                    const double& mABsq, const double& mACsq) const;
    /// Get
    std::complex<double> GetResVal(const double& mABsq, const double& mACsq,
                                   const int resnum) const {
        return ResAmp(resnum, mABsq, mACsq);
    }
    /// To add DalitzVeto to a DalitzModel
    void AddVeto(DalitzVeto* veto) {m_veto_v.push_back(veto);}
    /// Set amplitude modulus of the i'th DalitzPlotObject of DalitzModel
    void SetAmp(const unsigned resn, const double& a) {
        if (resn < m_res_v.size()) m_res_v[resn]->SetAmp(a);
    }
    /// Set amplitude phase of the i'th DalitzPlotObject of DalitzModel
    void SetPhase(const unsigned resn, const double& a) {
        if (resn < m_res_v.size()) m_res_v[resn]->SetPhase(a);
    }
    /// Set the list of resonances to be involved in the amplitude calculation
    bool SetRVec(const std::vector<unsigned>& rlist);
    /// Calculate Monte Carlo normalization integral
    int Norm(double* val, double* err, const uint64_t nc = 0);

    int UpdateNormMap(void);
    double NormElement(const int i, const int j) const;

 private:
    /// Chooses the correct combination of invariant masses squared and
    /// calculates resonance amplitude
    std::complex<double> GetResAmp(const DalitzPlotObject *res,
                                   const double& mAB,
                                   const double& mAC) const;
    std::vector<DalitzPlotObject*> m_res_v;
    std::vector<DalitzVeto*> m_veto_v;
    int use_subset;
    std::vector<int> m_rlist;
    std::vector<DStrip*> m_res_areas;

    std::vector<std::vector<double>> m_norm_map;
};

#endif  // SRC_DALITZMODEL_H_
