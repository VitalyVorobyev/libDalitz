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
    DalitzModel(double mm, double ma, double mb, double mc);
    /** Get resonance complex amplitude */
    auto ResAmp(uint32_t n, double mAB, double mAC) const;
    /** Check if the Dalitz plot point is vetoed */
    bool IsVetoed(double mAB, double mAC) const;
    /** To add DalitzResonance (or DalitzPlotObject) to a DalitzModel */
    void AddRes(DalitzPlotObject* res) {m_res_v.push_back(res);}
    /** Get pointer to i'th DalitzPlotObject of DalitzModel */
    const DalitzPlotObject* Res(uint32_t resn) {
        return resn < m_res_v.size() ? m_res_v[resn] : m_res_v[0];
    }
    ///
    double Norm(void) const override {return 1;}
    ///
    void SetParams(const std::vector<double>& pars) override {}
    /// Get vector of complex amplitudes for precalculation of the
    /// normalization integrals
    void GetResVals(std::vector<std::complex<double>>* resv,
                    double mABsq, double mACsq) const override;
    /** Get */
    std::complex<double> GetResVal(double mABsq, double mACsq, int resnum) const override;
    /** Add DalitzVeto to a DalitzModel */
    void AddVeto(DalitzVeto* veto) {m_veto_v.push_back(veto);}
    /** Set amplitude modulus of the i'th DalitzPlotObject of DalitzModel */
    void SetAmp(uint32_t resn, double a) {
        if (resn < m_res_v.size()) m_res_v[resn]->SetAmp(a);
    }
    /** Set amplitude phase of the i'th DalitzPlotObject of DalitzModel */
    void SetPhase(uint32_t resn, double a) {
        if (resn < m_res_v.size()) m_res_v[resn]->SetPhase(a);
    }
    /** Set the list of resonances to be involved in the amplitude calculation */
    bool SetRVec(const std::vector<uint32_t>& rlist);
    /** Calculate Monte Carlo normalization integral */
    int Norm(double* val, double* err, const uint64_t nc = 0);

    int UpdateNormMap(void);
    auto NormElement(int i, int j) const;
    /** Number of Dalitz plot objects in the model */
    auto NRes(void) const {return m_res_v.size();}

    void InitNames();

 private:
    /// Chooses the correct combination of invariant masses squared and
    /// calculates resonance amplitude
    auto GetResAmp(const DalitzPlotObject *res, double mAB, double mAC) const;
    std::vector<DalitzPlotObject*> m_res_v;
    std::vector<DalitzVeto*> m_veto_v;
    int use_subset;
    std::vector<int> m_rlist;
    std::vector<DStrip*> m_res_areas;
    std::vector<std::vector<double>> m_norm_map;
};

#endif  // SRC_DALITZMODEL_H_
