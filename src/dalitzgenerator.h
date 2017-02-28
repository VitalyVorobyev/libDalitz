/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzgenerator.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#ifndef SRC_DALITZGENERATOR_H_
#define SRC_DALITZGENERATOR_H_

#include <vector>
#include <string>
#include <cstdint>

#include "./absdalitzmodel.h"
#include "./randomdalitzpoint.h"

///
/// \brief The DalitzGenerator class for generation of Dalitz
/// plot distribution. Can generate events for any amplitude
/// defined with AbsDalitzModel descendant class.
///
class DalitzGenerator : public RandomDalitzPoint {
 public:
    ///
    /// \brief DalitzGenerator
    /// \param model
    ///
    explicit DalitzGenerator(const AbsDalitzModel* model);
    /**
     * @brief WriteDDist. Generate events and write it in text file
     * @param NEv. Number of events to generate
     * @param fname. Output text file name.
     * @return 0 if succeed
     */
    int WriteDDist(const uint64_t& NEv, const std::string& fname) const;
    ///
    /// \brief Generate method generates Dalitz distribution via
    /// Neumann method
    /// \param NEv --- number of events to be generated
    /// \param mABv --- vector of m^2(AB) to be filled
    /// \param mACv --- vector of m^2(AC) to be filled
    /// \param silent --- flag for suppression of terminal output.
    /// Suppression is done if flag is true (false by default).
    /// \return Returns 0 if worked properly
    ///
    int Generate(const uint64_t NEv, std::vector<double>* mABv,
                 std::vector<double>* mACv, const bool silent = false) const;
    ///
    /// \brief Generate method generates a single Dalitz plot point
    /// according to the PDF specified by model
    /// \param mABsq m^2(AB) value to be written in this variable
    /// \param mACsq m^2(AC) value to be written in this variable
    /// \return Returns 0 if worked properly
    ///
    int Generate(double* mABsq, double* mACsq,
                 std::uniform_real_distribution<double> *dist = nullptr) const;
    ///
    /// \brief SetMaxTries
    /// \param p
    ///
    void SetMaxTries(const uint64_t& p) {m_ntries = p;}
    ///
    /// \brief GetMaxTries
    /// \return
    ///
    uint64_t GetMaxTries(void) const {return m_ntries;}
    ///
    /// \brief SetNMajCounts
    /// \param p
    ///
    void SetNMajCounts(const unsigned p) {m_maj_counts = p;}
    ///
    /// \brief SetMajorant
    /// \param p
    ///
    void SetMajorant(const double& p);

 private:
    ///
    /// \brief CalcMajorant
    /// \return
    ///
    double CalcMajorant(void) const;
    ///
    /// \brief m_ntries
    ///
    uint64_t m_ntries;
    unsigned m_maj_counts;
    double m_maj;
    const AbsDalitzModel* m_model;

//    std::default_random_engine* ren;
};

#endif  // SRC_DALITZGENERATOR_H_
