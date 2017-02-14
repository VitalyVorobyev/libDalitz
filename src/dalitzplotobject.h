/** Copyright 2017 Vitaly Vorobyev
 ** @file dalitzplotobject.h
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **/

#ifndef SRC_DALITZPLOTOBJECT_H_
#define SRC_DALITZPLOTOBJECT_H_

#include <string>
#include <complex>

/**
 * @brief The DalitzPlotObject class. Abstract class for any object can
 * appears on a Dalitz diagram. The class contains name, complex amplitude
 * of the object and virtual method evaluate()
 */
class DalitzPlotObject {
 public:
    DalitzPlotObject(const std::string& name, const std::complex<double>& amp =
            std::complex<double>(0, 0));
    DalitzPlotObject(const std::string& name,
                     const double& a, const double& phi);
    virtual std::complex<double> evaluate(const double& mACsq,
                                          const double& mBCsq) const = 0;
    virtual int Path(void) const = 0;

    /// Modificators
    void SetName(const std::string& name) {m_name = name; return;}
    void SetCAmp(const std::complex<double>& amp);
    void SetAmp(const double& a);
    void SetPhase(const double& phi);

    /// Interface
    std::string Name(void) const {return m_name;}
    std::complex<double> CAmp(void) const {return m_amp;}
    double Amp(void) const {return m_a;}
    double Phase(void) const {return m_phi;}

 private:
    std::string m_name;
    std::complex<double> m_amp;
    double m_a;
    double m_phi;
};

#endif  // SRC_DALITZPLOTOBJECT_H_
