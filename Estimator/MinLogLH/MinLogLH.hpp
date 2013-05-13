//! Negative Log Likelihood-Estimator.
/*! \class MinLogLH
 * @file MinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef _MINLOGLH_HPP
#define _MINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Header
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

class MinLogLH : public Estimator {

public:

  virtual double controlParameter(ParameterList& minPar);
  static std::shared_ptr<ControlParameter> createInstance(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);
  static std::shared_ptr<ControlParameter> createInstance(std::shared_ptr<Amplitude>, std::shared_ptr<Data>, std::shared_ptr<Data>);

  /** Destructor */
  virtual ~MinLogLH();

protected:
  /// Default Constructor (0x0)
  MinLogLH(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);
  MinLogLH(std::shared_ptr<Amplitude>, std::shared_ptr<Data>, std::shared_ptr<Data>);

private:
  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<Data> pDIF_;
  std::shared_ptr<Data> pPHSP_;

};

#endif /* _MINLOGLH_HPP */