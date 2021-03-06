// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PHYSICS_NORMALIZATIONAMPLITUDEDECORATOR_HPP_
#define PHYSICS_NORMALIZATIONAMPLITUDEDECORATOR_HPP_

#include "Physics/Amplitude.hpp"

namespace ComPWA {
class Intensity;
namespace Tools {
class IntegrationStrategy;
}

namespace Physics {

class NormalizationAmplitudeDecorator : public NamedAmplitude {
public:
  NormalizationAmplitudeDecorator(
      const std::string &name, std::shared_ptr<NamedAmplitude> amplitude,
      std::shared_ptr<ComPWA::Tools::IntegrationStrategy> integrator);

  std::complex<double> evaluate(const ComPWA::DataPoint &point) const final;

  void updateParametersFrom(const ParameterList &list) final;
  void addUniqueParametersTo(ParameterList &list) final;
  void addFitParametersTo(std::vector<double> &FitParameters) final;

  std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const final;

  std::shared_ptr<const Amplitude> getUnnormalizedAmplitude() const;

private:
  bool checkParametersChanged() const;

  std::shared_ptr<Amplitude> UnnormalizedAmplitude;
  std::shared_ptr<Intensity> NormedAmplitude;

  double Normalization;
  std::vector<double> PreviousFitParameters;
  /// Phsp sample for numerical integration
  std::shared_ptr<ComPWA::Tools::IntegrationStrategy> Integrator;
};

} // namespace Physics
} // namespace ComPWA

#endif
