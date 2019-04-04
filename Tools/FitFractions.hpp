// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_FITFRACTIONS_HPP_
#define COMPWA_TOOLS_FITFRACTIONS_HPP_

#include <memory>
#include <vector>

#include "Core/ProgressBar.hpp"
#include "Data/DataSet.hpp"
#include "Physics/Amplitude.hpp"
#include "Physics/CoherentIntensity.hpp"
#include "Tools/Integration.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>

namespace ComPWA {
namespace Tools {

/// Print gsl_matrix
inline void gsl_matrix_print(const gsl_matrix *m) {
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      printf("%g ", gsl_matrix_get(m, i, j));
    }
    printf("\n");
  }
};

/// Convert ParameterList to gsl vector
inline gsl_vector *gsl_parameterList2Vec(const ParameterList &list) {
  size_t n = list.doubleParameters().size();
  gsl_vector *tmp = gsl_vector_alloc(n);
  size_t t = 0;
  for (size_t o = 0; o < n; o++) {
    auto outPar = list.doubleParameter(o);
    if (outPar->isFixed())
      continue;
    gsl_vector_set(tmp, t, outPar->value());
    t++;
  }
  // resize vector
  gsl_vector *vec = gsl_vector_alloc(t);
  for (unsigned int i = 0; i < vec->size; i++)
    gsl_vector_set(vec, i, gsl_vector_get(tmp, i));
  return vec;
};

/// Convert std::vector matrix to gsl matrix
inline gsl_matrix *
gsl_vecVec2Matrix(const std::vector<std::vector<double>> &m) {
  gsl_matrix *tmp = gsl_matrix_alloc(m.size(), m.at(0).size());
  for (size_t i = 0; i < tmp->size1; i++) {
    for (unsigned int j = 0; j < tmp->size2; j++) {
      gsl_matrix_set(tmp, i, j, m.at(i).at(j));
    }
  }
  return tmp;
};

/// Multivariate Gaussian using cholesky decomposition
/// A test application can be found at test/MultiVariateGaussianTestApp.cpp
///
/// \param rnd Random number generator
/// \param vecSize Size of data vector
/// \param in Mean value(s)
/// \param cov Covariance matrix
/// \param res Resulting Vector
inline void multivariateGaussian(const gsl_rng *rnd, const unsigned int vecSize,
                                 const gsl_vector *in, const gsl_matrix *cov,
                                 gsl_vector *res) {
  // Generate and fill temporary covariance matrix
  gsl_matrix *tmpM = gsl_matrix_alloc(vecSize, vecSize);
  gsl_matrix_memcpy(tmpM, cov);

  // Cholesky decomposition
  int status = gsl_linalg_cholesky_decomp(tmpM);
  if (status == GSL_EDOM)
    LOG(ERROR) << "Decomposition has failed!";

  // Compute vector of random gaussian variables
  for (unsigned int i = 0; i < vecSize; i++)
    gsl_vector_set(res, i, gsl_ran_ugaussian(rnd));

  // Debug
  //  gsl_matrix_print(cov);
  //  gsl_vector_print(in);
  //  gsl_vector_print(res);
  //  gsl_matrix_print(tmpM);

  // From the GNU GSL Documentation:
  // The function dtrmv compute the matrix-vector product x = op(A) x for the
  // triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans,
  // CblasTrans, CblasConjTrans. When Uplo is CblasUpper then the upper
  // triangle of A is used, and when Uplo is CblasLower then the lower
  // triangle of A is used. If Diag is CblasNonUnit then the diagonal of
  // the matrix is used, but if Diag is CblasUnit then the diagonal elements
  // of the matrix A are taken as unity and are not referenced.
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tmpM, res);

  gsl_vector_add(res, in);
  // free temporary object
  gsl_matrix_free(tmpM);
};

/// Calculates the fit fractions using the formula:
/// \f[
///  f_i = \frac{|c_i|^2 \int A_i A_i^*}{\int \sum c_l c_m^* A_l A_m^*}
/// \f]
/// The \f$c_i\f$ are the complex coefficient of the amplitudes \f$A_i\f$ and
/// the denominator is the integral over the whole intensity.
/// The integrals are performed via Monte Carlo integration method.
ComPWA::ParameterList calculateFitFractions(
    std::shared_ptr<const ComPWA::Physics::CoherentIntensity> intensity,
    std::shared_ptr<ComPWA::Data::DataSet> sample,
    const std::vector<std::string> &components = {}) {
  LOG(DEBUG) << "calculating fit fractions...";
  ComPWA::ParameterList FitFractionsList;

  // calculate denominator
  double IntegralDenominator = ComPWA::Tools::integrate(intensity, sample);

  // TODO: ultimately we want to have all combinations of amplitudes
  // A_i x A_j*
  // so loop over the amplitudes and in the second loop start with the same
  // index
  // -in that way we will only calculate one half of the matrix, plus the
  // diagonal
  // -then we need a complexconjugate amplitudedecorator, which we can
  // use here to perform the operation!

  // create list of numerators
  std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> Amps;
  for (auto const &AmpName : components) {
    for (auto x : intensity->getAmplitudes()) {
      if (x->getName() == AmpName) {
        Amps.push_back(x);
      }
    }
  }
  if (0 == Amps.size()) {
    for (auto x : intensity->getAmplitudes()) {
      Amps.push_back(x);
    }
  }

  // calculate nominators
  for (auto x : Amps) {
    auto ci = std::make_shared<ComPWA::Physics::CoherentIntensity>(
        "TempIntensity",
        std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>>{x});
    double IntegralNumerator = ComPWA::Tools::integrate(ci, sample);

    double fitfraction = IntegralNumerator / IntegralDenominator;
    LOG(TRACE) << "calculateFitFractions(): fit fraction for (" << x->getName()
               << ") is " << fitfraction;

    FitFractionsList.addParameter(
        std::make_shared<ComPWA::FitParameter>(x->getName(), fitfraction, 0.0));
  }
  LOG(DEBUG) << "finished fit fraction calculation!";
  return FitFractionsList;
}

/// Calculates the fit fractions with errors.
/// The errors are calculated by generating sets of fit parameters that are
/// smeared by a multidimensional gaussian using the covariance matrix of the
/// fit. For every set the fit fractions are calculated. From this distribution
/// the standard errors give the errors of the fit fractions. This can be a very
/// time consuming method.
inline ParameterList calculateFitFractionsWithSampledError(
    std::shared_ptr<const ComPWA::Physics::CoherentIntensity> CohIntensity,
    std::shared_ptr<ComPWA::Data::DataSet> Sample,
    const std::vector<std::vector<double>> &covariance, size_t nSets,
    const std::vector<std::string> &Components = {}) {

  ComPWA::ParameterList FitFractions =
      calculateFitFractions(CohIntensity, Sample, Components);

  LOG(INFO) << "calculateFitFractionsWithErrorSampling(): Calculating errors "
               "of fit fractions using "
            << nSets << " sets of parameters...";

  // in principle a copy of the intensity should be made to not interfere with
  // the constness
  auto intens =
      std::const_pointer_cast<ComPWA::Physics::CoherentIntensity>(CohIntensity);

  ParameterList parameters;
  intens->addUniqueParametersTo(parameters);

  if (nSets <= 0)
    return FitFractions;
  if (!parameters.doubleParameters().size())
    return FitFractions;

  size_t NumberOfFreeParameters(0);
  for (auto const &x : parameters.doubleParameters()) {
    if (!x->isFixed())
      ++NumberOfFreeParameters;
  }
  if (NumberOfFreeParameters != covariance.size())
    return FitFractions;

  // Check whether covariance matrix is set to zero.
  bool leave = true;
  for (unsigned int i = 0; i < covariance.size(); ++i)
    for (unsigned int j = 0; j < covariance.size(); ++j)
      if (covariance.at(i).at(j) != 0.)
        leave = false;
  if (leave)
    LOG(ERROR) << "CalcFractionError() | Covariance matrix is zero "
                  "(everywhere)! We skip further "
                  "calculation of fit fraction errors.";

  // Setting up random number generator
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gsl_rng *rnd = gsl_rng_alloc(T);

  // convert to GSL objects
  gsl_vector *gslFinalPar = gsl_parameterList2Vec(parameters);
  gsl_matrix *gslCov = gsl_vecVec2Matrix(covariance);
  gsl_matrix_print(gslCov); // DEBUG

  int nFreeParameter = covariance.size();

  ParameterList originalPar;
  originalPar.DeepCopy(parameters);

  std::vector<ParameterList> fracVect;
  ProgressBar bar(nSets);
  unsigned int i = 0;
  while (i < nSets) {
    bool error = 0;
    gsl_vector *gslNewPar = gsl_vector_alloc(nFreeParameter);
    // generate set of smeared parameters
    multivariateGaussian(rnd, nFreeParameter, gslFinalPar, gslCov, gslNewPar);
    // gsl_vector_print(gslNewPar); // Debug

    // deep copy of finalParameters
    ParameterList newPar;
    newPar.DeepCopy(parameters);

    std::size_t t = 0;
    for (std::size_t o = 0; o < newPar.doubleParameters().size(); o++) {
      std::shared_ptr<FitParameter> outPar = newPar.doubleParameter(o);
      if (outPar->isFixed())
        continue;
      // set floating values to smeared values
      try { // catch out-of-bound
        outPar->setValue(gslNewPar->data[t]);
      } catch (ParameterOutOfBound &ex) {
        error = 1;
      }
      t++;
    }
    // free vector
    gsl_vector_free(gslNewPar);
    if (error)
      continue; // skip this set if one parameter is out of bound

    // update amplitude with smeared parameters
    try {
      intens->updateParametersFrom(newPar);
    } catch (ParameterOutOfBound &ex) {
      continue;
    }
    fracVect.push_back(calculateFitFractions(intens, Sample, Components));
    bar.next();
    i++;
  }

  // free objects
  gsl_vector_free(gslFinalPar);
  gsl_matrix_free(gslCov);
  gsl_rng_free(rnd);

  // Calculate standard deviation
  for (unsigned int o = 0; o < FitFractions.doubleParameters().size(); ++o) {
    double mean = 0, sqSum = 0., stdev = 0;
    for (unsigned int i = 0; i < fracVect.size(); ++i) {
      double tmp = fracVect.at(i).doubleParameter(o)->value();
      mean += tmp;
      sqSum += tmp * tmp;
    }
    unsigned int s = fracVect.size();
    sqSum /= s;
    mean /= s;
    // Equivalent to RMS of the distribution
    stdev = std::sqrt(sqSum - mean * mean);
    FitFractions.doubleParameter(o)->setError(stdev);
  }

  intens->updateParametersFrom(originalPar);
  return FitFractions;
}

/// Calculates the fit fractions with errors via error propagation from the
/// covariance matrix. The gradients are calculated via numerical
/// differentiation:
/// \f[
/// f´(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
/// \f]
ComPWA::ParameterList calculateFitFractionsWithCovarianceErrorPropagation(
    std::shared_ptr<const ComPWA::Physics::CoherentIntensity> CohIntensity,
    std::shared_ptr<ComPWA::Data::DataSet> Sample,
    const std::vector<std::vector<double>> &CovarianceMatrix,
    const std::vector<std::string> &Components = {}) {

  ComPWA::ParameterList FitFractions =
      calculateFitFractions(CohIntensity, Sample, Components);

  // in principle a copy of the intensity should be made to not interfere with
  // the constness
  auto intens =
      std::const_pointer_cast<ComPWA::Physics::CoherentIntensity>(CohIntensity);

  ParameterList TempParameters;
  intens->addUniqueParametersTo(TempParameters);
  ParameterList PreviousParameters;
  PreviousParameters.DeepCopy(TempParameters);

  std::vector<ParameterList> Gradients;

  for (std::shared_ptr<FitParameter> FitPar :
       TempParameters.doubleParameters()) {
    if (FitPar->isFixed())
      continue;
    double TempValue = FitPar->value();

    double h = std::sqrt(std::numeric_limits<double>::epsilon()) * TempValue;

    FitPar->setValue(TempValue + h);
    intens->updateParametersFrom(TempParameters);
    ParameterList ff_ph = calculateFitFractions(intens, Sample, Components);

    FitPar->setValue(TempValue - h);
    intens->updateParametersFrom(TempParameters);
    ParameterList ff_mh = calculateFitFractions(intens, Sample, Components);

    size_t NumberOfFitFractions(ff_ph.doubleParameters().size());
    ParameterList FitFractionGradients;
    for (size_t i = 0; i < NumberOfFitFractions; ++i) {
      FitFractionGradients.addParameter(std::make_shared<ComPWA::FitParameter>(
          ff_ph.doubleParameter(i)->name(),
          (ff_ph.doubleParameter(i)->value() -
           ff_mh.doubleParameter(i)->value()) /
              (2. * h)));
    }
    Gradients.push_back(FitFractionGradients);
    FitPar->setValue(TempValue);
  }

  std::vector<double> FitFractionErrors(FitFractions.doubleParameters().size());

  for (unsigned int ffi = 0; ffi < FitFractionErrors.size(); ++ffi) {
    for (unsigned int par1 = 0; par1 < Gradients.size(); ++par1) {
      for (unsigned int par2 = 0; par2 < Gradients.size(); ++par2) {
        FitFractionErrors[ffi] +=
            CovarianceMatrix[par1][par2] *
            Gradients[par1].doubleParameter(ffi)->value() *
            Gradients[par2].doubleParameter(ffi)->value();
      }
    }
  }

  for (unsigned int ffi = 0; ffi < FitFractionErrors.size(); ++ffi) {
    FitFractions.doubleParameter(ffi)->setError(
        std::sqrt(FitFractionErrors[ffi]));
  }

  intens->updateParametersFrom(PreviousParameters);

  return FitFractions;
}

/// Check if intensity and components use same parameter set
/// i.e., parameter with same (unique) name, has same value
/// every parameter in components can be found in intensity
inline bool sharedParameters(
    const std::shared_ptr<const ComPWA::Intensity> intensity,
    const std::vector<std::shared_ptr<const ComPWA::Intensity>> &components) {

  ComPWA::ParameterList parameters;
  auto tempIntensity = std::const_pointer_cast<ComPWA::Intensity>(intensity);
  tempIntensity->addUniqueParametersTo(parameters);
  ComPWA::ParameterList componentParameters;
  for (std::size_t icomp = 0; icomp < components.size(); ++icomp) {
    auto tempComp = 
        std::const_pointer_cast<ComPWA::Intensity>(components.at(icomp));
    tempComp->addUniqueParametersTo(componentParameters);
  }

  auto const intensityFitPars = parameters.doubleParameters();
  auto const componentFitPars = componentParameters.doubleParameters();
  if (intensityFitPars.size() <= componentFitPars.size()) {
    LOG(INFO) << "ComPWA::Tools::sharedParameters(): "
        "global intensity and intensity of components have"
        " unshared parameters!";
    return false;
  }
  for (auto const & componentPar : componentFitPars) {
    std::string targetName = componentPar->name();
    double targetValue = componentPar->value();
    std::string name;
    double value;
    bool shared(false);
    for (auto const & intensityPar : intensityFitPars) {
      name = intensityPar->name();
      value = intensityPar->value();
      if (name == targetName) {
        if (abs(value - targetValue) < 1e-6) {
          shared = true;
          break;
        }
      } 
    }
    if (shared) continue;
    LOG(INFO) << "ComPWA::Tools::sharedParameters(): global intensity and "
        "intensity of components have unshared parameter: " << name << " = "
        << value << " vs " << targetName << " = " << targetValue << "!";
    return false;
  }
  return true;
}

inline ComPWA::ParameterList calculateFitFractions2(
    const std::shared_ptr<const ComPWA::Intensity> intensity,
    const std::shared_ptr<ComPWA::Data::DataSet> sample,
    const std::vector<std::string> & componentNames,
    const std::vector<std::shared_ptr<const ComPWA::Intensity>> &components) {
  LOG(DEBUG) << "calculating fit fractions...";
  ComPWA::ParameterList FitFractionsList;

  ////make sure components and intensity use same parameters
  if (!sharedParameters(intensity, components)) {
    LOG(INFO) << "calculateFitFractions2(): Found unshared parameters! Do not "
        "calculate fit fractions!";
    return FitFractionsList;
  }
  
  // calculate denominator
  double IntegralDenominator = ComPWA::Tools::integrate(intensity, sample);
  //calculate numerators and fit fractions
  for (std::size_t icomp = 0; icomp < componentNames.size(); ++icomp) {
    double IntegralNumerator = 
        ComPWA::Tools::integrate(components.at(icomp), sample);
    double FitFraction = IntegralNumerator / IntegralDenominator;
    LOG(TRACE) << "calculateFitFractions(): fit fraction for ("
        << componentNames.at(icomp) << ") is " << FitFraction;

    FitFractionsList.addParameter(
        std::make_shared<ComPWA::FitParameter>(componentNames.at(icomp),
        FitFraction, 0.0));
  }
  LOG(DEBUG) << "finished fit fraction calculation!";
  return FitFractionsList;
}

inline bool smearParameters(const ComPWA::ParameterList &parameters,
    const std::vector<std::vector<double>> &covariance, size_t nSets,
    std::vector<ComPWA::ParameterList> &vecNewPars) {
  vecNewPars.clear();

  if (nSets <= 0) return false;
  vecNewPars.reserve(nSets);

  if (!parameters.doubleParameters().size()) return false;

  std::size_t nFreeParameters(0);
  for (auto const &para : parameters.doubleParameters()) {
    if (!para->isFixed()) {
      ++nFreeParameters;
    }
  }
  if (nFreeParameters != covariance.size()) return false;

  // Check whether covariance matrix is set to zero
  bool isAllZero(true);
  for (std::size_t i = 0; i < covariance.size(); ++i) {
    for (std::size_t j = 0; j < covariance.size(); ++j) {
      if (covariance.at(i).at(j) != 0.) {
        isAllZero = false; 
        break;
      }
    }
  }
  if (isAllZero) {
    LOG(ERROR) << "smearParameters() | Covariance matrix is zero everywhere! "
        "We skip further smear calculation.";
    return false;
  }

  // Setting up random number generator
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gsl_rng *rnd = gsl_rng_alloc(T);

  // convert to GSL objects
  gsl_vector *gslFinalPar = gsl_parameterList2Vec(parameters);
  gsl_matrix *gslCov = gsl_vecVec2Matrix(covariance);
  gsl_matrix_print(gslCov); // DEBUG

  ProgressBar bar(nSets);
  std::size_t iSet = 0;
  while (iSet < nSets) {
    bool error = 0;
    gsl_vector *gslNewPar = gsl_vector_alloc(nFreeParameters);
    // generate set of smeared parameters
    multivariateGaussian(rnd, nFreeParameters, gslFinalPar, gslCov, gslNewPar);
    // gsl_vector_print(gslNewPar); // Debug

    // deep copy of finalParameters
    ParameterList newPar;
    newPar.DeepCopy(parameters);
    std::size_t i = 0;
    for (std::size_t j = 0; j < newPar.doubleParameters().size(); ++j) {
      std::shared_ptr<FitParameter> outPar = newPar.doubleParameter(j);
      if (outPar->isFixed()) continue;
      // set floating values to smeared values
      try { // catch out-of-bound
        outPar->setValue(gslNewPar->data[i]);
      } catch (ParameterOutOfBound &ex) {
        error = 1;
        break;
      }
      ++i;
    }
    // free vector
    gsl_vector_free(gslNewPar);
    if (error) continue;
    
    vecNewPars.push_back(newPar);
    bar.next();
    ++iSet;
  }
  return true;
}
    
inline ComPWA::ParameterList calculateFitFractionsWithSampledError2(
    const std::shared_ptr<const ComPWA::Intensity> intensity,
    const std::shared_ptr<ComPWA::Data::DataSet> sample,
    const std::vector<std::string> & componentNames,
    const std::vector<std::shared_ptr<const ComPWA::Intensity>> &components,
    const std::vector<std::vector<double>> &covariance, size_t nSets) {

  if (!sharedParameters(intensity, components)) {
    LOG(INFO) << "calculateFitFractionsWithSampledError2(): Found unshared "
        "parameters! Do not calculate fit fractions and errors!";
  }

  ComPWA::ParameterList calculateFitFractions2(
      const std::shared_ptr<const ComPWA::Intensity> intensity,
      const std::shared_ptr<ComPWA::Data::DataSet> sample,
      const std::vector<std::string> & componentNames,
      const std::vector<std::shared_ptr<const ComPWA::Intensity>> &components);

  ComPWA::ParameterList fitFractions = calculateFitFractions2(
      intensity, sample, componentNames, components);

  LOG(INFO) << "calculateFitFractionsWithSampledError(): Calculationg errors "
      "of fit fractions using " << nSets << " sets of parameters...";

  // copy original parameters, Note, parameters' ptr in tempIntensity 
  // now point to par in parameters(i.e., they depend on each other)
  // we need a deepcopy of parameters to get an independent copy of parameters
  ComPWA::ParameterList parameters;
  auto tempIntensity = std::const_pointer_cast<ComPWA::Intensity>(intensity);
  tempIntensity->addUniqueParametersTo(parameters);
  ParameterList originalParameters;
  originalParameters.DeepCopy(parameters);

  std::vector<ComPWA::ParameterList> componentParameters(
      components.size(), ComPWA::ParameterList());
  std::vector<ComPWA::ParameterList> originalComponentParameters(
      components.size(), ComPWA::ParameterList());
  std::vector<std::shared_ptr<ComPWA::Intensity>> tempComponents;
  tempComponents.reserve(components.size());
  for (std::size_t icomp = 0; icomp < components.size(); ++icomp) {
    auto tempComp =
        std::const_pointer_cast<ComPWA::Intensity>(components.at(icomp));
    tempComponents.push_back(tempComp);
    tempComponents.at(icomp)
        ->addUniqueParametersTo(componentParameters.at(icomp));

    originalComponentParameters.at(icomp)
        .DeepCopy(componentParameters.at(icomp));
  }
  
  std::vector<ComPWA::ParameterList> vecNewPars;

  bool smearParameters(const ComPWA::ParameterList &parameters,
    const std::vector<std::vector<double>> &covariance, size_t nSets,
    std::vector<ComPWA::ParameterList> &vecNewPars);

  bool smearOK = smearParameters(parameters, covariance, nSets, vecNewPars);
  if (!smearOK) return fitFractions;
  
  ProgressBar bar(nSets * componentNames.size());

  std::vector<ComPWA::ParameterList>
      fractionVec(nSets, ComPWA::ParameterList());
  std::vector<std::string> oneComponentNames(1, std::string());
  std::vector<std::shared_ptr<const ComPWA::Intensity>> oneComponent(1,
      std::shared_ptr<ComPWA::Intensity>());

  for (std::size_t iSet = 0; iSet < nSets; ++iSet) {
    tempIntensity->updateParametersFrom(vecNewPars.at(iSet));
    for (std::size_t icomp = 0; icomp < componentNames.size(); ++icomp) {
      tempComponents.at(icomp)->updateParametersFrom(vecNewPars.at(iSet));

      oneComponentNames.at(0) = componentNames.at(icomp);
      oneComponent.at(0) = tempComponents.at(icomp);
       
      ComPWA::ParameterList componentFraction = 
          calculateFitFractions2(tempIntensity, sample, oneComponentNames, 
          oneComponent);
      fractionVec.at(iSet).addParameter(componentFraction.doubleParameter(0));
      bar.next();
    }
  }


  //reset to original parameters
  tempIntensity->updateParametersFrom(originalParameters);
  for (std::size_t icomp = 0; icomp < components.size(); ++icomp) {
    tempComponents.at(icomp)->updateParametersFrom(
        originalComponentParameters.at(icomp));
  }

  // Calcualate standard deviation
  for (std::size_t icomp = 0; icomp < fitFractions.doubleParameters().size();
      ++icomp) {
    double mean = 0, sqSum = 0, stddev = 0;
    for (std::size_t iset = 0; iset < fractionVec.size(); ++iset) {
      double temp = fractionVec.at(iset).doubleParameter(icomp)->value();
      mean += temp;
      sqSum += temp * temp;
    }
    std::size_t nset = fractionVec.size();
    sqSum /= nset;
    mean /= nset;

    // Equivalent to RMS of the distribution
    stddev = std::sqrt(sqSum - mean * mean);
    fitFractions.doubleParameter(icomp)->setError(stddev);
  }

  return fitFractions;
}

} // namespace Tools
} // namespace ComPWA

#endif
