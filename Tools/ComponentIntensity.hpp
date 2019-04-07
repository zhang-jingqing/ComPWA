// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//===----------------------------------------------------------------------===//
////
//// \file
//// This file contains the declaration of functions used to get component
///  intensity which satisfies some requirements from the given intensity.
////
//===----------------------------------------------------------------------===//
#ifndef COMPWA_TOOLS_COMPONENT_INTENSITY_HPP_
#define COMPWA_TOOLS_COMPONENT_INTENSITY_HPP_

#include <memory>
#include <vector>
#include <string>
#include <utility>

#include "Physics/Amplitude.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"

namespace ComPWA {
namespace Tools {

/// Get compoent amplitudes from a coherent intensity. Returned intensity
/// satisfy the requirements on decay/products particles, orbital angular 
/// momentum and spin of helicity decays.
///
/// \param incoherentIntensity The global incohernt intensity.
/// \param componentName The name of the returned component intensity .
///
/// \param decList Names of mother and daughter particles of helicity decays.
//  Each vector<string> contains three elements: {mother, daughter1, daughter2}.
/// \param lRange Orbital angular momentum of two daughters of helicty decays.
/// Each pair<int, int> means a L range \f$(L_{min}, L_{max})\f$.
/// \param sRange Spin of two daughters of helicity decays.
/// Each pair<int, int> means a S range \f$(S_{min}, S_{max})\f$.
/// \param decayNames Splited amplitude name.
std::shared_ptr<ComPWA::Intensity> getComponentIntensity(
    const std::shared_ptr<ComPWA::Physics::IncoherentIntensity> 
        incoherentIntensity, const std::string &componentName,
    const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange);

/// Get compoent amplitudes from a coherent intensity. Returned amplitudes
/// satisfy the requirements on decay/products particles, orbital angular 
/// momentum and spin of helicity decays.
//
/// \param coherentIntensity The coherent intensity.
/// \param decList Names of mother and daughter particles of helicity decays.
//  Each vector<string> contains three elements: {mother, daughter1, daughter2}.
/// \param lRange Orbital angular momentum of two daughters of helicty decays.
/// Each pair<int, int> means a L range \f$(L_{min}, L_{max})\f$.
/// \param sRange Spin of two daughters of helicity decays.
/// Each pair<int, int> means a S range \f$(S_{min}, S_{max})\f$.
/// \param decayNames Splited amplitude name.
const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> 
    getAmplitudes(
    const std::shared_ptr<ComPWA::Physics::CoherentIntensity> coherentIntensity,
    const std::vector<std::vector<std::string>> & decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange);

/// Find if all decay requirements(decay and product particles' name, orbital
/// angular momentum ans spin of helicity decay nodes) are satisfied in the
/// splited amplitude name.
//
/// \param decList Names of mother and daughter particles of helicity decays.
//  Each vector<string> contains three elements: {mother, daughter1, daughter2}.
/// \param lRange Orbital angular momentum of two daughters of helicty decays.
/// Each pair<int, int> means a L range \f$(L_{min}, L_{max})\f$.
/// \param sRange Spin of two daughters of helicity decays.
/// Each pair<int, int> means a S range \f$(S_{min}, S_{max})\f$.
/// \param decayNames Splited amplitude name.
bool findAllDecayPatern(const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange,
    const std::vector<std::string> &decayNames);

/// Find if the decay requirements(decay and product particles' name, orbital
/// angular momentum ans spin of helicity decay nodes) are satisfied in the
/// splited amplitude name.
//
/// \param decay Names of mother and daughter particles of a helicity decay.
//  Elements are in the order mother, daughter1, daughter2.
/// \param lRange Orbital angular momentum of two daughters of a helicty decay.
/// The pair<int, int> means a L range \f$(L_{min}, L_{max})\f$.
/// \param sRange Spin of two daughters of a helicity decay.
/// The pair<int, int> means a S range \f$(S_{min}, S_{max})\f$.
/// \param decayNames Splited amplitude name.
bool findDecayPatern(const std::vector<std::string> &decay,
    const std::pair<int, int> &LRange, const std::pair<int, int> &SRange,
    const std::vector<std::string> &decayNames);

} // end namespace Tools
} //end namespace ComPWA

#endif
