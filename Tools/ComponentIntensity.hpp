// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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

//can not return a reference of a vector
std::shared_ptr<ComPWA::Intensity> getComponentIntensityFromIncoherentIntensity(
    const std::shared_ptr<ComPWA::Physics::IncoherentIntensity> 
        incoherentIntensity, const std::string &componentName,
    const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange);

//this could become a method of CoherentIntensity
const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> 
    getAmplitudesFromCoherentIntensity(
    const std::shared_ptr<ComPWA::Physics::CoherentIntensity> coherentIntensity,
    const std::string &componentName,
    const std::vector<std::vector<std::string>> & decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange);

bool findAllDecayPatern(const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange,
    const std::vector<std::string> &decayNames);

bool findDecayPatern(const std::vector<std::string> &decay,
    const std::pair<int, int> &LRange, const std::pair<int, int> &SRange,
    const std::vector<std::string> &decayNames);

} // end namespace Tools
} //end namespace ComPWA

#endif
