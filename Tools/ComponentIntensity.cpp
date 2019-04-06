#include "ComponentIntensity.hpp"
#include <boost/algorithm/string.hpp>

namespace ComPWA {
namespace Tools {

// I hope expert systems will not produce the incoherent_with_strength,
// (starts with incoherent)
// incoherent_with_strength
// incoherent
// coherent_0
// <<Moth, daugh1, daugh2>>, <(L_min, L_max)> <(S_min, S_max)>
// source intensity is IncoherentIntensity, not StrengthIncoherentIntensity
std::shared_ptr<ComPWA::Intensity> getComponentIntensityFromIncoherentIntensity(
    const std::shared_ptr<ComPWA::Physics::IncoherentIntensity> 
        incoherentIntensity, const std::string &componentName,
    const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange) {

  std::vector<std::shared_ptr<ComPWA::Intensity>> comps;
  const auto & intensities = incoherentIntensity->getIntensities();
  int coherentIndex = 0;
  for (const auto intensityPtr : intensities) {
    const auto coherentPtr = std::dynamic_pointer_cast<ComPWA::Physics
        ::CoherentIntensity>(intensityPtr);
    const auto amps = getAmplitudesFromCoherentIntensity(
        coherentPtr, componentName, decList, lRange, sRange);
    comps.push_back(std::make_shared<ComPWA::Physics::CoherentIntensity>(
        componentName + "_" + std::to_string(coherentIndex), amps));
    ++coherentIndex;
  }
  return std::make_shared<ComPWA::Physics::IncoherentIntensity>(componentName,
      comps);
}

const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>>
    getAmplitudesFromCoherentIntensity(
    const std::shared_ptr<ComPWA::Physics::CoherentIntensity> coherentIntensity,
    const std::string &componentName,
    const std::vector<std::vector<std::string>> & decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange) {

  std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>> componentAmps;
  const std::vector<std::shared_ptr<ComPWA::Physics::NamedAmplitude>>
      & allAmps = coherentIntensity->getAmplitudes();
  std::vector<std::string> decayNames;
  for (const auto & ampPtr : allAmps) {
    auto ampName = ampPtr->getName();
    decayNames.clear();
    boost::algorithm::split(decayNames, ampName, boost::is_any_of(";"));
    bool isTarget = findAllDecayPatern(decList, lRange, sRange, decayNames);
    if (!isTarget) continue;
    componentAmps.push_back(ampPtr); 
  }
  return componentAmps;
}

bool findAllDecayPatern(const std::vector<std::vector<std::string>> &decList,
    const std::vector<std::pair<int, int>> &lRange,
    const std::vector<std::pair<int, int>> &sRange,
    const std::vector<std::string> &decayNames) {

  // EpEm_-1_to_D*(2007)0_-1+D*(2007)0bar_-1_L_1.0_S_2.0;
  // EpEm_-1_to_D*(2007)0_0+D*(2007)0bar_-1_L_1.0_S_1.0;
  // EpEm_-1_to_D*(2007)0_1+D*(2007)0bar_1_L_3.0_S_2.0;D*(2007)0bar_1
  // _to_D0bar_0+pi0_0_L_1.0_S_0.0;D*(2007)0_1_to_D0_0+pi0_0_L_1.0_S_0.0;

  bool isTarget(true);

  for (std::size_t index = 0; index != decList.size(); ++index) {
    bool found = findDecayPatern(decList.at(index), lRange.at(index),
        sRange.at(index), decayNames);
    if (!found) {
      isTarget = false;
      break;
    }
  }
  return isTarget;
}

bool findDecayPatern(const std::vector<std::string> &decay,
    const std::pair<int, int> &LRange, const std::pair<int, int> &SRange,
    const std::vector<std::string> &decayNames) {

  //EpEm_-1_to_D*(2007)0_-1+D*(2007)0bar_-1_L_1.0_S_2.0;
  //Mother_?_to_daughter1_?+daughter2_?_L_lx_S_sx;
  //Mother_?_to_daughter2_?+daughter1_?_L_lx_S_sx;
  std::string mother = decay.at(0);
  std::string daughter1 = decay.at(1);
  std::string daughter2 = decay.at(2);

  int Lmin = LRange.first;
  int Lmax = LRange.second;
  int Smin = SRange.first;
  int Smax = SRange.second;

  bool found(false);
  std::vector<std::string> splitName;

  for (const auto & name : decayNames) {
    splitName.clear();
    boost::algorithm::split(splitName, name, boost::is_any_of("_+"));

    if (!mother.empty() && mother != splitName.at(0)) continue;
    if (!daughter1.empty() && daughter1 != splitName.at(3)
        && daughter1 != splitName.at(5)) continue;
    if (!daughter2.empty() && daughter2 != splitName.at(3)
        && daughter2 != splitName.at(5)) continue;
  
    int L = std::stoi(splitName.at(8));
    int S = std::stoi(splitName.at(10));

    if (Lmin >= 0 && L < Lmin) continue;
    if (Lmax >= 0 && L > Lmax) continue;
    if (Smin >= 0 && S < Smin) continue;
    if (Smax >= 0 && S > Smax) continue;

    found = true;
    break;
  }
  return found;
}

} // end namespace Tools
} //end namespace ComPWA

