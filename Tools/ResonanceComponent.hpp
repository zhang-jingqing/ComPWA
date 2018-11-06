#ifndef RESONANCE_COMPONENT_HPP
#define RESONANCE_COMPONENT_HPP

#include "Core/AmpIntensity.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"

namespace ComPWA{
namespace Tools{
  std::shared_ptr<AmpIntensity> ResonanceComponent(
      std::shared_ptr<ComPWA::Physics::IncoherentIntensity> incoIntensity,
      std::string name, std::string resName, std::string daug1Name = "", 
      std::string daug2Name = "", int L = -1, int S = -1);
  std::shared_ptr<AmpIntensity> ResonanceComponent(
      std::shared_ptr<ComPWA::Physics::CoherentIntensity> coIntensity,
      std::string name, std::string resName, std::string daug1Name = "", 
      std::string daug2Name = "", int L = -1, int S = -1);
}
}

#endif
