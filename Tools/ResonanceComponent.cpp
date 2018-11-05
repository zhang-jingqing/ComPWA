#include "ResonanceComponent.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/CoherentIntensity.hpp"

namespace ComPWA {
namespace Tools {

std::shared_ptr<AmpIntensity> ResonanceComponent(
    std::shared_ptr<ComPWA::Physics::IncoherentIntensity> incoIntensity,
    std::string name, std::string resName, std::string daug1Name,
    std::string daug2Name, int L, int S) {
  std::string lStr = "L_" + std::to_string(L) + ".0";
  if (L < 0) lStr = ""; 
  std::string sStr = "S_" + std::to_string(S) + ".0";
  if (S < 0) sStr = ""; 
  std::string lsStr = lStr + "_" + sStr;
  if (name == "") {
    name = resName + "_to_" + daug1Name + "+" + daug2Name;
    if (L >= 0) name += "_" + lStr;
    if (S >= 0) name += "_" + sStr;
  }
 
  auto icIn = std::make_shared<ComPWA::Physics::IncoherentIntensity>(*incoIntensity);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes
 
  bool found(false);
  for (auto ampIntensity : incoIntensity->intensities()) {
    auto coherentIntensity = 
        std::dynamic_pointer_cast<ComPWA::Physics::CoherentIntensity>(ampIntensity);
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try {
      comp = Tools::ResonanceComponent(coherentIntensity, name, resName, 
          daug1Name, daug2Name, L, S); 
    } catch (std::exception &ex) {
      continue;
    }   
    found = true;
    icIn->addIntensity(comp);
  }
 
  // Nothing found
  if (!found) {
    throw std::runtime_error(
    "Tools::ResonanceComponent(): | Error component " + name +
    " could not be found in IncoherentIntensity " + incoIntensity->name() + ".");
  }
 
  return icIn;
}
std::shared_ptr<AmpIntensity> ResonanceComponent(
    std::shared_ptr<ComPWA::Physics::CoherentIntensity> coIntensity,
    std::string name, std::string resName, std::string daug1Name,
    std::string daug2Name, int L, int S) {
  std::string lStr = "L_" + std::to_string(L) + ".0";
  if (L < 0) lStr = ""; 
  std::string sStr = "S_" + std::to_string(S) + ".0";
  if (S < 0) sStr = ""; 
  std::string lsStr = lStr + "_" + sStr;
  if (name == "") {
    name = resName + "_to_" + daug1Name + "+" + daug2Name;
    if (L >= 0) name += "_" + lStr;
    if (S >= 0) name += "_" + sStr;
  }

  auto icIn = std::make_shared<ComPWA::Physics::CoherentIntensity>(*coIntensity);
  icIn->setName(name);
  icIn->reset(); // delete all existing amplitudes
  bool found(false);
  for (auto seqAmp : coIntensity->amplitudes()) {
    std::string seqAmpName = seqAmp->name();
    bool foundComponent(false);
    std::size_t begPos = 0;
    std::size_t endPos = 0;
    while (begPos != std::string::npos) {
      endPos = seqAmpName.find(";", begPos);
      //endPos == std::string::npos in substr means until end of the string
      std::string subName = seqAmpName.substr(begPos, endPos);
      begPos = endPos;
      if (std::string::npos == subName.find(resName))
        continue;
      if (daug1Name != "" && std::string::npos == subName.find(daug1Name))
        continue;
      if (daug2Name != "" && std::string::npos == subName.find(daug2Name))
        continue;
      if (L >= 0 && std::string::npos == subName.find(lStr))
        continue; 
      if (S >= 0 && std::string::npos == subName.find(sStr))
        continue;
      foundComponent = true;
      break;
    }   
    if (!foundComponent)
      continue;
    found = true;
    icIn->addAmplitude(seqAmp);
  }

  // Nothing found, this could be true 
  // and it is not a error when this is called
  // in the IncoherentIntnesity::component()
  if (!found) {
    throw std::runtime_error(
        "Tools::ResonanceComponent() | Error component " + name +
        " could not be found in CoherentIntensity " + coIntensity->name() + ".");
  }

  return icIn;
}

}
}
