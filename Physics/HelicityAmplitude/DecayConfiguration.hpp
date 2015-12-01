//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_

#include <map>
#include <vector>

#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"

namespace HelicityFormalism {

struct DecayProductsInfo {
  std::vector<unsigned int> particle_indices_;
  boost::property_tree::ptree decay_strength_info_and_phase_;
};

typedef std::map<unsigned int, std::vector<DecayProductsInfo> > ParticleIndexDecayTree;


class DecayConfiguration {
  friend class DecayTreeFactory;

  std::vector<ParticleStateInfo> particles_;

  std::vector<ParticleIndexDecayTree> concrete_decay_trees_;

  ParticleIndexDecayTree current_concrete_decay_tree_;

public:
  DecayConfiguration();
  virtual ~DecayConfiguration();

  void addCurrentDecayTreeToList();

  void addDecayToCurrentDecayTree(const ParticleStateInfo& mother,
      const std::vector<ParticleStateInfo>& daughter_states,
      const boost::property_tree::ptree& decay_strength_info_and_phase);

  std::vector<unsigned int> addParticlesToList(
      const std::vector<ParticleStateInfo>& particle_list);

  unsigned int addParticleToList(ParticleStateInfo particle);

  void setRemainingParticleProperties(ParticleStateInfo& particle) const;
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYCONFIGURATION_HPP_ */
