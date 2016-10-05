//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "Core/Particle.hpp"
#include "Core/Event.hpp"

namespace ComPWA {

Event::Event():fWeight(1.),fName(""){

}

Event::Event(const std::string& name):fWeight(1.),fName(name),fFlavour(0),fCharge(0){

}

Event::Event(const double inWeight, const std::string& name="", const double inEff):fWeight(inWeight),fName(name),fFlavour(0),fCharge(0),fEff(inEff){

}

void Event::addParticle(Particle inParticle){
  fParticles.push_back(inParticle);
}

void Event::setParticleAt(const Particle &particle, unsigned int index) {
  fParticles[index] = particle;
}

void Event::reorderEvent(const Event &reference) {
  std::vector<unsigned int> temp;
  for(unsigned int i = 0; i < reference.getNParticles(); ++i) {
    auto iter = fParticles.begin() + i;
    int ref_pid = reference.getParticle(i).pid;
    if(ref_pid != iter->pid) {
      auto result = std::find_if(iter+1, fParticles.end(), [&] (const Particle& p) { return p.pid == ref_pid; });
      if(result == fParticles.end())
        throw std::runtime_error("Event::reorderEvent(): This should not happen...");
      std::iter_swap(iter, result);
    }
  }
}

Event::~Event() { /* nothing */	}

const Particle& Event::getParticle(const unsigned int id) const{
  if(id>=getNParticles()){
    //TODO Exception
    return Particle();
  }
  return fParticles[id];
}

} /* namespace ComPWA */
