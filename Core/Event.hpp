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
//! Internal container for event information.
/*! \class Event
 * @file Event.hpp
 * This class provides a internal container for event-based information. The
 * class provides a list of particles of the event.
*/

#ifndef _Event_HPP_
#define _Event_HPP_

#include <vector>
#include <string>
#include "Core/Particle.hpp"

class Event{

public:
  Event();

  Event(const std::string& name);

  Event(const double inWeight, const std::string& name);

  virtual void addParticle(Particle inParticle);

  virtual ~Event();

  virtual void inline setName(const std::string& name) { fName = name; }
  virtual const inline std::string& getName() { return fName; }
  virtual double inline getWeight() {return fWeight;};
  virtual void inline setWeight(double w) { fWeight=w;};
  virtual int inline getFlavour() {return fFlavour;};
  virtual void inline setFlavour(int fl) { fFlavour = fl;};
  virtual int inline getCharge() {return fCharge;};
  virtual void inline setCharge(int ch) { fCharge = ch;};

  virtual const inline unsigned int getNParticles() { return fParticles.size(); }
  virtual const Particle& getParticle(const unsigned int id);

protected:
  std::vector<Particle> fParticles;
  double fWeight;
  std::string fName;
  int fFlavour; //1 -> particle, 0 -> unknown, -1 anti-particle
  int fCharge;
  //Particle fParticleB;
  //TODO: other event info?

};

#endif