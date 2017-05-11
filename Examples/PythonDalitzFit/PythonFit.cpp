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
//		Peter Weidenkaff - adding functionality to generate set of
//events
//-------------------------------------------------------------------------------

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <ctime>
#include <numeric>
#include <time.h>

// Root header files go here
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TMath.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

//Core header files go here
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Efficiency.hpp"

// ComPWA header files go here
#include "DataReader/Data.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "DataReader/JakeReader/JakeReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"
#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"


#include <boost/filesystem.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>

#include "executables/Examples/PythonDalitzFit/PythonFit.hpp"

using namespace boost::log;

using namespace ComPWA;
using Physics::DPKinematics::DalitzKinematics;
using Physics::AmplitudeSum::AmpSumIntensity;
//using DataReader::RootReader::RootReader;
//using DataReader::JakeReader::JakeReader;
using Estimator::MinLogLH::MinLogLH;

PythonFit::PythonFit() {
  M = 3.096916; // GeV/c² (J/psi+)
  Br = 0.000093; // GeV/c² (width)
  m1 = 0.; // GeV/c² (gamma)
  m2 = 0.139570; // GeV/c² (pi)
  m3 = 0.139570; // GeV/c² (pi)
  //const Double_t c = 299792458.; // m/s
  PI = 3.14159; // m/s

  nFitEvents = 44000 - 1;
  nStartEvent = 0;

  int seed = 1234;

  bool useFctTree = false, resultGen = true;

  const char* pPath = getenv("COMPWA_DIR");
  std::string path = "";
  try {
      path = std::string(pPath);
  } catch (std::logic_error& ex) {
      BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
  }

  string output_filename("PythonFitResult.root");
  string output_directory("./executables/Examples/PythonDalitzFit/");
  string input_filename("3Parts-4vecs_1M_PDG.root");
  string input_directory("./executables/Examples/PythonDalitzFit/");

  // create correct output file name
  boost::filesystem::path input_file(input_filename);
  string input_file_path = input_directory + "/" + input_filename;

  Logging log("log", boost::log::trivial::debug); //initialize logging

  PhysConst::createInstance(path+"/Physics/");
  // initialize kinematics of decay
  DalitzKinematics* kin =
          dynamic_cast<DalitzKinematics*>(DalitzKinematics::createInstance(
                  "jpsi", "gamma", "pi0", "pi0"));

  std::string resoFile = path + "/executables/Examples/PythonDalitzFit/JPSI_ypipi.xml";

  auto fitAmpPtr = new AmpSumIntensity("amp",normStyle::none,
          std::shared_ptr<Efficiency>(new UnitEfficiency()), nFitEvents);
  fitAmpPtr->Configure(resoFile);
  std::shared_ptr<Amplitude> amps( fitAmpPtr );

  BOOST_LOG_TRIVIAL(info)<< "Load Modules";
  std::string file = path + "/executables/Examples/PythonDalitzFit/3Part-4vecs_1M_PDG.root";

  std::shared_ptr<RootReader> myReader(
      new RootReader(file, "data"));
  myReader->resetWeights(); //setting weights to 1
  std::shared_ptr<RootReader> myPHSPReader(
          new RootReader(file, "mc"));
  myPHSPReader->setEfficiency(shared_ptr<Efficiency>(new UnitEfficiency())); //setting efficiency to 1

  //ParameterList par;
  std::shared_ptr<ControlParameter> esti;
  amps->FillParameterList(par); //perfect startvalues
  esti = MinLogLH::createInstance(amps, myReader, myPHSPReader, nStartEvent,
          nFitEvents);
  MinLogLH* contrPar = dynamic_cast<MinLogLH*>(&*(esti->Instance()));
  std::shared_ptr<FunctionTree> tree;
  if (useFctTree) {
      contrPar->setUseFunctionTree(1);
      tree = contrPar->getTree();
  }

  std::shared_ptr<Optimizer::Optimizer> opti(new Optimizer::Minuit2::MinuitIF(esti, par));

}

/*PythonFit::PythonFit(std::shared_ptr<DataReader::Data> inD,
                       std::shared_ptr<Amplitude> inP,
                       std::shared_ptr<Optimizer::Optimizer> inO)
    : sampleData_(inD), amp_(inP), opti_(inO) {}

PythonFit::PythonFit(unsigned int size, std::shared_ptr<Amplitude> inP,
                       std::shared_ptr<Generator> gen)
    : amp_(inP), gen_(gen) {}*/

PythonFit::~PythonFit() {
  //LOG(debug) << "~RunManager: Last seed: " << gen_->getSeed();
}

std::shared_ptr<FitResult> PythonFit::startFit(ParameterList &inPar) {
  LOG(info) << "PythonFit::startFit() | Starting minimization.";

  // MINIMIZATION
  std::shared_ptr<FitResult> result = opti_->exec(inPar);

  LOG(info) << "PythonFit::startFit() | Minimization finished!"
               " Result = "
            << result->getResult() << ".";

  return result;
}

void PythonFit::StartFit() {
  LOG(info) << "PythonFit::StartFit() | call PythonFit::startFit";

  LOG(debug) << "PythonFit::StartFit() | test: member variable par has "
      << par.GetNParameter() << " parameter";

  LOG(debug) << "PythonFit::StartFit() | test: member variable sampleData_ has "
      << sampleData_->getNEvents() << " events";

  // MINIMIZATION
  std::shared_ptr<FitResult> result = startFit(par);

  return;
}

void PythonFit::setPhspSample(std::shared_ptr<Data> phsp,
                               std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && truePhsp->getNEvents() != phsp->getNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");
  samplePhsp_ = phsp;
  sampleTruePhsp_ = truePhsp;
}

void PythonFit::setTruePhspSample(std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && samplePhsp_ &&
      truePhsp->getNEvents() != samplePhsp_->getNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");

  sampleTruePhsp_ = truePhsp;
}

bool PythonFit::generate(int number) {
  LOG(info) << "RunManager::generate() | "
               "Generating "
            << number << " signal events!";

  return gen(number, gen_, amp_, sampleData_, samplePhsp_, sampleTruePhsp_);
}

bool PythonFit::generateBkg(int number) {
  LOG(info) << "RunManager::generateBkg() | "
               "Generating "
            << number << " background events!";

  return gen(number, gen_, ampBkg_, sampleBkg_, samplePhsp_, sampleTruePhsp_);
}

void PythonFit::SetAmplitudesData(
    std::vector<std::shared_ptr<Amplitude>> ampVec,
    std::vector<double> fraction,
    std::vector<std::shared_ptr<DataReader::Data>> dataVec) {
  _ampVec = ampVec;
  _fraction = fraction;
  _dataVec = dataVec;

  double sumFraction = std::accumulate(_fraction.begin(), _fraction.end(), 0.0);

  if (sumFraction > 1.0)
    throw std::runtime_error("RunManager::SetAmpltiudesData() | "
                             "Fractions sum larger 1.0!");

  if (_fraction.size() == _ampVec.size() - 1)
    _fraction.push_back(1 - sumFraction);

  // check size
  if (_fraction.size() != _ampVec.size())
    throw std::runtime_error("RunManager::SetAmpltiudesData() | "
                             "List of fractions (" +
                             std::to_string(_fraction.size()) +
                             ")"
                             " does not match with list of amplitudes"
                             "(" +
                             std::to_string(_ampVec.size()) + ")"
                                                              "!");

  if (_dataVec.size() != _ampVec.size())
    throw std::runtime_error("RunManager::SetAmpltiudesData() | "
                             "List of data samples (" +
                             std::to_string(_dataVec.size()) +
                             ")"
                             " does not match with list of amplitudes"
                             "(" +
                             std::to_string(_ampVec.size()) + ")"
                                                              "!");
}

void PythonFit::GenAmplitudesData(int nEvents) {
  if (nEvents <= 0)
    throw std::runtime_error(
        "RunManager::GenAmpltiudesData() | "
        "Number of events is small 0. Don't know how much to generate!");

  if (!_ampVec.size() || !_fraction.size() || !_dataVec.size())
    throw std::runtime_error(
        "RunManager::GenAmpltiudesData() | "
        "Amplitudes, fractions or data samples are not properly set!");

  if (_dataVec.size() != _ampVec.size() || _dataVec.size() != _fraction.size())
    throw std::runtime_error(
        "RunManager::GenAmpltiudesData() | "
        "Lists of data samples, fractions and amplitudes do "
        "not have the same size!");

  for (int i = 0; i < _ampVec.size(); ++i) {
    LOG(info) << "RunManager::GenAmplitudesData() | "
                 "Generating "
              << (int)(nEvents * _fraction.at(i)) << " events for amplitude "
              << _ampVec.at(i)->GetName();
    gen((int)(nEvents * _fraction.at(i)), gen_, _ampVec.at(i), _dataVec.at(i),
        samplePhsp_, sampleTruePhsp_);
  }
}

bool PythonFit::gen(int number, std::shared_ptr<Generator> gen,
                     std::shared_ptr<Amplitude> amp,
                     std::shared_ptr<DataReader::Data> data,
                     std::shared_ptr<DataReader::Data> phsp,
                     std::shared_ptr<DataReader::Data> phspTrue) {

  if (number == 0)
    return 0;

  // Doing some checks
  if (number < 0 && !phsp)
    throw std::runtime_error("RunManager: gen() negative number of events: " +
                             std::to_string((long double)number) +
                             ". And no phsp sample given!");
  if (!amp)
    throw std::runtime_error("RunManager::gen() | Amplitude not valid");
  if (!gen)
    throw std::runtime_error("RunManager::gen() | Generator not valid");
  if (!data)
    throw std::runtime_error("RunManager::gen() | Sample not valid");
  if (data->getNEvents() > 0)
    throw std::runtime_error("RunManager::gen() | Sample not empty!");
  if (phspTrue && !phsp)
    throw std::runtime_error("RunManager::gen() | We have a sample of true"
                             " phsp events, but no phsp sample!");
  if (phspTrue && phspTrue->getNEvents() != phsp->getNEvents())
    throw std::runtime_error(
        "RunManager::gen() | We have a sample of true "
        "phsp events, but the sample size doesn't match that one of "
        "the phsp sample!");

  double maxSampleWeight = 1.0;
  if (phsp)
    maxSampleWeight = phsp->getMaxWeight();
  if (phspTrue && phspTrue->getMaxWeight() > maxSampleWeight)
    maxSampleWeight = phspTrue->getMaxWeight();

  /* Maximum value for random number generation. We introduce an arbitrary
   * factor of 1.2 to make sure that the maximum value is never reached.
   */
  double generationMaxValue = 1.2 * amp->GetMaxVal(gen) * maxSampleWeight;

  unsigned int totalCalls = 0;
  //	unsigned int startTime = clock();
  double AMPpdf;

  unsigned int limit;
  unsigned int acceptedEvents = 0;
  if (phsp) {
    limit = phsp->getNEvents();
  } else
    limit = 100000000; // set large limit, should never be reached

  Event evt;     // event that we fill into generated sample
  Event evtTrue; // event that is used to evalutate amplitude
  progressBar bar(number);
  if (number <= 0)
    bar = progressBar(limit);
  for (unsigned int i = 0; i < limit; i++) {
    if (phsp && phspTrue) { // phsp and true sample is set
      evtTrue = phspTrue->getEvent(i);
      evt = phsp->getEvent(i);
    } else if (phsp) { // phsp sample is set
      evt = phsp->getEvent(i);
      evtTrue = evt;
    } else { // otherwise generate event
      gen->generate(evt);
      evtTrue = evt;
    }
    if (number <= 0)
      bar.nextEvent();

    // use reconstructed position for weights
    double weight = evt.getWeight();

    // use true position for amplitude value
    dataPoint point;
    try {
      point = dataPoint(evtTrue);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //		if(!Kinematics::instance()->isWithinPhsp(point)) continue;

    totalCalls++;
    double ampRnd = gen->getUniform() * generationMaxValue;
    ParameterList list;
    list = amp->intensity(point); // unfortunatly not thread safe
    AMPpdf = *list.GetDoubleParameter(0);

    if (generationMaxValue < (AMPpdf * weight))
      throw std::runtime_error(
          "RunManager::gen() | Error in HitMiss "
          "procedure: Maximum value of random number generation "
          "smaller then amplitude maximum!");

    if (ampRnd > (weight * AMPpdf))
      continue;

    // Fill event to sample
    /* reset weights: the weights are taken into account by hit and miss. The
     * resulting
     * sample is therefore unweighted */
    evt.setWeight(1.);     // reset weight
    evt.setEfficiency(1.); // reset weight
    data->pushEvent(evt);  // unfortunatly not thread safe

    // some statistics
    acceptedEvents++;
    if (number > 0)
      bar.nextEvent();
    // break if we have a sufficienct number of events
    if (acceptedEvents >= number)
      i = limit;
  }
  if (data->getNEvents() < number)
    LOG(error) << "RunManager::gen() | Not able to generate "
                  " "
               << number << " events. Phsp sample too small. Current size "
                            "of sample is now "
               << data->getNEvents();

  LOG(info) << "Efficiency of toy MC generation: "
            << (double)data->getNEvents() / totalCalls;

  return true;
}

bool PythonFit::generatePhsp(int number) {
  if (number == 0)
    return 0;
  if (!samplePhsp_)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "No phase-space sample set");
  if (samplePhsp_->getNEvents() > 0)
    throw std::runtime_error("RunManager: generatePhsp() | "
                             "Dataset not empty! abort!");

  LOG(info) << "Generating phase-space MC: [" << number << " events] ";

  progressBar bar(number);
  for (unsigned int i = 0; i < number; i++) {
    if (i > 0)
      i--;
    Event tmp;
    gen_->generate(tmp);
    double ampRnd = gen_->getUniform();
    if (ampRnd > tmp.getWeight())
      continue;

    /* Reset weights: weights are taken into account by hit&miss. The
     * resulting sample is therefore unweighted
     */
    tmp.setWeight(1.);

    tmp.setEfficiency(1.);
    i++;
    samplePhsp_->pushEvent(tmp); // unfortunatly not thread safe
    bar.nextEvent();
  }
  return true;
}