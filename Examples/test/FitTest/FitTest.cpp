// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#define BOOST_TEST_MODULE FitTest

#include <boost/test/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/ParticleList.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Generate.hpp"
#include "Tools/ParameterTools.hpp"
#include "Tools/Plotting/ROOT/RootPlotData.hpp"
#include "Tools/RootGenerator.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Optimizer/Minuit2/MinuitIF.hpp"

using namespace ComPWA;
using ComPWA::Optimizer::Minuit2::MinuitResult;
using ComPWA::Physics::IncoherentIntensity;
using ComPWA::Physics::HelicityFormalism::HelicityKinematics;

// We define an intensity model using a raw string literal. Currently, this is
// just a toy model without any physical meaning.
// (comments within the string are ignored!). This is convenient since we
// do not have to configure the build system to copy input files somewhere.
// In practise you may want to use a normal XML input file instead.
std::string amplitudeModel = R"####(
<Intensity Class='Incoherent' Name="jpsiGammaPiPi_inc">
  <Intensity Class='Coherent' Name="jpsiGammaPiPi">
    <Amplitude Class="SequentialPartialAmplitude" Name="f2(1270)">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_f2">
        <Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_f2">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>false</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name="f2ToPiPi">
        <DecayParticle Name="f2(1270)" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
    <Amplitude Class="SequentialPartialAmplitude" Name="myAmp">
      <Parameter Class='Double' Type="Magnitude"  Name="Magnitude_my">
        <Value>1.0</Value>
        <Min>-1.0</Min>
        <Max>2.0</Max>
        <Fix>true</Fix>
      </Parameter>
      <Parameter Class='Double' Type="Phase" Name="Phase_my`">
        <Value>0.0</Value>
        <Min>-100</Min>
        <Max>100</Max>
        <Fix>true</Fix>
      </Parameter>
      <PartialAmplitude Class="HelicityDecay" Name="MyResToPiPi">
        <DecayParticle Name="myRes" Helicity="0"/>
        <RecoilSystem FinalState="0" />
        <DecayProducts>
          <Particle Name="pi0" FinalState="1"  Helicity="0"/>
          <Particle Name="pi0" FinalState="2"  Helicity="0"/>
        </DecayProducts>
      </PartialAmplitude>
    </Amplitude>
  </Intensity>
</Intensity>
)####";

std::string myParticles = R"####(
<ParticleList>
  <Particle Name="f2(1270)">
    <Pid>225</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_f2(1270)">
      <Value>1.2755</Value>
      <Error>8.0E-04</Error>
      <Min>0.1</Min>
      <Max>2.0</Max>
      <Fix>false</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="2"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
      <Parameter Class='Double' Type="Width" Name="Width_f2(1270)">
        <Value>0.1867</Value>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_rho">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name="myRes">
    <Pid>999999</Pid>
    <Parameter Class='Double' Type="Mass" Name="Mass_myRes">
      <Value>2.0</Value>
      <Error>8.0E-04</Error>
      <Min>1.1</Min>
      <Max>4.0</Max>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class="Spin" Type="Spin" Value="1"/>
    <QuantumNumber Class="Int" Type="Charge" Value="0"/>
    <QuantumNumber Class="Int" Type="Parity" Value="+1"/>
    <QuantumNumber Class="Int" Type="Cparity" Value="+1"/>
    <DecayInfo Type="relativisticBreitWigner">
      <FormFactor Type="0" />
      <Parameter Class='Double' Type="Width" Name="Width_myRes">
        <Value>1.0</Value>
        <Min>0.1</Min>
        <Max>1.0</Max>
        <Fix>false</Fix>
      </Parameter>
      <Parameter Class='Double' Type="MesonRadius" Name="Radius_myRes">
        <Value>2.5</Value>
        <Fix>true</Fix>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

BOOST_AUTO_TEST_SUITE(FitTest)

BOOST_AUTO_TEST_CASE(HelicityDalitzFit) {

  // initialize logging
  // Logging log("DalitzFit-log.txt", boost::log::trivial::debug);
  // boost::test_tools::output_test_stream output;

  // List with all particle information needed
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, defaultParticleList);
  ReadParticles(partL, myParticles);

  //---------------------------------------------------
  // 1) Create Kinematics object
  //---------------------------------------------------
  std::vector<pid> initialState = {443};
  std::vector<pid> finalState = {22, 111, 111};
  auto kin =
      std::make_shared<HelicityKinematics>(partL, initialState, finalState);

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(partL, kin, 173);
  std::shared_ptr<ComPWA::Data::Data> phspSample =
      ComPWA::Tools::generatePhsp(100000, gen);

  //---------------------------------------------------
  // 3) Create intensity from pre-defined model
  //---------------------------------------------------
  // Read in model property_tree
  std::stringstream modelStream;
  modelStream << amplitudeModel;
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);

  // Construct intensity class from model string
  auto intens = std::make_shared<IncoherentIntensity>(
      partL, kin, modelTree.get_child("Intensity"));

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  auto phspPoints =
      std::make_shared<std::vector<DataPoint>>(phspSample->dataPoints(kin));
  intens->setPhspSample(phspPoints, phspPoints);

  //---------------------------------------------------
  // 4) Generate a data sample given intensity and kinematics
  //---------------------------------------------------
  gen->setSeed(1234);
  std::shared_ptr<ComPWA::Data::Data> sample =
      ComPWA::Tools::generate(1000, kin, gen, intens, phspSample);
  //---------------------------------------------------
  // 5) Fit the model to the data and print the result
  //---------------------------------------------------
  ParameterList FitParameters;
  intens->parameters(FitParameters);
  // Set start error of 0.05 for parameters, run Minos?
  setErrorOnParameterList(FitParameters, 0.05, false);
  auto esti = std::make_shared<Estimator::MinLogLH>(
      kin, intens, sample, phspSample, phspSample, 0, 0);

  esti->UseFunctionTree(true);
  //LOG(INFO) << esti->tree()->head()->print(25);

  auto minuitif = new Optimizer::Minuit2::MinuitIF(esti, FitParameters);
  minuitif->setUseHesse(true);

  // STARTING MINIMIZATION
  auto result =
      std::dynamic_pointer_cast<MinuitResult>(minuitif->exec(FitParameters));

  // Calculate fit fractions
  // std::vector<std::pair<std::string, std::string>> fitComponents;
  // fitComponents.push_back(
  //   std::pair<std::string, std::string>("myAmp", "jpsiGammaPiPi"));
  // fitComponents.push_back(
  // std::pair<std::string, std::string>("f2(1270)", "jpsiGammaPiPi"));

  // ParameterList fitFracs =
  // Tools::CalculateFitFractions(kin, intens, phspPoints, fitComponents);
  // A proper calculation of the fit fraction uncertainty requires
  // the uncertainty of the fit parameters propagated. We do this
  // using a numerical approach. Using the covariance matrix
  // 100 independend sets of fit parameters are generated and the fit fractions
  // are recalculated. In the end we take the RMS.
  // Tools::CalcFractionError(fitPar, result->GetCovarianceMatrix(), fitFracs,
  // kin, intens, phspPoints, 100, fitComponents);

  // result->SetFitFractions(fitFracs);
  // result->Print();

  //---------------------------------------------------
  // 5.1) Save the fit result
  //---------------------------------------------------
  // std::ofstream ofs("DalitzFit-fitResult.xml");
  // boost::archive::xml_oarchive oa(ofs);
  // oa << BOOST_SERIALIZATION_NVP(result);

  UpdateParticleList(partL, FitParameters);
  // boost::property_tree::ptree ptout;
  // ptout.add_child("ParticleList", SaveParticles(partL));
  // ptout.add_child("IncoherentIntensity", IncoherentIntensity::Save(intens));
  // boost::property_tree::xml_parser::write_xml("DalitzFit-Model.xml", ptout,
  // std::locale());
  //---------------------------------------------------
  // 6) Plot data sample and intensity
  //---------------------------------------------------
  // ComPWA::Tools::DalitzPlot pl(kin, "DalitzFit", 100);
  // pl.SetData(sample);
  // pl.SetPhspData(phspSample);
  // pl.SetFitAmp(intens, "", kBlue - 4);
  // pl.Plot();
  // LOG(info) << "Done";

  // output << result->finalLH();
  BOOST_CHECK_EQUAL(sample->numEvents(), 1000);
  BOOST_CHECK_CLOSE(result->finalLH(), -980, 5.); // 5% tolerance
  double sigma(2.0);
  auto fitpar = FindParameter("Magnitude_f2", FitParameters);
  BOOST_CHECK_GT((fitpar->value() + sigma * fitpar->error().second), 1.0);
  BOOST_CHECK_GT(1.0, (fitpar->value() - sigma * fitpar->error().first));
  fitpar = FindParameter("Phase_f2", FitParameters);
  BOOST_CHECK_GT((fitpar->value() + sigma * fitpar->error().second), 0.0);
  BOOST_CHECK_GT(0.0, (fitpar->value() - sigma * fitpar->error().first));
  fitpar = FindParameter("Mass_f2(1270)", FitParameters);
  BOOST_CHECK_GT((fitpar->value() + sigma * fitpar->error().second), 1.2755);
  BOOST_CHECK_GT(1.2755, (fitpar->value() - sigma * fitpar->error().first));
  fitpar = FindParameter("Width_myRes", FitParameters);
  BOOST_CHECK_GT((fitpar->value() + sigma * fitpar->error().second), 1.0);
  BOOST_CHECK_GT(1.0, (fitpar->value() - sigma * fitpar->error().first));
};
BOOST_AUTO_TEST_SUITE_END()