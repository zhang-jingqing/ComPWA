// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

// This can only be define once within the same library ?!
#define BOOST_TEST_MODULE ComponentIntensityTest

#include <vector>

#include "Core/Intensity.hpp"
#include "Core/Logging.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Properties.hpp"
#include "Core/FitParameter.hpp"
#include "Data/DataSet.hpp"
#include "Physics/Amplitude.hpp"
#include "Physics/Dynamics/RelativisticBreitWigner.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IntensityBuilderXML.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/ComponentIntensity.hpp"

#include <boost/foreach.hpp>
#include <boost/locale/utf.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/string.hpp>

using namespace ComPWA::Physics;
using namespace ComPWA::Tools;

BOOST_AUTO_TEST_SUITE(ToolsTest)

const std::string TestParticles = R"####(
<ParticleList>
  <Particle Name='pi0'>
    <Pid>111</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_pi0'>
      <Value>0.1349766</Value>
      <Error>0.0000006</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
  </Particle>
  <Particle Name='gamma'>
    <Pid>22</Pid>
    <Parameter Class='Double' Type='Mass' Name='mass_gamma'>
      <Value>0.</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
  </Particle>
  <Particle Name='jpsi'>
    <Pid>443</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_jpsi'>
      <Value>3.0969</Value>
      <Fix>true</Fix>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='1'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='-1'/>
    <QuantumNumber Class='Int' Type='Gparity' Value='1'/>
    <DecayInfo Type='nonResonant'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_jpsi'>
        <Value>0.0000929</Value>
        <Error>0.0000028</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_jpsi'>
        <Value>2.5</Value>
        <Fix>true</Fix>
        <Min>2.0</Min>
        <Max>3.0</Max>
      </Parameter>
    </DecayInfo>
  </Particle>
  <Particle Name='f0'>
    <Pid>6666</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_f0'>
      <Value>1.500</Value>
      <Fix>true</Fix>
      <Min>0.5</Min>
      <Max>1.5</Max>
      <Error>0.00012</Error>
    </Parameter>
    <QuantumNumber Class='Spin' Type='Spin' Value='0'/>
    <QuantumNumber Class='Int' Type='Charge' Value='0'/>
    <QuantumNumber Class='Int' Type='Parity' Value='1'/>
    <QuantumNumber Class='Int' Type='Cparity' Value='1'/>
    <DecayInfo Type='relativisticBreitWigner'>
      <FormFactor Type='1' />
      <Parameter Class='Double' Type='Width' Name='Width_f0'>
        <Value>0.050</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_f0'>
        <Value>1.5</Value>
        <Fix>true</Fix>
        <Min>1.0</Min>
        <Max>2.0</Max>
        <Error>0</Error>
      </Parameter>
    </DecayInfo>
  </Particle>
</ParticleList>
)####";

const std::string JpsiDecayKinematics = R"####(
<HelicityKinematics>
  <PhspVolume>1</PhspVolume>
  <InitialState>
    <Particle Name='jpsi' PositionIndex='0'/>
  </InitialState>
  <FinalState>
    <Particle Name='gamma' Id='0'/>
    <Particle Name='pi0' Id='1'/>
    <Particle Name='pi0' Id='2'/>
  </FinalState>
</HelicityKinematics>
)####";
//(1,-1)  (1,-1) (0)
//jpsi -> gamma f_0
//F_1,0 = F_-1,0
//F_1,0 = p1 (0 0 1 1|1 1)(1 1 0 0|1 1) a_01
//      + p2 (2 0 1 1|1 1)(1 1 0 0|1 1) a_21
//F_-1,0 = p1 (0 0 1 -1|1 -1)(1 -1 0 0|1 -1) a_01
//      + p2 (2 0 1 -1|1 -1)(1 -1 0 0|1 -1) a_21
//f_0 -> pi0 pi0
//F_0,0 = p * (0 0 0 0|0 0) a_0,0
const std::string JpsiDecayTree = R"####(
<Intensity Class='StrengthIntensity' Name='incoherent_with_strength'>
  <Parameter Class='Double' Type='Strength' Name='strength_incoherent'>
    <Value>1</Value>
    <Fix>True</Fix>
  </Parameter>
  <Intensity Class='IncoherentIntensity' Name='incoherent'>
    <!-- Helicity: (1, 1, 0)(0, 0, 0); -->
    <Intensity Class='CoherentIntensity' Name='coherent_0'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='1.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0_to_pi0_pi0'>
            <DecayParticle Name='f0' Helicity='0' />
            <DecayProducts>
              <Particle Name='pi0' FinalState='1' Helicity='0' />
              <Particle Name='pi0' FinalState='2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='0.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
              <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>

      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='2.0' S='1.0'>
              <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0_to_pi0_pi0'>
            <DecayParticle Name='f0' Helicity='0' />
            <DecayProducts>
              <Particle Name='pi0' FinalState='1' Helicity='0' />
              <Particle Name='pi0' FinalState='2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='0.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
              <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>
    </Intensity>

    <!-- Helicity: (-1, 1, 0)(0, 0, 0); -->
    <Intensity Class='CoherentIntensity' Name='coherent_1'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0'>
            <DecayParticle Name='jpsi' Helicity='-1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='1.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0_to_pi0_pi0'>
            <DecayParticle Name='f0' Helicity='0' />
            <DecayProducts>
              <Particle Name='pi0' FinalState='1' Helicity='0' />
              <Particle Name='pi0' FinalState='2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='0.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
              <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_0.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_0.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>

      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_-1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0'>
            <DecayParticle Name='jpsi' Helicity='-1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='2.0' S='1.0'>
              <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0_to_pi0_pi0'>
            <DecayParticle Name='f0' Helicity='0' />
            <DecayProducts>
              <Particle Name='pi0' FinalState='1' Helicity='0' />
              <Particle Name='pi0' FinalState='2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='0.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
              <ClebschGordan Type='s2s3' j1='0.0' m1='0.0' j2='0.0' m2='0.0' J='0.0' M='0.0'/>
            </CanonicalSum>
          </Amplitude>
        </Amplitude>
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0_L_2.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0_L_2.0_S_1.0;f0_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>
    </Intensity>

  </Intensity>
</Intensity>
)####";

//jpsi-> gamma f_0, f_0 -> pi0 pi0
//(a_01, a_02)(a_00)
//1. all f_0 contributions to total intensity(should be same)
//2. e.g., a_01 of f_0 (should be same to the one created from xml tree)
//3. extract component of a decay chain, e.g. a_01*a_00, or f_0 * a_00

BOOST_AUTO_TEST_CASE(ComponentIntensityTest) {
  ComPWA::Logging log("", "trace");

  LOG(INFO) << "Now check ComonentIntensity...";

  boost::property_tree::ptree modelTree;
  std::stringstream modelStream;
  ComPWA::Physics::IntensityBuilderXML Builder;

  //particle list
  modelStream << TestParticles;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, modelTree);

  //kinematics
  modelStream.clear();
  modelTree = boost::property_tree::ptree();
  modelStream << JpsiDecayKinematics;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);
  auto decayKin = Builder.createHelicityKinematics(partL,
      modelTree.get_child("HelicityKinematics"));

  //model intensity
  modelStream.clear();
  modelTree = boost::property_tree::ptree();
  modelStream << JpsiDecayTree;
  boost::property_tree::xml_parser::read_xml(modelStream, modelTree);
  auto modelIntensity = Builder.createIntensity(partL, decayKin,
      modelTree.get_child("Intensity"));
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(
      decayKin->getParticleStateTransitionKinematicsInfo(), 233);
  std::shared_ptr<ComPWA::Data::DataSet> sample =
      ComPWA::Tools::generatePhsp(20, gen);
  sample->convertEventsToDataPoints(decayKin);

  //to find components, 
  //for each decay node in a component, one need provide 
  //  std::vector<std::string>{mother, daughter1, daughter2},
  //  "" means no requirement
  //  std::pair<int, int>(Lmin, Lmax), std::pair<int, int>(Smin, Smax)
  //  negative L/S means no requirement
  //getComponent*() need component's decay paterns provided in below form:
  //std::vector<std::vector<std::string>>
  //std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>
  //getComponent*() also need (In)coherentIntensity as input

  auto incoherentIntensity = Builder.createIntensity(partL, decayKin,
      modelTree.get_child("Intensity.Intensity"));

  //ensure both modelIntensity and incoherentIntensity will use same set of
  //parameters(and parameter strengthIntensity should be fixed)
  ComPWA::ParameterList parameters;
  modelIntensity->addUniqueParametersTo(parameters);
  incoherentIntensity->updateParametersFrom(parameters);

  //check: get all contributions of f_0
  std::vector<std::vector<std::string>> subDecays(1,
      std::vector<std::string>({"f0", "", ""}));
  std::vector<std::pair<int, int>> subLRanges(1, std::pair<int, int>(-1, -1));
  std::vector<std::pair<int, int>> subSRanges(1, std::pair<int, int>(-1, -1));
  auto f0Intensity = getComponentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "f0", subDecays, subLRanges, subSRanges);
  
  for (const auto &point : sample->getDataPointList()) {
    double value1 = modelIntensity->evaluate(point);
    double value2 = f0Intensity->evaluate(point);
    LOG(INFO) << "< Jpsi -> gamma f_0, f_0 -> pi0 pi0:";
    LOG(INFO) << " total intensity value: " << value1;
    LOG(INFO) << " all f0 contributions : " << value2 << " >";
    BOOST_CHECK_EQUAL(value1, value2);
  }

  //check: get a_01 of jpsi -> gamma f_0 
  subDecays.at(0) = std::vector<std::string>({"jpsi", "gamma", "f0"});
  subLRanges.at(0) = std::pair<int, int>(0, 0);
  subSRanges.at(0).first = 1;
  subSRanges.at(0).second = 1;
  auto a01Intensity = getComponentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "a01", subDecays, subLRanges, subSRanges);
  //make a intensity tree with only a01
  boost::property_tree::ptree a01Tree(modelTree);
  for (auto &ch1 : a01Tree.get_child("Intensity.Intensity")) {
    //loop over coherent intensity
    if (ch1.first != "Intensity") continue;
    //loop over contents of coherent intensity
    boost::property_tree::ptree &node = ch1.second;
    for (auto it = node.begin(); it != node.end(); ) {
      if (it->first != "Amplitude") {
        ++it;
        continue;
      }
      std::string coeffAmpNames = it->second.get<std::string>("<xmlattr>.Name");
      std::vector<std::string> decayNames;
      boost::algorithm::split(decayNames, coeffAmpNames,
          boost::is_any_of(";"));
      LOG(INFO) << coeffAmpNames;
      bool found = findAllDecayPatern(subDecays, subLRanges, subSRanges,
          decayNames);
      if (found) {
        ++it;
      } else { //remove not a01 amplitudes
        it = node.erase(it);
        break;
      }
    }
  }
  //get a intensity only has a01 from this tree
  auto a01IntensityFromTree = Builder.createIntensity(partL, decayKin,
      a01Tree.get_child("Intensity"));
  a01Intensity->updateParametersFrom(parameters);

  for (const auto &point : sample->getDataPointList()) {
    double value1 = a01IntensityFromTree->evaluate(point);
    double value2 = a01Intensity->evaluate(point);
    LOG(INFO) << "< Jpsi -> gamma f_0, f_0 -> pi0 pi0, only a01 amplitude:";
    LOG(INFO) << " extracted from total intensity: " << value1;
    LOG(INFO) << " built from a01 model tree     : " << value2 << " >";
    BOOST_CHECK_EQUAL(value1, value2);
  }

  //check: get a_01 of jpsi -> gamma f_0 with a_00 of f_0 -> pi0 pi0
  subDecays.at(0) = std::vector<std::string>({"jpsi", "gamma", "f0"});
  subLRanges.at(0) = std::pair<int, int>(0, 1);
  subSRanges.at(0) = std::pair<int, int>(1, 1);
  //the second sub decay
  subDecays.push_back(std::vector<std::string>({"f0", "pi0", "pi0"}));
  subLRanges.push_back(std::pair<int, int>(0, 0));
  subSRanges.push_back(std::pair<int, int>(0, 0));
  auto a01_a00_intensity = getComponentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "a01xa00", subDecays, subLRanges, subSRanges);

  for (const auto &point : sample->getDataPointList()) {
    double value1 = a01Intensity->evaluate(point);
    double value2 = a01_a00_intensity->evaluate(point);
    BOOST_CHECK_EQUAL(value1, value2);
  }

  //check the find functions:
  // the LS combination is randomly chosen, may be not possible to happen
  std::vector<std::string> testAmpNames({
      "jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;",
      "jpsi_1_to_gamma_1+f0_0_L_0.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_0.0;",

      "jpsi_1_to_gamma_1+f2_0_L_2.0_S_1.0;f2_0_to_pi0_0+pi0_0_L_1.0_S_0.0;",

      "jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_3.0_S_2.0;",
      "jpsi_1_to_gamma_1+f0_0_L_2.0_S_1.0;f0_0_to_pi0_0+pi0_0_L_0.0_S_1.0;",
  });

  int f0Count_all = 0, f0Count_a21 = 0, f0Count_a21_a32 = 0;
  for (const auto &x : testAmpNames) {
    std::vector<std::string> decayNames;
    boost::algorithm::split(decayNames, x, boost::is_any_of(";"));

    //bool found = findAllDecayPatern(subDecays, 
    //    subLRanges, subSRanges, decayNames);
    //find all f0
    bool found = findAllDecayPatern(
        std::vector<std::vector<std::string>>(1,
            std::vector<std::string>({"f0", "", ""})),
        std::vector<std::pair<int, int>>(1, std::pair<int, int>(-1, -1)),
        std::vector<std::pair<int, int>>(1, std::pair<int, int>(-1, -1)),
        decayNames);
    if (found) f0Count_all++;

    found = false;
    //find a21_ jpis->gamma f0
    found = findAllDecayPatern(
        std::vector<std::vector<std::string>>(1,
            std::vector<std::string>({"jpsi", "gamma", "f0"})),
        std::vector<std::pair<int, int>>(1, std::pair<int, int>(2, 2)),
        std::vector<std::pair<int, int>>(1, std::pair<int, int>(1, 1)),
        decayNames);
    if (found) f0Count_a21++;
  
    found = false;
    //find a21 jpsi->gamma f0, a32 f0 -> pi0 pi0
    found = findAllDecayPatern(
        std::vector<std::vector<std::string>>({
            std::vector<std::string>({"jpsi", "gamma", "f0"}),
            std::vector<std::string>({"f0", "pi0", "pi0"})}),
        std::vector<std::pair<int, int>>({std::pair<int, int>(2, 2),
            std::pair<int, int>(3, 3)}),
        std::vector<std::pair<int, int>>({std::pair<int, int>(1, 1),
            std::pair<int, int>(2, 2)}),
        decayNames);
    if (found) f0Count_a21_a32++;
  }

  BOOST_CHECK_EQUAL(f0Count_all, 4);
  BOOST_CHECK_EQUAL(f0Count_a21, 2);
  BOOST_CHECK_EQUAL(f0Count_a21_a32, 1);
}


BOOST_AUTO_TEST_SUITE_END()
