#define BOOST_TEST_MODULE FitFractionsTest

#include <vector>
#include <string>

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
#include "Physics/IncoherentIntensity.hpp"
#include "Tools/Generate.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/ComponentIntensity.hpp"
#include "Tools/FitFractions.hpp"
#include "Tools/Integration.hpp"

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

  <Particle Name='f0clone'>
    <Pid>6666</Pid>
    <Parameter Class='Double' Type='Mass' Name='Mass_f0clone'>
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
      <Parameter Class='Double' Type='Width' Name='Width_f0clone'>
        <Value>0.050</Value>
        <Fix>true</Fix>
        <Min>0.0</Min>
        <Max>1.0</Max>
        <Error>0.00008</Error>
      </Parameter>
      <Parameter Class='Double' Type='MesonRadius' Name='Radius_f0clone'>
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

const std::string JpsiDecayTree = R"####(
<Intensity Class='StrengthIntensity' Name='incoherent_with_strength'>
  <Parameter Class='Double' Type='Strength' Name='strength_incoherent'>
    <Value>1</Value>
    <Fix>True</Fix>
  </Parameter>
  <Intensity Class='IncoherentIntensity' Name='incoherent'>
    <!-- Helicity: f0 -->
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

    <!-- Helicity: f0clone -->
    <!-- Helicity: (1, 1, 0)(0, 0, 0); -->
    <Intensity Class='CoherentIntensity' Name='coherent_2'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0clone_0_L_0.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0clone_0_L_0.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0clone'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0clone' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='1.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0clone_to_pi0_pi0'>
            <DecayParticle Name='f0clone' Helicity='0' />
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
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0clone_L_0.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0clone_L_0.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>

      <Amplitude Class='CoefficientAmplitude' Name='jpsi_1_to_gamma_1+f0clone_0_L_2.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_1_to_gamma_1+f0clone_0_L_2.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0clone'>
            <DecayParticle Name='jpsi' Helicity='+1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0clone' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='2.0' S='1.0'>
              <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0clone_to_pi0_pi0'>
            <DecayParticle Name='f0clone' Helicity='0' />
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
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0clone_L_2.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0clone_L_2.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>
    </Intensity>

    <!-- Helicity: (-1, 1, 0)(0, 0, 0); -->
    <Intensity Class='CoherentIntensity' Name='coherent_3'>
      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0clone_0_L_0.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_-1_to_gamma_1+f0clone_0_L_0.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0clone'>
            <DecayParticle Name='jpsi' Helicity='-1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0clone' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='0.0' S='1.0'>
              <ClebschGordan Type='LS' j1='0.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0clone_to_pi0_pi0'>
            <DecayParticle Name='f0clone' Helicity='0' />
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
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0clone_L_0.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0clone_L_0.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>

      <Amplitude Class='CoefficientAmplitude' Name='jpsi_-1_to_gamma_1+f0clone_0_L_2.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
        <Amplitude Class='SequentialAmplitude' Name='jpsi_-1_to_gamma_1+f0clone_0_L_2.0_S_1.0;f0clone_0_to_pi0_0+pi0_0_L_0.0_S_0.0;'>
          <Amplitude Class='HelicityDecay' Name='jpsi_to_gamma_f0clone'>
            <DecayParticle Name='jpsi' Helicity='-1' />
            <DecayProducts>
              <Particle Name='gamma' FinalState='0' Helicity='1' />
              <Particle Name='f0clone' FinalState='1 2' Helicity='0' />
            </DecayProducts>
            <CanonicalSum L='2.0' S='1.0'>
              <ClebschGordan Type='LS' j1='2.0' m1='0.0' j2='1.0' m2='1.0' J='1.0' M='1.0' />
              <ClebschGordan Type='s2s3' j1='1.0' m1='1.0' j2='0.0' m2='0.0' J='1.0' M='1.0' />
            </CanonicalSum>
          </Amplitude>
          <Amplitude Class='HelicityDecay' Name='f0clone_to_pi0_pi0'>
            <DecayParticle Name='f0clone' Helicity='0' />
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
        <Parameter Class='Double' Type='Magnitude' Name='Magnitude_jpsi_to_gamma+f0clone_L_2.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>1.0</Value>
          <Fix>False</Value>
        </Parameter>
        <Parameter Class='Double' Type='Phase' Name='Phase_jpsi_to_gamma+f0clone_L_2.0_S_1.0;f0clone_to_pi0+pi0_L_0.0_S_0.0;'>
          <Value>0.0</Value>
          <Fix>False</Value>
        </Parameter>
      </Amplitude>
    </Intensity>

  </Intensity>
</Intensity>
)####";

BOOST_AUTO_TEST_CASE(FitFractionsTest) {
  ComPWA::Logging log("", "trace");

  LOG(INFO) << "Now check FitFractions and ComponentIntensity...";

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

  //check: fit fractions of all f_0; should be 0.5
  std::vector<std::vector<std::string>> subDecays(1,
      std::vector<std::string>({"jpsi", "gamma", "f0"}));
  std::vector<std::pair<int, int>> subLRanges(1, std::pair<int, int>(-1, -1));
  std::vector<std::pair<int, int>> subSRanges(1, std::pair<int, int>(-1, -1));
  auto f0Intensity = getComponentIntensityFromIncoherentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "f0_all", subDecays, subLRanges, subSRanges);
  std::vector<std::shared_ptr<const ComPWA::Intensity>> componentIntensities;
  componentIntensities.push_back(f0Intensity);
  auto f0AllFF = ComPWA::Tools::calculateFitFractions(
      incoherentIntensity, sample, std::vector<std::string>(1, "f0_all"),
      componentIntensities);
  double f0_all_ff = f0AllFF.doubleParameter(0)->value();
  
  LOG(INFO) << "f0 (all contributions) fitfraction(should be 0.5): " << f0_all_ff;
  BOOST_CHECK_EQUAL(f0_all_ff, 0.5);

  subDecays.at(0) = std::vector<std::string>({"jpsi", "gamma", "f0clone"});
  auto f0cloneIntensity = getComponentIntensityFromIncoherentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "f0clone_all", subDecays, subLRanges, subSRanges);
  componentIntensities.at(0) = f0cloneIntensity;
  auto f0cloneAllFF = ComPWA::Tools::calculateFitFractions(
      incoherentIntensity, sample, std::vector<std::string>(1, "f0clone_all"),
      componentIntensities);
  double f0clone_all_ff = f0cloneAllFF.doubleParameter(0)->value();

  LOG(INFO) << "f0clone (all contributions) fitfraction(should be 0.5): "
      << f0clone_all_ff;
  BOOST_CHECK_EQUAL(f0clone_all_ff, 0.5);
      

  //check: get a_01 of jpsi -> gamma f_0 
  subDecays.at(0) = std::vector<std::string>({"jpsi", "gamma", "f0"});
  subLRanges.at(0) = std::pair<int, int>(0, 0);
  subSRanges.at(0) = std::pair<int, int>(1, 1);
  auto f0_a01Intensity = getComponentIntensityFromIncoherentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "f0_a01", subDecays, subLRanges, subSRanges);
  // get a_01 of jpsi -> gamma f0clone
  subDecays.at(0) = std::vector<std::string>({"jpsi", "gamma", "f0clone"});
  auto f0clone_a01Intensity = getComponentIntensityFromIncoherentIntensity(
      std::dynamic_pointer_cast<ComPWA::Physics::IncoherentIntensity>
      (incoherentIntensity),
      "f0clone_a01", subDecays, subLRanges, subSRanges);

  componentIntensities.at(0) = f0_a01Intensity;
  componentIntensities.push_back(f0clone_a01Intensity);
  auto a01_ff = ComPWA::Tools::calculateFitFractions(
      incoherentIntensity, sample, 
      std::vector<std::string>({"f0_a01", "f0clone_a01"}),
      componentIntensities);

 double f0_a01_ff = a01_ff.doubleParameter(0)->value();
 std::string f0_a01_name = a01_ff.doubleParameter(0)->name();
 double f0clone_a01_ff = a01_ff.doubleParameter(1)->value();
 std::string f0clone_a01_name = a01_ff.doubleParameter(1)->name();

 LOG(INFO) << "fit fraction of " << f0_a01_name << ":\t " << f0_a01_ff;
 LOG(INFO) << "fit fraction of " << f0clone_a01_name << ":\t " << f0_a01_ff;
 BOOST_CHECK_EQUAL(f0_a01_ff, f0clone_a01_ff);
 BOOST_CHECK_EQUAL(f0_a01_name, "f0_a01");
 BOOST_CHECK_EQUAL(f0clone_a01_name, "f0clone_a01");
}


BOOST_AUTO_TEST_SUITE_END()
