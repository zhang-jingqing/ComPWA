//! Application to generate J/Psi -> g pi pi.
/*!
 * @file HelicityCoherentAmplitudeGenApp.cpp
 * This application uses the HelicityAmplitude module and a
 * phase-space generator to generate a file with J/Psi -> gamma pi0 pi0 events.
 * Also uses the neat plotting tool to create a dalitz plot of the data.
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include <boost/filesystem.hpp>

#include "DataReader/RootGenerator/RootGenerator.hpp"
// Physics Interface header files go here
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/RunManager.hpp"
#include "Core/Efficiency.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Core/ProgressBar.hpp"

#include "Physics/DecayTree/DecayConfiguration.hpp"
#include "Physics/DecayTree/DecayXMLConfigReader.hpp"
#include "Physics/DecayTree/DecayTreeFactory.hpp"
#include "Physics/HelicityAmplitude/TopologyAmplitudeFactory.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"
#include "Physics/HelicityAmplitude/CoherentAmplitude.hpp"

/*#include "Tools/RootNeatPlotting/HelperFunctions.h"
 #include "Tools/RootNeatPlotting/plotting/PlotBundle.h"
 #include "Tools/RootNeatPlotting/plotting/Booky.h"
 #include "Tools/RootNeatPlotting/style/DataObjectStyle.h"
 #include "Tools/RootNeatPlotting/style/DefaultStyleSingleton.h"
 #include "Tools/RootNeatPlotting/style/xml-parser/XMLStyleConfigParser.h"*/
//#include "PWA/PlotData.hpp"
using namespace ComPWA;
using DataReader::Data;
using DataReader::RootReader::RootReader;
using DataReader::RootGenerator::RootGenerator;

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv) {
  ComPWA::Physics::HelicityFormalism::HelicityKinematics* kin =
      dynamic_cast<ComPWA::Physics::HelicityFormalism::HelicityKinematics*>(ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance());

  // we expect at least 1 argument
  if (argc > 1) {
    std::string input_config_file = argv[1];

    boost::filesystem::path output_path(input_config_file);
    if (!boost::filesystem::is_regular_file(output_path)) {
      std::runtime_error(
          "Specified xml config file does not exist! Please specify a valid url...");
    }
    output_path = boost::filesystem::absolute(output_path,
        boost::filesystem::current_path());

    std::cout << "using config file url " << input_config_file << std::endl;

    output_path = output_path.parent_path();

    //load resonances
    ComPWA::Physics::DecayTree::DecayConfiguration decay_configuration;
    ComPWA::Physics::DecayTree::DecayXMLConfigReader xml_reader(
        decay_configuration);
    xml_reader.readConfig(input_config_file);

    ComPWA::Physics::DecayTree::DecayTreeFactory decay_tree_factory(
        decay_configuration);

    std::vector<ComPWA::Physics::DecayTree::DecayTree> decay_trees =
        decay_tree_factory.createDecayTrees();

    std::cout << "created " << decay_trees.size() << " decay trees from "
        << input_config_file << " config file!" << std::endl;

    if (decay_trees.size() > 0) {
      ComPWA::Physics::HelicityFormalism::TopologyAmplitudeFactory topology_amp_factory;

      Event dummy_event = topology_amp_factory.createDummyEvent(decay_trees[0]);

      std::vector<ComPWA::Physics::DecayTree::DecayNode> leaves =
          decay_trees[0].getLeaves();
      std::vector<ComPWA::Physics::HelicityFormalism::ParticleStateInfo> fs_particles;
      for (auto iter = leaves.begin(); iter != leaves.end(); ++iter) {
        fs_particles.push_back(iter->state_info_);
      }

      ComPWA::Physics::HelicityFormalism::FinalStateParticleCombinatorics fsp_combinatorics;
      fsp_combinatorics.init(fs_particles, dummy_event);

      ComPWA::Physics::HelicityFormalism::HelicityKinematics* kinematics =
          (ComPWA::Physics::HelicityFormalism::HelicityKinematics*) ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance();

      std::vector<ComPWA::Physics::HelicityFormalism::TwoBodyDecayTopology> decay_topologies =
          topology_amp_factory.generateDecayTopologies(decay_trees);

      kinematics->setDecayTopologies(decay_topologies);
      kinematics->init(fsp_combinatorics);

      // The helicity amplitude tree factory sorts the decay trees based
      // on their topology, and then creates the list of amplitude trees from the
      // topology grouped decay trees

      std::vector<ComPWA::Physics::HelicityFormalism::TopologyAmplitude> topology_amplitudes =
          topology_amp_factory.generateTopologyAmplitudes(decay_trees);

      std::cout << "created " << topology_amplitudes.size()
          << " topology amplitudes from the decay trees!" << std::endl;

      std::shared_ptr<ComPWA::Physics::HelicityFormalism::CoherentAmplitude> amp(
          new ComPWA::Physics::HelicityFormalism::CoherentAmplitude(
              topology_amplitudes));

      ParameterList par;
      amp->copyParameterList(par);    //perfect startvalues
      std::cout<<"number of parameters: "<<par.GetNDouble()<<std::endl;
      unsigned int counter_free_params(0);
      for (auto const& param : par.GetDoubleParameters()) {
        std::cout << param->GetName() << " " << param->IsFixed() << std::endl;
        if(!param->IsFixed())
          ++counter_free_params;
      }
      std::cout<<"free params: "<<counter_free_params<<std::endl;

      //create dummy final state event to initialized the kinematics class
      unsigned int dataSize = 10000;

      std::shared_ptr<Data> data(new RootReader());
      std::shared_ptr<Data> phsp(new RootReader());
      std::shared_ptr<Generator> gen(new RootGenerator());

      RunManager run(dataSize, amp, gen);
      run.setGenerator(gen);
      //run.setData(data);
      //run.generate(dataSize);
      run.setPhspSample(phsp);
      run.generatePhsp(dataSize * 10);
      //std::cout << "Data size: " << data->getNEvents() << std::endl;
      //data->writeData(output_path.string() + "/data.root", "events");
      phsp->writeData(output_path.string() + "/phspdata.root", "events");
    }
  }
  return 0;
}
