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
  if (argc > 2) {
    int seed = 1234;
    boost::filesystem::path input_config_file(argv[1]);
    unsigned int num_events(std::atoi(argv[2]));
    boost::filesystem::path output_data_filename_template("data.root");
    boost::filesystem::path output_phspdata_filename_template("phspdata.root");

    boost::filesystem::path output_directory = boost::filesystem::absolute(
        input_config_file, boost::filesystem::current_path());
    output_directory = output_directory.parent_path();

    string output_file_suffix("");
    if (argc > 3)
      output_data_filename_template = argv[3];
    if (argc > 4)
      output_phspdata_filename_template = argv[4];
    if (argc > 5)
      output_directory = argv[5];
    if (argc > 6)
      seed = atoi(argv[6]);
    if (argc > 7)
      output_file_suffix = argv[7];

    if (!boost::filesystem::is_regular_file(input_config_file)) {
      std::runtime_error(
          "Specified xml config file does not exist! Please specify a valid url...");
    }

    std::cout << "using config file url " << input_config_file.string()
        << std::endl;

    //load resonances
    ComPWA::Physics::DecayTree::DecayConfiguration decay_configuration;
    ComPWA::Physics::DecayTree::DecayXMLConfigReader xml_reader(
        decay_configuration);
    xml_reader.readConfig(input_config_file.string());

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
      std::cout << "number of parameters: " << par.GetNDouble() << std::endl;
      unsigned int counter_free_params(0);
      for (auto const& param : par.GetDoubleParameters()) {
        std::cout << param->GetName() << " " << param->IsFixed() << std::endl;
        if (!param->IsFixed())
          ++counter_free_params;
      }
      std::cout << "free params: " << counter_free_params << std::endl;

      //create dummy final state event to initialized the kinematics class

      // create output filenames
      string output_data_filename(output_data_filename_template.string());
      string output_phspdata_filename(
          output_phspdata_filename_template.string());
      bool want_data_generation(true);
      bool want_phsp_generation(true);
      if (output_data_filename_template.string().compare("") == 0) {
        // no data generation wanted
        want_data_generation = false;
      }
      else {
        if (output_file_suffix.compare("") != 0) {
          // add output suffix
          std::stringstream ss;
          ss << output_data_filename_template.stem().string() << "_"
              << output_file_suffix
              << output_data_filename_template.extension().string();
          output_data_filename = ss.str();
        }
      }
      if (output_phspdata_filename_template.string().compare("") == 0) {
        want_phsp_generation = false;
      }
      else {
        if (output_file_suffix.compare("") != 0) {
          // add output suffix
          std::stringstream ss;
          ss << output_phspdata_filename_template.stem().string() << "_"
              << output_file_suffix
              << output_phspdata_filename_template.extension().string();
          output_phspdata_filename = ss.str();
        }
      }

      std::shared_ptr<Data> data(new RootReader());
      std::shared_ptr<Data> phsp(new RootReader());
      std::shared_ptr<Generator> gen(new RootGenerator());

      RunManager run(amp, gen);
      run.setGenerator(gen);
      if (want_data_generation) {
        run.setData(data);
        run.generate(num_events);
        data->writeData(output_directory.string() + "/" + output_data_filename,
            "events");
      }
      if (want_phsp_generation) {
        run.setPhspSample(phsp);
        if (want_data_generation)
          run.generatePhsp(num_events * 10);
        else
          run.generatePhsp(num_events);
        phsp->writeData(
            output_directory.string() + "/" + output_phspdata_filename,
            "events");
      }
      //std::cout << "Data size: " << data->getNEvents() << std::endl;

    }
  }
  return 0;
}
