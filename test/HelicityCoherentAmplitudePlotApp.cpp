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

#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"

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

struct HelicityAngles {
  double theta_cms;
  double phi_cms;
  double theta_hel;
  double phi_hel;
};

std::pair<double, double> helicityVec(const Vector4<double>& cm,
    const Vector4<double>& mother, const Vector4<double>& daughter) {
  Vector4<double> mother_cm(mother);
  mother_cm.Boost(cm);

  Vector4<double> mother_cm_rot(mother_cm);
  mother_cm_rot.RotateZ(-mother_cm.Phi());
  mother_cm_rot.RotateY(-mother_cm.Theta());

  Vector4<double> result(daughter);
  result.Boost(cm);
  result.RotateZ(-mother_cm.Phi());
  result.RotateY(-mother_cm.Theta());

  result.Boost(mother_cm_rot);

  return std::make_pair(result.CosTheta(), result.Phi());
}

HelicityAngles phiDiff(const TLorentzVector &lv_p1, const TLorentzVector &lv_p2,
    const TLorentzVector &lv_p3) {
  // define particle products of the two body decay
  Vector4<double> p1_4vector(lv_p1.E(), lv_p1.Px(), lv_p1.Py(), lv_p1.Pz());
  Vector4<double> p2_4vector(lv_p2.E(), lv_p2.Px(), lv_p2.Py(), lv_p2.Pz());
  Vector4<double> p3_4vector(lv_p3.E(), lv_p3.Px(), lv_p3.Py(), lv_p3.Pz());

  // define the two body state
  Vector4<double> decaying_state(p1_4vector + p2_4vector);
  // define mother state
  Vector4<double> mother(decaying_state + p3_4vector);

  HelicityAngles hel_angles;

  p1_4vector.Boost(decaying_state);

  decaying_state.Boost(mother);
  hel_angles.theta_cms = decaying_state.CosTheta();
  hel_angles.phi_cms = decaying_state.Phi();

  p1_4vector.Rotate(decaying_state.Phi(), decaying_state.Theta(),
      -decaying_state.Phi());

  hel_angles.theta_hel = p1_4vector.CosTheta();
  hel_angles.phi_hel = p1_4vector.Phi();

  return hel_angles;


  /*HelicityAngles ha;

  Vector4<double> mother2(std::sqrt(mother.M2() + 1.0), 0.0, 0.0, 1.0);

  std::pair<double, double> vals = helicityVec(mother2, mother, decaying_state);
  ha.theta_cms = vals.first;
  ha.phi_cms = vals.second;

  vals = helicityVec(mother, decaying_state, p1_4vector);
  ha.theta_hel = vals.first;
  ha.phi_hel = vals.second;
  return ha;*/
}

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv) {
  ComPWA::Physics::HelicityFormalism::HelicityKinematics* kin =
      dynamic_cast<ComPWA::Physics::HelicityFormalism::HelicityKinematics*>(ComPWA::Physics::HelicityFormalism::HelicityKinematics::createInstance());

  //load resonances
  std::string input_config_file("Physics/HelicityAmplitude/JPSI_ypipi.xml");
  if (argc > 1) {
    input_config_file = argv[1];
  }

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
            topology_amplitudes, decay_configuration.getBackgroundPart()));

    TFile pawian_file("Data.root", "READ");

    TTree *pawian_tree = (TTree*) pawian_file.Get("Tree");

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector p3;
    TLorentzVector *pp1(&p1);
    TLorentzVector *pp2(&p2);
    TLorentzVector *pp3(&p3);

    double weight;
    pawian_tree->SetBranchAddress("weigth", &weight);
    pawian_tree->SetBranchAddress("P1_P4", &pp1);
    pawian_tree->SetBranchAddress("P2_P4", &pp2);
    pawian_tree->SetBranchAddress("P3_P4", &pp3);

    TH2D dalitz_phsp("dalitz_events", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
    TH2D dalitz_pawian("dalitz_pawian", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
    TH2D dalitz_compwa("dalitz_compwa", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
    TH2D dalitz_reldiff("dalitz_reldiff", "", 100, 0.0, 10.0, 100, 0.0, 10.0);
    TH1D costhetagamma("cms_cos_theta_gamma", "", 100, -1, 1);
    TH1D costhetapi01("cms_cos_theta_pi01", "", 100, -1, 1);
    TH1D costhetapi02("cms_cos_theta_pi02", "", 100, -1, 1);
    TH1D phidiff_pawian("phidiff_pawian", "", 100, -3.146, 3.146);
    TH1D phidiff_compwa("phidiff_compwa", "", 100, -3.146, 3.146);
    TH1D phi_cms("phi_cms", "", 100, -3.146, 3.146);
    TH1D phi_hel("phi_hel", "", 100, -3.146, 3.146);
    TH1D costheta_cms("costheta_cms", "", 100, -1, 1);
    TH1D costheta_hel("costheta_hel", "", 100, -1, 1);
    TH2D phicorr("phicorr", "", 100, -3.146, 3.146, 100, -3.146, 3.146);
    /*TH2D hel_cos_theta_vs_phi_diff_pawian("hel_cos_theta_vs_phi_diff_pawian", "", 100, -3.146, 3.146, 100, -1, 1);
     TH2D hel_cos_theta_vs_phi_diff_compwa("hel_cos_theta_vs_phi_diff_compwa", "", 100, -3.146, 3.146, 100, -1, 1);
     TH1D phidiff_pawian("phidiff_pawian", "", 100, -3.146, 3.146);
     TH1D phidiff_compwa("phidiff_compwa", "", 100, -3.146, 3.146);*/

    std::cout << pawian_tree->GetEntries() << std::endl;
    unsigned int entries(20000);
    //unsigned int entries = pawian_tree->GetEntries();
    TGraph graph_pawian(entries);
    TGraph graph(entries);
    ComPWA::dataPoint data_point;
    for (unsigned int i = 0; i < entries; ++i) {
      pawian_tree->GetEntry(i);

      /*TLorentzVector cms = p1 + p2 + p3;
       p1.Boost(-cms.BoostVector());
       p2.Boost(-cms.BoostVector());
       p3.Boost(-cms.BoostVector());*/

      /* cms.Print();
       cms.Boost(-cms.BoostVector());
       std::cout<<"boost cms"<<std::endl;
       cms.Print();*/

      /* if (p1.M() > 1e-10 || p1.M() < 0.0) {
       p1.SetE(p1.Vect().Mag());
       if (p1.M() < 0.0)
       p1.SetE(p1.E() + std::numeric_limits<double>::epsilon());
       }*/


      ComPWA::Particle gamma(p1.Px(), p1.Py(), p1.Pz(), p1.E(), 22);
      ComPWA::Particle pi01(p2.Px(), p2.Py(), p2.Pz(), p2.E(), 111);
      ComPWA::Particle pi02(p3.Px(), p3.Py(), p3.Pz(), p3.E(), 11122);

      unsigned int pi0counter(0);
      for (unsigned int i = 0; i < dummy_event.getNParticles(); ++i) {
        //std::cout<<dummy_event.getParticle(i).pid<<std::endl;
        if (dummy_event.getParticle(i).pid == 22) {
          dummy_event.setParticleAt(gamma, i);
          //std::cout<<"asdf1\n";
        }
        else if (dummy_event.getParticle(i).pid == 111) {
          if (pi0counter == 0)
            dummy_event.setParticleAt(pi01, i);
          if (pi0counter == 1)
            dummy_event.setParticleAt(pi02, i);
          ++pi0counter;
          //std::cout<<"asdf2\n";
        }
        else if (dummy_event.getParticle(i).pid == 11122) {
          dummy_event.setParticleAt(pi02, i);
          //std::cout<<"asdf3\n";
        }
        else {
          std::cout << "wtf\n";
        }
      }

      dataPoint data_point(dummy_event);

      ParameterList list;
      list = amp->intensity(data_point);    //unfortunatly not thread safe
      double myval = *list.GetDoubleParameter(0);
      graph_pawian.SetPoint(i, i, weight);
      graph.SetPoint(i, i, myval);
      /*std::cout
       << parlist.GetDoubleParameter("coherent_amp")->GetValue() / weight
       << std::endl;*/

      TLorentzVector lv1(p1 + p2);
      TLorentzVector lv2(p2 + p3);

      dalitz_phsp.Fill(lv2.M2(), lv1.M2());
      dalitz_pawian.Fill(lv2.M2(), lv1.M2(), weight);
      dalitz_compwa.Fill(lv2.M2(), lv1.M2(), myval);
      dalitz_reldiff.Fill(lv2.M2(), lv1.M2(),
          (myval - weight) / (weight + myval));

      costhetagamma.Fill(p1.CosTheta(), (myval - weight) / (weight + myval));
      costhetapi01.Fill(p2.CosTheta(), (myval - weight) / (weight + myval));
      costhetapi02.Fill(p3.CosTheta(), (myval - weight) / (weight + myval));

      HelicityAngles ha = phiDiff(p1, p2, p3);

      double diff(ha.phi_cms - ha.phi_hel);
      if (diff < -TMath::Pi())
        diff += 2.0*TMath::Pi();
      else if (diff > TMath::Pi())
        diff -= 2.0*TMath::Pi();
      phidiff_pawian.Fill(diff, weight);
      phidiff_compwa.Fill(diff, myval);

      phicorr.Fill(ha.phi_cms, ha.phi_hel, weight);

      phi_cms.Fill(ha.phi_cms, myval);
      phi_hel.Fill(ha.phi_hel, myval);
      costheta_cms.Fill(ha.theta_cms, myval);
      costheta_hel.Fill(ha.theta_hel, myval);
    }
    TFile f("comparison.root", "RECREATE");
    graph.Write("compwa_weights");
    graph_pawian.Write("pawian_weights");
    dalitz_phsp.Write();
    dalitz_pawian.Write();
    dalitz_compwa.Write();
    dalitz_reldiff.Write();
    costhetagamma.Write();
    costhetapi01.Write();
    costhetapi02.Write();
    phidiff_pawian.Write();
    phidiff_compwa.Write();

    phicorr.Write();
    phi_cms.Write();
    phi_hel.Write();
    costheta_cms.Write();
    costheta_hel.Write();

  }
  return 0;
}
