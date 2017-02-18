//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------

#include <thread>
#include <future>
#include <functional>

#include "Core/DataPointStorage.hpp"
#include "Physics/HelicityAmplitude/CoherentAmplitude.hpp"
#include "Physics/HelicityAmplitude/HelicityKinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

CoherentAmplitude::CoherentAmplitude(
    const std::vector<TopologyAmplitude>& amplitude_trees,
    const boost::property_tree::ptree& background_part) :
    topology_amplitudes_(amplitude_trees), wasMaxAmplitudeValueCalculated_(
    false), maxAmplitudeValue_(0.0), background_part_(background_part) {

  init();
}

CoherentAmplitude::~CoherentAmplitude() {
}

void CoherentAmplitude::init() {
// TODO: I have to cast the kinematics instance to a helicity kinematics.
// This is pretty bad practice... I really don't get the point behind the
// singleton kinematics class...
// If you do not instantiate a different type of kinematics class beforehand
// you are fine, but otherwise we are in deep shit.
  HelicityKinematics* kinematics = (HelicityKinematics*) Kinematics::instance();

// initialize the kinematics class first
//  kinematics->init(event);
// now get the index lists that tell the topology amplitudes which
// data point variables to use in their evaluation
  data_point_index_lists_ =
      kinematics->getTopologyAmplitudeDataPointIndexLists();

  // and initialize parameter at position 0 in the result objects as result
  result_value_.reset(new DoubleParameter("coherent_amp"));
  result.AddParameter(result_value_);

  // ok this part is only for the function tree
  // now comes the actual hard part. constructing the amplitude from the parts

  // so the plan is to have evaluation lists as leaves. then those calculate two body decay amps
  // and from those individual building blocks we have to build the coherent/incoherent sum

  // ok we need to create two trees one for phsp and one for data
  // we give phasespace storage index 0 and data storage index 1
  // these trees are completely independent and are constructed independently
  unsigned int storage_index(0);
  constructSequentialDecayTreeNodes(storage_index);
  constructCoherentAmpTree(storage_index);
  storage_index = 1;
  clearMaps();
  constructSequentialDecayTreeNodes(storage_index);
  constructCoherentAmpTree(storage_index);
  //std::cout<<tree_[1]->print()<<std::endl;
}

void CoherentAmplitude::constructSequentialDecayTreeNodes(
    unsigned int storage_index) {
  // ok we will just go through the topology amplitudes and build up the basic building blocks
  unsigned int topology_index(0);
  for (auto const& decay_topology_amp : topology_amplitudes_) {
    // for each final state combination for this topology there exists a index list
    const std::vector<IndexList> evaluation_lists(
        data_point_index_lists_[topology_index]);
    ++topology_index;

    unsigned int seq_decay_amp_counter(0);
    for (auto const& sequential_decay : decay_topology_amp.getSequentialDecayList()) {
      std::vector<std::shared_ptr<TreeNode> > single_combinatoric_sequential_decay_nodes;

      std::set<unsigned int> seq_decay_info;

      for (unsigned int i = 0; i < evaluation_lists.size(); ++i) {

        // just connect all these sequential decay tree nodes with the appropriate
        unsigned int two_body_decay_index(0);

        std::vector<std::shared_ptr<TreeNode> > temp_two_body_decay_list;
        std::string seq_decay_node_name;

        for (auto const& full_two_body_decay_amp : sequential_decay.decay_amplitude) {
          std::stringstream eval_leaf_name;
          eval_leaf_name << "topology_" << topology_index << "_combination_"
              << i << "tbd_index_" << two_body_decay_index;
          IndexList temp_index_list;
          temp_index_list.push_back(evaluation_lists[i][two_body_decay_index]);
          std::shared_ptr<MultiUnsignedInteger> init_value(
              new MultiUnsignedInteger(eval_leaf_name.str(), temp_index_list));

          std::shared_ptr<TreeNode> eval_list_leaf(
              new TreeNode(eval_leaf_name.str(), init_value, nullptr, nullptr));

          std::string basename(full_two_body_decay_amp.name);
          basename = basename + "_" + eval_leaf_name.str();

          seq_decay_node_name = basename + "-&-" + seq_decay_node_name;

          auto tbd_psi = full_two_body_decay_amp.decay_spin_info_;
          auto result = std::find(particles_.begin(), particles_.end(),
              tbd_psi.first);
          unsigned int mother_index(result - particles_.begin());
          if (result == particles_.end()) {
            particles_.push_back(tbd_psi.first);
            mother_index = particles_.size() - 1;
          }
          result = std::find(particles_.begin(), particles_.end(),
              tbd_psi.second.first);
          unsigned int d1_index(result - particles_.begin());
          if (result == particles_.end()) {
            particles_.push_back(tbd_psi.second.first);
            d1_index = particles_.size() - 1;
          }
          result = std::find(particles_.begin(), particles_.end(),
              tbd_psi.second.second);
          unsigned int d2_index(result - particles_.begin());
          if (result == particles_.end()) {
            particles_.push_back(tbd_psi.second.second);
            d2_index = particles_.size() - 1;
          }

          // only do this once for the first evaluation list
          if (i == 0) {
            // ok just add it to the index decay tree, the indices should be unique
            seq_decay_info.insert(mother_index);
            seq_decay_info.insert(d1_index);
            seq_decay_info.insert(d2_index);
            /* seq_decay_info.unique_id_decay_tree_[mother_index].push_back(
             d1_index);
             seq_decay_info.unique_id_decay_tree_[mother_index].push_back(
             d2_index);*/
          }

          std::shared_ptr<MultiComplex> full_tbd_amp(
              new MultiComplex(basename + "_full_values",
                  std::vector<std::complex<double> >()));
          auto full_two_body_decay = full_two_body_decay_nodes_.insert(
              std::make_pair(basename + "_full",
                  std::shared_ptr<TreeNode>(
                      new TreeNode(basename + "_full", full_tbd_amp,
                          std::shared_ptr<MultAll>(
                              new MultAll(ParType::MCOMPLEX)), nullptr))));

          if (full_two_body_decay.second) {
            /*std::shared_ptr<TreeNode> strength(
             new TreeNode(basename + "_mag",
             full_two_body_decay_amp.strength_, nullptr, nullptr));

             std::shared_ptr<TreeNode> phase(
             new TreeNode(basename + "_phase",
             full_two_body_decay_amp.phase_, nullptr, nullptr));*/

            parameters_.AddParameter(full_two_body_decay_amp.strength_);
            parameters_.AddParameter(full_two_body_decay_amp.phase_);

            std::shared_ptr<ComplexParameter> temp_str_phase(
                new ComplexParameter(basename + "_str-phase_value",
                    std::complex<double>(0, 0)));
            std::shared_ptr<TreeNode> strength_and_phase(
                new TreeNode(basename + "_str-phase", temp_str_phase,
                    std::shared_ptr<Complexify>(
                        new Complexify(ParType::COMPLEX)), nullptr));

            //strength_and_phase->addChild(strength);
            //strength_and_phase->addChild(phase);
            tree_leaves_[strength_and_phase->getName()].push_back(
                std::make_pair(basename + "_mag",
                    full_two_body_decay_amp.strength_));
            tree_leaves_[strength_and_phase->getName()].push_back(
                std::make_pair(basename + "_phase",
                    full_two_body_decay_amp.phase_));

            full_two_body_decay.first->second->addChild(strength_and_phase);

            // create strategies
            auto tbdas = angular_part_strategies_.insert(
                std::make_pair(basename,
                    std::shared_ptr<TwoBodyDecayAngularStrategy>(
                        new TwoBodyDecayAngularStrategy(
                            full_two_body_decay_amp.angular_part_,
                            storage_index))));
            auto dfs = dynamical_part_strategies_.insert(
                std::make_pair(basename,
                    std::shared_ptr<
                        DynamicalFunctions::DynamicalFunctionStrategy>(
                        new DynamicalFunctions::DynamicalFunctionStrategy(
                            full_two_body_decay_amp.dynamical_part_,
                            storage_index))));

            std::shared_ptr<MultiComplex> initial_value(
                new MultiComplex("ang_part_values",
                    std::vector<std::complex<double> >()));
            auto two_body_decay_angular_part = angular_part_nodes_.insert(
                std::make_pair(basename,
                    std::shared_ptr<TreeNode>(
                        new TreeNode(basename + "_angular", initial_value,
                            tbdas.first->second, nullptr))));
            full_two_body_decay.first->second->addChild(
                two_body_decay_angular_part.first->second);

            std::shared_ptr<MultiComplex> initial_value2(
                new MultiComplex("dyn_part_values",
                    std::vector<std::complex<double> >()));
            auto two_body_decay_dynamical_part = dynamical_part_nodes_.insert(
                std::make_pair(basename,
                    std::shared_ptr<TreeNode>(
                        new TreeNode(basename + "_dynamic", initial_value2,
                            dfs.first->second, nullptr))));
            full_two_body_decay.first->second->addChild(
                two_body_decay_dynamical_part.first->second);

            if (two_body_decay_dynamical_part.second) {
              // we need to add the parameters of the dynamical function in order to trigger recalculation
              // they are not actually passed to the function because they have a copy of the same shared ptr
              // this concept has be be reworked... does not seem the best way
              /*std::vector<std::shared_ptr<TreeNode> > dynamic_part_leaves =
               createLeaves(
               full_two_body_decay_amp.dynamical_part_->getParameterList());
               for (auto leaf : dynamic_part_leaves) {
               two_body_decay_dynamical_part.first->second->addChild(leaf);
               }*/
              std::vector<std::shared_ptr<AbsParameter> > dynamic_part_leaves =
                  createLeaves(
                      full_two_body_decay_amp.dynamical_part_->getParameterList());
              auto& leaf_vec =
                  tree_leaves_[two_body_decay_dynamical_part.first->second->getName()];
              for (auto leaf : dynamic_part_leaves) {
                leaf_vec.push_back(std::make_pair(leaf->GetName(), leaf));
              }

              // add the evaluation lists so that final state combinatorics is taken care of as well
              two_body_decay_dynamical_part.first->second->addChild(
                  eval_list_leaf);
            }

            if (two_body_decay_angular_part.second) {
              // add the evaluation lists so that final state combinatorics is taken care of as well
              two_body_decay_angular_part.first->second->addChild(
                  eval_list_leaf);
            }
          }

          temp_two_body_decay_list.push_back(full_two_body_decay.first->second);

          ++two_body_decay_index;
        }

        // create tree node for sequential decay if it is not in storage yet
        std::shared_ptr<MultiComplex> seq_tbd_amp(
            new MultiComplex(seq_decay_node_name + "_values",
                std::vector<std::complex<double> >()));
        std::shared_ptr<TreeNode> seq_decay_node(
            new TreeNode(seq_decay_node_name, seq_tbd_amp,
                std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX)),
                nullptr));

        for (auto const& tbd_amp : temp_two_body_decay_list) {
          seq_decay_node->addChild(tbd_amp);
        }

        single_combinatoric_sequential_decay_nodes.push_back(seq_decay_node);
      }

      std::shared_ptr<TreeNode> combinatorics_seq_decay_amp_node;
      if (single_combinatoric_sequential_decay_nodes.size() == 1) {
        // just use this node
        combinatorics_seq_decay_amp_node =
            single_combinatoric_sequential_decay_nodes[0];
      }
      else {
        // add all nodes coherently
        std::stringstream node_label;

        for (auto node : single_combinatoric_sequential_decay_nodes) {
          node_label << "_" << node->getName();
        }

        std::shared_ptr<MultiComplex> combi_seq_amp_dummy_val(
            new MultiComplex("combinatorics" + node_label.str() + "_values",
                std::vector<std::complex<double> >()));
        combinatorics_seq_decay_amp_node.reset(
            new TreeNode("combinatorics" + node_label.str(),
                combi_seq_amp_dummy_val,
                std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX)),
                nullptr));

        for (auto node : single_combinatoric_sequential_decay_nodes) {
          combinatorics_seq_decay_amp_node->addChild(node);
        }

      }


      ++seq_decay_amp_counter;
      std::stringstream ss;
      ss<<seq_decay_amp_counter<<"_parity_factor";

      std::shared_ptr<DoubleParameter> parity_factor_value(
                      new DoubleParameter(ss.str()+"_value",
                          sequential_decay.factor));
      std::shared_ptr<TreeNode> parity_factor(
                      new TreeNode(ss.str(), parity_factor_value,
                          nullptr, nullptr));

      std::shared_ptr<MultiComplex> parity_corrected_amp_dummy_val(
                  new MultiComplex("parity_corrected_" + combinatorics_seq_decay_amp_node->getName() + "_values",
                      std::vector<std::complex<double> >()));
      std::shared_ptr<TreeNode> parity_corrected_amp(
                      new TreeNode("parity_corrected_" + combinatorics_seq_decay_amp_node->getName(), parity_corrected_amp_dummy_val,
                          std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX)), nullptr));
      // add the parity factor
      parity_corrected_amp->addChild(parity_factor);
      parity_corrected_amp->addChild(combinatorics_seq_decay_amp_node);

      // determine top node and final state of sequential decay
      /*for (auto const& decay_node : seq_decay_info.unique_id_decay_tree_) {
       bool is_top_node(true);
       bool is_final_state(true);
       for (auto const& test_decay_node : seq_decay_info.unique_id_decay_tree_) {
       if (std::find(test_decay_node.second.begin(),
       test_decay_node.second.end(), decay_node.first)
       != test_decay_node.second.end()) {
       is_top_node = false;
       break;
       }
       }
       if (is_top_node)
       seq_decay_info.top_node_ = decay_node.first;

       for (auto possible_final_state_unique_id : decay_node.second) {
       if (seq_decay_info.unique_id_decay_tree_.find(
       possible_final_state_unique_id)
       == seq_decay_info.unique_id_decay_tree_.end()) {
       seq_decay_info.final_state_.push_back(possible_final_state_unique_id);
       }
       }
       }*/

      auto insert_success = sequential_decay_amplitudes_map_.insert(
          std::make_pair(seq_decay_info,
              sequential_decay_amplitudes_vec_.size()));
      sequential_decay_amplitudes_vec_.push_back(
          parity_corrected_amp);

      // ok in principle we should never fail to create a new two body decay
      // safety check!
      if (!insert_success.second)
        std::runtime_error(
            "CoherentAmplitude::init(): error in the sequential decay tree node construction!");
    }
  }

  // prepare a phase space term
  boost::optional<boost::property_tree::ptree&> coherent_background_opt =
      background_part_.get_child_optional("coherent");
  if (coherent_background_opt.is_initialized()) {
    auto temp_coherent_part_config = background_part_.get_child("coherent");
    std::shared_ptr<DoubleParameter> coherent_phsp_mag_value(
        new DoubleParameter("coherent_phsp_mag_value",
            temp_coherent_part_config.get_child("strength").get<double>(
                "value")));
    coherent_phsp_mag_value->FixParameter(
        temp_coherent_part_config.get_child("strength").get<bool>("fix"));

    std::shared_ptr<DoubleParameter> coherent_phsp_phase_value(
        new DoubleParameter("coherent_phsp_phase_value",
            temp_coherent_part_config.get_child("phase").get<double>("value")));
    coherent_phsp_phase_value->FixParameter(
        temp_coherent_part_config.get_child("phase").get<bool>("fix"));

    parameters_.AddParameter(coherent_phsp_mag_value);
    parameters_.AddParameter(coherent_phsp_phase_value);

    auto result_coherent_phsp_mag_value = dynamical_parameter_nodes_.insert(
        std::make_pair(coherent_phsp_mag_value->GetName(),
            coherent_phsp_mag_value));
    auto result_coherent_phsp_phase_value = dynamical_parameter_nodes_.insert(
        std::make_pair(coherent_phsp_phase_value->GetName(),
            coherent_phsp_phase_value));

    std::shared_ptr<ComplexParameter> coherent_phsp_value(
        new ComplexParameter("coherent_phsp_value",
            std::complex<double>(0, 0)));
    coherent_phsp_.reset(
        new TreeNode("coherent_phsp", coherent_phsp_value,
            std::shared_ptr<Complexify>(new Complexify(ParType::COMPLEX)),
            nullptr));
    std::shared_ptr<AbsParameter> temp1 =
        result_coherent_phsp_mag_value.first->second;
    std::shared_ptr<AbsParameter> temp2 =
        result_coherent_phsp_phase_value.first->second;
    tree_leaves_[coherent_phsp_->getName()].push_back(
        std::make_pair("coherent_phsp_mag_value", temp1));
    tree_leaves_[coherent_phsp_->getName()].push_back(
        std::make_pair("coherent_phsp_phase_value", temp2));
  }
}

void CoherentAmplitude::constructCoherentAmpTree(unsigned int storage_index) {
  // now we build the intensity from that

  /* ok now when we construct the coherent amplitude we have to go through the sequential decay list
   and find suitable partners for the summation.
   so we need some kind of search algorithm that goes through the map and look for sequential decay trees
   that have the coherent index the same or different*/

  /* std::vector<std::pair<unsigned int, IndexList> > coherent_sum_pairings;

   // create the different pairings based on the coherency of the individual parts
   for (auto seq_amp_node : sequential_decay_amplitudes_map_) {
   IndexList partners = getListOfCoherentPartners(seq_amp_node);
   coherent_sum_pairings.push_back(
   std::make_pair(seq_amp_node.second, partners));
   }

   // construct the actual full tree amplitude
   std::shared_ptr<MultiComplex> coh_amp_val(
   new MultiComplex("coherent_sum_value",
   std::vector<std::complex<double> >()));
   std::shared_ptr<TreeNode> coherent_amp(
   new TreeNode("coherent_sum", coh_amp_val,
   std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX)), nullptr));

   unsigned int node_label_index(0);
   for (auto const& pairing : coherent_sum_pairings) {
   ++node_label_index;
   std::stringstream node_label;
   node_label << node_label_index;
   std::shared_ptr<MultiComplex> coh_amp_dummy_val(
   new MultiComplex("coherent_amp_part1_" + node_label.str() + "_values",
   std::vector<std::complex<double> >()));
   std::shared_ptr<TreeNode> coherent_amp_node_part1(
   new TreeNode("coherent_amp_part1_" + node_label.str(),
   coh_amp_dummy_val,
   std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX)), nullptr));
   for (auto const& node_to_conjugate : pairing.second) {
   coherent_amp_node_part1->addChild(
   sequential_decay_amplitudes_vec_[node_to_conjugate]);
   }

   std::shared_ptr<MultiComplex> coh_amp_dummy_val2(
   new MultiComplex(
   "coherent_amp_part1_conj_" + node_label.str() + "_values",
   std::vector<std::complex<double> >()));
   std::shared_ptr<TreeNode> coherent_amp_node_part1_conj(
   new TreeNode("coherent_amp_part1_conj_" + node_label.str(),
   coh_amp_dummy_val2,
   std::shared_ptr<ComplexConjugate>(
   new ComplexConjugate(ParType::MCOMPLEX)), nullptr));
   coherent_amp_node_part1_conj->addChild(coherent_amp_node_part1);

   std::shared_ptr<MultiComplex> coh_amp_dummy_val3(
   new MultiComplex("coherent_amp_" + node_label.str() + "_values",
   std::vector<std::complex<double> >()));
   std::shared_ptr<TreeNode> coherent_amp_node(
   new TreeNode("coherent_amp_" + node_label.str(), coh_amp_dummy_val3,
   std::shared_ptr<MultAll>(new MultAll(ParType::MCOMPLEX)), nullptr));
   coherent_amp_node->addChild(coherent_amp_node_part1_conj);
   coherent_amp_node->addChild(
   sequential_decay_amplitudes_vec_[pairing.first]);

   coherent_amp->addChild(coherent_amp_node);*/

  std::vector<IndexList> coherent_sum_parts;

  // create the different pairings based on the coherency of the individual parts
  for (auto seq_amp_node : sequential_decay_amplitudes_map_) {
    // getListofCoherentPartners gets index of all amplitudes that have same coherent topology
    // that also includes the amplitude itself that is used as reference because the loop goes
    // over all amplitudes again... make that nicer after handing in thesis TODO
    IndexList amplitude_group = getListOfCoherentPartners(seq_amp_node);
    auto result = std::find_if(coherent_sum_parts.begin(),
        coherent_sum_parts.end(), isIndexListContentEqual(amplitude_group));
    if (result == coherent_sum_parts.end())
      coherent_sum_parts.push_back(amplitude_group);
  }

  // construct the actual full tree amplitude
  std::shared_ptr<MultiDouble> coh_amp_val(
      new MultiDouble("coherent_sum_value", std::vector<double>()));
  std::shared_ptr<TreeNode> coherent_amp(
      new TreeNode("coherent_sum", coh_amp_val,
          std::shared_ptr<AddAll>(new AddAll(ParType::MDOUBLE)), nullptr));

  unsigned int node_label_index(0);
  for (auto const& amplitude_group : coherent_sum_parts) {
    ++node_label_index;
    std::stringstream node_label;
    node_label << node_label_index;
    std::shared_ptr<MultiComplex> coh_amp_dummy_val(
        new MultiComplex("coherent_amp_nodesum_" + node_label.str() + "_values",
            std::vector<std::complex<double> >()));
    std::shared_ptr<TreeNode> coherent_amp_nodesum(
        new TreeNode("coherent_amp_nodesum_" + node_label.str(),
            coh_amp_dummy_val,
            std::shared_ptr<AddAll>(new AddAll(ParType::MCOMPLEX)), nullptr));
    for (auto const& node : amplitude_group) {
      coherent_amp_nodesum->addChild(sequential_decay_amplitudes_vec_[node]);
    }

    boost::optional<boost::property_tree::ptree&> coherent_background_opt =
        background_part_.get_child_optional("coherent");
    if (coherent_background_opt.is_initialized()) {
      // add coherent background
      coherent_amp_nodesum->addChild(coherent_phsp_);
    }

    std::shared_ptr<MultiDouble> coh_amp_dummy_val2(
        new MultiDouble("coherent_amp_part_" + node_label.str() + "_values",
            std::vector<double>()));
    std::shared_ptr<TreeNode> coherent_amp_part(
        new TreeNode("coherent_amp_part_" + node_label.str(),
            coh_amp_dummy_val2,
            std::shared_ptr<AbsSquare>(new AbsSquare(ParType::MDOUBLE)),
            nullptr));
    coherent_amp_part->addChild(coherent_amp_nodesum);
    coherent_amp->addChild(coherent_amp_part);
  }

  /*  // construct the actual full tree amplitude
   std::shared_ptr<MultiDouble> real_coh_amp_val(
   new MultiDouble("real_coherent_sum_value", std::vector<double>()));
   std::shared_ptr<TreeNode> real_coherent_amp(
   new TreeNode("real_coherent_sum", real_coh_amp_val,
   std::shared_ptr<Real>(new Real(ParType::MCOMPLEX)), nullptr));
   real_coherent_amp->addChild(coherent_amp);*/

  tree_[storage_index] = std::shared_ptr<FunctionTree>(new FunctionTree());
  tree_[storage_index]->addHead(coherent_amp);

  // ok create leafs
  for (auto const& parent_leaves : tree_leaves_) {
    for (auto leaf : parent_leaves.second) {
      tree_[storage_index]->createLeaf(leaf.first, leaf.second,
          parent_leaves.first);
    }
  }

  tree_[storage_index]->fixLinks();
  tree_[storage_index]->sanityCheck();
  //tree_[storage_index]->forceRecalculate();
  //std::cout << tree_[storage_index]->print() << std::endl;
}

void CoherentAmplitude::clearMaps() {
// clear all of these containers
  particles_.clear();
  sequential_decay_amplitudes_map_.clear();
  sequential_decay_amplitudes_vec_.clear();

  angular_part_nodes_.clear();
  dynamical_part_nodes_.clear();
  //dynamical_parameter_nodes_.clear();
  full_two_body_decay_nodes_.clear();
  tree_leaves_.clear();

  angular_part_strategies_.clear();
  dynamical_part_strategies_.clear();
}

std::vector<IndexList> CoherentAmplitude::convertIndexLists(
    const std::vector<IndexList> evaluation_lists) const {
  std::vector<IndexList> converted_list;
  if (evaluation_lists.size() > 0) {
    converted_list.resize(evaluation_lists[0].size());
    for (auto const& eval_list : evaluation_lists) {
      for (unsigned int i = 0; i < eval_list.size(); ++i) {
        converted_list[i].push_back(eval_list[i]);
      }
    }
  }
  return converted_list;
}

std::vector<std::shared_ptr<AbsParameter> > CoherentAmplitude::createLeaves(
    const ParameterList& parameter_list) {
  std::vector<std::shared_ptr<AbsParameter> > leaves;

  auto const& parameters = parameter_list.GetDoubleParameters();
  for (auto const& param : parameters) {
    parameters_.AddParameter(param);
    auto result = dynamical_parameter_nodes_.insert(
        std::make_pair(param->GetName(), param));
    /*auto result = dynamical_parameter_nodes_.insert(
     std::make_pair(param->GetName(),
     std::shared_ptr<TreeNode>(
     new TreeNode(param->GetName(), param, nullptr, nullptr))));*/

    leaves.push_back(result.first->second);
  }

  return leaves;
}

IndexList CoherentAmplitude::getListOfCoherentPartners(
    const std::pair<std::set<unsigned int>, unsigned int>& seq_amp) const {
  IndexList coherent_partners;

  //find all elements with appropriate quantum numbers
  /*for (auto const& seq_amp_element : sequential_decay_amplitudes_map_) {
   if (isCoherentPartner(
   std::make_pair(seq_amp_element.first, seq_amp_element.first.top_node_),
   std::make_pair(seq_amp.first, seq_amp.first.top_node_)))
   coherent_partners.push_back(seq_amp_element.second);
   }
   return coherent_partners;*/

  // ok first remove all indices from the list that correspond to coherent particles
  std::set<unsigned int> incoherent_particles = removeCoherentParticleIndices(
      seq_amp.first);

  // loop over all sequential amplitudes
  for (auto const& seq_amp_element : sequential_decay_amplitudes_map_) {
    bool is_valid(true);
    // now we want to remove all coherent particles from the list
    std::set<unsigned int> incoherent_particles_other_seq_amp =
        removeCoherentParticleIndices(seq_amp_element.first);

    // then check if all incoherent particles are equal and thats it...
    if (incoherent_particles.size()
        == incoherent_particles_other_seq_amp.size()) {
      for (auto particle_index : incoherent_particles_other_seq_amp) {
        auto ref_partner = incoherent_particles.find(particle_index);
        if (ref_partner == incoherent_particles.end()) {
          is_valid = false;
          break;
        }
      }
      if (is_valid) {
        coherent_partners.push_back(seq_amp_element.second);
      }
    }
  }

  return coherent_partners;
}

std::set<unsigned int> CoherentAmplitude::removeCoherentParticleIndices(
    const std::set<unsigned int>& particles_indices) const {
  std::set<unsigned int> incoherent_particles;
  for (auto particle_index : particles_indices) {
    if (!particles_[particle_index].coherent) {
      incoherent_particles.insert(particle_index);
    }
  }
  return incoherent_particles;
}

/*bool CoherentAmplitude::isCoherentPartner(
 const std::pair<std::set<unsigned int>, unsigned int>& seq_amp_partner,
 const std::pair<std::set<unsigned int>, unsigned int>& seq_amp_ref) const {

 if (!compareCoherentParticleStates(particles_[seq_amp_partner.second],
 particles_[seq_amp_ref.second])) {
 return false;
 }

 auto result_ref = seq_amp_ref.first.unique_id_decay_tree_.find(
 seq_amp_ref.second);
 auto result_partner = seq_amp_partner.first.unique_id_decay_tree_.find(
 seq_amp_partner.second);
 if (result_ref != seq_amp_ref.first.unique_id_decay_tree_.end()
 && result_partner != seq_amp_partner.first.unique_id_decay_tree_.end()) {
 IndexList ref_index_list(result_ref->second);
 for (auto partner_index : result_partner->second) {
 auto found =
 std::find_if(ref_index_list.begin(), ref_index_list.end(),
 [&](unsigned int ref_index) {return particles_[partner_index].unique_id_ == particles_[ref_index].unique_id_;});
 if (found != ref_index_list.end()) {
 bool is_coherent_partner = isCoherentPartner(
 std::make_pair(seq_amp_partner.first, partner_index),
 std::make_pair(seq_amp_ref.first, *found));
 if (!is_coherent_partner)
 return false;
 ref_index_list.erase(found);
 }
 else {
 // ok we did not find a particle with same unique id...
 // but maybe we can find one with same pid
 auto found_same_pid =
 std::find_if(ref_index_list.begin(), ref_index_list.end(),
 [&](unsigned int ref_index) {return particles_[partner_index].pid_information_.particle_id_ == particles_[ref_index].pid_information_.particle_id_;});
 if (found_same_pid != ref_index_list.end()) {
 bool is_coherent_partner = isCoherentPartner(
 std::make_pair(seq_amp_partner.first, partner_index),
 std::make_pair(seq_amp_ref.first, *found_same_pid));
 if (!is_coherent_partner)
 return false;
 ref_index_list.erase(found_same_pid);
 }
 else {
 BOOST_LOG_TRIVIAL(debug)<<"CoherentAmplitude::isCoherentPartner() amplitudes are not matching from the topology";
 return false;
 }
 }
 }
 }

 return true;
 }*/

bool CoherentAmplitude::compareCoherentParticleStates(
    const ParticleStateInfo& ps, const ParticleStateInfo& ref) const {
  //if (ps.unique_id_ != ref.unique_id_)
  //  return false;
  //if (ps.pid_information_ != ref.pid_information_)
  //  return false;
  if (ps.coherent != ref.coherent)
    return false;
  else if (!ps.coherent) {
    if (ps.spin_information_ != ref.spin_information_)
      return false;
  }
  return true;
}

const double CoherentAmplitude::integral() {

}
const double CoherentAmplitude::integral(ParameterList& par) {

}
const double CoherentAmplitude::normalization() {

}
const double CoherentAmplitude::normalization(ParameterList& par) {

}
double CoherentAmplitude::getMaxVal(ParameterList& par,
    std::shared_ptr<Generator> gen) {
  setParameterList(par);
  return getMaxVal(gen);
}

double CoherentAmplitude::getMaxVal(std::shared_ptr<Generator> gen) {
  if (!wasMaxAmplitudeValueCalculated_) {
    unsigned int evaluations(10000);
    BOOST_LOG_TRIVIAL(info)<<"CoherentAmplitude::calcMaxVal() calculating amplitude max value with "
    <<evaluations<< " events...";
    HelicityKinematics* kin =
        dynamic_cast<HelicityKinematics*>(Kinematics::instance());

    /*std::vector<std::future<unsigned int> > future_list;
     std::vector<std::thread> thread_list;

     std::packaged_task<unsigned int()> task(
     std::bind(&IntegralStrategyGSL2D::determineOptimalCallNumber,
     integral_strategy.get(), divergence_model.get(),
     std::cref(int_ranges[i]), 1e-4));
     future_list.push_back(task.get_future());
     thread_list.push_back(std::thread(std::move(task)));

     // wait for futures and compute maximum number of calls

     for (auto& future : future_list) {
     unsigned int temp_calls = future.get();
     if (temp_calls > calls)
     calls = temp_calls;
     }

     // join all threads
     for (auto& thread : thread_list) {
     if (thread.joinable())
     thread.join();
     }*/

    Event event;
    //Event max_event;
    //dataPoint max_point;
    double maxVal = 0;
    for (unsigned int i = 0; i < evaluations; i++) {
      //create event
      gen->generate(event);
      dataPoint point;
      kin->eventToDataPoint(event, point);

      if (!kin->isWithinPhsp(point)) {
        if (i > 0)
          i--;
        continue;
      }    //only integrate over phase space

      intensity(point);
      double intens = result_value_->GetValue();

      if (intens > maxVal) {
        maxVal = intens;
        // if this event is the new max then we are maybe close to the actual maximum
        // so vary this event a little and see if that brings us even closer
        /*max_point = point;
         max_event = event;

         std::vector<Particle> aasdf;
         for (unsigned int i = 0; i < max_event.getNParticles(); ++i) {
         auto particle = max_event.getParticle(i);
         std::cout << particle.pid << ": (" << particle.E << ", " << particle.px
         << ", " << particle.py << ", " << particle.pz << ")" << std::endl;
         aasdf.push_back(particle);
         }
         std::cout << "inv masses: " << aasdf[0].invariantMass(aasdf[1]) << " "
         << aasdf[0].invariantMass(aasdf[2]) << " "
         << aasdf[1].invariantMass(aasdf[2]) << std::endl;

         intensity(max_point);
         std::cout << tree_.at(0)->print() << std::endl;*/

      }
      if (i % 1000 == 0) {
        std::cout << i << ": " << intens << std::endl;
        //std::cout<< tree_.at(0)->print() <<std::endl;
      }
    }

    maxAmplitudeValue_ = maxVal;
    wasMaxAmplitudeValueCalculated_ = true;

    BOOST_LOG_TRIVIAL(info)<<"CoherentAmplitude::calcMaxVal() calculated maximum of amplitude: "
    <<maxAmplitudeValue_;
  }

  return maxAmplitudeValue_;
}

bool CoherentAmplitude::hasTree() {
  return true;
}

std::shared_ptr<FunctionTree> CoherentAmplitude::getAmpTree(allMasses&,
    allMasses&, std::string label) {
  if (label.compare("data") == 0) {
    return tree_.at(1);
  }
  else
    return tree_.at(0);
}

const ParameterList& CoherentAmplitude::intensity(std::vector<double> point,
    ParameterList& par) {
  setParameterList(par);
  dataPoint dataP(point);
  return intensity(dataP);
}

const ParameterList& CoherentAmplitude::intensity(const dataPoint& point,
    ParameterList& par) {
  setParameterList(par);
  return intensity(point);
}

const ParameterList& CoherentAmplitude::intensityNoEff(const dataPoint& point) {
//  std::complex<double> coherent_amplitude = 0;
  double intensity(0.0);

  if (Kinematics::instance()->isWithinPhsp(point)) {
    /* for (unsigned int i = 0; i < topology_amplitudes_.size(); ++i) {
     // at first get the appropriate data evaluation index lists for this
     // topology which we pass to the topology amplitudes later
     const std::vector<IndexList> &topology_data_index_lists =
     data_point_index_lists_[i];
     for (unsigned int final_state_combination_index = 0;
     final_state_combination_index < topology_data_index_lists.size();
     ++final_state_combination_index) {
     coherent_amplitude += topology_amplitudes_[i].evaluate(point,
     topology_data_index_lists[final_state_combination_index]);
     }
     }

     // ok now we build the actual intensity

     for(auto const& sum_pairing : coherent_sum_pairings_) {
     std::complex<double> temp(0.0, 0.0);
     for(unsigned int conj_index : sum_pairing.second) {
     temp += std::conj(sequential_decay_amplitude_values_[conj_index]);
     }
     intensity += temp*sequential_decay_amplitude_values_[sum_pairing.first];
     }

     //intensity = std::pow(std::abs(coherent_amplitude), 2.0);
     }*/

    static bool first_time(true);
    if (first_time) {
      DataPointStorage::Instance().layoutDataStorageStructure(0, 1, point);
      first_time = false;
    }
    DataPointStorage::Instance().clearStorage();
    DataPointStorage::Instance().addDataPoint(0, point);

    tree_.at(0)->forceRecalculate();
    intensity = std::abs(
        ((MultiDouble*) tree_.at(0)->head()->getValue().get())->GetValue(0));

    if (intensity != intensity) {
      BOOST_LOG_TRIVIAL(error)<<"Intensity is not a number!!";
      intensity = 0;
    }
    result_value_->SetValue(intensity);
    return result;
  }
}

const ParameterList& CoherentAmplitude::intensity(const dataPoint& point) {
  intensityNoEff(point);
  if (0 != efficiency_.get()) {
    double eff = efficiency_->evaluate(point);
    result_value_->SetValue(result_value_->GetValue() * eff);
  }
  return result;
}

const bool CoherentAmplitude::fillStartParVec(ParameterList& outPar) {
  outPar = ParameterList(parameters_);
  return true;
}

void CoherentAmplitude::setParameterList(const ParameterList& par) {
  //parameters varied by Minimization algorithm
  if (par.GetNDouble() != parameters_.GetNDouble())
    throw std::runtime_error(
        "setParameterList(): size of parameter lists don't match");
  for (unsigned int i = 0; i < parameters_.GetNDouble(); i++) {
    std::shared_ptr<DoubleParameter> p = parameters_.GetDoubleParameter(i);
    if (p != par.GetDoubleParameter(i)) {
      //p->UpdateParameter(par.GetDoubleParameter(i));
      if (!p->IsFixed()) {
        p->SetValue(par.GetDoubleParameter(i)->GetValue());
        if (p->HasError())
          p->SetError(par.GetDoubleParameter(i)->GetError());
      }
    }
  }
  return;
}

bool CoherentAmplitude::copyParameterList(ParameterList& par) {
  par = parameters_;
  return true;
}

double CoherentAmplitude::getIntValue(std::string var1, double min1,
    double max1, std::string var2, double min2, double max2) {
  return 0.0;
}

void CoherentAmplitude::printAmps() {

}
void CoherentAmplitude::printFractions() {

}

Amplitude * CoherentAmplitude::Clone() {
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
