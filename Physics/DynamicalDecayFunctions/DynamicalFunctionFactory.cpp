//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#include "Physics/DynamicalDecayFunctions/DynamicalFunctionFactory.hpp"
#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/RelativisticBreitWigner.hpp"
#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/TopNodeConstantValue.hpp"

namespace DynamicalFunctions {

DynamicalFunctionFactory::DynamicalFunctionFactory() {
  // TODO Auto-generated constructor stub

}

DynamicalFunctionFactory::~DynamicalFunctionFactory() {
  // TODO Auto-generated destructor stub
}

std::shared_ptr<AbstractDynamicalFunction> DynamicalFunctionFactory::generateRelativisiticBreitWigner(
    const HelicityFormalism::TwoBodyDecayInformation& state_info,
    const ParameterList& external_parameters) {

  std::shared_ptr<RelativisticBreitWigner> rel_bw(
      new RelativisticBreitWigner(state_info.spin_info_.initial_state_.J_));
  rel_bw->initialiseParameters(state_info.dynamical_info_.initial_state_,
      external_parameters);

  dynamical_function_list_[state_info] = rel_bw;

  return rel_bw;
}

std::shared_ptr<AbstractDynamicalFunction> DynamicalFunctionFactory::generateDynamicalFunction(
    const HelicityFormalism::TwoBodyDecayInformation& state_info,
    const ParameterList& external_parameters) {

  // first check if we already have this dynamical function
  auto find_result = dynamical_function_list_.find(state_info);
  if (find_result != dynamical_function_list_.end()) {
    return find_result->second;
  }
  else {
    if (state_info.dynamical_info_.initial_state_.get<std::string>("type")
        == "relBW") {
      return generateRelativisiticBreitWigner(state_info, external_parameters);
    }
    else if (state_info.dynamical_info_.initial_state_.get<std::string>("type")
        == "topNode") {
      return std::shared_ptr<TopNodeConstantValue>(new TopNodeConstantValue());
    }
  }
}

} /* namespace DynamicalFunctions */