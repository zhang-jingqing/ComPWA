//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#include "Core/DataPointStorage.hpp"
#include "Core/Kinematics.hpp"
#include "Physics/DynamicalDecayFunctions/TwoBodyDecay/TopNodeConstantValue.hpp"
#include "Physics/DynamicalDecayFunctions/Kinematics.hpp"

namespace ComPWA {
namespace Physics {
namespace DynamicalFunctions {

TopNodeConstantValue::TopNodeConstantValue(const ParticleStateInfo& psi,
    const ExternalParameters& external_parameters) {
  // TODO Auto-generated constructor stub
  J_ = psi.spin_information_;
  resonance_mass_.reset(
      new DoubleParameter(psi.pid_information_.name_ + "_mass"));
  meson_radius_.reset(
      new DoubleParameter(psi.pid_information_.name_ + "_meson_radius"));

  index_cms_mass_squared_ = ComPWA::Kinematics::instance()->getVariableIndex(
      "cms_mass_squared");

  initialiseParameters(psi.dynamical_information_, external_parameters);
}

TopNodeConstantValue::~TopNodeConstantValue() {
  // TODO Auto-generated destructor stub
}

void TopNodeConstantValue::initialiseParameters(
    const boost::property_tree::ptree& parameter_info,
    const ExternalParameters& external_parameters) {
  resonance_mass_->SetValue(
      parameter_info.get_child("mass").get<double>("value"));
  if (parameter_info.get_child("mass").get<double>("fix"))
    resonance_mass_->SetParameterFixed();
  else
    resonance_mass_->SetParameterFree();

  resonance_mass_->SetMinValue(
      parameter_info.get_child("mass").get<double>("min"));
  resonance_mass_->SetMaxValue(
      parameter_info.get_child("mass").get<double>("max"));
  //resonance_mass_->SetUseBounds(false);

  auto meson_radius_pt = parameter_info.get_child("mesonRadius");
  if (meson_radius_pt.get_optional<double>("value")) {
    meson_radius_->SetValue(meson_radius_pt.get<double>("value"));
    if (meson_radius_pt.get<double>("fix"))
      meson_radius_->SetParameterFixed();
    else
      meson_radius_->SetParameterFree();
    meson_radius_->SetMinValue(meson_radius_pt.get<double>("min"));
    meson_radius_->SetMaxValue(meson_radius_pt.get<double>("max"));
  }
  else {
    meson_radius_->SetValue(parameter_info.get<double>("mesonRadius"));
    meson_radius_->SetParameterFixed();
  }

  // try to extract daughter masses from external parameters
  daughter1_mass_ = external_parameters.parameters_.GetDoubleParameter(
      external_parameters.parameter_name_mapping_.at("daughter1_mass"));
  daughter2_mass_ = external_parameters.parameters_.GetDoubleParameter(
      external_parameters.parameter_name_mapping_.at("daughter2_mass"));
}

double TopNodeConstantValue::phaseSpaceFactor(double two_body_mass) const {
  double ma = daughter1_mass_->GetValue();
  double mb = daughter2_mass_->GetValue();
  return std::sqrt(
      (1.0 - std::pow((ma + mb) / two_body_mass, 2))
          * (1.0 - std::pow((ma - mb) / two_body_mass, 2)));
}

double TopNodeConstantValue::qFactor(double two_body_mass) const {
  return 0.5 * phaseSpaceFactor(two_body_mass) * two_body_mass;
}

double TopNodeConstantValue::barrierTerm(double q, double q0) const {
  unsigned int J(J_.J_numerator_ / J_.J_denominator_);
  double z = std::pow(q * meson_radius_->GetValue(), 2);
  double z0 = std::pow(q0 * meson_radius_->GetValue(), 2);
  if (J == 0) {
    return 1;
  }
  else if (J == 1) {
    return std::sqrt((1.0 + z0) / (1.0 + z));
  }
  else if (J == 2) {
    return std::sqrt(
        (std::pow(z0 - 3.0, 2) + 9.0 * z0) / (std::pow(z - 3.0, 2) + 9.0 * z));
  }
  else {
    std::runtime_error(
        "RelativisticBreitWigner::barrierTerm: barrier factor only implemented for J <= 2");
  }
}

std::complex<double> TopNodeConstantValue::evaluate(const dataPoint& point,
    unsigned int evaluation_index) const {
  //return std::complex<double>(1.0, 0.0);
  double two_body_mass = point.getVal(evaluation_index + index_cms_mass_squared_);
  double q(qFactor (two_body_mass));
  double q0(qFactor(resonance_mass_->GetValue()));

  return barrierTerm(q, q0);
}

std::complex<double> TopNodeConstantValue::evaluate(unsigned int storage_index,
    unsigned int data_index, unsigned int evaluation_index) const {
  auto const& data_vec = DataPointStorage::Instance().getDataList(storage_index,
      evaluation_index + index_cms_mass_squared_);

  double two_body_mass = data_vec[data_index];

  double q(qFactor (two_body_mass));
  double q0(qFactor(resonance_mass_->GetValue()));

  return barrierTerm(q, q0);
  //return std::complex<double>(1.0, 0.0);
}

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */
