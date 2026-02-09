/*
 * Copyright (c) 2023 - 2026 Chair for Design Automation, TUM
 * Copyright (c) 2025 - 2026 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

//
// Created by cpsch on 11.12.2025.
//
#include "na/zoned/decomposer/NativeGateDecomposer.hpp"

#include "ir/operations/CompoundOperation.hpp"

#include <complex>
#include <vector>

namespace na::zoned {

auto NativeGateDecomposer::convertGateToQuaternion(
    std::reference_wrapper<const qc::Operation> op) -> Quaternion {
  assert(op.get().getNqubits() == 1);
  Quaternion quat{};
  if (op.get().getType() == qc::RZ || op.get().getType() == qc::P) {
    quat = {cos(op.get().getParameter().front() / 2), 0, 0,
            sin(op.get().getParameter().front() / 2)};
  } else if (op.get().getType() == qc::Z) {
    quat = {0, 0, 0, 1};
  } else if (op.get().getType() == qc::S) {
    quat = {cos(qc::PI_4), 0, 0, sin(qc::PI_4)};
  } else if (op.get().getType() == qc::Sdg) {
    quat = {cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)};
  } else if (op.get().getType() == qc::T) {
    quat = {cos(qc::PI_4 / 2), 0, 0, sin(qc::PI_4 / 2)};
  } else if (op.get().getType() == qc::Tdg) {
    quat = {cos(-qc::PI_4 / 2), 0, 0, sin(-qc::PI_4 / 2)};
  } else if (op.get().getType() == qc::U) {
    quat = combineQuaternions(
        combineQuaternions({cos(op.get().getParameter().at(1) / 2), 0, 0,
                            sin(op.get().getParameter().at(1) / 2)},
                           {cos(op.get().getParameter().front() / 2), 0,
                            sin(op.get().getParameter().front() / 2), 0}),
        {cos(op.get().getParameter().at(2) / 2), 0, 0,
         sin(op.get().getParameter().at(2) / 2)});
  } else if (op.get().getType() == qc::U2) {
    quat = combineQuaternions(
        combineQuaternions({cos(op.get().getParameter().front() / 2), 0, 0,
                            sin(op.get().getParameter().front() / 2)},
                           {cos(qc::PI_4), 0, sin(qc::PI_4), 0}),
        {cos(op.get().getParameter().at(1) / 2), 0, 0,
         sin(op.get().getParameter().at(1) / 2)});
  } else if (op.get().getType() == qc::RX) {
    quat = {cos(op.get().getParameter().front() / 2),
            sin(op.get().getParameter().front() / 2), 0, 0};
  } else if (op.get().getType() == qc::RY) {
    quat = {cos(op.get().getParameter().front() / 2), 0,
            sin(op.get().getParameter().front() / 2), 0};
  } else if (op.get().getType() == qc::H) {
    quat = combineQuaternions(
        combineQuaternions({1, 0, 0, 0}, {cos(qc::PI_4), 0, sin(qc::PI_4), 0}),
        {cos(qc::PI_2), 0, 0, sin(qc::PI_2)});
  } else if (op.get().getType() == qc::X) {
    quat = {0, 1, 0, 0};
  } else if (op.get().getType() == qc::Y) {
    quat = {0, 0, 1, 0};
  } else if (op.get().getType() == qc::Vdg) {
    quat = combineQuaternions(
        combineQuaternions({cos(qc::PI_4), 0, 0, sin(qc::PI_4)},
                           {cos(-qc::PI_4), 0, sin(-qc::PI_4), 0}),
        {cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)});
  } else if (op.get().getType() == qc::SX) {
    quat = combineQuaternions(
        combineQuaternions({cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)},
                           {cos(qc::PI_4), 0, sin(qc::PI_4), 0}),
        {cos(qc::PI_4), 0, 0, sin(qc::PI_4)});
  } else if (op.get().getType() == qc::SXdg || op.get().getType() == qc::V) {
    quat = combineQuaternions(
        combineQuaternions({cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)},
                           {cos(-qc::PI_4), 0, sin(-qc::PI_4), 0}),
        {cos(qc::PI_4), 0, 0, sin(qc::PI_4)});
  } else {
    // if the gate type is not recognized, an error is printed and the
    // gate is not included in the output.
    std::ostringstream oss;
    oss << "ERROR: Unsupported single-qubit gate: " << op.get().getType()
        << "\n";
    throw std::invalid_argument(oss.str());
  }
  return quat;
}

auto NativeGateDecomposer::combineQuaternions(const Quaternion& q1,
                                              const Quaternion& q2)
    -> Quaternion {
  Quaternion new_quat{};
  new_quat[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  new_quat[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  new_quat[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  new_quat[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
  return new_quat;
}

auto NativeGateDecomposer::getU3AnglesFromQuaternion(const Quaternion& quat)
    -> std::array<qc::fp, 3> {
  qc::fp theta;
  qc::fp phi;
  qc::fp lambda;
  if (std::fabs(quat[0]) > epsilon || std::fabs(quat[3]) > epsilon) {
    theta = 2. * std::atan2(std::sqrt(quat[2] * quat[2] + quat[1] * quat[1]),
                            std::sqrt(quat[0] * quat[0] + quat[3] * quat[3]));
    qc::fp alpha_1 = std::atan2(quat[3], quat[0]); // phi+ lambda
    if (std::fabs(quat[1]) > epsilon || std::fabs(quat[2]) > epsilon) {
      qc::fp alpha_2 = -1 * std::atan2(quat[1], quat[2]);
      phi = alpha_1 + alpha_2; // phi
      lambda = alpha_1 - alpha_2;
    } else {
      phi = 0;
      lambda = 2 * alpha_1;
    }
  } else {
    theta = qc::PI; // or § PI if sin(theta/2)=-1... Relevant? Or is the -1=
                    // exp(i*pi) just global phase
    if (std::fabs(quat[1]) > epsilon || std::fabs(quat[2]) > epsilon) {
      phi = 0;
      lambda = 2 * std::atan2(quat[1], quat[2]);
      // atan can give PI instead of 0 Problem?
    } else {
      // This should never happen! Exception??
      phi = 0.;
      lambda = 0.;
    }
  }
  return {theta, phi, lambda};
}

auto NativeGateDecomposer::calcThetaMax(const std::vector<StructU3>& layer)
    -> qc::fp {
  qc::fp theta_max = 0;
  for (auto gate : layer) {
    if (std::fabs(gate.angles[0]) > theta_max) {
      theta_max = std::fabs(gate.angles[0]);
    }
  }
  return theta_max;
}
auto NativeGateDecomposer::transformToU3(
    const std::vector<SingleQubitGateRefLayer>& layers) const
    -> std::vector<std::vector<StructU3>> {
  std::vector<std::vector<StructU3>> new_layers;
  for (const auto& layer : layers) {
    std::vector<std::vector<std::reference_wrapper<const qc::Operation>>> gates(
        this->nQubits_);
    std::vector<StructU3> new_layer;
    for (auto gate : layer) {
      gates[gate.get().getTargets().front()].push_back(gate);
    }

    for (auto qubit_gates : gates) {
      if (!qubit_gates.empty()) {
        std::array<qc::fp, 4> quat = convertGateToQuaternion(qubit_gates[0]);
        for (size_t i = 1; i < qubit_gates.size(); i++) {
          quat =
              combineQuaternions(quat, convertGateToQuaternion(qubit_gates[i]));
        }
        std::array<qc::fp, 3> angles = getU3AnglesFromQuaternion(quat);
        new_layer.emplace_back(
            StructU3{angles, qubit_gates[0].get().getTargets().front()});
      }
    }
    new_layers.push_back(new_layer);
  }
  return new_layers;
}
auto NativeGateDecomposer::getDecompositionAngles(
    const std::array<qc::fp, 3>& angles, qc::fp theta_max)
    -> std::array<qc::fp, 3> {
  qc::fp alpha, chi;
  // U3(theta,phi_min(phi),phi_plus(lambda)->Rz(gamma_minus)GR(theta_max/2,
  // PI_2)Rz(chi)GR(-theta_max/2,PI_2)RZ(gamma_plus)
  qc::fp sin_sq_diff = (sin(theta_max / 2) * sin(theta_max / 2) -
                        (sin(angles[0] / 2) * sin(angles[0] / 2)));
  if (std::fabs(sin_sq_diff) < epsilon) {
    chi = qc::PI;
    if (std::fabs(cos(theta_max / 2)) < epsilon) { // Periodicity covered?
      alpha = 0;
    } else {
      alpha = qc::PI_2;
    }
  } else {
    qc::fp kappa =
        std::sqrt((sin(angles[0] / 2) * sin(angles[0] / 2)) / sin_sq_diff);
    alpha = atan(cos(theta_max / 2) * kappa);
    chi = fmod(2 * atan(kappa), qc::TAU);
  }
  qc::fp beta = angles[0] < 0 ? -1 * qc::PI_2 : qc::PI_2;
  qc::fp gamma_plus = fmod(angles[2] - (alpha + beta), qc::TAU);
  qc::fp gamma_minus = fmod(angles[1] - (alpha - beta), qc::TAU);

  return {chi, gamma_minus, gamma_plus};
}

auto NativeGateDecomposer::decompose(
    const size_t nQubits,
    const std::vector<SingleQubitGateRefLayer>& singleQubitGateLayers)
    -> std::vector<SingleQubitGateLayer> {
  nQubits_ = nQubits;

  std::vector<std::vector<StructU3>> U3Layers =
      transformToU3(singleQubitGateLayers);
  std::vector<SingleQubitGateLayer> NewSingleQubitLayers =
      std::vector<SingleQubitGateLayer>{};

  for (const auto& layer : U3Layers) {
    qc::fp theta_max = calcThetaMax(layer);
    SingleQubitGateLayer FrontLayer;
    SingleQubitGateLayer MidLayer;
    SingleQubitGateLayer BackLayer;
    SingleQubitGateLayer NewLayer;

    for (auto gate : layer) {
      std::array<qc::fp, 3> decomp_angles =
          getDecompositionAngles(gate.angles, theta_max);

      // GR(theta_max/2, PI_2)==Global Y due to PI_2
      FrontLayer.emplace_back(std::make_unique<const qc::StandardOperation>(
          qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[1]})));

      MidLayer.emplace_back(std::make_unique<const qc::StandardOperation>(
          qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[0]})));

      BackLayer.emplace_back(std::make_unique<const qc::StandardOperation>(
          qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[2]})));
    } // gate::layer

    std::vector<std::unique_ptr<qc::Operation>> GR_plus;
    std::vector<std::unique_ptr<qc::Operation>> GR_minus;

    for (size_t i = 0; i < this->nQubits_; ++i) {
      GR_plus.emplace_back(std::make_unique<qc::StandardOperation>(
          i, qc::RY, std::initializer_list<qc::fp>{theta_max / 2}));
      GR_minus.emplace_back(std::make_unique<qc::StandardOperation>(
          i, qc::RY, std::initializer_list<qc::fp>{-1 * theta_max / 2}));
    }

    for (auto&& gate : FrontLayer) {
      NewLayer.push_back(std::move(gate));
    }

    NewLayer.emplace_back(
        std::move(std::make_unique<const qc::CompoundOperation>(
            qc::CompoundOperation(std::move(GR_plus), true))));

    for (auto&& gate : MidLayer) {
      NewLayer.push_back(std::move(gate));
    }
    NewLayer.emplace_back(
        std::move(std::make_unique<const qc::CompoundOperation>(
            qc::CompoundOperation(std::move(GR_minus), true))));

    for (auto&& gate : BackLayer) {
      NewLayer.push_back(std::move(gate));
    }
    NewSingleQubitLayers.push_back(std::move(NewLayer));
  } // layer::SingleQubitLayers
  return NewSingleQubitLayers;
}
} // namespace na::zoned
