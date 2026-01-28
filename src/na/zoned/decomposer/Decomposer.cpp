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
#include "na/zoned/decomposer/decomposer.hpp"

// #include "ir/operations/StandardOperation.hpp"
#include "ir/operations/CompoundOperation.hpp"

#include <complex>
#include <vector>

namespace na::zoned {
Decomposer::Decomposer(int n_qubits) { N_qubits = n_qubits; };

auto Decomposer::convert_gate_to_quaternion(
    std::reference_wrapper<const qc::Operation> op) -> std::array<qc::fp, 4> {
  assert(op.get().getNqubits() == 1);
  std::array<qc::fp, 4> quat{};
  // TODO: Phase?
  if (op.get().getType() == qc::RZ || op.get().getType() == qc::P) {
    quat = {cos(op.get().getParameter().front()), 0, 0,
            sin(op.get().getParameter().front())};
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
    quat = combine_quaternions(
        combine_quaternions({cos(op.get().getParameter().at(1) / 2), 0, 0,
                             sin(op.get().getParameter().at(1) / 2)},
                            {cos(op.get().getParameter().front() / 2), 0,
                             sin(op.get().getParameter().front() / 2), 0}),
        {cos(op.get().getParameter().at(2) / 2), 0, 0,
         sin(op.get().getParameter().at(2) / 2)});
  } else if (op.get().getType() == qc::U2) {
    quat = combine_quaternions(
        combine_quaternions({cos(op.get().getParameter().front() / 2), 0, 0,
                             sin(op.get().getParameter().front() / 2)},
                            {cos(qc::PI_2), 0, sin(qc::PI_2), 0}),
        {cos(op.get().getParameter().at(1) / 2), 0, 0,
         sin(op.get().getParameter().at(1) / 2)});
  } else if (op.get().getType() == qc::RX) {
    quat = {cos(op.get().getParameter().front() / 2),
            sin(op.get().getParameter().front() / 2), 0, 0};
  } else if (op.get().getType() == qc::RY) {
    quat = {cos(op.get().getParameter().front() / 2), 0,
            sin(op.get().getParameter().front() / 2), 0};
  } else if (op.get().getType() == qc::H) {
    quat = combine_quaternions(
        combine_quaternions({1, 0, 0, 0}, {cos(qc::PI_4), 0, sin(qc::PI_4), 0}),
        {cos(qc::PI_2), 0, 0, sin(qc::PI_2)});
  } else if (op.get().getType() == qc::X) {
    quat = {0, 1, 0, 0};
  } else if (op.get().getType() == qc::Y) {
    quat = {0, 1, 0, 0};
  } else if (op.get().getType() == qc::Vdg) {
    quat = combine_quaternions(
        combine_quaternions({cos(qc::PI_4), 0, 0, sin(qc::PI_4)},
                            {cos(-qc::PI_4), 0, sin(-qc::PI_4), 0}),
        {cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)});
  } else if (op.get().getType() == qc::SX) {
    quat = combine_quaternions(
        combine_quaternions({cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)},
                            {cos(qc::PI_4), 0, sin(qc::PI_4), 0}),
        {cos(qc::PI_4), 0, 0, sin(qc::PI_4)});
  } else if (op.get().getType() == qc::SXdg || op.get().getType() == qc::V) {
    quat = combine_quaternions(
        combine_quaternions({cos(-qc::PI_4), 0, 0, sin(-qc::PI_4)},
                            {cos(-qc::PI_4), 0, sin(-qc::PI_4), 0}),
        {cos(qc::PI_4), 0, 0, sin(qc::PI_4)});
  } else {
    // if the gate type is not recognized, an error is printed and the
    // gate is not included in the output.
    std::ostringstream oss;
    oss << "\033[1;31m[ERROR]\033[0m Unsupported single-qubit gate: "
        << op.get().getType() << "\n";
    throw std::invalid_argument(oss.str());
  }
  return quat;
}

auto Decomposer::combine_quaternions(const std::array<qc::fp, 4>& q1,
                                     const std::array<qc::fp, 4>& q2)
    -> std::array<qc::fp, 4> {
  std::array<qc::fp, 4> new_quat{};
  new_quat[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  new_quat[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  new_quat[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  new_quat[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
  return new_quat;
}

auto Decomposer::get_U3_angles_from_quaternion(
    const std::array<qc::fp, 4>& quat) -> std::array<qc::fp, 3> {
  // TODO: Is there a prescribed eps somewhere? Else where should THIS eps be
  // defined
  qc::fp theta;
  qc::fp phi;
  qc::fp lambda;
  if (abs(quat[0]) > epsilon || abs(quat[3]) > epsilon) {
    theta = 2. * std::atan2(std::sqrt(quat[2] * quat[2] + quat[1] * quat[1]),
                            std::sqrt(quat[0] * quat[0] + quat[3] * quat[3]));
    qc::fp alpha_1 = std::atan2(quat[3], quat[0]); // phi+ lambda
    if (abs(quat[1]) > epsilon || abs(quat[2]) > epsilon) {
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
    if (abs(quat[1]) > epsilon || abs(quat[2]) > epsilon) {
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

auto Decomposer::calc_theta_max(const std::vector<struct_U3>& layer) -> qc::fp {
  qc::fp theta_max = 0;
  for (auto gate : layer) {
    if (abs(gate.angles[0]) > theta_max) {
      theta_max = abs(gate.angles[0]);
    }
  }
  return theta_max;
}
auto Decomposer::transform_to_U3(
    const std::vector<SingleQubitGateRefLayer>& layers) const
    -> std::vector<std::vector<struct_U3>> {
  // auto u=struct_U3({0,0,0},0);
  std::vector<std::vector<struct_U3>> new_layers;
  for (const auto& layer : layers) {
    std::vector<std::vector<std::reference_wrapper<const qc::Operation>>> gates(
        this->N_qubits);
    std::vector<struct_U3> new_layer;
    for (auto gate : layer) {
      gates[gate.get().getTargets().front()].push_back(gate);
    }

    for (auto qubit_gates : gates) {
      if (!qubit_gates.empty()) {
        std::array<qc::fp, 4> quat = convert_gate_to_quaternion(qubit_gates[0]);
        for (auto i = 1; i < qubit_gates.size(); i++) {
          quat = combine_quaternions(
              quat, convert_gate_to_quaternion(qubit_gates[i]));
        }
        std::array<qc::fp, 3> angles = get_U3_angles_from_quaternion(quat);
        new_layer.emplace_back(
            struct_U3(angles, qubit_gates[0].get().getTargets().front()));
      }
    }
    new_layers.push_back(new_layer);
  }
  return new_layers;
}
auto Decomposer::get_decomposition_angles(const std::array<qc::fp, 3>& angles,
                                          qc::fp theta_max)
    -> std::array<qc::fp, 3> {
  qc::fp alpha, chi, beta;
  // U3(theta,phi_min(phi),phi_plus(lambda)->Rz(gamma_minus)GR(theta_max/2,
  // PI_2)Rz(chi)GR(-theta_max/2,PI_2)RZ(gamma_plus)
  if (abs(angles[0] - theta_max) < epsilon) {
    chi = qc::PI;
    if (abs(cos(theta_max / 2)) < epsilon) { // Periodicity covered?
      alpha = 0;
    } else {
      alpha = qc::PI_2;
    }
  } else {
    qc::fp kappa = sqrt((sin(angles[0] / 2) * sin(angles[0] / 2)) /
                        (sin(theta_max) * sin(theta_max) -
                         (sin(angles[0] / 2) * sin(angles[0] / 2))));
    alpha = atan(cos(theta_max / 2) * kappa);
    chi = fmod(2 * atan(kappa), qc::TAU);
  }
  beta = angles[0] < 0 ? -1 * qc::PI_2 : qc::PI_2;
  qc::fp gamma_plus = fmod(angles[2] - (alpha + beta), qc::TAU);
  qc::fp gamma_minus = fmod(angles[1] - (alpha - beta), qc::TAU);

  return {chi, gamma_minus, gamma_plus};
}

auto Decomposer::decompose(
    const std::vector<SingleQubitGateRefLayer>& singleQubitGateLayers) const
    -> std::vector<SingleQubitGateLayer> {

  std::vector<std::vector<struct_U3>> U3Layers =
      transform_to_U3(singleQubitGateLayers);
  std::vector<SingleQubitGateLayer> NewSingleQubitLayers =
      std::vector<SingleQubitGateLayer>{};

  for (const auto& layer : U3Layers) {
    qc::fp theta_max = calc_theta_max(layer);
    SingleQubitGateLayer FrontLayer;
    SingleQubitGateLayer MidLayer;
    SingleQubitGateLayer BackLayer;
    SingleQubitGateLayer NewLayer;

    for (auto gate : layer) {
      std::array<qc::fp, 3> decomp_angles =
          get_decomposition_angles(gate.angles, theta_max);

      // GR(theta_max/2, PI_2)==Global Y due to PI_2
      auto sop = qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[1]});
      std::unique_ptr<const qc::Operation> op =
          std::make_unique<const qc::StandardOperation>(sop);
      FrontLayer.emplace_back(std::move(op));

      sop = qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[0]});
      op = std::make_unique<const qc::StandardOperation>(sop);
      MidLayer.emplace_back(std::move(op));

      sop = qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[2]});
      op = std::make_unique<const qc::StandardOperation>(sop);
      BackLayer.emplace_back(std::move(op));
    } // gate::layer

    std::vector<std::unique_ptr<qc::Operation>> GR_plus;
    std::vector<std::unique_ptr<qc::Operation>> GR_minus;

    for (auto i = 0; i < this->N_qubits; ++i) {
      GR_plus.emplace_back(
          new qc::StandardOperation(i, qc::RY, {theta_max / 2}));
      GR_minus.emplace_back(
          new qc::StandardOperation(i, qc::RY, {-1 * theta_max / 2}));
    }

    for (auto&& gate : FrontLayer) {
      NewLayer.push_back(std::move(gate));
    }

    auto cop = qc::CompoundOperation(std::move(GR_plus), true);
    std::unique_ptr<const qc::Operation> ryp =
        std::make_unique<const qc::CompoundOperation>(cop);
    NewLayer.emplace_back(std::move(ryp));

    for (auto&& gate : MidLayer) {
      NewLayer.push_back(std::move(gate));
    }
    cop = qc::CompoundOperation(std::move(GR_minus), true);
    std::unique_ptr<const qc::Operation> rym =
        std::make_unique<const qc::CompoundOperation>(cop);
    NewLayer.emplace_back(std::move(rym));

    for (auto&& gate : BackLayer) {
      NewLayer.push_back(std::move(gate));
    }
    NewSingleQubitLayers.push_back(std::move(NewLayer));
  } // layer::SingleQubitLayers
  return NewSingleQubitLayers;
}
} // namespace na::zoned
