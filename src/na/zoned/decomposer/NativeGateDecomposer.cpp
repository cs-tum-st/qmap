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
      // WHat are operations with empty targets doing??
      if (!gate.get().getTargets().empty()) {
        gates[gate.get().getTargets().front()].push_back(gate);
      }
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
    qc::fp kappa = std::sqrt((sin(angles[0] / 2) * sin(angles[0] / 2)) /sin_sq_diff);
    alpha = atan(cos(theta_max / 2) * kappa);
    chi = fmod(2 * atan(kappa), qc::TAU);
  }
  qc::fp beta = angles[0] < 0 ? -1 * qc::PI_2 : qc::PI_2;
  qc::fp gamma_plus = fmod(angles[2] - (alpha + beta), qc::TAU);
  qc::fp gamma_minus = fmod(angles[1] - (alpha - beta), qc::TAU);

  return {chi, gamma_minus, gamma_plus};
}

auto NativeGateDecomposer::decompose(
    size_t nQubits,
    const std::pair<std::vector<SingleQubitGateRefLayer>,
                    std::vector<TwoQubitGateLayer>>& asap_schedule)
    -> std::pair<std::vector<SingleQubitGateLayer>,
                 std::vector<TwoQubitGateLayer>> {
  nQubits_ = nQubits;

  std::vector<std::vector<StructU3>> U3Layers =
      transformToU3(asap_schedule.first);
  std::vector<TwoQubitGateLayer> NewTwoQubitGateLayers = asap_schedule.second;
  if (theta_opt_schedule) {
    // TODO: Put in scheduling function here
    auto thetaOptSchedule =
        schedule(std::pair(U3Layers, NewTwoQubitGateLayers));
    U3Layers = thetaOptSchedule.first;
    NewTwoQubitGateLayers = thetaOptSchedule.second;
  }
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
      FrontLayer.emplace_back(std::make_unique<const qc::StandardOperation>(qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[1]})));

      MidLayer.emplace_back(std::make_unique<const qc::StandardOperation>(qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[0]})));

      BackLayer.emplace_back(std::make_unique<const qc::StandardOperation>(qc::StandardOperation(gate.qubit, qc::RZ, {decomp_angles[2]})));
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

    NewLayer.emplace_back(std::make_unique<const qc::CompoundOperation>(qc::CompoundOperation(std::move(GR_plus), true)));

    for (auto&& gate : MidLayer) {
      NewLayer.push_back(std::move(gate));
    }
    NewLayer.emplace_back(std::make_unique<const qc::CompoundOperation>(qc::CompoundOperation(std::move(GR_minus), true)));

    for (auto&& gate : BackLayer) {
      NewLayer.push_back(std::move(gate));
    }
    NewSingleQubitLayers.push_back(std::move(NewLayer));
  } // layer::SingleQubitLayers
  return {std::move(NewSingleQubitLayers), NewTwoQubitGateLayers};
}

auto NativeGateDecomposer::preprocess(
    const std::pair<std::vector<std::vector<StructU3>>,
                    std::vector<TwoQubitGateLayer>>& schedule)
    -> std::pair<std::vector<SingleQubitGateLayer>,
                 std::vector<SingleQubitGateLayer>> {
  std::vector<SingleQubitGateLayer> NewSingleQubitLayers =
      std::vector<SingleQubitGateLayer>{};
  SingleQubitGateLayer NewLayer;
  for (const auto& singleQubitLayer : schedule.first) {
    NewLayer.clear();
    for (auto gate : singleQubitLayer) {
      NewLayer.emplace_back(
          std::make_unique<const qc::StandardOperation>(qc::StandardOperation(
              (qc::Qubit)gate.qubit, qc::U,
              {gate.angles[0], gate.angles[1], gate.angles[2]})));
    }
    NewSingleQubitLayers.push_back(std::move((NewLayer)));
  }
  std::vector<SingleQubitGateLayer> NewTwoQubitLayers =
      std::vector<SingleQubitGateLayer>{};

  for (const auto& twoQubitLayer : schedule.second) {
    NewLayer.clear();
    for (auto gate : twoQubitLayer) {
      NewLayer.emplace_back(
          std::make_unique<const qc::StandardOperation>(qc::StandardOperation(
              qc::Control((qc::Qubit)gate[0]), (qc::Qubit)gate[1], qc::Z, {})));
    }
    NewTwoQubitLayers.push_back(std::move((NewLayer)));
  }
  return std::pair<std::vector<SingleQubitGateLayer>,
                   std::vector<SingleQubitGateLayer>>(
      std::move(NewSingleQubitLayers), std::move(NewTwoQubitLayers));
}

std::vector<size_t> NativeGateDecomposer::find_shortest_path(
    const DiGraph<std::pair<std::vector<std::size_t>,
                            std::vector<std::size_t>>>& subproblem_graph,
    const std::vector<std::size_t>& leaf_nodes) {
  std::set<std::size_t> leafs =
      std::set<std::size_t>(leaf_nodes.begin(), leaf_nodes.end());
  std::pair<std::vector<std::size_t>, double> leaf_path =
      shortest_path_to_start(subproblem_graph, 0, leafs);
  std::reverse(leaf_path.first.begin(), leaf_path.first.end());
  return leaf_path.first;
}
bool disjunct(const std::set<size_t>& set1, const std::set<size_t>& set2) {
  for (auto elem : set1) {
    if (set2.contains(elem)) {
      return false;
    }
  }
  return true;
}

bool disjunct(const std::set<qc::Qubit>& set1,
              const std::set<qc::Qubit>& set2) {
  for (auto elem : set1) {
    if (set2.contains(elem)) {
      return false;
    }
  }
  return true;
}

std::pair<std::vector<std::size_t>, double>
NativeGateDecomposer::shortest_path_to_start(
    const DiGraph<std::pair<std::vector<std::size_t>,
                            std::vector<std::size_t>>>& subproblem_graph,
    std::size_t curr_node, std::set<std::size_t> leaf_nodes) {
  std::vector<std::pair<std::vector<std::size_t>, double>> possible_paths = {};
  // Check if leaf nodes are reached
  for (auto node : subproblem_graph.get_adjacent(curr_node)) {
    if (leaf_nodes.contains(node.first)) {
      possible_paths.push_back({std::pair<std::vector<std::size_t>, double>(
          {node.first, curr_node}, node.second)});
    }
  }
  // Recursive Case
  if (possible_paths.empty()) {
    for (auto node : subproblem_graph.get_adjacent(curr_node)) {
      auto path =
          shortest_path_to_start(subproblem_graph, node.first, leaf_nodes);
      path.first.push_back(curr_node);
      path.second += node.second;
      possible_paths.push_back(path);
    }
  }
  // Base Case:
  // Choose shortest Paths
  auto min_length = possible_paths.at(0).first.size();
  std::vector<std::pair<std::vector<std::size_t>, double>> shortest_paths = {};
  for (auto path : possible_paths) {
    if (path.first.size() < min_length) {
      min_length = path.first.size();
    }
  }
  for (auto path : possible_paths) {
    if (path.first.size() == min_length) {
      ;
      shortest_paths.push_back(path);
    }
  }
  if (shortest_paths.size() == 1) {
    return shortest_paths.at(0);
  }
  // Find shortest path with minimal cost
  auto min_cost = shortest_paths.at(0).second;
  auto best_path = shortest_paths.at(0);
  for (auto path : shortest_paths) {
    if (path.second > min_cost) {
      min_cost = path.second;
      best_path = path;
    }
  }
  return best_path;
}
double NativeGateDecomposer::calc_cost(
    std::vector<std::vector<std::size_t>>::const_reference path,
    DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>&
        subproblem_graph) {
  double cost = 0;
  for (size_t node = 0; node < path.size() - 1; node++) {
    for (const auto [child, weight] :
         subproblem_graph.get_adjacent(path[node])) {
      if (child == path[node + 1]) {
        cost += weight;
        break;
      }
    }
  }
  return cost;
}
std::vector<std::size_t> NativeGateDecomposer::find_leaf_nodes(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const DiGraph<std::pair<std::vector<std::size_t>,
                            std::vector<std::size_t>>>& subproblem_graph) {
  std::vector<std::size_t> end_nodes = std::vector<std::size_t>{};
  for (size_t i = 0; i < subproblem_graph.get_N_Nodes(); i++) {
    if (subproblem_graph.get_adjacent(i).empty()) {
      end_nodes.push_back(i);
    }
  }
  return end_nodes;
}
std::vector<std::size_t> NativeGateDecomposer::sort_by_theta_dec(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<unsigned long long>& vector) {
  std::vector<std::size_t> sorted = {};
}
std::vector<std::size_t>
NativeGateDecomposer::remove_element(const std::vector<std::size_t>& vector,
                                     std::size_t node) {
  std::vector<std::size_t> new_vector = {};
  for (auto element : vector) {
    if (element != node) {
      new_vector.push_back(element);
    }
  }
  return new_vector;
}

// TODO: ADD last cond Bool??
std::vector<std::pair<std::array<std::vector<std::size_t>, 4>, qc::fp>>
NativeGateDecomposer::get_possible_moments(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<size_t>& v0_c,
    const std::array<std::vector<size_t>, 3>& v_new) {

  std::vector<size_t> v_p1_star = v_new[0];
  std::vector<size_t> v_p1_square = {};
  std::vector<size_t> v_c1_star = v_new[1];
  std::vector<size_t> v_c1_square = {};

  qc::fp v_c0_cost = max_theta(circuit, v0_c);
  qc::fp v_c1_cost = max_theta(circuit, v_new[1]);
  qc::fp orig_comb_cost = v_c0_cost + v_c1_cost;
  qc::fp new_v_C1_cost = std::max(v_c0_cost, v_c1_cost);

  std::array<std::vector<std::size_t>, 4> v_arg = {v0_c, v_new[0], v_new[1],
                                                   v_new[2]};
  std::vector<std::pair<std::array<std::vector<std::size_t>, 4>, qc::fp>> args =
      {std::pair(v_arg, v_c0_cost)};
  // Sort v_0C from highest to lowest theta
  std::vector<std::size_t> v_sort(v0_c);
  auto sort_by_theta = [&circuit](std::size_t a, std::size_t b) {
    return circuit.get_Node_Value(a)->get()->getParameter().front() >
           circuit.get_Node_Value(b)->get()->getParameter().front();
  };
  std::ranges::sort(v_sort, sort_by_theta);
  // TODO: Check Condition 1
  std::vector<std::pair<std::array<std::vector<std::size_t>, 2>,
                        std::pair<std::set<qc::Qubit>, qc::fp>>>
      potential_arg = {};
  auto prev_theta =
      circuit.get_Node_Value(v_sort[0])->get()->getParameter().front();
  auto this_theta = prev_theta;
  // TODO: This should be Qubits not Operations!!! change ARG
  std::set<qc::Qubit> mk_qubits = {
      circuit.get_Node_Value(v_sort[0])->get()->getUsedQubits()};
  for (auto i = 0; i < v_sort.size(); i++) {
    this_theta =
        circuit.get_Node_Value(v_sort[i])->get()->getParameter().front();
    if (this_theta != prev_theta) {
      std::vector<std::size_t> discarded = {v_sort.begin(),
                                            v_sort.begin() + (i - 1)};
      std::vector<std::size_t> kept = {v_sort.begin() + i, v_sort.end()};
      potential_arg.push_back(std::pair<std::array<std::vector<std::size_t>, 2>,
                                        std::pair<std::set<qc::Qubit>, qc::fp>>(
          {kept, discarded},
          std::pair<std::set<qc::Qubit>, qc::fp>(mk_qubits, this_theta)));
      prev_theta = this_theta;
      mk_qubits.clear();
    }
    mk_qubits.merge(circuit.get_Node_Value(v_sort[i])->get()->getUsedQubits());
  }
  std::vector<std::size_t> push_back_nodes = {};
  std::set<qc::Qubit> p_square_qubits = {};

  for (auto pot : potential_arg) {
    // TODO: Check Condition 2
    for (auto node : v_p1_star) {
      std::set<qc::Qubit> qubits =
          circuit.get_Node_Value(node)->get()->getUsedQubits();
      if (!disjunct(qubits, pot.second.first)) {
        push_back_nodes.emplace_back(node);
        p_square_qubits.merge(qubits);
      }
    }
    for (auto node : push_back_nodes) {
      v_p1_star = remove_element(v_p1_star, node);
      v_p1_square.emplace_back(node);
    }

    if (v_p1_star.empty()) {
      break;
    }
    //  TODO: Check Condition 3
    std::set<qc::Qubit> push_qubits = pot.second.first;
    push_qubits.merge(p_square_qubits);
    push_back_nodes.clear();

    for (auto node : v_c1_star) {
      std::set<qc::Qubit> qubits =
          circuit.get_Node_Value(node)->get()->getUsedQubits();
      if (!disjunct(qubits, push_qubits)) {
        push_back_nodes.emplace_back(node);
      }
    }
    for (auto node : push_back_nodes) {
      v_c1_star = remove_element(v_c1_star, node);
      v_c1_square.emplace_back(node);
    }

    if (v_c1_star.empty()) {
      break;
    }
    // TODO Check Condition 4 ?? Find out what it does!!
    bool check_last_cond = false;
    if ((!check_last_cond) ||
        (check_last_cond &&
         pot.second.second + new_v_C1_cost < orig_comb_cost)) {
      v_arg = {pot.first[0], v_p1_star, pot.first[1], v_new[2]};
      for (auto node : v_c1_star) {
        v_arg[2].push_back(node);
      }
      for (auto node : v_c1_square) {
        v_arg[3].push_back(node);
      }
      for (auto node : v_p1_square) {
        v_arg[3].push_back(node);
      }
      args.push_back(std::pair<std::array<std::vector<std::size_t>, 4>, qc::fp>(
          v_arg, pot.second.second));
    }
  }
  return args;
}
std::pair<std::vector<std::vector<NativeGateDecomposer::StructU3>>,
          std::vector<TwoQubitGateLayer>>
NativeGateDecomposer::postprocess(
    std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
        schedule) {
  std::vector<std::vector<NativeGateDecomposer::StructU3>> U3Layers = {};
  for (size_t i = 0; i < schedule.first.size(); i++) {
    U3Layers.emplace_back();
    for (size_t j = 0; j < schedule.first.at(i).size(); j++) {
      U3Layers.back().push_back(
          StructU3({schedule.first.at(i)[j]->getParameter()[0],
                    schedule.first.at(i)[j]->getParameter()[1],
                    schedule.first.at(i)[j]->getParameter()[2]},
                   schedule.first.at(i)[j]->getTargets()[0]));
    }
  }
}

NativeGateDecomposer::DiGraph<std::unique_ptr<const qc::Operation>>
NativeGateDecomposer::convert_circ_to_dag(
    std::pair<std::vector<SingleQubitGateLayer>,
              std::vector<SingleQubitGateLayer>>& qc) {
  NativeGateDecomposer::DiGraph<std::unique_ptr<const qc::Operation>> graph =
      DiGraph<std::unique_ptr<const qc::Operation>>();
  std::vector<std::vector<size_t>> qubit_paths =
      std::vector<std::vector<size_t>>{};
  // TODO:assert that One more sql exists than mql ??
  for (auto i = 0; i < qc.second.size(); ++i) {
    for (const auto& s : qc.first.at(i)) {
      size_t node = graph.add_Node(s->clone());
      for (auto t = 0; t < s->getNtargets(); t++) {
        qubit_paths.at(s->getTargets().at(t)).push_back(node);
      }
    }

    for (const auto& t : qc.second.at(i)) {
      size_t node = graph.add_Node(t->clone());
      qubit_paths.at(t->getControls().begin()->qubit).push_back(node);
      qubit_paths.at(t->getTargets().at(0)).push_back(node);
    }
  }
  for (auto s = 0; s < qc.first.end()->size(); ++s) {
    size_t node = graph.add_Node(qc.first.end()->at(s)->clone());

    for (auto t = 0; t < qc.first.end()->at(s)->getNtargets(); t++) {
      qubit_paths.at(qc.first.end()->at(s)->getTargets().at(t)).push_back(node);
    }
  }
  for (std::size_t i = 0; i < qubit_paths.size(); ++i) {
    for (std::size_t op = 0; op < qubit_paths.at(i).size(); ++op) {
      graph.add_Edge(i, qubit_paths.at(i).at(op));
    }
  }
  return graph;
}

qc::fp NativeGateDecomposer::max_theta(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<unsigned long long>& nodes) {
  qc::fp max_cost = 0;
  for (const auto node : nodes) {
    if ((*circuit.get_Node_Value(node))->getParameter()[0] >= max_cost) {
      max_cost = (*circuit.get_Node_Value(node))->getParameter()[0];
    }
  }
  return max_cost;
}

std::array<std::vector<size_t>, 3> NativeGateDecomposer::sift(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit, size_t n_qubits) {
  std::vector<size_t> v_p = std::vector<size_t>();
  std::vector<size_t> v_c = std::vector<size_t>();
  std::vector<size_t> v_r = std::vector<size_t>();

  std::set<size_t> removed = std::set<size_t>();

  for (auto node = 0; node < circuit.get_N_Nodes(); node++) {
    // TODO: Check for unique pointer issue
    auto op = circuit.get_Node_Value(node);
    std::set<size_t> op_qubits = std::set<size_t>();
    for (auto qubit : op->get()->getUsedQubits()) {
      op_qubits.insert(qubit);
    }
    if (removed.size() < n_qubits && disjunct(removed, op_qubits)) {
      if (op->get()->getNqubits() == 1) {
        v_c.push_back(node);
        removed.insert(op->get()->getTargets().at(0));
      } else {
        v_p.push_back(node);
        // Add something to make it only pick one 2-Qubit gate per Qubit per
        // moment???
      }
    } else {
      v_r.push_back(node);
    }
  }
  return std::array<std::vector<size_t>, 3>({v_p, v_c, v_r});
}

auto NativeGateDecomposer::build_schedule(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>&
        subproblem_graph) -> std::pair<std::vector<SingleQubitGateLayer>,
                                       std::vector<TwoQubitGateLayer>> {

  // TODO: Find Leaf Nodes of di_graph
  std::vector<std::size_t> leaf_nodes =
      find_leaf_nodes(circuit, subproblem_graph);
  // TODO: Find shortest path with minimal weight
  std::vector<std::size_t> minimal_path =
      find_shortest_path(subproblem_graph, leaf_nodes);
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
      schedule = std::pair<std::vector<SingleQubitGateLayer>,
                           std::vector<TwoQubitGateLayer>>{};
  // TODO:ADD in the check that makes sure SQGL's sandwich schedule: If 1st v_p
  // empty skip adding it . If not add empty SQGL
  for (auto i = 0; i < minimal_path.size(); i++) {
    for (auto j = 0;
         j < subproblem_graph.get_Node_Value(minimal_path[i])->second.size();
         ++j) {
      auto op = subproblem_graph.get_Node_Value(minimal_path[i])->second.at(j);
      schedule.first.at(i + 1).emplace_back(
          std::move(*circuit.get_Node_Value(op)));
    }

    for (auto j = 0;
         j < subproblem_graph.get_Node_Value(minimal_path[i])->first.size();
         ++j) {
      auto op = circuit.get_Node_Value(
          subproblem_graph.get_Node_Value(minimal_path[i])->first.at(j));
      schedule.second.at(i).emplace_back(
          std::array<qc::Qubit, 2>({op->get()->getControls().begin()->qubit,
                                    op->get()->getTargets()[0]}));
    }
  }
  return schedule;
}
size_t NativeGateDecomposer::add_node_to_sub_prob_graph(
    std::vector<size_t> v_p, std::vector<size_t> v_c, qc::fp cost,
    DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>&
        subproblem_graph,
    std::size_t prev_node) {
  std::size_t new_node = subproblem_graph.add_Node(std::pair(v_p, v_c));
  subproblem_graph.add_Edge(prev_node, new_node, cost);
  return new_node;
}

double NativeGateDecomposer::schedule_remaining(
    const std::array<std::vector<size_t>, 3>& v,
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
        subproblem_graph,
    size_t prev_node, size_t n_qubits,
    std::map<std::size_t, std::pair<std::size_t, std::array<double, 2>>>&
        memo) {
  double cost;
  // TODO: Check if subproblem has been computed
  std::size_t id = std::hash<std::array<std::vector<size_t>, 3>>{}(v);
  if (memo.contains(id)) {
    std::size_t sub_node = memo.at(id).first;
    double edge_weight = memo.at(id).second[1];
    subproblem_graph.add_Edge(prev_node, sub_node, edge_weight);
    return cost;
  }
  // TODO: Base Case-> V_rem is empty
  if (v[2].size() == 0) {
    cost = v[1].at(0);
    for (auto i = 0; i < v[1].size(); ++i) {
      // DO I need abs here???
      if (v[1].at(i) > cost) {
        cost = v[1].at(i);
      }
    }
    add_node_to_sub_prob_graph(v[0], v[1], cost, subproblem_graph, prev_node);
    return cost;
  }
  // TODO: Recursive Call
  auto v_new = sift(circuit, n_qubits);
  auto args = get_possible_moments(circuit, v[1], v_new);
  qc::fp temp_cost = 0;
  double min_cost = std::numeric_limits<double>::max();
  double min_weight = std::numeric_limits<double>::max();
  std::size_t min_node;
  for (const auto& arg : args) {
    auto new_node = add_node_to_sub_prob_graph(v[0], v_new[1], arg.second,
                                               subproblem_graph, prev_node);
    temp_cost = schedule_remaining(v_new, circuit, subproblem_graph, new_node,
                                   n_qubits, memo) +
                arg.second;
    if (temp_cost < min_cost) {
      min_cost = temp_cost;
      min_node = new_node;
      min_weight = arg.second;
    }
  }
  memo[id] = std::pair<std::size_t, std::array<double, 2>>(
      min_node, {min_cost, min_weight});
  return cost;
}

auto NativeGateDecomposer::schedule(
    const std::pair<std::vector<std::vector<StructU3>>,
                    std::vector<TwoQubitGateLayer>>& asap_schedule) const
    -> std::pair<std::vector<std::vector<StructU3>>,
                 std::vector<TwoQubitGateLayer>> {
  // TODO: Preprocessing-> Convert gates to combined U3 and CZ? Issue with
  // Single qubit layer/SIngle Qubit ReferenceLayer
  std::pair<std::vector<SingleQubitGateLayer>,
            std::vector<SingleQubitGateLayer>>
      qc_p = preprocess(asap_schedule);
  // TODO: Convert Circuit to DAG: How to handle the unique Pointer situation???
  DiGraph<std::unique_ptr<const qc::Operation>> circuit =
      convert_circ_to_dag(qc_p);
  // TODO: Get initial Moments( Nott does MQB THEN SQB!! SOl to get SQB MQB??)
  // v=(v_p,v_c,v_r)

  std::array<std::vector<size_t>, 3> v = sift(circuit, nQubits_);
  // TODO: Create Subproblem Graph
  DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>
      sub_prob_graph = DiGraph<
          std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>();
  // TODO: First Call of Recursive Function to create Schedule
  auto base_node = sub_prob_graph.add_Node(
      std::pair<std::vector<std::size_t>, std::vector<std::size_t>>({}, {}));
  std::map<std::size_t, std::pair<std::size_t, std::array<double, 2>>> memo =
      {};
  auto cost =
      schedule_remaining(v, circuit, sub_prob_graph, base_node, nQubits_, memo);

  // TODO: Create Schedule from Subproblem Graph
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
      final_circuit = build_schedule(circuit, sub_prob_graph);
  std::pair<std::vector<std::vector<StructU3>>, std::vector<TwoQubitGateLayer>>
      postprocessed = postprocess(final_circuit);
  return postprocessed;
}
} // namespace na::zoned
