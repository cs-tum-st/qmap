//
// Created by cpsch on 12.03.2026.
//
#include "na/zoned/scheduler/ThetaOptScheduler.hpp"

#include "dd/Approximation.hpp"
#include "ir/Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/decomposer/NativeGateDecomposer.hpp"
#include "na/zoned/scheduler/ASAPScheduler.hpp"
#include "spdlog/fmt/bundled/args.h"
#include "spdlog/fmt/bundled/printf.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na::zoned {

ThetaOptScheduler::ThetaOptScheduler(const Architecture& architecture,
                                     const Config& config)
    : architecture_(architecture), config_(config) {}

std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
ThetaOptScheduler::preprocess_qc(const qc::QuantumComputation& qc,
                                 const Architecture& architecture,
                                 const Config config) const {
  // TODO: Find out if only CZ's in Mulit qubit gate layers
  ASAPScheduler scheduler = ASAPScheduler(architecture, config);
  std::pair<std::vector<SingleQubitGateRefLayer>,
            std::vector<TwoQubitGateLayer>>
      asap_schedule = scheduler.schedule(qc);
  std::vector<std::vector<NativeGateDecomposer::StructU3>>
      singlequbitlayers_U3 =
          NativeGateDecomposer::transformToU3(asap_schedule.first);
  std::vector<SingleQubitGateLayer> NewSingleQubitLayers =
      std::vector<SingleQubitGateLayer>{};
  SingleQubitGateLayer NewLayer;
  for (auto layer : singlequbitlayers_U3) {
    NewLayer.clear();
    for (auto gate : layer) {
      NewLayer.emplace_back(
          std::make_unique<const qc::StandardOperation>(qc::StandardOperation(
              (qc::Qubit)gate.qubit, qc::U,
              {gate.angles[0], gate.angles[1], gate.angles[2]})));
    }
    NewSingleQubitLayers.push_back(std::move((NewLayer)));
  }
  return std::pair(std::move(NewSingleQubitLayers), asap_schedule.second);
}
std::vector<size_t> ThetaOptScheduler::find_shortest_path(
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
ThetaOptScheduler::shortest_path_to_start(
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
double ThetaOptScheduler::calc_cost(
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
std::vector<std::size_t> ThetaOptScheduler::find_leaf_nodes(
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
std::vector<std::size_t> ThetaOptScheduler::sort_by_theta_dec(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<unsigned long long>& vector) {
  std::vector<std::size_t> sorted = {};
}
std::vector<std::size_t>
ThetaOptScheduler::remove_element(const std::vector<std::size_t>& vector,
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
ThetaOptScheduler::get_possible_moments(
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
      v_arg = {pot.first[0], v_p1_star,pot.first[1],v_new[2]};
      for (auto node: v_c1_star) {
        v_arg[2].push_back(node);
      }
      for (auto node : v_c1_square) {v_arg[3].push_back(node);}
      for (auto node:v_p1_square) {v_arg[3].push_back(node);}
      args.push_back(std::pair<std::array<std::vector<std::size_t>, 4>, qc::fp>(v_arg,pot.second.second));
    }
  }
  return args;
}

ThetaOptScheduler::DiGraph<std::unique_ptr<const qc::Operation>>
ThetaOptScheduler::convert_circ_to_dag(
    const std::pair<std::vector<SingleQubitGateLayer>,
                    std::vector<TwoQubitGateLayer>>& qc) {
  ThetaOptScheduler::DiGraph<std::unique_ptr<const qc::Operation>> graph =
      DiGraph<std::unique_ptr<const qc::Operation>>();
  std::vector<std::vector<size_t>> qubit_paths =
      std::vector<std::vector<size_t>>{};
  // assert that One mor sql exists than mql
  for (auto i = 0; i < qc.first.size(); ++i) {
    for (auto s = 0; s < qc.first.at(i).size(); ++s) {
      size_t node = graph.add_Node(qc.first.at(i).at(s)->clone());
      for (auto t = 0; t < qc.first.at(i).at(s)->getNtargets(); t++) {
        qubit_paths.at(qc.first.at(i).at(s)->getTargets().at(t))
            .push_back(node);
      }
    }
    for (auto m = 0; m < qc.second.at(i).size(); ++m) {
      auto t = qc.second.at(i);
      // Create CZ & Make_unique
      size_t node = graph.add_Node(std::make_unique<qc::StandardOperation>(
          qc::StandardOperation(qc::Control(qc.second.at(i).at(m)[0]),
                                qc.second.at(i).at(m)[1], qc::Z, {})));
      qubit_paths.at(qc.second.at(i).at(m).at(0))
            .push_back(node);
      qubit_paths.at(qc.second.at(i).at(m).at(1)).push_back(node);
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

qc::fp ThetaOptScheduler::max_theta(
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

std::array<std::vector<size_t>, 3>
ThetaOptScheduler::sift(DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
                        size_t n_qubits) {
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

std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
ThetaOptScheduler::build_schedule(
    // TODO: FIgure ouot how to interface with scheduler
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>&
        subproblem_graph) {

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
          std::array<qc::Qubit, 2>({op->get()->getControls().begin()->qubit, op->get()->getTargets()[0]}));
    }
  }
  return schedule;
}
size_t ThetaOptScheduler::add_node_to_sub_prob_graph(
    std::vector<size_t> v_p, std::vector<size_t> v_c, qc::fp cost,
    DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>&
        subproblem_graph,
    std::size_t prev_node) {
  std::size_t new_node = subproblem_graph.add_Node(std::pair(v_p, v_c));
  subproblem_graph.add_Edge(prev_node, new_node, cost);
  return new_node;
}

double ThetaOptScheduler::schedule_remaining(
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
    temp_cost=schedule_remaining(v_new, circuit, subproblem_graph,
                                       new_node, n_qubits, memo) +arg.second;
    if (temp_cost<min_cost) {
      min_cost=temp_cost;
      min_node=new_node;
      min_weight=arg.second;
    }
  }
  memo[id]=std::pair<std::size_t,std::array<double,2>>(min_node,{min_cost,min_weight});
  return cost;
}

auto ThetaOptScheduler::schedule(const qc::QuantumComputation& qc) const
    -> std::pair<std::vector<SingleQubitGateRefLayer>,
                 std::vector<TwoQubitGateLayer>> {
  // TODO: Preprocessing-> Convert gates to combined U3 and CZ? Issue with
  // Single qubit layer/SIngle Qubit ReferenceLayer
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
      qc_p = preprocess_qc(qc, architecture_, config_);
  // TODO: Convert Circuit to DAG: How to handle the unique Pointer situation???
  // TODO: Circuit on NON_Unique pointers??
  DiGraph<std::unique_ptr<const qc::Operation>> circuit =
      convert_circ_to_dag(qc_p);
  // TODO: Get initial Moments( Nott does MQB THEN SQB!! SOl to get SQB MQB??)
  // v=(v_p,v_c,v_r)
  std::array<std::vector<size_t>, 3> v = sift(circuit, qc.size());
  // TODO: Create Subproblem Graph
  DiGraph<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>sub_prob_graph= DiGraph<
          std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>();
  // TODO: First Call of Recursive Function to create Schedule
  auto base_node = sub_prob_graph.add_Node(std::pair<std::vector<std::size_t>,std::vector<std::size_t>>({},{}));
  std::map<std::size_t,std::pair<std::size_t,std::array<double,2>>> memo={};
  auto cost=schedule_remaining(v, circuit, sub_prob_graph, base_node,
                                 qc.getNqubits(), memo);

  // TODO: Create Schedule from Subproblem Graph
  std::pair<std::vector<SingleQubitGateRefLayer>,
                 std::vector<TwoQubitGateLayer>> final_circuit=build_schedule(circuit,sub_prob_graph);
  return final_circuit;
}

} // namespace na::zoned