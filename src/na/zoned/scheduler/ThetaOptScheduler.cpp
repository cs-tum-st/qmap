//
// Created by cpsch on 12.03.2026.
//
#include "na/zoned/scheduler/ThetaOptScheduler.hpp"

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
                                 const ASAPScheduler::Config config) const {
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
      NewLayer.emplace_back(std::make_unique<const qc::StandardOperation>(
          qc::StandardOperation((qc::Qubit) gate.qubit, qc::U,{gate.angles[0],gate.angles[1],gate.angles[2]})));
    }
    NewSingleQubitLayers.push_back(std::move((NewLayer)));
  }
  return std::pair(std::move(NewSingleQubitLayers), asap_schedule.second);
}
std::vector<std::pair<std::array<std::vector<std::size_t>,4>,qc::fp>> ThetaOptScheduler::get_possible_moments(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<size_t>& v0_c,
    const std::array<std::vector<size_t>, 3>& v_new) {

    std::vector<size_t> v_p1_star=v_new[0];
    std::vector<size_t> v_p1_square=std::vector<size_t>();
    std::vector<size_t> v_c1_star=v_new[1];
    std::vector<size_t> v_c1_square=std::vector<size_t>();

    qc::fp v_c0_cost=max_theta(circuit, v0_c);
    qc::fp v_c1_cost=max_theta(circuit, v_new[1]);
    qc::fp orig_comb_cost=v_c0_cost+v_c1_cost;
    qc::fp new_v_C1_cost=std::max(v_c0_cost,v_c1_cost);

    std::array<std::vector<std::size_t>,4> v_arg={v0_c,v_new[0],v_new[1],v_new[2]};
    auto args={std::pair(v_arg,v_c0_cost)};
    //TODO: Check Condition 1
      //Sort v_0C from highest to lowest theta

    //TODO: Check Condition 2
    // TODO: Check Condition 3
    //TODO Check Condition 4 ?? Find out what it does!!

    return args;
}

ThetaOptScheduler::DiGraph<std::unique_ptr<const qc::Operation>>
ThetaOptScheduler::convert_circ_to_dag(
    const std::pair<std::vector<SingleQubitGateLayer>,
                    std::vector<TwoQubitGateLayer>>& qc) const {
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
      auto t=qc.second.at(i);
      //Create CZ & Make_unique
      size_t node = graph.add_Node(std::make_unique<qc::StandardOperation>(qc::StandardOperation(qc::Control(qc.second.at(i).at(m)[0]),qc.second.at(i).at(m)[1],qc::Z,{})));
      qubit_paths.at(qc.second.at(i).at(m).at(0))
            .push_back(node);
      qubit_paths.at(qc.second.at(i).at(m).at(1))
            .push_back(node);

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


bool disjunct(const std::unordered_set<size_t>& set1,
              const std::unordered_set<size_t>& set2) {
  for (auto elem: set1) {
    if (set2.contains(elem)){return false;}
  }
  return true;
}

qc::fp ThetaOptScheduler::max_theta(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    const std::vector<unsigned long long>& nodes) {
    qc::fp max_cost=0;
    for (const auto node : nodes) {
      if ((*circuit.get_Node_Value(node))->getParameter()[0]>=max_cost) {
        max_cost=(*circuit.get_Node_Value(node))->getParameter()[0];
      }
    }
    return max_cost;
}

std::array<std::vector<size_t>, 3> ThetaOptScheduler::sift(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    size_t n_qubits) {
  std::vector<size_t> v_p=std::vector<size_t>();
  std::vector<size_t> v_c=std::vector<size_t>();
  std::vector<size_t> v_r=std::vector<size_t>();

  std::unordered_set<size_t> removed=std::unordered_set<size_t>();

  for (auto node=0; node < circuit.get_N_Nodes(); node++) {
      //TODO: Check for unique pointer issue
    auto op=circuit.get_Node_Value(node);
    std::unordered_set<size_t> op_qubits=std::unordered_set<size_t>();
    for (auto qubit:op->get()->getUsedQubits()){op_qubits.insert(qubit);}
      if (removed.size()<n_qubits && disjunct(removed,op_qubits)) {
        if (op->get()->getNqubits()==1) {
          v_c.push_back(node);
          removed.insert(op->get()->getTargets().at(0));
        }else {
          v_p.push_back(node);
          //Add something to make it only pick one 2-Qubit gate per Qubit per moment???
        }
      }else {
        v_r.push_back(node);
      }
  }
  return std::array<std::vector<size_t>, 3>({v_p,v_c,v_r});
}

static std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
ThetaOptScheduler::build_schedule(
    DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
    DiGraph<std::pair<std::vector<std::size_t>,std::vector<std::size_t>>>&
        subproblem_graph) {

    //TODO: Find Leaf Nodes of di_graph
    std::vector<std::size_t> leaf_nodes;
    for (std::size_t i=0;i<subproblem_graph.get_N_Nodes();i++) {
      if (subproblem_graph.get_adjacent(i).empty()) {
        leaf_nodes.push_back(i);
      }
    }
   //TODO: Find shortest path to leaf nodes
  std::vector<std::size_t> start_nodes= find_start_nodes(subproblem_graph);
  std::vector<std::vector<std::size_t>> paths_to_leafs= find_shortest_paths(subproblem_graph,leaf_nodes,start_nodes);
  //TODO: Find shortest path with minimal weight
  std::vector<std::size_t> minimal_path=std::vector<std::size_t>();
  double min_cost= 100000000000000; //TODO: Find Max double
  for (auto path: paths_to_leafs) {
    double cost=calc_cost(subproblem_graph,cost);
    if (cost<min_cost) {min_cost=cost;}
  }
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>> schedule=std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>{};
  //TODO:ADD in the check that makes sure SQGL's sandwich schedule: If 1st v_p empty skip adding it . If not add empty SQGL
  for (auto i=0;i<minimal_path.size(); i++) {
    for (auto j=0;j<subproblem_graph.get_Node_Value(minimal_path[i])->second.size();++j) {
      auto op=subproblem_graph.get_Node_Value(minimal_path[i])->second.at(j);
      schedule.first.at(i+1).emplace_back(std::move(*circuit.get_Node_Value(op)));
    }

    for (auto j=0;j<subproblem_graph.get_Node_Value(minimal_path[i])->first.size();++j) {
      auto op=subproblem_graph.get_Node_Value(minimal_path[i])->first.at(j);
      schedule.second.at(i).emplace_back(std::move(*circuit.get_Node_Value(op)));
    }
  }
  return schedule;
}
size_t ThetaOptScheduler::add_node_to_sub_prob_graph(
    std::vector<size_t> v_p, std::vector<size_t> v_c, qc::fp cost,
    DiGraph<std::pair<std::vector<std::size_t>,std::vector<std::size_t>>>&
        subproblem_graph,std::size_t prev_node) {
    std::size_t new_node=subproblem_graph.add_Node(std::pair(v_p,v_c));
    subproblem_graph.add_Edge(prev_node,new_node,cost);
    return new_node;
}
double ThetaOptScheduler::schedule_remaining(
    const std::array<std::vector<size_t>, 3>& v,
     DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
     DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
        subproblem_graph,
    size_t prev_node, size_t n_qubits) {
    double cost;
    //TODO: Check if subproblem has been computed
    //auto id=?;//What to do for id??
    //if (memo.contains()) {//KEy set??
    //  cost=memo.at(id)[0];
    //  sub_node=memo.at(id)[1];
    //  edge_weight=memo.at(id)[2];
    //  subproblem_graph.add_edge(prev_node,sub_node,edge_weight);
    //   return cost;
    //  }
     //TODO: Base Case-> V_rem is empty
      if (v[2].size()==0) {
        cost=0;
        for (auto i=0; i <v[1].size();++i) {
          //DO I need abs here??? Should I set base value to - PI??
          if (std::fabs(v[1].at(i))>cost) {
            cost=std::fabs(v[1].at(i));
          }
        }
          add_node_to_sub_prob_graph(v[0],v[1],cost,subproblem_graph, prev_node);
        return cost;
      }
      //TODO: Recursive Call
      auto v_new=sift(circuit,n_qubits);
      auto args= get_possible_moments(circuit,v[1],v_new);
      std::vector<double> costs=std::vector<double>();
      for (const auto& arg:args) {
       auto new_node=add_node_to_sub_prob_graph(v[0],v_new[1],arg.second,subproblem_graph, prev_node);
       costs.push_back(schedule_remaining(v_new,circuit,subproblem_graph,new_node,n_qubits)+arg.second);
      }
      cost= *std::ranges::min_element(costs);
      //memo.add(id,)
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
  //TODO: Circuit on NON_Unique pointers??
  DiGraph<std::unique_ptr<const qc::Operation>> circuit =
      convert_circ_to_dag(qc_p);
  // TODO: Get initial Moments( Nott does MQB THEN SQB!! SOl to get SQB MQB??)
  // v=(v_p,v_c,v_r)
  std::array<std::vector<size_t>, 3> v=sift(circuit,qc.size());
  // TODO: Create Subproblem Graph
  DiGraph<std::pair<std::vector<std::size_t>,
  std::vector<std::size_t>>>sub_prob_graph= DiGraph<std::pair<std::vector<std::size_t>,
  std::vector<std::size_t>>> ();
  // TODO: First Call of Recursive Function to create Schedule
  auto base_node=sub_prob_graph.add_Node(std::pair<std::vector<std::size_t>,std::vector<std::size_t>>({},{}));
  auto cost=schedule_remaining(v, circuit, sub_prob_graph, base_node,qc.getNqubits());

  // TODO: Create Schedule from Subproblem Graph
  std::pair<std::vector<SingleQubitGateRefLayer>,
                 std::vector<TwoQubitGateLayer>> final_circuit=build_schedule(circuit,sub_prob_graph);
  return final_circuit;
}

} // namespace na::zoned