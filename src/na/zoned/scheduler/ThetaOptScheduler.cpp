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
          qc::StandardOperation::StandardOperation(gate.qubit, qc::U,
                                                   gate.angles)));
    }
    NewSingleQubitLayers.emplace_back(NewLayer);
  }
  return std::pair(NewSingleQubitLayers, asap_schedule.second);
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
      size_t node = graph.add_Node(qc.first.at(i).at(s));

      for (auto t = 0; t < qc.first.at(i).at(s)->getNtargets(); t++) {
        qubit_paths.at(qc.first.at(i).at(s)->getTargets().at(t))
            .push_back(node);
      }
    }
    for (auto m = 0; m < qc.second.at(i).size(); ++m) {
      size_t node = graph.add_Node(qc.second.at(i).at(m));

      for (auto t = 0; t < qc.first.at(i).at(m)->getNtargets(); t++) {
        qubit_paths.at(qc.first.at(i).at(m)->getTargets().at(t))
            .push_back(node);
      }
      for (auto c = 0; c < qc.first.at(i).at(m)->getNcontrols(); c++) {
        qubit_paths.at(qc.first.at(i).at(m)->getControls().at(m))
            .push_back(node);
      }
    }
  }
  for (auto s = 0; s < qc.first.end()->size(); ++s) {
    size_t node = graph.add_Node(qc.first.end()->at(s));

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

auto ThetaOptScheduler::schedule(const qc::QuantumComputation& qc) const
    -> std::pair<std::vector<SingleQubitGateRefLayer>,
                 std::vector<TwoQubitGateLayer>> {
  // TODO: Preprocessing-> Convert gates to combined U3 and CZ? Issue with
  // Single qubit layer/SIngle Qubit ReferenceLayer
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
      qc_p = preprocess_qc(qc, architecture_, config_);
  // TODO: Convert Circuit to DAG: How to handle the unique Pointer situation???
  DiGraph<std::unique_ptr<const qc::Operation>> circuit =
      convert_circ_to_dag(qc_p);
  // TODO: Get initial Moments( Nott does MQB THEN SQB!! SOl to get SQB MQB??)
  // SIFt IMPL?

  // TODO: Create Subproblem Graph

  // TODO: First Call of Recursive Function to create Schedule

  // TODO: Create Schedule from Subproblem Graph

  return std::pair(std::vector<SingleQubitGateRefLayer>(),
                   std::vector<TwoQubitGateLayer>());
}

} // namespace na::zoned