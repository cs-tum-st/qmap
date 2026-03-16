//
// Created by cpsch on 12.03.2026.
//
#pragma once

#include "ASAPScheduler.hpp"
#include "ir/QuantumComputation.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/Types.hpp"
#include "na/zoned/scheduler/SchedulerBase.hpp"

#include <functional>
#include <utility>
#include <vector>

namespace na::zoned {
/**
 * The class ThetaOptScheduler implements a scheduling strategy minimizing the
 * total global rotation angle theta over the circuit for the zoned neutral atom
 * compiler.
 */
class ThetaOptScheduler : public SchedulerBase {
  /// A reference to the zoned neutral atom architecture
  std::reference_wrapper<const Architecture> architecture_;
  /**
   * This value is calculated based on the architecture and indicates the
   * the entanglement zone.
   */
  size_t maxTwoQubitGateNumPerLayer_ = 0;

public:
  /// The configuration of the ASAPScheduler
  struct Config {
    /// The maximal share of traps that are used in the entanglement zone.
    double maxFillingFactor = 0.9;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Config, maxFillingFactor);
  };

  template <class T> class DiGraph {
    std::size_t nodes;
    std::vector<std::vector<std::size_t>> adjacencies_;
    std::vector<T> node_values_;

  public:
    DiGraph() {
      nodes = 0;
      adjacencies_ = std::vector<std::vector<std::size_t>>();
      node_values_ = std::vector<T>();
    }
    std::size_t add_Node(T node) {
      adjacencies_.push_back(std::vector<std::size_t>());
      node_values_.push_back(node);
      return nodes++;
    }
    bool add_Edge(std::size_t from, std::size_t to) {
      // TODO: Cycle Check
      if (from < nodes && to < nodes && from != to) {
        adjacencies_[from].push_back(to);
        return true;
      } else {
        return false;
      }
    }
    T get_Node_Value(std::size_t node) {
      if (node > nodes) {
        return node_values_.at(node);
      } else {
        std::ostringstream oss;
        oss << "ERROR: Node Number out of range: " << node << "\n";
        throw std::invalid_argument(oss.str());
      }
    }
    std::size_t get_N_Nodes() const { return nodes; }
    DiGraph<T> topographicSort(DiGraph<T> graph) {
      // TODO: Topographic Sort
      return graph;
    }
  };

private:
  /// The configuration of the ThetaOptScheduler
  Config config_;

public:
  /**
   * Create a new ThetaOptScheduler.
   * @param architecture is the architecture of the neutral atom system
   * @param config is the configuration for the scheduler
   */
  ThetaOptScheduler(const Architecture& architecture, const Config& config);
  DiGraph<std::unique_ptr<const qc::Operation>> convert_circ_to_dag(
      const std::pair<std::vector<SingleQubitGateLayer>,
                      std::vector<TwoQubitGateLayer>>& qc) const;
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
  preprocess_qc(const qc::QuantumComputation& qc,
                const Architecture& architecture,
                ASAPScheduler::Config config) const;
  /**
   * This function schedules the operations of a quantum computation.
   * @details
   * @param qc is the quantum computation
   * @return a pair of two vectors. The first vector contains the layers of
   * single-qubit operations. The second vector contains the layers of two-qubit
   * operations. A pair of qubits represents every two-qubit operation.
   */
  [[nodiscard]] auto schedule(const qc::QuantumComputation& qc) const
      -> std::pair<std::vector<SingleQubitGateRefLayer>,
                   std::vector<TwoQubitGateLayer>> override;
};
} // namespace na::zoned
