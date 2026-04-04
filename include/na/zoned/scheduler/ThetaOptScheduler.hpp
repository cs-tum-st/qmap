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
    // Add pedecessors???
    std::vector<std::vector<std::pair<std::size_t, double>>> adjacencies_;
    std::vector<T> node_values_;

  public:
    DiGraph() {
      nodes = 0;
      adjacencies_ = std::vector<std::vector<std::pair<std::size_t, double>>>();
      node_values_ = std::vector<T>();
    }
    std::size_t add_Node(T node) {
      adjacencies_.push_back(std::vector<std::pair<std::size_t, double>>());
      node_values_.push_back(std::move(node));
      return nodes++;
    }
    bool add_Edge(std::size_t from, std::size_t to, double weight) {
      if (from < nodes && to < nodes && from != to) {
        adjacencies_[from].push_back(
            std::pair<std::size_t, std::size_t>(to, weight));
        return true;
      } else {
        return false;
      }
    }

    bool add_Edge(std::size_t from, std::size_t to) {
      if (from < nodes && to < nodes && from != to) {
        adjacencies_[from].push_back(
            std::pair<std::size_t, std::size_t>(to, 1.0));
        return true;
      } else {
        return false;
      }
    }

    T* get_Node_Value(std::size_t node) {
      if (node < nodes) {
        return &node_values_[node];
        // return node_values_.at(node);
      } else {
        std::ostringstream oss;
        oss << "ERROR: Node Number out of range: " << node << "\n";
        throw std::invalid_argument(oss.str());
      }
    }
    std::size_t get_N_Nodes() const { return nodes; }

    std::vector<std::pair<std::size_t, double>>
    get_adjacent(std::size_t i) const {
      return adjacencies_.at(i);
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
  static DiGraph<std::unique_ptr<const qc::Operation>>
  convert_circ_to_dag(const std::pair<std::vector<SingleQubitGateLayer>,
                                      std::vector<TwoQubitGateLayer>>& qc);
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
  preprocess_qc(const qc::QuantumComputation& qc,
                const Architecture& architecture, Config config) const;
  static std::pair<std::vector<size_t>, double> shortest_path_to_start(
      const DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      unsigned long long curr_node, std::set<size_t> leaf_nodes);
  std::vector<size_t> find_shortest_path(
      const DiGraph<std::pair<std::vector<std::size_t>,
                              std::vector<std::size_t>>>& di_graph,
      const std::vector<std::size_t>& vector);
  static double
  calc_cost(std::vector<std::vector<size_t>>::const_reference path,
            DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
                subproblem_graph);
  static std::vector<std::size_t> find_leaf_nodes(
      DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
      const DiGraph<std::pair<std::vector<std::size_t>,
                              std::vector<std::size_t>>>& subproblem_graph);
  static std::vector<std::size_t>
  sort_by_theta_dec(DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
                    const std::vector<unsigned long long>& vector);
  static std::vector<std::size_t>
  remove_element(const std::vector<std::size_t>& vector, std::size_t node);
  static std::vector<std::pair<std::array<std::vector<size_t>, 4>, qc::fp>>
  get_possible_moments(DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
                       const std::vector<size_t>& v0_c,
                       const std::array<std::vector<size_t>, 3>& v_new);
  static qc::fp
  max_theta(DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
            const std::vector<unsigned long long>& nodes);
  static std::array<std::vector<size_t>, 3>
  sift(DiGraph<std::unique_ptr<const qc::Operation>>& circuit, size_t n_qubits);
  std::pair<std::vector<SingleQubitGateLayer>, std::vector<TwoQubitGateLayer>>
  build_schedule(DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
                 DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
                     subproblem_graph) const;

  static auto add_node_to_sub_prob_graph(
      std::vector<size_t> v_p, std::vector<size_t> v_c, qc::fp cost,
      DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      size_t prev_node) -> size_t;

  static double schedule_remaining(
      const std::array<std::vector<size_t>, 3>& v,
      DiGraph<std::unique_ptr<const qc::Operation>>& circuit,
      DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      size_t prev_node, size_t n_qubits,
      std::map<size_t, std::pair<size_t, std::array<double, 2>>>& memo);
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

template <> struct std::hash<std::array<std::vector<std::size_t>, 3>> {
  std::size_t operator()(
      const std::array<std::vector<std::size_t>, 3>& array) const noexcept {
    std::size_t seed = 0U;
    for (auto v : array) {
      for (auto node : v) {
        qc::hashCombine(seed, std::hash<std::size_t>{}(node));
      }
    }
    return seed;
  }
};
