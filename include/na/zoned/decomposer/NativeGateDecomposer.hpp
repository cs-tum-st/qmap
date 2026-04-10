/*
 * Copyright (c) 2023 - 2026 Chair for Design Automation, TUM
 * Copyright (c) 2025 - 2026 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "DecomposerBase.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/zoned/Compiler.hpp"
#include "na/zoned/Types.hpp"

#include <vector>

namespace na::zoned {

class NativeGateDecomposer : public DecomposerBase {

  /**
   * A quaternion is represented by an array of four `qc::fp` values `{q0, q1,
   * q2, q3}` denoting the components of the quaternion.
   */
  using Quaternion = std::array<qc::fp, 4>;
  size_t nQubits_ = 0;

  constexpr static qc::fp epsilon =
      std::numeric_limits<qc::fp>::epsilon() * 1024;

public:
  /// The configuration of the NativeGateDecomposer
  /// TODO:: ADD Theta_opt option (and check last?))
  struct Config {
    bool theta_opt_schedule = false;
    bool check_final_cond = false;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Config, theta_opt_schedule,
                                                check_final_cond);
  };

  /**
   * A minimal struct to store the parameters of a U3 gate along with the qubit
   * it acts on.
   */
  struct StructU3 {
    std::array<qc::fp, 3> angles;
    qc::Qubit qubit;
  };

private:
  /// The configuration of the NativeGateDecomposer
  Config config_;

public:
  /// Create a new NativeGateDecomposer.
  // TODO: Add in BOOL for scheduling here?
  NativeGateDecomposer(const Architecture& /* unused */,
                       const Config& /* unused */) {}

  /**
   * @brief Converts commonly used single qubit gates into their Quaternion
   * representation.
   * @details A single qubit gate R_v(phi) with rotation axis v=(v0,v1,v2)
   * and rotation angle phi can be represented as a quaternion:
   * @code quaternion(R_v(phi)) = (cos(phi/2) * I, v0 * sin(phi/2) * X, v1 *
   * sin(phi/2) * Y, v2 * sin(phi/2) * Z)@endcode with X, Y, Z Pauli Matrices.
   * @param op a reference_wrapper to the operation to be converted
   * @returns a quaternion.
   */
  static auto
  convertGateToQuaternion(std::reference_wrapper<const qc::Operation> op)
      -> Quaternion;
  /**
   * @brief Merges the quaternions representing two gates as in a matrix
   * multiplication of the gates.
   * @param q1 the first quaternion to be combined.
   * @param q2 the second quaternion to be combined.
   * @returns an quaternion.
   */
  static auto combineQuaternions(const Quaternion& q1, const Quaternion& q2)
      -> Quaternion;
  /**
   * @brief Calculates the values of the U3-gate parameters theta, phi, and
   * lambda.
   * @param quat is a quaternion representing a single qubit gate.
   * @returns an array of three `qc::fp` values `{theta, phi, lambda}` giving
   * the U3 gate angles.
   */
  static auto getU3AnglesFromQuaternion(const Quaternion& quat)
      -> std::array<qc::fp, 3>;

  /**
   * @brief Calculates the largest value of the U3-gate parameter theta from a
   * vector of operations.
   * @param layer is a vector of U3 parameters.
   * @returns the maximal value of theta in the given layer.
   */
  static auto calcThetaMax(const std::vector<StructU3>& layer) -> qc::fp;

  /**
   * @brief Takes a vector of SingleQubitGateLayers and, for each layer,
   * transforms all gates into U3 gates represented by `StructU3` objects.
   * @details It combines all gates acting on the same qubit into a single U3
   * gate.
   * @param layers is a std::vector of SingleQubitGateLayers of a scheduled
   * circuit.
   * @returns a vector of vectors of StructU3 objects representing the single
   * qubit gate layers.
   */
  [[nodiscard]] static auto
  transformToU3(const std::vector<SingleQubitGateRefLayer>& layers,
                size_t n_qubits) -> std::vector<std::vector<StructU3>>;
  /**
   * @brief Takes a vector of `qc::fp` representing the U3-gate angles of a
   * single-qubit gate and the maximal value of theta for the single qubit gate
   * layer and calculates the transversal decomposition angles as in Nottingham
   * et. al. 2024.
   * @param angles is a `std::array` of `qc::fp` representing (theta, phi,
   * lambda).
   * @param theta_max the maximal theta value of the single-qubit gate layer.
   * @returns an array of `qc::fp` values giving the angles (chi, gamma_minus,
   * gamma_plus).
   */
  auto static getDecompositionAngles(const std::array<qc::fp, 3>& angles,
                                     qc::fp theta_max) -> std::array<qc::fp, 3>;

  [[nodiscard]] auto
  decompose(size_t nQubits,
            const std::pair<std::vector<SingleQubitGateRefLayer>,
                            std::vector<TwoQubitGateLayer>>& asap_schedule)
      -> std::pair<std::vector<SingleQubitGateLayer>,
                   std::vector<TwoQubitGateLayer>> override;

  template <class T> class DiGraph {
    std::size_t nodes;
    std::vector<std::vector<std::pair<std::size_t, double>>> adjacencies_;
    std::vector<T> node_values_;

  public:
    DiGraph() {
      nodes = 0;
      adjacencies_ = std::vector<std::vector<std::pair<std::size_t, double>>>();
      node_values_ = std::vector<T>();
    }
    auto add_Node(T node) -> std::size_t {
      adjacencies_.emplace_back();
      node_values_.push_back(std::move(node));
      return nodes++;
    }
    auto add_Edge(std::size_t from, std::size_t to, double weight) -> bool {
      if (from < nodes && to < nodes && from != to) {
        adjacencies_[from].emplace_back(
            std::pair<std::size_t, std::size_t>(to, weight));
        return true;
      } else {
        return false;
      }
    }

    auto add_Edge(std::size_t from, std::size_t to) -> bool {
      if (from < nodes && to < nodes && from != to) {
        adjacencies_[from].emplace_back(
            std::pair<std::size_t, std::size_t>(to, 1.0));
        return true;
      } else {
        return false;
      }
    }

    auto get_Node_Value(std::size_t node) -> T {
      if (node < nodes) {
        return node_values_[node];
        // return node_values_.at(node);
      } else {
        std::ostringstream oss;
        oss << "ERROR: Node Number out of range: " << node << "\n";
        throw std::invalid_argument(oss.str());
      }
    }
    [[nodiscard]] auto get_N_Nodes() const -> std::size_t { return nodes; }

    [[nodiscard]] auto get_adjacent(std::size_t i) const
        -> std::vector<std::pair<std::size_t, double>> {
      return adjacencies_.at(i);
    }
  };

  static auto
  convert_circ_to_dag(const std::pair<std::vector<std::vector<StructU3>>,
                                      std::vector<TwoQubitGateLayer>>& qc,
                      size_t n_qubits)
      -> DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>;

  static auto shortest_path_to_start(
      const DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      unsigned long long curr_node, const std::set<size_t>& leaf_nodes)
      -> std::pair<std::vector<size_t>, double>;
  static auto find_shortest_path(
      const DiGraph<std::pair<std::vector<std::size_t>,
                              std::vector<std::size_t>>>& di_graph,
      const std::vector<std::size_t>& vector) -> std::vector<size_t>;
  static auto
  calc_cost(std::vector<std::vector<size_t>>::const_reference path,
            const DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
                subproblem_graph) -> double;
  static auto find_leaf_nodes(
      const DiGraph<std::pair<std::vector<std::size_t>,
                              std::vector<std::size_t>>>& subproblem_graph)
      -> std::vector<std::size_t>;
  static auto remove_element(const std::vector<std::size_t>& vector,
                             std::size_t node) -> std::vector<std::size_t>;
  static auto get_possible_moments(
      DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>& circuit,
      const std::vector<size_t>& v0_c,
      const std::array<std::vector<size_t>, 3>& v_new, bool check_final_cond)
      -> std::vector<std::pair<std::array<std::vector<size_t>, 4>, qc::fp>>;
  static auto
  max_theta(DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>& circuit,
            const std::vector<unsigned long long>& nodes) -> qc::fp;
  static auto
  sift(DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>& circuit,
       size_t n_qubits) -> std::array<std::vector<size_t>, 3>;

  static auto build_schedule(
      DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>& circuit,
      DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph) -> std::pair<std::vector<std::vector<StructU3>>,
                                         std::vector<TwoQubitGateLayer>>;

  static auto add_node_to_sub_prob_graph(
      const std::vector<size_t>& v_p, const std::vector<size_t>& v_c,
      qc::fp cost,
      DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      size_t prev_node) -> size_t;

  static auto schedule_remaining(
      const std::array<std::vector<size_t>, 3>& v,
      DiGraph<std::variant<StructU3, std::array<qc::Qubit, 2>>>& circuit,
      DiGraph<std::pair<std::vector<size_t>, std::vector<size_t>>>&
          subproblem_graph,
      size_t prev_node, size_t n_qubits, bool check_final_cond,
      std::map<size_t, std::pair<size_t, std::array<double, 2>>>& memo)
      -> double;
  /**
   * This function schedules the operations of a quantum computation.
   * @details
   * @param asap_schedule
   * @param qc is the quantum computation
   * @return a pair of two vectors. The first vector contains the layers of
   * single-qubit operations. The second vector contains the layers of two-qubit
   * operations. A pair of qubits represents every two-qubit operation.
   */
  [[nodiscard]] auto schedule_theta_opt(
      const std::pair<std::vector<std::vector<StructU3>>,
                      std::vector<TwoQubitGateLayer>>& asap_schedule,
      std::size_t n_qubits) const
      -> std::pair<std::vector<std::vector<StructU3>>,
                   std::vector<TwoQubitGateLayer>>;
};
} // namespace na::zoned

template <> struct std::hash<std::array<std::vector<std::size_t>, 3>> {
  auto operator()(const std::array<std::vector<std::size_t>, 3>& array)
      const noexcept -> std::size_t {
    std::size_t seed = 0U;
    for (const auto& v : array) {
      for (auto node : v) {
        qc::hashCombine(seed, std::hash<std::size_t>{}(node));
      }
    }
    return seed;
  }
};
