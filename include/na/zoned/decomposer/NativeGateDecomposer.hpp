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
#include "na/zoned/Types.hpp"

#include <vector>

namespace na::zoned {

class NativeGateDecomposer : public DecomposerBase {
  struct StructU3 {
    std::array<qc::fp, 3> angles;
    qc::Qubit qubit;
  };
  size_t nQubits_;

  constexpr static qc::fp epsilon =
      std::numeric_limits<qc::fp>::epsilon() * 1024;

public:
  /**
   * Converts commonly used single qubit gates into their Quaternion
   * representation
   * quaternion(R_v(phi))=(cos(phi/2)*I,v0*sin(phi/2)*X,v1*sin(phi/2)*Y,v2*sin(phi/2)*Z)
   * with X,Y,Z Pauli Matrices
   * @param op a reference_wrapper to the operation to be converted
   * @returns an array of four qc::fp values [q0,q1,q2,q3] denoting the
   * components of the quaternion
   */
  static auto
  convert_gate_to_quaternion(std::reference_wrapper<const qc::Operation> op)
      -> std::array<qc::fp, 4>;
  /**
   * Merges the quaternions representing two gates as in a Matrix multiplication
   * of the gates
   * @param q1 the first Quaternion to be combined
   * @param q2 the second Quaternion to be combined
   * @returns an array of four qc::fp values [q0,q1,q2,q3] denoting the
   * components of the combined quaternion
   */
  static auto combine_quaternions(const std::array<qc::fp, 4>& q1,
                                  const std::array<qc::fp, 4>& q2)
      -> std::array<qc::fp, 4>;
  /**
   * Calculates the values of the U3-Gate parameters theta, phi and lambda
   * @param quat is a quaternion representing a single qubit gate
   * @returns an array of three qc::fp values [theta, phi, lambda] giving the U§
   * gate angles
   */
  static auto get_U3_angles_from_quaternion(const std::array<qc::fp, 4>& quat)
      -> std::array<qc::fp, 3>;

  /**
   * Calculates the largest value of the U3-Gate parameter theta from a vector
   * of Operations. Fails when provided gates aren't all U3-Gates.
   * @param layer is a SingleQubitGateLayers of a scheduled
   * @returns the maximal value of theta in the given layer
   */
  static auto calc_theta_max(const std::vector<struct_U3>& layer) -> qc::fp;

  /**
   * Takes a vector of SingleQubitGateLayer's and, for each layer
   * , transforms all gates into U3 gates represented by struct_U3's
   * , combining  all gates acting on the same qubit into a single U3 gate
   * @param layers is a std::vector of SingleQubitGateLayers of a scheduled
   * circuit
   * @returns a vector of vectors of a struct_U3's representing the single qubit
   * gate layers
   */
  [[nodiscard]] auto
  transform_to_U3(const std::vector<SingleQubitGateRefLayer>& layers) const
      -> std::vector<std::vector<struct_U3>>;
  /**
   * Takes a vector of qc::fp's representing the U3-Gate angles of a single
   * qubit hate and the maximal value of theta for the single qubit gate layer
   * and calculates the transversal decomposition angles as in Nottingham et.
   * al. 2024
   * @param angles is a std::array of qc::fp representing (theta,phi, lambda)
   * @param theta_max the maximal theta value of the single-qubit qate layer
   * @returns an array of qc::fp values giving the angles (chi, gamma_minus,
   * gamma_plus)
   */
  auto static get_decomposition_angles(const std::array<qc::fp, 3>& angles,
                                       qc::fp theta_max)
      -> std::array<qc::fp, 3>;
  /**
   * Create a new Decomposer.
   * @param n_qubits is the number of qubits in the circuit to be decomposed
   */
  NativeGateDecomposer(int n_qubits);

  [[nodiscard]] auto decompose(
      const std::vector<SingleQubitGateRefLayer>& singleQubitGateLayers) const
      -> std::vector<SingleQubitGateLayer> override;
};

} // namespace na::zoned
