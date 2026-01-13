#pragma once

#include "ir/operations/StandardOperation.hpp"
#include "na/zoned/Types.hpp"

#include <vector>

namespace na::zoned {

  class Decomposer {
  public:
    /**
     * Converts commonly used single qubit gates into their Quaternion representation
     * quaternion(R_v(phi))=(cos(phi/2)*I,v0*sin(phi/2)*X,v1*sin(phi/2)*Y,v2*sin(phi/2)*Z) with X,Y,Z Pauli Matrices
     * @param op a reference_wrapper to the operation to be converted
     * @returns an array of four qc::fp values [q0,q1,q2,q3] denoting the components of the quaternion
     */
    static auto convert_gate_to_quaternion(std::reference_wrapper<const qc::Operation> op) -> std::array<qc::fp,4>;
    /**
     * Merges the quaternions representing two gates as in a Matrix multiplication of the gates
     * @param q1 the first Quaternion to be combined
     * @param q2 the second Quaternion to be combined
     * @returns an array of four qc::fp values [q0,q1,q2,q3] denoting the components of the combined quaternion
     */
    static auto combine_quaternions(const std::array<qc::fp,4>& q1,
                         const std::array<qc::fp,4>& q2) -> std::array<qc::fp,4>;
    /**
     * Calculates the values of the U3-Gate parameters theta, phi and lambda
     * @param quat is a quaternion representing a single qubit gate
     * @returns an array of three qc::fp values [theta, phi, lambda] giving the U§ gate angles
     */
    static auto get_U3_angles_from_quaternion(std::array<qc::fp,4> quat) -> std::array<qc::fp,3>;

    /**
     * Calculates the largest value of the U3-Gate parameter theta from a vector of Operations.
     * Fails when provided gates aren't all U3-Gates.
     * @param layer is a SingleQubitGateLayers of a scheduled
     */
    static auto calc_theta_max(const SingleQubitGateLayer& layer)->qc::fp;

    /**
     * Takes a vector of SingleQubitGateLayer's and, for each layer
     * , transforms all gates into U3 gates
     * , combining  all gates acting on the same qubit into a single U3 gate
     * @param layers is a std::vector of SingleQubitGateLayers of a scheduled circuit
     */
    [[nodiscard]] auto transform_to_U3(const std::vector<SingleQubitGateLayer>& layers) const
        -> std::vector<SingleQubitGateLayer>;

    /**
     * Create a new Decomposer.
     * @param n_qubits is the number of qubits in the circuit to be decomposed
     */
    Decomposer(int n_qubits);

    [[nodiscard]] auto
    decompose(const std::pair<std::vector<SingleQubitGateLayer>,
                              std::vector<TwoQubitGateLayer>>& schedule)
        -> std::pair<std::vector<SingleQubitGateLayer>,
                     std::vector<TwoQubitGateLayer>>;
    private:
    struct {
      std::array<qc::fp,3> angles;
      int qubit;
    } struct_U3;

      int N_qubits;

  };

} // namespace na::zoned