//
// Created by cpsch on 16.12.2025.
//

#include "ir/QuantumComputation.hpp"
#include "na/zoned/decomposer/Decomposer.hpp"
#include "na/zoned/scheduler/ASAPScheduler.hpp"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

namespace na::zoned {
constexpr std::string_view architectureJson = R"({
  "name": "asap_scheduler_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";


class DecomposerTest : public ::testing::Test {
protected:
  Decomposer decomposer=Decomposer(4);
  Architecture architecture;
  ASAPScheduler::Config config{.maxFillingFactor = .8};
  ASAPScheduler scheduler;
  DecomposerTest()
      : architecture(Architecture::fromJSONString(architectureJson)),
        scheduler(architecture, config) {}
};

qc::fp epsilon=1e-5;

TEST(Test,ThreeQuaternionCombiTest) {
  std::array<qc::fp,4> q1={cos(qc::PI_4),0,0,sin(qc::PI_4)};
  std::array<qc::fp,4> q2={cos(qc::PI_2),0,sin(qc::PI_2),0};
  std::array<qc::fp,4> q12=Decomposer::combine_quaternions(q1,q2);
  EXPECT_THAT(q12,::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
    ::testing::DoubleNear(-1*cos(qc::PI_4),epsilon),::testing::DoubleNear(cos(qc::PI_4),epsilon),::testing::DoubleNear(0,epsilon)));
  std::array<qc::fp,4> q3={cos(qc::PI_2),0,0,sin(qc::PI_2)};
  std::array<qc::fp,4> q13=Decomposer::combine_quaternions(q12,q3);
  EXPECT_THAT(q13,::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
    ::testing::DoubleNear(cos(qc::PI_4),epsilon),::testing::DoubleNear(cos(qc::PI_4),epsilon),::testing::DoubleNear(0,epsilon)));

}

TEST(Test,ThreeQuaternionU3Test) {
  std::array<qc::fp,4> q1={cos(qc::PI_2),0,0,sin(qc::PI_2)};
  std::array<qc::fp,4> q2={cos(qc::PI_4/2),0,sin(qc::PI_4/2),0};
  std::array<qc::fp,4> q12=Decomposer::combine_quaternions(q1,q2);
  EXPECT_THAT(q12,::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
    ::testing::DoubleNear(-1*sin(qc::PI_4/2),epsilon),::testing::DoubleNear(0,epsilon),::testing::DoubleNear(cos(qc::PI_4/2),epsilon)));
  std::array<qc::fp,4> q3={cos(qc::PI_4),0,0,sin(qc::PI_4)};
  std::array<qc::fp,4> q13=Decomposer::combine_quaternions(q12,q3);
  qc::fp r2=1/sqrt(2);
  EXPECT_THAT(q13,::testing::ElementsAre(::testing::DoubleNear(-r2*cos(qc::PI_4/2),epsilon),
   ::testing::DoubleNear(-r2*sin(qc::PI_4/2),epsilon),::testing::DoubleNear(r2*sin(qc::PI_4/2),epsilon),::testing::DoubleNear(r2*cos(qc::PI_4/2),epsilon)));

}

TEST(Test,SingleXGateAngleTest) {
  size_t n=1;
  const qc::Operation *op=new qc::StandardOperation(0,qc::X);

  qc::QuantumComputation qc(n);
  qc.rx(qc::PI, 0);
  std::array<qc::fp,4> q=Decomposer::convert_gate_to_quaternion(
      std::reference_wrapper<const qc::Operation>(*op));
  EXPECT_THAT(Decomposer::get_U3_angles_from_quaternion(q),::testing::ElementsAre(::testing::DoubleNear(qc::PI,epsilon),
    ::testing::DoubleNear(-1*qc::PI_2,epsilon),::testing::DoubleNear(qc::PI_2,epsilon)));
}

TEST(Test,SingleU3GateAngleTest) {
  size_t n=1;
  const qc::Operation *op=new qc::StandardOperation(0,qc::U,{qc::PI_4,qc::PI,qc::PI_2});
  std::array<qc::fp,4> q=Decomposer::convert_gate_to_quaternion(
      std::reference_wrapper<const qc::Operation>(*op));
  qc::fp r2=1/sqrt(2);
  EXPECT_THAT(q,::testing::ElementsAre(::testing::DoubleNear(-r2*cos(qc::PI_4/2),epsilon),
   ::testing::DoubleNear(-r2*sin(qc::PI_4/2),epsilon),::testing::DoubleNear(r2*sin(qc::PI_4/2),epsilon),::testing::DoubleNear(r2*cos(qc::PI_4/2),epsilon)));

  EXPECT_THAT(Decomposer::get_U3_angles_from_quaternion(q),::testing::ElementsAre(::testing::DoubleNear(qc::PI_4,epsilon),
    ::testing::DoubleNear(qc::PI,epsilon),::testing::DoubleNear(qc::PI_2,epsilon)));
}

TEST(Test, ThetaPiAngleTest) {
  size_t n=1;
  const qc::Operation *op=new qc::StandardOperation(0,qc::U,{qc::PI,qc::PI,qc::PI_2});
  std::array<qc::fp,4> q=Decomposer::convert_gate_to_quaternion(
      std::reference_wrapper<const qc::Operation>(*op));
  qc::fp r2=1/sqrt(2);
  EXPECT_THAT(q,::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
   ::testing::DoubleNear(-r2,epsilon),::testing::DoubleNear(r2,epsilon),::testing::DoubleNear(0,epsilon)));
  EXPECT_THAT(Decomposer::get_U3_angles_from_quaternion(q),::testing::ElementsAre(::testing::DoubleNear(qc::PI,epsilon),
    ::testing::DoubleNear(0,epsilon),::testing::DoubleNear(-1*qc::PI_2,epsilon)));

}

TEST(Test, ThetaZeroAngleTest) {
  size_t n=1;
  const qc::Operation *op=new qc::StandardOperation(0,qc::U,{0,qc::PI,qc::PI_2});
  std::array<qc::fp,4> q=Decomposer::convert_gate_to_quaternion(
      std::reference_wrapper<const qc::Operation>(*op));
  qc::fp r2=1/sqrt(2);
  EXPECT_THAT(q,::testing::ElementsAre(::testing::DoubleNear(-r2,epsilon),
   ::testing::DoubleNear(0,epsilon),::testing::DoubleNear(0,epsilon),::testing::DoubleNear(r2,epsilon)));

  EXPECT_THAT(Decomposer::get_U3_angles_from_quaternion(q),::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
    ::testing::DoubleNear(0,epsilon),::testing::DoubleNear(3*qc::PI_2,epsilon)));

}

//One or Two Decomposition Tests?


TEST(Test, RXDecompositionTest) {
  size_t n=1;
  std::array<qc::fp,3> rx={qc::PI,-qc::PI_2,qc::PI_2};
  //middle entry gamma_minus is -3/2 PI??
  EXPECT_THAT(Decomposer::get_decomposition_angles(rx,qc::PI),::testing::ElementsAre(::testing::DoubleNear(qc::PI,epsilon),
    ::testing::DoubleNear(qc::PI_2,epsilon),::testing::DoubleNear(qc::PI_2,epsilon)));

}

TEST(Test,U3DecompositionTest) {
  size_t n=1;
  std::array<qc::fp,3> u3={qc::PI_4,qc::PI,qc::PI_2};
  // gamm minus i - PI_" not PI_2 and gam_plus is 0 not PI!!
  EXPECT_THAT(Decomposer::get_decomposition_angles(u3,qc::PI_4),::testing::ElementsAre(::testing::DoubleNear(qc::PI,epsilon),
    ::testing::DoubleNear(qc::PI,epsilon),::testing::DoubleNear(-qc::PI_2,epsilon)));
}

TEST(Test,DoubleDecompositionTest) {
  size_t n=1;
  std::array<qc::fp,3> x1={qc::PI,-qc::PI_2,qc::PI_2};
  std::array<qc::fp,3> z2={0,0,qc::PI};
  //gamm_min is - 3/2 PI not 0 and gamm_plus is PI_2 not zero
  EXPECT_THAT(Decomposer::get_decomposition_angles(x1,qc::PI),::testing::ElementsAre(::testing::DoubleNear(qc::PI,epsilon),
    ::testing::DoubleNear(0,epsilon),::testing::DoubleNear(0,epsilon)));
  // gamm_min is 0 not PI_2 and gamm_plus is PI not PI_2
  EXPECT_THAT(Decomposer::get_decomposition_angles(z2,qc::PI),::testing::ElementsAre(::testing::DoubleNear(0,epsilon),
   ::testing::DoubleNear(qc::PI_2,epsilon),::testing::DoubleNear(qc::PI_2,epsilon)));

}

TEST_F(DecomposerTest, SingleRXGate) {
  //    ┌───────┐
  // q: ┤ Rx(π) ├
  //    └───────┘
  int n=1;
  qc::QuantumComputation qc(n);
  qc.rx(qc::PI, 0);
  Decomposer decomposer=Decomposer(n);
  const auto& sched =
      scheduler.schedule(qc);
  // auto decomp=decomposer.decompose(sched);
}

TEST_F(DecomposerTest, SingleU3Gate) {
  //    ┌───────────────┐
  // q:─┤ U3(π/4,π,π/2) ├─
  //    └───────────────┘
  size_t n=1;
  qc::QuantumComputation qc(n);
  qc.u(qc::PI_4,qc::PI,qc::PI_2,0);
  Decomposer decomposer=Decomposer(n);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);

}

TEST_F(DecomposerTest, TwoPauliGatesOneQubit) {
  //    ┌───────┐  ┌───────┐
  // q: ┤   X   ├──┤   Z   ├
  //    └───────┘  └───────┘
  size_t n=1;
  qc::QuantumComputation qc(n);
  qc.x(0);
  qc.z(0);
  Decomposer decomposer=Decomposer(n);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);

}

TEST_F(DecomposerTest, TwoPauliGatesTwoQubits) {
  //       ┌───────┐
  // q_0: ─┤   X   ├─
  //       └───────┘
  //       ┌───────┐
  // q_1: ─┤   Z   ├─
  //       └───────┘

  size_t n=2;
  qc::QuantumComputation qc(n);
  qc.x(0);
  qc.z(1);
  Decomposer decomposer=Decomposer(n);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);

}

TEST_F(DecomposerTest, TwoU3GatesTwoQubits) {
  //       ┌───────────────┐
  // q_0: ─┤ U3(π/4,π,π/2) ├─
  //       └───────────────┘
  //       ┌───────────────┐
  // q_1: ─┤ U3(π/2,π/2,π) ├─
  //       └───────────────┘
  size_t n=2;
  qc::QuantumComputation qc(n);
  qc.u(qc::PI_4,qc::PI,qc::PI_2,0);
  qc.u(qc::PI_2,qc::PI_2,qc::PI,1);
  Decomposer decomposer=Decomposer(n);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);

}


} // namespace na::zoned