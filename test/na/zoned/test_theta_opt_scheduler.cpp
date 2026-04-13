#include "ir/QuantumComputation.hpp"
#include "na/zoned/decomposer/NativeGateDecomposer.hpp"
#include "na/zoned/scheduler/ASAPScheduler.hpp"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

//
// Created by cpsch on 08.04.2026.
//
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

class ThetaOptTest : public ::testing::Test {
protected:
  Architecture architecture;
  ASAPScheduler::Config schedulerConfig{.maxFillingFactor = .8};
  ASAPScheduler scheduler;
  NativeGateDecomposer::Config decomposerConfig{.theta_opt_schedule = true,
                                                .check_final_cond = false};
  NativeGateDecomposer decomposer;
  ThetaOptTest()
      : architecture(Architecture::fromJSONString(architectureJson)),
        scheduler(architecture, schedulerConfig),
        decomposer(architecture, decomposerConfig) {}
};

constexpr static qc::fp epsilon = std::numeric_limits<qc::fp>::epsilon() * 1024;

TEST_F(ThetaOptTest, GraphTest) {
  // Circuit
  //        ┌───────┐       ┌───────┐
  //  q_0: ─┤   X   ├───■───┤   Z   ├─
  //        └───────┘   │   └───────┘
  //                    │   ┌───────┐
  //  q_1: ─────────────■───┤   X   ├─
  //                        └───────┘

  size_t n = 2;
  qc::QuantumComputation qc(n);
  qc.x(0);
  qc.cz(0, 1);
  qc.z(0);
  qc.x(1);

  auto schedule = scheduler.schedule(qc);
  auto one_qubit_gates = NativeGateDecomposer::transformToU3(schedule.first, n);
  auto graph = NativeGateDecomposer::convert_circ_to_dag(
      {one_qubit_gates, schedule.second}, n);
  EXPECT_EQ(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(0)).qubit,
      0);
  EXPECT_THAT(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(0)).angles,
      ::testing::ElementsAre(
          ::testing::DoubleNear(one_qubit_gates.front().front().angles.at(0),
                                epsilon),
          ::testing::DoubleNear(one_qubit_gates.front().front().angles.at(1),
                                epsilon),
          ::testing::DoubleNear(one_qubit_gates.front().front().angles.at(2),
                                epsilon)));
  EXPECT_THAT(graph.get_adjacent(0),
              ::testing::ElementsAre(std::pair<std::size_t, double>(1, 1.0)));

  auto tqg = std::get<std::array<qc::Qubit, 2>>(graph.get_Node_Value(1));
  EXPECT_THAT(tqg, ::testing::ElementsAre(static_cast<qc::Qubit>(0),
                                          static_cast<qc::Qubit>(1)));
  EXPECT_THAT(graph.get_adjacent(1),
              ::testing::ElementsAre(std::pair<std::size_t, double>(2, 1.0),
                                     std::pair<std::size_t, double>(3, 1.0)));

  EXPECT_EQ(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(2)).qubit,
      0);
  EXPECT_THAT(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(2)).angles,
      ::testing::ElementsAre(
          ::testing::DoubleNear(one_qubit_gates.at(1).front().angles.at(0),
                                epsilon),
          ::testing::DoubleNear(one_qubit_gates.at(1).front().angles.at(1), epsilon),
          ::testing::DoubleNear(one_qubit_gates.at(1).front().angles.at(2),
                                epsilon)));
  EXPECT_THAT(graph.get_adjacent(2), ::testing::IsEmpty());

  EXPECT_EQ(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(3)).qubit,
      1);
  EXPECT_THAT(
      std::get<NativeGateDecomposer::StructU3>(graph.get_Node_Value(3)).angles,
      ::testing::ElementsAre(
          ::testing::DoubleNear(one_qubit_gates.at(1).at(1).angles.at(0), epsilon),
          ::testing::DoubleNear(one_qubit_gates.at(1).at(1).angles.at(1),
                                epsilon),
          ::testing::DoubleNear(one_qubit_gates.at(1).at(1).angles.at(2), epsilon)));
  EXPECT_THAT(graph.get_adjacent(3), ::testing::IsEmpty());
}

TEST_F(ThetaOptTest, SiftTest) {
  // Circuit
  //         ┌───────┐           ┌───────┐       ┌───────┐
  //  q_0: ──┤   X   ├───■───────┤   Z   ├───■───┤   Y   ├─
  //         └───────┘   │       └───────┘   │   └───────┘
  //                     │       ┌───────┐   │
  //  q_1: ──────────────■───■───┤   X   ├───│─────────────
  //                         │   └───────┘   │
  //         ┌───────┐       │   ┌───────┐   │
  //  q_1: ──┤   X   ├───────■───┤   Y   ├───■─────────────
  //         └───────┘           └───────┘

  size_t n = 3;
  qc::QuantumComputation qc(n);
  qc.x(0);
  qc.x(2);
  qc.cz(0, 1);
  qc.cz(1, 2);
  qc.z(0);
  qc.x(1);
  qc.y(2);
  qc.cz(0, 2);
  qc.y(0);

  auto schedule = scheduler.schedule(qc);
  auto one_qubit_gates = NativeGateDecomposer::transformToU3(schedule.first, n);
  auto graph = NativeGateDecomposer::convert_circ_to_dag(
      {one_qubit_gates, schedule.second}, n);

  // Perform Sift (multiple times???)
  std::vector<std::size_t> v_start = {0};
  for (auto i = 0; i < graph.size(); ++i) {
    v_start.push_back(i);
  }
  std::array<std::vector<size_t>, 3> v =
      NativeGateDecomposer::sift(graph, v_start, n);

  EXPECT_THAT(v[0], ::testing::IsEmpty());
  EXPECT_THAT(v[1], ::testing::ElementsAre(0, 1));
  EXPECT_THAT(v[2], ::testing::ElementsAre(2, 3, 4, 5, 6, 7, 8));

  v = NativeGateDecomposer::sift(graph, v[2], n);

  EXPECT_THAT(v[0], ::testing::ElementsAre(2, 4));
  EXPECT_THAT(v[1], ::testing::ElementsAre(3, 5, 6));
  EXPECT_THAT(v[2], ::testing::ElementsAre(7, 8));

  v = NativeGateDecomposer::sift(graph, v[2], n);

  EXPECT_THAT(v[0], ::testing::ElementsAre(7));
  EXPECT_THAT(v[1], ::testing::ElementsAre(8));
  EXPECT_THAT(v[2], ::testing::IsEmpty());
}

TEST_F(ThetaOptTest, RecursionMemoTest) {
  // Circuit
}

TEST_F(ThetaOptTest, RecursionBaseTest) {
  // Circuit
}

TEST_F(ThetaOptTest, RecursionCallTest) {
  // Circuit
}

TEST_F(ThetaOptTest, ShortesPathTest) {
  // Subproblem graph!

}

TEST_F(ThetaOptTest, BuildScheduleTest) {
  // Circuit
}

TEST_F(ThetaOptTest, CompleteTestSmall) {
  // Circuit
}

TEST_F(ThetaOptTest, CompleteTestBig) {

}
} // namespace na::zoned