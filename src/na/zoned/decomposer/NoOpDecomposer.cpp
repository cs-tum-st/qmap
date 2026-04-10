/*
 * Copyright (c) 2023 - 2026 Chair for Design Automation, TUM
 * Copyright (c) 2025 - 2026 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "na/zoned/decomposer/NoOpDecomposer.hpp"

#include "ir/QuantumComputation.hpp"
#include "na/zoned/Architecture.hpp"

#include <utility>
#include <vector>

namespace na::zoned {
auto NoOpDecomposer::decompose(
    size_t /* unused */,
    const std::pair<std::vector<SingleQubitGateRefLayer>,
                    std::vector<TwoQubitGateLayer>>& asap_schedule)
    -> std::pair<std::vector<SingleQubitGateLayer>,
                 std::vector<TwoQubitGateLayer>> {
  std::vector<SingleQubitGateLayer> result;
  result.reserve(asap_schedule.first.size());
  for (const auto& layer : asap_schedule.first) {
    SingleQubitGateLayer newLayer;
    newLayer.reserve(layer.size());
    for (const auto& opRef : layer) {
      newLayer.emplace_back(opRef.get().clone());
    }
    result.emplace_back(std::move(newLayer));
  }
  return {std::move(result), asap_schedule.second};
}
} // namespace na::zoned
