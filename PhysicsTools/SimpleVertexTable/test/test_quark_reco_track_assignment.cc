#include "PhysicsTools/SimpleVertexTable/interface/QuarkRecoTrackAssignment.h"

#include <cassert>
#include <functional>
#include <limits>
#include <vector>

int main() {
  using qfrag::RecoGenCandidate;
  using qfrag::maximumCardinalityMinimumCostAssignment;

  // A greedy choice of track 0 -> gen 0 would strand track 1. The global
  // assignment finds the maximum-cardinality solution instead.
  auto greedyTrap = maximumCardinalityMinimumCostAssignment(
      2, 2, {{0, 0, 1}, {0, 1, 2}, {1, 0, 1}});
  assert((greedyTrap == std::vector<int>{1, 0}));

  // Among maximum-cardinality assignments, choose the lower total cost.
  auto minimumCost = maximumCardinalityMinimumCostAssignment(
      2, 2, {{0, 0, 1}, {0, 1, 8}, {1, 0, 8}, {1, 1, 1}});
  assert((minimumCost == std::vector<int>{0, 1}));

  // Stable indices resolve exact ties reproducibly.
  auto tied = maximumCardinalityMinimumCostAssignment(
      1, 2, {{0, 1, 3}, {0, 0, 3}});
  assert((tied == std::vector<int>{0}));

  // A shared generator can be used once only; the lower-index track wins a tie.
  auto duplicateTracks = maximumCardinalityMinimumCostAssignment(
      2, 1, {{1, 0, 2}, {0, 0, 2}});
  assert((duplicateTracks == std::vector<int>{0, -1}));

  auto noCandidates = maximumCardinalityMinimumCostAssignment(2, 2, {});
  assert((noCandidates == std::vector<int>{-1, -1}));

  // Exhaustively compare all 3x3 candidate graphs against enumeration. This
  // catches cardinality or min-cost regressions beyond the hand-picked trap.
  for (unsigned mask = 0; mask < (1u << 9); ++mask) {
    std::vector<RecoGenCandidate> graph;
    int costs[3][3] = {};
    for (int track = 0; track < 3; ++track) {
      for (int gen = 0; gen < 3; ++gen) {
        const unsigned bit = 1u << (track * 3 + gen);
        costs[track][gen] = 1 + static_cast<int>((mask + track * 3 + gen * 5) % 7);
        if (mask & bit) graph.push_back({track, gen, costs[track][gen]});
      }
    }
    int bestCardinality = -1;
    int bestCost = std::numeric_limits<int>::max();
    std::function<void(int, int, int, unsigned)> enumerate = [&](int track, int cardinality, int cost, unsigned used) {
      if (track == 3) {
        if (cardinality > bestCardinality || (cardinality == bestCardinality && cost < bestCost)) {
          bestCardinality = cardinality;
          bestCost = cost;
        }
        return;
      }
      enumerate(track + 1, cardinality, cost, used);
      for (int gen = 0; gen < 3; ++gen) {
        const unsigned bit = 1u << (track * 3 + gen);
        const unsigned genBit = 1u << gen;
        if ((mask & bit) && !(used & genBit))
          enumerate(track + 1, cardinality + 1, cost + costs[track][gen], used | genBit);
      }
    };
    enumerate(0, 0, 0, 0);
    const auto assignment = maximumCardinalityMinimumCostAssignment(3, 3, graph);
    int actualCardinality = 0, actualCost = 0;
    for (int track = 0; track < 3; ++track) {
      const int gen = assignment[track];
      if (gen >= 0) {
        ++actualCardinality;
        actualCost += costs[track][gen];
      }
    }
    assert(actualCardinality == bestCardinality);
    assert(actualCost == bestCost);
  }
  return 0;
}
