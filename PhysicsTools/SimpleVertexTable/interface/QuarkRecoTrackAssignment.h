#ifndef PhysicsTools_SimpleVertexTable_QuarkRecoTrackAssignment_h
#define PhysicsTools_SimpleVertexTable_QuarkRecoTrackAssignment_h

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

namespace qfrag {
  struct RecoGenCandidate {
    int track = -1;
    int gen = -1;
    std::int64_t cost = 0;
  };

  // Sparse successive-shortest-augmenting-path min-cost max-flow. Sending
  // flow until no augmenting path remains first maximizes cardinality; reduced
  // costs then minimize total candidate cost at that cardinality.
  inline std::vector<int> maximumCardinalityMinimumCostAssignment(
      int nTracks, int nGen, std::vector<RecoGenCandidate> candidates) {
    struct Edge {
      int to;
      int reverse;
      int capacity;
      std::int64_t cost;
    };
    const int source = 0;
    const int firstTrack = 1;
    const int firstGen = firstTrack + nTracks;
    const int sink = firstGen + nGen;
    const int nNodes = sink + 1;
    std::vector<std::vector<Edge>> graph(nNodes);
    auto addEdge = [&](int from, int to, int capacity, std::int64_t cost) {
      const int forward = graph[from].size();
      const int reverse = graph[to].size();
      graph[from].push_back({to, reverse, capacity, cost});
      graph[to].push_back({from, forward, 0, -cost});
    };

    for (int i = 0; i < nTracks; ++i) addEdge(source, firstTrack + i, 1, 0);
    for (int j = 0; j < nGen; ++j) addEdge(firstGen + j, sink, 1, 0);
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
      return std::tie(a.track, a.gen, a.cost) < std::tie(b.track, b.gen, b.cost);
    });
    for (const auto& candidate : candidates) {
      if (candidate.track >= 0 && candidate.track < nTracks && candidate.gen >= 0 && candidate.gen < nGen)
        addEdge(firstTrack + candidate.track, firstGen + candidate.gen, 1, candidate.cost);
    }

    constexpr auto infinity = std::numeric_limits<std::int64_t>::max() / 4;
    std::vector<std::int64_t> potential(nNodes, 0);
    while (true) {
      std::vector<std::int64_t> distance(nNodes, infinity);
      std::vector<std::pair<int, int>> previous(nNodes, {-1, -1});
      std::vector<bool> settled(nNodes, false);
      using QueueEntry = std::pair<std::int64_t, int>;
      std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
      distance[source] = 0;
      queue.push({0, source});
      while (!queue.empty()) {
        const auto [dist, node] = queue.top();
        queue.pop();
        if (settled[node] || dist != distance[node]) continue;
        settled[node] = true;
        for (int edgeIndex = 0; edgeIndex < static_cast<int>(graph[node].size()); ++edgeIndex) {
          const auto& edge = graph[node][edgeIndex];
          if (edge.capacity == 0 || settled[edge.to]) continue;
          const auto reduced = edge.cost + potential[node] - potential[edge.to];
          const auto next = dist + reduced;
          const auto predecessor = std::make_pair(node, edgeIndex);
          if (next < distance[edge.to] ||
              (next == distance[edge.to] &&
               (previous[edge.to].first < 0 || predecessor < previous[edge.to]))) {
            distance[edge.to] = next;
            previous[edge.to] = predecessor;
            queue.push({next, edge.to});
          }
        }
      }
      if (distance[sink] == infinity) break;
      for (int node = 0; node < nNodes; ++node)
        if (distance[node] != infinity) potential[node] += distance[node];
      for (int node = sink; node != source;) {
        const auto [parent, edgeIndex] = previous[node];
        auto& edge = graph[parent][edgeIndex];
        --edge.capacity;
        ++graph[node][edge.reverse].capacity;
        node = parent;
      }
    }

    std::vector<int> assignment(nTracks, -1);
    for (int track = 0; track < nTracks; ++track) {
      for (const auto& edge : graph[firstTrack + track]) {
        if (edge.to >= firstGen && edge.to < sink && edge.capacity == 0) {
          assignment[track] = edge.to - firstGen;
          break;
        }
      }
    }
    return assignment;
  }
}  // namespace qfrag

#endif
