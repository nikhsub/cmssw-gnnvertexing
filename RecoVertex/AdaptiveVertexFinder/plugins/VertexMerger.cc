#include <algorithm>
#include <memory>
#include <set>
#include <type_traits>

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

// ---------------------------------------------------------------------
// Trait used to pull a uniform list of track keys (TrackRef or
// CandidatePtr) out of either flavor of secondary vertex. This is what
// lets the same templated merger code look up entries in the edge-score
// ValueMaps regardless of which VTX/track collection flavor is in use,
// exactly mirroring how TemplatedInclusiveVertexFinder branches per VTX.
// ---------------------------------------------------------------------
template <class VTX>
struct VertexTrackKeys;

template <>
struct VertexTrackKeys<reco::Vertex> {
  using KeyType = reco::TrackRef;
  static std::vector<KeyType> get(const reco::Vertex &v) {
    std::vector<KeyType> out;
    out.reserve(v.tracksSize());
    for (auto it = v.tracks_begin(); it != v.tracks_end(); ++it) {
      out.push_back(it->castTo<reco::TrackRef>());
    }
    return out;
  }
};

template <>
struct VertexTrackKeys<reco::VertexCompositePtrCandidate> {
  using KeyType = reco::CandidatePtr;
  static std::vector<KeyType> get(const reco::VertexCompositePtrCandidate &v) { return v.daughterPtrVector(); }
};

template <class VTX>
class TemplatedVertexMerger : public edm::stream::EDProducer<> {
public:
  typedef std::vector<VTX> Product;
  typedef typename VertexTrackKeys<VTX>::KeyType TrackKey;

  TemplatedVertexMerger(const edm::ParameterSet &params);

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void produce(edm::Event &event, const edm::EventSetup &es) override;

private:
  // Edge-score-weighted analogue of vertexTools::computeSharedTracks(a, b):
  // for every track in `a`, contributes 1.0 if that same track is also in
  // `b`, otherwise contributes the best edge score connecting it to any
  // track that IS in `b` (0 if there is no such edge). Returns the sum
  // normalized by the number of tracks in `a`, so it is directly
  // comparable to the old fractional shared-tracks metric and to
  // maxFraction.
  double weightedSharedScore(const VTX &a,
                              const VTX &b,
                              const edm::ValueMap<std::vector<int>> &edgeIndices,
                              const edm::ValueMap<std::vector<float>> &edgeScores) const;

  edm::EDGetTokenT<Product> token_secondaryVertex;
  edm::EDGetTokenT<edm::ValueMap<std::vector<float>>> token_edgeScores;
  edm::EDGetTokenT<edm::ValueMap<std::vector<int>>> token_edgeIndices;
  double maxFraction;
  double minSignificance;
  bool useEdgeScore;
};

template <class VTX>
TemplatedVertexMerger<VTX>::TemplatedVertexMerger(const edm::ParameterSet &params)
    : maxFraction(params.getParameter<double>("maxFraction")),
      minSignificance(params.getParameter<double>("minSignificance")),
      useEdgeScore(params.getParameter<bool>("useEdgeScore")) {
  token_secondaryVertex = consumes<Product>(params.getParameter<edm::InputTag>("secondaryVertices"));
  token_edgeScores = consumes<edm::ValueMap<std::vector<float>>>(params.getParameter<edm::InputTag>("edgeScores"));
  token_edgeIndices = consumes<edm::ValueMap<std::vector<int>>>(params.getParameter<edm::InputTag>("edgeIndices"));
  produces<Product>();
}

template <class VTX>
double TemplatedVertexMerger<VTX>::weightedSharedScore(const VTX &a,
                                                         const VTX &b,
                                                         const edm::ValueMap<std::vector<int>> &edgeIndices,
                                                         const edm::ValueMap<std::vector<float>> &edgeScores) const {
  std::vector<TrackKey> tracksA = VertexTrackKeys<VTX>::get(a);
  std::vector<TrackKey> tracksB = VertexTrackKeys<VTX>::get(b);
  if (tracksA.empty())
    return 0.;

  std::set<size_t> keysB;
  for (auto const &tb : tracksB)
    keysB.insert(tb.key());

  double sum = 0.;
  for (auto const &ta : tracksA) {
    // Literally shared track: full weight, same as the old behavior.
    if (keysB.count(ta.key())) {
      sum += 1.0;
      continue;
    }
    // Otherwise fall back to the best edge score linking this track to
    // any track that is actually part of vertex b.
    if (!edgeIndices.contains(ta.id()) || !edgeScores.contains(ta.id()))
      continue;
    const std::vector<int> &idxs = edgeIndices[ta];
    const std::vector<float> &scores = edgeScores[ta];
    float best = 0.f;
    for (size_t k = 0; k < idxs.size() && k < scores.size(); ++k) {
      if (idxs[k] >= 0 && keysB.count(static_cast<size_t>(idxs[k]))) {
        best = std::max(best, scores[k]);
      }
    }
    sum += best;
  }
  return sum / tracksA.size();
}

template <class VTX>
void TemplatedVertexMerger<VTX>::produce(edm::Event &event, const edm::EventSetup &es) {
  using namespace reco;

  edm::Handle<Product> secondaryVertices;
  event.getByToken(token_secondaryVertex, secondaryVertices);

  edm::Handle<edm::ValueMap<std::vector<float>>> edgeScores;
  event.getByToken(token_edgeScores, edgeScores);
  edm::Handle<edm::ValueMap<std::vector<int>>> edgeIndices;
  event.getByToken(token_edgeIndices, edgeIndices);

  const bool haveEdgeInfo = useEdgeScore && edgeScores.isValid() && edgeIndices.isValid();

  VertexDistance3D dist;
  auto recoVertices = std::make_unique<Product>();
  for (typename Product::const_iterator sv = secondaryVertices->begin(); sv != secondaryVertices->end(); ++sv) {
    recoVertices->push_back(*sv);
  }
  for (typename Product::iterator sv = recoVertices->begin(); sv != recoVertices->end(); ++sv) {
    bool shared = false;
    VertexState s1(RecoVertex::convertPos(sv->position()), RecoVertex::convertError(sv->error()));
    for (typename Product::iterator sv2 = recoVertices->begin(); sv2 != recoVertices->end(); ++sv2) {
      if (sv - sv2 == 0)
        continue;
      VertexState s2(RecoVertex::convertPos(sv2->position()), RecoVertex::convertError(sv2->error()));

      double fr, frRev;
      if (haveEdgeInfo) {
        fr = weightedSharedScore(*sv2, *sv, *edgeIndices, *edgeScores);
        frRev = weightedSharedScore(*sv, *sv2, *edgeIndices, *edgeScores);
      } else {
        fr = vertexTools::computeSharedTracks(*sv2, *sv);
        frRev = vertexTools::computeSharedTracks(*sv, *sv2);
      }

      if (fr > maxFraction && dist.distance(s1, s2).significance() < minSignificance && fr >= frRev) {
        shared = true;
      }
    }
    if (shared) {
      sv = recoVertices->erase(sv) - 1;
    }
  }

  event.put(std::move(recoVertices));
}

template <class VTX>
void TemplatedVertexMerger<VTX>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("maxFraction", 0.7);
  desc.add<double>("minSignificance", 2);
  desc.add<edm::InputTag>("secondaryVertices", edm::InputTag("inclusiveVertexFinder"));
  desc.add<edm::InputTag>("edgeScores", edm::InputTag("dummyTrackValueMap", "edgeScores"));
  desc.add<edm::InputTag>("edgeIndices", edm::InputTag("dummyTrackValueMap", "edgeIndices"));
  desc.add<bool>("useEdgeScore", true);
  descriptions.addWithDefaultLabel(desc);
}

typedef TemplatedVertexMerger<reco::Vertex> VertexMerger;
typedef TemplatedVertexMerger<reco::VertexCompositePtrCandidate> CandidateVertexMerger;

DEFINE_FWK_MODULE(VertexMerger);
DEFINE_FWK_MODULE(CandidateVertexMerger);