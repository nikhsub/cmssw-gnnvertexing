// TrackGenMatcher
//
// Matches a reco::TrackCollection (e.g. the output of PATTrackAndVertexUnpacker,
// "unpackedTracksAndVertices") directly to a genParticle collection, without
// needing to first wrap the tracks into reco::Candidates (which the generic
// PAT MCMatcher template requires).
//
// reco::Track already exposes pt()/eta()/phi()/charge() directly, so matching
// is done by hand: for every (track, genParticle) pair passing the pdgId/status
// selection, deltaR and deltaPt/pt are computed; the best pairs (lowest deltaR)
// are assigned first, optionally resolving ambiguities so no track or gen
// particle is used twice.
//
// Output: edm::ValueMap<reco::GenParticleRef> keyed by the track collection
// (null Ref if unmatched), plus ValueMap<float> deltaR / deltaPtRel for the
// matched pair, for diagnostics/tuning.

#include <algorithm>
#include <cmath>
#include <vector>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"

class TrackGenMatcher : public edm::global::EDProducer<> {
public:
  explicit TrackGenMatcher(const edm::ParameterSet &cfg)
      : tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        genToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))),
        mcPdgId_(cfg.getParameter<std::vector<int>>("mcPdgId")),
        mcStatus_(cfg.getParameter<std::vector<int>>("mcStatus")),
        checkCharge_(cfg.getParameter<bool>("checkCharge")),
        maxDeltaR_(cfg.getParameter<double>("maxDeltaR")),
        maxDPtRel_(cfg.getParameter<double>("maxDPtRel")),
        resolveAmbiguities_(cfg.getParameter<bool>("resolveAmbiguities")) {
    produces<edm::ValueMap<reco::GenParticleRef>>("genMatch");
    produces<edm::ValueMap<int>>("genPartIdx");
    produces<edm::ValueMap<float>>("genMatchDeltaR");
    produces<edm::ValueMap<float>>("genMatchDPtRel");
  }

  ~TrackGenMatcher() override {}

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  const std::vector<int> mcPdgId_;      // abs(pdgId) allowed; empty = any
  const std::vector<int> mcStatus_;     // status allowed; empty = any
  const bool checkCharge_;
  const double maxDeltaR_;              // <=0 disables the cut
  const double maxDPtRel_;              // <=0 disables the cut
  const bool resolveAmbiguities_;       // forbid two tracks matching the same gen particle
};

void TrackGenMatcher::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {
  edm::Handle<reco::TrackCollection> tracks;
  evt.getByToken(tracksToken_, tracks);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genToken_, genParticles);

  const size_t nTrk = tracks->size();
  const size_t nGen = genParticles->size();

  std::vector<reco::GenParticleRef> bestMatch(nTrk);
  std::vector<float> bestDR(nTrk, -1.f);
  std::vector<float> bestDPtRel(nTrk, -1.f);

  struct Pair {
    float dr;
    size_t iTrk;
    size_t jGen;
  };
  std::vector<Pair> pairs;
  pairs.reserve(nTrk);

  for (size_t iTrk = 0; iTrk < nTrk; ++iTrk) {
    const reco::Track &trk = (*tracks)[iTrk];
    for (size_t jGen = 0; jGen < nGen; ++jGen) {
      const reco::GenParticle &gp = (*genParticles)[jGen];

      if (!mcStatus_.empty() &&
          std::find(mcStatus_.begin(), mcStatus_.end(), gp.status()) == mcStatus_.end())
        continue;
      if (!mcPdgId_.empty() &&
          std::find(mcPdgId_.begin(), mcPdgId_.end(), std::abs(gp.pdgId())) == mcPdgId_.end())
        continue;
      if (checkCharge_ && trk.charge() != gp.charge()) continue;
      if (gp.pt() <= 0) continue;

      const float dr = reco::deltaR(trk, gp);
      if (maxDeltaR_ > 0 && dr > maxDeltaR_) continue;

      const float dPtRel = std::abs(trk.pt() - gp.pt()) / gp.pt();
      if (maxDPtRel_ > 0 && dPtRel > maxDPtRel_) continue;

      pairs.push_back({dr, iTrk, jGen});
    }
  }

  // best (lowest deltaR) pairs assigned first
  std::sort(pairs.begin(), pairs.end(), [](const Pair &a, const Pair &b) { return a.dr < b.dr; });

  std::vector<bool> trkUsed(nTrk, false);
  std::vector<bool> genUsed(nGen, false);

  for (const auto &p : pairs) {
    if (trkUsed[p.iTrk]) continue;                       // a track keeps only its best match
    if (resolveAmbiguities_ && genUsed[p.jGen]) continue; // a gen particle used by only one track

    const reco::GenParticle &gp = (*genParticles)[p.jGen];
    bestMatch[p.iTrk] = reco::GenParticleRef(genParticles, p.jGen);
    bestDR[p.iTrk] = p.dr;
    bestDPtRel[p.iTrk] = std::abs((*tracks)[p.iTrk].pt() - gp.pt()) / gp.pt();

    trkUsed[p.iTrk] = true;
    if (resolveAmbiguities_) genUsed[p.jGen] = true;
  }

  auto genMatch = std::make_unique<edm::ValueMap<reco::GenParticleRef>>();
  edm::ValueMap<reco::GenParticleRef>::Filler matchFiller(*genMatch);
  matchFiller.insert(tracks, bestMatch.begin(), bestMatch.end());
  matchFiller.fill();

  std::vector<int> genPartIdx(nTrk, -1);
  for (size_t i = 0; i < nTrk; ++i) {
    if (bestMatch[i].isNonnull()) genPartIdx[i] = bestMatch[i].key();
  }
  auto genPartIdxMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler idxFiller(*genPartIdxMap);
  idxFiller.insert(tracks, genPartIdx.begin(), genPartIdx.end());
  idxFiller.fill();

  auto genMatchDR = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler drFiller(*genMatchDR);
  drFiller.insert(tracks, bestDR.begin(), bestDR.end());
  drFiller.fill();

  auto genMatchDPtRel = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler dptFiller(*genMatchDPtRel);
  dptFiller.insert(tracks, bestDPtRel.begin(), bestDPtRel.end());
  dptFiller.fill();

  evt.put(std::move(genMatch), "genMatch");
  evt.put(std::move(genPartIdxMap), "genPartIdx");
  evt.put(std::move(genMatchDR), "genMatchDeltaR");
  evt.put(std::move(genMatchDPtRel), "genMatchDPtRel");
}

DEFINE_FWK_MODULE(TrackGenMatcher);
