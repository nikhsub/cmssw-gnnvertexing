#include <cmath>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/AdaptiveVertexFinder/interface/TTHelpers.h"

// Computes, per track in `tracks` (e.g. unpackedTracksAndVertices) w.r.t. the
// leading vertex in `primaryVertices`, the SAME two quantities that
// TemplatedInclusiveVertexFinder::produce() cuts on when building its `tts`:
//
//   dz       = track.dz(pv.position())                                         (== "maximumLongitudinalImpactParameter" quantity)
//   timeSig  = |track.timeExt() - pv.t()| / sqrt(track.dtErrorExt()^2 + pv.covariance(3,3))
//              only defined when timing info is available, exactly as in IVF   (== "maximumTimeSignificance" quantity)
//
// Also computes whether the track would pass TracksClusteringFromDisplacedSeed's SEED
// selection, i.e. the window cut in clusters():
//
//   ip = IPTools::absoluteImpactParameter3D(tt, pv);
//   ip.first && ip.value() in [seedMin3DIPValue, seedMax3DIPValue]
//            && ip.significance() in [seedMin3DIPSignificance, seedMax3DIPSignificance]
//
// Both the raw values and pass/fail flags (against configurable thresholds,
// defaulting to IVF's own defaults) are produced as edm::ValueMap<...>, so they
// plug into SimpleTrackFlatTableProducer's externalVariables/ExtVar exactly like
// genPartIdx/SVscore already do in your track table.
//
// NOTE: minHits/minPt do NOT need this producer -- they only depend on the track
// itself and can be written directly as `Var(...)` expressions in the flat table
// (see the python snippet below).
class TrackVertexVars : public edm::stream::EDProducer<> {
public:
  explicit TrackVertexVars(const edm::ParameterSet &params);
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &event, const edm::EventSetup &es) override;

  edm::EDGetTokenT<std::vector<reco::Track>> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> token_trackBuilder;

  double maxLIP_;      // "maximumLongitudinalImpactParameter", same meaning/default as IVF
  double maxTimeSig_;  // "maximumTimeSignificance", same meaning/default as IVF

  // seed-selection window, same meaning/defaults as TracksClusteringFromDisplacedSeed's clusterizer PSet
  double seedMin3DIPValue_;
  double seedMax3DIPValue_;
  double seedMin3DIPSignificance_;
  double seedMax3DIPSignificance_;

  static constexpr float SENTINEL = -9999.f;  // used whenever a quantity can't be computed (no PV, no timing, invalid TT)
};

TrackVertexVars::TrackVertexVars(const edm::ParameterSet &params)
    : token_tracks(consumes<std::vector<reco::Track>>(params.getParameter<edm::InputTag>("tracks"))),
      token_primaryVertex(consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"))),
      token_beamSpot(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamSpot"))),
      token_trackBuilder(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxLIP_(params.getParameter<double>("maximumLongitudinalImpactParameter")),
      maxTimeSig_(params.getParameter<double>("maximumTimeSignificance")),
      seedMin3DIPValue_(params.getParameter<double>("seedMin3DIPValue")),
      seedMax3DIPValue_(params.getParameter<double>("seedMax3DIPValue")),
      seedMin3DIPSignificance_(params.getParameter<double>("seedMin3DIPSignificance")),
      seedMax3DIPSignificance_(params.getParameter<double>("seedMax3DIPSignificance")) {
  produces<edm::ValueMap<float>>("dz");
  produces<edm::ValueMap<float>>("timeSig");
  produces<edm::ValueMap<int>>("passLIP");
  produces<edm::ValueMap<int>>("passTimeSig");
  produces<edm::ValueMap<float>>("ip3dValue");
  produces<edm::ValueMap<float>>("ip3dSignificance");
  produces<edm::ValueMap<int>>("passSeed");
}

void TrackVertexVars::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  // same defaults as TemplatedInclusiveVertexFinder::fillDescriptions()
  desc.add<double>("maximumLongitudinalImpactParameter", 0.3);
  desc.add<double>("maximumTimeSignificance", 3.0);
  // same defaults as TracksClusteringFromDisplacedSeed's clusterizer PSet (seedMin3DIP*); override
  // these to match whatever you actually pass to the clusterizer PSet in your IVF config
  desc.add<double>("seedMin3DIPValue", 0.005);
  desc.add<double>("seedMax3DIPValue", 9999.0);
  desc.add<double>("seedMin3DIPSignificance", 1.2);
  desc.add<double>("seedMax3DIPSignificance", 9999.0);
  descriptions.addWithDefaultLabel(desc);
}

void TrackVertexVars::produce(edm::Event &event, const edm::EventSetup &es) {
  using namespace reco;

  edm::Handle<std::vector<Track>> tracks;
  event.getByToken(token_tracks, tracks);

  edm::Handle<VertexCollection> primaryVertices;
  event.getByToken(token_primaryVertex, primaryVertices);

  edm::Handle<BeamSpot> beamSpot;
  event.getByToken(token_beamSpot, beamSpot);

  edm::ESHandle<TransientTrackBuilder> trackBuilder = es.getHandle(token_trackBuilder);

  size_t nTrk = tracks->size();
  std::vector<float> dz(nTrk, SENTINEL);
  std::vector<float> timeSig(nTrk, SENTINEL);
  std::vector<int> passLIP(nTrk, 0);
  std::vector<int> passTimeSig(nTrk, 1);  // IVF does not cut when timing is unavailable -> treat as "pass" by default, same as IVF's `continue` logic (no continue = pass)
  std::vector<float> ip3dValue(nTrk, SENTINEL);
  std::vector<float> ip3dSignificance(nTrk, SENTINEL);
  std::vector<int> passSeed(nTrk, 0);

  if (!primaryVertices->empty()) {
    const Vertex &pv = (*primaryVertices)[0];

    for (size_t i = 0; i < nTrk; ++i) {
      const Track &trk = (*tracks)[i];

      // dz: identical call to what IVF does (there, on tt.track(); here, directly on the reco::Track -- same value)
      float trkDz = trk.dz(pv.position());
      dz[i] = trkDz;
      passLIP[i] = (std::abs(trkDz) <= maxLIP_) ? 1 : 0;

      // needs a TransientTrack for timeExt()/dtErrorExt() and for the IP tools below,
      // built the same way IVF builds tts, then setBeamSpot the same way IVF does
      // before it hands tracks to the clusterizer's seed selection
      TransientTrack tt(tthelpers::buildTT(tracks, trackBuilder, i));
      if (!tt.isValid())
        continue;
      tt.setBeamSpot(*beamSpot);

      // timeSig
      if (edm::isFinite(tt.timeExt()) && pv.covariance(3, 3) > 0.) {
        auto tError = std::sqrt(std::pow(tt.dtErrorExt(), 2) + pv.covariance(3, 3));
        float sig = std::abs(tt.timeExt() - pv.t()) / tError;
        timeSig[i] = sig;
        passTimeSig[i] = (sig <= maxTimeSig_) ? 1 : 0;
      }
      // else: timeSig stays SENTINEL and passTimeSig stays 1 (no timing -> IVF does not reject on this cut)

      // seed selection: same window cut as TracksClusteringFromDisplacedSeed::clusters()
      std::pair<bool, Measurement1D> ip = IPTools::absoluteImpactParameter3D(tt, pv);
      if (ip.first) {
        float ipVal = ip.second.value();
        float ipSig = ip.second.significance();
        ip3dValue[i] = ipVal;
        ip3dSignificance[i] = ipSig;
        passSeed[i] = (ipVal >= seedMin3DIPValue_ && ipVal <= seedMax3DIPValue_ &&
                       ipSig >= seedMin3DIPSignificance_ && ipSig <= seedMax3DIPSignificance_)
                          ? 1
                          : 0;
      }
      // else: ip3dValue/ip3dSignificance stay SENTINEL, passSeed stays 0 (same as IVF: ip.first==false -> not a seed)
    }
  }

  auto dzMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler dzFiller(*dzMap);
  dzFiller.insert(tracks, dz.begin(), dz.end());
  dzFiller.fill();
  event.put(std::move(dzMap), "dz");

  auto timeSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler timeSigFiller(*timeSigMap);
  timeSigFiller.insert(tracks, timeSig.begin(), timeSig.end());
  timeSigFiller.fill();
  event.put(std::move(timeSigMap), "timeSig");

  auto passLIPMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passLIPFiller(*passLIPMap);
  passLIPFiller.insert(tracks, passLIP.begin(), passLIP.end());
  passLIPFiller.fill();
  event.put(std::move(passLIPMap), "passLIP");

  auto passTimeSigMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passTimeSigFiller(*passTimeSigMap);
  passTimeSigFiller.insert(tracks, passTimeSig.begin(), passTimeSig.end());
  passTimeSigFiller.fill();
  event.put(std::move(passTimeSigMap), "passTimeSig");

  auto ip3dValueMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip3dValueFiller(*ip3dValueMap);
  ip3dValueFiller.insert(tracks, ip3dValue.begin(), ip3dValue.end());
  ip3dValueFiller.fill();
  event.put(std::move(ip3dValueMap), "ip3dValue");

  auto ip3dSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip3dSigFiller(*ip3dSigMap);
  ip3dSigFiller.insert(tracks, ip3dSignificance.begin(), ip3dSignificance.end());
  ip3dSigFiller.fill();
  event.put(std::move(ip3dSigMap), "ip3dSignificance");

  auto passSeedMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passSeedFiller(*passSeedMap);
  passSeedFiller.insert(tracks, passSeed.begin(), passSeed.end());
  passSeedFiller.fill();
  event.put(std::move(passSeedMap), "passSeed");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackVertexVars);
