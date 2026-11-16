// Physics-only track-pair (edge) feature producer.
//
// This is the truth-independent half of what used to be TrackGVFeaturesProducer's
// "trackPair" table: everything needed to compute dca/dcaSig/cptopv/pvtoPCA/dotprod/
// pairMom/pairInvMass/pt/eta/deltaR/deltaPhi/deltaEta, with NO dependency on gvTable/
// gvDaughtersTable/genPartIdx (those are generator-truth quantities that only exist
// in simulated MC events and would make this producer unusable on real data).
//
// The truth-dependent half (per-track hadidx/label, per-pair edgeLabel) lives in
// TrackPairLabelProducer, which consumes THIS producer's "trackPair" table (by
// trk1Idx/trk2Idx) to add an aligned edgeLabel column, and only needs to run in the
// training/ntuple-production sequence, never in general reconstruction.
//
// modelInference.cc reads this producer's "trackPair" table directly (by column name)
// instead of recomputing any of this inline, so there is exactly one implementation
// of the edge-feature math, shared between training-ntuple production and inference.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/SimpleVertexTable/interface/TwoTrackVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

class TrackPairFeaturesProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackPairFeaturesProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

  // Cut thresholds -- defaults match exactly what TrackGVFeaturesProducer used at
  // training time. Kept configurable rather than hardcoded so this stays the single
  // source of truth; if these ever change, everything downstream (labeling, training
  // ntuples, inference) picks it up automatically as long as they share one cfi.
  double maxDR_;
  double maxDEta_;
  double maxDPhi_;
  double maxPairEta_;
  double maxInvMass_;
  double maxDca_;
  bool enableTwoTrackVertexFit_;
  bool enableVertexFitTiming_;
  std::unique_ptr<TwoTrackVertexFitter> twoTrackVertexFitter_;

  static constexpr float PION_MASS = 0.13957018f;
};

TrackPairFeaturesProducer::TrackPairFeaturesProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxDR_(iConfig.getParameter<double>("maxDR")),
      maxDEta_(iConfig.getParameter<double>("maxDEta")),
      maxDPhi_(iConfig.getParameter<double>("maxDPhi")),
      maxPairEta_(iConfig.getParameter<double>("maxPairEta")),
      maxInvMass_(iConfig.getParameter<double>("maxInvMass")),
      maxDca_(iConfig.getParameter<double>("maxDca")),
      enableTwoTrackVertexFit_(iConfig.getParameter<bool>("enableTwoTrackVertexFit")),
      enableVertexFitTiming_(iConfig.getParameter<bool>("enableVertexFitTiming")),
      twoTrackVertexFitter_(enableTwoTrackVertexFit_ ? std::make_unique<TwoTrackVertexFitter>() : nullptr) {
  produces<nanoaod::FlatTable>("trackPair");
  if (enableTwoTrackVertexFit_) produces<nanoaod::FlatTable>("vertexFitSummary");
}

void TrackPairFeaturesProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<double>("maxDR", 0.8);
  desc.add<double>("maxDEta", 0.4);
  desc.add<double>("maxDPhi", 0.5);
  desc.add<double>("maxPairEta", 2.5);
  desc.add<double>("maxInvMass", 5.0);
  desc.add<double>("maxDca", 0.4);
  desc.add<bool>("enableTwoTrackVertexFit", false);
  desc.add<bool>("enableVertexFitTiming", false);
  descriptions.addWithDefaultLabel(desc);
}

void TrackPairFeaturesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using Clock = std::chrono::steady_clock;
  const auto totalStart = enableVertexFitTiming_ ? Clock::now() : Clock::time_point{};
  edm::Handle<std::vector<reco::Track>> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  edm::Handle<std::vector<reco::Vertex>> pvH;
  iEvent.getByToken(pvToken_, pvH);
  const reco::Vertex& pv = pvH->front();

  const auto& theB = iSetup.getData(ttbToken_);

  size_t nTrk = tracks.size();
  std::vector<reco::TransientTrack> t_trks(nTrk);
  const auto transientBuildStart = enableVertexFitTiming_ ? Clock::now() : Clock::time_point{};
  for (size_t i = 0; i < nTrk; ++i) t_trks[i] = theB.build(tracks[i]);
  const auto transientBuildEnd = enableVertexFitTiming_ ? Clock::now() : Clock::time_point{};

  std::vector<int> trk1Idx, trk2Idx;
  std::vector<float> outDeltaR, outDeltaPhi, outDeltaEta, outDca, outDcaSig, outCptopv, outPvtoPCA1,
      outPvtoPCA2, outDotprod1, outDotprod2, outPairMom, outPairInvMass, outPt, outEta;

  VertexDistance3D vdist3D;

  std::vector<TwoTrackVertexFitter::Result> fitResults;
  std::vector<reco::TransientTrack> fitterTrackScratch;
  if (enableTwoTrackVertexFit_) {
    fitterTrackScratch.reserve(2);
    fitResults.reserve(nTrk);
  }
  size_t attemptedFits = 0, successfulFits = 0, invalidCovarianceFits = 0, invalidTrackFits = 0, failedFits = 0;
  double fitTimeMs = 0.;
  std::unordered_map<std::uint64_t, TwoTrackVertexFitter::Result> fitCache;
  if (enableTwoTrackVertexFit_) fitCache.reserve(nTrk * 2);

  // Full ordered pairs (i,j) with i != j -- NOT i<j. dotprod1/pvtoPCA1 treat i as
  // "seed" and dotprod2/pvtoPCA2 treat j as "track", so (i,j) and (j,i) are genuinely
  // different feature vectors, not a symmetric edge to deduplicate.
  for (size_t i = 0; i < nTrk; ++i) {
    if (!t_trks[i].isValid()) continue;
    for (size_t j = 0; j < nTrk; ++j) {
      if (i == j) continue;
      if (!t_trks[j].isValid()) continue;

      float dr = reco::deltaR(tracks[i].eta(), tracks[i].phi(), tracks[j].eta(), tracks[j].phi());
      float dEta = std::abs(tracks[i].eta() - tracks[j].eta());
      float dPhi = reco::deltaPhi(tracks[i].phi(), tracks[j].phi());
      if (dr > maxDR_) continue;
      if (dEta > maxDEta_) continue;
      if (std::abs(dPhi) > maxDPhi_) continue;

      float e1 = std::sqrt(tracks[i].p() * tracks[i].p() + PION_MASS * PION_MASS);
      float e2 = std::sqrt(tracks[j].p() * tracks[j].p() + PION_MASS * PION_MASS);
      float sumPx = tracks[i].px() + tracks[j].px();
      float sumPy = tracks[i].py() + tracks[j].py();
      float sumPz = tracks[i].pz() + tracks[j].pz();
      float pt_pair = std::sqrt(sumPx * sumPx + sumPy * sumPy);
      float compositeP = std::sqrt(sumPx * sumPx + sumPy * sumPy + sumPz * sumPz);
      float eta_pair;
      if (pt_pair < 1e-6f) {
        eta_pair = (sumPz > 0) ? 9999.f : -9999.f;
      } else {
        eta_pair = 0.5f * std::log((compositeP + sumPz) / (compositeP - sumPz));
      }
      if (std::abs(eta_pair) > maxPairEta_) continue;

      float sumE = e1 + e2;
      float invMass = std::sqrt(std::max(0.f, sumE * sumE - (sumPx * sumPx + sumPy * sumPy + sumPz * sumPz)));
      if (invMass > maxInvMass_) continue;

      float dca_val = -1.f, dcaSig_val = -1.f, cptopv_val = -1.f;
      float pvToPCA1_val = -1.f, pvToPCA2_val = -1.f, dot1_val = -999.f, dot2_val = -999.f, pairMom_val = -1.f;

      TwoTrackMinimumDistance minDist;
      if (minDist.calculate(t_trks[i].impactPointState(), t_trks[j].impactPointState())) {
        auto m = vdist3D.distance(
            VertexState(minDist.points().second, t_trks[i].impactPointState().cartesianError().position()),
            VertexState(minDist.points().first, t_trks[j].impactPointState().cartesianError().position()));
        dca_val = m.value();
        if (m.error() > 0) dcaSig_val = m.value() / m.error();

        GlobalPoint cp(minDist.crossingPoint());
        GlobalPoint pvp(pv.x(), pv.y(), pv.z());
        GlobalPoint seedPCA = minDist.points().second;
        GlobalPoint trackPCA = minDist.points().first;

        pvToPCA1_val = (seedPCA - pvp).mag();
        pvToPCA2_val = (trackPCA - pvp).mag();
        cptopv_val = (cp - pvp).mag();
        dot2_val = (trackPCA - pvp).unit().dot(t_trks[j].impactPointState().globalDirection().unit());
        dot1_val = (seedPCA - pvp).unit().dot(t_trks[i].impactPointState().globalDirection().unit());

        GlobalVector pairMomentum((Basic3DVector<float>)(t_trks[i].track().momentum() + t_trks[j].track().momentum()));
        pairMom_val = pairMomentum.mag();
      }

      if (dca_val > maxDca_) continue;

      trk1Idx.push_back(static_cast<int>(i));
      trk2Idx.push_back(static_cast<int>(j));
      outDeltaR.push_back(dr);
      outDeltaPhi.push_back(dPhi);
      outDeltaEta.push_back(dEta);
      outDca.push_back(dca_val);
      outDcaSig.push_back(dcaSig_val);
      outCptopv.push_back(cptopv_val);
      outPvtoPCA1.push_back(pvToPCA1_val);
      outPvtoPCA2.push_back(pvToPCA2_val);
      outDotprod1.push_back(dot1_val);
      outDotprod2.push_back(dot2_val);
      outPairMom.push_back(pairMom_val);
      outPairInvMass.push_back(invMass);
      outPt.push_back(pt_pair);
      outEta.push_back(eta_pair);

      if (enableTwoTrackVertexFit_) {
        const auto low = static_cast<std::uint32_t>(std::min(i, j));
        const auto high = static_cast<std::uint32_t>(std::max(i, j));
        const std::uint64_t key = (static_cast<std::uint64_t>(low) << 32) | high;
        const auto cached = fitCache.find(key);
        if (cached != fitCache.end()) {
          // The historical graph stores ordered (i,j) and (j,i) rows because
          // some pair features are endpoint-specific. The vertex fit itself is
          // symmetric, so reuse one physical fit while preserving both rows.
          fitResults.push_back(cached->second);
        } else {
          const auto fitStart = enableVertexFitTiming_ ? Clock::now() : Clock::time_point{};
          const reco::Vertex* pv = (!pvH->empty() && !pvH->front().isFake()) ? &pvH->front() : nullptr;
          auto result = twoTrackVertexFitter_->fit(t_trks[i], t_trks[j], tracks[i], tracks[j], pv, fitterTrackScratch);
          const auto fitEnd = enableVertexFitTiming_ ? Clock::now() : Clock::time_point{};
          if (enableVertexFitTiming_)
            fitTimeMs += std::chrono::duration<double, std::milli>(fitEnd - fitStart).count();
          ++attemptedFits;
          if (result.fitValid) ++successfulFits;
          else if (result.status == 2) ++invalidTrackFits;
          else ++failedFits;
          if (result.fitValid && !result.covarianceValid) ++invalidCovarianceFits;
          fitCache.emplace(key, result);
          fitResults.push_back(std::move(result));
        }
      }
    }
  }

  auto pairTable = std::make_unique<nanoaod::FlatTable>(trk1Idx.size(), "trackPair", false);
  pairTable->addColumn<int>("trk1Idx", trk1Idx, "index of first track in the pair (into the source track collection)");
  pairTable->addColumn<int>("trk2Idx", trk2Idx, "index of second track in the pair");
  pairTable->addColumn<float>("deltaR", outDeltaR, "deltaR between the two tracks");
  pairTable->addColumn<float>("deltaPhi", outDeltaPhi, "deltaPhi between the two tracks");
  pairTable->addColumn<float>("deltaEta", outDeltaEta, "abs(eta1-eta2)");
  pairTable->addColumn<float>("dca", outDca, "distance of closest approach");
  pairTable->addColumn<float>("dcaSig", outDcaSig, "DCA significance");
  pairTable->addColumn<float>("cptopv", outCptopv, "distance from crossing point to PV");
  pairTable->addColumn<float>("pvtoPCA1", outPvtoPCA1, "distance PV to track1 PCA");
  pairTable->addColumn<float>("pvtoPCA2", outPvtoPCA2, "distance PV to track2 PCA");
  pairTable->addColumn<float>("dotprod1", outDotprod1, "direction dot product, track1 side");
  pairTable->addColumn<float>("dotprod2", outDotprod2, "direction dot product, track2 side");
  pairTable->addColumn<float>("pairMom", outPairMom, "magnitude of the pair momentum sum");
  pairTable->addColumn<float>("pairInvMass", outPairInvMass, "invariant mass, pion mass hypothesis");
  pairTable->addColumn<float>("pt", outPt, "transverse momentum of the pair");
  pairTable->addColumn<float>("eta", outEta, "eta of the pair");
  if (enableTwoTrackVertexFit_) {
    std::vector<int> fitValid, fitStatus, covarianceValid, pvValid, lxyValid, signedLxyValid, l3dValid,
                     pointing2DValid, pointing3DValid, perpendicularValid;
    std::vector<float> chi2, ndof, chi2Ndof, probability, vx, vy, vz, covXX, covXY, covXZ,
                       covYY, covYZ, covZZ, sigmaX, sigmaY, sigmaZ, lxy, lxyErr, lxySig,
                       signedLxy, signedLxySig, l3d, l3dErr, l3dSig, cos2D, cos3D, perp, perpSig;
    const size_t n = fitResults.size();
    fitValid.reserve(n); fitStatus.reserve(n); covarianceValid.reserve(n); pvValid.reserve(n);
    lxyValid.reserve(n); signedLxyValid.reserve(n); l3dValid.reserve(n); pointing2DValid.reserve(n);
    pointing3DValid.reserve(n); perpendicularValid.reserve(n); chi2.reserve(n); ndof.reserve(n);
    chi2Ndof.reserve(n); probability.reserve(n); vx.reserve(n); vy.reserve(n); vz.reserve(n);
    covXX.reserve(n); covXY.reserve(n); covXZ.reserve(n); covYY.reserve(n); covYZ.reserve(n); covZZ.reserve(n);
    sigmaX.reserve(n); sigmaY.reserve(n); sigmaZ.reserve(n); lxy.reserve(n); lxyErr.reserve(n); lxySig.reserve(n);
    signedLxy.reserve(n); signedLxySig.reserve(n); l3d.reserve(n); l3dErr.reserve(n); l3dSig.reserve(n);
    cos2D.reserve(n); cos3D.reserve(n); perp.reserve(n); perpSig.reserve(n);
    for (const auto& r : fitResults) {
      fitValid.push_back(r.fitValid); fitStatus.push_back(r.status); covarianceValid.push_back(r.covarianceValid);
      pvValid.push_back(r.pvValid); lxyValid.push_back(r.lxyValid); signedLxyValid.push_back(r.signedLxyValid);
      l3dValid.push_back(r.l3dValid); pointing2DValid.push_back(r.pointing2DValid);
      pointing3DValid.push_back(r.pointing3DValid); perpendicularValid.push_back(r.perpendicularValid);
      chi2.push_back(r.chi2); ndof.push_back(r.ndof); chi2Ndof.push_back(r.chi2Ndof); probability.push_back(r.probability);
      vx.push_back(r.x); vy.push_back(r.y); vz.push_back(r.z);
      covXX.push_back(r.covariance[0]); covXY.push_back(r.covariance[1]); covXZ.push_back(r.covariance[2]);
      covYY.push_back(r.covariance[3]); covYZ.push_back(r.covariance[4]); covZZ.push_back(r.covariance[5]);
      sigmaX.push_back(r.sigmaX); sigmaY.push_back(r.sigmaY); sigmaZ.push_back(r.sigmaZ);
      lxy.push_back(r.lxy); lxyErr.push_back(r.lxyError); lxySig.push_back(r.lxySignificance);
      signedLxy.push_back(r.signedLxy); signedLxySig.push_back(r.signedLxySignificance);
      l3d.push_back(r.l3d); l3dErr.push_back(r.l3dError); l3dSig.push_back(r.l3dSignificance);
      cos2D.push_back(r.cosPointing2D); cos3D.push_back(r.cosPointing3D);
      perp.push_back(r.perpendicular); perpSig.push_back(r.perpendicularSignificance);
    }
    pairTable->addColumn<int>("vtxFitValid", fitValid, "1 for valid Kalman two-track vertex fit");
    pairTable->addColumn<int>("vtxFitStatus", fitStatus, "1 success, 2 invalid track state, 3 invalid fit, 4 exception/numerics");
    pairTable->addColumn<int>("vtxCovarianceValid", covarianceValid, "1 when the fitted position covariance is positive semidefinite");
    pairTable->addColumn<float>("vtxChi2", chi2, "Kalman vertex-fit chi2");
    pairTable->addColumn<float>("vtxNdof", ndof, "Kalman vertex-fit degrees of freedom");
    pairTable->addColumn<float>("vtxChi2Ndof", chi2Ndof, "Kalman vertex-fit chi2/ndof");
    pairTable->addColumn<float>("vtxProb", probability, "chi-square fit probability");
    pairTable->addColumn<float>("vtxX", vx, "fitted vertex x [cm]"); pairTable->addColumn<float>("vtxY", vy, "fitted vertex y [cm]");
    pairTable->addColumn<float>("vtxZ", vz, "fitted vertex z [cm]");
    pairTable->addColumn<float>("vtxCovXX", covXX, "fitted covariance xx [cm^2]"); pairTable->addColumn<float>("vtxCovXY", covXY, "fitted covariance xy [cm^2]");
    pairTable->addColumn<float>("vtxCovXZ", covXZ, "fitted covariance xz [cm^2]"); pairTable->addColumn<float>("vtxCovYY", covYY, "fitted covariance yy [cm^2]");
    pairTable->addColumn<float>("vtxCovYZ", covYZ, "fitted covariance yz [cm^2]"); pairTable->addColumn<float>("vtxCovZZ", covZZ, "fitted covariance zz [cm^2]");
    pairTable->addColumn<float>("vtxSigmaX", sigmaX, "fitted x uncertainty [cm]"); pairTable->addColumn<float>("vtxSigmaY", sigmaY, "fitted y uncertainty [cm]");
    pairTable->addColumn<float>("vtxSigmaZ", sigmaZ, "fitted z uncertainty [cm]");
    pairTable->addColumn<int>("vtxPVValid", pvValid, "1 when a usable leading reconstructed PV exists");
    pairTable->addColumn<int>("vtxLxyValid", lxyValid, "1 when Lxy and its uncertainty are valid");
    pairTable->addColumn<float>("vtxLxy", lxy, "PV-to-fit transverse displacement [cm]"); pairTable->addColumn<float>("vtxLxyErr", lxyErr, "Lxy uncertainty [cm]");
    pairTable->addColumn<float>("vtxLxySig", lxySig, "Lxy significance");
    pairTable->addColumn<int>("vtxSignedLxyValid", signedLxyValid, "1 when signed Lxy and significance are valid");
    pairTable->addColumn<float>("vtxSignedLxy", signedLxy, "signed projected transverse displacement [cm]");
    pairTable->addColumn<float>("vtxSignedLxySig", signedLxySig, "signed Lxy significance");
    pairTable->addColumn<int>("vtxL3DValid", l3dValid, "1 when L3D and its uncertainty are valid");
    pairTable->addColumn<float>("vtxL3D", l3d, "PV-to-fit 3D displacement [cm]"); pairTable->addColumn<float>("vtxL3DErr", l3dErr, "L3D uncertainty [cm]");
    pairTable->addColumn<float>("vtxL3DSig", l3dSig, "L3D significance");
    pairTable->addColumn<int>("vtxPointing2DValid", pointing2DValid, "1 when 2D pointing cosine is valid");
    pairTable->addColumn<float>("vtxCosPointing2D", cos2D, "cosine between transverse displacement and pair momentum");
    pairTable->addColumn<int>("vtxPointing3DValid", pointing3DValid, "1 when 3D pointing cosine is valid");
    pairTable->addColumn<float>("vtxCosPointing3D", cos3D, "cosine between 3D displacement and pair momentum");
    pairTable->addColumn<int>("vtxPerpendicularValid", perpendicularValid, "1 when transverse pointing residual and significance are valid");
    pairTable->addColumn<float>("vtxPerpDisplacement", perp, "transverse displacement perpendicular to pair pT [cm]");
    pairTable->addColumn<float>("vtxPerpDisplacementSig", perpSig, "perpendicular displacement significance");
  }
  iEvent.put(std::move(pairTable), "trackPair");

  if (enableTwoTrackVertexFit_) {
    auto summary = std::make_unique<nanoaod::FlatTable>(1, "vertexFitSummary", true);
    std::vector<int> tracksCount{static_cast<int>(nTrk)}, pairCount{static_cast<int>(trk1Idx.size())};
    std::vector<int> attempted{static_cast<int>(attemptedFits)}, successful{static_cast<int>(successfulFits)};
    std::vector<int> covarianceInvalid{static_cast<int>(invalidCovarianceFits)};
    std::vector<int> invalidTracks{static_cast<int>(invalidTrackFits)}, failed{static_cast<int>(failedFits)};
    std::vector<float> buildMs{enableVertexFitTiming_ ? static_cast<float>(std::chrono::duration<double, std::milli>(transientBuildEnd - transientBuildStart).count()) : -1.f};
    std::vector<float> fitMs{enableVertexFitTiming_ ? static_cast<float>(fitTimeMs) : -1.f};
    std::vector<float> totalMs{enableVertexFitTiming_ ? static_cast<float>(std::chrono::duration<double, std::milli>(Clock::now() - totalStart).count()) : -1.f};
    summary->addColumn<int>("nTracks", tracksCount, "reconstructed tracks in event");
    summary->addColumn<int>("nSelectedPairs", pairCount, "rows in the original selected trackPair table");
    summary->addColumn<int>("nFitAttempts", attempted, "number of attempted two-track fits");
    summary->addColumn<int>("nFitSuccess", successful, "number of valid two-track fits");
    summary->addColumn<int>("nCovarianceInvalid", covarianceInvalid, "valid fits whose position covariance failed the PSD check");
    summary->addColumn<int>("nInvalidTrackState", invalidTracks, "fits skipped due to invalid transient-track state");
    summary->addColumn<int>("nFitFailure", failed, "invalid fit results, exceptions, or invalid fit numerics");
    summary->addColumn<float>("transientTrackBuildMs", buildMs, "transient-track construction time; -1 when timing disabled");
    summary->addColumn<float>("vertexFitMs", fitMs, "cumulative fitter-call time; -1 when timing disabled");
    summary->addColumn<float>("producerTotalMs", totalMs, "producer wall time; -1 when timing disabled");
    iEvent.put(std::move(summary), "vertexFitSummary");
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackPairFeaturesProducer);
