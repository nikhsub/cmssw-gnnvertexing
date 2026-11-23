#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "PhysicsTools/SimpleVertexTable/interface/TwoTrackVertexFitter.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
  struct Candidate {
    int first = -1;
    int second = -1;
    int edge12 = -1;
    int edge21 = -1;
    float logit12 = 0.f;
    float logit21 = 0.f;
    float score12 = 0.f;
    float score21 = 0.f;
    float meanLogit = 0.f;
    float score = 0.f;
    float pt = -1.f;
    float eta = -1.f;
    float phi = -1.f;
    float mass = -1.f;
    float deltaR = -1.f;
    float dca = -1.f;
    float dcaSig = -1.f;
    float cptopv = -1.f;
    TwoTrackVertexFitter::Result fit;
  };

  std::uint64_t pairKey(unsigned int first, unsigned int second) {
    return (static_cast<std::uint64_t>(first) << 32) | second;
  }

  float sigmoid(float x) {
    if (x >= 0.f) return 1.f / (1.f + std::exp(-x));
    const float e = std::exp(x);
    return e / (1.f + e);
  }
}

class TrackGNNSeedProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackGNNSeedProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> trackPairToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> seederPairToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  unsigned int maxSeeds_;
  TwoTrackVertexFitter fitter_;
};

TrackGNNSeedProducer::TrackGNNSeedProducer(const edm::ParameterSet& config)
    : tracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("tracks"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("primaryVertices"))),
      trackPairToken_(consumes<nanoaod::FlatTable>(config.getParameter<edm::InputTag>("trackPairSrc"))),
      seederPairToken_(consumes<nanoaod::FlatTable>(config.getParameter<edm::InputTag>("seederPairSrc"))),
      ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxSeeds_(config.getParameter<unsigned int>("maxSeeds")) {
  produces<nanoaod::FlatTable>("Seed");
  produces<nanoaod::FlatTable>("SeedSummary");
}

void TrackGNNSeedProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("filteredTracks"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("filteredTracks"));
  desc.add<edm::InputTag>("trackPairSrc", edm::InputTag("trackPairFeatures", "trackPair"));
  desc.add<edm::InputTag>("seederPairSrc", edm::InputTag("trackGNNSeederInference", "trackPair"));
  desc.add<unsigned int>("maxSeeds", 128);
  descriptions.add("trackGNNSeedProducer", desc);
}

void TrackGNNSeedProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<reco::TrackCollection> tracks;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<nanoaod::FlatTable> pairs;
  edm::Handle<nanoaod::FlatTable> predictions;
  event.getByToken(tracksToken_, tracks);
  event.getByToken(primaryVerticesToken_, vertices);
  event.getByToken(trackPairToken_, pairs);
  event.getByToken(seederPairToken_, predictions);
  if (!tracks.isValid() || !pairs.isValid() || !predictions.isValid()) {
    throw cms::Exception("ProductNotFound") << "TrackGNNSeedProducer requires tracks, the original trackPair table, and its seeder score extension";
  }
  if (pairs->size() != predictions->size()) {
    throw cms::Exception("TableSizeMismatch") << "TrackGNNSeedProducer input table row mismatch: pair rows=" << pairs->size()
                                                << ", prediction rows=" << predictions->size();
  }

  const auto trk1 = pairs->columnData<int>(pairs->columnIndex("trk1Idx"));
  const auto trk2 = pairs->columnData<int>(pairs->columnIndex("trk2Idx"));
  const auto pairPt = pairs->columnData<float>(pairs->columnIndex("pt"));
  const auto pairEta = pairs->columnData<float>(pairs->columnIndex("eta"));
  const auto pairMass = pairs->columnData<float>(pairs->columnIndex("pairInvMass"));
  const auto pairDR = pairs->columnData<float>(pairs->columnIndex("deltaR"));
  const auto pairDca = pairs->columnData<float>(pairs->columnIndex("dca"));
  const auto pairDcaSig = pairs->columnData<float>(pairs->columnIndex("dcaSig"));
  const auto pairCptopv = pairs->columnData<float>(pairs->columnIndex("cptopv"));
  const auto logits = predictions->columnData<float>(predictions->columnIndex("seederLogit"));
  const auto scoreValid = predictions->columnData<int>(predictions->columnIndex("seederValid"));

  const size_t nTracks = tracks->size();
  const size_t nEdges = pairs->size();
  std::vector<Candidate> candidates;
  candidates.reserve(nEdges / 2);
  std::unordered_map<std::uint64_t, size_t> candidateByPair;
  candidateByPair.reserve(nEdges);
  std::vector<int> directedRows12, directedRows21;
  std::vector<float> directedLogit12, directedLogit21;
  std::vector<int> seen12, seen21;
  size_t nInvalidPairIndices = 0;
  size_t nInvalidScores = 0;
  size_t nDuplicateDirectedEdges = 0;

  // Preserve first-observed physical-pair order. TrackPairFeaturesProducer's
  // nested (i,j) loops emit the canonical i<j row first, in lexicographic
  // physical-pair order; filtering that directed stream to i<j has the same
  // enumeration. The lookup map only resolves pairs and never breaks ties.
  for (size_t row = 0; row < nEdges; ++row) {
    const int a = trk1[row], b = trk2[row];
    if (a < 0 || b < 0 || static_cast<size_t>(a) >= nTracks || static_cast<size_t>(b) >= nTracks || a == b) {
      ++nInvalidPairIndices;
      continue;
    }
    const bool rowScoreValid = scoreValid[row] == 1 && std::isfinite(logits[row]);
    if (!rowScoreValid) ++nInvalidScores;
    const int first = std::min(a, b);
    const int second = std::max(a, b);
    const std::uint64_t key = pairKey(first, second);
    auto found = candidateByPair.find(key);
    size_t index;
    if (found == candidateByPair.end()) {
      index = candidates.size();
      candidateByPair.emplace(key, index);
      Candidate c;
      c.first = first;
      c.second = second;
      const auto& p1 = (*tracks)[first];
      const auto& p2 = (*tracks)[second];
      c.phi = std::atan2(p1.py() + p2.py(), p1.px() + p2.px());
      candidates.push_back(c);
      directedRows12.push_back(-1); directedRows21.push_back(-1);
      directedLogit12.push_back(0.f); directedLogit21.push_back(0.f);
      seen12.push_back(0); seen21.push_back(0);
    } else {
      index = found->second;
    }
    const bool is12 = (a == first && b == second);
    int& seen = is12 ? seen12[index] : seen21[index];
    if (seen) {
      ++nDuplicateDirectedEdges;
      continue;  // deterministic policy: retain the first row in each direction
    }
    seen = 1;
    if (is12) directedRows12[index] = static_cast<int>(row);
    else directedRows21[index] = static_cast<int>(row);
    if (rowScoreValid) {
      if (is12) directedLogit12[index] = logits[row];
      else directedLogit21[index] = logits[row];
    } else {
      // Keep a separate invalid marker via a nonfinite value; score sentinels
      // in the inference extension are finite, so they cannot be mistaken as valid.
      if (is12) directedLogit12[index] = std::numeric_limits<float>::quiet_NaN();
      else directedLogit21[index] = std::numeric_limits<float>::quiet_NaN();
    }
  }

  size_t nMissingReciprocal = 0;
  std::vector<Candidate> validPairs;
  validPairs.reserve(candidates.size());
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (directedRows12[i] < 0 || directedRows21[i] < 0) {
      ++nMissingReciprocal;
      continue;
    }
    const float z12 = directedLogit12[i], z21 = directedLogit21[i];
    if (!std::isfinite(z12) || !std::isfinite(z21)) continue;
    Candidate c = candidates[i];
    c.edge12 = directedRows12[i]; c.edge21 = directedRows21[i];
    // Pair-level quantities are read from the canonical i->j row, independent
    // of which direction appeared first in the input table.
    const size_t canonicalRow = static_cast<size_t>(c.edge12);
    c.pt = pairPt[canonicalRow]; c.eta = pairEta[canonicalRow]; c.mass = pairMass[canonicalRow];
    c.deltaR = pairDR[canonicalRow]; c.dca = pairDca[canonicalRow];
    c.dcaSig = pairDcaSig[canonicalRow]; c.cptopv = pairCptopv[canonicalRow];
    c.logit12 = z12; c.logit21 = z21;
    c.score12 = sigmoid(z12); c.score21 = sigmoid(z21);
    c.meanLogit = 0.5f * (z12 + z21);
    c.score = sigmoid(c.meanLogit);
    if (!std::isfinite(c.score) || !std::isfinite(c.meanLogit)) continue;
    validPairs.push_back(c);
  }
  const size_t nUniquePairs = candidates.size();
  const size_t nValidScoredPairs = validPairs.size();
  std::stable_sort(validPairs.begin(), validPairs.end(), [](const Candidate& a, const Candidate& b) {
    return a.meanLogit > b.meanLogit;
  });
  const size_t nSeeds = std::min<size_t>(maxSeeds_, validPairs.size());
  validPairs.resize(nSeeds);

  const auto& builder = setup.getData(ttbToken_);
  std::vector<reco::TransientTrack> transientTracks(nTracks);
  std::vector<unsigned char> transientBuilt(nTracks, 0);
  std::vector<reco::TransientTrack> scratch;
  scratch.reserve(2);
  const reco::Vertex* pv = (vertices.isValid() && !vertices->empty() && !vertices->front().isFake())
                               ? &vertices->front()
                               : nullptr;
  size_t nFitSuccess = 0;
  size_t nCovarianceInvalid = 0;
  for (auto& seed : validPairs) {
    for (int index : {seed.first, seed.second}) {
      if (!transientBuilt[index]) {
        transientTracks[index] = builder.build((*tracks)[index]);
        transientBuilt[index] = 1;
      }
    }
    seed.fit = fitter_.fit(transientTracks[seed.first], transientTracks[seed.second],
                           (*tracks)[seed.first], (*tracks)[seed.second], pv, scratch);
    if (seed.fit.fitValid) {
      ++nFitSuccess;
      if (!seed.fit.covarianceValid) ++nCovarianceInvalid;
    }
  }

  std::vector<int> rank, seedTrk1, seedTrk2, edgeIdx12, edgeIdx21;
  std::vector<float> logit12, logit21, score12, score21, meanLogit, score;
  std::vector<float> outPt, outEta, outPhi, outMass, deltaR, dca, dcaSig, cptopv;
  std::vector<int> fitValid, fitStatus, covarianceValid, pvValid, lxyValid, signedLxyValid, l3dValid;
  std::vector<int> pointing2DValid, pointing3DValid, perpendicularValid;
  std::vector<float> chi2, ndof, chi2Ndof, probability, vx, vy, vz, covXX, covXY, covXZ, covYY, covYZ, covZZ;
  std::vector<float> sigmaX, sigmaY, sigmaZ, lxy, lxyErr, lxySig, signedLxy, signedLxySig;
  std::vector<float> l3d, l3dErr, l3dSig, cos2D, cos3D, perp, perpSig;
  for (size_t i = 0; i < nSeeds; ++i) {
    const auto& c = validPairs[i];
    const auto& r = c.fit;
    rank.push_back(static_cast<int>(i + 1)); seedTrk1.push_back(c.first); seedTrk2.push_back(c.second);
    edgeIdx12.push_back(c.edge12); edgeIdx21.push_back(c.edge21);
    logit12.push_back(c.logit12); logit21.push_back(c.logit21); score12.push_back(c.score12); score21.push_back(c.score21);
    meanLogit.push_back(c.meanLogit); score.push_back(c.score);
    outPt.push_back(c.pt); outEta.push_back(c.eta); outPhi.push_back(c.phi); outMass.push_back(c.mass);
    deltaR.push_back(c.deltaR); dca.push_back(c.dca); dcaSig.push_back(c.dcaSig); cptopv.push_back(c.cptopv);
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

  auto seedTable = std::make_unique<nanoaod::FlatTable>(nSeeds, "Seed", false);
#define ADD_SEED_COLUMN(TYPE, NAME, DATA, DOC) seedTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_SEED_COLUMN(int, "rank", rank, "stable descending mean-logit rank from 1; ties retain first physical-pair occurrence");
  ADD_SEED_COLUMN(int, "trk1Idx", seedTrk1, "lower filteredTracks index in this physical pair");
  ADD_SEED_COLUMN(int, "trk2Idx", seedTrk2, "higher filteredTracks index in this physical pair");
  ADD_SEED_COLUMN(int, "edgeIdx12", edgeIdx12, "original trackPair row for canonical trk1Idx to trk2Idx");
  ADD_SEED_COLUMN(int, "edgeIdx21", edgeIdx21, "original trackPair row for reciprocal direction");
  ADD_SEED_COLUMN(float, "logit12", logit12, "raw directed ONNX logit for edgeIdx12");
  ADD_SEED_COLUMN(float, "logit21", logit21, "raw directed ONNX logit for edgeIdx21");
  ADD_SEED_COLUMN(float, "score12", score12, "sigmoid of logit12");
  ADD_SEED_COLUMN(float, "score21", score21, "sigmoid of logit21");
  ADD_SEED_COLUMN(float, "meanLogit", meanLogit, "mean of reciprocal directed logits");
  ADD_SEED_COLUMN(float, "score", score, "sigmoid(meanLogit); ranking uses meanLogit");
  ADD_SEED_COLUMN(float, "pt", outPt, "pair transverse momentum from canonical trackPair row [GeV]");
  ADD_SEED_COLUMN(float, "eta", outEta, "pair pseudorapidity from canonical trackPair row");
  ADD_SEED_COLUMN(float, "phi", outPhi, "pair azimuth from summed canonical constituent momenta [rad]");
  ADD_SEED_COLUMN(float, "mass", outMass, "pair invariant mass from canonical trackPair row, pion mass hypothesis [GeV]");
  ADD_SEED_COLUMN(float, "deltaR", deltaR, "angular separation from canonical trackPair row");
  ADD_SEED_COLUMN(float, "dca", dca, "distance of closest approach from canonical trackPair row [cm]");
  ADD_SEED_COLUMN(float, "dcaSig", dcaSig, "DCA significance from canonical trackPair row");
  ADD_SEED_COLUMN(float, "cptopv", cptopv, "crossing-point distance from primary vertex from canonical trackPair row [cm]");
  ADD_SEED_COLUMN(int, "vtxFitValid", fitValid, "1 for valid selected-pair Kalman fit");
  ADD_SEED_COLUMN(int, "vtxFitStatus", fitStatus, "1 success, 2 invalid track, 3 invalid fitter result, 4 exception/numerics; 0 default");
  ADD_SEED_COLUMN(int, "vtxCovarianceValid", covarianceValid, "1 when fitted position covariance passes PSD validation");
  ADD_SEED_COLUMN(int, "vtxPVValid", pvValid, "1 when a usable leading reconstructed PV exists");
  ADD_SEED_COLUMN(int, "vtxLxyValid", lxyValid, "1 when Lxy and uncertainty are valid");
  ADD_SEED_COLUMN(int, "vtxSignedLxyValid", signedLxyValid, "1 when signed Lxy and significance are valid");
  ADD_SEED_COLUMN(int, "vtxL3DValid", l3dValid, "1 when 3D displacement and uncertainty are valid");
  ADD_SEED_COLUMN(int, "vtxPointing2DValid", pointing2DValid, "1 when transverse pointing cosine is valid");
  ADD_SEED_COLUMN(int, "vtxPointing3DValid", pointing3DValid, "1 when 3D pointing cosine is valid");
  ADD_SEED_COLUMN(int, "vtxPerpendicularValid", perpendicularValid, "1 when transverse perpendicular displacement is valid");
  ADD_SEED_COLUMN(float, "vtxChi2", chi2, "Kalman vertex-fit chi2");
  ADD_SEED_COLUMN(float, "vtxNdof", ndof, "Kalman vertex-fit degrees of freedom");
  ADD_SEED_COLUMN(float, "vtxChi2Ndof", chi2Ndof, "Kalman vertex-fit chi2/ndof");
  ADD_SEED_COLUMN(float, "vtxProb", probability, "chi-square fit probability");
  ADD_SEED_COLUMN(float, "vtxX", vx, "fitted vertex x [cm]"); ADD_SEED_COLUMN(float, "vtxY", vy, "fitted vertex y [cm]"); ADD_SEED_COLUMN(float, "vtxZ", vz, "fitted vertex z [cm]");
  ADD_SEED_COLUMN(float, "vtxCovXX", covXX, "fitted covariance xx [cm^2]"); ADD_SEED_COLUMN(float, "vtxCovXY", covXY, "fitted covariance xy [cm^2]"); ADD_SEED_COLUMN(float, "vtxCovXZ", covXZ, "fitted covariance xz [cm^2]");
  ADD_SEED_COLUMN(float, "vtxCovYY", covYY, "fitted covariance yy [cm^2]"); ADD_SEED_COLUMN(float, "vtxCovYZ", covYZ, "fitted covariance yz [cm^2]"); ADD_SEED_COLUMN(float, "vtxCovZZ", covZZ, "fitted covariance zz [cm^2]");
  ADD_SEED_COLUMN(float, "vtxSigmaX", sigmaX, "fitted x uncertainty [cm]"); ADD_SEED_COLUMN(float, "vtxSigmaY", sigmaY, "fitted y uncertainty [cm]"); ADD_SEED_COLUMN(float, "vtxSigmaZ", sigmaZ, "fitted z uncertainty [cm]");
  ADD_SEED_COLUMN(float, "vtxLxy", lxy, "PV-to-fit transverse displacement [cm]"); ADD_SEED_COLUMN(float, "vtxLxyErr", lxyErr, "Lxy uncertainty [cm]"); ADD_SEED_COLUMN(float, "vtxLxySig", lxySig, "Lxy significance");
  ADD_SEED_COLUMN(float, "vtxSignedLxy", signedLxy, "signed projected transverse displacement [cm]"); ADD_SEED_COLUMN(float, "vtxSignedLxySig", signedLxySig, "signed Lxy significance");
  ADD_SEED_COLUMN(float, "vtxL3D", l3d, "PV-to-fit 3D displacement [cm]"); ADD_SEED_COLUMN(float, "vtxL3DErr", l3dErr, "L3D uncertainty [cm]"); ADD_SEED_COLUMN(float, "vtxL3DSig", l3dSig, "L3D significance");
  ADD_SEED_COLUMN(float, "vtxCosPointing2D", cos2D, "cosine between transverse displacement and pair momentum");
  ADD_SEED_COLUMN(float, "vtxCosPointing3D", cos3D, "cosine between 3D displacement and pair momentum");
  ADD_SEED_COLUMN(float, "vtxPerpDisplacement", perp, "transverse displacement perpendicular to pair pT [cm]");
  ADD_SEED_COLUMN(float, "vtxPerpDisplacementSig", perpSig, "perpendicular displacement significance");
#undef ADD_SEED_COLUMN
  event.put(std::move(seedTable), "Seed");

  auto summary = std::make_unique<nanoaod::FlatTable>(1, "SeedSummary", true);
  const auto putCount = [&summary](const std::string& name, int value, const std::string& doc) {
    summary->addColumn<int>(name, std::vector<int>{value}, doc);
  };
  putCount("nTracks", static_cast<int>(nTracks), "number of rows in filteredTracks");
  putCount("nDirectedEdges", static_cast<int>(nEdges), "rows in the original directed trackPair table");
  putCount("nUniquePairs", static_cast<int>(nUniquePairs), "unique canonical pairs with valid, non-self track indices");
  putCount("nValidScoredPairs", static_cast<int>(nValidScoredPairs), "unique pairs with reciprocal rows and two finite valid logits");
  putCount("nSeeds", static_cast<int>(nSeeds), "selected top-K valid scored physical pairs");
  putCount("nFitAttempts", static_cast<int>(nSeeds), "number of selected pairs passed to the two-track fitter");
  putCount("nFitSuccess", static_cast<int>(nFitSuccess), "selected fits with fitValid=1");
  putCount("nFitFailure", static_cast<int>(nSeeds - nFitSuccess), "selected fits with fitValid=0; failed seeds remain in Seed");
  putCount("nCovarianceInvalid", static_cast<int>(nCovarianceInvalid), "fit-valid selected rows with covarianceValid=0");
  putCount("nMissingReciprocal", static_cast<int>(nMissingReciprocal), "unique valid-index pairs missing one directed row in trackPair");
  putCount("nInvalidPairIndices", static_cast<int>(nInvalidPairIndices), "directed rows with invalid, out-of-range, or self-pair indices");
  putCount("nInvalidScores", static_cast<int>(nInvalidScores), "valid-index directed rows with seederValid=0 or nonfinite logit, including duplicate rows");
  putCount("nDuplicateDirectedEdges", static_cast<int>(nDuplicateDirectedEdges), "extra rows duplicating an already-seen directed physical pair; first row retained");
  event.put(std::move(summary), "SeedSummary");
}

DEFINE_FWK_MODULE(TrackGNNSeedProducer);
