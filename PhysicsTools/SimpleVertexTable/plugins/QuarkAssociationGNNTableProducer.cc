// Joins the finalized Stage-2 reco->QGenCharged assignment to the historical
// TrackGVFeaturesProducer output. That producer remains the single source for
// AssociationGNN track selection, pair features, candidate edges and legacy
// hadron labels; this module adds QFrag-specific truth joins and pair targets.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <cmath>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace {
  template <typename T>
  std::span<const T> column(const nanoaod::FlatTable& table, const std::string& name) {
    const int idx = table.columnIndex(name);
    if (idx < 0)
      throw cms::Exception("ProductNotFound") << "FlatTable is missing required column '" << name << "'";
    return table.columnData<T>(static_cast<unsigned int>(idx));
  }

  bool isBHadron(int pdg) {
    switch (std::abs(pdg)) {
      case 511: case 521: case 531: case 541: case 5122: case 5132: case 5232:
      case 5332: case 5142: case 5242: case 5342: case 5512: case 5532:
      case 5542: case 5554: return true;
      default: return false;
    }
  }

  bool isDHadron(int pdg) {
    switch (std::abs(pdg)) {
      case 411: case 421: case 431: case 4122: case 4232: case 4132:
      case 4332: case 4412: case 4422: case 4432: case 4444: return true;
      default: return false;
    }
  }

  // -1 unresolved, 0 other, 1 B, 2 prompt D, 3 D from B, 4 D origin unknown.
  int hadronClass(int pdg, int charmOrigin) {
    if (isBHadron(pdg)) return 1;
    if (!isDHadron(pdg)) return 0;
    if (charmOrigin == 1) return 2;
    if (charmOrigin == 2) return 3;
    return 4;
  }
}

class QuarkAssociationGNNTableProducer : public edm::stream::EDProducer<> {
public:
  explicit QuarkAssociationGNNTableProducer(const edm::ParameterSet& cfg)
      : recoToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("recoTrackTable"))),
        chargedToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("genChargedTable"))),
        hadronToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("hadronTable"))),
        pairToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("pairTable"))) {
    produces<nanoaod::FlatTable>("QGNNTrack");
    produces<nanoaod::FlatTable>("QGNNPairTruth");
    produces<nanoaod::FlatTable>("QGNNHadron");
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& reco = event.get(recoToken_);
    const auto& charged = event.get(chargedToken_);
    const auto& hadrons = event.get(hadronToken_);
    const auto& pairs = event.get(pairToken_);

    const auto rTrackIdx = column<int>(reco, "trackIdx");
    const auto rGenIdx = column<int>(reco, "genChargedIdx");
    const auto rPt = column<float>(reco, "pt");
    const auto rEta = column<float>(reco, "eta");
    const auto gHadronIdx = column<int>(charged, "hadronIdx");
    const auto gQuarkIdx = column<int>(charged, "quarkIdx");
    const auto gQuality = column<int>(charged, "associationQuality");
    const auto gMembership = column<int>(charged, "membershipType");
    const auto gSourceCertain = column<int>(charged, "fragmentationSourceCertain");
    const auto hPdgId = column<int>(hadrons, "pdgId");
    const auto hCharmOrigin = column<int>(hadrons, "charmOrigin");
    const auto edgeSrc = column<int>(pairs, "trk1Idx");
    const auto edgeDst = column<int>(pairs, "trk2Idx");
    const size_t nTracks = reco.nRows(), nGen = charged.nRows(), nHadrons = hadrons.nRows(), nPairs = pairs.nRows();
    if (rTrackIdx.size() != nTracks || rGenIdx.size() != nTracks || rPt.size() != nTracks || rEta.size() != nTracks ||
        gHadronIdx.size() != nGen || gQuarkIdx.size() != nGen || gQuality.size() != nGen || gMembership.size() != nGen ||
        gSourceCertain.size() != nGen || hPdgId.size() != nHadrons || hCharmOrigin.size() != nHadrons ||
        edgeSrc.size() != nPairs || edgeDst.size() != nPairs)
      throw cms::Exception("LogicError") << "QFrag GNN tables have inconsistent row counts";

    // One QGNNTrack row for every QRecoTrack row, in exactly the unpacker order.
    std::vector<int> trackIdx, genIdx, hadronIdx, quarkIdx, quality, membership, sourceCertain;
    std::vector<int> hadClass, hadTruthValid, isB, isPromptD, isDFromB;
    std::vector<int> nAcceptedDecayTracks(nHadrons, 0);
    trackIdx.reserve(nTracks); genIdx.reserve(nTracks); hadronIdx.reserve(nTracks); quarkIdx.reserve(nTracks);
    quality.reserve(nTracks); membership.reserve(nTracks); sourceCertain.reserve(nTracks);
    hadClass.reserve(nTracks); hadTruthValid.reserve(nTracks); isB.reserve(nTracks); isPromptD.reserve(nTracks); isDFromB.reserve(nTracks);
    for (size_t i = 0; i < nTracks; ++i) {
      const int ig = rGenIdx[i];
      const bool validGen = ig >= 0 && static_cast<size_t>(ig) < nGen;
      const int ih = validGen ? gHadronIdx[ig] : -1;
      const bool validHadron = ih >= 0 && static_cast<size_t>(ih) < nHadrons;
      const bool knownDecay = validGen && validHadron && gMembership[ig] == 1;
      const int cls = knownDecay ? hadronClass(hPdgId[ih], hCharmOrigin[ih]) : -1;
      trackIdx.push_back(rTrackIdx[i]); genIdx.push_back(ig);
      hadronIdx.push_back(validHadron ? ih : -1);
      quarkIdx.push_back(validGen ? gQuarkIdx[ig] : -1);
      quality.push_back(validGen ? gQuality[ig] : 0);
      membership.push_back(validGen ? gMembership[ig] : 0);
      sourceCertain.push_back(validGen ? gSourceCertain[ig] : 0);
      hadClass.push_back(cls); hadTruthValid.push_back(knownDecay ? 1 : 0);
      isB.push_back(knownDecay && cls == 1 ? 1 : 0);
      isPromptD.push_back(knownDecay && cls == 2 ? 1 : 0);
      isDFromB.push_back(knownDecay && cls == 3 ? 1 : 0);
      // The Stage-3 accepted-decay-track definition; each reco row is unique
      // because the Stage-2 assignment is one-to-one.
      if (knownDecay && rPt[i] >= 0.4f && std::abs(rEta[i]) <= 2.5f)
        ++nAcceptedDecayTracks[ih];
    }

    auto trackTable = std::make_unique<nanoaod::FlatTable>(nTracks, "QGNNTrack", false, false);
    trackTable->addColumn<int>("trackIdx", trackIdx, "event-local row in unpackedTracksAndVertices, QRecoTrack and base track table");
    trackTable->addColumn<int>("genChargedIdx", genIdx, "QGenCharged row, or -1 when unmatched");
    trackTable->addColumn<int>("hadronIdx", hadronIdx, "event-local QHadron row, or -1");
    trackTable->addColumn<int>("quarkIdx", quarkIdx, "event-local QFrag row, or -1");
    trackTable->addColumn<int>("associationQuality", quality, "Stage-1 ancestry quality: 0 missing, 1 unique, 2 ambiguous, 3 unknown");
    trackTable->addColumn<int>("membershipType", membership, "0 unassigned, 1 heavy-hadron decay, 2 additional candidate, 3 ambiguous");
    trackTable->addColumn<int>("fragmentationSourceCertain", sourceCertain, "Stage-1 physical-source certainty; additional uncertain candidates remain 0");
    trackTable->addColumn<int>("hadronClass", hadClass, "-1 unresolved, 0 other, 1 B, 2 prompt D, 3 D from B, 4 D origin unknown");
    trackTable->addColumn<int>("hadronTruthValid", hadTruthValid, "1 only for membershipType=1 association to a valid QHadron; quark quality is retained separately");
    trackTable->addColumn<int>("isB", isB, "B decay-track target; valid only when hadronTruthValid=1");
    trackTable->addColumn<int>("isPromptD", isPromptD, "prompt-D decay-track target; valid only when hadronTruthValid=1");
    trackTable->addColumn<int>("isDFromB", isDFromB, "D-from-B decay-track target; valid only when hadronTruthValid=1");
    event.put(std::move(trackTable), "QGNNTrack");

    std::vector<int> hadIdx, hadN, hadCategory, hadClassOut;
    for (size_t i = 0; i < nHadrons; ++i) {
      const int n = nAcceptedDecayTracks[i];
      hadIdx.push_back(static_cast<int>(i)); hadN.push_back(n);
      hadCategory.push_back(n == 0 ? 0 : n == 1 ? 1 : n == 2 ? 2 : 3);
      hadClassOut.push_back(hadronClass(hPdgId[i], hCharmOrigin[i]));
    }
    auto hadronTable = std::make_unique<nanoaod::FlatTable>(nHadrons, "QGNNHadron", false, false);
    hadronTable->addColumn<int>("hadronIdx", hadIdx, "row in QHadron");
    hadronTable->addColumn<int>("nAcceptedDecayTracks", hadN, "distinct matched membershipType=1 decay tracks, pT>=0.4 GeV and |eta|<=2.5");
    hadronTable->addColumn<int>("nCategory", hadCategory, "0=N0, 1=N1, 2=N2, 3=N3+");
    hadronTable->addColumn<int>("hadronClass", hadClassOut, "-1 unresolved, 0 other, 1 B, 2 prompt D, 3 D from B, 4 D origin unknown");
    event.put(std::move(hadronTable), "QGNNHadron");

    // This table follows the exact row order of TrackGVFeaturesProducer's
    // one-row-per-physical-pair trackPair table. Legacy edgeLabel remains there.
    std::vector<int> pairIdx, hadronTruth, sameHadron, sameB, samePromptD, sameDFromB;
    std::vector<int> quarkTruth, sameQuark;
    for (size_t i = 0; i < nPairs; ++i) {
      const int ti = edgeSrc[i], tj = edgeDst[i];
      if (ti < 0 || tj < 0 || static_cast<size_t>(ti) >= nTracks || static_cast<size_t>(tj) >= nTracks)
        throw cms::Exception("LogicError") << "trackPair contains an out-of-range event-local track index";
      const int gi = rGenIdx[ti], gj = rGenIdx[tj];
      const bool validGi = gi >= 0 && static_cast<size_t>(gi) < nGen;
      const bool validGj = gj >= 0 && static_cast<size_t>(gj) < nGen;
      const int hi = validGi ? gHadronIdx[gi] : -1, hj = validGj ? gHadronIdx[gj] : -1;
      const bool knownHi = validGi && gMembership[gi] == 1 && hi >= 0 && static_cast<size_t>(hi) < nHadrons;
      const bool knownHj = validGj && gMembership[gj] == 1 && hj >= 0 && static_cast<size_t>(hj) < nHadrons;
      const bool validH = knownHi && knownHj, sameH = validH && hi == hj;
      const int ci = knownHi ? hadronClass(hPdgId[hi], hCharmOrigin[hi]) : -1;
      const int qi = validGi ? gQuarkIdx[gi] : -1, qj = validGj ? gQuarkIdx[gj] : -1;
      const bool validQ = validGi && validGj && gQuality[gi] == 1 && gQuality[gj] == 1 &&
                          gSourceCertain[gi] == 1 && gSourceCertain[gj] == 1 && qi >= 0 && qj >= 0;
      pairIdx.push_back(static_cast<int>(i));
      hadronTruth.push_back(validH ? 1 : 0); sameHadron.push_back(validH ? (sameH ? 1 : 0) : -1);
      sameB.push_back(validH ? (sameH && ci == 1 ? 1 : 0) : -1);
      samePromptD.push_back(validH ? (sameH && ci == 2 ? 1 : 0) : -1);
      sameDFromB.push_back(validH ? (sameH && ci == 3 ? 1 : 0) : -1);
      quarkTruth.push_back(validQ ? 1 : 0); sameQuark.push_back(validQ ? (qi == qj ? 1 : 0) : -1);
    }
    auto pairTruthTable = std::make_unique<nanoaod::FlatTable>(nPairs, "QGNNPairTruth", false, false);
    pairTruthTable->addColumn<int>("pairIdx", pairIdx, "row in trackPair / historical AssociationGNN physical-edge table");
    pairTruthTable->addColumn<int>("hadronTruthValid", hadronTruth, "both endpoints have membershipType=1 association to valid QHadron rows");
    pairTruthTable->addColumn<int>("sameHadron", sameHadron, "1 same QHadron, 0 different known QHadrons, -1 unknown");
    pairTruthTable->addColumn<int>("sameB", sameB, "same B-hadron pair label; -1 when truth invalid");
    pairTruthTable->addColumn<int>("samePromptD", samePromptD, "same prompt-D pair label; -1 when truth invalid");
    pairTruthTable->addColumn<int>("sameDFromB", sameDFromB, "same D-from-B pair label; -1 when truth invalid");
    pairTruthTable->addColumn<int>("quarkTruthValid", quarkTruth, "both endpoints have unique, source-certain QFrag membership");
    pairTruthTable->addColumn<int>("sameQuark", sameQuark, "same QFrag label; -1 when certainty is insufficient");
    event.put(std::move(pairTruthTable), "QGNNPairTruth");
  }

private:
  edm::EDGetTokenT<nanoaod::FlatTable> recoToken_, chargedToken_, hadronToken_, pairToken_;
};

DEFINE_FWK_MODULE(QuarkAssociationGNNTableProducer);
