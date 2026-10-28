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
#include <limits>
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
        quarkToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("quarkTable"))),
        legacyHadronToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("legacyHadronTable"))),
        pairToken_(consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("pairTable"))),
        produceFullValidationTables_(cfg.getParameter<bool>("produceFullValidationTables")),
        saveFragmentationObservables_(cfg.getUntrackedParameter<bool>("saveFragmentationObservables", true)) {
    if (produceFullValidationTables_) {
      produces<nanoaod::FlatTable>("QGNNTrack");
      produces<nanoaod::FlatTable>("QGNNPairTruth");
      produces<nanoaod::FlatTable>("QGNNHadron");
    }
    produces<nanoaod::FlatTable>("QTrainTrack");
    produces<nanoaod::FlatTable>("QTrainPair");
    produces<nanoaod::FlatTable>("QTrainHadron");
    produces<nanoaod::FlatTable>("QTrainQuark");
    produces<nanoaod::FlatTable>("QTrainLegacyHadron");
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& reco = event.get(recoToken_);
    const auto& charged = event.get(chargedToken_);
    const auto& hadrons = event.get(hadronToken_);
    const auto& quarks = event.get(quarkToken_);
    const auto& legacyHadrons = event.get(legacyHadronToken_);
    const auto& pairs = event.get(pairToken_);

    const auto rTrackIdx = column<int>(reco, "trackIdx");
    const auto rGenIdx = column<int>(reco, "genChargedIdx");
    const auto rPt = column<float>(reco, "pt");
    const auto rEta = column<float>(reco, "eta");
    const auto rMatchDR = column<float>(reco, "matchDR");
    const auto rMatchRelPt = column<float>(reco, "matchRelPt");
    const auto rMatchStatus = column<int>(reco, "matchStatus");
    const auto rNCandidates = column<int>(reco, "nGenCandidates");
    const auto gHadronIdx = column<int>(charged, "hadronIdx");
    const auto gQuarkIdx = column<int>(charged, "quarkIdx");
    const auto gQuality = column<int>(charged, "associationQuality");
    const auto gMembership = column<int>(charged, "membershipType");
    const auto gSourceCertain = column<int>(charged, "fragmentationSourceCertain");
    const auto hPdgId = column<int>(hadrons, "pdgId");
    const auto hCharmOrigin = column<int>(hadrons, "charmOrigin");
    const auto hPt = column<float>(hadrons, "pt");
    const auto hEta = column<float>(hadrons, "eta");
    const auto hPhi = column<float>(hadrons, "phi");
    const auto hParentBIdx = column<int>(hadrons, "parentBHadronIdx");
    const auto hQuarkIdx = column<int>(hadrons, "quarkIdx");
    const auto hAssocQuality = column<int>(hadrons, "associationQuality");
    const auto hSourceCertain = column<int>(hadrons, "fragmentationSourceCertain");
    std::span<const float> hPreHadPt,hPreHadEta,hPreHadPhi;
    std::span<const int> hPreHadValid;
    if (saveFragmentationObservables_) {
      hPreHadPt = column<float>(hadrons, "preHadPt");
      hPreHadEta = column<float>(hadrons, "preHadEta");
      hPreHadPhi = column<float>(hadrons, "preHadPhi");
      hPreHadValid = column<int>(hadrons, "preHadKinematicsValid");
    }
    const auto qPdgId = column<int>(quarks, "pdgId");
    const auto qPx = column<float>(quarks, "px"), qPy = column<float>(quarks, "py"), qPz = column<float>(quarks, "pz"), qEnergy = column<float>(quarks, "energy");
    const auto qOriginCategory = column<int>(quarks, "originCategory"), qOriginValid = column<int>(quarks, "originCategoryValid");
    const auto legacyPdgId = column<int>(legacyHadrons, "Hadron_pdgId");
    const auto legacyIsB = column<int>(legacyHadrons, "isB");
    const auto legacyIsD = column<int>(legacyHadrons, "isD");
    const auto legacyFromHF = column<int>(legacyHadrons, "fromHF");
    const auto edgeSrc = column<int>(pairs, "trk1Idx");
    const auto edgeDst = column<int>(pairs, "trk2Idx");
    const size_t nTracks = reco.nRows(), nGen = charged.nRows(), nHadrons = hadrons.nRows(), nQuarks = quarks.nRows(), nPairs = pairs.nRows();
    if (rTrackIdx.size() != nTracks || rGenIdx.size() != nTracks || rPt.size() != nTracks || rEta.size() != nTracks ||
        rMatchDR.size() != nTracks || rMatchRelPt.size() != nTracks || rMatchStatus.size() != nTracks || rNCandidates.size() != nTracks ||
        gHadronIdx.size() != nGen || gQuarkIdx.size() != nGen || gQuality.size() != nGen || gMembership.size() != nGen ||
        gSourceCertain.size() != nGen || hPdgId.size() != nHadrons || hCharmOrigin.size() != nHadrons ||
        hPt.size() != nHadrons || hEta.size() != nHadrons || hPhi.size() != nHadrons ||
        hParentBIdx.size() != nHadrons || hQuarkIdx.size() != nHadrons || hAssocQuality.size() != nHadrons || hSourceCertain.size() != nHadrons ||
        (saveFragmentationObservables_ && (hPreHadPt.size() != nHadrons || hPreHadEta.size() != nHadrons || hPreHadPhi.size() != nHadrons || hPreHadValid.size() != nHadrons)) ||
        qPdgId.size() != nQuarks || qPx.size() != nQuarks || qPy.size() != nQuarks || qPz.size() != nQuarks || qEnergy.size() != nQuarks ||
        qOriginCategory.size() != nQuarks || qOriginValid.size() != nQuarks || legacyPdgId.size() != legacyHadrons.nRows() ||
        legacyIsB.size() != legacyHadrons.nRows() || legacyIsD.size() != legacyHadrons.nRows() ||
        legacyFromHF.size() != legacyHadrons.nRows() ||
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

    if (produceFullValidationTables_) {
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
      trackTable->addColumn<int>("isDFromB", isDFromB, "D-from-B target; valid only when hadronTruthValid=1");
      event.put(std::move(trackTable), "QGNNTrack");
    }

    std::vector<int> hadIdx, hadN, hadCategory, hadClassOut;
    for (size_t i = 0; i < nHadrons; ++i) {
      const int n = nAcceptedDecayTracks[i];
      hadIdx.push_back(static_cast<int>(i)); hadN.push_back(n);
      hadCategory.push_back(n == 0 ? 0 : n == 1 ? 1 : n == 2 ? 2 : 3);
      hadClassOut.push_back(hadronClass(hPdgId[i], hCharmOrigin[i]));
    }
    if (produceFullValidationTables_) {
      auto hadronTable = std::make_unique<nanoaod::FlatTable>(nHadrons, "QGNNHadron", false, false);
      hadronTable->addColumn<int>("hadronIdx", hadIdx, "row in QHadron");
      hadronTable->addColumn<int>("nAcceptedDecayTracks", hadN, "distinct matched membershipType=1 decay tracks, pT>=0.4 GeV and |eta|<=2.5");
      hadronTable->addColumn<int>("nCategory", hadCategory, "0=N0, 1=N1, 2=N2, 3=N3+");
      hadronTable->addColumn<int>("hadronClass", hadClassOut, "-1 unresolved, 0 other, 1 B, 2 prompt D, 3 D from B, 4 D origin unknown");
      event.put(std::move(hadronTable), "QGNNHadron");
    }

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
      if (produceFullValidationTables_) pairIdx.push_back(static_cast<int>(i));
      hadronTruth.push_back(validH ? 1 : 0); sameHadron.push_back(validH ? (sameH ? 1 : 0) : -1);
      sameB.push_back(validH ? (sameH && ci == 1 ? 1 : 0) : -1);
      samePromptD.push_back(validH ? (sameH && ci == 2 ? 1 : 0) : -1);
      sameDFromB.push_back(validH ? (sameH && ci == 3 ? 1 : 0) : -1);
      quarkTruth.push_back(validQ ? 1 : 0); sameQuark.push_back(validQ ? (qi == qj ? 1 : 0) : -1);
    }
    if (produceFullValidationTables_) {
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

    // Compact training tables contain only truth targets and references. The
    // original track/trackPair feature tables remain the single source of all
    // reconstructed inputs, so this adds no duplicate node or edge features.
    std::vector<int> trainTrackIdx, trainHadronIdx, trainQuarkIdx, trainMatchStatus, trainNCandidates;
    std::vector<int> trainQuality, trainMembership, trainSourceCertain, trainHadronClass;
    std::vector<int> trainHadronValid, trainIsB, trainIsPromptD, trainIsDFromB, trainIsDUnknown;
    std::vector<int> trainQuarkWeakValid, trainQuarkStrictValid;
    std::vector<float> trainMatchDR, trainMatchRelPt;
    for (size_t i = 0; i < nTracks; ++i) {
      const int ig = genIdx[i];
      const bool validGen = ig >= 0 && static_cast<size_t>(ig) < nGen;
      const int iq = validGen ? gQuarkIdx[ig] : -1;
      const bool validQuark = iq >= 0 && static_cast<size_t>(iq) < nQuarks;
      const int ih = hadronIdx[i];
      const bool validHad = ih >= 0 && static_cast<size_t>(ih) < nHadrons;
      const bool uniqueLineage = validGen && validQuark && quality[i] == 1;
      const bool weakLineage = uniqueLineage && (membership[i] == 1 || membership[i] == 2);
      const bool strictLineage = weakLineage && sourceCertain[i] == 1;
      trainTrackIdx.push_back(rTrackIdx[i]);
      trainHadronIdx.push_back(validHad ? ih : -1);
      trainQuarkIdx.push_back(validQuark ? iq : -1);
      trainMatchStatus.push_back(rMatchStatus[i]); trainMatchDR.push_back(rMatchDR[i]);
      trainMatchRelPt.push_back(rMatchRelPt[i]); trainNCandidates.push_back(rNCandidates[i]);
      trainQuality.push_back(quality[i]); trainMembership.push_back(membership[i]);
      trainSourceCertain.push_back(sourceCertain[i]); trainHadronClass.push_back(hadClass[i]);
      trainHadronValid.push_back(hadTruthValid[i]);
      trainIsB.push_back(isB[i]); trainIsPromptD.push_back(isPromptD[i]);
      trainIsDFromB.push_back(isDFromB[i]); trainIsDUnknown.push_back(hadTruthValid[i] && hadClass[i] == 4 ? 1 : 0);
      trainQuarkWeakValid.push_back(weakLineage ? 1 : 0);
      trainQuarkStrictValid.push_back(strictLineage ? 1 : 0);
    }
    auto trainTrackTable = std::make_unique<nanoaod::FlatTable>(nTracks, "QTrainTrack", false, false);
    trainTrackTable->addColumn<int>("trackIdx", trainTrackIdx, "row in the shared track feature table");
    trainTrackTable->addColumn<int>("hadronIdx", trainHadronIdx, "row in QTrainHadron, or -1");
    trainTrackTable->addColumn<int>("quarkIdx", trainQuarkIdx, "row in QTrainQuark, or -1");
    trainTrackTable->addColumn<int>("matchStatus", trainMatchStatus, "Stage-2 status: 0 no candidate, 1 unopposed, 2 competing, 3 lost one-to-one conflict");
    trainTrackTable->addColumn<float>("matchDR", trainMatchDR, "Stage-2 reco-to-generator deltaR, or -1");
    trainTrackTable->addColumn<float>("matchRelPt", trainMatchRelPt, "Stage-2 relative-pT residual, or -1");
    trainTrackTable->addColumn<int>("nGenCandidates", trainNCandidates, "admissible Stage-2 generator candidates before global assignment");
    trainTrackTable->addColumn<int>("associationQuality", trainQuality, "0 missing, 1 unique, 2 ambiguous, 3 unknown/incomplete ancestry");
    trainTrackTable->addColumn<int>("membershipType", trainMembership, "0 unassigned, 1 heavy-hadron decay, 2 additional candidate, 3 ambiguous");
    trainTrackTable->addColumn<int>("fragmentationSourceCertain", trainSourceCertain, "physical fragmentation-source certainty");
    trainTrackTable->addColumn<int>("quarkWeakValid", trainQuarkWeakValid, "unique QFrag ancestry for heavy-decay or additional candidate membership");
    trainTrackTable->addColumn<int>("quarkStrictValid", trainQuarkStrictValid, "unique QFrag ancestry with source certainty");
    trainTrackTable->addColumn<int>("hadronClass", trainHadronClass, "-1 unresolved, 0 other, 1 B, 2 prompt D, 3 D from B, 4 D origin unknown");
    trainTrackTable->addColumn<int>("hadronTruthValid", trainHadronValid, "valid heavy-hadron decay-track label mask");
    trainTrackTable->addColumn<int>("isB", trainIsB, "B decay-track label");
    trainTrackTable->addColumn<int>("isPromptD", trainIsPromptD, "prompt-D decay-track label");
    trainTrackTable->addColumn<int>("isDFromB", trainIsDFromB, "D-from-B decay-track label");
    trainTrackTable->addColumn<int>("isDUnknown", trainIsDUnknown, "D-origin-unknown decay-track label");
    event.put(std::move(trainTrackTable), "QTrainTrack");

    std::vector<float> zProj(nHadrons, std::numeric_limits<float>::quiet_NaN());
    std::vector<float> zT(nHadrons, std::numeric_limits<float>::quiet_NaN());
    std::vector<float> deltaR(nHadrons, std::numeric_limits<float>::quiet_NaN());
    std::vector<float> kTrel(nHadrons, std::numeric_limits<float>::quiet_NaN());
    std::vector<int> auxInputsValid(nHadrons, 0), auxPromptValid(nHadrons, 0), bFragmentationTargetValid(nHadrons, 0);
    std::vector<int> trainHadIdx, trainHadPdg, trainHadCharmOrigin, trainHadParentB;
    std::vector<int> trainHadQuarkIdx, trainHadAssocQuality, trainHadSourceCertain, trainHadPreValid;
    std::vector<float> trainHadPt, trainHadEta, trainHadPhi;
    std::vector<float> trainHadPrePt, trainHadPreEta, trainHadPrePhi;
    trainHadIdx.reserve(nHadrons); trainHadPdg.reserve(nHadrons); trainHadCharmOrigin.reserve(nHadrons);
    trainHadParentB.reserve(nHadrons); trainHadQuarkIdx.reserve(nHadrons); trainHadAssocQuality.reserve(nHadrons);
    trainHadSourceCertain.reserve(nHadrons); trainHadPreValid.reserve(nHadrons);
    trainHadPt.reserve(nHadrons); trainHadEta.reserve(nHadrons); trainHadPhi.reserve(nHadrons);
    if (saveFragmentationObservables_) { trainHadPrePt.reserve(nHadrons); trainHadPreEta.reserve(nHadrons); trainHadPrePhi.reserve(nHadrons); }
    constexpr double pi = 3.14159265358979323846;
    for (size_t ih = 0; ih < nHadrons; ++ih) {
      const int iq = hQuarkIdx[ih];
      const bool validQ = iq >= 0 && static_cast<size_t>(iq) < nQuarks;
      if (saveFragmentationObservables_) {
        const bool validAncestry = validQ && hPreHadValid[ih] == 1 && hAssocQuality[ih] == 1;
        const bool isD = isDHadron(hPdgId[ih]);
        const bool isB = isBHadron(hPdgId[ih]);
        const bool dInputsEligible = isD && validAncestry;
        const bool bInputsEligible = isB && validAncestry && hSourceCertain[ih] == 1 && std::abs(qPdgId[iq]) == 5;
        bool finiteInputs = dInputsEligible || bInputsEligible;
        const double pdPt = hPt[ih], pdEta = hEta[ih], pdPhi = hPhi[ih];
        const double pqPt = hPreHadPt[ih], pqEta = hPreHadEta[ih], pqPhi = hPreHadPhi[ih];
        finiteInputs = finiteInputs && std::isfinite(pdPt) && std::isfinite(pdEta) && std::isfinite(pdPhi) &&
                       std::isfinite(pqPt) && std::isfinite(pqEta) && std::isfinite(pqPhi) && pdPt > 0. && pqPt > 0.;
        if (finiteInputs) {
        const double pdx = pdPt * std::cos(pdPhi), pdy = pdPt * std::sin(pdPhi), pdz = pdPt * std::sinh(pdEta);
        const double pqx = pqPt * std::cos(pqPhi), pqy = pqPt * std::sin(pqPhi), pqz = pqPt * std::sinh(pqEta);
        const double q2 = pqx*pqx + pqy*pqy + pqz*pqz, qnorm = std::sqrt(q2);
        const double dot = pdx*pqx + pdy*pqy + pdz*pqz;
        const double cx = pdy*pqz - pdz*pqy, cy = pdz*pqx - pdx*pqz, cz = pdx*pqy - pdy*pqx;
        const double dphi = std::remainder(pdPhi - pqPhi, 2. * pi);
        const double dr = std::hypot(pdEta - pqEta, dphi);
        const double proj = dot / q2, trans = pdPt / pqPt, kt = std::sqrt(cx*cx + cy*cy + cz*cz) / qnorm;
        if (qnorm > 0. && std::isfinite(qnorm) && std::isfinite(pdx) && std::isfinite(pdy) && std::isfinite(pdz) &&
            std::isfinite(pqx) && std::isfinite(pqy) && std::isfinite(pqz) && std::isfinite(proj) &&
            std::isfinite(trans) && std::isfinite(dr) && std::isfinite(kt)) {
          zProj[ih] = static_cast<float>(proj); zT[ih] = static_cast<float>(trans);
          deltaR[ih] = static_cast<float>(dr); kTrel[ih] = static_cast<float>(kt);
          if (dInputsEligible) {
            // Preserve the existing D-only meaning of these two masks.
            auxInputsValid[ih] = 1;
            auxPromptValid[ih] = hCharmOrigin[ih] == 1 && std::abs(qPdgId[iq]) == 4 ? 1 : 0;
          }
          if (bInputsEligible) bFragmentationTargetValid[ih] = 1;
        }
        }
      }
      trainHadIdx.push_back(static_cast<int>(ih)); trainHadPdg.push_back(hPdgId[ih]);
      trainHadCharmOrigin.push_back(hCharmOrigin[ih]); trainHadParentB.push_back(hParentBIdx[ih]);
      trainHadQuarkIdx.push_back(validQ ? iq : -1); trainHadAssocQuality.push_back(hAssocQuality[ih]);
      trainHadSourceCertain.push_back(hSourceCertain[ih]);
      trainHadPt.push_back(hPt[ih]); trainHadEta.push_back(hEta[ih]); trainHadPhi.push_back(hPhi[ih]);
      if (saveFragmentationObservables_) {
        trainHadPreValid.push_back(hPreHadValid[ih]);
        trainHadPrePt.push_back(hPreHadPt[ih]); trainHadPreEta.push_back(hPreHadEta[ih]); trainHadPrePhi.push_back(hPreHadPhi[ih]);
      }
    }
    auto trainHadronTable = std::make_unique<nanoaod::FlatTable>(nHadrons, "QTrainHadron", false, false);
    trainHadronTable->addColumn<int>("hadronIdx", trainHadIdx, "event-local row; QTrainTrack.hadronIdx references this table");
    trainHadronTable->addColumn<int>("pdgId", trainHadPdg, "signed hadron PDG identifier");
    trainHadronTable->addColumn<int>("charmOrigin", trainHadCharmOrigin, "0 not D, 1 prompt D, 2 D from B, 3 unknown D origin");
    trainHadronTable->addColumn<float>("pt", trainHadPt, "generator hadron pT [GeV]");
    trainHadronTable->addColumn<float>("eta", trainHadEta, "generator hadron eta");
    trainHadronTable->addColumn<float>("phi", trainHadPhi, "generator hadron phi [rad]");
    trainHadronTable->addColumn<int>("parentBHadronIdx", trainHadParentB, "QTrainHadron row of parent B, or -1");
    trainHadronTable->addColumn<int>("quarkIdx", trainHadQuarkIdx, "QTrainQuark row, or -1");
    trainHadronTable->addColumn<int>("associationQuality", trainHadAssocQuality, "hadron-to-quark ancestry quality");
    trainHadronTable->addColumn<int>("fragmentationSourceCertain", trainHadSourceCertain, "physical fragmentation-source certainty");
    trainHadronTable->addColumn<int>("nAcceptedDecayTracks", hadN, "matched heavy-decay tracks with pT>=0.4 GeV and |eta|<=2.5");
    trainHadronTable->addColumn<int>("nCategory", hadCategory, "0=N0, 1=N1, 2=N2, 3=N3+");
    if (saveFragmentationObservables_) {
      trainHadronTable->addColumn<float>("preHadPt", trainHadPrePt, "pre-hadronization quark proxy pT [GeV]");
      trainHadronTable->addColumn<float>("preHadEta", trainHadPreEta, "pre-hadronization quark proxy eta");
      trainHadronTable->addColumn<float>("preHadPhi", trainHadPrePhi, "pre-hadronization quark proxy phi [rad]");
      trainHadronTable->addColumn<int>("preHadKinematicsValid", trainHadPreValid, "pre-hadronization kinematics-valid flag");
      trainHadronTable->addColumn<float>("zProj", zProj, "dot(pH,pQ)/|pQ|^2; D or B training-only target, NaN when unavailable");
      trainHadronTable->addColumn<float>("zT", zT, "pT(H)/pT(pre-hadronization quark); D or B target, NaN when unavailable");
      trainHadronTable->addColumn<float>("deltaR_DQuark", deltaR, "deltaR(H, pre-hadronization quark) for D or B, NaN when unavailable");
      trainHadronTable->addColumn<float>("kTrel", kTrel, "|pH cross pQ|/|pQ| [GeV]; D or B target, NaN when unavailable");
      trainHadronTable->addColumn<int>("observableInputsValid", auxInputsValid, "D-only observable-input mask, preserved for compatibility");
      trainHadronTable->addColumn<int>("auxiliaryTargetValid", auxPromptValid, "prompt D + unique charm QFrag + valid finite pre-hadronization momentum");
      trainHadronTable->addColumn<int>("bFragmentationTargetValid", bFragmentationTargetValid, "B + unique |pdgId|=5 QFrag + source-certain association + valid finite nonzero proxy kinematics");
    }
    const auto copyHadronFloat = [&](const std::string& name) {
      const auto values = column<float>(hadrons, name);
      return std::vector<float>(values.begin(), values.end());
    };
    const int generatorKinematicsColumn = hadrons.columnIndex("p");
    if (generatorKinematicsColumn >= 0) {
      trainHadronTable->addColumn<float>("p", copyHadronFloat("p"), "generator total momentum [GeV]");
      trainHadronTable->addColumn<float>("px", copyHadronFloat("px"), "generator px [GeV]");
      trainHadronTable->addColumn<float>("py", copyHadronFloat("py"), "generator py [GeV]");
      trainHadronTable->addColumn<float>("pz", copyHadronFloat("pz"), "generator pz [GeV]");
      trainHadronTable->addColumn<float>("energy", copyHadronFloat("energy"), "generator energy [GeV]");
      for (const std::string name : {"status", "isPrompt", "isHardProcess", "statusFlagLastCopy"}) {
        const auto values = column<int>(hadrons, name);
        trainHadronTable->addColumn<int>(name, std::vector<int>(values.begin(), values.end()), "selected hadron generator status metadata");
      }
    }
    const int prodXColumn = hadrons.columnIndex("prodX");
    const int decayXColumn = hadrons.columnIndex("decayX");
    if (generatorKinematicsColumn >= 0 && prodXColumn >= 0 && decayXColumn >= 0) {
      const auto copyIntColumn = [&](const std::string& name) {
        const auto values = column<int>(hadrons, name);
        return std::vector<int>(values.begin(), values.end());
      };
      trainHadronTable->addColumn<float>("prodX", copyHadronFloat("prodX"), "generator hadron production x [cm]");
      trainHadronTable->addColumn<float>("prodY", copyHadronFloat("prodY"), "generator hadron production y [cm]");
      trainHadronTable->addColumn<float>("prodZ", copyHadronFloat("prodZ"), "generator hadron production z [cm]");
      trainHadronTable->addColumn<int>("prodVertexValid", copyIntColumn("prodVertexValid"), "generator production-position validity");
      trainHadronTable->addColumn<float>("decayX", copyHadronFloat("decayX"), "generator decay x [cm]");
      trainHadronTable->addColumn<float>("decayY", copyHadronFloat("decayY"), "generator decay y [cm]");
      trainHadronTable->addColumn<float>("decayZ", copyHadronFloat("decayZ"), "generator decay z [cm]");
      trainHadronTable->addColumn<int>("decayVertexValid", copyIntColumn("decayVertexValid"), "generator decay-position validity");
      trainHadronTable->addColumn<int>("decayVertexResolved", copyIntColumn("decayVertexResolved"), "decay position uniquely resolved from immediate daughter vertices");
      trainHadronTable->addColumn<int>("decayVertexSource", copyIntColumn("decayVertexSource"), "1 immediate daughter production vertex");
      if (hadrons.columnIndex("trueLxy") >= 0) {
        trainHadronTable->addColumn<float>("trueLxy", copyHadronFloat("trueLxy"), "true transverse decay length [cm]");
        trainHadronTable->addColumn<float>("trueL3D", copyHadronFloat("trueL3D"), "true 3D decay length [cm]");
        trainHadronTable->addColumn<int>("trueDecayLengthValid", copyIntColumn("trueDecayLengthValid"), "true decay-length validity");
        trainHadronTable->addColumn<float>("trueCt", copyHadronFloat("trueCt"), "projected true proper decay length [cm]");
        trainHadronTable->addColumn<int>("trueCtValid", copyIntColumn("trueCtValid"), "true proper decay-length validity");
      }
    }
    event.put(std::move(trainHadronTable), "QTrainHadron");

    std::vector<int> trainQuarkPdg, trainQuarkOrigin, trainQuarkOriginValid;
    std::vector<float> trainQuarkPx, trainQuarkPy, trainQuarkPz, trainQuarkEnergy;
    for (size_t iq = 0; iq < nQuarks; ++iq) {
      trainQuarkPdg.push_back(qPdgId[iq]);
      trainQuarkPx.push_back(qPx[iq]); trainQuarkPy.push_back(qPy[iq]); trainQuarkPz.push_back(qPz[iq]); trainQuarkEnergy.push_back(qEnergy[iq]);
      trainQuarkOrigin.push_back(qOriginCategory[iq]); trainQuarkOriginValid.push_back(qOriginValid[iq]);
    }
    auto trainQuarkTable = std::make_unique<nanoaod::FlatTable>(nQuarks, "QTrainQuark", false, false);
    trainQuarkTable->addColumn<int>("pdgId", trainQuarkPdg, "signed pre-hadronization quark PDG identifier");
    trainQuarkTable->addColumn<float>("px", trainQuarkPx, "pre-hadronization quark px [GeV]");
    trainQuarkTable->addColumn<float>("py", trainQuarkPy, "pre-hadronization quark py [GeV]");
    trainQuarkTable->addColumn<float>("pz", trainQuarkPz, "pre-hadronization quark pz [GeV]");
    trainQuarkTable->addColumn<float>("energy", trainQuarkEnergy, "pre-hadronization quark energy [GeV]");
    trainQuarkTable->addColumn<int>("originCategory", trainQuarkOrigin, "0 unavailable/ambiguous hard-copy ancestry, 1 unique hard-process copy retained");
    trainQuarkTable->addColumn<int>("originCategoryValid", trainQuarkOriginValid, "origin-category validity mask");
    event.put(std::move(trainQuarkTable), "QTrainQuark");

    std::vector<int> legacyHadIdx, legacyPdg, legacyB, legacyD, legacyHF;
    for (size_t i = 0; i < legacyHadrons.nRows(); ++i) {
      legacyHadIdx.push_back(static_cast<int>(i)); legacyPdg.push_back(legacyPdgId[i]);
      legacyB.push_back(legacyIsB[i]); legacyD.push_back(legacyIsD[i]); legacyHF.push_back(legacyFromHF[i]);
    }
    auto trainLegacyHadronTable = std::make_unique<nanoaod::FlatTable>(legacyHadrons.nRows(), "QTrainLegacyHadron", false, false);
    trainLegacyHadronTable->addColumn<int>("hadronIdx", legacyHadIdx, "row in legacy GVTable; indexed by track.hadidx");
    trainLegacyHadronTable->addColumn<int>("pdgId", legacyPdg, "legacy GV hadron PDG identifier");
    trainLegacyHadronTable->addColumn<int>("isB", legacyB, "legacy GV B-hadron flag used by track.label");
    trainLegacyHadronTable->addColumn<int>("isD", legacyD, "legacy GV D-hadron flag");
    trainLegacyHadronTable->addColumn<int>("fromHF", legacyHF, "legacy GV heavy-flavor ancestry flag used by track.label");
    event.put(std::move(trainLegacyHadronTable), "QTrainLegacyHadron");

    std::vector<int> weakQuarkTruth, sameQuarkWeak;
    weakQuarkTruth.reserve(nPairs); sameQuarkWeak.reserve(nPairs);
    for (size_t i = 0; i < nPairs; ++i) {
      const int ti = edgeSrc[i], tj = edgeDst[i];
      const int gi = rGenIdx[ti], gj = rGenIdx[tj];
      const bool vi = gi >= 0 && static_cast<size_t>(gi) < nGen;
      const bool vj = gj >= 0 && static_cast<size_t>(gj) < nGen;
      const int qi = vi ? gQuarkIdx[gi] : -1, qj = vj ? gQuarkIdx[gj] : -1;
      const bool validQi = qi >= 0 && static_cast<size_t>(qi) < nQuarks;
      const bool validQj = qj >= 0 && static_cast<size_t>(qj) < nQuarks;
      const bool weak = vi && vj && validQi && validQj && gQuality[gi] == 1 && gQuality[gj] == 1 &&
                        (gMembership[gi] == 1 || gMembership[gi] == 2) &&
                        (gMembership[gj] == 1 || gMembership[gj] == 2);
      weakQuarkTruth.push_back(weak ? 1 : 0);
      sameQuarkWeak.push_back(weak ? (qi == qj ? 1 : 0) : -1);
    }
    auto trainPairTable = std::make_unique<nanoaod::FlatTable>(nPairs, "QTrainPair", false, false);
    trainPairTable->addColumn<int>("hadronTruthValid", hadronTruth, "both endpoints have valid heavy-hadron decay assignments");
    trainPairTable->addColumn<int>("sameHadron", sameHadron, "1 same known QHadron, 0 different known QHadrons, -1 unknown");
    trainPairTable->addColumn<int>("sameB", sameB, "same-B-hadron pair target; -1 when masked");
    trainPairTable->addColumn<int>("samePromptD", samePromptD, "same-prompt-D pair target; -1 when masked");
    trainPairTable->addColumn<int>("sameDFromB", sameDFromB, "same-D-from-B pair target; -1 when masked");
    trainPairTable->addColumn<int>("quarkTruthValid", quarkTruth, "both endpoints have unique source-certain quark ancestry");
    trainPairTable->addColumn<int>("sameQuark", sameQuark, "strict same-quark target; -1 when source certainty is insufficient");
    trainPairTable->addColumn<int>("quarkWeakTruthValid", weakQuarkTruth, "both endpoints have unique ancestry; uncertain additional candidates are weak-only");
    trainPairTable->addColumn<int>("sameQuarkWeak", sameQuarkWeak, "weak same-quark candidate target; -1 when ancestry is unresolved");
    event.put(std::move(trainPairTable), "QTrainPair");
  }

private:
  edm::EDGetTokenT<nanoaod::FlatTable> recoToken_, chargedToken_, hadronToken_, quarkToken_, legacyHadronToken_, pairToken_;
  bool produceFullValidationTables_;
  bool saveFragmentationObservables_;
};

DEFINE_FWK_MODULE(QuarkAssociationGNNTableProducer);
