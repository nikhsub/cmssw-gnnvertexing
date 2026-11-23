#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace {
  constexpr float kInvalid = -9999.f;

  struct TrackMapping {
    int unpacked = -1;
    int filtered = -1;
    int invalidReference = 0;
  };

  template <typename Collection>
  TrackMapping mapToFilteredTrack(const edm::Handle<Collection>& source,
                                  size_t sourceIndex,
                                  const edm::Handle<edm::ValueMap<int>>& sourceToUnpacked,
                                  const edm::Handle<reco::TrackCollection>& unpackedTracks,
                                  const edm::Handle<edm::ValueMap<int>>& unpackedToFiltered,
                                  size_t nFilteredTracks) {
    TrackMapping result;
    if (!source.isValid() || sourceIndex >= source->size()) {
      result.invalidReference = 1;
      return result;
    }
    const edm::Ref<Collection> sourceRef(source, sourceIndex);
    if (!sourceToUnpacked.isValid() || !sourceToUnpacked->contains(sourceRef.id())) {
      result.invalidReference = 1;
      return result;
    }
    result.unpacked = (*sourceToUnpacked)[sourceRef];
    if (result.unpacked < 0) return result;  // absent from the unpacked track collection by selection
    if (!unpackedTracks.isValid() || static_cast<size_t>(result.unpacked) >= unpackedTracks->size()) {
      result.unpacked = -1;
      result.invalidReference = 1;
      return result;
    }
    const reco::TrackRef unpackedRef(unpackedTracks, static_cast<unsigned int>(result.unpacked));
    if (!unpackedToFiltered.isValid() || !unpackedToFiltered->contains(unpackedRef.id())) {
      result.invalidReference = 1;
      return result;
    }
    result.filtered = (*unpackedToFiltered)[unpackedRef];
    if (result.filtered >= 0 && static_cast<size_t>(result.filtered) >= nFilteredTracks) {
      result.filtered = -1;
      result.invalidReference = 1;
    }
    return result;
  }

  float finiteOrInvalid(double value) {
    return std::isfinite(value) ? static_cast<float>(value) : kInvalid;
  }
}

class TrackGNNPFContextProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackGNNPFContextProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostToken_;
  edm::EDGetTokenT<reco::TrackCollection> unpackedTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> filteredTracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> packedTrackIndexToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> lostTrackIndexToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> filteredTrackIndexToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> embeddingToken_;
  double minChargedPt_;
  double minNeutralPt_;
  double maxPFAbsEta_;
};

TrackGNNPFContextProducer::TrackGNNPFContextProducer(const edm::ParameterSet& config)
    : packedToken_(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("packedCandidates"))),
      lostToken_(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("lostTracks"))),
      unpackedTracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("unpackedTracks"))),
      filteredTracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("filteredTracks"))),
      packedTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("packedTrackIndex"))),
      lostTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("lostTrackIndex"))),
      filteredTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("filteredTrackIndex"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("primaryVertices"))),
      embeddingToken_(mayConsume<nanoaod::FlatTable>(config.getParameter<edm::InputTag>("seederEmbeddingSrc"))),
      minChargedPt_(config.getParameter<double>("minChargedPt")),
      minNeutralPt_(config.getParameter<double>("minNeutralPt")),
      maxPFAbsEta_(config.getParameter<double>("maxPFAbsEta")) {
  produces<nanoaod::FlatTable>("CtxChargedPF");
  produces<nanoaod::FlatTable>("CtxNeutralPF");
  produces<nanoaod::FlatTable>("CtxPFSummary");
}

void TrackGNNPFContextProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("unpackedTracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("filteredTracks", edm::InputTag("filteredTracks"));
  desc.add<edm::InputTag>("packedTrackIndex", edm::InputTag("unpackedTracksAndVertices", "packedTrackIndex"));
  desc.add<edm::InputTag>("lostTrackIndex", edm::InputTag("unpackedTracksAndVertices", "lostTrackIndex"));
  desc.add<edm::InputTag>("filteredTrackIndex", edm::InputTag("filteredTracks", "filteredIndex"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("seederEmbeddingSrc", edm::InputTag("trackGNNSeederInference", "SeederEmbedding"));
  desc.add<double>("minChargedPt", 0.1);
  desc.add<double>("minNeutralPt", 0.1);
  desc.add<double>("maxPFAbsEta", 3.0);
  descriptions.add("trackGNNPFContextProducer", desc);
}

void TrackGNNPFContextProducer::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<pat::PackedCandidateCollection> packed;
  edm::Handle<pat::PackedCandidateCollection> lost;
  edm::Handle<reco::TrackCollection> unpackedTracks, filteredTracks;
  edm::Handle<edm::ValueMap<int>> packedTrackIndex, lostTrackIndex, filteredTrackIndex;
  edm::Handle<reco::VertexCollection> primaryVertices;
  edm::Handle<nanoaod::FlatTable> embeddings;
  event.getByToken(packedToken_, packed);
  event.getByToken(lostToken_, lost);
  event.getByToken(unpackedTracksToken_, unpackedTracks);
  event.getByToken(filteredTracksToken_, filteredTracks);
  event.getByToken(packedTrackIndexToken_, packedTrackIndex);
  event.getByToken(lostTrackIndexToken_, lostTrackIndex);
  event.getByToken(filteredTrackIndexToken_, filteredTrackIndex);
  event.getByToken(primaryVerticesToken_, primaryVertices);
  event.getByToken(embeddingToken_, embeddings);
  if (!packed.isValid() || !unpackedTracks.isValid() || !filteredTracks.isValid() ||
      !packedTrackIndex.isValid() || !filteredTrackIndex.isValid() || !primaryVertices.isValid()) {
    throw cms::Exception("ProductNotFound")
        << "TrackGNNPFContextProducer requires packed candidates, track mappings, and PVs";
  }
  if (embeddings.isValid() && embeddings->size() != filteredTracks->size()) {
    throw cms::Exception("TableSizeMismatch") << "SeederEmbedding has " << embeddings->size()
                                                << " rows for " << filteredTracks->size() << " filtered tracks";
  }
  const int* embeddingValid = nullptr;
  if (embeddings.isValid()) {
    const auto embeddingTrackIdx = embeddings->columnData<int>(embeddings->columnIndex("trackIdx"));
    embeddingValid = embeddings->columnData<int>(embeddings->columnIndex("valid"));
    for (size_t i = 0; i < embeddings->size(); ++i)
      if (embeddingTrackIdx[i] != static_cast<int>(i))
        throw cms::Exception("InvalidReference") << "SeederEmbedding trackIdx is not aligned with filteredTracks rows";
  }

  std::vector<int> chargedSourceId, chargedSourceIdx, chargedUnpackedIdx, chargedFilteredIdx;
  std::vector<int> chargedEmbeddingValid, chargedCharge, chargedPdgId, chargedHasTrackDetails;
  std::vector<int> chargedTrackHighPurity, chargedLostInnerHits, chargedNumberOfHits;
  std::vector<int> chargedNumberOfPixelHits, chargedNumberOfStripHits, chargedFromPV;
  std::vector<int> chargedPVAssociationQuality, chargedAssociatedPVIdx, chargedAssociatedPVValid;
  std::vector<int> chargedIPValid, chargedIPErrorValid, chargedIPReference;
  std::vector<float> chargedPt, chargedEta, chargedPhi, chargedMass, chargedEnergy;
  std::vector<float> chargedPx, chargedPy, chargedPz, chargedPuppiWeight, chargedPuppiWeightNoLep;
  std::vector<float> chargedDxy, chargedDz, chargedDxyError, chargedDzError;
  size_t invalidReferences = 0;
  size_t chargedMapped = 0;
  const auto appendCharged = [&](const pat::PackedCandidate& candidate,
                                 int sourceCollectionId,
                                 int sourceIndex,
                                 const TrackMapping& mapping,
                                 std::vector<int>& sourceIds,
                                 std::vector<int>& sourceIndices,
                                 std::vector<int>& unpackedIndices,
                                 std::vector<int>& filteredIndices,
                                 std::vector<int>& embeddingValids,
                                 std::vector<int>& charges,
                                 std::vector<int>& pdgIds,
                                 std::vector<int>& hasTrackDetails,
                                 std::vector<int>& highPurity,
                                 std::vector<int>& lostInnerHits,
                                 std::vector<int>& numberOfHits,
                                 std::vector<int>& numberOfPixelHits,
                                 std::vector<int>& numberOfStripHits,
                                 std::vector<int>& fromPV,
                                 std::vector<int>& pvAssociationQuality,
                                 std::vector<int>& associatedPVIdx,
                                 std::vector<int>& associatedPVValid,
                                 std::vector<int>& ipValid,
                                 std::vector<int>& ipErrorValid,
                                 std::vector<int>& ipReference,
                                 std::vector<float>& pt,
                                 std::vector<float>& eta,
                                 std::vector<float>& phi,
                                 std::vector<float>& mass,
                                 std::vector<float>& energy,
                                 std::vector<float>& px,
                                 std::vector<float>& py,
                                 std::vector<float>& pz,
                                 std::vector<float>& puppiWeight,
                                 std::vector<float>& puppiWeightNoLep,
                                 std::vector<float>& dxy,
                                 std::vector<float>& dz,
                                 std::vector<float>& dxyError,
                                 std::vector<float>& dzError) {
    sourceIds.push_back(sourceCollectionId);
    sourceIndices.push_back(sourceIndex);
    unpackedIndices.push_back(mapping.unpacked);
    filteredIndices.push_back(mapping.filtered);
    const bool embeddingOK = embeddingValid && mapping.filtered >= 0 &&
                              embeddingValid[static_cast<size_t>(mapping.filtered)] == 1;
    embeddingValids.push_back(embeddingOK ? 1 : 0);
    if (mapping.filtered >= 0) ++chargedMapped;
    invalidReferences += static_cast<size_t>(mapping.invalidReference);
    charges.push_back(candidate.charge());
    pdgIds.push_back(candidate.pdgId());
    const bool hasDetails = candidate.hasTrackDetails();
    hasTrackDetails.push_back(hasDetails ? 1 : 0);
    highPurity.push_back(hasDetails ? (candidate.trackHighPurity() ? 1 : 0) : -1);
    lostInnerHits.push_back(hasDetails ? static_cast<int>(candidate.lostInnerHits()) : -1);
    numberOfHits.push_back(hasDetails ? candidate.numberOfHits() : -1);
    numberOfPixelHits.push_back(hasDetails ? candidate.numberOfPixelHits() : -1);
    numberOfStripHits.push_back(hasDetails ? candidate.numberOfHits() - candidate.numberOfPixelHits() : -1);
    pvAssociationQuality.push_back(static_cast<int>(candidate.pvAssociationQuality()));
    int assocIdx = -1;
    bool assocValid = false;
    const reco::VertexRef vertexRef = candidate.vertexRef();
    if (vertexRef.isNonnull() && vertexRef.isAvailable() && vertexRef.id() == primaryVertices.id() &&
        vertexRef.key() < primaryVertices->size()) {
      assocIdx = static_cast<int>(vertexRef.key());
      assocValid = true;
    }
    associatedPVIdx.push_back(assocIdx);
    associatedPVValid.push_back(assocValid ? 1 : 0);
    fromPV.push_back(assocValid && !primaryVertices->empty() ? static_cast<int>(candidate.fromPV(0)) : -1);

    float dxyValue = kInvalid, dzValue = kInvalid;
    bool impactValid = false;
    int impactReference = -1;
    if (hasDetails && !primaryVertices->empty()) {
      dxyValue = candidate.dxy(primaryVertices->front().position());
      dzValue = candidate.dz(primaryVertices->front().position());
      impactValid = std::isfinite(dxyValue) && std::isfinite(dzValue);
      if (impactValid) impactReference = 0;  // both values are explicitly evaluated relative to MiniAOD PV0
    }
    const float dxyErr = candidate.dxyError();
    const float dzErr = candidate.dzError();
    const bool impactErrorValid = hasDetails && std::isfinite(dxyErr) && std::isfinite(dzErr) &&
                                  dxyErr >= 0.f && dzErr >= 0.f;
    ipValid.push_back(impactValid ? 1 : 0);
    ipErrorValid.push_back(impactErrorValid ? 1 : 0);
    ipReference.push_back(impactReference);
    dxy.push_back(impactValid ? dxyValue : kInvalid);
    dz.push_back(impactValid ? dzValue : kInvalid);
    dxyError.push_back(impactErrorValid ? dxyErr : kInvalid);
    dzError.push_back(impactErrorValid ? dzErr : kInvalid);
    pt.push_back(finiteOrInvalid(candidate.pt()));
    eta.push_back(finiteOrInvalid(candidate.eta()));
    phi.push_back(finiteOrInvalid(candidate.phi()));
    mass.push_back(finiteOrInvalid(candidate.mass()));
    energy.push_back(finiteOrInvalid(candidate.energy()));
    px.push_back(finiteOrInvalid(candidate.px()));
    py.push_back(finiteOrInvalid(candidate.py()));
    pz.push_back(finiteOrInvalid(candidate.pz()));
    puppiWeight.push_back(finiteOrInvalid(candidate.puppiWeight()));
    puppiWeightNoLep.push_back(finiteOrInvalid(candidate.puppiWeightNoLep()));
  };

  for (size_t i = 0; i < packed->size(); ++i) {
    const auto& candidate = (*packed)[i];
    if (candidate.charge() == 0 || !std::isfinite(candidate.pt()) || !std::isfinite(candidate.eta()) ||
        candidate.pt() < minChargedPt_ || std::abs(candidate.eta()) > maxPFAbsEta_) continue;
    const auto mapping = mapToFilteredTrack(packed, i, packedTrackIndex, unpackedTracks,
                                            filteredTrackIndex, filteredTracks->size());
    appendCharged(candidate, 0, static_cast<int>(i), mapping,
                  chargedSourceId, chargedSourceIdx, chargedUnpackedIdx, chargedFilteredIdx,
                  chargedEmbeddingValid, chargedCharge, chargedPdgId, chargedHasTrackDetails,
                  chargedTrackHighPurity, chargedLostInnerHits, chargedNumberOfHits,
                  chargedNumberOfPixelHits, chargedNumberOfStripHits, chargedFromPV,
                  chargedPVAssociationQuality, chargedAssociatedPVIdx, chargedAssociatedPVValid,
                  chargedIPValid, chargedIPErrorValid, chargedIPReference, chargedPt, chargedEta,
                  chargedPhi, chargedMass, chargedEnergy, chargedPx, chargedPy, chargedPz,
                  chargedPuppiWeight, chargedPuppiWeightNoLep, chargedDxy, chargedDz,
                  chargedDxyError, chargedDzError);
  }
  if (lost.isValid()) {
    for (size_t i = 0; i < lost->size(); ++i) {
      const auto& candidate = (*lost)[i];
      if (candidate.charge() == 0 || !std::isfinite(candidate.pt()) || !std::isfinite(candidate.eta()) ||
          candidate.pt() < minChargedPt_ || std::abs(candidate.eta()) > maxPFAbsEta_) continue;
      const auto mapping = mapToFilteredTrack(lost, i, lostTrackIndex, unpackedTracks,
                                              filteredTrackIndex, filteredTracks->size());
      appendCharged(candidate, 1, static_cast<int>(i), mapping,
                    chargedSourceId, chargedSourceIdx, chargedUnpackedIdx, chargedFilteredIdx,
                    chargedEmbeddingValid, chargedCharge, chargedPdgId, chargedHasTrackDetails,
                    chargedTrackHighPurity, chargedLostInnerHits, chargedNumberOfHits,
                    chargedNumberOfPixelHits, chargedNumberOfStripHits, chargedFromPV,
                    chargedPVAssociationQuality, chargedAssociatedPVIdx, chargedAssociatedPVValid,
                    chargedIPValid, chargedIPErrorValid, chargedIPReference, chargedPt, chargedEta,
                    chargedPhi, chargedMass, chargedEnergy, chargedPx, chargedPy, chargedPz,
                    chargedPuppiWeight, chargedPuppiWeightNoLep, chargedDxy, chargedDz,
                    chargedDxyError, chargedDzError);
    }
  }

  auto chargedTable = std::make_unique<nanoaod::FlatTable>(chargedPt.size(), "CtxChargedPF", false);
#define ADD_CHARGED(TYPE, NAME, DATA, DOC) chargedTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_CHARGED(int, "sourceCollectionId", chargedSourceId, "0=packedPFCandidates, 1=lostTracks");
  ADD_CHARGED(int, "sourceIdx", chargedSourceIdx, "index in the original source collection");
  ADD_CHARGED(int, "unpackedTrackIdx", chargedUnpackedIdx, "exact ValueMap mapping to unpackedTracksAndVertices; -1 if not unpacked");
  ADD_CHARGED(int, "filteredTrackIdx", chargedFilteredIdx, "exact ValueMap mapping to filteredTracks; -1 if absent from GNN graph");
  ADD_CHARGED(int, "embeddingValid", chargedEmbeddingValid, "SeederEmbedding_valid for filteredTrackIdx; 0 when unmapped/invalid");
  ADD_CHARGED(float, "pt", chargedPt, "raw candidate transverse momentum [GeV]");
  ADD_CHARGED(float, "eta", chargedEta, "raw candidate pseudorapidity");
  ADD_CHARGED(float, "phi", chargedPhi, "raw candidate azimuth [rad]");
  ADD_CHARGED(float, "mass", chargedMass, "raw candidate mass [GeV]");
  ADD_CHARGED(float, "energy", chargedEnergy, "raw candidate energy [GeV]");
  ADD_CHARGED(float, "px", chargedPx, "raw candidate px [GeV]");
  ADD_CHARGED(float, "py", chargedPy, "raw candidate py [GeV]");
  ADD_CHARGED(float, "pz", chargedPz, "raw candidate pz [GeV]");
  ADD_CHARGED(int, "charge", chargedCharge, "candidate electric charge");
  ADD_CHARGED(int, "pdgId", chargedPdgId, "candidate PDG identifier");
  ADD_CHARGED(float, "puppiWeight", chargedPuppiWeight, "raw PUPPI weight");
  ADD_CHARGED(float, "puppiWeightNoLep", chargedPuppiWeightNoLep, "raw PUPPI weight with leptons removed");
  ADD_CHARGED(int, "hasTrackDetails", chargedHasTrackDetails, "PackedCandidate has stored track details");
  ADD_CHARGED(int, "trackHighPurity", chargedTrackHighPurity, "packed track high-purity bit; -1 if unavailable");
  ADD_CHARGED(int, "lostInnerHits", chargedLostInnerHits, "PackedCandidate lostInnerHits enum; -1 if unavailable");
  ADD_CHARGED(int, "numberOfHits", chargedNumberOfHits, "packed valid hit count; -1 if unavailable");
  ADD_CHARGED(int, "numberOfPixelHits", chargedNumberOfPixelHits, "packed valid pixel-hit count; -1 if unavailable");
  ADD_CHARGED(int, "numberOfStripHits", chargedNumberOfStripHits, "numberOfHits-numberOfPixelHits; -1 if unavailable");
  ADD_CHARGED(int, "fromPV", chargedFromPV, "PackedCandidate fromPV classification relative to PV0; -1 if unavailable");
  ADD_CHARGED(int, "pvAssociationQuality", chargedPVAssociationQuality, "PackedCandidate PVAssociationQuality enum");
  ADD_CHARGED(int, "associatedPVIdx", chargedAssociatedPVIdx, "original MiniAOD PV index from vertexRef; -1 if unavailable");
  ADD_CHARGED(int, "associatedPVValid", chargedAssociatedPVValid, "1 if vertexRef points into offlineSlimmedPrimaryVertices");
  ADD_CHARGED(float, "dxy", chargedDxy, "transverse impact parameter explicitly relative to MiniAOD PV0 [cm]");
  ADD_CHARGED(float, "dz", chargedDz, "longitudinal impact parameter explicitly relative to MiniAOD PV0 [cm]");
  ADD_CHARGED(float, "dxyError", chargedDxyError, "packed transverse impact parameter uncertainty [cm]");
  ADD_CHARGED(float, "dzError", chargedDzError, "packed longitudinal impact parameter uncertainty [cm]");
  ADD_CHARGED(int, "ipValid", chargedIPValid, "1 if PV0-relative dxy and dz are finite");
  ADD_CHARGED(int, "ipErrorValid", chargedIPErrorValid, "1 if packed dxyError and dzError are finite/nonnegative");
  ADD_CHARGED(int, "ipReference", chargedIPReference, "0=offlineSlimmedPrimaryVertices PV0, -1=unavailable");
#undef ADD_CHARGED
  event.put(std::move(chargedTable), "CtxChargedPF");

  std::vector<int> neutralSourceIdx, neutralPdgId, neutralCharge, neutralCaloFractionValid;
  std::vector<float> neutralPt, neutralEta, neutralPhi, neutralMass, neutralEnergy;
  std::vector<float> neutralPx, neutralPy, neutralPz, neutralPuppiWeight, neutralPuppiWeightNoLep;
  std::vector<float> neutralCaloFraction, neutralHcalFraction, neutralRawCaloFraction, neutralRawHcalFraction;
  for (size_t i = 0; i < packed->size(); ++i) {
    const auto& candidate = (*packed)[i];
    if (candidate.charge() != 0 || !std::isfinite(candidate.pt()) || !std::isfinite(candidate.eta()) ||
        candidate.pt() < minNeutralPt_ || std::abs(candidate.eta()) > maxPFAbsEta_) continue;
    neutralSourceIdx.push_back(static_cast<int>(i));
    neutralPt.push_back(finiteOrInvalid(candidate.pt())); neutralEta.push_back(finiteOrInvalid(candidate.eta()));
    neutralPhi.push_back(finiteOrInvalid(candidate.phi())); neutralMass.push_back(finiteOrInvalid(candidate.mass()));
    neutralEnergy.push_back(finiteOrInvalid(candidate.energy())); neutralPx.push_back(finiteOrInvalid(candidate.px()));
    neutralPy.push_back(finiteOrInvalid(candidate.py())); neutralPz.push_back(finiteOrInvalid(candidate.pz()));
    neutralPdgId.push_back(candidate.pdgId()); neutralCharge.push_back(candidate.charge());
    neutralPuppiWeight.push_back(finiteOrInvalid(candidate.puppiWeight()));
    neutralPuppiWeightNoLep.push_back(finiteOrInvalid(candidate.puppiWeightNoLep()));
    const bool fractionsValid = std::isfinite(candidate.caloFraction()) && std::isfinite(candidate.hcalFraction()) &&
                                std::isfinite(candidate.rawCaloFraction()) && std::isfinite(candidate.rawHcalFraction()) &&
                                candidate.energy() > 0.;
    neutralCaloFractionValid.push_back(fractionsValid ? 1 : 0);
    neutralCaloFraction.push_back(fractionsValid ? candidate.caloFraction() : kInvalid);
    neutralHcalFraction.push_back(fractionsValid ? candidate.hcalFraction() : kInvalid);
    neutralRawCaloFraction.push_back(fractionsValid ? candidate.rawCaloFraction() : kInvalid);
    neutralRawHcalFraction.push_back(fractionsValid ? candidate.rawHcalFraction() : kInvalid);
  }
  auto neutralTable = std::make_unique<nanoaod::FlatTable>(neutralPt.size(), "CtxNeutralPF", false);
#define ADD_NEUTRAL(TYPE, NAME, DATA, DOC) neutralTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_NEUTRAL(int, "sourceIdx", neutralSourceIdx, "index in packedPFCandidates");
  ADD_NEUTRAL(float, "pt", neutralPt, "raw candidate transverse momentum [GeV]");
  ADD_NEUTRAL(float, "eta", neutralEta, "raw candidate pseudorapidity");
  ADD_NEUTRAL(float, "phi", neutralPhi, "raw candidate azimuth [rad]");
  ADD_NEUTRAL(float, "mass", neutralMass, "raw candidate mass [GeV]");
  ADD_NEUTRAL(float, "energy", neutralEnergy, "raw candidate energy [GeV]");
  ADD_NEUTRAL(float, "px", neutralPx, "raw candidate px [GeV]");
  ADD_NEUTRAL(float, "py", neutralPy, "raw candidate py [GeV]");
  ADD_NEUTRAL(float, "pz", neutralPz, "raw candidate pz [GeV]");
  ADD_NEUTRAL(int, "pdgId", neutralPdgId, "candidate PDG identifier (e.g. 22 photon, 130 neutral hadron)");
  ADD_NEUTRAL(int, "charge", neutralCharge, "candidate electric charge; selected rows are neutral");
  ADD_NEUTRAL(float, "puppiWeight", neutralPuppiWeight, "raw PUPPI weight");
  ADD_NEUTRAL(float, "puppiWeightNoLep", neutralPuppiWeightNoLep, "raw PUPPI weight with leptons removed");
  ADD_NEUTRAL(float, "caloFraction", neutralCaloFraction, "PackedCandidate ECAL+HCAL energy fraction; -9999 if invalid");
  ADD_NEUTRAL(float, "hcalFraction", neutralHcalFraction, "PackedCandidate HCAL energy fraction; -9999 if invalid");
  ADD_NEUTRAL(float, "rawCaloFraction", neutralRawCaloFraction, "raw PackedCandidate ECAL+HCAL fraction; -9999 if invalid");
  ADD_NEUTRAL(float, "rawHcalFraction", neutralRawHcalFraction, "raw PackedCandidate HCAL fraction; -9999 if invalid");
  ADD_NEUTRAL(int, "caloFractionsValid", neutralCaloFractionValid, "1 if packed calorimeter fractions are finite and candidate energy is positive");
#undef ADD_NEUTRAL
  event.put(std::move(neutralTable), "CtxNeutralPF");

  auto summary = std::make_unique<nanoaod::FlatTable>(1, "CtxPFSummary", true);
  const int nLostInput = lost.isValid() ? static_cast<int>(lost->size()) : 0;
  const int nMapped = static_cast<int>(chargedMapped);
  const int nCharged = static_cast<int>(chargedPt.size());
  summary->addColumn<int>("nPackedInput", std::vector<int>{static_cast<int>(packed->size())}, "input packedPFCandidates size before selection");
  summary->addColumn<int>("nLostInput", std::vector<int>{nLostInput}, "input lostTracks size before selection; zero if missing");
  summary->addColumn<int>("nUnpackedTracks", std::vector<int>{static_cast<int>(unpackedTracks->size())}, "unpacked track collection size");
  summary->addColumn<int>("nFilteredTracks", std::vector<int>{static_cast<int>(filteredTracks->size())}, "filtered GNN track collection size");
  summary->addColumn<int>("nChargedPF", std::vector<int>{nCharged}, "rows in CtxChargedPF");
  summary->addColumn<int>("nNeutralPF", std::vector<int>{static_cast<int>(neutralPt.size())}, "rows in CtxNeutralPF");
  summary->addColumn<int>("nChargedMappedToFilteredTracks", std::vector<int>{nMapped}, "charged PF rows with filteredTrackIdx >= 0");
  summary->addColumn<int>("nChargedUnmapped", std::vector<int>{nCharged - nMapped}, "charged PF rows with filteredTrackIdx < 0");
  summary->addColumn<int>("nInvalidReferences", std::vector<int>{static_cast<int>(invalidReferences)}, "invalid/missing ValueMap references or out-of-range indices; expected -1 selection exclusions are not errors");
  event.put(std::move(summary), "CtxPFSummary");
}

DEFINE_FWK_MODULE(TrackGNNPFContextProducer);
