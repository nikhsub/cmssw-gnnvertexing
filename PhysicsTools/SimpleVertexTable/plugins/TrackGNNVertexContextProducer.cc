#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace {
  constexpr float kInvalid = -9999.f;

  bool finite3(double x, double y, double z) {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
  }

  bool covarianceValid(const reco::Vertex::CovarianceMatrix& cov) {
    const double xx = cov(0, 0), xy = cov(0, 1), xz = cov(0, 2);
    const double yy = cov(1, 1), yz = cov(1, 2), zz = cov(2, 2);
    if (!std::isfinite(xx) || !std::isfinite(xy) || !std::isfinite(xz) ||
        !std::isfinite(yy) || !std::isfinite(yz) || !std::isfinite(zz)) return false;
    const double scale = std::max({std::abs(xx), std::abs(yy), std::abs(zz), 1.e-30});
    const double tol2 = 1.e-10 * scale * scale;
    const double tol3 = tol2 * scale;
    const double minorXY = xx * yy - xy * xy;
    const double minorXZ = xx * zz - xz * xz;
    const double minorYZ = yy * zz - yz * yz;
    const double determinant = xx * yy * zz + 2. * xy * xz * yz - xx * yz * yz - yy * xz * xz - zz * xy * xy;
    return xx >= 0. && yy >= 0. && zz >= 0. && minorXY >= -tol2 && minorXZ >= -tol2 &&
           minorYZ >= -tol2 && determinant >= -tol3;
  }

  float matrixValue(const reco::Vertex::CovarianceMatrix& cov, int i, int j) {
    const double value = cov(i, j);
    return std::isfinite(value) ? static_cast<float>(value) : kInvalid;
  }

  double clampCosine(double value) {
    return std::max(-1., std::min(1., value));
  }

  struct SourceMapping {
    int collectionId = -1;
    int sourceIndex = -1;
    int unpackedIndex = -1;
    int filteredIndex = -1;
    int sourceIdentityValid = 0;
    int chargedPFIndex = -1;
    int chargedPFMappingValid = 0;
    int filteredTrackMappingValid = 0;
  };

  std::uint64_t sourceKey(int collectionId, int sourceIndex) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(collectionId)) << 32) |
           static_cast<std::uint32_t>(sourceIndex);
  }
}

class TrackGNNVertexContextProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackGNNVertexContextProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secondaryVerticesToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandidatesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> unpackedTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> filteredTracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> packedTrackIndexToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> lostTrackIndexToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> filteredTrackIndexToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> chargedPFToken_;
};

TrackGNNVertexContextProducer::TrackGNNVertexContextProducer(const edm::ParameterSet& config)
    : secondaryVerticesToken_(consumes<reco::VertexCompositePtrCandidateCollection>(config.getParameter<edm::InputTag>("secondaryVertices"))),
      primaryVerticesToken_(consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("primaryVertices"))),
      packedCandidatesToken_(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("packedCandidates"))),
      lostTracksToken_(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("lostTracks"))),
      unpackedTracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("unpackedTracks"))),
      filteredTracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("filteredTracks"))),
      packedTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("packedTrackIndex"))),
      lostTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("lostTrackIndex"))),
      filteredTrackIndexToken_(consumes<edm::ValueMap<int>>(config.getParameter<edm::InputTag>("filteredTrackIndex"))),
      chargedPFToken_(consumes<nanoaod::FlatTable>(config.getParameter<edm::InputTag>("chargedPF"))) {
  produces<nanoaod::FlatTable>("CtxPV");
  produces<nanoaod::FlatTable>("CtxSV");
  produces<nanoaod::FlatTable>("CtxSVConstituent");
}

void TrackGNNVertexContextProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("secondaryVertices", edm::InputTag("slimmedSecondaryVertices"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("packedCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("unpackedTracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("filteredTracks", edm::InputTag("filteredTracks"));
  desc.add<edm::InputTag>("packedTrackIndex", edm::InputTag("unpackedTracksAndVertices", "packedTrackIndex"));
  desc.add<edm::InputTag>("lostTrackIndex", edm::InputTag("unpackedTracksAndVertices", "lostTrackIndex"));
  desc.add<edm::InputTag>("filteredTrackIndex", edm::InputTag("filteredTracks", "filteredIndex"));
  desc.add<edm::InputTag>("chargedPF", edm::InputTag("trackGNNPFContext", "CtxChargedPF"));
  descriptions.add("trackGNNVertexContextProducer", desc);
}

void TrackGNNVertexContextProducer::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVertices;
  edm::Handle<reco::VertexCollection> primaryVertices;
  edm::Handle<pat::PackedCandidateCollection> packedCandidates, lostTracks;
  edm::Handle<reco::TrackCollection> unpackedTracks, filteredTracks;
  edm::Handle<edm::ValueMap<int>> packedTrackIndex, lostTrackIndex, filteredTrackIndex;
  edm::Handle<nanoaod::FlatTable> chargedPF;
  event.getByToken(secondaryVerticesToken_, secondaryVertices);
  event.getByToken(primaryVerticesToken_, primaryVertices);
  event.getByToken(packedCandidatesToken_, packedCandidates);
  event.getByToken(lostTracksToken_, lostTracks);
  event.getByToken(unpackedTracksToken_, unpackedTracks);
  event.getByToken(filteredTracksToken_, filteredTracks);
  event.getByToken(packedTrackIndexToken_, packedTrackIndex);
  event.getByToken(lostTrackIndexToken_, lostTrackIndex);
  event.getByToken(filteredTrackIndexToken_, filteredTrackIndex);
  event.getByToken(chargedPFToken_, chargedPF);
  if (!secondaryVertices.isValid() || !primaryVertices.isValid() || !packedCandidates.isValid() ||
      !unpackedTracks.isValid() || !filteredTracks.isValid() || !packedTrackIndex.isValid() ||
      !filteredTrackIndex.isValid() || !chargedPF.isValid()) {
    throw cms::Exception("ProductNotFound") << "TrackGNNVertexContextProducer requires source PV/SV collections, track maps, and CtxChargedPF";
  }
  auto chargedPFSourceId = chargedPF->columnData<int>(chargedPF->columnIndex("sourceCollectionId"));
  auto chargedPFSourceIdx = chargedPF->columnData<int>(chargedPF->columnIndex("sourceIdx"));
  std::map<std::uint64_t, int> chargedPFIndexBySource;
  for (size_t i = 0; i < chargedPF->size(); ++i) {
    const std::uint64_t key = sourceKey(chargedPFSourceId[i], chargedPFSourceIdx[i]);
    if (!chargedPFIndexBySource.emplace(key, static_cast<int>(i)).second)
      throw cms::Exception("DuplicatePFIdentity") << "CtxChargedPF contains duplicate source identity "
                                                   << chargedPFSourceId[i] << ":" << chargedPFSourceIdx[i];
  }

  std::vector<int> pvSourceIdx, pvIsFake, pvNTracks, pvPositionValid, pvCovarianceValid, pvChi2NdofValid;
  std::vector<float> pvX, pvY, pvZ, pvChi2, pvNdof, pvChi2Ndof;
  std::vector<float> pvCovXX, pvCovXY, pvCovXZ, pvCovYY, pvCovYZ, pvCovZZ;
  for (size_t i = 0; i < primaryVertices->size(); ++i) {
    const auto& pv = (*primaryVertices)[i];
    const auto& cov = pv.covariance();
    const bool positionOK = finite3(pv.x(), pv.y(), pv.z());
    const bool covOK = covarianceValid(cov);
    const bool chi2OK = std::isfinite(pv.chi2());
    const bool ndofOK = std::isfinite(pv.ndof());
    const bool chi2NdofOK = chi2OK && ndofOK && pv.ndof() > 0.;
    pvSourceIdx.push_back(static_cast<int>(i)); pvIsFake.push_back(pv.isFake() ? 1 : 0);
    pvNTracks.push_back(static_cast<int>(pv.tracksSize()));
    pvPositionValid.push_back(positionOK ? 1 : 0); pvCovarianceValid.push_back(covOK ? 1 : 0);
    pvChi2NdofValid.push_back(chi2NdofOK ? 1 : 0);
    pvX.push_back(positionOK ? pv.x() : kInvalid); pvY.push_back(positionOK ? pv.y() : kInvalid); pvZ.push_back(positionOK ? pv.z() : kInvalid);
    pvChi2.push_back(chi2OK ? pv.chi2() : kInvalid); pvNdof.push_back(ndofOK ? pv.ndof() : kInvalid);
    pvChi2Ndof.push_back(chi2NdofOK ? pv.chi2() / pv.ndof() : kInvalid);
    pvCovXX.push_back(matrixValue(cov, 0, 0)); pvCovXY.push_back(matrixValue(cov, 0, 1)); pvCovXZ.push_back(matrixValue(cov, 0, 2));
    pvCovYY.push_back(matrixValue(cov, 1, 1)); pvCovYZ.push_back(matrixValue(cov, 1, 2)); pvCovZZ.push_back(matrixValue(cov, 2, 2));
  }
  auto pvTable = std::make_unique<nanoaod::FlatTable>(primaryVertices->size(), "CtxPV", false);
#define ADD_PV(TYPE, NAME, DATA, DOC) pvTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_PV(int, "sourceIdx", pvSourceIdx, "original row in offlineSlimmedPrimaryVertices; PV0 is sourceIdx 0");
  ADD_PV(float, "x", pvX, "PV x position [cm]"); ADD_PV(float, "y", pvY, "PV y position [cm]"); ADD_PV(float, "z", pvZ, "PV z position [cm]");
  ADD_PV(float, "chi2", pvChi2, "PV fit chi2"); ADD_PV(float, "ndof", pvNdof, "PV fit degrees of freedom");
  ADD_PV(float, "chi2Ndof", pvChi2Ndof, "PV chi2/ndof; -9999 if ndof is nonpositive");
  ADD_PV(int, "isFake", pvIsFake, "reco::Vertex isFake flag"); ADD_PV(int, "nTracks", pvNTracks, "number of referenced PV tracks");
  ADD_PV(float, "covXX", pvCovXX, "PV covariance xx [cm^2]"); ADD_PV(float, "covXY", pvCovXY, "PV covariance xy [cm^2]"); ADD_PV(float, "covXZ", pvCovXZ, "PV covariance xz [cm^2]");
  ADD_PV(float, "covYY", pvCovYY, "PV covariance yy [cm^2]"); ADD_PV(float, "covYZ", pvCovYZ, "PV covariance yz [cm^2]"); ADD_PV(float, "covZZ", pvCovZZ, "PV covariance zz [cm^2]");
  ADD_PV(int, "positionValid", pvPositionValid, "1 if PV coordinates are finite");
  ADD_PV(int, "covarianceValid", pvCovarianceValid, "1 if PV covariance is finite and positive semidefinite");
  ADD_PV(int, "chi2NdofValid", pvChi2NdofValid, "1 if chi2/ndof can be computed");
#undef ADD_PV
  event.put(std::move(pvTable), "CtxPV");

  std::vector<int> svSourceIdx, svNDaughters, svNMappedTracks, svCovarianceValid, svChi2NdofValid;
  std::vector<int> svPVGeometryValid, svLxyValid, svLxyErrValid, svLxySigValid;
  std::vector<int> svL3DValid, svL3DErrValid, svL3DSigValid, svCosPointing2DValid, svCosPointing3DValid;
  std::vector<float> svPt, svEta, svPhi, svMass, svEnergy, svX, svY, svZ, svChi2, svNdof, svChi2Ndof;
  std::vector<float> svCovXX, svCovXY, svCovXZ, svCovYY, svCovYZ, svCovZZ;
  std::vector<float> svLxy, svLxyErr, svLxySig, svL3D, svL3DErr, svL3DSig, svCosPointing2D, svCosPointing3D;
  std::vector<int> constituentSVIdx, constituentSourceSVIdx, constituentDaughterIdx;
  std::vector<int> constituentSourceCollectionId, constituentSourceIdx, constituentUnpackedTrackIdx;
  std::vector<int> constituentChargedPFIdx, constituentFilteredTrackIdx;
  std::vector<int> constituentMappingValid, constituentChargedPFMappingValid, constituentFilteredTrackMappingValid;
  std::vector<float> constituentPt, constituentEta, constituentPhi;
  std::vector<int> constituentCharge;
  VertexDistanceXY distanceXY;
  VertexDistance3D distance3D;
  const bool pv0Usable = !primaryVertices->empty() && !primaryVertices->front().isFake() &&
                         finite3(primaryVertices->front().x(), primaryVertices->front().y(), primaryVertices->front().z()) &&
                         covarianceValid(primaryVertices->front().covariance());
  for (size_t isv = 0; isv < secondaryVertices->size(); ++isv) {
    const auto& sv = (*secondaryVertices)[isv];
    const auto cov = sv.vertexCovariance();
    const bool svCovOK = covarianceValid(cov);
    const bool positionOK = finite3(sv.vx(), sv.vy(), sv.vz());
    const bool chi2OK = std::isfinite(sv.vertexChi2());
    const bool ndofOK = std::isfinite(sv.vertexNdof());
    const bool chi2NdofOK = chi2OK && ndofOK && sv.vertexNdof() > 0.;
    int nMappedTracks = 0;
    svSourceIdx.push_back(static_cast<int>(isv));
    svNDaughters.push_back(static_cast<int>(sv.numberOfDaughters()));
    svPt.push_back(std::isfinite(sv.pt()) ? sv.pt() : kInvalid);
    svEta.push_back(std::isfinite(sv.eta()) ? sv.eta() : kInvalid);
    svPhi.push_back(std::isfinite(sv.phi()) ? sv.phi() : kInvalid);
    svMass.push_back(std::isfinite(sv.mass()) ? sv.mass() : kInvalid);
    svEnergy.push_back(std::isfinite(sv.energy()) ? sv.energy() : kInvalid);
    svX.push_back(positionOK ? sv.vx() : kInvalid); svY.push_back(positionOK ? sv.vy() : kInvalid); svZ.push_back(positionOK ? sv.vz() : kInvalid);
    svChi2.push_back(chi2OK ? sv.vertexChi2() : kInvalid); svNdof.push_back(ndofOK ? sv.vertexNdof() : kInvalid);
    svChi2Ndof.push_back(chi2NdofOK ? sv.vertexChi2() / sv.vertexNdof() : kInvalid);
    svCovarianceValid.push_back(svCovOK ? 1 : 0); svChi2NdofValid.push_back(chi2NdofOK ? 1 : 0);
    svCovXX.push_back(matrixValue(cov, 0, 0)); svCovXY.push_back(matrixValue(cov, 0, 1)); svCovXZ.push_back(matrixValue(cov, 0, 2));
    svCovYY.push_back(matrixValue(cov, 1, 1)); svCovYZ.push_back(matrixValue(cov, 1, 2)); svCovZZ.push_back(matrixValue(cov, 2, 2));

    float lxyValue = kInvalid, lxyError = kInvalid, lxySignificance = kInvalid;
    float l3dValue = kInvalid, l3dError = kInvalid, l3dSignificance = kInvalid;
    float cos2d = kInvalid, cos3d = kInvalid;
    int lxyOK = 0, lxyErrOK = 0, lxySigOK = 0, l3dOK = 0, l3dErrOK = 0, l3dSigOK = 0, cos2dOK = 0, cos3dOK = 0;
    const bool geometryInputsValid = pv0Usable && svCovOK && positionOK;
    if (geometryInputsValid) {
      reco::Vertex recoSV(sv.vertex(), cov, sv.vertexChi2(), sv.vertexNdof(), 0);
      const auto measurementXY = distanceXY.distance(primaryVertices->front(), recoSV);
      const auto measurement3D = distance3D.distance(primaryVertices->front(), recoSV);
      if (std::isfinite(measurementXY.value())) {
        lxyValue = measurementXY.value(); lxyOK = 1;
      }
      if (std::isfinite(measurementXY.error()) && measurementXY.error() > 0.) {
        lxyError = measurementXY.error(); lxyErrOK = 1;
      }
      if (lxyOK && lxyErrOK && std::isfinite(measurementXY.significance())) {
        lxySignificance = measurementXY.significance(); lxySigOK = 1;
      }
      if (std::isfinite(measurement3D.value())) {
        l3dValue = measurement3D.value(); l3dOK = 1;
      }
      if (std::isfinite(measurement3D.error()) && measurement3D.error() > 0.) {
        l3dError = measurement3D.error(); l3dErrOK = 1;
      }
      if (l3dOK && l3dErrOK && std::isfinite(measurement3D.significance())) {
        l3dSignificance = measurement3D.significance(); l3dSigOK = 1;
      }
      const double dx = sv.vx() - primaryVertices->front().x();
      const double dy = sv.vy() - primaryVertices->front().y();
      const double dz = sv.vz() - primaryVertices->front().z();
      const double px = sv.px(), py = sv.py(), pz = sv.pz();
      const double pt = std::hypot(px, py);
      const double p = std::sqrt(px * px + py * py + pz * pz);
      if (lxyOK && lxyValue > 0.f && pt > 0. && std::isfinite(pt)) {
        const double value = clampCosine((dx * px + dy * py) / (lxyValue * pt));
        if (std::isfinite(value)) { cos2d = value; cos2dOK = 1; }
      }
      if (l3dOK && l3dValue > 0.f && p > 0. && std::isfinite(p)) {
        const double value = clampCosine((dx * px + dy * py + dz * pz) / (l3dValue * p));
        if (std::isfinite(value)) { cos3d = value; cos3dOK = 1; }
      }
    }
    svPVGeometryValid.push_back(geometryInputsValid ? 1 : 0);
    svLxyValid.push_back(lxyOK); svLxyErrValid.push_back(lxyErrOK); svLxySigValid.push_back(lxySigOK);
    svL3DValid.push_back(l3dOK); svL3DErrValid.push_back(l3dErrOK); svL3DSigValid.push_back(l3dSigOK);
    svCosPointing2DValid.push_back(cos2dOK); svCosPointing3DValid.push_back(cos3dOK);
    svLxy.push_back(lxyValue); svLxyErr.push_back(lxyError); svLxySig.push_back(lxySignificance);
    svL3D.push_back(l3dValue); svL3DErr.push_back(l3dError); svL3DSig.push_back(l3dSignificance);
    svCosPointing2D.push_back(cos2d); svCosPointing3D.push_back(cos3d);

    for (size_t daughterIndex = 0; daughterIndex < sv.numberOfDaughters(); ++daughterIndex) {
      const auto daughterPtr = sv.daughterPtr(daughterIndex);
      const reco::Candidate* daughter = daughterPtr.isNonnull() && daughterPtr.isAvailable() ? daughterPtr.get() : nullptr;
      SourceMapping mapping;
      if (daughterPtr.isNonnull()) {
        mapping.sourceIndex = static_cast<int>(daughterPtr.key());
        if (daughterPtr.id() == packedCandidates.id()) mapping.collectionId = 0;
        else if (lostTracks.isValid() && daughterPtr.id() == lostTracks.id()) mapping.collectionId = 1;
        if (mapping.collectionId >= 0 && mapping.sourceIndex >= 0) mapping.sourceIdentityValid = 1;
      }
      if (mapping.collectionId == 0 && packedTrackIndex->contains(daughterPtr.id())) {
        mapping.unpackedIndex = (*packedTrackIndex)[daughterPtr];
      } else if (mapping.collectionId == 1 && lostTrackIndex.isValid() && lostTrackIndex->contains(daughterPtr.id())) {
        mapping.unpackedIndex = (*lostTrackIndex)[daughterPtr];
      }
      if (mapping.unpackedIndex >= 0 && static_cast<size_t>(mapping.unpackedIndex) < unpackedTracks->size()) {
        const reco::TrackRef unpackedRef(unpackedTracks, static_cast<unsigned int>(mapping.unpackedIndex));
        if (filteredTrackIndex->contains(unpackedRef.id())) {
          mapping.filteredIndex = (*filteredTrackIndex)[unpackedRef];
          if (mapping.filteredIndex >= 0 && static_cast<size_t>(mapping.filteredIndex) >= filteredTracks->size())
            mapping.filteredIndex = -1;
        }
      }
      if (mapping.collectionId >= 0 && mapping.sourceIndex >= 0) {
        const auto found = chargedPFIndexBySource.find(sourceKey(mapping.collectionId, mapping.sourceIndex));
        if (found != chargedPFIndexBySource.end()) {
          mapping.chargedPFIndex = found->second;
          mapping.chargedPFMappingValid = 1;
        }
      }
      mapping.filteredTrackMappingValid = mapping.filteredIndex >= 0 ? 1 : 0;
      if (mapping.filteredIndex >= 0) ++nMappedTracks;
      constituentSVIdx.push_back(static_cast<int>(isv));
      constituentSourceSVIdx.push_back(static_cast<int>(isv));
      constituentDaughterIdx.push_back(static_cast<int>(daughterIndex));
      constituentSourceCollectionId.push_back(mapping.collectionId);
      constituentSourceIdx.push_back(mapping.sourceIndex);
      constituentUnpackedTrackIdx.push_back(mapping.unpackedIndex);
      constituentChargedPFIdx.push_back(mapping.chargedPFIndex);
      constituentFilteredTrackIdx.push_back(mapping.filteredIndex);
      constituentMappingValid.push_back(mapping.sourceIdentityValid);
      constituentChargedPFMappingValid.push_back(mapping.chargedPFMappingValid);
      constituentFilteredTrackMappingValid.push_back(mapping.filteredTrackMappingValid);
      constituentPt.push_back(daughter && std::isfinite(daughter->pt()) ? daughter->pt() : kInvalid);
      constituentEta.push_back(daughter && std::isfinite(daughter->eta()) ? daughter->eta() : kInvalid);
      constituentPhi.push_back(daughter && std::isfinite(daughter->phi()) ? daughter->phi() : kInvalid);
      constituentCharge.push_back(daughter ? daughter->charge() : 0);
    }
    svNMappedTracks.push_back(nMappedTracks);
  }

  auto svTable = std::make_unique<nanoaod::FlatTable>(secondaryVertices->size(), "CtxSV", false);
#define ADD_SV(TYPE, NAME, DATA, DOC) svTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_SV(int, "sourceIdx", svSourceIdx, "original row in slimmedSecondaryVertices; no SV selection is applied");
  ADD_SV(float, "pt", svPt, "source SV transverse momentum [GeV]"); ADD_SV(float, "eta", svEta, "source SV pseudorapidity");
  ADD_SV(float, "phi", svPhi, "source SV azimuth [rad]"); ADD_SV(float, "mass", svMass, "source SV mass [GeV]"); ADD_SV(float, "energy", svEnergy, "source SV energy [GeV]");
  ADD_SV(float, "x", svX, "source SV x position [cm]"); ADD_SV(float, "y", svY, "source SV y position [cm]"); ADD_SV(float, "z", svZ, "source SV z position [cm]");
  ADD_SV(float, "chi2", svChi2, "source SV fit chi2"); ADD_SV(float, "ndof", svNdof, "source SV fit degrees of freedom");
  ADD_SV(float, "chi2Ndof", svChi2Ndof, "source SV chi2/ndof; -9999 if ndof is nonpositive");
  ADD_SV(int, "nDaughters", svNDaughters, "all source SV daughter pointers, including unmapped daughters");
  ADD_SV(int, "nMappedTracks", svNMappedTracks, "constituents with filteredTrackIdx >= 0");
  ADD_SV(float, "covXX", svCovXX, "source SV covariance xx [cm^2]"); ADD_SV(float, "covXY", svCovXY, "source SV covariance xy [cm^2]"); ADD_SV(float, "covXZ", svCovXZ, "source SV covariance xz [cm^2]");
  ADD_SV(float, "covYY", svCovYY, "source SV covariance yy [cm^2]"); ADD_SV(float, "covYZ", svCovYZ, "source SV covariance yz [cm^2]"); ADD_SV(float, "covZZ", svCovZZ, "source SV covariance zz [cm^2]");
  ADD_SV(int, "covarianceValid", svCovarianceValid, "1 if source covariance is finite and positive semidefinite");
  ADD_SV(int, "chi2NdofValid", svChi2NdofValid, "1 if chi2/ndof can be computed");
  ADD_SV(int, "pvGeometryValid", svPVGeometryValid, "1 if nonfake PV0 and both PV/SV covariances permit standard distance calculations");
  ADD_SV(float, "lxy", svLxy, "VertexDistanceXY PV0-to-SV distance [cm]");
  ADD_SV(float, "lxyErr", svLxyErr, "VertexDistanceXY uncertainty [cm]");
  ADD_SV(float, "lxySig", svLxySig, "VertexDistanceXY significance");
  ADD_SV(int, "lxyValid", svLxyValid, "1 if the VertexDistanceXY value is finite");
  ADD_SV(int, "lxyErrValid", svLxyErrValid, "1 if the VertexDistanceXY uncertainty is finite and positive");
  ADD_SV(int, "lxySigValid", svLxySigValid, "1 if the VertexDistanceXY significance is finite");
  ADD_SV(float, "l3d", svL3D, "VertexDistance3D PV0-to-SV distance [cm]");
  ADD_SV(float, "l3dErr", svL3DErr, "VertexDistance3D uncertainty [cm]");
  ADD_SV(float, "l3dSig", svL3DSig, "VertexDistance3D significance");
  ADD_SV(int, "l3dValid", svL3DValid, "1 if the VertexDistance3D value is finite");
  ADD_SV(int, "l3dErrValid", svL3DErrValid, "1 if the VertexDistance3D uncertainty is finite and positive");
  ADD_SV(int, "l3dSigValid", svL3DSigValid, "1 if the VertexDistance3D significance is finite");
  ADD_SV(float, "cosPointing2D", svCosPointing2D, "cosine of PV0-to-SV transverse displacement and source SV pT");
  ADD_SV(float, "cosPointing3D", svCosPointing3D, "cosine of PV0-to-SV 3D displacement and source SV momentum");
  ADD_SV(int, "cosPointing2DValid", svCosPointing2DValid, "1 if transverse pointing cosine is defined");
  ADD_SV(int, "cosPointing3DValid", svCosPointing3DValid, "1 if 3D pointing cosine is defined");
#undef ADD_SV
  event.put(std::move(svTable), "CtxSV");

  auto constituentTable = std::make_unique<nanoaod::FlatTable>(constituentSVIdx.size(), "CtxSVConstituent", false);
#define ADD_CONSTITUENT(TYPE, NAME, DATA, DOC) constituentTable->addColumn<TYPE>(NAME, DATA, DOC)
  ADD_CONSTITUENT(int, "svIdx", constituentSVIdx, "event-local CtxSV row; equal to sourceSVIdx because no rows are filtered");
  ADD_CONSTITUENT(int, "sourceSVIdx", constituentSourceSVIdx, "original slimmedSecondaryVertices row");
  ADD_CONSTITUENT(int, "daughterIdx", constituentDaughterIdx, "original daughter ordinal within the source SV");
  ADD_CONSTITUENT(int, "sourceCollectionId", constituentSourceCollectionId, "0=packedPFCandidates, 1=lostTracks, -1=unknown source");
  ADD_CONSTITUENT(int, "sourceIdx", constituentSourceIdx, "original source candidate key; -1 for null Ptr");
  ADD_CONSTITUENT(int, "unpackedTrackIdx", constituentUnpackedTrackIdx, "exact source-to-unpackedTracksAndVertices ValueMap result; -1 if not unpacked");
  ADD_CONSTITUENT(int, "chargedPFIdx", constituentChargedPFIdx, "CtxChargedPF row; -1 if unselected, neutral, or unmapped");
  ADD_CONSTITUENT(int, "filteredTrackIdx", constituentFilteredTrackIdx, "filteredTracks row from source ValueMaps; -1 if not in graph");
  ADD_CONSTITUENT(int, "mappingValid", constituentMappingValid, "1 if original source collection and source key are identified");
  ADD_CONSTITUENT(int, "chargedPFMappingValid", constituentChargedPFMappingValid, "1 if chargedPFIdx references the exact source row");
  ADD_CONSTITUENT(int, "filteredTrackMappingValid", constituentFilteredTrackMappingValid, "1 if filteredTrackIdx references a filtered track row");
  ADD_CONSTITUENT(float, "pt", constituentPt, "source daughter transverse momentum [GeV]");
  ADD_CONSTITUENT(float, "eta", constituentEta, "source daughter pseudorapidity");
  ADD_CONSTITUENT(float, "phi", constituentPhi, "source daughter azimuth [rad]");
  ADD_CONSTITUENT(int, "charge", constituentCharge, "source daughter charge");
#undef ADD_CONSTITUENT
  event.put(std::move(constituentTable), "CtxSVConstituent");
}

DEFINE_FWK_MODULE(TrackGNNVertexContextProducer);
