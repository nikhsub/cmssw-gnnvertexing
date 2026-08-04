#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
  constexpr float kPionMass = 0.13957039f;
  constexpr unsigned int kNTrackFeatures = 16;
  constexpr unsigned int kNEdgeFeatures = 10;
  constexpr unsigned int kNGlobalFeatures = 24;
  constexpr unsigned int kNClasses = 4;

  float safeRatio(float numerator, float denominator) {
    return denominator > 0.f ? numerator / denominator : 0.f;
  }

  template <typename T> bool containsName(const std::vector<T>& names, const std::string& name) {
    return std::find(names.begin(), names.end(), name) != names.end();
  }
}

class SVGraphVertexGNNInferenceProducer : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit SVGraphVertexGNNInferenceProducer(const edm::ParameterSet&, const cms::Ort::ONNXRuntime*);
  void produce(edm::Event&, const edm::EventSetup&) override;

  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const cms::Ort::ONNXRuntime*) {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  struct GraphTrackInfo {
    int globalIdx = -1;
    const reco::Track* trk = nullptr;
    int isFitConstituent = 0;
    int isConeConstituent = 0;
    float dRToSVAxis = -1.f;
    float dzPV0 = 0.f;
    float dzSigPV0 = 0.f;
    float absDzPV0 = 0.f;
    float dxyPV0 = 0.f;
    float dxySigPV0 = 0.f;
    int closestPVIdx = -1;
    float dzClosestPV = 999.f;
    int isClosestPV0 = 0;
  };

  struct Candidate {
    int svIdx = -1;
    std::vector<float> trk, edg, eidx, trkValid, edgeValid, glb;
    bool valid = false;
    bool truncatedTracks = false;
    bool truncatedEdges = false;
  };

  GraphTrackInfo makeGraphTrackInfo(const reco::Track*, int, int, int, float, const std::vector<reco::Vertex>&) const;
  Candidate makeCandidate(const reco::Vertex&,
                          const std::vector<reco::Vertex>&,
                          int,
                          const edm::Handle<reco::TrackCollection>&,
                          const edm::ValueMap<int>&,
                          const TransientTrackBuilder&) const;

  edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_, pvToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> globalTrackIdxToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  std::vector<std::string> inputNames_;
  std::vector<std::string> outputNames_;
  unsigned int maxTracks_, maxEdges_;
  double dlenSigMin_, nearbyTrackDR_, nearbyTrackPtMin_, nearbyTrackEtaMax_;
  int maxExtraTracks_;
  bool includeNearbyTracks_, requireExtraHighPurity_, debug_;
};

SVGraphVertexGNNInferenceProducer::SVGraphVertexGNNInferenceProducer(const edm::ParameterSet& iConfig,
                                                                     const cms::Ort::ONNXRuntime* cache) :
    svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
    pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
    globalTrackIdxToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("globalTrackIdxMap"))),
    ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    inputNames_({"trk", "edg", "eidx", "trk_valid", "edge_valid", "glb"}),
    outputNames_({"flavor_logits", "flavor_probs", "reliability_logits", "reliability_score", "valid_graph_mask"}),
    maxTracks_(iConfig.getParameter<unsigned int>("maxTracks")),
    maxEdges_(iConfig.getParameter<unsigned int>("maxEdges")),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    nearbyTrackDR_(iConfig.getParameter<double>("nearbyTrackDR")),
    nearbyTrackPtMin_(iConfig.getParameter<double>("nearbyTrackPtMin")),
    nearbyTrackEtaMax_(iConfig.getParameter<double>("nearbyTrackEtaMax")),
    maxExtraTracks_(iConfig.getParameter<int>("maxExtraTracks")),
    includeNearbyTracks_(iConfig.getParameter<bool>("includeNearbyTracks")),
    requireExtraHighPurity_(iConfig.getParameter<bool>("requireExtraHighPurity")),
    debug_(iConfig.getUntrackedParameter<bool>("debug", false)) {
  if (maxTracks_ == 0 || maxEdges_ == 0)
    throw cms::Exception("Configuration") << "maxTracks and maxEdges must be positive";
  const auto& modelOutputs = cache->getOutputNames();
  for (const auto& outputName : outputNames_) {
    if (!containsName(modelOutputs, outputName))
      throw cms::Exception("Configuration") << "GraphVertexGNN model is missing required output '" << outputName << "'";
    (void)cache->getOutputShape(outputName);
  }
  produces<nanoaod::FlatTable>("SVGraphVertexGNN");
}

std::unique_ptr<cms::Ort::ONNXRuntime>
SVGraphVertexGNNInferenceProducer::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  return std::make_unique<cms::Ort::ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void SVGraphVertexGNNInferenceProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("myFinalInclusiveSecondaryVertices"));
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("globalTrackIdxMap", edm::InputTag("dummyValueMap", "globalTrackIdxMap"));
  desc.add<edm::FileInPath>("model_path");
  desc.add<unsigned int>("maxTracks", 32);
  desc.add<unsigned int>("maxEdges", 128);
  desc.add<double>("dlenSigMin", 3.0);
  desc.add<bool>("includeNearbyTracks", false);
  desc.add<double>("nearbyTrackDR", 0.3);
  desc.add<double>("nearbyTrackPtMin", 0.8);
  desc.add<double>("nearbyTrackEtaMax", 2.5);
  desc.add<int>("maxExtraTracks", 4);
  desc.add<bool>("requireExtraHighPurity", false);
  desc.addUntracked<bool>("debug", false);
  descriptions.addWithDefaultLabel(desc);
}

SVGraphVertexGNNInferenceProducer::GraphTrackInfo
SVGraphVertexGNNInferenceProducer::makeGraphTrackInfo(const reco::Track* trk,
                                                      int globalIdx,
                                                      int isFitConstituent,
                                                      int isConeConstituent,
                                                      float dRToSVAxis,
                                                      const std::vector<reco::Vertex>& pvs) const {
  GraphTrackInfo info;
  info.globalIdx = globalIdx;
  info.trk = trk;
  info.isFitConstituent = isFitConstituent;
  info.isConeConstituent = isConeConstituent;
  info.dRToSVAxis = dRToSVAxis;
  if (!trk || pvs.empty())
    return info;
  const auto& pv0 = pvs.front();
  const float dzPV0 = (trk->dz(pv0.position()));
  info.dzPV0 = dzPV0;
  info.absDzPV0 = std::abs(dzPV0);
  info.dzSigPV0 = safeRatio(dzPV0, trk->dzError());
  const float dxyPV0 = (trk->dxy(pv0.position()));
  info.dxyPV0 = dxyPV0;
  info.dxySigPV0 = safeRatio(dxyPV0, trk->dxyError());
  for (size_t ipv = 0; ipv < pvs.size(); ++ipv) {
    const float absDz = std::abs((trk->dz(pvs[ipv].position())));
    if (absDz < info.dzClosestPV) {
      info.dzClosestPV = absDz;
      info.closestPVIdx = static_cast<int>(ipv);
    }
  }
  info.isClosestPV0 = (info.closestPVIdx == 0 ? 1 : 0);
  return info;
}

SVGraphVertexGNNInferenceProducer::Candidate
SVGraphVertexGNNInferenceProducer::makeCandidate(const reco::Vertex& sv,
                                                 const std::vector<reco::Vertex>& pvs,
                                                 int svIdx,
                                                 const edm::Handle<reco::TrackCollection>& tracksHandle,
                                                 const edm::ValueMap<int>& globalIdxMap,
                                                 const TransientTrackBuilder& ttBuilder) const {
  const auto& tracks = *tracksHandle;
  Candidate cand;
  cand.svIdx = svIdx;
  cand.trk.assign(maxTracks_ * kNTrackFeatures, 0.f);
  cand.edg.assign(maxEdges_ * kNEdgeFeatures, 0.f);
  cand.eidx.assign(2 * maxEdges_, 0.f);
  cand.trkValid.assign(maxTracks_, 0.f);
  cand.edgeValid.assign(maxEdges_, 0.f);
  cand.glb.assign(kNGlobalFeatures, 0.f);
  const auto& pv0 = pvs.front();
  VertexDistance3D vdist;
  const Measurement1D dl =
      vdist.distance(pv0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
  if (!(dl.value() > 0. && dl.significance() > dlenSigMin_))
    return cand;

  TVector3 flightVec(sv.x() - pv0.x(), sv.y() - pv0.y(), sv.z() - pv0.z());
  float svAxisEta = 0.f;
  float svAxisPhi = 0.f;
  if (flightVec.Mag() > 0.f) {
    svAxisEta = flightVec.Eta();
    svAxisPhi = flightVec.Phi();
  }
  std::vector<GraphTrackInfo> fitTracks, graphTracks;
  std::unordered_set<int> usedGlobalIdx;
  TLorentzVector fitP4;
  for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
    const edm::RefToBase<reco::Track>& trkRefBase = *it;
    if (trkRefBase.isNull())
      continue;
    reco::TrackRef trkRef = trkRefBase.castTo<reco::TrackRef>();
    if (trkRef.isNull() || trkRef.id() != tracksHandle.id())
      continue;
    const int globalIdx = globalIdxMap[trkRef];
    if (globalIdx < 0)
      continue;
    const reco::Track* trk = trkRef.get();
    const float dRToAxis =
        reco::deltaR(trk->eta(), trk->phi(), static_cast<double>(svAxisEta), static_cast<double>(svAxisPhi));
    auto info = makeGraphTrackInfo(trk, globalIdx, 1, 0, dRToAxis, pvs);
    fitTracks.push_back(info);
    graphTracks.push_back(info);
    usedGlobalIdx.insert(globalIdx);
    TLorentzVector trkP4;
    trkP4.SetPtEtaPhiM(trk->pt(), trk->eta(), trk->phi(), kPionMass);
    fitP4 += trkP4;
  }
  if (fitTracks.size() < 2)
    return cand;
  if (!(std::isfinite(svAxisEta) && std::isfinite(svAxisPhi)) && fitP4.Pt() > 0.f) {
    svAxisEta = fitP4.Eta();
    svAxisPhi = fitP4.Phi();
  }
  if (includeNearbyTracks_) {
    std::vector<std::pair<float, GraphTrackInfo>> extraCandidates;
    for (size_t itrk = 0; itrk < tracks.size(); ++itrk) {
      reco::TrackRef trkRef(tracksHandle, itrk);
      if (trkRef.isNull())
        continue;
      const int globalIdx = globalIdxMap[trkRef];
      if (globalIdx < 0 || usedGlobalIdx.count(globalIdx))
        continue;
      const reco::Track* trk = trkRef.get();
      if (trk->pt() < nearbyTrackPtMin_ || std::abs(trk->eta()) > nearbyTrackEtaMax_)
        continue;
      if (requireExtraHighPurity_ && !trk->quality(reco::TrackBase::highPurity))
        continue;
      const float dRToAxis =
          reco::deltaR(trk->eta(), trk->phi(), static_cast<double>(svAxisEta), static_cast<double>(svAxisPhi));
      if (dRToAxis > nearbyTrackDR_)
        continue;
      extraCandidates.emplace_back(dRToAxis, makeGraphTrackInfo(trk, globalIdx, 0, 1, dRToAxis, pvs));
    }
    std::sort(extraCandidates.begin(), extraCandidates.end(), [](const auto& a, const auto& b) {
      return a.first < b.first;
    });
    int nAdded = 0;
    for (const auto& extra : extraCandidates) {
      if (maxExtraTracks_ >= 0 && nAdded >= maxExtraTracks_)
        break;
      graphTracks.push_back(extra.second);
      ++nAdded;
    }
  }
  if (graphTracks.size() < 2)
    return cand;
  cand.truncatedTracks = graphTracks.size() > maxTracks_;
  const unsigned int nTracks = std::min<unsigned int>(graphTracks.size(), maxTracks_);

  int nFitClosestPV0 = 0;
  std::unordered_set<int> fitDistinctPVs;
  std::unordered_map<int, int> fitPVCounts;
  for (const auto& info : fitTracks) {
    if (info.closestPVIdx >= 0) {
      fitDistinctPVs.insert(info.closestPVIdx);
      ++fitPVCounts[info.closestPVIdx];
    }
    if (info.isClosestPV0)
      ++nFitClosestPV0;
  }
  int fitMajorityPVCount = 0;
  for (const auto& kv : fitPVCounts)
    if (kv.second > fitMajorityPVCount)
      fitMajorityPVCount = kv.second;
  int svClosestPVIdx = -1;
  float svDzClosestPV = 999.f;
  for (size_t ipv = 0; ipv < pvs.size(); ++ipv) {
    const float absDz = std::abs((sv.z() - pvs[ipv].z()));
    if (absDz < svDzClosestPV) {
      svDzClosestPV = absDz;
      svClosestPVIdx = static_cast<int>(ipv);
    }
  }

  cand.glb[0] = (fitP4.Pt());
  cand.glb[1] = (fitP4.Eta());
  cand.glb[2] = (fitP4.Phi());
  cand.glb[3] = (sv.x());
  cand.glb[4] = (sv.y());
  cand.glb[5] = (sv.z());
  cand.glb[6] = static_cast<float>(graphTracks.size());
  cand.glb[7] = static_cast<float>(graphTracks.size());
  cand.glb[8] = static_cast<float>(fitTracks.size());
  cand.glb[9] = static_cast<float>(graphTracks.size() - fitTracks.size());
  cand.glb[10] = svClosestPVIdx == 0 ? 1.f : 0.f;
  cand.glb[11] = static_cast<float>(fitDistinctPVs.size());
  cand.glb[12] = (fitP4.M());
  cand.glb[13] = (dl.value());
  cand.glb[14] = (dl.significance());
  cand.glb[15] = (sv.chi2());
  cand.glb[16] = (sv.ndof());
  cand.glb[17] = (sv.z() - pv0.z());
  cand.glb[18] = std::abs(cand.glb[17]);
  cand.glb[19] = (std::hypot(sv.x() - pv0.x(), sv.y() - pv0.y()));
  cand.glb[20] = svDzClosestPV;
  cand.glb[21] = safeRatio(static_cast<float>(nFitClosestPV0), static_cast<float>(fitTracks.size()));
  cand.glb[22] = 1.f - cand.glb[21];
  cand.glb[23] = safeRatio(static_cast<float>(fitMajorityPVCount), static_cast<float>(fitTracks.size()));

  for (unsigned int i = 0; i < nTracks; ++i) {
    const auto& info = graphTracks[i];
    const auto& trk = *info.trk;
    const unsigned int off = i * kNTrackFeatures;
    cand.trk[off + 0] = static_cast<float>(info.isFitConstituent);
    cand.trk[off + 1] = static_cast<float>(info.isConeConstituent);
    cand.trk[off + 2] = static_cast<float>(info.isClosestPV0);
    cand.trk[off + 3] = (info.dRToSVAxis);
    cand.trk[off + 4] = (trk.pt());
    cand.trk[off + 5] = (trk.eta());
    cand.trk[off + 6] = (trk.phi());
    cand.trk[off + 7] = (trk.p());
    cand.trk[off + 8] = info.dzPV0;
    cand.trk[off + 9] = info.dzSigPV0;
    cand.trk[off + 10] = info.absDzPV0;
    cand.trk[off + 11] = info.dxyPV0;
    cand.trk[off + 12] = info.dxySigPV0;
    cand.trk[off + 13] = info.dzClosestPV;
    cand.trk[off + 14] = static_cast<float>(trk.charge());
    cand.trk[off + 15] = static_cast<float>(trk.numberOfValidHits());
    cand.trkValid[i] = 1.f;
  }

  unsigned int edge = 0;
  for (unsigned int i = 0; i < nTracks; ++i) {
    TLorentzVector p4i;
    p4i.SetPtEtaPhiM(graphTracks[i].trk->pt(), graphTracks[i].trk->eta(), graphTracks[i].trk->phi(), kPionMass);
    for (unsigned int j = i + 1; j < nTracks; ++j) {
      if (edge >= maxEdges_) {
        cand.truncatedEdges = true;
        break;
      }
      TLorentzVector p4j;
      p4j.SetPtEtaPhiM(graphTracks[j].trk->pt(), graphTracks[j].trk->eta(), graphTracks[j].trk->phi(), kPionMass);
      float dca = -1.f,
            dcaSig = -1.f,
            cpToPv = -1.f,
            pvToPcaSrc = -1.f,
            pvToPcaDst = -1.f,
            dotprodSrc = -999.f,
            dotprodDst = -999.f,
            pairMom = -1.f;
      const auto ttrkI = ttBuilder.build(*graphTracks[i].trk);
      const auto ttrkJ = ttBuilder.build(*graphTracks[j].trk);
      TwoTrackMinimumDistance minDist;
      if (ttrkI.isValid() && ttrkJ.isValid() && minDist.calculate(ttrkI.impactPointState(), ttrkJ.impactPointState())) {
        VertexDistance3D distanceComputer;
        auto m = distanceComputer.distance(VertexState(minDist.points().second, ttrkI.impactPointState().cartesianError().position()),
            VertexState(minDist.points().first, ttrkJ.impactPointState().cartesianError().position()));
        dca = (m.value());
        if (m.error() > 0.f)
          dcaSig = (m.value() / m.error());
        GlobalPoint pvp(pv0.position().x(), pv0.position().y(), pv0.position().z());
        GlobalPoint cp(minDist.crossingPoint());
        GlobalPoint srcPCA = minDist.points().second;
        GlobalPoint dstPCA = minDist.points().first;
        cpToPv = ((cp - pvp).mag());
        pvToPcaSrc = ((srcPCA - pvp).mag());
        pvToPcaDst = ((dstPCA - pvp).mag());
        dotprodSrc = ((srcPCA - pvp).unit().dot(ttrkI.impactPointState().globalDirection().unit()));
        dotprodDst = ((dstPCA - pvp).unit().dot(ttrkJ.impactPointState().globalDirection().unit()));
        GlobalVector pairMomentum((Basic3DVector<float>)(ttrkI.track().momentum() + ttrkJ.track().momentum()));
        pairMom = (pairMomentum.mag());
      }
      const unsigned int off = edge * kNEdgeFeatures;
      cand.edg[off + 0] = (p4i.DeltaR(p4j));
      cand.edg[off + 1] = ((p4i + p4j).M());
      cand.edg[off + 2] = dca;
      cand.edg[off + 3] = dcaSig;
      cand.edg[off + 4] = cpToPv;
      cand.edg[off + 5] = pvToPcaSrc;
      cand.edg[off + 6] = pvToPcaDst;
      cand.edg[off + 7] = dotprodSrc;
      cand.edg[off + 8] = dotprodDst;
      cand.edg[off + 9] = pairMom;
      cand.eidx[edge] = static_cast<float>(i);
      cand.eidx[maxEdges_ + edge] = static_cast<float>(j);
      cand.edgeValid[edge] = 1.f;
      ++edge;
    }
    if (cand.truncatedEdges)
      break;
  }
  cand.valid = nTracks >= 2 && edge > 0;
  return cand;
}

void SVGraphVertexGNNInferenceProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto svs = iEvent.getHandle(svToken_);
  const auto pvs = iEvent.getHandle(pvToken_);
  const auto tracks = iEvent.getHandle(tracksToken_);
  const auto globalIdxMap = iEvent.getHandle(globalTrackIdxToken_);
  std::vector<Candidate> candidates;
  if (svs.isValid() && pvs.isValid() && !pvs->empty() && tracks.isValid() && !tracks->empty() && globalIdxMap.isValid()) {
    const auto& ttBuilder = iSetup.getData(ttbToken_);
    int keptSvIdx = 0;
    candidates.reserve(svs->size());

    for (std::size_t svIdx = 0; svIdx < svs->size(); ++svIdx) {
      auto cand = makeCandidate(
          svs->at(svIdx),
          *pvs,
          static_cast<int>(svIdx),
          tracks,
          *globalIdxMap,
          ttBuilder
      );
    
      if (cand.valid) {
        candidates.emplace_back(std::move(cand));
      }
}
    
  }
  const std::size_t nCand = candidates.size();
  std::vector<int> outSvIdx;
  outSvIdx.reserve(nCand);
  std::vector<int> outValidGraph(nCand, 0), outPredClass(nCand, -1);
  std::vector<float> outProbB(nCand, -1.f),
                     outProbDPrompt(nCand, -1.f),
                     outProbDFromB(nCand, -1.f),
                     outProbOther(nCand, -1.f),
                     outReliabilityScore(nCand, -1.f),
                     outReliabilityLogit(nCand, -999.f);
  for (const auto& cand : candidates)
    outSvIdx.push_back(cand.svIdx);
  if (!candidates.empty()) {
    const int64_t batchSize = static_cast<int64_t>(nCand);
    cms::Ort::FloatArrays inputValues(6);
    for (auto& values : inputValues)
      values.reserve(batchSize * maxEdges_ * kNEdgeFeatures);
    for (const auto& cand : candidates) {
      inputValues[0].insert(inputValues[0].end(), cand.trk.begin(), cand.trk.end());
      inputValues[1].insert(inputValues[1].end(), cand.edg.begin(), cand.edg.end());
      inputValues[2].insert(inputValues[2].end(), cand.eidx.begin(), cand.eidx.end());
      inputValues[3].insert(inputValues[3].end(), cand.trkValid.begin(), cand.trkValid.end());
      inputValues[4].insert(inputValues[4].end(), cand.edgeValid.begin(), cand.edgeValid.end());
      inputValues[5].insert(inputValues[5].end(), cand.glb.begin(), cand.glb.end());
    }
    const std::vector<std::vector<int64_t>> inputShapes = {{batchSize, static_cast<int64_t>(maxTracks_), static_cast<int64_t>(kNTrackFeatures)},
        {batchSize, static_cast<int64_t>(maxEdges_), static_cast<int64_t>(kNEdgeFeatures)},
        {batchSize, 2, static_cast<int64_t>(maxEdges_)},
        {batchSize, static_cast<int64_t>(maxTracks_)},
        {batchSize, static_cast<int64_t>(maxEdges_)},
        {batchSize, static_cast<int64_t>(kNGlobalFeatures)}};
    const auto outputs = globalCache()->run(inputNames_, inputValues, inputShapes, outputNames_, batchSize);
    if (outputs.size() != 5)
      throw cms::Exception("RuntimeError") << "GraphVertexGNN returned " << outputs.size() << " outputs, expected 5";
    const auto& probs = outputs[1];
    const auto& relLogits = outputs[2];
    const auto& relScores = outputs[3];
    const auto& validMask = outputs[4];
    if (probs.size() != nCand * kNClasses || relLogits.size() != nCand || relScores.size() != nCand ||
        validMask.size() != nCand)
      throw cms::Exception("RuntimeError") << "GraphVertexGNN output shape mismatch";
    for (std::size_t i = 0; i < nCand; ++i) {
      const std::size_t base = i * kNClasses;
      outValidGraph[i] = validMask[i] > 0.5f ? 1 : 0;
      outProbB[i] = probs[base];
      outProbDPrompt[i] = probs[base + 1];
      outProbDFromB[i] = probs[base + 2];
      outProbOther[i] = probs[base + 3];
      outReliabilityLogit[i] = relLogits[i];
      outReliabilityScore[i] = relScores[i];
      auto begin = probs.begin() + base;
      outPredClass[i] = static_cast<int>(std::distance(begin, std::max_element(begin, begin + kNClasses)));
    }
    if (debug_)
      edm::LogInfo("SVGraphVertexGNNInferenceProducer") << "Ran GraphVertexGNN on " << nCand << " vertices";
  }
  auto table = std::make_unique<nanoaod::FlatTable>(outSvIdx.size(), "SVGraphVertexGNN", false);
  table->addColumn<int>("svIdx", outSvIdx, "Original secondary-vertex index in the input collection");
  table->addColumn<int>("validGraph", outValidGraph, "ONNX valid graph mask");
  table->addColumn<int>("predClass", outPredClass, "Argmax flavor class: 0=isB, 1=isD_prompt, 2=isD_fromB, 3=isOther");
  table->addColumn<float>("prob_isB", outProbB, "GraphVertexGNN probability for class isB");
  table->addColumn<float>("prob_isD_prompt", outProbDPrompt, "GraphVertexGNN probability for class isD_prompt");
  table->addColumn<float>("prob_isD_fromB", outProbDFromB, "GraphVertexGNN probability for class isD_fromB");
  table->addColumn<float>("prob_isOther", outProbOther, "GraphVertexGNN probability for class isOther");
  table->addColumn<float>("reliability", outReliabilityScore, "GraphVertexGNN reliability score");
  table->addColumn<float>("reliabilityLogit", outReliabilityLogit, "GraphVertexGNN reliability logit");
  iEvent.put(std::move(table), "SVGraphVertexGNN");
}

DEFINE_FWK_MODULE(SVGraphVertexGNNInferenceProducer);
