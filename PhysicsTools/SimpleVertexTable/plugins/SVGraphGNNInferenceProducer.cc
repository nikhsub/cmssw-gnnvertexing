#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Provenance/interface/ProductID.h"
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  constexpr float kPionMass = 0.13957039f;
  constexpr unsigned int kNTrackFeatures = 8;
  constexpr unsigned int kNEdgeFeatures = 10;
  constexpr unsigned int kNGlobalFeatures = 8;
  constexpr unsigned int kNClasses = 4;

  float finiteOrZero(float value) { return std::isfinite(value) ? value : 0.f; }

  float safeRatio(float numerator, float denominator) {
    if (denominator == 0.f || !std::isfinite(denominator)) {
      return 0.f;
    }
    return finiteOrZero(numerator / denominator);
  }

  template <typename T>
  bool containsName(const std::vector<T>& names, const std::string& name) {
    return std::find(names.begin(), names.end(), name) != names.end();
  }
}  // namespace

class SVGraphGNNInferenceProducer : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit SVGraphGNNInferenceProducer(const edm::ParameterSet&, const cms::Ort::ONNXRuntime*);
  void produce(edm::Event&, const edm::EventSetup&) override;

  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const cms::Ort::ONNXRuntime*) {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  struct Candidate {
    int svIdx = -1;
    std::vector<float> trk;
    std::vector<float> edg;
    std::vector<float> eidx;
    std::vector<float> trkValid;
    std::vector<float> edgeValid;
    std::vector<float> glb;
    bool valid = false;
    bool truncatedTracks = false;
    bool truncatedEdges = false;
  };

  Candidate makeCandidate(const reco::Vertex& sv,
                          const reco::Vertex& pv0,
                          int svIdx,
                          const edm::ProductID& trackProductId,
                          const TransientTrackBuilder& ttBuilder) const;

  edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  std::vector<std::string> inputNames_;
  std::vector<std::string> outputNames_;
  unsigned int maxTracks_;
  unsigned int maxEdges_;
  double dlenSigMin_;
  bool debug_;
};

SVGraphGNNInferenceProducer::SVGraphGNNInferenceProducer(const edm::ParameterSet& iConfig,
                                                         const cms::Ort::ONNXRuntime* cache)
    : svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
      ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      inputNames_({"trk", "edg", "eidx", "trk_valid", "edge_valid", "glb"}),
      outputNames_({"vertex_logits", "vertex_probs", "valid_graph_mask"}),
      maxTracks_(iConfig.getParameter<unsigned int>("maxTracks")),
      maxEdges_(iConfig.getParameter<unsigned int>("maxEdges")),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
      debug_(iConfig.getUntrackedParameter<bool>("debug", false)) {
  if (maxTracks_ == 0 || maxEdges_ == 0) {
    throw cms::Exception("Configuration") << "SVGraphGNNInferenceProducer requires positive maxTracks and maxEdges";
  }

  const auto& modelOutputs = cache->getOutputNames();
  for (const auto& outputName : outputNames_) {
    if (!containsName(modelOutputs, outputName)) {
      throw cms::Exception("Configuration") << "GraphVertexGNN model is missing required output '" << outputName << "'";
    }
    (void)cache->getOutputShape(outputName);
  }

  edm::LogInfo("SVGraphGNNInferenceProducer") << "Configured GraphVertexGNN inference with maxTracks=" << maxTracks_
                                              << ", maxEdges=" << maxEdges_ << ", outputs=" << modelOutputs.size();

  produces<nanoaod::FlatTable>("SVGraphGNN");
}

std::unique_ptr<cms::Ort::ONNXRuntime> SVGraphGNNInferenceProducer::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  return std::make_unique<cms::Ort::ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void SVGraphGNNInferenceProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("myFinalInclusiveSecondaryVertices"));
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::FileInPath>("model_path", edm::FileInPath("PhysicsTools/data/GraphVertexGNN.onnx"));
  desc.add<unsigned int>("maxTracks", 16);
  desc.add<unsigned int>("maxEdges", 128);
  desc.add<double>("dlenSigMin", 3.0);
  desc.addUntracked<bool>("debug", false);
  descriptions.addWithDefaultLabel(desc);
}

SVGraphGNNInferenceProducer::Candidate SVGraphGNNInferenceProducer::makeCandidate(
    const reco::Vertex& sv,
    const reco::Vertex& pv0,
    int svIdx,
    const edm::ProductID& trackProductId,
    const TransientTrackBuilder& ttBuilder) const {
  Candidate cand;
  cand.svIdx = svIdx;
  cand.trk.assign(maxTracks_ * kNTrackFeatures, 0.f);
  cand.edg.assign(maxEdges_ * kNEdgeFeatures, 0.f);
  cand.eidx.assign(2 * maxEdges_, 0.f);
  cand.trkValid.assign(maxTracks_, 0.f);
  cand.edgeValid.assign(maxEdges_, 0.f);
  cand.glb.assign(kNGlobalFeatures, 0.f);

  VertexDistance3D vdist;
  const Measurement1D dl =
      vdist.distance(pv0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
  if (!(dl.value() > 0. && dl.significance() > dlenSigMin_)) {
    return cand;
  }

  std::vector<const reco::Track*> tracks;
  tracks.reserve(sv.tracksSize());
  TLorentzVector svP4;
  for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
    const edm::RefToBase<reco::Track>& trkRefBase = *it;
    if (trkRefBase.isNull()) {
      continue;
    }
    reco::TrackRef trkRef = trkRefBase.castTo<reco::TrackRef>();
    if (trkRef.isNull() || trkRef.id() != trackProductId) {
      continue;
    }
    const reco::Track* trk = trkRef.get();
    tracks.push_back(trk);
    TLorentzVector trkP4;
    trkP4.SetPtEtaPhiM(trk->pt(), trk->eta(), trk->phi(), kPionMass);
    svP4 += trkP4;
  }

  if (tracks.size() < 2) {
    return cand;
  }

  cand.truncatedTracks = tracks.size() > maxTracks_;
  const unsigned int nTracks = std::min<unsigned int>(tracks.size(), maxTracks_);

  for (unsigned int i = 0; i < nTracks; ++i) {
    const auto& trk = *tracks[i];
    const unsigned int offset = i * kNTrackFeatures;
    cand.trk[offset + 0] = finiteOrZero(trk.pt());
    cand.trk[offset + 1] = finiteOrZero(trk.eta());
    cand.trk[offset + 2] = finiteOrZero(trk.phi());
    cand.trk[offset + 3] = finiteOrZero(trk.p());
    const float dz = finiteOrZero(trk.dz(pv0.position()));
    cand.trk[offset + 4] = dz;
    cand.trk[offset + 5] = safeRatio(dz, trk.dzError());
    cand.trk[offset + 6] = static_cast<float>(trk.charge());
    cand.trk[offset + 7] = static_cast<float>(trk.numberOfValidHits());
    cand.trkValid[i] = 1.f;
  }

  cand.glb[0] = finiteOrZero(svP4.Pt());
  cand.glb[1] = finiteOrZero(svP4.Eta());
  cand.glb[2] = static_cast<float>(tracks.size());
  cand.glb[3] = finiteOrZero(svP4.M());
  cand.glb[4] = finiteOrZero(dl.value());
  cand.glb[5] = finiteOrZero(dl.significance());
  cand.glb[6] = finiteOrZero(sv.chi2());
  cand.glb[7] = finiteOrZero(sv.ndof());

  unsigned int edge = 0;
  for (unsigned int i = 0; i < nTracks; ++i) {
    TLorentzVector p4i;
    p4i.SetPtEtaPhiM(tracks[i]->pt(), tracks[i]->eta(), tracks[i]->phi(), kPionMass);
    for (unsigned int j = i + 1; j < nTracks; ++j) {
      if (edge >= maxEdges_) {
        cand.truncatedEdges = true;
        break;
      }

      TLorentzVector p4j;
      p4j.SetPtEtaPhiM(tracks[j]->pt(), tracks[j]->eta(), tracks[j]->phi(), kPionMass);

      float dca = -1.f;
      float dcaSig = -1.f;
      float cpToPv = -1.f;
      float pvToPcaSrc = -1.f;
      float pvToPcaDst = -1.f;
      float dotprodSrc = -999.f;
      float dotprodDst = -999.f;
      float pairMom = -1.f;

      const auto ttrkI = ttBuilder.build(*tracks[i]);
      const auto ttrkJ = ttBuilder.build(*tracks[j]);
      TwoTrackMinimumDistance minDist;
      if (ttrkI.isValid() && ttrkJ.isValid() && minDist.calculate(ttrkI.impactPointState(), ttrkJ.impactPointState())) {
        VertexDistance3D distanceComputer;
        auto m = distanceComputer.distance(
            VertexState(minDist.points().second, ttrkI.impactPointState().cartesianError().position()),
            VertexState(minDist.points().first, ttrkJ.impactPointState().cartesianError().position()));
        dca = finiteOrZero(m.value());
        if (m.error() > 0.f) {
          dcaSig = finiteOrZero(m.value() / m.error());
        }
        GlobalPoint pvp(pv0.position().x(), pv0.position().y(), pv0.position().z());
        GlobalPoint cp(minDist.crossingPoint());
        GlobalPoint srcPCA = minDist.points().second;
        GlobalPoint dstPCA = minDist.points().first;
        cpToPv = finiteOrZero((cp - pvp).mag());
        pvToPcaSrc = finiteOrZero((srcPCA - pvp).mag());
        pvToPcaDst = finiteOrZero((dstPCA - pvp).mag());
        dotprodSrc = finiteOrZero((srcPCA - pvp).unit().dot(ttrkI.impactPointState().globalDirection().unit()));
        dotprodDst = finiteOrZero((dstPCA - pvp).unit().dot(ttrkJ.impactPointState().globalDirection().unit()));
        GlobalVector pairMomentum((Basic3DVector<float>)(ttrkI.track().momentum() + ttrkJ.track().momentum()));
        pairMom = finiteOrZero(pairMomentum.mag());
      }

      const unsigned int offset = edge * kNEdgeFeatures;
      cand.edg[offset + 0] = finiteOrZero(p4i.DeltaR(p4j));
      cand.edg[offset + 1] = finiteOrZero((p4i + p4j).M());
      cand.edg[offset + 2] = dca;
      cand.edg[offset + 3] = dcaSig;
      cand.edg[offset + 4] = cpToPv;
      cand.edg[offset + 5] = pvToPcaSrc;
      cand.edg[offset + 6] = pvToPcaDst;
      cand.edg[offset + 7] = dotprodSrc;
      cand.edg[offset + 8] = dotprodDst;
      cand.edg[offset + 9] = pairMom;
      cand.eidx[edge] = static_cast<float>(i);
      cand.eidx[maxEdges_ + edge] = static_cast<float>(j);
      cand.edgeValid[edge] = 1.f;
      ++edge;
    }
    if (cand.truncatedEdges) {
      break;
    }
  }

  cand.valid = nTracks >= 2 && edge > 0;
  return cand;
}

void SVGraphGNNInferenceProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto svs = iEvent.getHandle(svToken_);
  const auto pvs = iEvent.getHandle(pvToken_);
  const auto tracks = iEvent.getHandle(tracksToken_);

  std::vector<int> outSvIdx;
  std::vector<int> outValidGraph;
  std::vector<int> outPredClass;
  std::vector<float> outProbB;
  std::vector<float> outProbDPrompt;
  std::vector<float> outProbDFromB;
  std::vector<float> outProbOther;
  std::vector<float> outLogitB;
  std::vector<float> outLogitDPrompt;
  std::vector<float> outLogitDFromB;
  std::vector<float> outLogitOther;

  std::vector<Candidate> candidates;
  if (svs.isValid() && pvs.isValid() && !pvs->empty() && tracks.isValid()) {
    const auto& pv0 = pvs->front();
    const auto& ttBuilder = iSetup.getData(ttbToken_);
    candidates.reserve(svs->size());
    int svIdx = 0;
    for (const auto& sv : *svs) {
      Candidate cand = makeCandidate(sv, pv0, svIdx, tracks.id(), ttBuilder);
      if (cand.valid) {
        candidates.emplace_back(std::move(cand));
      }
      ++svIdx;
    }
  }

  outSvIdx.reserve(candidates.size());
  outValidGraph.assign(candidates.size(), 0);
  outPredClass.assign(candidates.size(), -1);
  outProbB.assign(candidates.size(), -1.f);
  outProbDPrompt.assign(candidates.size(), -1.f);
  outProbDFromB.assign(candidates.size(), -1.f);
  outProbOther.assign(candidates.size(), -1.f);
  outLogitB.assign(candidates.size(), -999.f);
  outLogitDPrompt.assign(candidates.size(), -999.f);
  outLogitDFromB.assign(candidates.size(), -999.f);
  outLogitOther.assign(candidates.size(), -999.f);
  for (const auto& cand : candidates) {
    outSvIdx.push_back(cand.svIdx);
  }

  if (!candidates.empty()) {
    const int64_t batchSize = candidates.size();
    cms::Ort::FloatArrays inputValues(6);
    inputValues[0].reserve(batchSize * maxTracks_ * kNTrackFeatures);
    inputValues[1].reserve(batchSize * maxEdges_ * kNEdgeFeatures);
    inputValues[2].reserve(batchSize * 2 * maxEdges_);
    inputValues[3].reserve(batchSize * maxTracks_);
    inputValues[4].reserve(batchSize * maxEdges_);
    inputValues[5].reserve(batchSize * kNGlobalFeatures);

    unsigned int nTruncatedTracks = 0;
    unsigned int nTruncatedEdges = 0;
    for (const auto& cand : candidates) {
      inputValues[0].insert(inputValues[0].end(), cand.trk.begin(), cand.trk.end());
      inputValues[1].insert(inputValues[1].end(), cand.edg.begin(), cand.edg.end());
      inputValues[2].insert(inputValues[2].end(), cand.eidx.begin(), cand.eidx.end());
      inputValues[3].insert(inputValues[3].end(), cand.trkValid.begin(), cand.trkValid.end());
      inputValues[4].insert(inputValues[4].end(), cand.edgeValid.begin(), cand.edgeValid.end());
      inputValues[5].insert(inputValues[5].end(), cand.glb.begin(), cand.glb.end());
      nTruncatedTracks += cand.truncatedTracks;
      nTruncatedEdges += cand.truncatedEdges;
    }

    const std::vector<std::vector<int64_t>> inputShapes = {
        {batchSize, static_cast<int64_t>(maxTracks_), static_cast<int64_t>(kNTrackFeatures)},
        {batchSize, static_cast<int64_t>(maxEdges_), static_cast<int64_t>(kNEdgeFeatures)},
        {batchSize, 2, static_cast<int64_t>(maxEdges_)},
        {batchSize, static_cast<int64_t>(maxTracks_)},
        {batchSize, static_cast<int64_t>(maxEdges_)},
        {batchSize, static_cast<int64_t>(kNGlobalFeatures)}};

    const auto outputs = globalCache()->run(inputNames_, inputValues, inputShapes, outputNames_, batchSize);
    if (outputs.size() != 3) {
      throw cms::Exception("RuntimeError") << "GraphVertexGNN returned " << outputs.size() << " outputs, expected 3";
    }
    const auto& logits = outputs[0];
    const auto& probs = outputs[1];
    const auto& validMask = outputs[2];
    if (validMask.size() != candidates.size()) {
      throw cms::Exception("RuntimeError")
          << "GraphVertexGNN valid_graph_mask has length " << validMask.size() << ", expected " << candidates.size();
    }

    std::size_t compressedRow = 0;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      const bool valid = validMask[i] > 0.5f;
      outValidGraph[i] = valid ? 1 : 0;
      if (!valid) {
        continue;
      }
      const std::size_t base = compressedRow * kNClasses;
      if (base + kNClasses > probs.size() || base + kNClasses > logits.size()) {
        throw cms::Exception("RuntimeError") << "GraphVertexGNN compressed output row count is inconsistent with "
                                             << "valid_graph_mask";
      }
      outLogitB[i] = logits[base + 0];
      outLogitDPrompt[i] = logits[base + 1];
      outLogitDFromB[i] = logits[base + 2];
      outLogitOther[i] = logits[base + 3];
      outProbB[i] = probs[base + 0];
      outProbDPrompt[i] = probs[base + 1];
      outProbDFromB[i] = probs[base + 2];
      outProbOther[i] = probs[base + 3];
      const auto begin = probs.begin() + base;
      outPredClass[i] = static_cast<int>(std::distance(begin, std::max_element(begin, begin + kNClasses)));
      ++compressedRow;
    }
    if (compressedRow * kNClasses != probs.size()) {
      throw cms::Exception("RuntimeError") << "GraphVertexGNN produced " << probs.size() / kNClasses
                                           << " probability rows, but valid_graph_mask selected " << compressedRow;
    }

    if (debug_) {
      edm::LogInfo("SVGraphGNNInferenceProducer")
          << "Processed " << candidates.size() << " candidates, valid ONNX rows=" << compressedRow
          << ", track-truncated=" << nTruncatedTracks << ", edge-truncated=" << nTruncatedEdges;
    }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(outSvIdx.size(), "SVGraphGNN", false);
  table->addColumn<int>("svIdx", outSvIdx, "Original secondary-vertex index in the input collection");
  table->addColumn<int>(
      "validGraph", outValidGraph, "GraphVertexGNN valid_graph_mask mapped to the original candidate");
  table->addColumn<int>("predClass", outPredClass, "Argmax class: 0=isB, 1=isD_prompt, 2=isD_fromB, 3=isOther");
  table->addColumn<float>("prob_isB", outProbB, "GraphVertexGNN probability for class isB");
  table->addColumn<float>("prob_isD_prompt", outProbDPrompt, "GraphVertexGNN probability for class isD_prompt");
  table->addColumn<float>("prob_isD_fromB", outProbDFromB, "GraphVertexGNN probability for class isD_fromB");
  table->addColumn<float>("prob_isOther", outProbOther, "GraphVertexGNN probability for class isOther");
  table->addColumn<float>("logit_isB", outLogitB, "GraphVertexGNN logit for class isB");
  table->addColumn<float>("logit_isD_prompt", outLogitDPrompt, "GraphVertexGNN logit for class isD_prompt");
  table->addColumn<float>("logit_isD_fromB", outLogitDFromB, "GraphVertexGNN logit for class isD_fromB");
  table->addColumn<float>("logit_isOther", outLogitOther, "GraphVertexGNN logit for class isOther");
  iEvent.put(std::move(table), "SVGraphGNN");
}

DEFINE_FWK_MODULE(SVGraphGNNInferenceProducer);
