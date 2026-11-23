#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
using namespace cms::Ort;
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>


class TrackGNNSeederInferenceProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
    explicit TrackGNNSeederInferenceProducer(const edm::ParameterSet&, const ONNXRuntime*);
    void produce(edm::Event &iEvent,const edm::EventSetup &iSetup) override;
    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const ONNXRuntime*);

private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvsToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> passSeedToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> dzToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> dzSignificanceToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ip2dSignificanceToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ip2dValueToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ip3dSignificanceToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ip3dValueToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> trackPairToken_;  // TrackPairFeaturesProducer's output
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;

  bool saveEmbeddings_;

};

TrackGNNSeederInferenceProducer::TrackGNNSeederInferenceProducer(const edm::ParameterSet& iConfig, const ONNXRuntime*):
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  pvsToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
  passSeedToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("passSeedSrc"))),
  dzToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dzSrc"))),
  dzSignificanceToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dzSignificanceSrc"))),
  ip2dSignificanceToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ip2dSignificanceSrc"))),
  ip2dValueToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ip2dValueSrc"))),
  ip3dSignificanceToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ip3dSignificanceSrc"))),
  ip3dValueToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ip3dValueSrc"))),
  trackPairToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("trackPairSrc"))),
  theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  saveEmbeddings_(iConfig.getParameter<bool>("saveEmbeddings"))
{
    produces<nanoaod::FlatTable>("SeederNode");
    if (saveEmbeddings_) produces<nanoaod::FlatTable>("SeederEmbedding");
    produces<nanoaod::FlatTable>("trackPair");  // extension of TrackPairFeaturesProducer table

}
std::unique_ptr<ONNXRuntime> TrackGNNSeederInferenceProducer::initializeGlobalCache(const edm::ParameterSet &iConfig)
{
    return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void TrackGNNSeederInferenceProducer::globalEndJob(const ONNXRuntime*) {}

void TrackGNNSeederInferenceProducer::produce(edm::Event& iEvent,const edm::EventSetup &iSetup)
{

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    edm::Handle<std::vector<reco::Vertex>> pvs;
    iEvent.getByToken(pvsToken_, pvs);
    edm::Handle<edm::ValueMap<int>> passSeedMap;
    iEvent.getByToken(passSeedToken_, passSeedMap);
    edm::Handle<edm::ValueMap<float>> dzMap;
    iEvent.getByToken(dzToken_, dzMap);
    edm::Handle<edm::ValueMap<float>> dzSignificanceMap;
    iEvent.getByToken(dzSignificanceToken_, dzSignificanceMap);
    edm::Handle<edm::ValueMap<float>> ip2dSignificanceMap;
    iEvent.getByToken(ip2dSignificanceToken_, ip2dSignificanceMap);
    edm::Handle<edm::ValueMap<float>> ip2dValueMap;
    iEvent.getByToken(ip2dValueToken_, ip2dValueMap);
    edm::Handle<edm::ValueMap<float>> ip3dSignificanceMap;
    iEvent.getByToken(ip3dSignificanceToken_, ip3dSignificanceMap);
    edm::Handle<edm::ValueMap<float>> ip3dValueMap;
    iEvent.getByToken(ip3dValueToken_, ip3dValueMap);
    edm::Handle<nanoaod::FlatTable> pairTable;
    iEvent.getByToken(trackPairToken_, pairTable);


    if (!tracks.isValid() || !pairTable.isValid() || !passSeedMap.isValid() || !dzMap.isValid() ||
        !dzSignificanceMap.isValid() || !ip2dSignificanceMap.isValid() || !ip2dValueMap.isValid() ||
        !ip3dSignificanceMap.isValid() || !ip3dValueMap.isValid()) {
        throw cms::Exception("ProductNotFound")
            << "TrackGNNSeederInferenceProducer requires tracks, trackPair, and all configured TrackVertexVars ValueMaps";
    }

    const size_t nPairs = pairTable->size();
    const size_t num_tracks = tracks->size();
    std::vector<int> nodeValid(num_tracks, 0), embeddingValid(num_tracks, 0);
    std::vector<int> originalTrackIdx(num_tracks);
    std::vector<float> nodeLogits(num_tracks, -1.f);
    std::vector<float> nodeScores(num_tracks, -1.f);
    std::vector<float> nodeEmbeddings(num_tracks * 128, -1.f);
    std::vector<int> edgeValid(nPairs, 0);
    std::vector<float> edgeLogits(nPairs, -1.f), edgeProbabilities(nPairs, -1.f);
    for (size_t i = 0; i < num_tracks; ++i) originalTrackIdx[i] = static_cast<int>(i);

    auto putTables = [&]() {
        auto nodeTable = std::make_unique<nanoaod::FlatTable>(num_tracks, "SeederNode", false);
        nodeTable->addColumn<int>("trackIdx", originalTrackIdx, "row index in the original filteredTracks collection");
        nodeTable->addColumn<int>("valid", nodeValid, "1 if finite node model outputs were produced, else 0");
        nodeTable->addColumn<float>("logit", nodeLogits,
                                    "raw binary logit for accepted generator-level hadron-decay membership; -1 when invalid");
        nodeTable->addColumn<float>("score", nodeScores,
                                    "sigmoid binary hadron-association score; not prompt-D-specific or necessarily calibrated; -1 when invalid");
        iEvent.put(std::move(nodeTable), "SeederNode");

        if (saveEmbeddings_) {
            auto embeddingTable = std::make_unique<nanoaod::FlatTable>(num_tracks, "SeederEmbedding", false);
            embeddingTable->addColumn<int>("trackIdx", originalTrackIdx, "row index in the original filteredTracks collection");
            embeddingTable->addColumn<int>("valid", embeddingValid, "1 if all 128 embedding values are finite, else 0");
            for (int k = 0; k < 128; ++k) {
                std::vector<float> values(num_tracks);
                for (size_t i = 0; i < num_tracks; ++i) values[i] = nodeEmbeddings[i * 128 + k];
                embeddingTable->addColumn<float>("emb" + std::to_string(k), values, "ONNX node embedding component; -1 when invalid");
            }
            iEvent.put(std::move(embeddingTable), "SeederEmbedding");
        }

        auto edgeTable = std::make_unique<nanoaod::FlatTable>(nPairs, "trackPair", false, true);
        edgeTable->addColumn<float>("seederLogit", edgeLogits, "raw ONNX logit for this directed trackPair row; -1 when invalid");
        edgeTable->addColumn<float>("seederScore", edgeProbabilities, "sigmoid of directed edge logit; -1 when invalid");
        edgeTable->addColumn<int>("seederValid", edgeValid, "1 if this directed edge was evaluated with a finite output");
        iEvent.put(std::move(edgeTable), "trackPair");
    };

    if (!pvs.isValid() || pvs->empty()) {
        putTables();
        return;
    }
    const auto& ttBuilder = iSetup.getData(theTTBToken);

    // Node inputs follow the optimized checkpoint's HDF5 feature order exactly.
    std::vector<float> trk_isHighPurity, trk_charge, trk_nHits, trk_nPixelHits, trk_nStripHits;
    std::vector<float> trk_passSeed, trk_chi2, trk_energy, trk_eta, trk_phi, trk_pt;
    std::vector<float> trk_px, trk_py, trk_pz, trk_qdotp, trk_qoverp;
    std::vector<float> trk_dz, trk_dzSignificance;
    std::vector<float> trk_ip2dSignificance, trk_ip2dValue, trk_ip3dSignificance, trk_ip3dValue;

    std::vector<int> origToNode(num_tracks, -1); // length of tracks. -1 for non valid ttrk. otherwise counts 0,1, 2, -1, 3, ..
    std::vector<int> nodeToOrig;
    nodeToOrig.reserve(num_tracks);
    int node = 0;

    // fill the vectors for each valid transient track
    for (size_t i = 0; i < tracks->size(); ++i) {
        const auto& trk = (*tracks)[i];
        reco::TransientTrack ttrk = ttBuilder.build(trk);

        if (!(ttrk.isValid())) continue;

        origToNode[i] = node;
        nodeToOrig.push_back((int)i);

        const reco::TrackRef trackRef(tracks, i);
        trk_isHighPurity.push_back(trk.quality(reco::TrackBase::highPurity) ? 1.f : 0.f);
        trk_charge.push_back(trk.charge());
        trk_nHits.push_back(trk.hitPattern().numberOfValidHits());
        trk_nPixelHits.push_back(trk.hitPattern().numberOfValidPixelHits());
        trk_nStripHits.push_back(trk.hitPattern().numberOfValidStripHits());
        trk_passSeed.push_back(static_cast<float>((*passSeedMap)[trackRef]));
        trk_chi2.push_back(trk.chi2());
        trk_energy.push_back(std::sqrt(trk.p() * trk.p() + 0.13957 * 0.13957));
        trk_eta.push_back(trk.eta());
        trk_phi.push_back(trk.phi());
        trk_pt.push_back(trk.pt());
        trk_px.push_back(trk.px());
        trk_py.push_back(trk.py());
        trk_pz.push_back(trk.pz());
        trk_qdotp.push_back(trk.charge() * trk.p());
        trk_qoverp.push_back(trk.qoverp());
        trk_dz.push_back((*dzMap)[trackRef]);
        trk_dzSignificance.push_back((*dzSignificanceMap)[trackRef]);
        trk_ip2dSignificance.push_back((*ip2dSignificanceMap)[trackRef]);
        trk_ip2dValue.push_back((*ip2dValueMap)[trackRef]);
        trk_ip3dSignificance.push_back((*ip3dSignificanceMap)[trackRef]);
        trk_ip3dValue.push_back((*ip3dValueMap)[trackRef]);

        node++;
    }

    // Edge (track-pair) features now come straight from TrackPairFeaturesProducer's
    // "trackPair" FlatTable -- read by column name, so this is immune to whatever
    // internal column order that producer happens to use. That producer is the
    // single shared implementation of the pair loop/cuts/geometry, used both for
    // training-ntuple production and here, so there is no risk of the two drifting
    // apart the way the old inline recomputation could.
    auto tp_trk1Idx = pairTable->columnData<int>(pairTable->columnIndex("trk1Idx"));
    auto tp_trk2Idx = pairTable->columnData<int>(pairTable->columnIndex("trk2Idx"));
    auto tp_deltaEta = pairTable->columnData<float>(pairTable->columnIndex("deltaEta"));
    auto tp_deltaPhi = pairTable->columnData<float>(pairTable->columnIndex("deltaPhi"));
    auto tp_pt = pairTable->columnData<float>(pairTable->columnIndex("pt"));
    auto tp_eta = pairTable->columnData<float>(pairTable->columnIndex("eta"));
    auto tp_dca = pairTable->columnData<float>(pairTable->columnIndex("dca"));
    auto tp_deltaR = pairTable->columnData<float>(pairTable->columnIndex("deltaR"));
    auto tp_dcaSig = pairTable->columnData<float>(pairTable->columnIndex("dcaSig"));
    auto tp_cptopv = pairTable->columnData<float>(pairTable->columnIndex("cptopv"));
    auto tp_pvtoPCA1 = pairTable->columnData<float>(pairTable->columnIndex("pvtoPCA1"));
    auto tp_pvtoPCA2 = pairTable->columnData<float>(pairTable->columnIndex("pvtoPCA2"));
    auto tp_dotprod1 = pairTable->columnData<float>(pairTable->columnIndex("dotprod1"));
    auto tp_dotprod2 = pairTable->columnData<float>(pairTable->columnIndex("dotprod2"));
    auto tp_pairMom = pairTable->columnData<float>(pairTable->columnIndex("pairMom"));
    auto tp_pairInvMass = pairTable->columnData<float>(pairTable->columnIndex("pairInvMass"));

    std::vector<int> trk_idxs, trk_jdxs;
    std::vector<size_t> kept_row_idx;  // maps each kept entry back to its row in pairTable
    std::vector<float> pair_deltaEta, pair_deltaPhi, pair_pt, pair_eta;
    std::vector<float> deltaR, dca, dca_sig, cptopv, pvtoPCA_j, pvtoPCA_i, dotprod_j, dotprod_i, pair_mom, pair_invmass;

    for (size_t k = 0; k < nPairs; ++k) {
        int oi = tp_trk1Idx[k];
        int oj = tp_trk2Idx[k];
        // Defensive: TrackPairFeaturesProducer uses the same isValid() criterion to
        // decide which tracks can appear in a pair, so these should always resolve,
        // but skip cleanly rather than index with -1 if that ever changes upstream.
        if (oi < 0 || oj < 0 || static_cast<size_t>(oi) >= origToNode.size() || static_cast<size_t>(oj) >= origToNode.size()) continue;
        if (origToNode[oi] < 0 || origToNode[oj] < 0) continue;

        trk_idxs.push_back(oi);
        trk_jdxs.push_back(oj);
        kept_row_idx.push_back(k);
        pair_deltaEta.push_back(tp_deltaEta[k]);
        pair_deltaPhi.push_back(tp_deltaPhi[k]);
        pair_pt.push_back(tp_pt[k]);
        pair_eta.push_back(tp_eta[k]);
        deltaR.push_back(tp_deltaR[k]);
        dca.push_back(tp_dca[k]);
        dca_sig.push_back(tp_dcaSig[k]);
        cptopv.push_back(tp_cptopv[k]);
        pvtoPCA_i.push_back(tp_pvtoPCA1[k]);
        pvtoPCA_j.push_back(tp_pvtoPCA2[k]);
        dotprod_i.push_back(tp_dotprod1[k]);
        dotprod_j.push_back(tp_dotprod2[k]);
        pair_mom.push_back(tp_pairMom[k]);
        pair_invmass.push_back(tp_pairInvMass[k]);
    }


    // Assemble raw node/edge feature rows in the exact optimized-checkpoint HDF5
    // order. The ONNX graph embeds normalization; do not normalize here.
    std::vector<std::vector<float>> track_features;
    static constexpr int N_NODE_FEATURES = 22;
    static constexpr int N_EDGE_FEATURES = 14;

    for (int ni = 0; ni < node; ++ni) {
        // [0..21] = isHighPurity, charge, nHits, pixelHits, stripHits,
        // passSeed, chi2, energy, eta, phi, pt, px, py, pz, qdotp,
        // qoverp, dz, dzSignificance, ip2dSignificance, ip2dValue,
        // ip3dSignificance, ip3dValue.
        std::vector<float> features = {
            trk_isHighPurity[ni],
            trk_charge[ni],
            trk_nHits[ni],
            trk_nPixelHits[ni],
            trk_nStripHits[ni],
            trk_passSeed[ni],
            trk_chi2[ni],
            trk_energy[ni],
            trk_eta[ni],
            trk_phi[ni],
            trk_pt[ni],
            trk_px[ni],
            trk_py[ni],
            trk_pz[ni],
            trk_qdotp[ni],
            trk_qoverp[ni],
            trk_dz[ni],
            trk_dzSignificance[ni],
            trk_ip2dSignificance[ni],
            trk_ip2dValue[ni],
            trk_ip3dSignificance[ni],
            trk_ip3dValue[ni],
        };

        // Match the HDF5 feature builder component-by-component. Keep finite
        // sentinel values (for example -9999) unchanged.
        for (float& value : features) {
            if (!std::isfinite(value)) value = 0.f;
        }
        track_features.push_back(std::move(features));
    }

    // Build edge_index / per-edge feature rows from the (already column-name-matched)
    // pair data assembled above.
    std::vector<int64_t> edge_i, edge_j;
    std::vector<std::vector<float>> edge_features;
    for (size_t idx = 0; idx < trk_idxs.size(); ++idx) {
        int oi = trk_idxs[idx];
        int oj = trk_jdxs[idx];

        int ni = origToNode[oi];
        int nj = origToNode[oj];

        edge_i.push_back(ni);
        edge_j.push_back(nj);

        // [0..13] = deltaR, deltaPhi, deltaEta, dca, dcaSig, cptopv,
        // pvtoPCA1, pvtoPCA2, dotprod1, dotprod2, pairMom, pairInvMass,
        // pt, eta.
        std::vector<float> features = {
            deltaR[idx],
            pair_deltaPhi[idx],
            pair_deltaEta[idx],
            dca[idx],
            dca_sig[idx],
            cptopv[idx],
            pvtoPCA_i[idx],
            pvtoPCA_j[idx],
            dotprod_i[idx],
            dotprod_j[idx],
            pair_mom[idx],
            pair_invmass[idx],
            pair_pt[idx],
            pair_eta[idx],
        };
        // Match the HDF5 builder's per-component nonfinite replacement.
        for (float& value : features) {
            if (!std::isfinite(value)) value = 0.f;
        }
        edge_features.push_back(std::move(features));
     }

    if (track_features.empty() || edge_features.empty()) {
        putTables();
        return;
    }

    std::vector<float> x_in_flat;
    x_in_flat.reserve(track_features.size() * N_NODE_FEATURES);
    for (const auto& feat : track_features) x_in_flat.insert(x_in_flat.end(), feat.begin(), feat.end());

    std::vector<float> edge_index_flat_f;
    edge_index_flat_f.reserve(edge_i.size() * 2);
    for (int64_t idx : edge_i) edge_index_flat_f.push_back(static_cast<float>(idx));
    for (int64_t idx : edge_j) edge_index_flat_f.push_back(static_cast<float>(idx));

    std::vector<float> edge_attr_flat;
    edge_attr_flat.reserve(edge_features.size() * N_EDGE_FEATURES);
    for (const auto& feat : edge_features) edge_attr_flat.insert(edge_attr_flat.end(), feat.begin(), feat.end());

    const auto allFinite = [](const std::vector<float>& values) {
        return std::all_of(values.begin(), values.end(), [](float value) { return std::isfinite(value); });
    };
    if (!allFinite(x_in_flat) || !allFinite(edge_index_flat_f) || !allFinite(edge_attr_flat)) {
        putTables();
        return;
    }

    std::vector<std::string> input_names = {"x_in", "edge_index", "edge_attr"};
    std::vector<std::string> output_names = {"node_logits", "edge_logits", "node_embeddings"};
    std::vector<std::vector<int64_t>> input_shapes = {
        {1, static_cast<int64_t>(track_features.size()), N_NODE_FEATURES},
        {1, 2, static_cast<int64_t>(edge_i.size())},
        {1, static_cast<int64_t>(edge_features.size()), N_EDGE_FEATURES}
    };
    std::vector<std::vector<float>> data = {x_in_flat, edge_index_flat_f, edge_attr_flat};
    const auto output = globalCache()->run(input_names, data, input_shapes, output_names);
    const size_t expectedNodeLogits = track_features.size();
    const size_t expectedEdgeLogits = edge_i.size();
    const size_t expectedEmbeddings = track_features.size() * 128;
    const bool outputSizesValid = output.size() == 3 && output[0].size() == expectedNodeLogits &&
                                  output[1].size() == expectedEdgeLogits && output[2].size() == expectedEmbeddings;
    if (!outputSizesValid) {
        std::ostringstream details;
        details << "TrackGNNSeederInferenceProducer output size mismatch: internalNodes=" << track_features.size()
                << ", directedEdges=" << edge_i.size() << ", outputCount=" << output.size()
                << ", expected element counts=[node_logits:" << expectedNodeLogits
                << ", edge_logits:" << expectedEdgeLogits << ", node_embeddings:" << expectedEmbeddings << "]"
                << ", actual element counts=[";
        for (size_t i = 0; i < output.size(); ++i) {
            if (i) details << ", ";
            details << (i < output_names.size() ? output_names[i] : "unknown") << ":" << output[i].size();
        }
        details << "]";
        for (const auto& name : output_names) {
            const auto& shape = globalCache()->getOutputShape(name);
            details << "; model " << name << " shape=[";
            for (size_t dim = 0; dim < shape.size(); ++dim) {
                if (dim) details << ",";
                details << shape[dim];
            }
            details << "]";
        }
        throw cms::Exception("InvalidModelOutput") << details.str();
    }

    const auto& rawNodeLogits = output[0];
    const auto& rawEdgeLogits = output[1];
    const auto& rawEmbeddings = output[2];
    const auto sigmoid = [](float value) {
        if (value >= 0.f) return 1.f / (1.f + std::exp(-value));
        const float expValue = std::exp(value);
        return expValue / (1.f + expValue);
    };

    for (size_t ni = 0; ni < nodeToOrig.size(); ++ni) {
        const size_t original = static_cast<size_t>(nodeToOrig[ni]);
        const float logit = rawNodeLogits[ni];
        const float score = sigmoid(logit);
        const bool valid = std::isfinite(logit) && std::isfinite(score);
        bool embeddingIsValid = true;
        for (int k = 0; k < 128; ++k) {
            if (!std::isfinite(rawEmbeddings[ni * 128 + k])) embeddingIsValid = false;
        }
        if (valid) {
            nodeValid[original] = 1;
            nodeLogits[original] = logit;
            nodeScores[original] = score;
        }
        if (embeddingIsValid) {
            embeddingValid[original] = 1;
            for (int k = 0; k < 128; ++k) nodeEmbeddings[original * 128 + k] = rawEmbeddings[ni * 128 + k];
        }
    }

    for (size_t e = 0; e < kept_row_idx.size(); ++e) {
        const float logit = rawEdgeLogits[e];
        if (!std::isfinite(logit)) continue;
        const size_t originalPairRow = kept_row_idx[e];
        if (originalPairRow >= nPairs) continue;
        edgeLogits[originalPairRow] = logit;
        edgeProbabilities[originalPairRow] = sigmoid(logit);
        edgeValid[originalPairRow] = std::isfinite(edgeProbabilities[originalPairRow]) ? 1 : 0;
        if (!edgeValid[originalPairRow]) {
            edgeLogits[originalPairRow] = -1.f;
            edgeProbabilities[originalPairRow] = -1.f;
        }
    }
    putTables();
}

DEFINE_FWK_MODULE(TrackGNNSeederInferenceProducer);
