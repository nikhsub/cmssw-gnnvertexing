#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
using namespace cms::Ort;
#include <iostream>
#include <cmath>
#include <limits>
#include <array>
#include <string>
#include "DataFormats/Common/interface/Wrapper.h"
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
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronsToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> trackPairToken_;  // TrackPairFeaturesProducer's output
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;

  // Same style/meaning as TrackVertexVars' lepton-match window. IMPORTANT: the *training*
  // config (trackVertexVars_cfi in the python fragment) sets these to 0.001 / 0.2, which
  // differ from TrackVertexVars.cc's own fillDescriptions() defaults (0.01 / 0.1). Match
  // whatever was actually used in the ntuple production for these values.
  double dRMatchMax_;
  double relPtMatchMax_;
  bool saveEmbeddings_;

  static constexpr float PION_MASS = 0.13957018f;
  static constexpr int N_NODE_CLASSES = 6;  // sv_pred's output width (GraphEventGNN.sv_pred -> 6-way)
};

TrackGNNSeederInferenceProducer::TrackGNNSeederInferenceProducer(const edm::ParameterSet& iConfig, const ONNXRuntime*):
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  pvsToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  trackPairToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("trackPairSrc"))),
  theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  dRMatchMax_(iConfig.getParameter<double>("dRMatchMax")),
  relPtMatchMax_(iConfig.getParameter<double>("relPtMatchMax")),
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
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotToken_, beamSpot);
    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);
    edm::Handle<std::vector<pat::Electron>> electrons;
    iEvent.getByToken(electronsToken_, electrons);
    edm::Handle<nanoaod::FlatTable> pairTable;
    iEvent.getByToken(trackPairToken_, pairTable);


    if (!tracks.isValid() || !pairTable.isValid() || !beamSpot.isValid() || !muons.isValid() || !electrons.isValid()) {
        throw cms::Exception("ProductNotFound") << "TrackGNNSeederInferenceProducer requires valid track, pair, beam-spot, muon, and electron inputs";
    }

    const size_t nPairs = pairTable->size();
    const size_t num_tracks = tracks->size();
    std::vector<int> nodeValid(num_tracks, 0), embeddingValid(num_tracks, 0);
    std::vector<int> originalTrackIdx(num_tracks);
    std::vector<float> nodeLogits(num_tracks * N_NODE_CLASSES, -1.f);
    std::vector<float> nodeProbabilities(num_tracks * N_NODE_CLASSES, -1.f);
    std::vector<float> nodeEmbeddings(num_tracks * 128, -1.f);
    std::vector<int> edgeValid(nPairs, 0);
    std::vector<float> edgeLogits(nPairs, -1.f), edgeProbabilities(nPairs, -1.f);
    for (size_t i = 0; i < num_tracks; ++i) originalTrackIdx[i] = static_cast<int>(i);

    auto putTables = [&]() {
        auto nodeTable = std::make_unique<nanoaod::FlatTable>(num_tracks, "SeederNode", false);
        nodeTable->addColumn<int>("trackIdx", originalTrackIdx, "row index in the original filteredTracks collection");
        nodeTable->addColumn<int>("valid", nodeValid, "1 if finite node model outputs were produced, else 0");
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            std::vector<float> logits(num_tracks), probabilities(num_tracks);
            for (size_t i = 0; i < num_tracks; ++i) {
                logits[i] = nodeLogits[i * N_NODE_CLASSES + k];
                probabilities[i] = nodeProbabilities[i * N_NODE_CLASSES + k];
            }
            nodeTable->addColumn<float>("logit" + std::to_string(k), logits, "raw ONNX node class logit; -1 when invalid");
            nodeTable->addColumn<float>("prob" + std::to_string(k), probabilities, "six-class softmax probability; -1 when invalid");
        }
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
    const auto& PV0 = pvs->front();
    const auto& ttBuilder = iSetup.getData(theTTBToken);


    // retrieve information of tracks. Everything here mirrors TrackVertexVars.cc exactly
    // (same PV, same setBeamSpot call, same IPTools entry points, same dz-significance
    // formula) so that inference features match what the model was trained on.
    std::vector<float> trk_pt, trk_eta, trk_phi;
    std::vector<float> trk_dz, trk_dzSignificance;
    std::vector<float> trk_ip2dValue, trk_ip2dSignificance, trk_ip3dValue, trk_ip3dSignificance;
    std::vector<float> trk_p, trk_charge, trk_energy, trk_chi2;
    std::vector<float> trk_numberOfValidHits, trk_numberOfValidPixelHits, trk_numberOfValidStripHits;
    std::vector<float> trk_px, trk_py, trk_pz, trk_qdotp, trk_qoverp;
    std::vector<float> trk_isElectron, trk_isMuon;

    std::vector<float> track_SVscore;

    track_SVscore.reserve(num_tracks);

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

        trk_pt.push_back(trk.pt());
        trk_eta.push_back(trk.eta());
        trk_phi.push_back(trk.phi());
        trk_p.push_back(trk.p());
        trk_charge.push_back(trk.charge());
        trk_numberOfValidHits.push_back(trk.numberOfValidHits());
        trk_numberOfValidPixelHits.push_back(trk.hitPattern().numberOfValidPixelHits());
        trk_numberOfValidStripHits.push_back(trk.hitPattern().numberOfValidStripHits());
        trk_energy.push_back(std::sqrt(trk.p() * trk.p() + PION_MASS * PION_MASS));
        trk_chi2.push_back(trk.chi2());
        trk_px.push_back(trk.px());
        trk_py.push_back(trk.py());
        trk_pz.push_back(trk.pz());
        trk_qdotp.push_back(trk.charge() * trk.p());
        trk_qoverp.push_back(trk.qoverp());

        // dz / dz significance -- identical formula to TrackVertexVars.cc
        float dzVal = trk.dz(PV0.position());
        float dzErr = std::sqrt(std::pow(trk.dzError(), 2) + PV0.covariance(2, 2));
        trk_dz.push_back(dzVal);
        trk_dzSignificance.push_back((dzErr > 0.f && std::isfinite(dzErr)) ? dzVal / dzErr : -9999.f);

        // IP2D / IP3D -- same IPTools entry points as TrackVertexVars.cc, same
        // setBeamSpot call beforehand (this was missing before and changes the result)
        ttrk.setBeamSpot(*beamSpot);

        float ip2dVal = -9999.f, ip2dSig = -9999.f;
        std::pair<bool, Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(ttrk, PV0);
        if (ip2d.first) {
            if (edm::isFinite(ip2d.second.value())) ip2dVal = ip2d.second.value();
            if (edm::isFinite(ip2d.second.significance())) ip2dSig = ip2d.second.significance();
        }
        trk_ip2dValue.push_back(ip2dVal);
        trk_ip2dSignificance.push_back(ip2dSig);

        float ip3dVal = -9999.f, ip3dSig = -9999.f;
        std::pair<bool, Measurement1D> ip3d = IPTools::absoluteImpactParameter3D(ttrk, PV0);
        if (ip3d.first) {
            if (edm::isFinite(ip3d.second.value())) ip3dVal = ip3d.second.value();
            if (edm::isFinite(ip3d.second.significance())) ip3dSig = ip3d.second.significance();
        }
        trk_ip3dValue.push_back(ip3dVal);
        trk_ip3dSignificance.push_back(ip3dSig);

        node++;
    }

    // Track-lepton matching (isMuon / isElectron) -- same best-dR-within-relPt-window
    // logic as TrackVertexVars.cc, using the SAME dRMatchMax/relPtMatchMax values as
    // the training-time trackVertexVars configuration (pass 0.001 / 0.2 in the cfi,
    // not TrackVertexVars.cc's own fillDescriptions defaults of 0.01 / 0.1).
    trk_isMuon.assign(node, 0.f);
    trk_isElectron.assign(node, 0.f);
    for (int ni = 0; ni < node; ++ni) {
        int oi = nodeToOrig[ni];
        const auto& trk = (*tracks)[oi];
        float trkEta = trk.eta();
        float trkPhi = trk.phi();
        float trkPt = trk.pt();

        float bestDR_mu = std::numeric_limits<float>::infinity();
        for (const auto& mu : *muons) {
            float ptRatio = std::abs(trkPt - mu.pt()) / mu.pt();
            if (ptRatio > relPtMatchMax_) continue;
            float dr = reco::deltaR(trkEta, trkPhi, mu.eta(), mu.phi());
            if (dr >= dRMatchMax_ || dr >= bestDR_mu) continue;
            bestDR_mu = dr;
            trk_isMuon[ni] = 1.f;
        }

        float bestDR_e = std::numeric_limits<float>::infinity();
        for (const auto& ele : *electrons) {
            float ptRatio = std::abs(trkPt - ele.pt()) / ele.pt();
            if (ptRatio > relPtMatchMax_) continue;
            float dr = reco::deltaR(trkEta, trkPhi, ele.eta(), ele.phi());
            if (dr >= dRMatchMax_ || dr >= bestDR_e) continue;
            bestDR_e = dr;
            trk_isElectron[ni] = 1.f;
        }
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


    // === Assemble the exact raw 22 node features used by the established training path.
    // The ONNX graph embeds the training normalization; do not normalize here.
    std::vector<std::vector<float>> track_features;
    static constexpr int N_NODE_FEATURES = 22;
    static constexpr int N_EDGE_FEATURES = 14;

    for (int ni = 0; ni < node; ++ni) {
        std::vector<float> features = {
            trk_pt[ni],
            trk_eta[ni],
            trk_phi[ni],
            trk_dz[ni],
            trk_dzSignificance[ni],
            trk_ip2dSignificance[ni],
            trk_ip2dValue[ni],
            trk_ip3dSignificance[ni],
            trk_ip3dValue[ni],
            trk_numberOfValidHits[ni],
            trk_charge[ni],
            trk_energy[ni],
            trk_chi2[ni],
            trk_numberOfValidStripHits[ni],
            trk_numberOfValidPixelHits[ni],
            trk_px[ni],
            trk_py[ni],
            trk_pz[ni],
            trk_qdotp[ni],
            trk_qoverp[ni],
            trk_isElectron[ni],
            trk_isMuon[ni],
        };

        bool has_nan = false;
        for (float val : features) {
            if (!std::isfinite(val)) {
                has_nan = true;
                break;
            }
        }

        if (has_nan) {
            // Whole-track invalid sentinel. This is the -999.0 that GraphEventGNN.forward()
            // checks on feature[0] (track_pt) to mask the node with the model's learned
            // dummy token -- confirm this matches whatever DatasetConstructorTask pads
            // invalid/missing tracks with in cpf_arr, since it's a DIFFERENT convention
            // from the per-quantity -9999.f sentinel used above for individual IP/dz
            // failures on an otherwise-valid track.
            features.assign(N_NODE_FEATURES, -999.0f);
        }

        track_features.push_back(features);
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

         edge_features.push_back({
         pair_deltaEta[idx],
         pair_deltaPhi[idx],
         pair_pt[idx],
         pair_eta[idx],
         dca[idx],
         deltaR[idx],
         dca_sig[idx],
         cptopv[idx],
         pvtoPCA_i[idx],
         pvtoPCA_j[idx],
         dotprod_i[idx],
         dotprod_j[idx],
         pair_mom[idx],
         pair_invmass[idx]
         });
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
    if (output.size() != 3 || output[0].size() != track_features.size() * N_NODE_CLASSES ||
        output[1].size() != edge_i.size() || output[2].size() != track_features.size() * 128) {
        throw cms::Exception("InvalidModelOutput") << "TrackGNNSeederInferenceProducer got unexpected output tensor count or length";
    }

    const auto& rawNodeLogits = output[0];
    const auto& rawEdgeLogits = output[1];
    const auto& rawEmbeddings = output[2];
    auto softmaxNode = [&](size_t node) -> std::array<float, N_NODE_CLASSES> {
        std::array<float, N_NODE_CLASSES> probabilities{};
        float maximum = -std::numeric_limits<float>::infinity();
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            const float value = rawNodeLogits[node * N_NODE_CLASSES + k];
            if (!std::isfinite(value)) return probabilities;
            maximum = std::max(maximum, value);
        }
        float denominator = 0.f;
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            probabilities[k] = std::exp(rawNodeLogits[node * N_NODE_CLASSES + k] - maximum);
            denominator += probabilities[k];
        }
        if (!(denominator > 0.f) || !std::isfinite(denominator)) return {};
        for (float& value : probabilities) value /= denominator;
        return probabilities;
    };
    const auto sigmoid = [](float value) {
        if (value >= 0.f) return 1.f / (1.f + std::exp(-value));
        const float expValue = std::exp(value);
        return expValue / (1.f + expValue);
    };

    for (size_t ni = 0; ni < nodeToOrig.size(); ++ni) {
        const size_t original = static_cast<size_t>(nodeToOrig[ni]);
        const auto probabilities = softmaxNode(ni);
        bool valid = true;
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            if (!std::isfinite(rawNodeLogits[ni * N_NODE_CLASSES + k]) || !std::isfinite(probabilities[k])) valid = false;
        }
        bool embeddingIsValid = true;
        for (int k = 0; k < 128; ++k) {
            if (!std::isfinite(rawEmbeddings[ni * 128 + k])) embeddingIsValid = false;
        }
        if (valid) {
            nodeValid[original] = 1;
            for (int k = 0; k < N_NODE_CLASSES; ++k) {
                nodeLogits[original * N_NODE_CLASSES + k] = rawNodeLogits[ni * N_NODE_CLASSES + k];
                nodeProbabilities[original * N_NODE_CLASSES + k] = probabilities[k];
            }
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
