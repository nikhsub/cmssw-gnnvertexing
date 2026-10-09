#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Math/SMatrix.h"
#include "Math/SVector.h"
typedef ROOT::Math::SVector<double, 3> Vector3D;
typedef reco::Vertex::CovarianceMatrix CovMatrix;

// ---------------------------------------------------------------------------
// GVMatchProducer
//
// Consumes the plain vector<T> products published by GenVertexProducer
// (hadronGVx/y/z, hadronEta/Phi, daughterPt/Eta/Phi/Charge/GVidx, nGV) plus
// the reco-level collections (primary vertices, secondary vertices, tracks),
// and performs:
//   - GV <-> SV matching (distance-significance based, with track-overlap
//     confirmation), producing SVGVMatchTable / SVGVtrkMatchTable and a
//     "GV"-extension table (Hadron_SVIdx, SV_distanceSig, minDistNotMatched,
//     nDaughtersMatchedToTracks, maxDaughterTrkDeltaR)
//   - GVDaughters <-> reconstructed-track matching, producing a
//     "GVDaughters"-extension table (trkIdx, isTrkMatched, trkDeltaR, trkDPtRel)
//
// Row order for both extension tables matches GenVertexProducer's GVTable/
// GVDaughtersTable exactly, since they are built from the same vectors,
// untouched, in the same order.
// ---------------------------------------------------------------------------

class GVMatchProducer : public edm::stream::EDProducer<> {
public:
    explicit GVMatchProducer(const edm::ParameterSet&);
    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    std::vector<std::vector<float>> computeDistanceMatrix(const std::vector<float>& SV_x,
                                                            const std::vector<float>& SV_y,
                                                            const std::vector<float>& SV_z,
                                                            std::vector<CovMatrix> SV_cov,
                                                            const std::vector<float>& Hadron_GVx,
                                                            const std::vector<float>& Hadron_GVy,
                                                            const std::vector<float>& Hadron_GVz);
    void printDistanceMatrix(const std::vector<std::vector<float>>& distances);

    std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>,
               std::vector<int>>
    matchHadronsToSV(std::vector<std::vector<float>> distances, const std::vector<float>& SVtrk_pt,
                      const std::vector<float>& SVtrk_eta, const std::vector<float>& SVtrk_phi,
                      const std::vector<int>& SVtrk_SVidx, const std::vector<float>& Daughters_pt,
                      const std::vector<float>& Daughters_eta, const std::vector<float>& Daughters_phi,
                      const std::vector<int>& Daughters_GVidx, const std::vector<float>& SV_eta,
                      const std::vector<float>& SV_phi, const std::vector<float>& GV_eta,
                      const std::vector<float>& GV_phi, int n_Hadrons, int nRequiredCommonTracks, double dR_max,
                      double relPt_max, bool doubleMatching, int doubleMatching_nRequiredCommonTracks,
                      double doubleMatching_maxSignificance, double doubleMatching_dR_max,
                      double doubleMatching_relPt_max);

    // Matches gen daughters (by kinematics passed in) to reconstructed tracks.
    std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>> matchDaughtersToTracks(
        const std::vector<float>& Daughters_pt, const std::vector<float>& Daughters_eta,
        const std::vector<float>& Daughters_phi, const std::vector<int>& Daughters_charge,
        const std::vector<reco::Track>& tracks, double maxDeltaR, double maxDPtRel, bool checkCharge,
        bool resolveAmbiguities) const;

    const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
    edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;

    // Inputs from GenVertexProducer
    edm::EDGetTokenT<std::vector<float>> hadronGVxToken_, hadronGVyToken_, hadronGVzToken_;
    edm::EDGetTokenT<std::vector<float>> hadronEtaToken_, hadronPhiToken_;
    edm::EDGetTokenT<std::vector<float>> daughterPtToken_, daughterEtaToken_, daughterPhiToken_;
    edm::EDGetTokenT<std::vector<int>> daughterChargeToken_, daughterGVidxToken_;
    edm::EDGetTokenT<int> nGVToken_;

    int nRequiredCommonTracks_;
    double dlenSigMin_;
    double dR_max_;
    double relPt_max_;
    bool doubleMatching_;
    int doubleMatching_nRequiredCommonTracks_;
    double doubleMatching_maxSignificance_;
    double doubleMatching_dR_max_;
    double doubleMatching_relPt_max_;

    double trkMinPt_;
    double trkMaxEta_;
    double trkMaxDeltaR_;
    double trkMaxDPtRel_;
    bool trkCheckCharge_;
    bool trkResolveAmbiguities_;
};

GVMatchProducer::GVMatchProducer(const edm::ParameterSet& iConfig)
    : pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
      tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      hadronGVxToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("hadronGVx"))),
      hadronGVyToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("hadronGVy"))),
      hadronGVzToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("hadronGVz"))),
      hadronEtaToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("hadronEta"))),
      hadronPhiToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("hadronPhi"))),
      daughterPtToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("daughterPt"))),
      daughterEtaToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("daughterEta"))),
      daughterPhiToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("daughterPhi"))),
      daughterChargeToken_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("daughterCharge"))),
      daughterGVidxToken_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("daughterGVidx"))),
      nGVToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nGV"))),
      nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
      dR_max_(iConfig.getParameter<double>("dR_max")),
      relPt_max_(iConfig.getParameter<double>("relPt_max")),
      doubleMatching_(iConfig.getParameter<bool>("doubleMatching")),
      doubleMatching_nRequiredCommonTracks_(iConfig.getParameter<int>("doubleMatching_nRequiredCommonTracks")),
      doubleMatching_maxSignificance_(iConfig.getParameter<double>("doubleMatching_maxSignificance")),
      doubleMatching_dR_max_(iConfig.getParameter<double>("doubleMatching_dR_max")),
      doubleMatching_relPt_max_(iConfig.getParameter<double>("doubleMatching_relPt_max")),
      trkMinPt_(iConfig.getParameter<double>("trkMinPt")),
      trkMaxEta_(iConfig.getParameter<double>("trkMaxEta")),
      trkMaxDeltaR_(iConfig.getParameter<double>("trkMaxDeltaR")),
      trkMaxDPtRel_(iConfig.getParameter<double>("trkMaxDPtRel")),
      trkCheckCharge_(iConfig.getParameter<bool>("trkCheckCharge")),
      trkResolveAmbiguities_(iConfig.getParameter<bool>("trkResolveAmbiguities")) {
    produces<nanoaod::FlatTable>("SVGVMatchTable");
    produces<nanoaod::FlatTable>("SVGVtrkMatchTable");
    produces<nanoaod::FlatTable>("GVMatchExtTable");        // extension=true, table name "GV"
    produces<nanoaod::FlatTable>("GVDaughtersMatchTable");  // extension=true, table name "GVDaughters"
}

void GVMatchProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
    edm::Handle<std::vector<reco::Vertex>> svHandle;
    iEvent.getByToken(svToken_, svHandle);
    auto pvsIn = iEvent.getHandle(pvs_);
    edm::Handle<std::vector<reco::Track>> tracksHandle;
    iEvent.getByToken(tracksToken_, tracksHandle);

    const auto& secondaryVertices = svHandle;
    const auto& tracks = *tracksHandle;

    edm::Handle<std::vector<float>> hadronGVxH, hadronGVyH, hadronGVzH, hadronEtaH, hadronPhiH;
    iEvent.getByToken(hadronGVxToken_, hadronGVxH);
    iEvent.getByToken(hadronGVyToken_, hadronGVyH);
    iEvent.getByToken(hadronGVzToken_, hadronGVzH);
    iEvent.getByToken(hadronEtaToken_, hadronEtaH);
    iEvent.getByToken(hadronPhiToken_, hadronPhiH);

    edm::Handle<std::vector<float>> daughterPtH, daughterEtaH, daughterPhiH;
    edm::Handle<std::vector<int>> daughterChargeH, daughterGVidxH;
    iEvent.getByToken(daughterPtToken_, daughterPtH);
    iEvent.getByToken(daughterEtaToken_, daughterEtaH);
    iEvent.getByToken(daughterPhiToken_, daughterPhiH);
    iEvent.getByToken(daughterChargeToken_, daughterChargeH);
    iEvent.getByToken(daughterGVidxToken_, daughterGVidxH);

    edm::Handle<int> nGVH;
    iEvent.getByToken(nGVToken_, nGVH);
    const int ngv = *nGVH;

    const std::vector<float>& Hadron_GVx = *hadronGVxH;
    const std::vector<float>& Hadron_GVy = *hadronGVyH;
    const std::vector<float>& Hadron_GVz = *hadronGVzH;
    const std::vector<float>& Hadron_eta = *hadronEtaH;
    const std::vector<float>& Hadron_phi = *hadronPhiH;

    const std::vector<float>& Daughters_pt = *daughterPtH;
    const std::vector<float>& Daughters_eta = *daughterEtaH;
    const std::vector<float>& Daughters_phi = *daughterPhiH;
    const std::vector<int>& Daughters_charge = *daughterChargeH;
    const std::vector<int>& Daughters_GVidx = *daughterGVidxH;

    VertexDistance3D vdist;
    const auto& PV0 = pvsIn->front();

    // save coordinates of SV (will be used for matching with GV)
    std::vector<float> SV_x, SV_y, SV_z, SV_eta, SV_phi;
    std::vector<CovMatrix> SV_cov;
    for (auto const& sv : *secondaryVertices) {
        Measurement1D dl =
            vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
        if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
            SV_x.push_back(sv.x());
            SV_y.push_back(sv.y());
            SV_z.push_back(sv.z());

            // Get Eta and Phi from tracks
            TLorentzVector p4s_SV = TLorentzVector(0, 0, 0, 0);
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                TLorentzVector p4;
                p4.SetPtEtaPhiM(trkRef->pt(), trkRef->eta(), trkRef->phi(), 0.13957039);
                p4s_SV += p4;
            }
            SV_eta.push_back(p4s_SV.Eta());
            SV_phi.push_back(p4s_SV.Phi());
            SV_cov.push_back(sv.covariance());
        }
    }

    // Filling tracks from reco SV
    std::vector<float> SVtrk_pt, SVtrk_eta, SVtrk_phi;
    std::vector<int> SVtrk_SVidx;
    int SV_index = 0;
    for (const auto& sv : *secondaryVertices) {
        Measurement1D dl =
            vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
        if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                if (trkRef.isNull()) continue;
                SVtrk_pt.push_back(trkRef->pt());
                SVtrk_eta.push_back(trkRef->eta());
                SVtrk_phi.push_back(trkRef->phi());
                SVtrk_SVidx.push_back(SV_index);
            }
            SV_index++;
        }
    }

    // Compute matrix of distances between SV and GV
    auto distances = computeDistanceMatrix(SV_x, SV_y, SV_z, SV_cov, Hadron_GVx, Hadron_GVy, Hadron_GVz);

    std::vector<int> Hadron_SVIdx(ngv, -1);
    std::vector<float> Hadron_SVDistance(ngv, -1);
    std::vector<float> Hadron_minDistNotMatched(ngv, 999.f);

    // perform matching based on distance matrix and track-to-daughter matching
    auto result = matchHadronsToSV(distances, SVtrk_pt, SVtrk_eta, SVtrk_phi, SVtrk_SVidx, Daughters_pt, Daughters_eta,
                                    Daughters_phi, Daughters_GVidx, SV_eta, SV_phi, Hadron_eta, Hadron_phi, ngv,
                                    nRequiredCommonTracks_, dR_max_, relPt_max_, doubleMatching_,
                                    doubleMatching_nRequiredCommonTracks_, doubleMatching_maxSignificance_,
                                    doubleMatching_dR_max_, doubleMatching_relPt_max_);

    Hadron_SVIdx = std::get<0>(result);
    Hadron_SVDistance = std::get<1>(result);
    Hadron_minDistNotMatched = std::get<2>(result);
    std::vector<int> SVtrk_isMatched = std::get<3>(result);
    std::vector<int> SVtrk_GVIdx = std::get<4>(result);
    std::vector<int> SVtrk_daughterIdx = std::get<5>(result);

    auto svTrkGVTable = std::make_unique<nanoaod::FlatTable>(SVtrk_pt.size(), "mySVtrks", false, true);
    svTrkGVTable->addColumn<int>("isMatched", SVtrk_isMatched, "1 if track is matched to a genParticle daughter of the matched GV");
    svTrkGVTable->addColumn<int>("GVIdx", SVtrk_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
    svTrkGVTable->addColumn<int>("daughterIdx", SVtrk_daughterIdx,
                                 "Index into GVDaughters table of the matched gen daughter, -1 if unmatched");
    iEvent.put(std::move(svTrkGVTable), "SVGVtrkMatchTable");

    // Build per-SV matched flag / GV index (same ordering as SV_x, since same dlenSig cut)
    std::vector<int> SV_GVIdx(SV_x.size(), -1);
    std::vector<int> SV_isMatched(SV_x.size(), 0);
    for (int had = 0; had < ngv; ++had) {
        int sv = Hadron_SVIdx[had];
        if (sv >= 0) {
            SV_GVIdx[sv] = had;
            SV_isMatched[sv] = 1;
        }
    }

    auto svGVTable = std::make_unique<nanoaod::FlatTable>(SV_x.size(), "mySV", false, true);
    svGVTable->addColumn<int>("GVIdx", SV_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
    svGVTable->addColumn<int>("isMatched", SV_isMatched, "1 if SV matched to a GV");
    iEvent.put(std::move(svGVTable), "SVGVMatchTable");

    // match GVDaughters (gen daughters) to reconstructed tracks
    auto trkMatchResult = matchDaughtersToTracks(Daughters_pt, Daughters_eta, Daughters_phi, Daughters_charge, tracks,
                                                  trkMaxDeltaR_, trkMaxDPtRel_, trkCheckCharge_, trkResolveAmbiguities_);

    std::vector<int> Daughters_trkIdx = std::get<0>(trkMatchResult);
    std::vector<int> Daughters_isTrkMatched = std::get<1>(trkMatchResult);
    std::vector<float> Daughters_trkDeltaR = std::get<2>(trkMatchResult);
    std::vector<float> Daughters_trkDPtRel = std::get<3>(trkMatchResult);

    // per-GV count of daughters matched to a track
    std::vector<int> GV_nDaughtersMatchedToTracks(ngv, 0);
    for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
        if (!Daughters_isTrkMatched[i]) continue;
        int gvIdx = Daughters_GVidx[i];
        if (gvIdx >= 0 && gvIdx < ngv) {
            GV_nDaughtersMatchedToTracks[gvIdx]++;
        }
    }

    // per-GV max trkDeltaR among its (matched) daughters
    std::vector<float> GV_maxDaughterTrkDeltaR(ngv, -1.f);
    for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
        if (!Daughters_isTrkMatched[i]) continue;
        int gvIdx = Daughters_GVidx[i];
        if (gvIdx >= 0 && gvIdx < ngv) {
            if (Daughters_trkDeltaR[i] > GV_maxDaughterTrkDeltaR[gvIdx]) {
                GV_maxDaughterTrkDeltaR[gvIdx] = Daughters_trkDeltaR[i];
            }
        }
    }

    // "GV"-extension table: matching-derived per-GV columns
    auto gvExtTable = std::make_unique<nanoaod::FlatTable>(ngv, "GV", false, true);
    gvExtTable->addColumn<int>("Hadron_SVIdx", Hadron_SVIdx, "SVIdx");
    gvExtTable->addColumn<float>("SV_distanceSig", Hadron_SVDistance, "SV_distanceSig");
    gvExtTable->addColumn<float>("minDistNotMatched", Hadron_minDistNotMatched, "Minimum distance to SV among unmatched hadrons");
    gvExtTable->addColumn<int>("nDaughtersMatchedToTracks", GV_nDaughtersMatchedToTracks,
                              "Number of GVDaughters of this GV matched to a reconstructed track");
    gvExtTable->addColumn<float>("maxDaughterTrkDeltaR", GV_maxDaughterTrkDeltaR,
                                 "Max trkDeltaR among this GV's track-matched GVDaughters (not GVDirectDaughters); -1 if none matched");
    iEvent.put(std::move(gvExtTable), "GVMatchExtTable");

    // "GVDaughters"-extension table: track-matching columns per daughter
    auto dauExtTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(), "GVDaughters", false, true);
    dauExtTable->addColumn<int>("trkIdx", Daughters_trkIdx, "Index of matched track in the input track collection, -1 if unmatched");
    dauExtTable->addColumn<int>("isTrkMatched", Daughters_isTrkMatched, "1 if daughter matched to a reconstructed track");
    dauExtTable->addColumn<float>("trkDeltaR", Daughters_trkDeltaR, "deltaR to matched track, -1 if unmatched");
    dauExtTable->addColumn<float>("trkDPtRel", Daughters_trkDPtRel, "relative pT difference to matched track, -1 if unmatched");
    iEvent.put(std::move(dauExtTable), "GVDaughtersMatchTable");
}

std::vector<std::vector<float>> GVMatchProducer::computeDistanceMatrix(const std::vector<float>& SV_x,
                                                                        const std::vector<float>& SV_y,
                                                                        const std::vector<float>& SV_z,
                                                                        std::vector<CovMatrix> SV_cov,
                                                                        const std::vector<float>& Hadron_GVx,
                                                                        const std::vector<float>& Hadron_GVy,
                                                                        const std::vector<float>& Hadron_GVz) {
    // Returns matrix of "chi2-like" distances between SV and GV (nSV x nHadron)
    size_t nSV = SV_x.size();
    size_t nHadron = Hadron_GVx.size();

    std::vector<std::vector<float>> distances(nSV, std::vector<float>(nHadron, 999.0));

    for (size_t i = 0; i < nSV; ++i) {
        CovMatrix covInv = SV_cov[i];
        covInv.Invert();
        for (size_t j = 0; j < nHadron; ++j) {
            float dx = SV_x[i] - Hadron_GVx[j];
            float dy = SV_y[i] - Hadron_GVy[j];
            float dz = SV_z[i] - Hadron_GVz[j];
            float chi2 = dx * (covInv(0, 0) * dx + covInv(0, 1) * dy + covInv(0, 2) * dz) +
                         dy * (covInv(1, 0) * dx + covInv(1, 1) * dy + covInv(1, 2) * dz) +
                         dz * (covInv(2, 0) * dx + covInv(2, 1) * dy + covInv(2, 2) * dz);

            float dist = std::sqrt(chi2);
            distances[i][j] = dist;
        }
    }

    return distances;
}

void GVMatchProducer::printDistanceMatrix(const std::vector<std::vector<float>>& distances) {
    size_t nSV = distances.size();
    if (nSV == 0) return;

    size_t nGV = distances[0].size();

    std::cout << "\n Distance Matrix (SV rows x GV cols):\n\n";

    std::cout << std::setw(8) << "SV/GV";
    for (size_t j = 0; j < nGV; ++j) {
        std::cout << std::setw(10) << "GV[" + std::to_string(j) + "]";
    }
    std::cout << "\n";

    for (size_t i = 0; i < nSV; ++i) {
        std::cout << std::setw(8) << "SV[" + std::to_string(i) + "]";
        for (size_t j = 0; j < nGV; ++j) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(3) << distances[i][j];
        }
        std::cout << "\n";
    }
}

static float computeDR_SV_Had(int bestSV, int bestHad, const std::vector<float>& SV_eta,
                               const std::vector<float>& SV_phi, const std::vector<float>& GV_eta,
                               const std::vector<float>& GV_phi) {
    float dEta = SV_eta[bestSV] - GV_eta[bestHad];
    float dPhi = SV_phi[bestSV] - GV_phi[bestHad];

    while (dPhi > M_PI) dPhi -= 2.0 * M_PI;
    while (dPhi < -M_PI) dPhi += 2.0 * M_PI;

    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>,
           std::vector<int>>
GVMatchProducer::matchHadronsToSV(std::vector<std::vector<float>> distances, const std::vector<float>& SVtrk_pt,
                                   const std::vector<float>& SVtrk_eta, const std::vector<float>& SVtrk_phi,
                                   const std::vector<int>& SVtrk_SVidx, const std::vector<float>& Daughters_pt,
                                   const std::vector<float>& Daughters_eta, const std::vector<float>& Daughters_phi,
                                   const std::vector<int>& Daughters_GVidx, const std::vector<float>& SV_eta,
                                   const std::vector<float>& SV_phi, const std::vector<float>& GV_eta,
                                   const std::vector<float>& GV_phi, int n_Hadrons, int nRequiredCommonTracks,
                                   double dR_max, double relPt_max, bool doubleMatching,
                                   int doubleMatching_nRequiredCommonTracks, double doubleMatching_maxSignificance,
                                   double doubleMatching_dR_max, double doubleMatching_relPt_max) {
    size_t nSV = distances.size();
    std::vector<int> Hadron_SVIdx(n_Hadrons, -1);
    std::vector<float> Hadron_SVDistance(n_Hadrons, -1);
    std::vector<int> SVtrk_isMatched(SVtrk_pt.size(), 0);
    std::vector<int> SVtrk_GVIdx(SVtrk_pt.size(), -1);
    std::vector<int> SVtrk_daughterIdx(SVtrk_pt.size(), -1);
    std::vector<size_t> svTrackIdxs_fromBestSV;
    std::vector<std::vector<float>> distancesOriginal = distances;

    while (true) {
        float minDist = 999.0;
        int bestSV = -1;
        int bestHad = -1;

        for (size_t sv = 0; sv < nSV; ++sv) {
            for (int had = 0; had < n_Hadrons; ++had) {
                if (distances[sv][had] < minDist) {
                    minDist = distances[sv][had];
                    bestSV = sv;
                    bestHad = had;
                }
            }
        }
        if (minDist >= 997.0) break;  // done

        svTrackIdxs_fromBestSV.clear();
        for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
            if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > trkMinPt_ && std::fabs(SVtrk_eta[i]) < trkMaxEta_) {
                svTrackIdxs_fromBestSV.push_back(i);
            }
        }

        std::vector<size_t> GenDaughtersIdxs_fromBestHad;
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (Daughters_GVidx[i] == bestHad) {
                GenDaughtersIdxs_fromBestHad.push_back(i);
            }
        }

        int common = 0;
        std::vector<std::pair<size_t, size_t>> matchedTrackToDaughter;
        for (size_t iSV : svTrackIdxs_fromBestSV) {
            bool trackMatched = false;
            size_t matchedDauIdx = 0;
            for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / Daughters_pt[iHad];
                if (dR < dR_max && relPt < relPt_max) {
                    trackMatched = true;
                    matchedDauIdx = iHad;
                    break;
                }
            }
            if (trackMatched) {
                ++common;
                matchedTrackToDaughter.push_back({iSV, matchedDauIdx});
            }
        }

        if (common >= nRequiredCommonTracks) {
            Hadron_SVIdx[bestHad] = bestSV;
            Hadron_SVDistance[bestHad] = minDist;
            for (auto& [iSV, iHad] : matchedTrackToDaughter) {
                SVtrk_isMatched[iSV] = 1;
                SVtrk_GVIdx[iSV] = bestHad;
                SVtrk_daughterIdx[iSV] = iHad;
            }
            for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 1000.0;
            for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 1000.0;
        } else {
            distances[bestSV][bestHad] = 998.0;
        }
    }

    if (doubleMatching) {
        while (true) {
            float minDist = 999.0;
            int bestSV = -1;
            int bestHad = -1;
            for (size_t sv = 0; sv < nSV; ++sv) {
                for (int had = 0; had < n_Hadrons; ++had) {
                    if (distances[sv][had] < minDist) {
                        minDist = distances[sv][had];
                        bestSV = sv;
                        bestHad = had;
                    }
                }
            }

            if (minDist >= 998.5) break;  // done
            float dR_SV_Had = computeDR_SV_Had(bestSV, bestHad, SV_eta, SV_phi, GV_eta, GV_phi);
            (void)dR_SV_Had;  // computed for potential use in the cut below (kept for parity with original)

            svTrackIdxs_fromBestSV.clear();
            for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
                if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > trkMinPt_ && std::fabs(SVtrk_eta[i]) < trkMaxEta_) {
                    svTrackIdxs_fromBestSV.push_back(i);
                }
            }

            std::vector<size_t> GenDaughtersIdxs_fromBestHad;
            for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
                if (Daughters_GVidx[i] == bestHad) {
                    GenDaughtersIdxs_fromBestHad.push_back(i);
                }
            }

            int common = 0;
            for (size_t iSV : svTrackIdxs_fromBestSV) {
                for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                    float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                    float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / Daughters_pt[iHad];
                    if (dR < doubleMatching_dR_max && relPt < doubleMatching_relPt_max) {
                        ++common;
                        if (common >= 1) break;
                    }
                }
                if (common >= doubleMatching_nRequiredCommonTracks) break;
            }

            if (common >= doubleMatching_nRequiredCommonTracks &&
                distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance) {
                Hadron_SVIdx[bestHad] = bestSV;
                Hadron_SVDistance[bestHad] = -distancesOriginal[bestSV][bestHad];
                for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 1000.0;
                for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 1000.0;
            } else {
                distances[bestSV][bestHad] = 999.0;
            }
        }
    }

    // check the minDist not matched using original matrix (not modified by matching procedure)
    std::vector<float> minDistNotMatched(n_Hadrons, 999.);
    for (int had = 0; had < n_Hadrons; ++had) {
        for (size_t sv = 0; sv < nSV; ++sv) {
            if (distancesOriginal[sv][had] < minDistNotMatched[had] && Hadron_SVIdx[had] == -1) {
                minDistNotMatched[had] = distancesOriginal[sv][had];
            }
        }
    }
    return std::make_tuple(Hadron_SVIdx, Hadron_SVDistance, minDistNotMatched, SVtrk_isMatched, SVtrk_GVIdx,
                            SVtrk_daughterIdx);
}

// Matches gen-level GVDaughters to reconstructed tracks by deltaR and relative pT,
// following the same TrackGenMatcher-style logic (optional charge check, optional
// one-to-one "resolveAmbiguities" assignment).
std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>>
GVMatchProducer::matchDaughtersToTracks(const std::vector<float>& Daughters_pt, const std::vector<float>& Daughters_eta,
                                        const std::vector<float>& Daughters_phi, const std::vector<int>& Daughters_charge,
                                        const std::vector<reco::Track>& tracks, double maxDeltaR, double maxDPtRel,
                                        bool checkCharge, bool resolveAmbiguities) const {
    size_t nDau = Daughters_pt.size();
    size_t nTrk = tracks.size();

    std::vector<int> trkIdx(nDau, -1);
    std::vector<int> isMatched(nDau, 0);
    std::vector<float> matchDeltaR(nDau, -1.f);
    std::vector<float> matchDPtRel(nDau, -1.f);

    struct Pair {
        float dR;
        size_t dau;
        size_t trk;
    };
    std::vector<Pair> pairs;
    pairs.reserve(nDau);

    for (size_t i = 0; i < nDau; ++i) {
        for (size_t j = 0; j < nTrk; ++j) {
            if (checkCharge && tracks[j].charge() != Daughters_charge[i]) continue;

            float dR = deltaR(Daughters_eta[i], Daughters_phi[i], tracks[j].eta(), tracks[j].phi());
            if (dR >= maxDeltaR) continue;

            float dPtRel = std::fabs(tracks[j].pt() - Daughters_pt[i]) / Daughters_pt[i];
            if (dPtRel >= maxDPtRel) continue;

            pairs.push_back({dR, i, j});
        }
    }

    std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.dR < b.dR; });

    std::vector<bool> dauUsed(nDau, false);
    std::vector<bool> trkUsed(nTrk, false);
    for (const auto& p : pairs) {
        if (dauUsed[p.dau]) continue;
        if (resolveAmbiguities && trkUsed[p.trk]) continue;

        trkIdx[p.dau] = static_cast<int>(p.trk);
        isMatched[p.dau] = 1;
        matchDeltaR[p.dau] = p.dR;
        matchDPtRel[p.dau] = std::fabs(tracks[p.trk].pt() - Daughters_pt[p.dau]) / Daughters_pt[p.dau];

        dauUsed[p.dau] = true;
        if (resolveAmbiguities) trkUsed[p.trk] = true;
    }

    return std::make_tuple(trkIdx, isMatched, matchDeltaR, matchDPtRel);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GVMatchProducer);
