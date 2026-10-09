#include "PhysicsTools/SimpleVertexTable/interface/trackGVFeaturesProducer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include <limits>
#include <cmath>
#include <algorithm>

TrackGVFeaturesProducer::TrackGVFeaturesProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      //svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
      gvTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvTable"))),
      gvDauTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvDaughtersTable"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      trkPtCut_(iConfig.getParameter<double>("trkPtCut")),
      dRMatchMax_(iConfig.getParameter<double>("dRMatchMax")),
      relPtMatchMax_(iConfig.getParameter<double>("relPtMatchMax"))
{
    produces<nanoaod::FlatTable>("track");      // extension of trackTable
    produces<nanoaod::FlatTable>("trackPair");  // standalone edge table
    //produces<nanoaod::FlatTable>("GV");         // extension of GenVertexProducer's GV table
}

void TrackGVFeaturesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<std::vector<reco::Track>> tracksH;
    iEvent.getByToken(tracksToken_, tracksH);
    const auto& tracks = *tracksH;

    edm::Handle<std::vector<reco::Vertex>> pvH;
    iEvent.getByToken(pvToken_, pvH);
    const reco::Vertex& pv = pvH->front();


    edm::Handle<nanoaod::FlatTable> gvTable;
    iEvent.getByToken(gvTableToken_, gvTable);
    edm::Handle<nanoaod::FlatTable> gvDauTable;
    iEvent.getByToken(gvDauTableToken_, gvDauTable);

    const auto& theB = iSetup.getData(ttbToken_);

    // ---- Read GV / GVDaughters columns straight out of the upstream FlatTables ----
    //size_t nGV = gvTable->size();
    auto GV_isB = gvTable->columnData<int>(gvTable->columnIndex("isB"));
    auto GV_isD = gvTable->columnData<int>(gvTable->columnIndex("isD"));
    auto GV_fromHF = gvTable->columnData<int>(gvTable->columnIndex("fromHF"));
    //auto GV_pdgClass = gvTable->columnData<int>(gvTable->columnIndex("pdgClass"));

    size_t nDau = gvDauTable->size();
    auto Dau_pt = gvDauTable->columnData<float>(gvDauTable->columnIndex("pt"));
    auto Dau_eta = gvDauTable->columnData<float>(gvDauTable->columnIndex("eta"));
    auto Dau_phi = gvDauTable->columnData<float>(gvDauTable->columnIndex("phi"));
    auto Dau_hadIdx = gvDauTable->columnData<int>(gvDauTable->columnIndex("hadronIndex"));


    // Determine which track is from PV and which is from other PVs
    std::vector<int> trk_vtxIdx(tracks.size(), -1);
    for (size_t ipv = 0; ipv < pvH->size(); ++ipv) {
        const reco::Vertex& vtx = (*pvH)[ipv];
        for (auto it = vtx.tracks_begin(); it != vtx.tracks_end(); ++it) {
            const reco::TrackBaseRef& trkRef = *it;
            if (trkRef.isNull()) continue;
            size_t idx = trkRef.key();   // index into the same track collection `tracks` was built from
            if (idx < trk_vtxIdx.size()) {
                trk_vtxIdx[idx] = static_cast<int>(ipv);
            }
        }
    }

    // Track - Daughters Matching
    size_t nTrk = tracks.size();
    std::vector<int>   trk_hadidx(nTrk, -1);
    std::vector<int>   trk_label(nTrk, 6);  
    std::vector<float> trk_delr(nTrk, std::numeric_limits<float>::infinity());
    std::vector<float> trk_ptrat(nTrk, std::numeric_limits<float>::quiet_NaN());
    // match trks to genParticles
    // match is the best dR provided that dPt/pT< ptRatio_cut
    for (size_t it = 0; it < nTrk; ++it) {
        //if (tracks[it].pt() < trkPtCut_) continue;
        int bestDau = -1;
	if (tracks[it].pt() >= trkPtCut_) {
        	//float bestPtRatio = std::numeric_limits<float>::quiet_NaN();
        	float bestDR = std::numeric_limits<float>::infinity();
        	for (size_t id = 0; id < nDau; ++id) {
        	    if (Dau_pt[id] <= 0.f) continue;
        	    float ptRatio = std::abs(tracks[it].pt() - Dau_pt[id])/ Dau_pt[id];
        	    if (ptRatio > relPtMatchMax_) continue;
        	    float dr = reco::deltaR(tracks[it].eta(), tracks[it].phi(), Dau_eta[id], Dau_phi[id]);
        	    if (dr >= dRMatchMax_ || dr >= bestDR) continue;
        	    bestDR = dr; 
        	    bestDau = static_cast<int>(id); //bestPtRatio = ptRatio;
        	}
	}
        if (bestDau >= 0) {
            // Matched to a daughter
            trk_hadidx[it] = Dau_hadIdx[bestDau];
	    // Track truth labels:
	    //
	    // 0 : Primary     - associated with reconstructed PV0
	    // 1 : Pileup      - associated with reconstructed PV index > 0
	    // 2 : FromB       - daughter of B hadron
	    // 3 : FromBC      - daughter in B -> C / heavy-flavor chain
	    // 4 : FromC       - daughter of prompt C hadron
	    // 5 : OtherSec    - other truth-matched secondary
	    // 6 : Unassociated - not associated with any reconstructed PV

            if (trk_hadidx[it] == -1) {
                trk_label[it] = 0;       // matched to a daughter, but daughter has no hadron? impossible 
                std::cout << "  -> label=0 (hadidx == -1, daughter has no hadron)" << std::endl;
            } else {
                if (GV_isB[trk_hadidx[it]] == 1) {
                    trk_label[it] = 2;       // fromB
                } else if (GV_fromHF[trk_hadidx[it]] == 1) {
                    trk_label[it] = 3;       // fromBC (covers isD&&fromHF, and fromHF alone)
                } else if (GV_isD[trk_hadidx[it]] == 1) {
                    trk_label[it] = 4;       // fromC
                } else {
                    trk_label[it] = 5;       // OtherSecondary
                }
            }

        } else {
		// Track is not truth-matched to any GV daughter.
	       // Classify using reconstructed primary-vertex association.
	   
	       if (trk_vtxIdx[it] == 0) {
	           trk_label[it] = 0;   // Primary / PV0
	       }
	       else if (trk_vtxIdx[it] > 0) {
	           trk_label[it] = 1;   // Pileup / non-PV0
	       }
	       else {
	           trk_label[it] = 6;   // Unassociated with any reconstructed PV
	       }
	   }
        

    }

    auto trackTableOut = std::make_unique<nanoaod::FlatTable>(nTrk, "track", false, true);  // extension=true
    trackTableOut->addColumn<int>("hadidx", trk_hadidx, "matched GV daughter's hadron index, -1 if unmatched");
    trackTableOut->addColumn<int>("label", trk_label, "truth label (currently mirrors flav; adjust as needed)");
    iEvent.put(std::move(trackTableOut), "track");

    // ---- Track-pair (edge) features ----
    std::vector<int> trk1Idx, trk2Idx;
    std::vector<float> edgeDeltaR, edgeDeltaEta, edgeDeltaPhi, edgeDca, edgeDcaSig, edgeCptopv, edgePvtoPCA1, edgePvtoPCA2,
                       edgeDotprod1, edgeDotprod2, edgePairMom, edgePairInvMass, edgeLabel, edgePair_pt, edgePair_eta;

    std::vector<reco::TransientTrack> t_trks(nTrk);
    for (size_t i = 0; i < nTrk; ++i) t_trks[i] = theB.build(tracks[i]);

    const float PION_MASS = 0.13957018f;
    VertexDistance3D vdist3D;

    for (size_t i = 0; i < nTrk; ++i) {
        if (!t_trks[i].isValid()) continue;
        for (size_t j = i + 1; j < nTrk; ++j) {
            if (!t_trks[j].isValid()) continue;

            float dr = reco::deltaR(tracks[i].eta(), tracks[i].phi(), tracks[j].eta(), tracks[j].phi());
            float dEta = abs(tracks[i].eta() - tracks[j].eta());
            float dPhi = reco::deltaPhi(tracks[i].phi(), tracks[j].phi());
            if (dr>1.5) continue;
            //if (dEta>0.4) continue;
            //if (abs(dPhi)>0.5) continue;

            float e1 = std::sqrt(tracks[i].p() * tracks[i].p() + PION_MASS * PION_MASS);
            float e2 = std::sqrt(tracks[j].p() * tracks[j].p() + PION_MASS * PION_MASS);
            float sumPx = tracks[i].px() + tracks[j].px();
            float sumPy = tracks[i].py() + tracks[j].py();
            float pt_t1t2 = std::sqrt(sumPx*sumPx + sumPy*sumPy);
            float sumPz = tracks[i].pz() + tracks[j].pz();
            float compositeP   = std::sqrt(sumPx*sumPx + sumPy*sumPy + sumPz*sumPz);
            float eta_t1t2 = 0.f;
            if (pt_t1t2 < 1e-6f) {
                eta_t1t2 = (sumPz > 0) ? 9999.f
                                            : -9999.f;
            } else {
                eta_t1t2 = 0.5f * std::log((compositeP + sumPz) / (compositeP - sumPz));
            }
            //if (abs(eta_t1t2) > 2.5) continue;
            float sumE  = e1 + e2;
            float invMass = std::sqrt(std::max(0.f, sumE * sumE - (sumPx*sumPx + sumPy*sumPy + sumPz*sumPz)));
            if (invMass>5.0) continue;

            float dcaVal = -1.f, dcaSig = -1.f, cptopvVal = -1.f;
            float pvToPCA1 = -1.f, pvToPCA2 = -1.f, dot1 = -999.f, dot2 = -999.f, pairMomMag = -1.f;

            TwoTrackMinimumDistance minDist;
            if (minDist.calculate(t_trks[i].impactPointState(), t_trks[j].impactPointState())) {
                auto m = vdist3D.distance(
                    VertexState(minDist.points().second, t_trks[i].impactPointState().cartesianError().position()),
                    VertexState(minDist.points().first,  t_trks[j].impactPointState().cartesianError().position()));
                dcaVal = m.value();
                if (m.error() > 0) dcaSig = m.value() / m.error();

                GlobalPoint cp(minDist.crossingPoint());
                GlobalPoint pvp(pv.x(), pv.y(), pv.z());
                GlobalPoint seedPCA  = minDist.points().second;
                GlobalPoint trackPCA = minDist.points().first;

                pvToPCA1 = (seedPCA - pvp).mag();
                pvToPCA2 = (trackPCA - pvp).mag();
                cptopvVal = (cp - pvp).mag();
                dot2 = (trackPCA - pvp).unit().dot(t_trks[j].impactPointState().globalDirection().unit());
                dot1 = (seedPCA - pvp).unit().dot(t_trks[i].impactPointState().globalDirection().unit());

                GlobalVector pairMomentum((Basic3DVector<float>)(t_trks[i].track().momentum() + t_trks[j].track().momentum()));
                pairMomMag = pairMomentum.mag();
            }
            if (dcaVal>0.1) continue;

            trk1Idx.push_back(static_cast<int>(i));
            trk2Idx.push_back(static_cast<int>(j));
            edgeDeltaR.push_back(dr);
            edgeDeltaEta.push_back(dEta);
            edgeDeltaPhi.push_back(dPhi);
            edgeDca.push_back(dcaVal);
            edgeDcaSig.push_back(dcaSig);
            edgeCptopv.push_back(cptopvVal);
            edgePvtoPCA1.push_back(pvToPCA1);
            edgePvtoPCA2.push_back(pvToPCA2);
            edgeDotprod1.push_back(dot1);
            edgeDotprod2.push_back(dot2);
            edgePairMom.push_back(pairMomMag);
            edgePairInvMass.push_back(invMass);
            edgePair_pt.push_back(pt_t1t2);
            edgePair_eta.push_back(eta_t1t2);

            bool sameHadron = trk_hadidx[i] >= 0 && trk_hadidx[i] == trk_hadidx[j];
            edgeLabel.push_back(sameHadron ? 1.f : 0.f);
        }
    }

    auto pairTable = std::make_unique<nanoaod::FlatTable>(trk1Idx.size(), "trackPair", false);  // standalone
    pairTable->addColumn<int>("trk1Idx", trk1Idx, "index of first track in the pair");
    pairTable->addColumn<int>("trk2Idx", trk2Idx, "index of second track in the pair");
    pairTable->addColumn<float>("deltaR", edgeDeltaR, "deltaR between the two tracks");
    pairTable->addColumn<float>("deltaPhi", edgeDeltaPhi, "deltaPhi between the two tracks");
    pairTable->addColumn<float>("deltaEta", edgeDeltaEta, "deltaEta between the two tracks");
    pairTable->addColumn<float>("dca", edgeDca, "distance of closest approach");
    pairTable->addColumn<float>("dcaSig", edgeDcaSig, "DCA significance");
    pairTable->addColumn<float>("cptopv", edgeCptopv, "distance from crossing point to PV");
    pairTable->addColumn<float>("pvtoPCA1", edgePvtoPCA1, "distance PV to track1 PCA");
    pairTable->addColumn<float>("pvtoPCA2", edgePvtoPCA2, "distance PV to track2 PCA");
    pairTable->addColumn<float>("dotprod1", edgeDotprod1, "direction dot product, track1 side");
    pairTable->addColumn<float>("dotprod2", edgeDotprod2, "direction dot product, track2 side");
    pairTable->addColumn<float>("pairMom", edgePairMom, "magnitude of the pair momentum sum");
    pairTable->addColumn<float>("pairInvMass", edgePairInvMass, "invariant mass, pion mass hypothesis");
    pairTable->addColumn<float>("pt", edgePair_pt, "transverse momentum of the pair");
    pairTable->addColumn<float>("eta", edgePair_eta, "eta of the pair");
    pairTable->addColumn<float>("edgeLabel", edgeLabel, "1 if both tracks truth-matched to the same GV, else 0");
    iEvent.put(std::move(pairTable), "trackPair");


}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackGVFeaturesProducer);
