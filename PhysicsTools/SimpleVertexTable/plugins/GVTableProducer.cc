#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <unordered_map>
#include "TLorentzVector.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <optional>
#include <numeric>
#include <memory>
#include <iostream>
#include <iomanip>

#include "Math/SMatrix.h"
#include "Math/SVector.h"
//#include "Math/SMatrixFunctions.h"
//typedef ROOT::Math::SMatrix<float, 3, 3> Matrix3x3;
//typedef ROOT::Math::SVector<float, 3> Vector3;
typedef ROOT::Math::SVector<double,3> Vector3D;
typedef reco::Vertex::CovarianceMatrix CovMatrix;
class GenVertexProducer : public edm::stream::EDProducer<> {
public:
    explicit GenVertexProducer(const edm::ParameterSet&);
    void produce(edm::Event&, const edm::EventSetup&) override;

private:

    int checkPDG(int abs_pdg) const;

    std::optional<std::tuple<float, float, float>>isAncestor(const reco::Candidate* mother,const reco::Candidate* daughter) const;
    bool hasHFAncestor(const reco::Candidate* hadron) const;
    bool hasHFDescendant(const reco::Candidate* hadron) const;
    bool hasBHadronAncestor(const reco::Candidate* hadron) const;
    int getDaughterOriginLabelNoPU(const reco::Candidate* daughter) const;

    std::vector<std::vector<float>> computeDistanceMatrix(
                    const std::vector<float>& SV_x,const std::vector<float>& SV_y,const std::vector<float>& SV_z,
                    std::vector<CovMatrix> SV_cov,
                    const std::vector<float>& Hadron_GVx,const std::vector<float>& Hadron_GVy,const std::vector<float>& Hadron_GVz);
    void printDistanceMatrix(const std::vector<std::vector<float>>& distances);
    std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>, std::vector<int>>  matchHadronsToSV(
                                                                        std::vector<std::vector<float>> distances,
                                                                        const std::vector<float>& SVtrk_pt, const std::vector<float>& SVtrk_eta, const std::vector<float>& SVtrk_phi, const std::vector<int>& SVtrk_SVidx,
                                                                        const std::vector<float>& Daughters_pt,     //genparticles
                                                                        const std::vector<float>& Daughters_eta,  //genparticles
                                                                        const std::vector<float>& Daughters_phi,  //genparticles
                                                                        const std::vector<int>& Daughters_GVidx, // hadron index per daughter
                                                                        const std::vector<float>& SV_eta,
                                                                        const std::vector<float>& SV_phi,
                                                                        const std::vector<float>& GV_eta,
                                                                        const std::vector<float>& GV_phi,
                                                                        int n_Hadrons,
                                                                        int nRequiredCommonTracks,
                                                                        double dR_max,
                                                                        double relPt_max,
                                                                        bool doubleMatching,
                                                                        int doubleMatching_nRequiredCommonTracks,
                                                                        double doubleMatching_maxSignificance,
                                                                        double doubleMatching_dR_max,
                                                                        double doubleMatching_relPt_max
                                                                    );

    // NEW: matches GVDaughters (gen daughters) to reconstructed tracks.
    // Returns, indexed the same way as the Daughters_* vectors:
    //   trkIdx      : index into the `tracks` collection of the matched track, -1 if unmatched
    //   isMatched   : 1 if matched, 0 otherwise
    //   matchDeltaR : deltaR to the matched track, -1 if unmatched
    //   matchDPtRel : |pt_trk - pt_dau| / pt_dau of the matched track, -1 if unmatched
    std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>> matchDaughtersToTracks(
                                                                        const std::vector<float>& Daughters_pt,
                                                                        const std::vector<float>& Daughters_eta,
                                                                        const std::vector<float>& Daughters_phi,
                                                                        const std::vector<int>& Daughters_charge,
                                                                        const std::vector<reco::Track>& tracks,
                                                                        double maxDeltaR,
                                                                        double maxDPtRel,
                                                                        bool checkCharge,
                                                                        bool resolveAmbiguities) const;
	
    int findMatchingPrunedHadron(
	    const reco::Candidate* mergedHadron,
	    const reco::GenParticleCollection& prunedParticles
	) const;

    const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
    edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;  // NEW: general tracks used for GVDaughters matching
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
    int nRequiredCommonTracks_;
    double hadPt_min_;
    double dlenSigMin_;
    double dR_max_;
    double relPt_max_;
    bool doubleMatching_;
    int doubleMatching_nRequiredCommonTracks_;
    double doubleMatching_maxSignificance_;
    double doubleMatching_dR_max_;
    double doubleMatching_relPt_max_;

    // NEW: GVDaughters <-> tracks matching configuration
    double trkMaxDeltaR_;
    double trkMaxDPtRel_;
    bool trkCheckCharge_;
    bool trkResolveAmbiguities_;
};


GenVertexProducer::GenVertexProducer(const edm::ParameterSet& iConfig):
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
        genToken_(
        consumes<edm::View<reco::Candidate>>(
            iConfig.getParameter<edm::InputTag>("genParticles")
        )
    ),
    
    prunedGenToken_(
        consumes<reco::GenParticleCollection>(
            iConfig.getParameter<edm::InputTag>("prunedGenParticles")
        )
    ),
    svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
    tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
    jetFlavourInfosToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetFlavourInfos"))),
    nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
    hadPt_min_(iConfig.getParameter<double>("hadPt_min")),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    dR_max_(iConfig.getParameter<double>("dR_max")),
    relPt_max_(iConfig.getParameter<double>("relPt_max")),
    doubleMatching_(iConfig.getParameter<bool>("doubleMatching")),
    doubleMatching_nRequiredCommonTracks_(iConfig.getParameter<int>("doubleMatching_nRequiredCommonTracks")),
    doubleMatching_maxSignificance_(iConfig.getParameter<double>("doubleMatching_maxSignificance")),
    doubleMatching_dR_max_(iConfig.getParameter<double>("doubleMatching_dR_max")),
    doubleMatching_relPt_max_(iConfig.getParameter<double>("doubleMatching_relPt_max")),
    trkMaxDeltaR_(iConfig.getParameter<double>("trkMaxDeltaR")),
    trkMaxDPtRel_(iConfig.getParameter<double>("trkMaxDPtRel")),
    trkCheckCharge_(iConfig.getParameter<bool>("trkCheckCharge")),
    trkResolveAmbiguities_(iConfig.getParameter<bool>("trkResolveAmbiguities"))
{
    produces<nanoaod::FlatTable>("GVTable");
    produces<nanoaod::FlatTable>("rejectedGVTable");
    produces<nanoaod::FlatTable>("GVDaughtersTable");
    produces<nanoaod::FlatTable>("GVDirectDaughters");
    produces<nanoaod::FlatTable>("SVGVMatchTable");
    produces<nanoaod::FlatTable>("SVGVtrkMatchTable"); 

}


void GenVertexProducer::produce(edm::Event& iEvent,
             const edm::EventSetup&) 
    {
	edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosHandle;
	iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfosHandle);
        edm::Handle<edm::View<reco::Candidate>> genHandle;
        iEvent.getByToken(genToken_, genHandle);
	edm::Handle<reco::GenParticleCollection> prunedGenHandle;
	iEvent.getByToken(prunedGenToken_, prunedGenHandle);
        edm::Handle<std::vector<reco::Vertex>> svHandle;
        iEvent.getByToken(svToken_, svHandle);
        auto pvsIn = iEvent.getHandle(pvs_);
        edm::Handle<std::vector<reco::Track>> tracksHandle;   // NEW
        iEvent.getByToken(tracksToken_, tracksHandle);

        if (!genHandle.isValid() ||!prunedGenHandle.isValid()|| !svHandle.isValid() || !pvsIn.isValid() ||
            pvsIn->empty() || !tracksHandle.isValid()) {
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "GV", false), "GVTable");
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "RejectedGV", false), "rejectedGVTable");
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "GVDaughters", false), "GVDaughtersTable");
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "GVDirectDaughters", false), "GVDirectDaughters");
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "mySV", false, true), "SVGVMatchTable");
            iEvent.put(std::make_unique<nanoaod::FlatTable>(0, "mySVtrks", false, true), "SVGVtrkMatchTable");
            return;
        }

        const auto& genParticles = genHandle;
	const auto& prunedGenParticles = *prunedGenHandle;
        const auto& secondaryVertices = svHandle;
        const auto& tracks = *tracksHandle;                   // NEW
	const bool hasJetFlavourInfos = jetFlavourInfosHandle.isValid();

	struct GenJetMatchInfo {
	    int jetIdx = -1;
	
	    float pt = -1.f;
	    float eta = 0.f;
	    float phi = 0.f;
	    float mass = -1.f;
	
	    int hadronFlavour = 0;
	    int partonFlavour = 0;
	    int nCHadrons = 0;
	    int nBHadrons = 0;
	};

        // Output vectors
        std::vector<float> Hadron_pt, Hadron_eta, Hadron_phi;
	std::vector<int> Hadron_mergedGenPartIdx;
	std::vector<int> Hadron_prunedGenPartIdx;
	std::vector<int> Hadron_genJetIdx;
	std::vector<int> Hadron_hasGenJet;
	
	std::vector<float> Hadron_genJetPt;
	std::vector<float> Hadron_genJetEta;
	std::vector<float> Hadron_genJetPhi;
	std::vector<float> Hadron_genJetMass;
	std::vector<float> Hadron_genJetDeltaR;
	
	std::vector<int> Hadron_genJetHadronFlavour;
	std::vector<int> Hadron_genJetPartonFlavour;
	std::vector<int> Hadron_genJetNCHadrons;
	std::vector<int> Hadron_genJetNBHadrons;
        std::vector<float> SV_x, SV_y, SV_z, SV_eta, SV_phi;
        std::vector<CovMatrix> SV_cov;
        std::vector<float> Hadron_GVx, Hadron_GVy, Hadron_GVz;
        std::vector<float>  Hadron_GVx_i, Hadron_GVy_i, Hadron_GVz_i;
        std::vector<int> Hadron_pdgId;
        std::vector<int> Hadron_pdgClass, Hadron_isB, Hadron_isD, Hadron_isBtoD;
        std::vector<int> Hadron_fromHF, Hadron_toHF;
        std::vector<int> Hadron_nDauNoRecognizedSecondaryAncestor;
        std::vector<int> Hadron_nDauFromB, Hadron_nDauFromBC, Hadron_nDauFromC;
        std::vector<int> Hadron_nDauOtherSecondary, Hadron_nDauOriginUnknown;
        std::vector<int> GV_nDaughters; // NEW: total number of GVDaughters per GV (denominator for nDaughtersMatchedToTracks)
        std::vector<float> GV_maxDaughterPairDeltaR; // NEW: max pairwise deltaR among a GV's own GVDaughters (no tracks involved)
        std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
        std::vector<float> Daughters_vx, Daughters_vy, Daughters_vz;
        std::vector<int> Daughters_charge, Daughters_GVidx, Daughters_pdgId;
        std::vector<int> Daughters_originLabel;
        
        std::vector<float> allHadron_GVx, allHadron_GVy, allHadron_GVz;
        std::vector<float>  allHadron_GVx_i, allHadron_GVy_i, allHadron_GVz_i;
        std::vector<int> allHadron_pdgId;
	std::vector<int> allHadron_isB, allHadron_isD, allHadron_isBtoD;
        std::vector<float> allHadron_pt, allHadron_eta, allHadron_phi;
	std::vector<int> allHadron_mergedGenPartIdx;
	std::vector<int> allHadron_prunedGenPartIdx;
	std::vector<int> allHadron_genJetIdx;
	std::vector<int> allHadron_hasGenJet;
	
	std::vector<float> allHadron_genJetPt;
	std::vector<float> allHadron_genJetEta;
	std::vector<float> allHadron_genJetPhi;
	std::vector<float> allHadron_genJetDeltaR;

	std::vector<int> Hadron_hasPrunedGenMatch;
	std::vector<int> allHadron_hasPrunedGenMatch;

        std::vector<float> directDaughters_pt, directDaughters_eta, directDaughters_phi;
        std::vector<int> directDaughters_charge, directDaughters_GVidx, directDaughters_pdgId;
        VertexDistance3D vdist;
        const auto& PV0 = pvsIn->front();

	std::unordered_map<unsigned int, GenJetMatchInfo> genParticleToJet;

	if (hasJetFlavourInfos) {
	    for (size_t iJet = 0; iJet < jetFlavourInfosHandle->size(); ++iJet) {
	        const auto jetAndFlavour = (*jetFlavourInfosHandle)[iJet];
	
	        const auto& jetRef = jetAndFlavour.first;
	        const auto& flavourInfo = jetAndFlavour.second;
	
	        if (jetRef.isNull()) {
	            continue;
	        }
	
	        GenJetMatchInfo info;
		info.jetIdx = static_cast<int>(jetRef.key());
	
	        info.pt = jetRef->pt();
	        info.eta = jetRef->eta();
	        info.phi = jetRef->phi();
	        info.mass = jetRef->mass();
	
	        info.hadronFlavour = flavourInfo.getHadronFlavour();
	        info.partonFlavour = flavourInfo.getPartonFlavour();
	        info.nCHadrons =
	            static_cast<int>(flavourInfo.getcHadrons().size());
	        info.nBHadrons =
	            static_cast<int>(flavourInfo.getbHadrons().size());
	
	        for (const auto& cHadronRef : flavourInfo.getcHadrons()) {
	            if (cHadronRef.isNull()) {
	                continue;
	            }
	
	            const unsigned int genPartIdx = cHadronRef.key();
	
	            const auto existing = genParticleToJet.find(genPartIdx);
	
	            if (existing == genParticleToJet.end()) {
	                genParticleToJet.emplace(genPartIdx, info);
	            } else {
	                // This should normally not happen. Retain the harder jet
	                // deterministically if the same hadron appears more than once.
	                if (info.pt > existing->second.pt) {
	                    existing->second = info;
	                }
	            }
	        }
	
	        // Also map B hadrons, since your GV collection contains B GVs.
	        for (const auto& bHadronRef : flavourInfo.getbHadrons()) {
	            if (bHadronRef.isNull()) {
	                continue;
	            }
	
	            const unsigned int genPartIdx = bHadronRef.key();
	
	            const auto existing = genParticleToJet.find(genPartIdx);
	
	            if (existing == genParticleToJet.end()) {
	                genParticleToJet.emplace(genPartIdx, info);
	            } else if (info.pt > existing->second.pt) {
	                existing->second = info;
	            }
	        }
	    }
	}

        // save coordinates of SV (will be used for matching with GV)
        for (auto const& sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                SV_x.push_back(sv.x());
                SV_y.push_back(sv.y());
                SV_z.push_back(sv.z());

                // Get Eta and Phi from tracks
                TLorentzVector p4s_SV = TLorentzVector(0,0,0,0);
                for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                    const edm::RefToBase<reco::Track>& trkRef = *it;
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                    p4s_SV += p4;
                }
                SV_eta.push_back(p4s_SV.Eta());
                SV_phi.push_back(p4s_SV.Phi());
                SV_cov.push_back(sv.covariance());
            }
        }

        // Filling Hadrons and Daughters vectors
        int ngv=0;
        int nRejectedGV=0;
        int ngv_b=0, ngv_d=0, ngv_s=0, ngv_tau=0;
        for(size_t i=0; i<genParticles->size(); ++i){
            const reco::Candidate* hadron = &(*genParticles)[i];
            //std::cout<<"Hadron "<<i<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
            if(!(hadron->pt()>hadPt_min_ && std::abs(hadron->eta())<2.5)) continue;

            int hadPDG = checkPDG(std::abs(hadron->pdgId())); // 1: Beauty, 2: Charmed, 3: Strange,  4: Tau,  0: Else
            if(hadPDG==0) continue;
            //     code here
            //    
            //    
            //    
            //    
            //  

                


            //  Collect stable charged daughters
            std::vector<float> temp_pt, temp_eta, temp_phi, temp_vx, temp_vy, temp_vz; // kinematics of gen daughters of the hadron in the loop
            std::vector<int> temp_charge, temp_GVidx, temp_flav, temp_pdgId, temp_originLabel;
            int nPack=0;
            float vx=std::numeric_limits<float>::quiet_NaN();
            float vy=std::numeric_limits<float>::quiet_NaN();
            float vz=std::numeric_limits<float>::quiet_NaN();

            for(size_t j=0; j<genParticles->size(); ++j){
                const reco::Candidate* dau = &(*genParticles)[j];
                if(dau==hadron) continue;
                if(!(dau->status()==1 && dau->charge()!=0 && dau->pt()>0.4 && std::abs(dau->eta())<2.5)) continue;

                auto GV = isAncestor(hadron,dau); //takes the x,y,z of the daughter (decay point of the hadron) if daughters otherwise return nan
                if(GV.has_value()){
                    std::tie(vx,vy,vz) = *GV;
                    if(!std::isnan(vx)){
                        nPack++;
                        temp_pt.push_back(dau->pt());
                        temp_eta.push_back(dau->eta());
                        temp_phi.push_back(dau->phi());
                        temp_vx.push_back(dau->vx());
                        temp_vy.push_back(dau->vy());
                        temp_vz.push_back(dau->vz());
                        temp_charge.push_back(dau->charge());
                        temp_pdgId.push_back(dau->pdgId());
                        temp_GVidx.push_back(ngv); // hadron index
                        temp_originLabel.push_back(getDaughterOriginLabelNoPU(dau));
                    }
                }
            }
            // If has more than 2 good daughters, the Hadron is Good, we found a GV:
            if(nPack>=1){
                // Save hadron info
                //std::cout<<"Found hadron "<<ngv<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
                Hadron_pt.push_back(hadron->pt());
                Hadron_eta.push_back(hadron->eta());
                Hadron_phi.push_back(hadron->phi());
                Hadron_pdgId.push_back(hadron->pdgId());
                Hadron_pdgClass.push_back(hadPDG);

		const int mergedGenPartIdx = static_cast<int>(i);

		const int prunedGenPartIdx =
		    findMatchingPrunedHadron(
		        hadron,
		        prunedGenParticles
		    );
		
		Hadron_mergedGenPartIdx.push_back(mergedGenPartIdx);
		Hadron_prunedGenPartIdx.push_back(prunedGenPartIdx);

		Hadron_hasPrunedGenMatch.push_back(prunedGenPartIdx >= 0 ? 1 : 0);
		
		auto jetMatch = genParticleToJet.end();
		
		if (prunedGenPartIdx >= 0) {
		    jetMatch = genParticleToJet.find(
		        static_cast<unsigned int>(prunedGenPartIdx)
		    );
		}
		
		if (jetMatch != genParticleToJet.end()) {
		    const auto& info = jetMatch->second;
		
		    Hadron_genJetIdx.push_back(info.jetIdx);
		    Hadron_hasGenJet.push_back(1);
		
		    Hadron_genJetPt.push_back(info.pt);
		    Hadron_genJetEta.push_back(info.eta);
		    Hadron_genJetPhi.push_back(info.phi);
		    Hadron_genJetMass.push_back(info.mass);
		
		    Hadron_genJetDeltaR.push_back(
		        deltaR(
		            static_cast<float>(hadron->eta()),
		            static_cast<float>(hadron->phi()),
		            info.eta,
		            info.phi
		        )
		    );
		
		    Hadron_genJetHadronFlavour.push_back(info.hadronFlavour);
		    Hadron_genJetPartonFlavour.push_back(info.partonFlavour);
		    Hadron_genJetNCHadrons.push_back(info.nCHadrons);
		    Hadron_genJetNBHadrons.push_back(info.nBHadrons);

		} else {
		    Hadron_genJetIdx.push_back(-1);
		    Hadron_hasGenJet.push_back(0);
		
		    Hadron_genJetPt.push_back(-1.f);
		    Hadron_genJetEta.push_back(0.f);
		    Hadron_genJetPhi.push_back(0.f);
		    Hadron_genJetMass.push_back(-1.f);
		    Hadron_genJetDeltaR.push_back(-1.f);
		
		    Hadron_genJetHadronFlavour.push_back(0);
		    Hadron_genJetPartonFlavour.push_back(0);
		    Hadron_genJetNCHadrons.push_back(0);
		    Hadron_genJetNBHadrons.push_back(0);
		}
                
                

                // Save GenVertex
                ngv++;
                if(hadPDG==1) {
                    ngv_b++;
                    Hadron_isB.push_back(1);
                    }
                else{
                    Hadron_isB.push_back(0);
                }
                if(hadPDG==2) {
                    ngv_d++;
                    Hadron_isD.push_back(1);
                    }
                else{
                    Hadron_isD.push_back(0);
                }
                if(hadPDG==3) ngv_s++;
                if(hadPDG==4) ngv_tau++;
                Hadron_fromHF.push_back(hasHFAncestor(hadron) ? 1 : 0);
                Hadron_toHF.push_back(hasHFDescendant(hadron) ? 1 : 0);
                Hadron_isBtoD.push_back(hadPDG == 2 && hasBHadronAncestor(hadron) ? 1 : 0);

                int nNoRecognizedSecondaryAncestor = 0;
                int nFromB = 0;
                int nFromBC = 0;
                int nFromC = 0;
                int nOtherSecondary = 0;
                int nOriginUnknown = 0;
                for (const int label : temp_originLabel) {
                    if (label == 0) ++nNoRecognizedSecondaryAncestor;
                    else if (label == 2) ++nFromB;
                    else if (label == 3) ++nFromBC;
                    else if (label == 4) ++nFromC;
                    else if (label == 5) ++nOtherSecondary;
                    else ++nOriginUnknown;
                }
                Hadron_nDauNoRecognizedSecondaryAncestor.push_back(nNoRecognizedSecondaryAncestor);
                Hadron_nDauFromB.push_back(nFromB);
                Hadron_nDauFromBC.push_back(nFromBC);
                Hadron_nDauFromC.push_back(nFromC);
                Hadron_nDauOtherSecondary.push_back(nOtherSecondary);
                Hadron_nDauOriginUnknown.push_back(nOriginUnknown);
                GV_nDaughters.push_back(nPack);

                // NEW: max pairwise deltaR among this GV's own daughters (gen-level, no tracks involved)
                float maxPairDR = 0.f;
                for (size_t a = 0; a < temp_eta.size(); ++a) {
                    for (size_t b = a + 1; b < temp_eta.size(); ++b) {
                        float dR = deltaR(temp_eta[a], temp_phi[a], temp_eta[b], temp_phi[b]);
                        if (dR > maxPairDR) maxPairDR = dR;
                    }
                }
                GV_maxDaughterPairDeltaR.push_back(maxPairDR);
                
                Hadron_GVx.push_back(vx);               // point of decay of the hadron
                Hadron_GVy.push_back(vy);               // point of decay of the hadron
                Hadron_GVz.push_back(vz);               // point of decay of the hadron
                Hadron_GVx_i.push_back(hadron->vx());   // point of origin of the hadron
                Hadron_GVy_i.push_back(hadron->vy());   // point of origin of the hadron
                Hadron_GVz_i.push_back(hadron->vz());   // point of origin of the hadron

                // Save daughters
                Daughters_pt.insert(Daughters_pt.end(), temp_pt.begin(), temp_pt.end());
                Daughters_vx.insert(Daughters_vx.end(), temp_vx.begin(), temp_vx.end());
                Daughters_vy.insert(Daughters_vy.end(), temp_vy.begin(), temp_vy.end());
                Daughters_vz.insert(Daughters_vz.end(), temp_vz.begin(), temp_vz.end());
                Daughters_eta.insert(Daughters_eta.end(), temp_eta.begin(), temp_eta.end());
                Daughters_phi.insert(Daughters_phi.end(), temp_phi.begin(), temp_phi.end());
                Daughters_charge.insert(Daughters_charge.end(), temp_charge.begin(), temp_charge.end());
                Daughters_pdgId.insert(Daughters_pdgId.end(), temp_pdgId.begin(), temp_pdgId.end());
                Daughters_GVidx.insert(Daughters_GVidx.end(), temp_GVidx.begin(), temp_GVidx.end());
                Daughters_originLabel.insert(Daughters_originLabel.end(), temp_originLabel.begin(), temp_originLabel.end());
                for(size_t j=0; j<genParticles->size(); ++j){   
                    const reco::Candidate* dau = &(*genParticles)[j];
                    if(dau==hadron) continue;
                    if (dau->numberOfMothers() > 0){
                        const reco::Candidate* mother = dau->mother(0);
                        if (mother == hadron){
                            directDaughters_pt.push_back(dau->pt());
                            directDaughters_eta.push_back(dau->eta());
                            directDaughters_phi.push_back(dau->phi());
                            directDaughters_charge.push_back(dau->charge());
                            directDaughters_pdgId.push_back(dau->pdgId());
                            directDaughters_GVidx.push_back(ngv-1);
                        }
                    }
                }
            }
            else{
                //fallback to save the decay point of the hadron if it has no daughters, but is still a rejected GV
                if (std::isnan(vx) && hadron->numberOfDaughters() > 0) {
                    const reco::Candidate* directDau = hadron->daughter(0);
                    vx = directDau->vx();
                    vy = directDau->vy();
                    vz = directDau->vz();
                }

		if(hadPDG==1) {
                    allHadron_isB.push_back(1);
                    }
                else{
                    allHadron_isB.push_back(0);
                }
                if(hadPDG==2) {
                    allHadron_isD.push_back(1);
                }
                else{
                    allHadron_isD.push_back(0);
                }
                allHadron_isBtoD.push_back(hadPDG == 2 && hasBHadronAncestor(hadron) ? 1 : 0);
                nRejectedGV++;
                allHadron_pt.push_back(hadron->pt());
                allHadron_eta.push_back(hadron->eta());
                allHadron_phi.push_back(hadron->phi());
                allHadron_pdgId.push_back(hadron->pdgId());
                allHadron_GVx.push_back(vx);               // point of decay of the hadron
                allHadron_GVy.push_back(vy);               // point of decay of the hadron
                allHadron_GVz.push_back(vz);               // point of decay of the hadron
                allHadron_GVx_i.push_back(hadron->vx());   // point of origin of the hadron
                allHadron_GVy_i.push_back(hadron->vy());   // point of origin of the hadron
                allHadron_GVz_i.push_back(hadron->vz());   // point of origin of the hadron
	
		const int mergedGenPartIdx = static_cast<int>(i);
		
		const int prunedGenPartIdx =
		    findMatchingPrunedHadron(
		        hadron,
		        prunedGenParticles
		    );
		
		allHadron_mergedGenPartIdx.push_back(
		    mergedGenPartIdx
		);
		
		allHadron_prunedGenPartIdx.push_back(prunedGenPartIdx);

		allHadron_hasPrunedGenMatch.push_back(prunedGenPartIdx >= 0 ? 1 : 0);
		
		auto jetMatch = genParticleToJet.end();
		
		if (prunedGenPartIdx >= 0) {
		    jetMatch = genParticleToJet.find(
		        static_cast<unsigned int>(prunedGenPartIdx)
		    );
		}
		
                if (jetMatch != genParticleToJet.end()) {
                    const auto& info = jetMatch->second;

                    allHadron_genJetIdx.push_back(info.jetIdx);
                    allHadron_hasGenJet.push_back(1);

                    allHadron_genJetPt.push_back(info.pt);
                    allHadron_genJetEta.push_back(info.eta);
                    allHadron_genJetPhi.push_back(info.phi);

                    allHadron_genJetDeltaR.push_back(
                        deltaR(
                            static_cast<float>(hadron->eta()),
                            static_cast<float>(hadron->phi()),
                            info.eta,
                            info.phi
                        )
                    );

                } else {
                    allHadron_genJetIdx.push_back(-1);
                    allHadron_hasGenJet.push_back(0);

                    allHadron_genJetPt.push_back(-1.f);
                    allHadron_genJetEta.push_back(0.f);
                    allHadron_genJetPhi.push_back(0.f);
                    allHadron_genJetDeltaR.push_back(-1.f);
                }
            }
        }


        // Filling tracks from reco SV
        std::vector<float> SVtrk_pt, SVtrk_eta, SVtrk_phi;
        std::vector<int> SVtrk_SVidx;
        int SV_index=0;
        for (const auto &sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                if (trkRef.isNull()) continue;
                TLorentzVector p4;
                p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                SVtrk_pt.push_back(trkRef->pt());
                SVtrk_eta.push_back(trkRef->eta());
                SVtrk_phi.push_back(trkRef->phi());
                SVtrk_SVidx.push_back(SV_index); 
            }
            SV_index++;
            }
        }

        
        // Compute matrix of distances between SV and GV
        auto distances = computeDistanceMatrix(SV_x, SV_y, SV_z, SV_cov,Hadron_GVx, Hadron_GVy, Hadron_GVz);
        //printDistanceMatrix(distances);
        
        std::vector<int> Hadron_SVIdx(ngv, -1); // 
        std::vector<float> Hadron_SVDistance(ngv, -1); // 
        std::vector<float> Hadron_minDistNotMatched(ngv, 999.f);

        // perform matching based on distance matrix and track-to-daughter matching
        auto result = matchHadronsToSV(distances,SVtrk_pt, SVtrk_eta, SVtrk_phi, SVtrk_SVidx,Daughters_pt, Daughters_eta, Daughters_phi, Daughters_GVidx,
                                        SV_eta,SV_phi,Hadron_eta,Hadron_phi,ngv,nRequiredCommonTracks_,dR_max_,relPt_max_,doubleMatching_ , doubleMatching_nRequiredCommonTracks_, doubleMatching_maxSignificance_ , doubleMatching_dR_max_, doubleMatching_relPt_max_  );
        
        Hadron_SVIdx             = std::get<0>(result);
        Hadron_SVDistance        = std::get<1>(result);
        Hadron_minDistNotMatched = std::get<2>(result);
        std::vector<int> SVtrk_isMatched = std::get<3>(result);
        std::vector<int> SVtrk_GVIdx     = std::get<4>(result);
        std::vector<int> SVtrk_daughterIdx = std::get<5>(result);  // NEW

        // ---------------------------------------------------------------------
        // NEW DIAGNOSTIC: characterize the spatially closest SV candidate for
        // every accepted GV BEFORE the greedy one-to-one assignment is applied.
        //
        // Important:
        //   * `distances` in this scope is still the original unmodified matrix,
        //     because matchHadronsToSV() receives it by value.
        //   * These quantities are diagnostic only; they do not alter matching.
        //   * Daughter matching here is one-to-one in the daughter index so that
        //     one GV daughter cannot be counted multiple times for the same SV.
        // ---------------------------------------------------------------------
        std::vector<int> GV_bestCandidateSVIdx(ngv, -1);
        std::vector<float> GV_bestCandidateDistanceSig(ngv, -1.f);
        std::vector<int> GV_bestCandidateNCommonTracks(ngv, 0);
        std::vector<int> GV_bestCandidateNMatchedDaughters(ngv, 0);
        std::vector<float> GV_bestCandidateMaxTrackDeltaR(ngv, -1.f);
        std::vector<float> GV_bestCandidateMaxRelPtDiff(ngv, -1.f);
        std::vector<float> GV_bestCandidateDeltaR(ngv, -1.f);
        std::vector<int> GV_bestCandidateWasFinalMatch(ngv, 0);

        for (int had = 0; had < ngv; ++had) {
            float bestDist = 999.f;
            int bestSV = -1;

            // Find the spatially closest reconstructed SV for this GV using
            // the original covariance-weighted distance-significance matrix.
            for (size_t sv = 0; sv < distances.size(); ++sv) {
                if (had >= static_cast<int>(distances[sv].size())) continue;
                const float dist = distances[sv][had];
                if (std::isfinite(dist) && dist < bestDist && dist < 997.f) {
                    bestDist = dist;
                    bestSV = static_cast<int>(sv);
                }
            }

            if (bestSV < 0) continue;

            GV_bestCandidateSVIdx[had] = bestSV;
            GV_bestCandidateDistanceSig[had] = bestDist;

            if (bestSV < static_cast<int>(SV_eta.size()) &&
                had < static_cast<int>(Hadron_eta.size())) {
                GV_bestCandidateDeltaR[had] =
                    deltaR(SV_eta[bestSV], SV_phi[bestSV],
                           Hadron_eta[had], Hadron_phi[had]);
            }

            // Collect the selected reco tracks belonging to this candidate SV.
            std::vector<size_t> candidateSVTracks;
            for (size_t itrk = 0; itrk < SVtrk_SVidx.size(); ++itrk) {
                if (SVtrk_SVidx[itrk] == bestSV &&
                    SVtrk_pt[itrk] > 0.4 &&
                    std::fabs(SVtrk_eta[itrk]) < 2.5) {
                    candidateSVTracks.push_back(itrk);
                }
            }

            // Collect this GV's selected gen daughters.
            std::vector<size_t> candidateGVDaughters;
            for (size_t idau = 0; idau < Daughters_GVidx.size(); ++idau) {
                if (Daughters_GVidx[idau] == had) {
                    candidateGVDaughters.push_back(idau);
                }
            }

            // One-to-one daughter bookkeeping for the diagnostic. Each SV
            // track can match at most one daughter, and each daughter can be
            // used at most once.
            std::vector<bool> daughterUsed(candidateGVDaughters.size(), false);

            int nCommon = 0;
            float maxMatchedDR = -1.f;
            float maxMatchedRelPt = -1.f;

            for (size_t itrk : candidateSVTracks) {
                int bestLocalDau = -1;
                float bestLocalDR = std::numeric_limits<float>::max();
                float bestLocalRelPt = -1.f;

                for (size_t ilocal = 0; ilocal < candidateGVDaughters.size(); ++ilocal) {
                    if (daughterUsed[ilocal]) continue;

                    const size_t idau = candidateGVDaughters[ilocal];
                    const float dR = deltaR(
                        SVtrk_eta[itrk], SVtrk_phi[itrk],
                        Daughters_eta[idau], Daughters_phi[idau]);
                    const float relPt =
                        std::fabs(SVtrk_pt[itrk] - Daughters_pt[idau]) /
                        std::max(Daughters_pt[idau], 1e-6f);

                    if (dR < dR_max_ && relPt < relPt_max_ && dR < bestLocalDR) {
                        bestLocalDau = static_cast<int>(ilocal);
                        bestLocalDR = dR;
                        bestLocalRelPt = relPt;
                    }
                }

                if (bestLocalDau >= 0) {
                    daughterUsed[bestLocalDau] = true;
                    ++nCommon;
                    maxMatchedDR = std::max(maxMatchedDR, bestLocalDR);
                    maxMatchedRelPt = std::max(maxMatchedRelPt, bestLocalRelPt);
                }
            }

            GV_bestCandidateNCommonTracks[had] = nCommon;
            GV_bestCandidateNMatchedDaughters[had] = nCommon;
            GV_bestCandidateMaxTrackDeltaR[had] = maxMatchedDR;
            GV_bestCandidateMaxRelPtDiff[had] = maxMatchedRelPt;

            if (Hadron_SVIdx[had] >= 0 && Hadron_SVIdx[had] == bestSV) {
                GV_bestCandidateWasFinalMatch[had] = 1;
            }
        }
        auto svTrkGVTable = std::make_unique<nanoaod::FlatTable>(SVtrk_pt.size(), "mySVtrks", false, true);
        svTrkGVTable->addColumn<int>("isMatched", SVtrk_isMatched, "1 if track is matched to a genParticle daughter of the matched GV");
        svTrkGVTable->addColumn<int>("GVIdx", SVtrk_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
        svTrkGVTable->addColumn<int>("daughterIdx", SVtrk_daughterIdx, "Index into GVDaughters table of the matched gen daughter, -1 if unmatched"); // NEW

        iEvent.put(std::move(svTrkGVTable), "SVGVtrkMatchTable");

        // Build the complete per-SV truth record from the matching performed above.
        const size_t nSV = SV_x.size();
        std::vector<int> SV_GVIdx(nSV, -1);
        std::vector<int> SV_isMatched(nSV, 0);
        std::vector<int> SV_nMatchedGV(nSV, 0);
        std::vector<int> SV_nCommonTracks(nSV, 0);
        std::vector<int> SV_nMatchedDaughters(nSV, 0);
        std::vector<float> SV_bestMatchScore(nSV, -1.f);
        std::vector<float> SV_bestMatchDistanceSig(nSV, -1.f);

        std::vector<int> SV_truthPdgClass(nSV, 0), SV_truthPdgId(nSV, 0);
        std::vector<int> SV_truthIsB(nSV, 0), SV_truthIsD(nSV, 0), SV_truthIsBtoD(nSV, 0);
        std::vector<int> SV_nOriginNoRecognizedSecondaryAncestor(nSV, 0);
        std::vector<int> SV_nOriginFromB(nSV, 0), SV_nOriginFromBC(nSV, 0), SV_nOriginFromC(nSV, 0);
        std::vector<int> SV_nOriginOtherSecondary(nSV, 0), SV_nOriginUnknown(nSV, 0);
        std::vector<float> SV_fracOriginNoRecognizedSecondaryAncestor(nSV, -1.f);
        std::vector<float> SV_fracOriginHF(nSV, -1.f);
        std::vector<int> SV_dominantOriginLabel(nSV, 9);

        std::vector<int> SV_matchedGV_nDauNoRecognizedSecondaryAncestor(nSV, -1);
        std::vector<int> SV_matchedGV_nDauFromB(nSV, -1), SV_matchedGV_nDauFromBC(nSV, -1);
        std::vector<int> SV_matchedGV_nDauFromC(nSV, -1), SV_matchedGV_nDauOtherSecondary(nSV, -1);
        std::vector<int> SV_matchedGV_nDauOriginUnknown(nSV, -1);

        for (int had = 0; had < ngv; ++had) {
            const int sv = Hadron_SVIdx[had];
            if (sv < 0 || sv >= static_cast<int>(nSV)) continue;
            SV_GVIdx[sv] = had;
            SV_isMatched[sv] = 1;
            SV_nMatchedGV[sv] = 1;
            SV_bestMatchDistanceSig[sv] = Hadron_SVDistance[had];
            SV_truthPdgClass[sv] = Hadron_pdgClass[had];
            SV_truthPdgId[sv] = Hadron_pdgId[had];
            SV_truthIsB[sv] = Hadron_isB[had];
            SV_truthIsD[sv] = Hadron_isD[had];
            SV_truthIsBtoD[sv] = Hadron_isBtoD[had];
            SV_matchedGV_nDauNoRecognizedSecondaryAncestor[sv] = Hadron_nDauNoRecognizedSecondaryAncestor[had];
            SV_matchedGV_nDauFromB[sv] = Hadron_nDauFromB[had];
            SV_matchedGV_nDauFromBC[sv] = Hadron_nDauFromBC[had];
            SV_matchedGV_nDauFromC[sv] = Hadron_nDauFromC[had];
            SV_matchedGV_nDauOtherSecondary[sv] = Hadron_nDauOtherSecondary[had];
            SV_matchedGV_nDauOriginUnknown[sv] = Hadron_nDauOriginUnknown[had];
        }

        // Derive matched-track and matched-daughter summaries from the exact associations
        // produced by matchHadronsToSV. This guarantees a single matching source of truth.
        for (size_t itrk = 0; itrk < SVtrk_isMatched.size(); ++itrk) {
            if (!SVtrk_isMatched[itrk]) continue;
            const int sv = SVtrk_SVidx[itrk];
            const int gv = SVtrk_GVIdx[itrk];
            const int dau = SVtrk_daughterIdx[itrk];
            if (sv < 0 || sv >= static_cast<int>(nSV) ||
                gv < 0 || gv >= ngv ||
                dau < 0 || dau >= static_cast<int>(Daughters_pt.size())) continue;

            ++SV_nCommonTracks[sv];
            ++SV_nMatchedDaughters[sv];
            const float dR = deltaR(SVtrk_eta[itrk], SVtrk_phi[itrk], Daughters_eta[dau], Daughters_phi[dau]);
            const float relPt = std::fabs(SVtrk_pt[itrk] - Daughters_pt[dau]) /
                                std::max(Daughters_pt[dau], 1e-6f);
            if (SV_bestMatchScore[sv] < 0.f) SV_bestMatchScore[sv] = 0.f;
            SV_bestMatchScore[sv] += std::max(0.f, 1.f - dR / static_cast<float>(dR_max_)) +
                                     std::max(0.f, 1.f - relPt / static_cast<float>(relPt_max_));

            const int label = Daughters_originLabel[dau];
            if (label == 0) ++SV_nOriginNoRecognizedSecondaryAncestor[sv];
            else if (label == 2) ++SV_nOriginFromB[sv];
            else if (label == 3) ++SV_nOriginFromBC[sv];
            else if (label == 4) ++SV_nOriginFromC[sv];
            else if (label == 5) ++SV_nOriginOtherSecondary[sv];
            else ++SV_nOriginUnknown[sv];
        }

        for (size_t sv = 0; sv < nSV; ++sv) {
            const int nLab = SV_nMatchedDaughters[sv];
            if (nLab <= 0) continue;
            SV_fracOriginNoRecognizedSecondaryAncestor[sv] =
                static_cast<float>(SV_nOriginNoRecognizedSecondaryAncestor[sv]) / nLab;
            const int nHF = SV_nOriginFromB[sv] + SV_nOriginFromBC[sv] + SV_nOriginFromC[sv];
            SV_fracOriginHF[sv] = static_cast<float>(nHF) / nLab;
            const std::vector<std::pair<int,int>> counts = {
                {0, SV_nOriginNoRecognizedSecondaryAncestor[sv]},
                {2, SV_nOriginFromB[sv]}, {3, SV_nOriginFromBC[sv]},
                {4, SV_nOriginFromC[sv]}, {5, SV_nOriginOtherSecondary[sv]},
                {9, SV_nOriginUnknown[sv]}};
            SV_dominantOriginLabel[sv] = std::max_element(
                counts.begin(), counts.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; })->first;
        }

        auto svGVTable = std::make_unique<nanoaod::FlatTable>(nSV, "mySV", false, true);
        svGVTable->addColumn<int>("GVIdx", SV_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
        svGVTable->addColumn<int>("isMatched", SV_isMatched, "1 if SV matched to a GV");
        svGVTable->addColumn<int>("nMatchedGV", SV_nMatchedGV, "One if this SV was matched to a GV");
        svGVTable->addColumn<int>("nCommonTracks", SV_nCommonTracks, "Number of matched SV tracks");
        svGVTable->addColumn<int>("nMatchedDaughters", SV_nMatchedDaughters, "Number of matched GV daughters");
        svGVTable->addColumn<float>("bestMatchScore", SV_bestMatchScore, "Track-sharing score from GenVertexProducer matching");
        svGVTable->addColumn<float>("bestMatchDistanceSig", SV_bestMatchDistanceSig, "Matched SV-GV distance significance");
        svGVTable->addColumn<int>("truth_pdgClass", SV_truthPdgClass, "Truth class of matched GV");
        svGVTable->addColumn<int>("truth_pdgId", SV_truthPdgId, "PDG id of matched GV hadron");
        svGVTable->addColumn<int>("truth_isB", SV_truthIsB, "Matched GV is a B hadron");
        svGVTable->addColumn<int>("truth_isD", SV_truthIsD, "Matched GV is a D hadron");
        svGVTable->addColumn<int>("truth_isBtoD", SV_truthIsBtoD, "Matched GV is D from B");
       svGVTable->addColumn<int>("nOriginNoRecognizedSecondaryAncestor", SV_nOriginNoRecognizedSecondaryAncestor, "Matched daughters with originLabel 0");
       svGVTable->addColumn<int>("nOriginFromB", SV_nOriginFromB, "Matched daughters with B ancestry");
       svGVTable->addColumn<int>("nOriginFromBC", SV_nOriginFromBC, "Matched daughters with B and C ancestry");
       svGVTable->addColumn<int>("nOriginFromC", SV_nOriginFromC, "Matched daughters with C ancestry");
       svGVTable->addColumn<int>("nOriginOtherSecondary", SV_nOriginOtherSecondary, "Matched daughters from other secondary sources");
       svGVTable->addColumn<int>("nOriginUnknown", SV_nOriginUnknown, "Matched daughters with unknown origin");
       svGVTable->addColumn<float>("fracOriginNoRecognizedSecondaryAncestor", SV_fracOriginNoRecognizedSecondaryAncestor, "Fraction of matched daughters with originLabel 0");
       svGVTable->addColumn<float>("fracOriginHF", SV_fracOriginHF, "Fraction of matched daughters with B/BC/C labels");
       svGVTable->addColumn<int>("dominantOriginLabel", SV_dominantOriginLabel, "Dominant matched-daughter origin label");
       svGVTable->addColumn<int>("matchedGV_nDauNoRecognizedSecondaryAncestor", SV_matchedGV_nDauNoRecognizedSecondaryAncestor, "Matched GV count for originLabel 0");
       svGVTable->addColumn<int>("matchedGV_nDauFromB", SV_matchedGV_nDauFromB, "Matched GV daughter count from B");
       svGVTable->addColumn<int>("matchedGV_nDauFromBC", SV_matchedGV_nDauFromBC, "Matched GV daughter count from BC");
       svGVTable->addColumn<int>("matchedGV_nDauFromC", SV_matchedGV_nDauFromC, "Matched GV daughter count from C");
       svGVTable->addColumn<int>("matchedGV_nDauOtherSecondary", SV_matchedGV_nDauOtherSecondary, "Matched GV other-secondary daughter count");
       svGVTable->addColumn<int>("matchedGV_nDauOriginUnknown", SV_matchedGV_nDauOriginUnknown, "Matched GV unknown-origin daughter count");
        iEvent.put(std::move(svGVTable), "SVGVMatchTable");

        // NEW: match GVDaughters (gen daughters) to reconstructed tracks
        auto trkMatchResult = matchDaughtersToTracks(
            Daughters_pt, Daughters_eta, Daughters_phi, Daughters_charge, tracks,
            trkMaxDeltaR_, trkMaxDPtRel_, trkCheckCharge_, trkResolveAmbiguities_);

        std::vector<int>   Daughters_trkIdx       = std::get<0>(trkMatchResult);
        std::vector<int>   Daughters_isTrkMatched = std::get<1>(trkMatchResult);
        std::vector<float> Daughters_trkDeltaR     = std::get<2>(trkMatchResult);
        std::vector<float> Daughters_trkDPtRel     = std::get<3>(trkMatchResult);

        // NEW: per-GV count of daughters matched to a track
        std::vector<int> GV_nDaughtersMatchedToTracks(ngv, 0);
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (!Daughters_isTrkMatched[i]) continue;
            int gvIdx = Daughters_GVidx[i];
            if (gvIdx >= 0 && gvIdx < ngv) {
                GV_nDaughtersMatchedToTracks[gvIdx]++;
            }
        }

        // NEW: per-GV max trkDeltaR among its (matched) daughters (GVDaughters, not GVDirectDaughters).
        // Stays at -1 for a GV that has no track-matched daughters.
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

        //  Build FlatTables 
        auto rejectedGVTable = std::make_unique<nanoaod::FlatTable>(nRejectedGV,"RejectedGV",false);
        rejectedGVTable->addColumn<float>("pt",allHadron_pt,"Rejected Hadron pt");
        rejectedGVTable->addColumn<float>("eta",allHadron_eta,"Rejected Hadron eta");
        rejectedGVTable->addColumn<float>("phi",allHadron_phi,"Rejected Hadron phi");
        rejectedGVTable->addColumn<float>("x",allHadron_GVx,"Rejected GV x");
        rejectedGVTable->addColumn<float>("y",allHadron_GVy,"Rejected GV y");
        rejectedGVTable->addColumn<float>("z",allHadron_GVz,"Rejected GV z");
        rejectedGVTable->addColumn<float>("x_i",allHadron_GVx_i,"Born x coordinate of Rejected GV ");
        rejectedGVTable->addColumn<float>("y_i",allHadron_GVy_i,"Born y coordinate of Rejected GV ");
        rejectedGVTable->addColumn<float>("z_i",allHadron_GVz_i,"Born z coordinate of Rejected GV ");
	rejectedGVTable->addColumn<int>("isB",allHadron_isB,"isB");
        rejectedGVTable->addColumn<int>("isD",allHadron_isD,"isD");
        rejectedGVTable->addColumn<int>("isBtoD",allHadron_isBtoD,"D hadron has a B-hadron ancestor");
	rejectedGVTable->addColumn<int>(
	    "mergedGenPartIdx",
	    allHadron_mergedGenPartIdx,
	    "Index of this hadron in mergedGenParticles"
	);
	
	rejectedGVTable->addColumn<int>(
	    "prunedGenPartIdx",
	    allHadron_prunedGenPartIdx,
	    "Index of the corresponding hadron in prunedGenParticles; -1 if unmatched"
	);
	
	rejectedGVTable->addColumn<int>(
	    "genJetIdx",
	    allHadron_genJetIdx,
	    "Ghost-associated generator-jet index; -1 if absent"
	);
	
	rejectedGVTable->addColumn<int>(
	    "hasGenJet",
	    allHadron_hasGenJet,
	    "One if the rejected GV hadron is associated with a generator jet"
	);
	
	rejectedGVTable->addColumn<float>(
	    "genJetPt",
	    allHadron_genJetPt,
	    "Associated generator-jet pT; -1 if absent"
	);

	rejectedGVTable->addColumn<float>(
	    "genJetEta",
	    allHadron_genJetEta,
	    "Eta of associated generator jet"
	);
	
	rejectedGVTable->addColumn<float>(
	    "genJetPhi",
	    allHadron_genJetPhi,
	    "Phi of associated generator jet"
	);
	
	rejectedGVTable->addColumn<float>(
	    "genJetDeltaR",
	    allHadron_genJetDeltaR,
	    "DeltaR to associated generator jet; -1 if absent"
	);
	
	rejectedGVTable->addColumn<int>(
	    "hasPrunedGenMatch",
	    allHadron_hasPrunedGenMatch,
	    "One if merged rejected hadron was matched to prunedGenParticles"
	);

        auto gvTable = std::make_unique<nanoaod::FlatTable>(ngv,"GV",false);
        gvTable->addColumn<float>("pt",Hadron_pt,"Hadron pt");
        gvTable->addColumn<float>("eta",Hadron_eta,"Hadron eta");
        gvTable->addColumn<float>("phi",Hadron_phi,"Hadron phi");
        gvTable->addColumn<float>("x",Hadron_GVx,"GV x");
        gvTable->addColumn<float>("y",Hadron_GVy,"GV y");
        gvTable->addColumn<float>("z",Hadron_GVz,"GV z");
        gvTable->addColumn<float>("x_i",Hadron_GVx_i,"Born x coordinate of GV ");
        gvTable->addColumn<float>("y_i",Hadron_GVy_i,"Born y coordinate of GV ");
        gvTable->addColumn<float>("z_i",Hadron_GVz_i,"Born z coordinate of GV ");
        gvTable->addColumn<int>("Hadron_SVIdx",Hadron_SVIdx,"SVIdx");
        gvTable->addColumn<int>("Hadron_pdgId",Hadron_pdgId,"Hadron_pdgId");
        gvTable->addColumn<float>("SV_distanceSig",Hadron_SVDistance,"SV_distanceSig");
        // new class
        gvTable->addColumn<int>("isB",Hadron_isB,"isB");
        gvTable->addColumn<int>("isD",Hadron_isD,"isD");
        gvTable->addColumn<int>("isBtoD",Hadron_isBtoD,"D hadron has a B-hadron ancestor");
        gvTable->addColumn<int>("fromHF",Hadron_fromHF,"1 if an ancestor in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
        gvTable->addColumn<int>("toHF",Hadron_toHF,"1 if a descendant in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
        gvTable->addColumn<int>("pdgClass",Hadron_pdgClass,"Hadron class: 1 B, 2 D, 3 strange, 4 tau");
        gvTable->addColumn<float>("minDistNotMatched",Hadron_minDistNotMatched,"Minimum distance to SV among unmatched hadrons");
        gvTable->addColumn<int>("nDaughters",GV_nDaughters,"Total number of GVDaughters belonging to this GV (denominator for nDaughtersMatchedToTracks)"); // NEW
        gvTable->addColumn<float>("maxDaughterPairDeltaR",GV_maxDaughterPairDeltaR,"Max pairwise deltaR among this GV's own GVDaughters (gen-level, no tracks involved)"); // NEW
        gvTable->addColumn<int>("nDaughtersMatchedToTracks",GV_nDaughtersMatchedToTracks,"Number of GVDaughters of this GV matched to a reconstructed track"); // NEW
        gvTable->addColumn<float>("maxDaughterTrkDeltaR",GV_maxDaughterTrkDeltaR,"Max trkDeltaR among this GV's track-matched GVDaughters (not GVDirectDaughters); -1 if none matched");

        // NEW: nearest-SV candidate diagnostics. These are computed from the
        // original distance matrix and do not affect the final greedy match.
        //gvTable->addColumn<int>(
        //    "bestCandidateSVIdx",
        //    GV_bestCandidateSVIdx,
        //    "Index of spatially closest reconstructed SV before one-to-one assignment; -1 if none"
        //);
        //gvTable->addColumn<float>(
        //    "bestCandidateDistanceSig",
        //    GV_bestCandidateDistanceSig,
        //    "Covariance-weighted distance significance to the spatially closest reconstructed SV; -1 if none"
        //);
        //gvTable->addColumn<int>(
        //    "bestCandidateNCommonTracks",
        //    GV_bestCandidateNCommonTracks,
        //    "Number of distinct daughter-track matches between this GV and its spatially closest reconstructed SV"
        //);
        //gvTable->addColumn<int>(
        //    "bestCandidateNMatchedDaughters",
        //    GV_bestCandidateNMatchedDaughters,
        //    "Number of distinct GV daughters matched to tracks in the spatially closest reconstructed SV"
        //);
        //gvTable->addColumn<float>(
        //    "bestCandidateMaxTrackDeltaR",
        //    GV_bestCandidateMaxTrackDeltaR,
        //    "Maximum daughter-track deltaR among diagnostic matches to the spatially closest reconstructed SV; -1 if none"
        //);
        //gvTable->addColumn<float>(
        //    "bestCandidateMaxRelPtDiff",
        //    GV_bestCandidateMaxRelPtDiff,
        //    "Maximum relative pT difference among diagnostic daughter-track matches to the spatially closest reconstructed SV; -1 if none"
        //);
        //gvTable->addColumn<float>(
        //    "bestCandidateDeltaR",
        //    GV_bestCandidateDeltaR,
        //    "DeltaR between GV hadron direction and the spatially closest reconstructed SV direction; -1 if none"
        //);
        //gvTable->addColumn<int>(
        //    "bestCandidateWasFinalMatch",
        //    GV_bestCandidateWasFinalMatch,
        //    "One if the spatially closest reconstructed SV is also the final assigned SV"
        //);
        //gvTable->addColumn<int>("nDauNoRecognizedSecondaryAncestor", Hadron_nDauNoRecognizedSecondaryAncestor, "Selected daughters with no recognized secondary ancestor");
        //gvTable->addColumn<int>("nDauFromB", Hadron_nDauFromB, "Selected daughters with B ancestry");
        //gvTable->addColumn<int>("nDauFromBC", Hadron_nDauFromBC, "Selected daughters with B and C ancestry");
        //gvTable->addColumn<int>("nDauFromC", Hadron_nDauFromC, "Selected daughters with C ancestry");
        //gvTable->addColumn<int>("nDauOtherSecondary", Hadron_nDauOtherSecondary, "Selected daughters from strange/tau/conversion-like ancestry");
        //gvTable->addColumn<int>("nDauOriginUnknown", Hadron_nDauOriginUnknown, "Selected daughters with unknown origin");	
	gvTable->addColumn<int>(
	    "mergedGenPartIdx",
	    Hadron_mergedGenPartIdx,
	    "Index of this GV hadron in mergedGenParticles"
	);
	
	gvTable->addColumn<int>(
	    "prunedGenPartIdx",
	    Hadron_prunedGenPartIdx,
	    "Index of the corresponding hadron in prunedGenParticles; -1 if unmatched"
	);
	
	gvTable->addColumn<int>(
	    "genJetIdx",
	    Hadron_genJetIdx,
	    "Index of the ghost-associated generator jet; -1 if no associated jet"
	);
	
	gvTable->addColumn<int>(
	    "hasGenJet",
	    Hadron_hasGenJet,
	    "One if this hadron is ghost-associated with a generator jet"
	);
	
	gvTable->addColumn<float>(
	    "genJetPt",
	    Hadron_genJetPt,
	    "pT of the ghost-associated generator jet; -1 if absent"
	);
	
	gvTable->addColumn<float>(
	    "genJetEta",
	    Hadron_genJetEta,
	    "Eta of the ghost-associated generator jet"
	);
	
	gvTable->addColumn<float>(
	    "genJetPhi",
	    Hadron_genJetPhi,
	    "Phi of the ghost-associated generator jet"
	);
	
	gvTable->addColumn<float>(
	    "genJetMass",
	    Hadron_genJetMass,
	    "Mass of the ghost-associated generator jet; -1 if absent"
	);
	
	gvTable->addColumn<float>(
	    "genJetDeltaR",
	    Hadron_genJetDeltaR,
	    "DeltaR between the GV hadron and its ghost-associated generator jet"
	);
	
	gvTable->addColumn<int>(
	    "genJetHadronFlavour",
	    Hadron_genJetHadronFlavour,
	    "Hadron flavour assigned to the associated generator jet"
	);
	
	gvTable->addColumn<int>(
	    "genJetPartonFlavour",
	    Hadron_genJetPartonFlavour,
	    "Parton flavour assigned to the associated generator jet"
	);
	
	gvTable->addColumn<int>(
	    "genJetNCHadrons",
	    Hadron_genJetNCHadrons,
	    "Number of ghost-associated charm hadrons in the generator jet"
	);
	
	gvTable->addColumn<int>(
	    "genJetNBHadrons",
	    Hadron_genJetNBHadrons,
	    "Number of ghost-associated bottom hadrons in the generator jet"
	);
	gvTable->addColumn<int>(
	    "hasPrunedGenMatch",
	    Hadron_hasPrunedGenMatch,
	    "One if merged GV hadron was matched to prunedGenParticles"
	);
        
        //

        auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(),"GVDaughters",false);
        dauTable->addColumn<float>("pt",Daughters_pt,"Daughter pt");
        dauTable->addColumn<float>("eta",Daughters_eta,"Daughter eta");
        dauTable->addColumn<float>("phi",Daughters_phi,"Daughter phi");
        dauTable->addColumn<float>("vx", Daughters_vx, "Daughters_vx");
        dauTable->addColumn<float>("vy", Daughters_vy, "Daughters_vy");
        dauTable->addColumn<float>("vz", Daughters_vz, "Daughters_vz");
        dauTable->addColumn<int>("charge",Daughters_charge,"Daughter charge");
        dauTable->addColumn<int>("pdgId",Daughters_pdgId,"Daughter pdgId");
        dauTable->addColumn<int>("hadronIndex",Daughters_GVidx,"Hadron index");
        dauTable->addColumn<int>("originLabel", Daughters_originLabel, "0 no recognized secondary ancestor, 2 fromB, 3 fromBC, 4 fromC, 5 otherSecondary, 9 unknown");
        dauTable->addColumn<int>("trkIdx",Daughters_trkIdx,"Index of matched track in the input track collection, -1 if unmatched");     // NEW
        dauTable->addColumn<int>("isTrkMatched",Daughters_isTrkMatched,"1 if daughter matched to a reconstructed track");                // NEW
        dauTable->addColumn<float>("trkDeltaR",Daughters_trkDeltaR,"deltaR to matched track, -1 if unmatched");                          // NEW
        dauTable->addColumn<float>("trkDPtRel",Daughters_trkDPtRel,"relative pT difference to matched track, -1 if unmatched");          // NEW

        auto directdauTable = std::make_unique<nanoaod::FlatTable>(directDaughters_pt.size(),"GVDirectDaughters",false);
        directdauTable->addColumn<float>("pt",directDaughters_pt,"Daughter pt");
        directdauTable->addColumn<float>("eta",directDaughters_eta,"Daughter eta");
        directdauTable->addColumn<float>("phi",directDaughters_phi,"Daughter phi");
        directdauTable->addColumn<int>("charge",directDaughters_charge,"Daughter charge");
        directdauTable->addColumn<int>("pdgId",directDaughters_pdgId,"Daughter pdgId");
        directdauTable->addColumn<int>("hadronIndex",directDaughters_GVidx,"Hadron index");


        //dauTable->addColumn<int>("hadronFlav",Daughters_flav,"Hadron flavor");


        //
        iEvent.put(std::move(gvTable),"GVTable");
        iEvent.put(std::move(rejectedGVTable),"rejectedGVTable");
        iEvent.put(std::move(dauTable),"GVDaughtersTable");
        iEvent.put(std::move(directdauTable),"GVDirectDaughters");

    }






//  checkPDG() 
int GenVertexProducer::checkPDG(int abs_pdg) const {
    std::vector<int> pdgList_B = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    std::vector<int> pdgList_D = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    std::vector<int> pdgList_S = {3122,3222,3212,3312,3322,3334};
    std::vector<int> pdgList_Tau = {15};

    if(std::find(pdgList_B.begin(),pdgList_B.end(),abs_pdg)!=pdgList_B.end()) return 1;
    if(std::find(pdgList_D.begin(),pdgList_D.end(),abs_pdg)!=pdgList_D.end()) return 2;
    if(std::find(pdgList_S.begin(),pdgList_S.end(),abs_pdg)!=pdgList_S.end()) return 3;
    if(std::find(pdgList_Tau.begin(),pdgList_Tau.end(),abs_pdg)!=pdgList_Tau.end()) return 4;
    return 0;
}


bool GenVertexProducer::hasBHadronAncestor(const reco::Candidate* hadron) const {
    if (hadron == nullptr) return false;
    const reco::Candidate* current = hadron;
    int guard = 0;
    while (current != nullptr && current->numberOfMothers() > 0 && guard++ < 100) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == nullptr || mother == current) break;
        if (checkPDG(std::abs(mother->pdgId())) == 1) return true;
        current = mother;
    }
    return false;
}

int GenVertexProducer::getDaughterOriginLabelNoPU(const reco::Candidate* daughter) const {
    if (daughter == nullptr) return 9;
    static const std::unordered_set<int> pdgSetB = {
        521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    static const std::unordered_set<int> pdgSetC = {
        411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    static const std::unordered_set<int> pdgSetOther = {
        310,130,3122,3222,3212,3312,3322,3334};
    bool foundB = false, foundC = false, foundOther = false;
    const reco::Candidate* current = daughter;
    int guard = 0;
    while (current != nullptr && current->numberOfMothers() > 0 && guard++ < 100) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == nullptr || mother == current) break;
        const int pdg = std::abs(mother->pdgId());
        foundB |= pdgSetB.count(pdg);
        foundC |= pdgSetC.count(pdg);
        foundOther |= pdgSetOther.count(pdg) || pdg == 15 || pdg == 22;
        current = mother;
    }
    if (foundB && foundC) return 3;
    if (foundB) return 2;
    if (foundC) return 4;
    if (foundOther) return 5;
    return 0;
}

//  hasHFAncestor() 
// Walks upstream from `hadron` (excluding itself). Returns true if any
// ancestor along the mother(0) chain is itself an HF/long-lived particle
// (B, D, S baryon, or Tau, per the same PDG lists as checkPDG()).
bool GenVertexProducer::hasHFAncestor(const reco::Candidate* hadron) const {
    const reco::Candidate* current = hadron;
    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == nullptr || mother == current) break;
        if (checkPDG(std::abs(mother->pdgId())) != 0) {
            return true;
        }
        current = mother;
    }
    return false;
}

//  hasHFDescendant() 
// Walks the full daughter tree of `hadron` (all daughters, not just stable
// charged ones). Returns true if any descendant at any depth is itself an
// HF/long-lived particle (B, D, S baryon, or Tau).
bool GenVertexProducer::hasHFDescendant(const reco::Candidate* hadron) const {
    for (size_t i = 0; i < hadron->numberOfDaughters(); ++i) {
        const reco::Candidate* dau = hadron->daughter(i);
        if (checkPDG(std::abs(dau->pdgId())) != 0) {
            return true;
        }
        if (hasHFDescendant(dau)) {
            return true;
        }
    }
    return false;
}

//  isAncestor() 
std::optional<std::tuple<float, float, float>> GenVertexProducer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle) const
    {
    std::vector<int> pdgList_B = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    std::vector<int> pdgList_D = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    std::vector<int> pdgList_S = {3122,3222,3212,3312,3322,3334};
    std::vector<int> pdgList_Tau = {15};
    std::unordered_set<int> pdgSet_D(pdgList_D.begin(), pdgList_D.end());
    std::unordered_set<int> pdgSet_B(pdgList_B.begin(), pdgList_B.end());
    std::unordered_set<int> pdgSet_S(pdgList_S.begin(), pdgList_S.end());
    std::unordered_set<int> pdgSet_Tau(pdgList_Tau.begin(), pdgList_Tau.end());
    const reco::Candidate* current = particle;
    //const reco::Candidate* child = nullptr;

    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == ancestor) {
            // Found the ancestor; return the vertex of the current particle (i.e., the direct daughter)
            return std::make_optional(std::make_tuple(current->vx(), current->vy(), current->vz()));
        }
        int mother_pdg = std::abs(mother->pdgId());
        if (pdgSet_B.count(mother_pdg) || pdgSet_D.count(mother_pdg) || pdgSet_S.count(mother_pdg) || pdgSet_Tau.count(mother_pdg)) break;
        current = mother;
    }

    // If we reached here, the ancestor was not found in the chain
    return std::nullopt;
}


std::vector<std::vector<float>> GenVertexProducer::computeDistanceMatrix(
                const std::vector<float>& SV_x,
                const std::vector<float>& SV_y,
                const std::vector<float>& SV_z,
                std::vector<CovMatrix> SV_cov,
                const std::vector<float>& Hadron_GVx,
                const std::vector<float>& Hadron_GVy,
                const std::vector<float>& Hadron_GVz) {
    // computeDistanceMatrix
    // Returns
    // distances = Matrix of Euclidean distances between SV and GV
    
    
    size_t nSV = SV_x.size();
    size_t nHadron = Hadron_GVx.size();
    
    // 2D vector initialized to 999 nSV x nHadron
    std::vector<std::vector<float>> distances(nSV, std::vector<float>(nHadron, 999.0));

    for (size_t i = 0; i < nSV; ++i) {
        //CovMatrix covInv = SV_cov[i].Inverse();
        CovMatrix covInv = SV_cov[i];
        covInv.Invert();
        for (size_t j = 0; j < nHadron; ++j) {
            float dx = SV_x[i] - Hadron_GVx[j];
            float dy = SV_y[i] - Hadron_GVy[j];
            float dz = SV_z[i] - Hadron_GVz[j];
            float chi2 =
              dx * (covInv(0,0) * dx +
                    covInv(0,1) * dy +
                    covInv(0,2) * dz)

            + dy * (covInv(1,0) * dx +
                    covInv(1,1) * dy +
                    covInv(1,2) * dz)

            + dz * (covInv(2,0) * dx +
                    covInv(2,1) * dy +
                    covInv(2,2) * dz);

            float dist = 999.0f;
            if (std::isfinite(chi2) && chi2 >= 0.0f) dist = std::sqrt(chi2);
            distances[i][j] = dist;
        }
    }

    return distances;
}


void GenVertexProducer::printDistanceMatrix(
    const std::vector<std::vector<float>>& distances
) {
    size_t nSV = distances.size();
    if (nSV == 0) return;

    size_t nGV = distances[0].size();

    std::cout << "\n Distance Matrix (SV rows × GV cols):\n\n";

    // Print header row
    std::cout << std::setw(8) << "SV/GV";
    for (size_t j = 0; j < nGV; ++j) {
        std::cout << std::setw(10) << "GV[" + std::to_string(j) + "]";
    }
    std::cout << "\n";

    // Print matrix values
    for (size_t i = 0; i < nSV; ++i) {
        std::cout << std::setw(8) << "SV[" + std::to_string(i) + "]";
        for (size_t j = 0; j < nGV; ++j) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(3) << distances[i][j];
        }
        std::cout << "\n";
    }
}
float computeDR_SV_Had(int bestSV,
                       int bestHad,
                       const std::vector<float>& SV_eta,
                       const std::vector<float>& SV_phi,
                       const std::vector<float>& GV_eta,
                       const std::vector<float>& GV_phi)
{
    float dEta = SV_eta[bestSV] - GV_eta[bestHad];

    float dPhi = SV_phi[bestSV] - GV_phi[bestHad];

    // wrap phi into [-pi, pi]
    while (dPhi > M_PI)  dPhi -= 2.0 * M_PI;
    while (dPhi < -M_PI) dPhi += 2.0 * M_PI;

    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>, std::vector<int>>  GenVertexProducer::matchHadronsToSV(
    std::vector<std::vector<float>> distances,
    const std::vector<float>& SVtrk_pt,
    const std::vector<float>& SVtrk_eta,
    const std::vector<float>& SVtrk_phi,
    const std::vector<int>& SVtrk_SVidx,
    const std::vector<float>& Daughters_pt,     //genparticles
    const std::vector<float>& Daughters_eta,  //genparticles
    const std::vector<float>& Daughters_phi,  //genparticles
    const std::vector<int>& Daughters_GVidx, // hadron index per daughter
    const std::vector<float>& SV_eta,
    const std::vector<float>& SV_phi,
    const std::vector<float>& GV_eta,
    const std::vector<float>& GV_phi,
    int n_Hadrons,
    int nRequiredCommonTracks,
    double dR_max,
    double relPt_max,
    bool doubleMatching,
    int doubleMatching_nRequiredCommonTracks,
    double doubleMatching_maxSignificance,
    double doubleMatching_dR_max,
    double doubleMatching_relPt_max
) {
    // Returns
    // Hadron_SVIdx = Array of lenght = Hadron_pt.size() with index of the SV
    size_t nSV = distances.size(); //distances is nSV x nHadron matrix (the first dimension is nSV)
    std::vector<int> Hadron_SVIdx(n_Hadrons, -1); // Output
    std::vector<float> Hadron_SVDistance(n_Hadrons, -1); // Output
    std::vector<int> SVtrk_isMatched(SVtrk_pt.size(), 0);   // NEW
    std::vector<int> SVtrk_GVIdx(SVtrk_pt.size(), -1);      // NEW
    std::vector<int> SVtrk_daughterIdx(SVtrk_pt.size(), -1);
    std::vector<size_t> svTrackIdxs_fromBestSV;
    std::vector<std::vector<float>> distancesOriginal = distances;
    while (true) {
        float minDist = 999.0;
        int bestSV = -1;
        int bestHad = -1;

        // Find minimum distance in current matrix
        // store in 
        // - minDist
        // - bestSV
        // - bestHad
        for (size_t sv = 0; sv < nSV; ++sv) {
            for (int had = 0; had < n_Hadrons; ++had) {
                if (distances[sv][had] < minDist) {
                    minDist = distances[sv][had];
                    bestSV = sv;
                    bestHad = had;
                }
            }
        }
        //std::cout<<"\n Considering best pair: SV["<<bestSV<<"] and Hadron["<<bestHad<<"] with distance "<<minDist<<std::endl;
        if (minDist >= 997.0) break;  // done

        // Select tracks from SV
        // svTrackIdxs_fromBestSV is initialized every time
        // it stores the index of the tracks which originate from SV candidate in this loop
        
        svTrackIdxs_fromBestSV.clear();
        for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
            // among all tracks from all SV, select those from the candidate SV
            //std::cout<<" Track index "<<i<<" SVtrk_SVidx: "<<SVtrk_SVidx[i]<<" SVtrk_pt: "<<SVtrk_pt[i]<<" Best SV :"<<bestSV<<std::endl;
            if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > 0.4 && std::fabs(SVtrk_eta[i]) < 2.5) {
                svTrackIdxs_fromBestSV.push_back(i);
            }
        }

        // Select daughters of Hadron
        std::vector<size_t> GenDaughtersIdxs_fromBestHad;
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (Daughters_GVidx[i] == bestHad) {
                GenDaughtersIdxs_fromBestHad.push_back(i);
            }
        }

        // Match logic: check for 1 (2) or more matched tracks by ΔR & dPt/pT
        //int nRequiredCommonTracks = 1;
        int common = 0;
        //std::vector<size_t> matchedTrackIdxs; // global SVtrk indices matched in this SV/Hadron trial
        std::vector<std::pair<size_t,size_t>> matchedTrackToDaughter; // 
        for (size_t iSV : svTrackIdxs_fromBestSV) {
            bool trackMatched = false;
            size_t matchedDauIdx = 0;
            //std::cout<<" Checking SVtrack index "<<iSV<<std::endl;
            for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                //std::cout<<" Checking GVDaughters index "<<iHad<<std::endl;
                float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / std::max(Daughters_pt[iHad], 1e-6f);
                //std::cout<<"Comparing SV track (pt: "<<SVtrk_pt[iSV]<<", eta: "<<SVtrk_eta[iSV]<<", phi: "<<SVtrk_phi[iSV]<<") with Daughter (pt: "<<Daughters_pt[iHad]<<", eta: "<<Daughters_eta[iHad]<<", phi: "<<Daughters_phi[iHad]<<") => dR: "<<dR<<", relPt: "<<relPt<<std::endl;
                if (dR < dR_max && relPt < relPt_max) {
                    
                    trackMatched = true;
                    matchedDauIdx = iHad;
                    break;
                    //std::cout<<"  -> Matched! Common tracks: "<<common<<std::endl;
                    //if (common >= nRequiredCommonTracks) break; // break the iHad cycle
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




    //printDistanceMatrix(distances);
    if (doubleMatching){

        while (true){
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
            //check whether deltaR bewteen bestSV and bestHad is less than doublematching_dR_max
            float dR_SV_Had = computeDR_SV_Had(bestSV, bestHad, SV_eta, SV_phi, GV_eta, GV_phi);
            

            svTrackIdxs_fromBestSV.clear();
            for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
                // among all tracks from all SV, select those from the candidate SV
                //std::cout<<" Track index "<<i<<" SVtrk_SVidx: "<<SVtrk_SVidx[i]<<" SVtrk_pt: "<<SVtrk_pt[i]<<" Best SV :"<<bestSV<<std::endl;
                if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > 0.4 && std::fabs(SVtrk_eta[i]) < 2.5) {
                    svTrackIdxs_fromBestSV.push_back(i);
                }
            }

            // Select daughters of Hadron
            std::vector<size_t> GenDaughtersIdxs_fromBestHad;
            for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
                if (Daughters_GVidx[i] == bestHad) {
                    GenDaughtersIdxs_fromBestHad.push_back(i);
                }
            }

            // Match logic: check for 1 (2) or more matched tracks by ΔR & dPt/pT
            //int nRequiredCommonTracks = 1;
            int common = 0;
            std::vector<std::pair<size_t,size_t>> matchedTrackToDaughter;
            for (size_t iSV : svTrackIdxs_fromBestSV) {
                for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                    const float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                    const float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) /
                                        std::max(Daughters_pt[iHad], 1e-6f);
                    if (dR < doubleMatching_dR_max && relPt < doubleMatching_relPt_max) {
                        ++common;
                        matchedTrackToDaughter.push_back({iSV, iHad});
                        break;
                    }
                }
                if (common >= doubleMatching_nRequiredCommonTracks) break;
            }
            //std::cout<<"\n [DoubleMatching] Pair: SV["<<bestSV<<"] and Hadron["<<bestHad<<"] with distance "<<distancesOriginal[bestSV][bestHad]<<" and dR "<<dR_SV_Had<<" and common tracks "<<common<<std::endl;
            if (common >= doubleMatching_nRequiredCommonTracks && distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance ) {
            //if (common >= 1 && distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance && dR_SV_Had < doubleMatching_dR_max ) {

                Hadron_SVIdx[bestHad] = bestSV;
                Hadron_SVDistance[bestHad]= -distancesOriginal[bestSV][bestHad];
                for (const auto& match : matchedTrackToDaughter) {
                    SVtrk_isMatched[match.first] = 1;
                    SVtrk_GVIdx[match.first] = bestHad;
                    SVtrk_daughterIdx[match.first] = static_cast<int>(match.second);
                }
                for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 1000.0;  // remove SV row (MATCHED)
                for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 1000.0;   // remove Hadron column (MATCHED)
                //std::cout << "[V] Matched Hadron[" << bestHad << "] to SV[" << bestSV << "] (distance = " << distancesOriginal[bestSV][bestHad] << ", common tracks = " << common << ")\n";
            } else {
                distances[bestSV][bestHad] = 999.0;  // exclude this pair
            }



        }

    }

    // check the minDist not mathced using original matrix (not modified by matching procedure)
    std::vector<float> minDistNotMatched(n_Hadrons, 999.);
    for (int had = 0; had < n_Hadrons; ++had) {
        for (size_t sv = 0; sv < nSV; ++sv) {
            if (distancesOriginal[sv][had] < minDistNotMatched[had] && Hadron_SVIdx[had] == -1) {
                minDistNotMatched[had] = distancesOriginal[sv][had];
            }
        }
    }
    return std::make_tuple(Hadron_SVIdx, Hadron_SVDistance, minDistNotMatched, SVtrk_isMatched, SVtrk_GVIdx, SVtrk_daughterIdx);
    }


// NEW: matchDaughtersToTracks()
// Matches gen-level GVDaughters to reconstructed tracks by deltaR and relative pT,
// following the same TrackGenMatcher-style logic (optional charge check, optional
// one-to-one "resolveAmbiguities" assignment).
//
// Strategy: build the list of all (daughter, track) pairs that pass the deltaR/relPt
// (and, if requested, charge) cuts, sort them by ascending deltaR, then greedily assign
// pairs. When resolveAmbiguities is true, once a track is used it can't be reused by
// another daughter (one-to-one matching, analogous to TrackGenMatcher's resolveAmbiguities).
// When false, several daughters may be matched to the same track.
std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>>
GenVertexProducer::matchDaughtersToTracks(
    const std::vector<float>& Daughters_pt,
    const std::vector<float>& Daughters_eta,
    const std::vector<float>& Daughters_phi,
    const std::vector<int>& Daughters_charge,
    const std::vector<reco::Track>& tracks,
    double maxDeltaR,
    double maxDPtRel,
    bool checkCharge,
    bool resolveAmbiguities) const
{
    size_t nDau = Daughters_pt.size();
    size_t nTrk = tracks.size();

    std::vector<int>   trkIdx(nDau, -1);
    std::vector<int>   isMatched(nDau, 0);
    std::vector<float> matchDeltaR(nDau, -1.f);
    std::vector<float> matchDPtRel(nDau, -1.f);

    struct Pair { float dR; size_t dau; size_t trk; };
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

    // best (smallest deltaR) candidates get assigned first
    std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.dR < b.dR; });

    std::vector<bool> dauUsed(nDau, false);
    std::vector<bool> trkUsed(nTrk, false);
    for (const auto& p : pairs) {
        if (dauUsed[p.dau]) continue;                        // each daughter is matched at most once
        if (resolveAmbiguities && trkUsed[p.trk]) continue;   // one-to-one if requested

        trkIdx[p.dau]      = static_cast<int>(p.trk);
        isMatched[p.dau]   = 1;
        matchDeltaR[p.dau] = p.dR;
        matchDPtRel[p.dau] = std::fabs(tracks[p.trk].pt() - Daughters_pt[p.dau]) / Daughters_pt[p.dau];

        dauUsed[p.dau] = true;
        if (resolveAmbiguities) trkUsed[p.trk] = true;
    }

    return std::make_tuple(trkIdx, isMatched, matchDeltaR, matchDPtRel);
}

int GenVertexProducer::findMatchingPrunedHadron(
    const reco::Candidate* mergedHadron,
    const reco::GenParticleCollection& prunedParticles
) const {
    if (mergedHadron == nullptr) {
        return -1;
    }

    int bestIdx = -1;
    float bestScore = std::numeric_limits<float>::max();

    for (size_t i = 0; i < prunedParticles.size(); ++i) {
        const reco::GenParticle& pruned = prunedParticles[i];

        if (pruned.pdgId() != mergedHadron->pdgId()) {
            continue;
        }

        if (pruned.status() != mergedHadron->status()) {
            continue;
        }

        const float dR = deltaR(
            static_cast<float>(mergedHadron->eta()),
            static_cast<float>(mergedHadron->phi()),
            static_cast<float>(pruned.eta()),
            static_cast<float>(pruned.phi())
        );

        const float relPt =
            std::abs(
                static_cast<float>(mergedHadron->pt()) -
                static_cast<float>(pruned.pt())
            ) /
            std::max(
                static_cast<float>(mergedHadron->pt()),
                1.0e-6f
            );

        const float dx =
            static_cast<float>(mergedHadron->vx() - pruned.vx());
        const float dy =
            static_cast<float>(mergedHadron->vy() - pruned.vy());
        const float dz =
            static_cast<float>(mergedHadron->vz() - pruned.vz());

        const float vertexDistance =
            std::sqrt(dx * dx + dy * dy + dz * dz);

        // MergedGenParticleProducer copies pruned particles, so these
        // differences should normally be extremely small or exactly zero.
        if (dR > 1.0e-5f) {
            continue;
        }

        if (relPt > 1.0e-5f) {
            continue;
        }

        if (vertexDistance > 1.0e-5f) {
            continue;
        }

        const float score = dR + relPt + vertexDistance;

        if (score < bestScore) {
            bestScore = score;
            bestIdx = static_cast<int>(i);
        }
    }

    return bestIdx;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenVertexProducer);
