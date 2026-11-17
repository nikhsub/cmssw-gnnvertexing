#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

class SVTableProducer : public edm::global::EDProducer<> {
public:
    explicit SVTableProducer(const edm::ParameterSet &iConfig);
    void produce(edm::StreamID,
                 edm::Event &iEvent,
                 const edm::EventSetup &iSetup) const override;

private:
    edm::EDGetTokenT<std::vector<reco::Vertex>> svToken;
    edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> sourceSVToken_;
    edm::EDGetTokenT<reco::TrackCollection> unpackedTracksToken_;
    edm::EDGetTokenT<edm::ValueMap<int>> packedTrackIndexToken_;
    edm::EDGetTokenT<edm::ValueMap<int>> lostTrackIndexToken_;
    edm::EDGetTokenT<edm::ValueMap<int>> filteredTrackIndexToken_;
    bool hasSourceConstituents_;
    bool hasFilteredTrackIndex_;
    double dlenSigMin_;  
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
    //edm::EDGetTokenT<edm::ValueMap<float>> svscoreToken_;
};

SVTableProducer::SVTableProducer(const edm::ParameterSet &iConfig): 
    svToken(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    hasSourceConstituents_(iConfig.existsAs<edm::InputTag>("sourceSecondaryVertices")),
    hasFilteredTrackIndex_(iConfig.existsAs<edm::InputTag>("filteredTrackIndex")),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
    //svscoreToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pvSrc")))
{
    if (hasSourceConstituents_) {
        sourceSVToken_ = consumes<std::vector<reco::VertexCompositePtrCandidate>>(
            iConfig.getParameter<edm::InputTag>("sourceSecondaryVertices"));
        unpackedTracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("unpackedTracks"));
        packedTrackIndexToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("packedTrackIndex"));
        lostTrackIndexToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("lostTrackIndex"));
    }
    if (hasFilteredTrackIndex_)
        filteredTrackIndexToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("filteredTrackIndex"));
    produces<nanoaod::FlatTable>("SVTable");
    produces<nanoaod::FlatTable>("SVtrksTable");
    produces<nanoaod::FlatTable>("PVTable");   // <-- new
}

void SVTableProducer::produce(edm::StreamID,
                 edm::Event &iEvent,
                 const edm::EventSetup &iSetup) const 
    {
        edm::Handle<std::vector<reco::Vertex>> svs;
        iEvent.getByToken(svToken, svs);
        auto pvsIn = iEvent.getHandle(pvs_);
        edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> sourceSVs;
        edm::Handle<reco::TrackCollection> unpackedTracks;
        edm::Handle<edm::ValueMap<int>> packedTrackIndex, lostTrackIndex;
        if (hasSourceConstituents_) {
            iEvent.getByToken(sourceSVToken_, sourceSVs);
            iEvent.getByToken(unpackedTracksToken_, unpackedTracks);
            iEvent.getByToken(packedTrackIndexToken_, packedTrackIndex);
            iEvent.getByToken(lostTrackIndexToken_, lostTrackIndex);
            if (sourceSVs->size() != svs->size())
                throw cms::Exception("MismatchedSVCollections")
                    << "sourceSecondaryVertices and reco SV collections have different sizes";
        }
        edm::Handle<edm::ValueMap<int>> filteredTrackIndex;
        if (hasFilteredTrackIndex_)
            iEvent.getByToken(filteredTrackIndexToken_, filteredTrackIndex);
        VertexDistance3D vdist;
        //edm::Handle<edm::ValueMap<float>> svscoreHandle;
        //iEvent.getByToken(svscoreToken_, svscoreHandle);

        //std::cout<<"[DEBUG] SVTableProducer " << "Number of SVs: " << svs->size()<< "\n";
        //for (size_t i = 0; i < svs->size(); ++i) {
        //    const auto &sv = svs->at(i);
        //    //std::cout<<"SVTableProducer\n" 
        //    //    << "SV " << i 
        //    //    << ": x=" << sv.x() 
        //    //    << " y=" << sv.y() 
        //    //    << " z=" << sv.z()
        //    //    << " nTracks=" << sv.tracksSize()<<"\n";
        //}

        const auto& PV0 = pvsIn->front();
        // --- PV table (leading vertex only) ---
        auto pv_table = std::make_unique<nanoaod::FlatTable>(1, "myPV", true); // true = singleton

        std::vector<float> pv_x{static_cast<float>(PV0.x())};
        std::vector<float> pv_y{static_cast<float>(PV0.y())};
        std::vector<float> pv_z{static_cast<float>(PV0.z())};
        std::vector<float> pv_chi2{static_cast<float>(PV0.chi2())};
        std::vector<float> pv_ndof{static_cast<float>(PV0.ndof())};
        std::vector<float> pv_rho{static_cast<float>(PV0.position().Rho())};
        std::vector<int> pv_ntracks{static_cast<int>(PV0.tracksSize())};
        std::vector<bool> pv_isFake{PV0.isFake()};

        pv_table->addColumn<float>("x", pv_x, "X position of leading PV");
        pv_table->addColumn<float>("y", pv_y, "Y position of leading PV");
        pv_table->addColumn<float>("z", pv_z, "Z position of leading PV");
        pv_table->addColumn<float>("chi2", pv_chi2, "Chi2 of leading PV fit");
        pv_table->addColumn<float>("ndof", pv_ndof, "Degrees of freedom of leading PV fit");
        pv_table->addColumn<float>("rho", pv_rho, "Transverse position of leading PV");
        pv_table->addColumn<int>("nTracks", pv_ntracks, "Number of tracks in leading PV");
        pv_table->addColumn<bool>("isFake", pv_isFake, "Is fake PV flag");

        iEvent.put(std::move(pv_table), "PVTable");
        unsigned int nSVtracks = 0;
        int nTrksCurrentSV = 0;
        unsigned int nSV_cutdlen = 0;
        for (size_t isv = 0; isv < svs->size(); ++isv) {
            const auto& sv = (*svs)[isv];
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                nTrksCurrentSV = 0;
                if (hasSourceConstituents_) {
                    nTrksCurrentSV = (*sourceSVs)[isv].numberOfDaughters();
                } else {
                    for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                        const edm::RefToBase<reco::Track>& trkRef = *it;
                        if (trkRef.isNull()) continue;
                        nTrksCurrentSV++;
                    }
                }
                if (nTrksCurrentSV < 2) continue; // skip SVs with less than 2 high-weight tracks
                nSVtracks += nTrksCurrentSV;
                nSV_cutdlen+=1;
            }
        }
        auto table = std::make_unique<nanoaod::FlatTable>(static_cast<unsigned int>(nSV_cutdlen), "mySV", false);
        auto trk_table = std::make_unique<nanoaod::FlatTable>(static_cast<unsigned int>(nSVtracks), "mySVtrks", false);


        std::vector<float> x, y, z, chi2, ndof, pt, eta, phi, mass, dlen, dlenSig;
        std::vector<float> trk_pt, trk_eta, trk_phi, trk_weight;// trk_SVscore;
        std::vector<int> trk_trackIdx, trk_sourceCollectionId, trk_sourceIdx;
        std::vector<float> trk_ip_z, trk_ip_z_sig, trk_ip2d, trk_ip3d, trk_ip2d_sig, trk_ip3d_sig, trk_p, trk_charge, trk_numberOfValidHits, trk_numberOfValidPixelHits, trk_numberOfValidStripHits;
        // vector pair [num_tracks * num_tracks / 2]
        // trk_i 
        // trk_j 
        // deltaR
        // dca
        // dca sig 
        // pvtoPCA_i
        // pvtoPCA_j
        // dotprod_i
        // dotprod_j
        // pair_mom
        // pair_invmass

        std::vector<int> nTracks, nMappedTracks;
        std::vector<int> trk_SVidx;

        int nTrksPerSV = 0;
        //std::cout<<svs->size()<<" SVs to process\n";
        for (size_t isv = 0; isv < svs->size(); ++isv) {
            const auto& sv = (*svs)[isv];
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                if (hasSourceConstituents_) {
                    nTrksPerSV = (*sourceSVs)[isv].numberOfDaughters();
                } else {
                    nTrksPerSV = 0;
                    for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                        const edm::RefToBase<reco::Track>& trkRef = *it;
                        if (!trkRef.isNull()) ++nTrksPerSV;
                    }
                }
                if (nTrksPerSV < 2) continue;

                x.push_back(sv.x());
                y.push_back(sv.y());
                z.push_back(sv.z());
                chi2.push_back(sv.chi2());
                ndof.push_back(sv.ndof());
                nTracks.push_back(nTrksPerSV);
                int nMappedThisSV = 0;
                dlen.push_back(dl.value());
                dlenSig.push_back(dl.significance());
                //std::cout<<"dl significance is "<< dl.significance()<<std::endl;


                TLorentzVector p4s_SV;
                if (hasSourceConstituents_) {
                    const auto& sourceSV = (*sourceSVs)[isv];
                    for (size_t idau = 0; idau < sourceSV.numberOfDaughters(); ++idau) {
                        const auto& daughterPtr = sourceSV.daughterPtr(idau);
                        const reco::Candidate* daughter = daughterPtr.isNonnull() ? daughterPtr.get() : nullptr;
                        const double dauPt = daughter ? daughter->pt() : -1.0;
                        const double dauEta = daughter ? daughter->eta() : -1.0;
                        const double dauPhi = daughter ? daughter->phi() : -1.0;

                        int sourceCollectionId = -1;
                        int sourceIndex = daughterPtr.isNonnull() ? static_cast<int>(daughterPtr.key()) : -1;
                        int unpackedIndex = -1;
                        if (daughterPtr.isNonnull() && packedTrackIndex.isValid() &&
                            packedTrackIndex->contains(daughterPtr.id())) {
                            sourceCollectionId = 0;  // packedPFCandidates
                            unpackedIndex = (*packedTrackIndex)[daughterPtr];
                        } else if (daughterPtr.isNonnull() && lostTrackIndex.isValid() &&
                                   lostTrackIndex->contains(daughterPtr.id())) {
                            sourceCollectionId = 1;  // lostTracks
                            unpackedIndex = (*lostTrackIndex)[daughterPtr];
                        }

                        int trackIdx = -1;
                        if (unpackedIndex >= 0 && unpackedTracks.isValid() && filteredTrackIndex.isValid()) {
                            reco::TrackRef unpackedRef(unpackedTracks, static_cast<unsigned int>(unpackedIndex));
                            if (filteredTrackIndex->contains(unpackedRef.id()))
                                trackIdx = (*filteredTrackIndex)[unpackedRef];
                        }
                        if (trackIdx >= 0) ++nMappedThisSV;

                        const auto* packedDaughter = dynamic_cast<const pat::PackedCandidate*>(daughter);
                        const bool hasTrackDetails = packedDaughter && packedDaughter->hasTrackDetails();
                        reco::Track sourceTrack;
                        if (hasTrackDetails) sourceTrack = packedDaughter->pseudoTrack();

                        trk_pt.push_back(dauPt);
                        trk_eta.push_back(dauEta);
                        trk_phi.push_back(dauPhi);
                        trk_weight.push_back(1.0f);
                        trk_trackIdx.push_back(trackIdx);
                        trk_sourceCollectionId.push_back(sourceCollectionId);
                        trk_sourceIdx.push_back(sourceIndex);
                        trk_SVidx.push_back(x.size() - 1);

                        const double momentum = daughter ? daughter->p() : -1.0;
                        const int charge = daughter ? daughter->charge() : 0;
                        trk_p.push_back(momentum);
                        trk_charge.push_back(charge);
                        trk_numberOfValidHits.push_back(hasTrackDetails ? sourceTrack.numberOfValidHits() : -1.0f);
                        trk_numberOfValidPixelHits.push_back(hasTrackDetails ? sourceTrack.hitPattern().numberOfValidPixelHits() : -1.0f);
                        trk_numberOfValidStripHits.push_back(hasTrackDetails ? sourceTrack.hitPattern().numberOfValidStripHits() : -1.0f);

                        if (daughter) {
                            TLorentzVector p4;
                            p4.SetPtEtaPhiM(dauPt, dauEta, dauPhi, 0.13957039);
                            p4s_SV += p4;
                        }

                        if (hasTrackDetails) {
                            const double ipZ = sourceTrack.dz(PV0.position());
                            trk_ip_z.push_back(ipZ);
                            trk_ip_z_sig.push_back(sourceTrack.dzError() > 0.0 ? ipZ / sourceTrack.dzError() : -1.0);

                            GlobalVector direction(1, 0, 0);
                            const auto& ttBuilder = iSetup.getData(theTTBToken);
                            reco::TransientTrack transientTrack = ttBuilder.build(sourceTrack);
                            auto ip2d = IPTools::signedTransverseImpactParameter(transientTrack, direction, PV0).second;
                            auto ip3d = IPTools::signedImpactParameter3D(transientTrack, direction, PV0).second;
                            trk_ip2d.push_back(ip2d.value());
                            trk_ip3d.push_back(ip3d.value());
                            trk_ip2d_sig.push_back(ip2d.significance());
                            trk_ip3d_sig.push_back(ip3d.significance());
                        } else {
                            trk_ip_z.push_back(-1.0f);
                            trk_ip_z_sig.push_back(-1.0f);
                            trk_ip2d.push_back(-1.0f);
                            trk_ip3d.push_back(-1.0f);
                            trk_ip2d_sig.push_back(-1.0f);
                            trk_ip3d_sig.push_back(-1.0f);
                        }
                    }
                } else {
                for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                    const edm::RefToBase<reco::Track>& trkRef = *it;
                    if (trkRef.isNull()) continue;
                    double w = sv.trackWeight(trkRef);
                    //if (w < 0.5) continue;



                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                    trk_pt.push_back(trkRef->pt());
                    trk_weight.push_back(w);
                    trk_eta.push_back(trkRef->eta());
                    trk_phi.push_back(trkRef->phi());
                    int trackIdx = -1;
                    if (filteredTrackIndex.isValid() && filteredTrackIndex->contains(trkRef.id()))
                        trackIdx = (*filteredTrackIndex)[trkRef];
                    trk_trackIdx.push_back(trackIdx);
                    trk_sourceCollectionId.push_back(-1);
                    trk_sourceIdx.push_back(-1);
                    if (trackIdx >= 0) ++nMappedThisSV;
                    //edm::Ref<TrackCollection> mapRef = trackMap[trkRef];
                    //if (mapRef.isNonnull()) {
                    //    trk_SVscore.push_back((*svscoreHandle)[mapRef]);
                    //} else {
                    //    trk_SVscore.push_back(-1.0);
                    //}
                    trk_SVidx.push_back(x.size()-1);

                    // for GNN model:
                    //trk_ip_z.push_back(trkRef->ip_z())
                    //trk_ip_z_sig.push_back(trkRef->ip_z_sig())
                    //trk_ip2d.push_back(trkRef->ip2d())
                    //trk_ip3d.push_back(trkRef->ip3d())
                    //trk_ip2d_sig.push_back(trkRef->ip2d_sig())
                    //trk_ip3d_sig.push_back(trkRef->ip3d_sig())
                    trk_p.push_back(trkRef->p());
                    trk_charge.push_back(trkRef->charge());
                    trk_numberOfValidHits.push_back(trkRef->numberOfValidHits());
                    trk_numberOfValidPixelHits.push_back(trkRef->hitPattern().numberOfValidPixelHits());
                    trk_numberOfValidStripHits.push_back(trkRef->hitPattern().numberOfValidStripHits());


                    double ip_z = trkRef->dz(PV0.position());
                    double ip_z_sig = ip_z / trkRef->dzError();
                    trk_ip_z.push_back(ip_z);
                    trk_ip_z_sig.push_back(ip_z_sig);



                    GlobalVector direction(1,0,0);
                    direction = direction.unit();
                    
                    // Building Transient Track
                    const auto& ttBuilder = iSetup.getData(theTTBToken);
                    reco::TransientTrack ttrk = ttBuilder.build(*trkRef);
                    auto ip2d_val = IPTools::signedTransverseImpactParameter(ttrk, direction, PV0).second;
                    auto ip3d_val = IPTools::signedImpactParameter3D(ttrk, direction, PV0).second;

                    trk_ip2d.push_back(ip2d_val.value());
                    trk_ip3d.push_back(ip3d_val.value());
                    trk_ip2d_sig.push_back(ip2d_val.significance());
                    trk_ip3d_sig.push_back(ip3d_val.significance());

                    p4s_SV += p4;
                }
                }
                nMappedTracks.push_back(nMappedThisSV);
                pt.push_back(p4s_SV.Pt());
                eta.push_back(p4s_SV.Eta());
                phi.push_back(p4s_SV.Phi());
                mass.push_back(p4s_SV.M());


                //std::cout<<"SV: x=" << sv.x() 
                //    << " y=" << sv.y() 
                //    << " z=" << sv.z()
                //    << " nTracks=" << sv.tracksSize()
                //    << " pt=" << p4s_SV.Pt();
            }
        }


        table->addColumn<float>("x", x, "X position of SV");
        table->addColumn<float>("y", y, "Y position of SV");
        table->addColumn<float>("z", z, "Z position of SV");
        table->addColumn<float>("dlen", dlen, "dlen of SV");
        table->addColumn<float>("dlenSig", dlenSig, "dlenSig of SV");
        //std::cout<<"[DEBUG] SVTableProducer added x,y,z columns\n";
        table->addColumn<float>("chi2", chi2, "Chi2 of vertex fit");
        table->addColumn<float>("ndof", ndof, "Degrees of freedom of vertex fit");
        table->addColumn<int>("nTracks", nTracks, "Number of tracks in SV");
        table->addColumn<float>("pt", pt, "pt");
        table->addColumn<float>("eta", eta, "eta ");
        table->addColumn<float>("phi", phi, "phi");
        table->addColumn<float>("mass", mass, "mass");
        table->addColumn<int>("nMappedTracks", nMappedTracks,
                              "Number of SV constituents with a valid row in filteredTracks/trackTable");
        trk_table->addColumn<float>("trk_pt", trk_pt, "trk_pt");
        trk_table->addColumn<int>("trk_trackIdx", trk_trackIdx,
                                  "Event-local row index in filteredTracks/trackTable; -1 when absent");
        trk_table->addColumn<int>("trk_sourceCollectionId", trk_sourceCollectionId,
                                  "Original daughter collection: 0=packedPFCandidates, 1=lostTracks, -1=other/unknown");
        trk_table->addColumn<int>("trk_sourceIdx", trk_sourceIdx,
                                  "Original daughter key; -1 when the source pointer is null");
        trk_table->addColumn<float>("trk_weight", trk_weight, "trk_weight");
        trk_table->addColumn<float>("trk_eta", trk_eta, "trk_eta ");
        trk_table->addColumn<float>("trk_phi", trk_phi, "trk_phi");
        //trk_table->addColumn<float>("trk_SVscore", trk_SVscore, "trk_SVscore");


        trk_table->addColumn<float>("trk_p", trk_p, "trk_p");
        trk_table->addColumn<float>("trk_charge", trk_charge, "trk_charge ");
        trk_table->addColumn<float>("trk_numberOfValidHits", trk_numberOfValidHits, "trk_numberOfValidHits");
        trk_table->addColumn<float>("trk_numberOfValidPixelHits", trk_numberOfValidPixelHits, "trk_numberOfValidPixelHits");
        trk_table->addColumn<float>("trk_numberOfValidStripHits", trk_numberOfValidStripHits, "trk_numberOfValidStripHits ");
        trk_table->addColumn<float>("trk_ip_z", trk_ip_z, "trk_ip_z ");

        trk_table->addColumn<float>("trk_ip_z_sig", trk_ip_z_sig, "trk_ip_z_sig");
        trk_table->addColumn<float>("trk_ip2d", trk_ip2d, "trk_ip2d");
        trk_table->addColumn<float>("trk_ip3d", trk_ip3d, "trk_ip3d");
        trk_table->addColumn<float>("trk_ip2d_sig", trk_ip2d_sig, "trk_ip2d_sig");
        trk_table->addColumn<float>("trk_ip3d_sig", trk_ip3d_sig, "trk_ip3d_sig");


        trk_table->addColumn<int>("trk_SVidx", trk_SVidx, "trk_SVidx");

        iEvent.put(std::move(table), "SVTable");
        iEvent.put(std::move(trk_table), "SVtrksTable");
    } 


    
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVTableProducer);
