#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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
    double dlenSigMin_;  
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
    //edm::EDGetTokenT<edm::ValueMap<float>> svscoreToken_;
};

SVTableProducer::SVTableProducer(const edm::ParameterSet &iConfig): 
    svToken(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
    //svscoreToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pvSrc")))
{
    produces<nanoaod::FlatTable>("SVTable");
    produces<nanoaod::FlatTable>("SVtrksTable");
}

void SVTableProducer::produce(edm::StreamID,
                 edm::Event &iEvent,
                 const edm::EventSetup &iSetup) const 
    {
        edm::Handle<std::vector<reco::Vertex>> svs;
        iEvent.getByToken(svToken, svs);
        auto pvsIn = iEvent.getHandle(pvs_);
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
        unsigned int nSVtracks = 0;
        int nTrksCurrentSV = 0;
        unsigned int nSV_cutdlen = 0;
        for (auto const& sv : *svs) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                nTrksCurrentSV = 0;
                for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                    const edm::RefToBase<reco::Track>& trkRef = *it;
                    if (trkRef.isNull()) continue;
                    double w = sv.trackWeight(trkRef);   // <-- weight comes from vertex
                    //if (w < 0.5) continue;               // skip low-weight tracks
                    nTrksCurrentSV++;
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

        std::vector<int> nTracks;
        std::vector<int> trk_SVidx;

        int nTrksPerSV = 0;
        //std::cout<<svs->size()<<" SVs to process\n";
        for (const auto &sv : *svs) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                nTrksPerSV = 0;
                // First count tracks with weight >= 0.5
                for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                    const edm::RefToBase<reco::Track>& trkRef = *it;
                    if (trkRef.isNull()) continue;
                    double w = sv.trackWeight(trkRef);   // <-- weight comes from vertex
                    //if (w < 0.5) continue;               // skip low-weight tracks
                    nTrksPerSV++;
                }
                if (nTrksPerSV < 2) continue; // skip SVs with less than 2 high-weight tracks

                x.push_back(sv.x());
                y.push_back(sv.y());
                z.push_back(sv.z());
                chi2.push_back(sv.chi2());
                ndof.push_back(sv.ndof());
                nTracks.push_back(nTrksPerSV);
                dlen.push_back(dl.value());
                dlenSig.push_back(dl.significance());
                //std::cout<<"dl significance is "<< dl.significance()<<std::endl;


                TLorentzVector p4s_SV;
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
        trk_table->addColumn<float>("trk_pt", trk_pt, "trk_pt");
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
