#include <string>
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
    std::string tableName_;
    std::string trackTableName_;
    std::string pvTableName_;
    bool applySelection_;
    bool savePV_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
    //edm::EDGetTokenT<edm::ValueMap<float>> svscoreToken_;
};

SVTableProducer::SVTableProducer(const edm::ParameterSet &iConfig):
    svToken(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    tableName_(iConfig.getParameter<std::string>("tableName")),
    trackTableName_(iConfig.getParameter<std::string>("trackTableName")),
    pvTableName_(iConfig.getParameter<std::string>("pvTableName")),
    applySelection_(iConfig.getParameter<bool>("applySelection")),
    savePV_(iConfig.getParameter<bool>("savePV")),
    theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
{
    produces<nanoaod::FlatTable>("SVTable");
    produces<nanoaod::FlatTable>("SVtrksTable");
    if (savePV_) {
        produces<nanoaod::FlatTable>("PVTable");
    }
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
	        // --- Optional PV table (leading vertex only) ---
        if (savePV_) {
            auto pv_table = std::make_unique<nanoaod::FlatTable>(1, pvTableName_, true);

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
        }

        unsigned int nSVtracks = 0;
        unsigned int nSVSelected = 0;

        for (const auto& sv : *svs) {
            Measurement1D dl = vdist.distance(
                PV0,
                VertexState(
                    RecoVertex::convertPos(sv.position()),
                    RecoVertex::convertError(sv.error())
                )
            );

            unsigned int nNonNullTracks = 0;
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                if (trkRef.isNull()) continue;
                ++nNonNullTracks;
            }

            const bool acceptVertex =
                !applySelection_ ||
                (
                    dl.value() > 0.0 &&
                    dl.significance() > dlenSigMin_ &&
                    nNonNullTracks >= 2
                );

            if (!acceptVertex) continue;

            nSVtracks += nNonNullTracks;
            ++nSVSelected;
        }

        auto table = std::make_unique<nanoaod::FlatTable>(
            nSVSelected, tableName_, false
        );
        auto trk_table = std::make_unique<nanoaod::FlatTable>(
            nSVtracks, trackTableName_, false
        );


        std::vector<float> x, y, z, chi2, ndof, pt, eta, phi, mass, dlen, dlenSig;
	std::vector<float> covXX, covXY, covXZ, covYY, covYZ, covZZ;
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
        std::vector<int> svIdx;
        std::vector<int> trk_SVidx;
        std::vector<int> trk_globalIdx;
        std::vector<int> trk_localIdx_inSV;

        int nTrksPerSV = 0;
        for (size_t inputSVIdx = 0; inputSVIdx < svs->size(); ++inputSVIdx) {
            const auto& sv = svs->at(inputSVIdx);

            Measurement1D dl = vdist.distance(
                PV0,
                VertexState(
                    RecoVertex::convertPos(sv.position()),
                    RecoVertex::convertError(sv.error())
                )
            );

            nTrksPerSV = 0;
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                if (trkRef.isNull()) continue;
                ++nTrksPerSV;
            }

            const bool acceptVertex =
                !applySelection_ ||
                (
                    dl.value() > 0.0 &&
                    dl.significance() > dlenSigMin_ &&
                    nTrksPerSV >= 2
                );

            if (!acceptVertex) continue;

            svIdx.push_back(static_cast<int>(inputSVIdx));

            x.push_back(sv.x());
                y.push_back(sv.y());
                z.push_back(sv.z());

		const auto& cov = sv.covariance();
		covXX.push_back(cov(0,0));
		covXY.push_back(cov(0,1));
		covXZ.push_back(cov(0,2));
		covYY.push_back(cov(1,1));
		covYZ.push_back(cov(1,2));
		covZZ.push_back(cov(2,2));

                chi2.push_back(sv.chi2());
                ndof.push_back(sv.ndof());
                nTracks.push_back(nTrksPerSV);
                dlen.push_back(dl.value());
                dlenSig.push_back(dl.significance());
                //std::cout<<"dl significance is "<< dl.significance()<<std::endl;


                TLorentzVector p4s_SV;
                int localIdxInSV = 0;
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
                    trk_localIdx_inSV.push_back(localIdxInSV);
                    trk_globalIdx.push_back(static_cast<int>(trkRef.key()));

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
                    ++localIdxInSV;
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


        table->addColumn<int>("svIdx", svIdx, "Index in input reco::Vertex collection");
        table->addColumn<float>("x", x, "X position of SV");
        table->addColumn<float>("y", y, "Y position of SV");
        table->addColumn<float>("z", z, "Z position of SV");

	table->addColumn<float>("covXX", covXX, "SV covariance matrix element xx");
	table->addColumn<float>("covXY", covXY, "SV covariance matrix element xy");
	table->addColumn<float>("covXZ", covXZ, "SV covariance matrix element xz");
	table->addColumn<float>("covYY", covYY, "SV covariance matrix element yy");
	table->addColumn<float>("covYZ", covYZ, "SV covariance matrix element yz");
	table->addColumn<float>("covZZ", covZZ, "SV covariance matrix element zz");

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
        trk_table->addColumn<int>("trk_localIdx_inSV", trk_localIdx_inSV, "Track index local to each SV (convenience only)");
        trk_table->addColumn<int>("trk_globalIdx", trk_globalIdx, "Authoritative global index in upstream unpacked track collection");

        iEvent.put(std::move(table), "SVTable");
        iEvent.put(std::move(trk_table), "SVtrksTable");
    } 


    
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVTableProducer);
