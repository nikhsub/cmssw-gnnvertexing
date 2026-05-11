#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TLorentzVector.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>
    std::vector<int> pdgList_B = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    std::vector<int> pdgList_D = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    std::vector<int> pdgList_S = {3122,3222,3212,3312,3322,3334};
    std::vector<int> pdgList_Tau = {15};
    std::unordered_set<int> pdgSet_D(pdgList_D.begin(), pdgList_D.end());
    std::unordered_set<int> pdgSet_B(pdgList_B.begin(), pdgList_B.end());
    std::unordered_set<int> pdgSet_S(pdgList_S.begin(), pdgList_S.end());
    std::unordered_set<int> pdgSet_Tau(pdgList_Tau.begin(), pdgList_Tau.end());
class GenVertexCandidateProducer : public edm::stream::EDProducer<> {
public:
    explicit GenVertexCandidateProducer(const edm::ParameterSet&);

    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    int checkPDG(int abs_pdg) const;

    std::optional<std::tuple<float, float, float>>isAncestor(const reco::Candidate* mother,const reco::Candidate* daughter) const;

    std::vector<std::vector<float>> computeDistanceMatrix(
                    const std::vector<float>& SV_x,const std::vector<float>& SV_y,const std::vector<float>& SV_z,
                    const std::vector<float>& Hadron_GVx,const std::vector<float>& Hadron_GVy,const std::vector<float>& Hadron_GVz);
    void printDistanceMatrix(const std::vector<std::vector<float>>& distances);
    std::pair<std::vector<int>, std::vector<float>>  matchHadronsToSV(
                std::vector<std::vector<float>> distances,
                const std::vector<float>& SVtrk_pt, const std::vector<float>& SVtrk_eta, const std::vector<float>& SVtrk_phi, const std::vector<int>& SVtrk_SVidx,
                const std::vector<float>& Daughters_pt,     //genparticles
                const std::vector<float>& Daughters_eta,  //genparticles
                const std::vector<float>& Daughters_phi,  //genparticles
                const std::vector<int>& Daughters_GVidx, // hadron index per daughter
                int n_Hadrons,
                int nRequiredCommonTracks,
                double dR_max,
                double relPt_max
            );

    const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> svToken_;
    int nRequiredCommonTracks_;
    double dlenSigMin_;
    double dR_max_;
    double relPt_max_;
};


GenVertexCandidateProducer::GenVertexCandidateProducer(const edm::ParameterSet& iConfig):
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    genToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    svToken_(consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
    nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    dR_max_(iConfig.getParameter<double>("dR_max")),
    relPt_max_(iConfig.getParameter<double>("relPt_max"))
{
    produces<nanoaod::FlatTable>("GVTable");
    produces<nanoaod::FlatTable>("GVDaughtersTable");
    produces<nanoaod::FlatTable>("SVDaughtersTable");
}


void GenVertexCandidateProducer::produce(edm::Event& iEvent,
             const edm::EventSetup&) 
    {
        edm::Handle<edm::View<reco::Candidate>> genHandle;
        iEvent.getByToken(genToken_, genHandle);
        edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> svHandle;
        iEvent.getByToken(svToken_, svHandle);
        auto pvsIn = iEvent.getHandle(pvs_);

        



        const auto& genParticles = genHandle;
        const auto& secondaryVertices = svHandle;

        // Output vectors
        std::vector<float> Hadron_pt, Hadron_eta, Hadron_phi;
        std::vector<float> SV_x, SV_y, SV_z;
        std::vector<float> Hadron_GVx, Hadron_GVy, Hadron_GVz;
        std::vector<float>  Hadron_GVx_i, Hadron_GVy_i, Hadron_GVz_i;
        std::vector<int> Hadron_pdgClass, Hadron_isB, Hadron_isD;
        std::vector<int> Hadron_pdgId;
        std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
        std::vector<int> Daughters_charge, Daughters_GVidx;
        VertexDistance3D vdist;

        //const reco::VertexCompositePtrCandidate* vcand;
        


        const auto& PV0 = pvsIn->front();
        // save coordinates of SV (will be used for matching with GV)
        for (auto const& sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
            SV_x.push_back(sv.position().x());
            SV_y.push_back(sv.position().y());
            SV_z.push_back(sv.position().z());
            
            
            //const auto& vtxPos = vcand->vertex();
            //const auto& cov3 = vcand->vertexCovariance();
            }
        }

        // Filling Hadrons and Daughters vectors
        int ngv=0;
        int ngv_b=0, ngv_d=0, ngv_s=0, ngv_tau=0;
        for(size_t i=0; i<genParticles->size(); ++i){
            const reco::Candidate* hadron = &(*genParticles)[i];
            //std::cout<<"Hadron "<<i<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
            if(!(hadron->pt()>10 && std::abs(hadron->eta())<2.5)) continue;

            int hadPDG = checkPDG(std::abs(hadron->pdgId())); // 1: Beauty, 2: Charmed, 3: Strange,  4: Tau,  0: Else
            if(hadPDG==0) continue;

            //  Collect stable charged daughters
            std::vector<float> temp_pt, temp_eta, temp_phi; // kinematics of gen daughters of the hadron in the loop
            std::vector<int> temp_charge, temp_GVidx, temp_flav;
            int nPack=0;
            float vx=std::numeric_limits<float>::quiet_NaN();
            float vy=std::numeric_limits<float>::quiet_NaN();
            float vz=std::numeric_limits<float>::quiet_NaN();

            for(size_t j=0; j<genParticles->size(); ++j){
                const reco::Candidate* dau = &(*genParticles)[j];
                if(dau==hadron) continue;
                if(!(dau->status()==1 && dau->charge()!=0 && dau->pt()>0.8 && std::abs(dau->eta())<2.5)) continue;

                auto GV = isAncestor(hadron,dau); //takes the x,y,z of the daughter (decay point of the hadron)
                if(GV.has_value()){
                    std::tie(vx,vy,vz) = *GV;
                    if(!std::isnan(vx)){
                        nPack++;
                        temp_pt.push_back(dau->pt());
                        temp_eta.push_back(dau->eta());
                        temp_phi.push_back(dau->phi());
                        temp_charge.push_back(dau->charge());
                        temp_GVidx.push_back(ngv); // hadron index
                        //temp_flav.push_back(hadPDG);
                    }
                }
            }
            // If has more than 2 good daughters, the Hadron is Good, we found a GV:
            if(nPack>=2){
                // Save hadron info
                //std::cout<<"Found hadron "<<ngv<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
                Hadron_pt.push_back(hadron->pt());
                Hadron_eta.push_back(hadron->eta());
                Hadron_phi.push_back(hadron->phi());
                Hadron_pdgId.push_back(hadron->pdgId());

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
                if(hadPDG==3) {
                    ngv_s++;
                    }
                if(hadPDG==4) {
                    ngv_tau++;
                    }
                Hadron_GVx.push_back(vx);               // point of decay of the hadron
                Hadron_GVy.push_back(vy);               // point of decay of the hadron
                Hadron_GVz.push_back(vz);               // point of decay of the hadron
                Hadron_GVx_i.push_back(hadron->vx());   // point of origin of the hadron
                Hadron_GVy_i.push_back(hadron->vy());   // point of origin of the hadron
                Hadron_GVz_i.push_back(hadron->vz());   // point of origin of the hadron

                // Save daughters
                Daughters_pt.insert(Daughters_pt.end(), temp_pt.begin(), temp_pt.end());
                Daughters_eta.insert(Daughters_eta.end(), temp_eta.begin(), temp_eta.end());
                Daughters_phi.insert(Daughters_phi.end(), temp_phi.begin(), temp_phi.end());
                Daughters_charge.insert(Daughters_charge.end(), temp_charge.begin(), temp_charge.end());
                Daughters_GVidx.insert(Daughters_GVidx.end(), temp_GVidx.begin(), temp_GVidx.end());
            }
        }


        // Filling tracks from reco SV
        std::vector<float> SVtrk_pt, SVtrk_eta, SVtrk_phi;
        std::vector<int> SVtrk_SVidx;
        
        
        int SV_index=0;
        for (const auto &sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
        
        


            for (size_t i = 0; i < sv.numberOfDaughters(); ++i) {
                 const reco::Candidate* dau = sv.daughter(i);
                 if(dau->charge()==0) continue; // only charged tracks


                //TLorentzVector p4;
                //p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                SVtrk_pt.push_back(dau->pt());
                SVtrk_eta.push_back(dau->eta());
                SVtrk_phi.push_back(dau->phi());
                SVtrk_SVidx.push_back(SV_index); 
            }
            SV_index++;
            }
        }
        
        // Compute matrix of distances between SV and GV
        auto distances = computeDistanceMatrix(SV_x, SV_y, SV_z, Hadron_GVx, Hadron_GVy, Hadron_GVz);
        //printDistanceMatrix(distances);
        
        std::vector<int> Hadron_SVIdx(ngv, -1); // 
        std::vector<float> Hadron_SVDistance(ngv, -1); // 

        // perform matching based on distance matrix and track-to-daughter matching
        auto result = matchHadronsToSV(
            distances,
            SVtrk_pt, SVtrk_eta, SVtrk_phi, SVtrk_SVidx,
            Daughters_pt, Daughters_eta, Daughters_phi, Daughters_GVidx,
            ngv,
            nRequiredCommonTracks_,
            dR_max_,
            relPt_max_
            );

        Hadron_SVIdx = result.first;
        Hadron_SVDistance = result.second;

        //  Build FlatTables 

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
        // new class
        gvTable->addColumn<int>("isB",Hadron_isB,"isB");
        gvTable->addColumn<int>("isD",Hadron_isD,"isD");
        gvTable->addColumn<int>("pdgClass",Hadron_pdgId,"pdgClass");
        
        
        //

        auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(),"GVDaughters",false);
        dauTable->addColumn<float>("pt",Daughters_pt,"Daughter pt");
        dauTable->addColumn<float>("eta",Daughters_eta,"Daughter eta");
        dauTable->addColumn<float>("phi",Daughters_phi,"Daughter phi");
        dauTable->addColumn<int>("charge",Daughters_charge,"Daughter charge");
        dauTable->addColumn<int>("hadronIndex",Daughters_GVidx,"Hadron index");


        auto svdauTable = std::make_unique<nanoaod::FlatTable>(SVtrk_pt.size(),"SVDaughters",false);
        svdauTable->addColumn<float>("pt",SVtrk_pt,"Daughter pt");
        svdauTable->addColumn<float>("eta",SVtrk_eta,"Daughter eta");
        svdauTable->addColumn<float>("phi",SVtrk_phi,"Daughter phi");
        //svdauTable->addColumn<int>("charge",SVtrk_charge,"Daughter charge");
        svdauTable->addColumn<int>("SVIdx",SVtrk_SVidx,"Hadron index");


        //
        iEvent.put(std::move(gvTable),"GVTable");
        iEvent.put(std::move(dauTable),"GVDaughtersTable");
        iEvent.put(std::move(svdauTable),"SVDaughtersTable");
    }






//  checkPDG() 
int GenVertexCandidateProducer::checkPDG(int abs_pdg) const {

    if(std::find(pdgList_B.begin(),pdgList_B.end(),abs_pdg)!=pdgList_B.end()) return 1;
    if(std::find(pdgList_D.begin(),pdgList_D.end(),abs_pdg)!=pdgList_D.end()) return 2;
    if(std::find(pdgList_S.begin(),pdgList_S.end(),abs_pdg)!=pdgList_S.end()) return 3;
    if(std::find(pdgList_Tau.begin(),pdgList_Tau.end(),abs_pdg)!=pdgList_Tau.end()) return 4;
    return 0;
}



//  isAncestor() 
std::optional<std::tuple<float, float, float>> GenVertexCandidateProducer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle) const
// ancestor: Candidate Object of the Mother
// particle: Candidate object of the particle that can be daughter
    {
    const reco::Candidate* current = particle;

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


std::vector<std::vector<float>> GenVertexCandidateProducer::computeDistanceMatrix(
                const std::vector<float>& SV_x,
                const std::vector<float>& SV_y,
                const std::vector<float>& SV_z,
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
        for (size_t j = 0; j < nHadron; ++j) {
            float dx = SV_x[i] - Hadron_GVx[j];
            float dy = SV_y[i] - Hadron_GVy[j];
            float dz = SV_z[i] - Hadron_GVz[j];
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            distances[i][j] = dist;
        }
    }

    return distances;
}


void GenVertexCandidateProducer::printDistanceMatrix(
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


std::pair<std::vector<int>, std::vector<float>>  GenVertexCandidateProducer::matchHadronsToSV(
    std::vector<std::vector<float>> distances,
    const std::vector<float>& SVtrk_pt,
    const std::vector<float>& SVtrk_eta,
    const std::vector<float>& SVtrk_phi,
    const std::vector<int>& SVtrk_SVidx,
    const std::vector<float>& Daughters_pt,     //genparticles
    const std::vector<float>& Daughters_eta,  //genparticles
    const std::vector<float>& Daughters_phi,  //genparticles
    const std::vector<int>& Daughters_GVidx, // hadron index per daughter
    int n_Hadrons,
    int nRequiredCommonTracks,
    double dR_max,
    double relPt_max
) {
    // Returns
    // Hadron_SVIdx = Array of lenght = Hadron_pt.size() with index of the SV
    size_t nSV = distances.size(); //distances is nSV x nHadron matrix (the first dimension is nSV)
    std::vector<int> Hadron_SVIdx(n_Hadrons, -1); // Output
    std::vector<float> Hadron_SVDistance(n_Hadrons, -1); // Output
    std::vector<size_t> svTrackIdxs_fromBestSV;

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
            if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > 0.8 && std::fabs(SVtrk_eta[i]) < 2.5) {
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
        for (size_t iSV : svTrackIdxs_fromBestSV) {
            //std::cout<<" Checking SVtrack index "<<iSV<<std::endl;
            for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                //std::cout<<" Checking GVDaughters index "<<iHad<<std::endl;
                float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / Daughters_pt[iHad];
                //std::cout<<"Comparing SV track (pt: "<<SVtrk_pt[iSV]<<", eta: "<<SVtrk_eta[iSV]<<", phi: "<<SVtrk_phi[iSV]<<") with Daughter (pt: "<<Daughters_pt[iHad]<<", eta: "<<Daughters_eta[iHad]<<", phi: "<<Daughters_phi[iHad]<<") => dR: "<<dR<<", relPt: "<<relPt<<std::endl;
                if (dR < dR_max && relPt < relPt_max) {
                    ++common;
                    //std::cout<<"  -> Matched! Common tracks: "<<common<<std::endl;
                    if (common >= nRequiredCommonTracks) break; // break the iHad cycle
                }
            }
            if (common >= nRequiredCommonTracks) break; // break the iSV cycle
        }

        if (common >= nRequiredCommonTracks) {
            Hadron_SVIdx[bestHad] = bestSV;
            Hadron_SVDistance[bestHad]= minDist;
            for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 998.0;  // remove SV row
            for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 998.0;   // remove Hadron column
            //std::cout << "[V] Matched Hadron[" << bestHad << "] to SV[" << bestSV << "] (distance = " << minDist << ", common tracks = " << common << ")\n";
        } else {
            distances[bestSV][bestHad] = 998.0;  // exclude this pair
        }
    }

    return std::make_pair(Hadron_SVIdx, Hadron_SVDistance);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenVertexCandidateProducer);

