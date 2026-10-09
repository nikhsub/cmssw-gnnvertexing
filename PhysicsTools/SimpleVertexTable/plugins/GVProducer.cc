#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <optional>

// ---------------------------------------------------------------------------
// GVProducer
//
// Gen-level only. Walks genParticles, finds "good" B/D/S/Tau hadrons (>=2
// stable charged daughters passing kinematic cuts) and builds:
//   - GVTable            : one row per accepted hadron/GenVertex
//   - rejectedGVTable     : one row per hadron that failed the >=2 daughter cut
//   - GVDaughtersTable    : one row per stable charged daughter used to build a GV
//   - GVDirectDaughters   : one row per direct daughter of an accepted hadron
//
// It also republishes, as plain vector<T>/scalar EDM products, everything a
// downstream matching producer (GVMatchProducer) needs to redo the
// SV/track matching without re-walking the gen tree itself. Row order/
// indexing of these vectors matches GVTable/GVDaughtersTable exactly, so
// GVIdx-style indices stay valid across the two modules.
// ---------------------------------------------------------------------------

int packGenStatusFlags(const reco::GenStatusFlags& sf) {
    return  sf.isPrompt()                                <<  0 |
            sf.isDecayedLeptonHadron()                    <<  1 |
            sf.isTauDecayProduct()                        <<  2 |
            sf.isPromptTauDecayProduct()                  <<  3 |
            sf.isDirectTauDecayProduct()                  <<  4 |
            sf.isDirectPromptTauDecayProduct()            <<  5 |
            sf.isDirectHadronDecayProduct()                <<  6 |
            sf.isHardProcess()                             <<  7 |
            sf.fromHardProcess()                           <<  8 |
            sf.isHardProcessTauDecayProduct()              <<  9 |
            sf.isDirectHardProcessTauDecayProduct()        << 10 |
            sf.fromHardProcessBeforeFSR()                  << 11 |
            sf.isFirstCopy()                                << 12 |
            sf.isLastCopy()                                 << 13 |
            sf.isLastCopyBeforeFSR()                        << 14;
}

class GVProducer : public edm::stream::EDProducer<> {
public:
    explicit GVProducer(const edm::ParameterSet&);
    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    int checkPDG(int abs_pdg) const;
    std::optional<std::tuple<float, float, float>> isAncestor(const reco::Candidate* mother,
                                                                const reco::Candidate* daughter) const;
    bool hasHFAncestor(const reco::Candidate* hadron) const;
    bool hasHFDescendant(const reco::Candidate* hadron) const;

    edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
    double minHadronPt_;
    double minDaughterPt_;
    double maxHadronEta_;
    double maxDaughterEta_;
};

GVProducer::GVProducer(const edm::ParameterSet& iConfig)
    : genToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    minHadronPt_(iConfig.getParameter<double>("minHadronPt")),
      minDaughterPt_(iConfig.getParameter<double>("minDaughterPt")),
      maxHadronEta_(iConfig.getParameter<double>("maxHadronEta")),
      maxDaughterEta_(iConfig.getParameter<double>("maxDaughterEta")) {
    produces<nanoaod::FlatTable>("GVTable");
    produces<nanoaod::FlatTable>("rejectedGVTable");
    produces<nanoaod::FlatTable>("GVDaughtersTable");
    produces<nanoaod::FlatTable>("GVDirectDaughters");

    // Raw vectors handed off to GVMatchProducer. Row order matches the
    // corresponding FlatTable exactly (same fill loop, nothing re-sorted).
    produces<std::vector<float>>("hadronGVx");
    produces<std::vector<float>>("hadronGVy");
    produces<std::vector<float>>("hadronGVz");
    produces<std::vector<float>>("hadronEta");
    produces<std::vector<float>>("hadronPhi");

    produces<std::vector<float>>("daughterPt");
    produces<std::vector<float>>("daughterEta");
    produces<std::vector<float>>("daughterPhi");
    produces<std::vector<int>>("daughterCharge");
    produces<std::vector<int>>("daughterGVidx");

    produces<int>("nGV");
}

void GVProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
    edm::Handle<edm::View<reco::Candidate>> genHandle;
    iEvent.getByToken(genToken_, genHandle);
    const auto& genParticles = genHandle;

    // Output vectors
    std::vector<float> Hadron_pt, Hadron_eta, Hadron_phi;
    std::vector<float> Hadron_GVx, Hadron_GVy, Hadron_GVz;
    std::vector<float> Hadron_GVx_i, Hadron_GVy_i, Hadron_GVz_i;
    std::vector<int> Hadron_pdgId;
    std::vector<int> Hadron_pdgClass, Hadron_isB, Hadron_isD;
    std::vector<int> Hadron_fromHF, Hadron_toHF;
    std::vector<int> GV_nDaughters;
    std::vector<float> GV_maxDaughterPairDeltaR;
    std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
    std::vector<float> Daughters_vx, Daughters_vy, Daughters_vz;
    std::vector<int> Daughters_charge, Daughters_GVidx, Daughters_pdgId, Daughters_statusFlags;

    std::vector<float> allHadron_GVx, allHadron_GVy, allHadron_GVz;
    std::vector<float> allHadron_GVx_i, allHadron_GVy_i, allHadron_GVz_i;
    std::vector<int> allHadron_pdgId;
    std::vector<float> allHadron_pt, allHadron_eta, allHadron_phi;

    std::vector<float> directDaughters_pt, directDaughters_eta, directDaughters_phi;
    std::vector<int> directDaughters_charge, directDaughters_GVidx, directDaughters_pdgId;

    // Filling Hadrons and Daughters vectors
    int ngv = 0;
    int nRejectedGV = 0;
    int ngv_b = 0, ngv_d = 0, ngv_s = 0, ngv_tau = 0;
    for (size_t i = 0; i < genParticles->size(); ++i) {
        const reco::Candidate* hadron = &(*genParticles)[i];
        if (!(hadron->pt() > minHadronPt_ && std::abs(hadron->eta()) < maxHadronEta_)) continue;

        int hadPDG = checkPDG(std::abs(hadron->pdgId()));  // 1: Beauty, 2: Charmed, 3: Strange, 4: Tau, 0: Else
        if (hadPDG == 0) continue;

        // Collect stable charged daughters
        std::vector<float> temp_pt, temp_eta, temp_phi, temp_vx, temp_vy, temp_vz;
        std::vector<int> temp_charge, temp_GVidx, temp_flav, temp_pdgId, temp_statusFlags;
        int nPack = 0;
        float vx = std::numeric_limits<float>::quiet_NaN();
        float vy = std::numeric_limits<float>::quiet_NaN();
        float vz = std::numeric_limits<float>::quiet_NaN();

        for (size_t j = 0; j < genParticles->size(); ++j) {
            const reco::Candidate* dau = &(*genParticles)[j];
            if (dau == hadron) continue;
            if (!(dau->status() == 1 && dau->charge() != 0 && dau->pt() > minDaughterPt_ && std::abs(dau->eta()) < maxDaughterEta_)) continue;

            auto GV = isAncestor(hadron, dau);
            if (GV.has_value()) {
                std::tie(vx, vy, vz) = *GV;
                if (!std::isnan(vx)) {
                    nPack++;
                    temp_pt.push_back(dau->pt());
                    temp_eta.push_back(dau->eta());
                    temp_phi.push_back(dau->phi());
                    temp_vx.push_back(dau->vx());
                    temp_vy.push_back(dau->vy());
                    temp_vz.push_back(dau->vz());
                    temp_charge.push_back(dau->charge());
                    temp_pdgId.push_back(dau->pdgId());
                    const reco::GenParticle* genDau = dynamic_cast<const reco::GenParticle*>(dau);
                    if (genDau) {
                        temp_statusFlags.push_back(packGenStatusFlags(genDau->statusFlags()));
                    } else {
                        // handle non-GenParticle daughters, e.g. push a default/empty flags object
                    }
                    temp_GVidx.push_back(ngv);
                }
            }
        }

        // If has more than 2 good daughters, the Hadron is Good, we found a GV:
        if (nPack >= 2) {
            Hadron_pt.push_back(hadron->pt());
            Hadron_eta.push_back(hadron->eta());
            Hadron_phi.push_back(hadron->phi());
            Hadron_pdgId.push_back(hadron->pdgId());
            Hadron_pdgClass.push_back(hadPDG);

            ngv++;
            if (hadPDG == 1) {
                ngv_b++;
                Hadron_isB.push_back(1);
            } else {
                Hadron_isB.push_back(0);
            }
            if (hadPDG == 2) {
                ngv_d++;
                Hadron_isD.push_back(1);
            } else {
                Hadron_isD.push_back(0);
            }
            if (hadPDG == 3) ngv_s++;
            if (hadPDG == 4) ngv_tau++;
            Hadron_fromHF.push_back(hasHFAncestor(hadron) ? 1 : 0);
            Hadron_toHF.push_back(hasHFDescendant(hadron) ? 1 : 0);
            GV_nDaughters.push_back(nPack);

            // max pairwise deltaR among this GV's own daughters (gen-level, no tracks involved)
            float maxPairDR = 0.f;
            for (size_t a = 0; a < temp_eta.size(); ++a) {
                for (size_t b = a + 1; b < temp_eta.size(); ++b) {
                    float dR = deltaR(temp_eta[a], temp_phi[a], temp_eta[b], temp_phi[b]);
                    if (dR > maxPairDR) maxPairDR = dR;
                }
            }
            GV_maxDaughterPairDeltaR.push_back(maxPairDR);

            Hadron_GVx.push_back(vx);
            Hadron_GVy.push_back(vy);
            Hadron_GVz.push_back(vz);
            Hadron_GVx_i.push_back(hadron->vx());
            Hadron_GVy_i.push_back(hadron->vy());
            Hadron_GVz_i.push_back(hadron->vz());

            Daughters_pt.insert(Daughters_pt.end(), temp_pt.begin(), temp_pt.end());
            Daughters_vx.insert(Daughters_vx.end(), temp_vx.begin(), temp_vx.end());
            Daughters_vy.insert(Daughters_vy.end(), temp_vy.begin(), temp_vy.end());
            Daughters_vz.insert(Daughters_vz.end(), temp_vz.begin(), temp_vz.end());
            Daughters_eta.insert(Daughters_eta.end(), temp_eta.begin(), temp_eta.end());
            Daughters_phi.insert(Daughters_phi.end(), temp_phi.begin(), temp_phi.end());
            Daughters_charge.insert(Daughters_charge.end(), temp_charge.begin(), temp_charge.end());
            Daughters_pdgId.insert(Daughters_pdgId.end(), temp_pdgId.begin(), temp_pdgId.end());
            Daughters_statusFlags.insert(Daughters_statusFlags.end(), temp_statusFlags.begin(), temp_statusFlags.end());
            Daughters_GVidx.insert(Daughters_GVidx.end(), temp_GVidx.begin(), temp_GVidx.end());

            for (size_t j = 0; j < genParticles->size(); ++j) {
                const reco::Candidate* dau = &(*genParticles)[j];
                if (dau == hadron) continue;
                if (dau->numberOfMothers() > 0) {
                    const reco::Candidate* mother = dau->mother(0);
                    if (mother == hadron) {
                        directDaughters_pt.push_back(dau->pt());
                        directDaughters_eta.push_back(dau->eta());
                        directDaughters_phi.push_back(dau->phi());
                        directDaughters_charge.push_back(dau->charge());
                        directDaughters_pdgId.push_back(dau->pdgId());
                        directDaughters_GVidx.push_back(ngv - 1);
                    }
                }
            }
        } else {
            // fallback to save the decay point of the hadron if it has no daughters, but is still a rejected GV
            if (std::isnan(vx) && hadron->numberOfDaughters() > 0) {
                const reco::Candidate* directDau = hadron->daughter(0);
                vx = directDau->vx();
                vy = directDau->vy();
                vz = directDau->vz();
            }
            nRejectedGV++;
            allHadron_pt.push_back(hadron->pt());
            allHadron_eta.push_back(hadron->eta());
            allHadron_phi.push_back(hadron->phi());
            allHadron_pdgId.push_back(hadron->pdgId());
            allHadron_GVx.push_back(vx);
            allHadron_GVy.push_back(vy);
            allHadron_GVz.push_back(vz);
            allHadron_GVx_i.push_back(hadron->vx());
            allHadron_GVy_i.push_back(hadron->vy());
            allHadron_GVz_i.push_back(hadron->vz());
        }
    }

    // Build FlatTables
    auto rejectedGVTable = std::make_unique<nanoaod::FlatTable>(nRejectedGV, "RejectedGV", false);
    rejectedGVTable->addColumn<float>("pt", allHadron_pt, "Rejected Hadron pt");
    rejectedGVTable->addColumn<float>("eta", allHadron_eta, "Rejected Hadron eta");
    rejectedGVTable->addColumn<float>("phi", allHadron_phi, "Rejected Hadron phi");
    rejectedGVTable->addColumn<float>("x", allHadron_GVx, "Rejected GV x");
    rejectedGVTable->addColumn<float>("y", allHadron_GVy, "Rejected GV y");
    rejectedGVTable->addColumn<float>("z", allHadron_GVz, "Rejected GV z");
    rejectedGVTable->addColumn<float>("x_i", allHadron_GVx_i, "Born x coordinate of Rejected GV ");
    rejectedGVTable->addColumn<float>("y_i", allHadron_GVy_i, "Born y coordinate of Rejected GV ");
    rejectedGVTable->addColumn<float>("z_i", allHadron_GVz_i, "Born z coordinate of Rejected GV ");

    auto gvTable = std::make_unique<nanoaod::FlatTable>(ngv, "GV", false);
    gvTable->addColumn<float>("pt", Hadron_pt, "Hadron pt");
    gvTable->addColumn<float>("eta", Hadron_eta, "Hadron eta");
    gvTable->addColumn<float>("phi", Hadron_phi, "Hadron phi");
    gvTable->addColumn<float>("x", Hadron_GVx, "GV x");
    gvTable->addColumn<float>("y", Hadron_GVy, "GV y");
    gvTable->addColumn<float>("z", Hadron_GVz, "GV z");
    gvTable->addColumn<float>("x_i", Hadron_GVx_i, "Born x coordinate of GV ");
    gvTable->addColumn<float>("y_i", Hadron_GVy_i, "Born y coordinate of GV ");
    gvTable->addColumn<float>("z_i", Hadron_GVz_i, "Born z coordinate of GV ");
    gvTable->addColumn<int>("Hadron_pdgId", Hadron_pdgId, "Hadron_pdgId");
    gvTable->addColumn<int>("isB", Hadron_isB, "isB");
    gvTable->addColumn<int>("isD", Hadron_isD, "isD");
    gvTable->addColumn<int>("fromHF", Hadron_fromHF,"1 if an ancestor in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
    gvTable->addColumn<int>("toHF", Hadron_toHF,"1 if a descendant in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
    gvTable->addColumn<int>("pdgClass", Hadron_pdgId, "pdgClass");
    gvTable->addColumn<int>("nDaughters", GV_nDaughters,"Total number of GVDaughters belonging to this GV");
    gvTable->addColumn<float>("maxDaughterPairDeltaR", GV_maxDaughterPairDeltaR,"Max pairwise deltaR among this GV's own GVDaughters (gen-level, no tracks involved)");

    auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(), "GVDaughters", false);
    dauTable->addColumn<float>("pt", Daughters_pt, "Daughter pt");
    dauTable->addColumn<float>("eta", Daughters_eta, "Daughter eta");
    dauTable->addColumn<float>("phi", Daughters_phi, "Daughter phi");
    dauTable->addColumn<float>("vx", Daughters_vx, "Daughters_vx");
    dauTable->addColumn<float>("vy", Daughters_vy, "Daughters_vy");
    dauTable->addColumn<float>("vz", Daughters_vz, "Daughters_vz");
    dauTable->addColumn<int>("charge", Daughters_charge, "Daughter charge");
    dauTable->addColumn<int>("pdgId", Daughters_pdgId, "Daughter pdgId");
    dauTable->addColumn<int>("statusFlags", Daughters_statusFlags, "statusFlags");
    dauTable->addColumn<int>("hadronIndex", Daughters_GVidx, "Hadron index");

    auto directdauTable = std::make_unique<nanoaod::FlatTable>(directDaughters_pt.size(), "GVDirectDaughters", false);
    directdauTable->addColumn<float>("pt", directDaughters_pt, "Daughter pt");
    directdauTable->addColumn<float>("eta", directDaughters_eta, "Daughter eta");
    directdauTable->addColumn<float>("phi", directDaughters_phi, "Daughter phi");
    directdauTable->addColumn<int>("charge", directDaughters_charge, "Daughter charge");
    directdauTable->addColumn<int>("pdgId", directDaughters_pdgId, "Daughter pdgId");
    directdauTable->addColumn<int>("hadronIndex", directDaughters_GVidx, "Hadron index");

    iEvent.put(std::move(gvTable), "GVTable");
    iEvent.put(std::move(rejectedGVTable), "rejectedGVTable");
    iEvent.put(std::move(dauTable), "GVDaughtersTable");
    iEvent.put(std::move(directdauTable), "GVDirectDaughters");

    // Raw vectors for GVMatchProducer -- same order as gvTable/dauTable above.
    iEvent.put(std::make_unique<std::vector<float>>(Hadron_GVx), "hadronGVx");
    iEvent.put(std::make_unique<std::vector<float>>(Hadron_GVy), "hadronGVy");
    iEvent.put(std::make_unique<std::vector<float>>(Hadron_GVz), "hadronGVz");
    iEvent.put(std::make_unique<std::vector<float>>(Hadron_eta), "hadronEta");
    iEvent.put(std::make_unique<std::vector<float>>(Hadron_phi), "hadronPhi");

    iEvent.put(std::make_unique<std::vector<float>>(Daughters_pt), "daughterPt");
    iEvent.put(std::make_unique<std::vector<float>>(Daughters_eta), "daughterEta");
    iEvent.put(std::make_unique<std::vector<float>>(Daughters_phi), "daughterPhi");
    iEvent.put(std::make_unique<std::vector<int>>(Daughters_charge), "daughterCharge");
    iEvent.put(std::make_unique<std::vector<int>>(Daughters_GVidx), "daughterGVidx");

    iEvent.put(std::make_unique<int>(ngv), "nGV");
}

//  checkPDG()
int GVProducer::checkPDG(int abs_pdg) const {
    std::vector<int> pdgList_B = {521, 511, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242, 5342, 5512, 5532, 5542, 5554};
    std::vector<int> pdgList_D = {411, 421, 431, 4122, 4232, 4132, 4332, 4412, 4422, 4432, 4444};
    std::vector<int> pdgList_S = {3122, 3222, 3212, 3312, 3322, 3334};
    std::vector<int> pdgList_Tau = {15};

    if (std::find(pdgList_B.begin(), pdgList_B.end(), abs_pdg) != pdgList_B.end()) return 1;
    if (std::find(pdgList_D.begin(), pdgList_D.end(), abs_pdg) != pdgList_D.end()) return 2;
    if (std::find(pdgList_S.begin(), pdgList_S.end(), abs_pdg) != pdgList_S.end()) return 3;
    if (std::find(pdgList_Tau.begin(), pdgList_Tau.end(), abs_pdg) != pdgList_Tau.end()) return 4;
    return 0;
}

//  hasHFAncestor()
bool GVProducer::hasHFAncestor(const reco::Candidate* hadron) const {
    const reco::Candidate* current = hadron;
    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (checkPDG(std::abs(mother->pdgId())) != 0) {
            return true;
        }
        current = mother;
    }
    return false;
}

//  hasHFDescendant()
bool GVProducer::hasHFDescendant(const reco::Candidate* hadron) const {
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
std::optional<std::tuple<float, float, float>> GVProducer::isAncestor(const reco::Candidate* ancestor,
                                                                              const reco::Candidate* particle) const {
    std::vector<int> pdgList_B = {521, 511, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242, 5342, 5512, 5532, 5542, 5554};
    std::vector<int> pdgList_D = {411, 421, 431, 4122, 4232, 4132, 4332, 4412, 4422, 4432, 4444};
    std::vector<int> pdgList_S = {3122, 3222, 3212, 3312, 3322, 3334};
    std::vector<int> pdgList_Tau = {15};
    std::unordered_set<int> pdgSet_D(pdgList_D.begin(), pdgList_D.end());
    std::unordered_set<int> pdgSet_B(pdgList_B.begin(), pdgList_B.end());
    std::unordered_set<int> pdgSet_S(pdgList_S.begin(), pdgList_S.end());
    std::unordered_set<int> pdgSet_Tau(pdgList_Tau.begin(), pdgList_Tau.end());
    const reco::Candidate* current = particle;

    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == ancestor) {
            // Found the ancestor; return the vertex of the current particle (i.e., the direct daughter)
            return std::make_optional(std::make_tuple(current->vx(), current->vy(), current->vz()));
        }
        int mother_pdg = std::abs(mother->pdgId());
        if (pdgSet_B.count(mother_pdg) || pdgSet_D.count(mother_pdg) || pdgSet_S.count(mother_pdg) ||
            pdgSet_Tau.count(mother_pdg))
            break;
        current = mother;
    }

    // If we reached here, the ancestor was not found in the chain
    return std::nullopt;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GVProducer);
