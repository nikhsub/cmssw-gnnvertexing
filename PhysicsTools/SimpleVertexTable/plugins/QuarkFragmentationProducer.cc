// Generator-only MiniAOD truth producer for quark-fragmentation studies.
//
// Stage 1 deliberately makes no angular or momentum based assignments.  An
// object receives a quark-system index only when an unambiguous retained
// mother-chain path reaches one selected pruned quark system.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
  constexpr float kInvalid = -999.f;
  enum AssociationType { kNoAssociation = 0, kRetainedAncestry = 1 };
  enum AssociationQuality { kMissing = 0, kUnique = 1, kAmbiguous = 2 };
  enum MembershipType { kUnassigned = 0, kHeavyHadronDecay = 1, kRetainedFragmentation = 2, kAmbiguousMembership = 3 };
  enum CharmOrigin { kCharmNotApplicable = 0, kPromptCharm = 1, kCharmFromB = 2, kCharmUnknown = 3 };

  bool isQuark(int id) { return std::abs(id) >= 1 && std::abs(id) <= 5; }
  int flavour(int id) { return std::abs(id); }
  bool isBHadron(int id) {
    const int a = std::abs(id);
    static const std::unordered_set<int> ids = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    return ids.count(a) != 0;
  }
  bool isDHadron(int id) {
    const int a = std::abs(id);
    static const std::unordered_set<int> ids = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    return ids.count(a) != 0;
  }
  bool isHeavyHadron(int id) { return isBHadron(id) || isDHadron(id); }
  int hadronFlavour(int id) { return isBHadron(id) ? 5 : (isDHadron(id) ? 4 : 0); }

  struct Ancestors {
    std::vector<int> indices;
    bool incomplete = false;
  };
}

class QuarkFragmentationProducer : public edm::stream::EDProducer<> {
public:
  explicit QuarkFragmentationProducer(const edm::ParameterSet& cfg)
      : prunedToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
        packedToken_(consumes<pat::PackedGenParticleCollection>(cfg.getParameter<edm::InputTag>("packedGenParticles"))) {
    produces<nanoaod::FlatTable>("QFragTable");
    produces<nanoaod::FlatTable>("QHadronTable");
    produces<nanoaod::FlatTable>("QGenChargedTable");
  }
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  Ancestors ancestors(int start, const reco::GenParticleCollection& particles,
                      const std::unordered_map<const reco::Candidate*, int>& index) const;
  bool isFinalCopy(int index, const reco::GenParticleCollection& particles,
                   const std::unordered_map<const reco::Candidate*, int>& lookup) const;
  std::vector<int> samePdgComponent(int start, const reco::GenParticleCollection& particles,
                                    const std::unordered_map<const reco::Candidate*, int>& lookup) const;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedToken_;
};

Ancestors QuarkFragmentationProducer::ancestors(
    int start, const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup) const {
  Ancestors result;
  if (start < 0 || static_cast<size_t>(start) >= particles.size()) { result.incomplete = true; return result; }
  std::queue<int> todo;
  std::unordered_set<int> seen;
  todo.push(start);
  seen.insert(start);
  while (!todo.empty()) {
    const int current = todo.front(); todo.pop();
    const auto& particle = particles[current];
    for (size_t im = 0; im < particle.numberOfMothers(); ++im) {
      const reco::Candidate* mother = particle.mother(im);
      if (!mother) { result.incomplete = true; continue; }
      const auto found = lookup.find(mother);
      if (found == lookup.end()) { result.incomplete = true; continue; }
      if (seen.insert(found->second).second) {
        result.indices.push_back(found->second);
        todo.push(found->second);
      }
    }
  }
  return result;
}

bool QuarkFragmentationProducer::isFinalCopy(
    int index, const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup) const {
  const auto& particle = particles[index];
  for (size_t id = 0; id < particle.numberOfDaughters(); ++id) {
    const reco::Candidate* daughter = particle.daughter(id);
    if (!daughter || daughter == &particle) continue;
    if (daughter->pdgId() == particle.pdgId() && lookup.count(daughter)) return false;
  }
  return true;
}

std::vector<int> QuarkFragmentationProducer::samePdgComponent(
    int start, const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup) const {
  std::vector<int> out;
  std::queue<int> todo;
  std::unordered_set<int> seen;
  const int id = particles[start].pdgId();
  todo.push(start); seen.insert(start);
  while (!todo.empty()) {
    const int current = todo.front(); todo.pop(); out.push_back(current);
    const auto& p = particles[current];
    auto add = [&](const reco::Candidate* other) {
      const auto found = lookup.find(other);
      if (other && found != lookup.end() && other->pdgId() == id && seen.insert(found->second).second) todo.push(found->second);
    };
    for (size_t im = 0; im < p.numberOfMothers(); ++im) add(p.mother(im));
    for (size_t idau = 0; idau < p.numberOfDaughters(); ++idau) add(p.daughter(idau));
  }
  return out;
}

void QuarkFragmentationProducer::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<reco::GenParticleCollection> pruned;
  edm::Handle<pat::PackedGenParticleCollection> packed;
  event.getByToken(prunedToken_, pruned);
  event.getByToken(packedToken_, packed);

  std::unordered_map<const reco::Candidate*, int> prunedIndex;
  for (size_t i = 0; i < pruned->size(); ++i) prunedIndex.emplace(&(*pruned)[i], static_cast<int>(i));

  // A system is the terminal retained member of a connected same-PDG copy chain.
  // The representative is the terminal member; hardQuarkIdx stores a hard-process
  // member of that chain when exactly one is retained.
  std::vector<int> qRep, qHard, qLast;
  std::unordered_map<int, int> prunedToQFrag;
  std::unordered_set<int> handledQuarks;
  for (size_t i = 0; i < pruned->size(); ++i) {
    if (!isQuark((*pruned)[i].pdgId()) || handledQuarks.count(i)) continue;
    const auto component = samePdgComponent(i, *pruned, prunedIndex);
    for (const int member : component) handledQuarks.insert(member);
    std::vector<int> terminal, hard;
    for (const int member : component) {
      if (isFinalCopy(member, *pruned, prunedIndex)) terminal.push_back(member);
      if ((*pruned)[member].statusFlags().isHardProcess()) hard.push_back(member);
    }
    if (terminal.size() != 1) continue;  // preserves ambiguous split histories without inventing a system
    const int row = static_cast<int>(qRep.size());
    const int representative = terminal.front();
    qRep.push_back(representative);
    qLast.push_back(representative);
    qHard.push_back(hard.size() == 1 ? hard.front() : -1);
    for (const int member : component) prunedToQFrag[member] = row;
  }

  std::vector<int> hPruned;
  for (size_t i = 0; i < pruned->size(); ++i)
    if (isHeavyHadron((*pruned)[i].pdgId()) && isFinalCopy(i, *pruned, prunedIndex)) hPruned.push_back(i);
  std::unordered_map<int, int> prunedToHadron;
  for (size_t i = 0; i < hPruned.size(); ++i) prunedToHadron[hPruned[i]] = static_cast<int>(i);

  std::vector<int> hadQuark(hPruned.size(), -1), hadAssocQ(hPruned.size(), kMissing), hadFromB(hPruned.size(), 0), hadParentB(hPruned.size(), -1);
  std::vector<int> hadCharmOrigin(hPruned.size(), kCharmNotApplicable);
  for (size_t ih = 0; ih < hPruned.size(); ++ih) {
    const auto ancestry = ancestors(hPruned[ih], *pruned, prunedIndex);
    std::unordered_set<int> systems, bAncestors;
    for (const int ancestor : ancestry.indices) {
      const auto q = prunedToQFrag.find(ancestor); if (q != prunedToQFrag.end()) systems.insert(q->second);
      if (isBHadron((*pruned)[ancestor].pdgId())) {
        hadFromB[ih] = 1;
        const auto b = prunedToHadron.find(ancestor); if (b != prunedToHadron.end()) bAncestors.insert(b->second);
      }
    }
    if (systems.size() == 1) { hadQuark[ih] = *systems.begin(); hadAssocQ[ih] = kUnique; }
    else if (systems.size() > 1) hadAssocQ[ih] = kAmbiguous;
    if (bAncestors.size() == 1) hadParentB[ih] = *bAncestors.begin();
    if (isDHadron((*pruned)[hPruned[ih]].pdgId())) {
      // Absence of a retained B is not evidence of prompt production.  The
      // prompt class requires MiniAOD's explicit generator status flag.
      hadCharmOrigin[ih] = hadFromB[ih] ? kCharmFromB
                                        : ((*pruned)[hPruned[ih]].statusFlags().isPrompt() ? kPromptCharm : kCharmUnknown);
    }
  }

  std::vector<int> qNHadrons(qRep.size(), 0), qNCharged(qRep.size(), 0);
  std::vector<float> qChargedPt(qRep.size(), 0.f), qChargedPx(qRep.size(), 0.f), qChargedPy(qRep.size(), 0.f), qChargedPz(qRep.size(), 0.f), qChargedMass(qRep.size(), 0.f);
  for (const int q : hadQuark) if (q >= 0) ++qNHadrons[q];

  std::vector<int> cPacked, cPdg, cCharge, cQuark, cHadron, cLastPruned, cHFDecay, cDirectHFDecay, cFragmentation, cFromB, cFromPromptD, cFromBtoD, cMembership, cAmbiguous;
  std::vector<float> cPt, cEta, cPhi, cMass, cDR, cZ, cJT;
  for (size_t ip = 0; ip < packed->size(); ++ip) {
    const auto& p = (*packed)[ip];
    if (p.status() != 1 || p.charge() == 0) continue;
    int last = -1;
    const auto& ref = p.lastPrunedRef();
    if (ref.isNonnull() && ref.id() == pruned.id()) last = static_cast<int>(ref.key());
    Ancestors ancestry;
    if (last >= 0) ancestry = ancestors(last, *pruned, prunedIndex); else ancestry.incomplete = true;
    std::unordered_set<int> systems, hfHadrons;
    bool hasB = false, hasD = false, hasPromptD = false, hasBtoD = false;
    for (const int ancestor : ancestry.indices) {
      const auto qs = prunedToQFrag.find(ancestor); if (qs != prunedToQFrag.end()) systems.insert(qs->second);
      const auto hh = prunedToHadron.find(ancestor);
      if (hh != prunedToHadron.end()) hfHadrons.insert(hh->second);
      if (isBHadron((*pruned)[ancestor].pdgId())) hasB = true;
      if (isDHadron((*pruned)[ancestor].pdgId())) {
        hasD = true;
        const auto d = prunedToHadron.find(ancestor);
        if (d != prunedToHadron.end() && hadCharmOrigin[d->second] == kPromptCharm) hasPromptD = true;
      }
    }
    hasBtoD = hasB && hasD;
    int qidx = -1, hidx = -1, quality = kMissing;
    if (systems.size() == 1) { qidx = *systems.begin(); quality = kUnique; }
    else if (systems.size() > 1) quality = kAmbiguous;
    if (hfHadrons.size() == 1) hidx = *hfHadrons.begin();
    const bool hfDecay = hidx >= 0;
    const bool directHF = p.statusFlags().isDirectHadronDecayProduct();
    int membership = kUnassigned;
    if (quality == kAmbiguous || hfHadrons.size() > 1) membership = kAmbiguousMembership;
    else if (hfDecay) membership = kHeavyHadronDecay;
    else if (qidx >= 0) membership = kRetainedFragmentation;
    float dr = kInvalid, z = kInvalid, jt = kInvalid;
    if (qidx >= 0) {
      const auto& q = (*pruned)[qRep[qidx]];
      const float q2 = q.px()*q.px() + q.py()*q.py() + q.pz()*q.pz();
      if (q2 > 0.f) {
        const float dot = p.px()*q.px() + p.py()*q.py() + p.pz()*q.pz();
        z = dot/q2;
        jt = static_cast<float>(std::sqrt(std::max(0.0, p.p()*p.p() - dot*dot/q2)));
        dr = reco::deltaR(p.eta(), p.phi(), q.eta(), q.phi());
      }
      ++qNCharged[qidx]; qChargedPt[qidx] += p.pt(); qChargedPx[qidx] += p.px(); qChargedPy[qidx] += p.py(); qChargedPz[qidx] += p.pz(); qChargedMass[qidx] += p.mass();
    }
    cPacked.push_back(ip); cPdg.push_back(p.pdgId()); cCharge.push_back(p.charge()); cPt.push_back(p.pt()); cEta.push_back(p.eta()); cPhi.push_back(p.phi()); cMass.push_back(p.mass());
    cQuark.push_back(qidx); cHadron.push_back(hidx); cLastPruned.push_back(last); cHFDecay.push_back(hfDecay); cDirectHFDecay.push_back(directHF); cFragmentation.push_back(membership == kRetainedFragmentation);
    cFromB.push_back(hasB); cFromPromptD.push_back(hasPromptD); cFromBtoD.push_back(hasBtoD); cMembership.push_back(membership); cAmbiguous.push_back(quality == kAmbiguous || hfHadrons.size() > 1); cDR.push_back(dr); cZ.push_back(z); cJT.push_back(jt);
  }

  auto qTable = std::make_unique<nanoaod::FlatTable>(qRep.size(), "QFrag", false, false);
  std::vector<int> qPdg, qFlavour, qHardPdg, qIsHard, qFromHard, qIsLast, qIsLastBefore, qProd, qQuality;
  std::vector<float> qPt, qEta, qPhi, qMass, qEnergy, qHardPt, qHardEta, qHardPhi, qHardEnergy, qZCharged;
  for (size_t iq = 0; iq < qRep.size(); ++iq) {
    const auto& q = (*pruned)[qRep[iq]]; const auto& flags = q.statusFlags();
    qPdg.push_back(q.pdgId()); qFlavour.push_back(flavour(q.pdgId())); qIsHard.push_back(flags.isHardProcess()); qFromHard.push_back(flags.fromHardProcess()); qIsLast.push_back(flags.isLastCopy()); qIsLastBefore.push_back(flags.isLastCopyBeforeFSR());
    qPt.push_back(q.pt()); qEta.push_back(q.eta()); qPhi.push_back(q.phi()); qMass.push_back(q.mass()); qEnergy.push_back(q.energy());
    if (qHard[iq] >= 0) { const auto& h = (*pruned)[qHard[iq]]; qHardPt.push_back(h.pt()); qHardEta.push_back(h.eta()); qHardPhi.push_back(h.phi()); qHardEnergy.push_back(h.energy()); }
    else { qHardPt.push_back(kInvalid); qHardEta.push_back(kInvalid); qHardPhi.push_back(kInvalid); qHardEnergy.push_back(kInvalid); }
    qProd.push_back(flags.isHardProcess() ? 1 : (flags.fromHardProcess() ? 2 : 0)); qQuality.push_back(qHard[iq] >= 0 || flags.isLastCopy() ? kUnique : kMissing);
    const float q2 = q.p()*q.p(); qZCharged.push_back(q2 > 0 ? (qChargedPx[iq]*q.px()+qChargedPy[iq]*q.py()+qChargedPz[iq]*q.pz())/q2 : kInvalid);
  }
  qTable->addColumn<int>("pdgId", qPdg, "signed PDG identifier"); qTable->addColumn<int>("flavour", qFlavour, "absolute quark flavour"); qTable->addColumn<int>("prunedGenIdx", qRep, "representative terminal pruned quark index"); qTable->addColumn<int>("hardQuarkIdx", qHard, "unique hard-process same-copy-chain index, or -1"); qTable->addColumn<int>("lastQuarkIdx", qLast, "terminal same-copy-chain index");
  qTable->addColumn<float>("pt", qPt, "representative transverse momentum"); qTable->addColumn<float>("eta", qEta, "representative eta"); qTable->addColumn<float>("phi", qPhi, "representative phi"); qTable->addColumn<float>("mass", qMass, "representative mass"); qTable->addColumn<float>("energy", qEnergy, "representative energy"); qTable->addColumn<float>("hardPt", qHardPt, "hard-copy pt or -999"); qTable->addColumn<float>("hardEta", qHardEta, "hard-copy eta or -999"); qTable->addColumn<float>("hardPhi", qHardPhi, "hard-copy phi or -999"); qTable->addColumn<float>("hardEnergy", qHardEnergy, "hard-copy energy or -999");
  qTable->addColumn<int>("isHardProcess", qIsHard, "status flag"); qTable->addColumn<int>("fromHardProcess", qFromHard, "status flag"); qTable->addColumn<int>("isLastCopy", qIsLast, "status flag"); qTable->addColumn<int>("isLastCopyBeforeFSR", qIsLastBefore, "status flag"); qTable->addColumn<int>("productionType", qProd, "0 unknown, 1 hard process, 2 from hard process"); qTable->addColumn<int>("assignmentQuality", qQuality, "0 missing, 1 retained unique"); qTable->addColumn<int>("nHeavyHadrons", qNHadrons, "uniquely associated selected B/D hadrons"); qTable->addColumn<int>("nGenCharged", qNCharged, "uniquely associated packed charged particles"); qTable->addColumn<float>("chargedPtSum", qChargedPt, "scalar charged pt sum"); qTable->addColumn<float>("chargedPx", qChargedPx, "charged px sum"); qTable->addColumn<float>("chargedPy", qChargedPy, "charged py sum"); qTable->addColumn<float>("chargedPz", qChargedPz, "charged pz sum"); qTable->addColumn<float>("chargedMass", qChargedMass, "sum of charged particle masses"); qTable->addColumn<float>("zCharged", qZCharged, "charged momentum projection or -999");

  auto hTable = std::make_unique<nanoaod::FlatTable>(hPruned.size(), "QHadron", false, false);
  std::vector<int> hPdg, hFlavor, hPrompt, hLeading, hNCharged(hPruned.size(), 0), hNMatched(hPruned.size(), -1), hType(hPruned.size(), kNoAssociation); std::vector<float> hPt, hEta, hPhi, hMass, hEnergy, hZ(hPruned.size(), kInvalid), hDR(hPruned.size(), kInvalid);
  for (const int h : cHadron) if (h >= 0) ++hNCharged[h];
  for (size_t ih = 0; ih < hPruned.size(); ++ih) { const auto& h = (*pruned)[hPruned[ih]]; hPdg.push_back(h.pdgId()); hFlavor.push_back(hadronFlavour(h.pdgId())); hPrompt.push_back(hadCharmOrigin[ih] == kPromptCharm); hLeading.push_back(0); hPt.push_back(h.pt()); hEta.push_back(h.eta()); hPhi.push_back(h.phi()); hMass.push_back(h.mass()); hEnergy.push_back(h.energy()); if (hadQuark[ih] >= 0) { const auto& q = (*pruned)[qRep[hadQuark[ih]]]; const float q2=q.p()*q.p(); if(q2>0){hZ[ih]=(h.px()*q.px()+h.py()*q.py()+h.pz()*q.pz())/q2; hDR[ih]=reco::deltaR(h.eta(),h.phi(),q.eta(),q.phi()); hType[ih]=kRetainedAncestry;} } }
  hTable->addColumn<int>("quarkIdx", hadQuark, "QFrag row, or -1"); hTable->addColumn<int>("prunedGenIdx", hPruned, "pruned generator index"); hTable->addColumn<int>("pdgId", hPdg, "signed PDG identifier"); hTable->addColumn<int>("flavour", hFlavor, "5 B, 4 D"); hTable->addColumn<float>("pt",hPt,"pt"); hTable->addColumn<float>("eta",hEta,"eta"); hTable->addColumn<float>("phi",hPhi,"phi"); hTable->addColumn<float>("mass",hMass,"mass"); hTable->addColumn<float>("energy",hEnergy,"energy"); hTable->addColumn<int>("isPromptCharm",hPrompt,"1 confirmed prompt D, 0 otherwise"); hTable->addColumn<int>("fromB",hadFromB,"retained B-hadron ancestry found"); hTable->addColumn<int>("parentBHadronIdx",hadParentB,"nearest unique selected B QHadron row, or -1"); hTable->addColumn<int>("isLeadingHFHadron",hLeading,"reserved; 0 in Stage 1"); hTable->addColumn<float>("zHadron",hZ,"quark projection or -999"); hTable->addColumn<float>("deltaRQuark",hDR,"quark deltaR or -999"); hTable->addColumn<int>("nChargedDaughters",hNCharged,"packed charged particles with this unique nearest heavy-hadron ancestor"); hTable->addColumn<int>("nMatchedDaughters",hNMatched,"-1: no reco matching in Stage 1"); hTable->addColumn<int>("associationType",hType,"0 none, 1 retained ancestry"); hTable->addColumn<int>("associationQuality",hadAssocQ,"0 missing, 1 unique, 2 ambiguous"); hTable->addColumn<int>("charmOrigin",hadCharmOrigin,"0 n/a, 1 prompt, 2 from B, 3 unknown");

  auto cTable = std::make_unique<nanoaod::FlatTable>(cPacked.size(), "QGenCharged", false, false);
  cTable->addColumn<int>("packedGenIdx",cPacked,"packedGenParticles index"); cTable->addColumn<int>("pdgId",cPdg,"signed PDG identifier"); cTable->addColumn<int>("charge",cCharge,"electric charge"); cTable->addColumn<float>("pt",cPt,"pt"); cTable->addColumn<float>("eta",cEta,"eta"); cTable->addColumn<float>("phi",cPhi,"phi"); cTable->addColumn<float>("mass",cMass,"mass"); cTable->addColumn<int>("quarkIdx",cQuark,"QFrag row or -1"); cTable->addColumn<int>("hadronIdx",cHadron,"QHadron row or -1"); cTable->addColumn<int>("lastPrunedIdx",cLastPruned,"lastPrunedRef key or -1"); cTable->addColumn<int>("isHFDecay",cHFDecay,"unique selected heavy-hadron ancestor"); cTable->addColumn<int>("isDirectHFDecay",cDirectHFDecay,"generator direct-hadron-decay flag"); cTable->addColumn<int>("isFragmentation",cFragmentation,"retained quark ancestry outside selected HF decay"); cTable->addColumn<int>("isFromB",cFromB,"retained B ancestor"); cTable->addColumn<int>("isFromPromptD",cFromPromptD,"retained confirmed-prompt D ancestor"); cTable->addColumn<int>("isFromBtoD",cFromBtoD,"both B and D retained ancestors"); cTable->addColumn<int>("membershipType",cMembership,"0 unassigned, 1 HF decay, 2 retained fragmentation, 3 ambiguous"); cTable->addColumn<int>("isAmbiguous",cAmbiguous,"multiple retained systems or HF ancestors"); cTable->addColumn<float>("deltaRQuark",cDR,"quark deltaR or -999"); cTable->addColumn<float>("zTrackQuark",cZ,"quark projection or -999"); cTable->addColumn<float>("jTQuark",cJT,"transverse momentum to quark or -999");
  event.put(std::move(qTable), "QFragTable"); event.put(std::move(hTable), "QHadronTable"); event.put(std::move(cTable), "QGenChargedTable");
}

DEFINE_FWK_MODULE(QuarkFragmentationProducer);
