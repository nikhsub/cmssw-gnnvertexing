// Generator-only MiniAOD truth producer for quark-fragmentation studies.
//
// Stage 1 deliberately makes no angular or momentum based assignments.  An
// object receives a quark-system index only when an unambiguous retained
// mother-chain path reaches one selected pruned quark system.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
  constexpr float kInvalid = -999.f;
  enum AssociationType { kNoAssociation = 0, kRetainedAncestry = 1 };
  enum AssociationQuality { kMissing = 0, kUnique = 1, kAmbiguous = 2 };
  enum MembershipType { kUnassigned = 0, kHeavyHadronDecay = 1, kRetainedFragmentation = 2, kAmbiguousMembership = 3 };
  enum CharmOrigin { kCharmNotApplicable = 0, kPromptCharm = 1, kCharmFromB = 2, kCharmUnknown = 3 };
  // Diagnostic reason codes describe current algorithm decisions only.
  enum AssociationReason { kReasonUniqueSystem = 1, kReasonMultipleSystems = 2, kReasonNoQuarkAncestor = 3, kReasonNoSelectedSystem = 4 };
  enum OriginReason { kOriginNotD = 0, kOriginRetainedB = 1, kOriginPromptFlag = 2, kOriginNoRetainedEvidence = 3, kOriginDepthGuard = 4, kOriginMissingLink = 5 };
  enum MembershipReason { kReasonHFDecay = 1, kReasonRetainedFragmentation = 2, kReasonAmbiguousMembership = 3, kReasonUnassigned = 4 };
  enum RepresentativeReason { kReasonUniqueTerminalCopy = 1 };

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
    unsigned int maxDepth = 0;
    unsigned int multipleMotherNodes = 0;
    bool incomplete = false;  // null or non-pruned reference encountered
    bool truncated = false;   // traversal depth guard reached
    bool cycle = false;
  };
}

class QuarkFragmentationProducer : public edm::stream::EDProducer<> {
public:
  explicit QuarkFragmentationProducer(const edm::ParameterSet& cfg)
      : prunedToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
        packedToken_(consumes<pat::PackedGenParticleCollection>(cfg.getParameter<edm::InputTag>("packedGenParticles"))),
        debugAncestry_(cfg.getParameter<bool>("debugAncestry")),
        debugOnlyAmbiguous_(cfg.getParameter<bool>("debugOnlyAmbiguous")),
        debugOnlyD_(cfg.getParameter<bool>("debugOnlyD")),
        debugMaxEvents_(cfg.getParameter<unsigned int>("debugMaxEvents")),
        debugMaxHadronsPerEvent_(cfg.getParameter<unsigned int>("debugMaxHadronsPerEvent")),
        debugMaxDepth_(cfg.getParameter<unsigned int>("debugMaxDepth")) {
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
  void logAncestry(const edm::Event& event, int hadronRow, int hadronIndex,
                   const reco::GenParticleCollection& particles,
                   const std::unordered_map<const reco::Candidate*, int>& lookup,
                   const std::vector<int>& candidateCopies, const Ancestors& ancestry, int quarkRow,
                   int associationQuality, int associationReason, int charmOrigin,
                   int originReason, unsigned int depthLimit) const;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedToken_;
  bool debugAncestry_;
  bool debugOnlyAmbiguous_;
  bool debugOnlyD_;
  unsigned int debugMaxEvents_;
  unsigned int debugMaxHadronsPerEvent_;
  unsigned int debugMaxDepth_;
  unsigned int debugEventsSeen_ = 0;
};

Ancestors QuarkFragmentationProducer::ancestors(
    int start, const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup) const {
  Ancestors result;
  if (start < 0 || static_cast<size_t>(start) >= particles.size()) { result.incomplete = true; return result; }
  constexpr unsigned int kTraversalDepthGuard = 256;
  std::queue<std::pair<int, unsigned int>> todo;
  std::unordered_set<int> seen;
  todo.push({start, 0});
  seen.insert(start);
  while (!todo.empty()) {
    const int current = todo.front().first;
    const unsigned int depth = todo.front().second;
    todo.pop();
    result.maxDepth = std::max(result.maxDepth, depth);
    const auto& particle = particles[current];
    if (particle.numberOfMothers() > 1) ++result.multipleMotherNodes;
    for (size_t im = 0; im < particle.numberOfMothers(); ++im) {
      const reco::Candidate* mother = particle.mother(im);
      if (!mother) { result.incomplete = true; continue; }
      const auto found = lookup.find(mother);
      if (found == lookup.end()) { result.incomplete = true; continue; }
      if (seen.insert(found->second).second) {
        result.indices.push_back(found->second);
        todo.push({found->second, depth + 1});
      }
    }
  }
  std::unordered_set<int> path;
  std::unordered_set<int> completed;
  std::function<void(int, unsigned int)> detectCycle = [&](int current, unsigned int depth) {
    if (depth >= kTraversalDepthGuard) { result.truncated = true; return; }
    if (completed.count(current)) return;
    path.insert(current);
    const auto& particle = particles[current];
    for (size_t im = 0; im < particle.numberOfMothers(); ++im) {
      const auto found = lookup.find(particle.mother(im));
      if (found == lookup.end()) continue;
      if (path.count(found->second)) { result.cycle = true; continue; }
      detectCycle(found->second, depth + 1);
    }
    path.erase(current);
    completed.insert(current);
  };
  detectCycle(start, 0);
  return result;
}

void QuarkFragmentationProducer::logAncestry(
    const edm::Event& event, int hadronRow, int hadronIndex,
    const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup,
    const std::vector<int>& candidateCopies, const Ancestors& ancestry, int quarkRow,
    int associationQuality, int associationReason, int charmOrigin,
    int originReason, unsigned int depthLimit) const {
  const auto& h = particles[hadronIndex];
  auto log = edm::LogInfo("QFragAncestry");
  log << "Event " << event.id() << " QHadron=" << hadronRow << " prunedGenIdx=" << hadronIndex
      << " pdgId=" << h.pdgId() << " status=" << h.status() << " pt=" << h.pt()
      << " eta=" << h.eta() << " phi=" << h.phi() << "\nMother tree:";
  std::unordered_set<int> path;
  std::unordered_set<int> visited;
  std::function<void(int, unsigned int)> visit = [&](int index, unsigned int depth) {
    const auto& p = particles[index];
    log << "\n" << std::string(depth * 2, ' ') << "prunedGenIdx=" << index << " pdgId=" << p.pdgId()
        << " status=" << p.status() << " flags=" << p.statusFlags().flags_ << " pt=" << p.pt()
        << " eta=" << p.eta() << " phi=" << p.phi() << " nMothers=" << p.numberOfMothers();
    if (depth >= depthLimit) { log << " [debug depth limit]"; return; }
    if (visited.count(index)) { log << " [shared ancestry node]"; return; }
    visited.insert(index);
    path.insert(index);
    for (size_t im = 0; im < p.numberOfMothers(); ++im) {
      const auto* mother = p.mother(im);
      const auto found = lookup.find(mother);
      if (!mother || found == lookup.end()) { log << "\n" << std::string((depth + 1) * 2, ' ') << "[null or unavailable mother reference]"; continue; }
      if (path.count(found->second)) { log << "\n" << std::string((depth + 1) * 2, ' ') << "[cycle to prunedGenIdx=" << found->second << "]"; continue; }
      visit(found->second, depth + 1);
    }
    path.erase(index);
  };
  visit(hadronIndex, 0);
  unsigned int nCharm = 0, nBottom = 0, nLight = 0;
  log << "\nCandidate quarks:";
  for (int idx : candidateCopies) {
    const auto& q = particles[idx];
    const int f = flavour(q.pdgId());
    if (f == 4) ++nCharm; else if (f == 5) ++nBottom; else ++nLight;
    log << "\n  prunedGenIdx=" << idx << " pdgId=" << q.pdgId() << " flags=" << q.statusFlags().flags_
        << " pt=" << q.pt() << " eta=" << q.eta() << " phi=" << q.phi();
  }
  log << "\nnumberOfCandidateQuarkCopies=" << candidateCopies.size() << " charmCopies=" << nCharm
      << " bottomCopies=" << nBottom << " lightCopies=" << nLight
      << " distinctIdentifiableLineages=unknown nAvailableMothers=" << h.numberOfMothers()
      << " ancestryDepth=" << ancestry.maxDepth << " hasMultipleMotherBranches=" << (ancestry.multipleMotherNodes > 0)
      << " ancestryTruncated=" << ancestry.truncated << " ancestryIncomplete=" << ancestry.incomplete
      << " ancestryCycleDetected=" << ancestry.cycle
      << " quarkIdx=" << quarkRow << " associationQuality=" << associationQuality
      << " associationReason=" << associationReason << " charmOrigin=" << charmOrigin
      << " originReason=" << originReason;
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
  std::vector<int> qRep, qHard, qLast, qFirst, qCopies, qMothers, qSameMother, qSameDaughter, qCopyDepth, qCopyAmbiguous, qSelectionReason;
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
    std::unordered_set<int> componentSet(component.begin(), component.end());
    std::vector<int> roots;
    bool branching = false;
    for (const int member : component) {
      const auto& p = (*pruned)[member];
      bool sameMother = false;
      for (size_t im = 0; im < p.numberOfMothers(); ++im) {
        const auto* m = p.mother(im);
        const auto mi = prunedIndex.find(m);
        if (m && mi != prunedIndex.end() && m->pdgId() == p.pdgId() && componentSet.count(mi->second)) sameMother = true;
      }
      if (!sameMother) roots.push_back(member);
      unsigned int sameDaughters = 0;
      for (size_t id = 0; id < p.numberOfDaughters(); ++id) {
        const auto* d = p.daughter(id);
        const auto di = prunedIndex.find(d);
        if (d && di != prunedIndex.end() && d->pdgId() == p.pdgId() && componentSet.count(di->second)) ++sameDaughters;
      }
      if (sameDaughters > 1) branching = true;
    }
    qFirst.push_back(roots.size() == 1 ? roots.front() : -1);
    qCopies.push_back(component.size()); qMothers.push_back((*pruned)[representative].numberOfMothers());
    const auto& rep = (*pruned)[representative];
    bool repSameMother = false, repSameDaughter = false;
    for (size_t im = 0; im < rep.numberOfMothers(); ++im) { const auto* m=rep.mother(im); repSameMother |= m && m->pdgId()==rep.pdgId(); }
    for (size_t id = 0; id < rep.numberOfDaughters(); ++id) { const auto* d=rep.daughter(id); repSameDaughter |= d && d->pdgId()==rep.pdgId(); }
    qSameMother.push_back(repSameMother); qSameDaughter.push_back(repSameDaughter);
    qCopyDepth.push_back(roots.size() == 1 && !branching ? static_cast<int>(component.size() - 1) : -1);
    qCopyAmbiguous.push_back(branching || roots.size() != 1);
    qSelectionReason.push_back(kReasonUniqueTerminalCopy);  // exactly one terminal same-PDG copy in the retained component
    for (const int member : component) prunedToQFrag[member] = row;
  }

  std::vector<int> hPruned;
  for (size_t i = 0; i < pruned->size(); ++i)
    if (isHeavyHadron((*pruned)[i].pdgId()) && isFinalCopy(i, *pruned, prunedIndex)) hPruned.push_back(i);
  std::unordered_map<int, int> prunedToHadron;
  for (size_t i = 0; i < hPruned.size(); ++i) prunedToHadron[hPruned[i]] = static_cast<int>(i);

  std::vector<int> hadQuark(hPruned.size(), -1), hadAssocQ(hPruned.size(), kMissing), hadFromB(hPruned.size(), 0), hadParentB(hPruned.size(), -1);
  std::vector<int> hadCharmOrigin(hPruned.size(), kCharmNotApplicable);
  std::vector<int> hNCandidates(hPruned.size(), 0), hNCandidateC(hPruned.size(), 0), hNCandidateB(hPruned.size(), 0), hNCandidateLight(hPruned.size(), 0), hNDistinctLineages(hPruned.size(), -1), hNCandidateCopies(hPruned.size(), 0);
  std::vector<int> hMultiCandidates(hPruned.size(), 0), hConflictingFlavours(hPruned.size(), 0), hNMothers(hPruned.size(), 0), hHasB(hPruned.size(), 0), hHasD(hPruned.size(), 0), hHasCQuark(hPruned.size(), 0), hHasBQuark(hPruned.size(), 0), hMultipleBranches(hPruned.size(), 0), hAncestryDepth(hPruned.size(), 0), hAncestryTruncated(hPruned.size(), 0), hAncestryIncomplete(hPruned.size(), 0), hAncestryExhausted(hPruned.size(), 0), hCycle(hPruned.size(), 0), hAssocReason(hPruned.size(), kReasonNoQuarkAncestor), hOriginReason(hPruned.size(), kOriginNotD), hCandidateSystems(hPruned.size(), 0);
  std::vector<std::vector<int>> hadCandidateCopies(hPruned.size());
  ++debugEventsSeen_;
  unsigned int debuggedHadrons = 0;
  for (size_t ih = 0; ih < hPruned.size(); ++ih) {
    const auto ancestry = ancestors(hPruned[ih], *pruned, prunedIndex);
    std::unordered_set<int> systems, bAncestors;
    std::unordered_set<int> candidateFlavours;
    for (const int ancestor : ancestry.indices) {
      const auto q = prunedToQFrag.find(ancestor); if (q != prunedToQFrag.end()) systems.insert(q->second);
      if (isQuark((*pruned)[ancestor].pdgId())) {
        hadCandidateCopies[ih].push_back(ancestor);
      }
      if (isBHadron((*pruned)[ancestor].pdgId())) {
        hadFromB[ih] = 1;
        hHasB[ih] = 1;
        const auto b = prunedToHadron.find(ancestor); if (b != prunedToHadron.end()) bAncestors.insert(b->second);
      }
      if (isDHadron((*pruned)[ancestor].pdgId())) hHasD[ih] = 1;
    }
    for (const int system : systems) candidateFlavours.insert(flavour((*pruned)[qRep[system]].pdgId()));
    hNCandidates[ih] = systems.size(); hNCandidateCopies[ih] = hadCandidateCopies[ih].size();
    for (const int f : candidateFlavours) {
      if (f == 4) ++hNCandidateC[ih];
      else if (f == 5) ++hNCandidateB[ih];
      else ++hNCandidateLight[ih];
    }
    hMultiCandidates[ih] = hNCandidates[ih] > 1; hConflictingFlavours[ih] = candidateFlavours.size() > 1;
    hCandidateSystems[ih] = systems.size(); hNMothers[ih] = (*pruned)[hPruned[ih]].numberOfMothers();
    hMultipleBranches[ih] = ancestry.multipleMotherNodes > 0; hAncestryDepth[ih] = ancestry.maxDepth;
    hAncestryTruncated[ih] = ancestry.truncated; hAncestryIncomplete[ih] = ancestry.incomplete;
    hAncestryExhausted[ih] = !ancestry.truncated && !ancestry.incomplete; hCycle[ih] = ancestry.cycle;
    for (int idx : hadCandidateCopies[ih]) {
      const int f = flavour((*pruned)[idx].pdgId()); hHasCQuark[ih] |= f == 4; hHasBQuark[ih] |= f == 5;
    }
    if (systems.size() == 1) { hadQuark[ih] = *systems.begin(); hadAssocQ[ih] = kUnique; }
    else if (systems.size() > 1) hadAssocQ[ih] = kAmbiguous;
    hAssocReason[ih] = systems.size() == 1 ? kReasonUniqueSystem
                           : (systems.size() > 1 ? kReasonMultipleSystems
                              : (hadCandidateCopies[ih].empty() ? kReasonNoQuarkAncestor : kReasonNoSelectedSystem));
    if (bAncestors.size() == 1) hadParentB[ih] = *bAncestors.begin();
    if (isDHadron((*pruned)[hPruned[ih]].pdgId())) {
      // Absence of a retained B is not evidence of prompt production.  The
      // prompt class requires MiniAOD's explicit generator status flag.
      hadCharmOrigin[ih] = hadFromB[ih] ? kCharmFromB
                                        : ((*pruned)[hPruned[ih]].statusFlags().isPrompt() ? kPromptCharm : kCharmUnknown);
      hOriginReason[ih] = hadFromB[ih] ? kOriginRetainedB
                           : ((*pruned)[hPruned[ih]].statusFlags().isPrompt() ? kOriginPromptFlag
                              : (ancestry.truncated ? kOriginDepthGuard : (ancestry.incomplete ? kOriginMissingLink : kOriginNoRetainedEvidence)));
    }
    if (debugAncestry_ && debugEventsSeen_ <= debugMaxEvents_ && debuggedHadrons < debugMaxHadronsPerEvent_ &&
        (!debugOnlyAmbiguous_ || hadAssocQ[ih] == kAmbiguous) &&
        (!debugOnlyD_ || isDHadron((*pruned)[hPruned[ih]].pdgId()))) {
      logAncestry(event, ih, hPruned[ih], *pruned, prunedIndex, hadCandidateCopies[ih], ancestry, hadQuark[ih],
                  hadAssocQ[ih], hAssocReason[ih], hadCharmOrigin[ih], hOriginReason[ih], debugMaxDepth_);
      ++debuggedHadrons;
    }
  }

  std::vector<int> qNHadrons(qRep.size(), 0), qNCharged(qRep.size(), 0);
  std::vector<float> qChargedPt(qRep.size(), 0.f), qChargedPx(qRep.size(), 0.f), qChargedPy(qRep.size(), 0.f), qChargedPz(qRep.size(), 0.f), qChargedMass(qRep.size(), 0.f);
  for (const int q : hadQuark) if (q >= 0) ++qNHadrons[q];

  std::vector<int> cPacked, cPdg, cCharge, cQuark, cHadron, cLastPruned, cLastPdg, cLastStatus, cLastSelf, cHFDecay, cDirectHFDecay, cFragmentation, cFromB, cFromPromptD, cFromBtoD, cMembership, cAmbiguous, cNCandidateSystems, cNCandidateHadrons, cHasBAncestor, cHasDAncestor, cMultipleBranches, cAncestryTruncated, cMembershipReason;
  std::vector<float> cPt, cEta, cPhi, cMass, cDR, cZ, cJT;
  for (size_t ip = 0; ip < packed->size(); ++ip) {
    const auto& p = (*packed)[ip];
    if (p.status() != 1 || p.charge() == 0) continue;
    int last = -1;
    const auto& ref = p.lastPrunedRef();
    if (ref.isNonnull() && ref.isAvailable() && ref.id() == pruned.id()) last = static_cast<int>(ref.key());
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
    // PackedGenParticle::motherRef() treats a stored status-1 pruned reference
    // as the packed particle's own pruned representation and steps past it.
    const bool lastSelf = last >= 0 && (*pruned)[last].status() == 1;
    cPacked.push_back(ip); cPdg.push_back(p.pdgId()); cCharge.push_back(p.charge()); cPt.push_back(p.pt()); cEta.push_back(p.eta()); cPhi.push_back(p.phi()); cMass.push_back(p.mass());
    cQuark.push_back(qidx); cHadron.push_back(hidx); cLastPruned.push_back(last); cHFDecay.push_back(hfDecay); cDirectHFDecay.push_back(directHF); cFragmentation.push_back(membership == kRetainedFragmentation);
    cLastPdg.push_back(last >= 0 ? (*pruned)[last].pdgId() : -999); cLastStatus.push_back(last >= 0 ? (*pruned)[last].status() : -1); cLastSelf.push_back(lastSelf);
    cFromB.push_back(hasB); cFromPromptD.push_back(hasPromptD); cFromBtoD.push_back(hasBtoD); cMembership.push_back(membership); cAmbiguous.push_back(quality == kAmbiguous || hfHadrons.size() > 1);
    cNCandidateSystems.push_back(systems.size()); cNCandidateHadrons.push_back(hfHadrons.size()); cHasBAncestor.push_back(hasB); cHasDAncestor.push_back(hasD); cMultipleBranches.push_back(ancestry.multipleMotherNodes > 0); cAncestryTruncated.push_back(ancestry.truncated);
    cMembershipReason.push_back(membership == kHeavyHadronDecay ? kReasonHFDecay
                                  : (membership == kRetainedFragmentation ? kReasonRetainedFragmentation
                                     : (membership == kAmbiguousMembership ? kReasonAmbiguousMembership : kReasonUnassigned)));
    cDR.push_back(dr); cZ.push_back(z); cJT.push_back(jt);
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
  qTable->addColumn<int>("firstQuarkIdx", qFirst, "unique retained same-PDG chain root, or -1"); qTable->addColumn<int>("nRetainedQuarkCopies", qCopies, "number of pruned records in the connected same-PDG component"); qTable->addColumn<int>("nAvailableMothers", qMothers, "mother count on representative record"); qTable->addColumn<int>("hasSameFlavourMother", qSameMother, "representative has a same-PDG mother"); qTable->addColumn<int>("hasSameFlavourDaughter", qSameDaughter, "representative has a same-PDG daughter"); qTable->addColumn<int>("copyChainDepth", qCopyDepth, "component size minus one for an unbranched chain; -1 if ambiguous"); qTable->addColumn<int>("copyChainAmbiguous", qCopyAmbiguous, "branching or non-unique retained root"); qTable->addColumn<int>("representativeSelectionReason", qSelectionReason, "1: unique terminal record in connected same-PDG component"); qTable->addColumn<int>("nDistinctQuarkLineages", std::vector<int>(qRep.size(), -1), "-1: physical lineage cannot be established from retained copies alone");

  auto hTable = std::make_unique<nanoaod::FlatTable>(hPruned.size(), "QHadron", false, false);
  std::vector<int> hPdg, hFlavor, hPrompt, hLeading, hNCharged(hPruned.size(), 0), hNMatched(hPruned.size(), -1), hType(hPruned.size(), kNoAssociation); std::vector<float> hPt, hEta, hPhi, hMass, hEnergy, hZ(hPruned.size(), kInvalid), hDR(hPruned.size(), kInvalid);
  for (const int h : cHadron) if (h >= 0) ++hNCharged[h];
  for (size_t ih = 0; ih < hPruned.size(); ++ih) { const auto& h = (*pruned)[hPruned[ih]]; hPdg.push_back(h.pdgId()); hFlavor.push_back(hadronFlavour(h.pdgId())); hPrompt.push_back(hadCharmOrigin[ih] == kPromptCharm); hLeading.push_back(0); hPt.push_back(h.pt()); hEta.push_back(h.eta()); hPhi.push_back(h.phi()); hMass.push_back(h.mass()); hEnergy.push_back(h.energy()); if (hadQuark[ih] >= 0) { const auto& q = (*pruned)[qRep[hadQuark[ih]]]; const float q2=q.p()*q.p(); if(q2>0){hZ[ih]=(h.px()*q.px()+h.py()*q.py()+h.pz()*q.pz())/q2; hDR[ih]=reco::deltaR(h.eta(),h.phi(),q.eta(),q.phi()); hType[ih]=kRetainedAncestry;} } }
  hTable->addColumn<int>("quarkIdx", hadQuark, "QFrag row, or -1"); hTable->addColumn<int>("prunedGenIdx", hPruned, "pruned generator index"); hTable->addColumn<int>("pdgId", hPdg, "signed PDG identifier"); hTable->addColumn<int>("flavour", hFlavor, "5 B, 4 D"); hTable->addColumn<float>("pt",hPt,"pt"); hTable->addColumn<float>("eta",hEta,"eta"); hTable->addColumn<float>("phi",hPhi,"phi"); hTable->addColumn<float>("mass",hMass,"mass"); hTable->addColumn<float>("energy",hEnergy,"energy"); hTable->addColumn<int>("isPromptCharm",hPrompt,"1 confirmed prompt D, 0 otherwise"); hTable->addColumn<int>("fromB",hadFromB,"retained B-hadron ancestry found"); hTable->addColumn<int>("parentBHadronIdx",hadParentB,"nearest unique selected B QHadron row, or -1"); hTable->addColumn<int>("isLeadingHFHadron",hLeading,"reserved; 0 in Stage 1"); hTable->addColumn<float>("zHadron",hZ,"quark projection or -999"); hTable->addColumn<float>("deltaRQuark",hDR,"quark deltaR or -999"); hTable->addColumn<int>("nChargedDaughters",hNCharged,"packed charged particles with this unique nearest heavy-hadron ancestor"); hTable->addColumn<int>("nMatchedDaughters",hNMatched,"-1: no reco matching in Stage 1"); hTable->addColumn<int>("associationType",hType,"0 none, 1 retained ancestry"); hTable->addColumn<int>("associationQuality",hadAssocQ,"0 missing, 1 unique, 2 ambiguous"); hTable->addColumn<int>("charmOrigin",hadCharmOrigin,"0 n/a, 1 prompt, 2 from B, 3 unknown");
  hTable->addColumn<int>("nCandidateQuarks",hNCandidates,"number of existing QFrag systems reached in retained ancestry"); hTable->addColumn<int>("nCandidateCharm",hNCandidateC,"candidate QFrag systems with charm flavour"); hTable->addColumn<int>("nCandidateBottom",hNCandidateB,"candidate QFrag systems with bottom flavour"); hTable->addColumn<int>("nCandidateLight",hNCandidateLight,"candidate QFrag systems with u/d/s flavour"); hTable->addColumn<int>("nDistinctQuarkLineages",hNDistinctLineages,"-1: physical lineages are not identifiable reliably"); hTable->addColumn<int>("nCandidateQuarkCopies",hNCandidateCopies,"number of distinct retained ancestor quark records before QFrag copy grouping"); hTable->addColumn<int>("nCandidateQuarkSystems",hCandidateSystems,"number of existing QFrag copy components among candidate ancestors"); hTable->addColumn<int>("hasMultipleQuarkCandidates",hMultiCandidates,"more than one QFrag system reached"); hTable->addColumn<int>("hasConflictingQuarkFlavours",hConflictingFlavours,"candidate systems have multiple absolute flavours");
  hTable->addColumn<int>("nAvailableMothers",hNMothers,"direct mother count of this hadron record"); hTable->addColumn<int>("hasBHadronAncestor",hHasB,"retained B-hadron ancestor"); hTable->addColumn<int>("hasCHadronAncestor",hHasD,"retained D-hadron ancestor"); hTable->addColumn<int>("hasCharmQuarkAncestor",hHasCQuark,"retained charm quark ancestor"); hTable->addColumn<int>("hasBottomQuarkAncestor",hHasBQuark,"retained bottom quark ancestor"); hTable->addColumn<int>("hasMultipleMotherBranches",hMultipleBranches,"an ancestor record has more than one mother"); hTable->addColumn<int>("ancestryDepth",hAncestryDepth,"maximum discovered breadth-first retained ancestry depth"); hTable->addColumn<int>("ancestryTruncated",hAncestryTruncated,"depth guard reached"); hTable->addColumn<int>("ancestryIncomplete",hAncestryIncomplete,"null or non-pruned mother reference encountered"); hTable->addColumn<int>("ancestryExhausted",hAncestryExhausted,"retained traversal ended without detected truncation/reference gap; not proof physical history is complete"); hTable->addColumn<int>("ancestryCycleDetected",hCycle,"cycle found in retained mother graph");
  hTable->addColumn<int>("associationReason",hAssocReason,"1 unique retained QFrag component, 2 multiple components, 3 no quark ancestor, 4 quark ancestor(s) found but none mapped to a selected QFrag row"); hTable->addColumn<int>("originReason",hOriginReason,"0 not D, 1 retained B ancestor, 2 explicit prompt status flag, 3 no retained B or prompt flag, 4 traversal depth-limited, 5 missing retained link");

  auto cTable = std::make_unique<nanoaod::FlatTable>(cPacked.size(), "QGenCharged", false, false);
  cTable->addColumn<int>("packedGenIdx",cPacked,"packedGenParticles index"); cTable->addColumn<int>("pdgId",cPdg,"signed PDG identifier"); cTable->addColumn<int>("charge",cCharge,"electric charge"); cTable->addColumn<float>("pt",cPt,"pt"); cTable->addColumn<float>("eta",cEta,"eta"); cTable->addColumn<float>("phi",cPhi,"phi"); cTable->addColumn<float>("mass",cMass,"mass"); cTable->addColumn<int>("quarkIdx",cQuark,"QFrag row or -1"); cTable->addColumn<int>("hadronIdx",cHadron,"QHadron row or -1"); cTable->addColumn<int>("lastPrunedIdx",cLastPruned,"lastPrunedRef key or -1"); cTable->addColumn<int>("isHFDecay",cHFDecay,"unique selected heavy-hadron ancestor"); cTable->addColumn<int>("isDirectHFDecay",cDirectHFDecay,"generator direct-hadron-decay flag"); cTable->addColumn<int>("isFragmentation",cFragmentation,"retained quark ancestry outside selected HF decay"); cTable->addColumn<int>("isFromB",cFromB,"retained B ancestor"); cTable->addColumn<int>("isFromPromptD",cFromPromptD,"retained confirmed-prompt D ancestor"); cTable->addColumn<int>("isFromBtoD",cFromBtoD,"both B and D retained ancestors"); cTable->addColumn<int>("membershipType",cMembership,"0 unassigned, 1 HF decay, 2 retained fragmentation, 3 ambiguous"); cTable->addColumn<int>("isAmbiguous",cAmbiguous,"multiple retained systems or HF ancestors"); cTable->addColumn<float>("deltaRQuark",cDR,"quark deltaR or -999"); cTable->addColumn<float>("zTrackQuark",cZ,"quark projection or -999"); cTable->addColumn<float>("jTQuark",cJT,"transverse momentum to quark or -999");
  cTable->addColumn<int>("nCandidateQuarkSystems",cNCandidateSystems,"retained QFrag components encountered"); cTable->addColumn<int>("nCandidateHeavyHadrons",cNCandidateHadrons,"selected B/D hadron ancestors"); cTable->addColumn<int>("hasBHadronAncestor",cHasBAncestor,"retained B ancestor"); cTable->addColumn<int>("hasCHadronAncestor",cHasDAncestor,"retained D ancestor"); cTable->addColumn<int>("hasMultipleMotherBranches",cMultipleBranches,"an ancestor record has multiple mothers"); cTable->addColumn<int>("ancestryTruncated",cAncestryTruncated,"depth guard reached"); cTable->addColumn<int>("membershipReason",cMembershipReason,"1 unique selected HF hadron, 2 unique quark system outside HF decay, 3 multiple candidates, 4 no candidate"); cTable->addColumn<int>("lastPrunedPdgId",cLastPdg,"last pruned PDG id or -999"); cTable->addColumn<int>("lastPrunedStatus",cLastStatus,"last pruned status or -1"); cTable->addColumn<int>("lastPrunedIsSelf",cLastSelf,"1 when retained reference is status-1, the case PackedGenParticle::motherRef identifies as its self-copy");
  event.put(std::move(qTable), "QFragTable"); event.put(std::move(hTable), "QHadronTable"); event.put(std::move(cTable), "QGenChargedTable");
}

DEFINE_FWK_MODULE(QuarkFragmentationProducer);
