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

#include <cmath>
#include <array>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
  constexpr float kInvalid = -999.f;
  constexpr double kMinMomentumGeV = 1.e-6;
  constexpr double kMinTransverseFraction = 1.e-6;
  enum AssociationQuality { kMissing = 0, kUnique = 1, kAmbiguous = 2 };
  enum LineageQuality { kV2Missing = 0, kV2Unique = 1, kV2Ambiguous = 2, kV2Unknown = 3 };
  enum MembershipType { kUnassigned = 0, kHeavyHadronDecay = 1, kRetainedFragmentation = 2, kAmbiguousMembership = 3 };
  enum CharmOrigin { kCharmNotApplicable = 0, kPromptCharm = 1, kCharmFromB = 2, kCharmUnknown = 3 };
  enum HeavyContent { kContentUnsupported = 0, kContentCharm = 4, kContentAntiCharm = -4,
                      kContentBottom = 5, kContentAntiBottom = -5, kContentMultiHeavy = 6 };

  bool isQuark(int id) { return std::abs(id) >= 1 && std::abs(id) <= 5; }
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

  // Explicitly cover the open-heavy species selected by isBHadron/isDHadron.
  // Unsupported and multi-heavy species stay distinct; PDG sign alone is not
  // used because meson and baryon conventions differ.
  int heavyContent(int id) {
    const int a = std::abs(id), sign = id < 0 ? -1 : 1;
    if (a == 541 || a == 5142 || a == 5242 || a == 5342 || a == 5512 || a == 5532 || a == 5542 || a == 5554 || a == 4412 || a == 4422 || a == 4432 || a == 4444)
      return kContentMultiHeavy;
    if (a == 411 || a == 421 || a == 431) return sign > 0 ? kContentCharm : kContentAntiCharm;
    if (a == 4122 || a == 4132 || a == 4232 || a == 4332) return sign > 0 ? kContentCharm : kContentAntiCharm;
    if (a == 511 || a == 521 || a == 531) return sign > 0 ? kContentAntiBottom : kContentBottom;
    if (a == 5122 || a == 5132 || a == 5232 || a == 5332) return sign > 0 ? kContentBottom : kContentAntiBottom;
    return kContentUnsupported;
  }

  struct QuarkLineage {
    std::vector<int> copies;
  };

  struct CopyPath {
    std::vector<int> ordered;
    bool unique = false;
  };

  CopyPath analyzeCopyPath(const QuarkLineage& lineage, const std::unordered_set<int>& ancestryQuarks,
                           const reco::GenParticleCollection& particles,
                           const std::unordered_map<const reco::Candidate*, int>& lookup) {
    CopyPath result;
    std::unordered_set<int> nodes;
    for (const int copy : lineage.copies) if (ancestryQuarks.count(copy)) nodes.insert(copy);
    if (nodes.empty()) return result;
    std::unordered_map<int, std::vector<int>> daughters;
    std::unordered_map<int, int> nMothers;
    for (const int node : nodes) {
      const auto& particle = particles[node];
      std::unordered_set<int> uniqueMothers;
      for (size_t im = 0; im < particle.numberOfMothers(); ++im) {
        const auto found = lookup.find(particle.mother(im));
        if (particle.mother(im) && found != lookup.end() && nodes.count(found->second) &&
            particles[found->second].pdgId() == particle.pdgId() && uniqueMothers.insert(found->second).second) {
          ++nMothers[node];
          daughters[found->second].push_back(node);
        }
      }
    }
    int roots = 0, terminals = 0;
    bool branching = false, multiMother = false;
    for (const int node : nodes) {
      roots += nMothers[node] == 0;
      terminals += daughters[node].empty();
      branching |= daughters[node].size() > 1;
      multiMother |= nMothers[node] > 1;
    }
    if (roots != 1 || terminals != 1 || branching || multiMother) return result;
    int current = -1;
    for (const int node : nodes) if (nMothers[node] == 0) current = node;
    std::unordered_set<int> seen;
    while (current >= 0 && seen.insert(current).second) {
      result.ordered.push_back(current);
      if (daughters[current].empty()) break;
      current = daughters[current].front();
    }
    if (result.ordered.size() != nodes.size() || result.ordered.empty() || daughters[result.ordered.back()].size() != 0) {
      result.ordered.clear();
      return result;
    }
    result.unique = true;
    return result;
  }

  bool isHadronRecord(int pdgId) {
    const int id = std::abs(pdgId);
    return id >= 100 && id < 1000000000;
  }

  struct Ancestors {
    std::vector<int> indices;
    bool incomplete = false;  // null or non-pruned reference encountered
    bool truncated = false;   // traversal depth guard reached
  };
}

class QuarkFragmentationProducer : public edm::stream::EDProducer<> {
public:
  explicit QuarkFragmentationProducer(const edm::ParameterSet& cfg)
      : prunedToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
        packedToken_(consumes<pat::PackedGenParticleCollection>(cfg.getParameter<edm::InputTag>("packedGenParticles"))),
        saveDecayVertexTruth_(cfg.getUntrackedParameter<bool>("saveDecayVertexTruth", false)),
        saveFragmentationObservables_(cfg.getUntrackedParameter<bool>("saveFragmentationObservables", true)),
        minimalQuarkTruth_(cfg.getUntrackedParameter<bool>("minimalQuarkTruth", false)),
        saveDerivedDecayTruth_(cfg.getUntrackedParameter<bool>("saveDerivedDecayTruth", true)) {
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
  bool saveDecayVertexTruth_;
  bool saveFragmentationObservables_;
  bool minimalQuarkTruth_;
  bool saveDerivedDecayTruth_;
};

Ancestors QuarkFragmentationProducer::ancestors(
    int start, const reco::GenParticleCollection& particles,
    const std::unordered_map<const reco::Candidate*, int>& lookup) const {
  Ancestors result;
  if (start < 0 || static_cast<size_t>(start) >= particles.size()) { result.incomplete = true; return result; }
  constexpr unsigned int kTraversalDepthGuard = 256;
  std::queue<int> todo;
  std::unordered_set<int> seen;
  todo.push(start);
  seen.insert(start);
  while (!todo.empty()) {
    const int current = todo.front();
    todo.pop();
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
      if (path.count(found->second)) continue;
      detectCycle(found->second, depth + 1);
    }
    path.erase(current);
    completed.insert(current);
  };
  detectCycle(start, 0);
  return result;
}

namespace {
  struct QuarkFrontier {
    std::unordered_set<int> lineages;
    std::unordered_set<int> qfragSystems;
    bool incomplete = false;
  };

  // For additional charged particles, stop each retained mother branch at its
  // first quark-copy lineage. This removes upstream production quarks reached
  // only by continuing through an already encountered outgoing quark, while
  // retaining distinct quark co-mothers on separate branches.
  QuarkFrontier firstQuarkFrontier(
      int start,
      const reco::GenParticleCollection& particles,
      const std::unordered_map<const reco::Candidate*, int>& prunedIndex,
      const std::unordered_map<int, int>& prunedToLineage,
      const std::unordered_map<int, int>& prunedToQFrag) {
    QuarkFrontier result;
    if (start < 0 || static_cast<size_t>(start) >= particles.size()) {
      result.incomplete = true;
      return result;
    }
    std::queue<int> todo;
    std::unordered_set<int> seen;
    todo.push(start);
    seen.insert(start);
    while (!todo.empty()) {
      const int current = todo.front();
      todo.pop();
      const auto& particle = particles[current];
      if (isQuark(particle.pdgId())) {
        const auto lineage = prunedToLineage.find(current);
        if (lineage == prunedToLineage.end()) {
          result.incomplete = true;
        } else {
          result.lineages.insert(lineage->second);
          const auto qfrag = prunedToQFrag.find(current);
          if (qfrag != prunedToQFrag.end()) result.qfragSystems.insert(qfrag->second);
        }
        continue;
      }
      for (size_t im = 0; im < particle.numberOfMothers(); ++im) {
        const reco::Candidate* mother = particle.mother(im);
        if (!mother) {
          result.incomplete = true;
          continue;
        }
        const auto found = prunedIndex.find(mother);
        if (found == prunedIndex.end()) {
          result.incomplete = true;
          continue;
        }
        if (seen.insert(found->second).second) todo.push(found->second);
      }
    }
    return result;
  }
}  // namespace

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

  // Connected same-signed-PDG copies represent one showered quark system, so
  // retaining one terminal representative avoids counting copies as systems.
  // hardPrunedIdx is stored only when exactly one copy is flagged hard-process.
  std::vector<int> qRep, qHard;
  std::unordered_map<int, int> prunedToQFrag;
  std::unordered_map<int, int> prunedToLineage;
  std::vector<QuarkLineage> quarkLineages;
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
    QuarkLineage lineage;
    lineage.copies = component;
    const int lineageId = static_cast<int>(quarkLineages.size());
    quarkLineages.push_back(lineage);
    for (const int member : component) prunedToLineage[member] = lineageId;
    if (terminal.size() != 1) continue;
    const int row = static_cast<int>(qRep.size());
    const int representative = terminal.front();
    qRep.push_back(representative);
    qHard.push_back(hard.size() == 1 ? hard.front() : -1);
    for (const int member : component) prunedToQFrag[member] = row;
  }

  std::vector<int> hPruned;
  for (size_t i = 0; i < pruned->size(); ++i)
    if (isHeavyHadron((*pruned)[i].pdgId()) && isFinalCopy(i, *pruned, prunedIndex)) hPruned.push_back(i);
  std::unordered_map<int, int> prunedToHadron;
  for (size_t i = 0; i < hPruned.size(); ++i) prunedToHadron[hPruned[i]] = static_cast<int>(i);

  std::vector<int> hadFromB(hPruned.size(), 0), hadParentB(hPruned.size(), -1);
  std::vector<int> hadCharmOrigin(hPruned.size(), kCharmNotApplicable);
  std::vector<int> v2QFrag(hPruned.size(), -1), v2Quality(hPruned.size(), kV2Missing), v2SourceCertain(hPruned.size(), 0);
  std::vector<int> v2PreHad(hPruned.size(), -1), v2KinValid(hPruned.size(), 0);
  std::vector<float> v2PreHadPt(hPruned.size(), kInvalid), v2PreHadEta(hPruned.size(), kInvalid), v2PreHadPhi(hPruned.size(), kInvalid);
  std::vector<float> hZProj(hPruned.size(), kInvalid), hZT(hPruned.size(), kInvalid),
      hDeltaRQuark(hPruned.size(), kInvalid), hKTrel(hPruned.size(), kInvalid);
  std::vector<int> hObservableInputsValid(hPruned.size(), 0), hAuxiliaryTargetValid(hPruned.size(), 0),
      hBFragmentationTargetValid(hPruned.size(), 0);
  for (size_t ih = 0; ih < hPruned.size(); ++ih) {
    const auto ancestry = ancestors(hPruned[ih], *pruned, prunedIndex);
    std::unordered_set<int> bAncestors;
    for (const int ancestor : ancestry.indices) {
      if (!isBHadron((*pruned)[ancestor].pdgId())) continue;
      hadFromB[ih] = 1;
      const auto b = prunedToHadron.find(ancestor);
      if (b != prunedToHadron.end()) bAncestors.insert(b->second);
    }
    if (bAncestors.size() == 1) hadParentB[ih] = *bAncestors.begin();
    const int hadronId = (*pruned)[hPruned[ih]].pdgId();
    const bool isD = isDHadron(hadronId);
    if (isD) {
      // Missing retained B ancestry is not evidence of prompt charm; prompt
      // origin requires MiniAOD's explicit status flag.
      const bool explicitPrompt = (*pruned)[hPruned[ih]].statusFlags().isPrompt();
      hadCharmOrigin[ih] = hadFromB[ih] ? kCharmFromB : (explicitPrompt ? kPromptCharm : kCharmUnknown);
    }

    // Match the compatible heavy-quark lineage while retaining unknown D origin.
    int v2SourcePruned = -1;
    int v2Content = kContentUnsupported;
    bool v2CanSearch = true;
    if (isD) {
      if (hadFromB[ih]) {
        if (bAncestors.size() == 1) {
          // A B-to-D hadronization source is the parent bottom system.
          const int bRow = *bAncestors.begin();
          v2SourceCertain[ih] = 1;
          v2SourcePruned = hPruned[bRow];
        } else {
          v2SourceCertain[ih] = 0;
          v2CanSearch = false;  // positive B evidence must not fall back to charm
          v2Quality[ih] = kV2Unknown;
        }
      } else {
        v2SourcePruned = hPruned[ih];
        if (hadCharmOrigin[ih] == kPromptCharm) {
          v2SourceCertain[ih] = 1;
        } else {
          v2SourceCertain[ih] = 0;
        }
      }
    } else if (isBHadron(hadronId)) {
      v2SourceCertain[ih] = 1;
      v2SourcePruned = hPruned[ih];
    }
    if (v2CanSearch && v2SourcePruned >= 0) {
      v2Content = heavyContent((*pruned)[v2SourcePruned].pdgId());
      if (v2Content == kContentUnsupported || v2Content == kContentMultiHeavy) {
        v2Quality[ih] = kV2Unknown;
      } else {
        const Ancestors v2Ancestry = ancestors(v2SourcePruned, *pruned, prunedIndex);
        std::unordered_set<int> ancestryQuarks, candidateLineages, compatibleLineages;
        for (const int ancestor : v2Ancestry.indices) {
          if (!isQuark((*pruned)[ancestor].pdgId())) continue;
          ancestryQuarks.insert(ancestor);
          const auto lineageIt = prunedToLineage.find(ancestor);
          if (lineageIt == prunedToLineage.end()) continue;
          candidateLineages.insert(lineageIt->second);
          // Keep only the signed heavy-quark flavour compatible with the
          // selected B or D hadron's valence heavy flavour.
          if ((*pruned)[ancestor].pdgId() == v2Content) compatibleLineages.insert(lineageIt->second);
        }
        std::vector<std::pair<int, CopyPath>> compatiblePaths;
        for (const int lineageId : compatibleLineages) {
          auto path = analyzeCopyPath(quarkLineages[lineageId], ancestryQuarks, *pruned, prunedIndex);
          compatiblePaths.emplace_back(lineageId, std::move(path));
        }
        if (compatibleLineages.size() > 1) {
          v2Quality[ih] = kV2Ambiguous;
        } else if (compatibleLineages.empty()) {
          const bool incomplete = v2Ancestry.incomplete || v2Ancestry.truncated;
          v2Quality[ih] = incomplete ? kV2Unknown : kV2Missing;
        } else {
          const CopyPath& path = compatiblePaths.front().second;
          if (!path.unique) {
            v2Quality[ih] = kV2Ambiguous;
          } else if (v2Ancestry.incomplete || v2Ancestry.truncated) {
            v2Quality[ih] = kV2Unknown;
          } else {
            v2Quality[ih] = kV2Unique;
          }
          if (!path.ordered.empty()) {
            std::unordered_set<int> sourceAncestryNodes(v2Ancestry.indices.begin(), v2Ancestry.indices.end());
            sourceAncestryNodes.insert(v2SourcePruned);
            std::vector<int> hadronAdjacent;
            for (const int copy : path.ordered) {
              bool hasHadronDaughter = false;
              for (const int descendant : sourceAncestryNodes) {
                if (!isHadronRecord((*pruned)[descendant].pdgId())) continue;
                const auto& hadron = (*pruned)[descendant];
                for (size_t im = 0; im < hadron.numberOfMothers(); ++im) {
                  const auto mother = prunedIndex.find(hadron.mother(im));
                  if (hadron.mother(im) && mother != prunedIndex.end() && mother->second == copy) hasHadronDaughter = true;
                }
              }
              if (hasHadronDaughter) hadronAdjacent.push_back(copy);
            }
            if (hadronAdjacent.size() == 1) {
              v2PreHad[ih] = hadronAdjacent.front();
            } else if (hadronAdjacent.size() > 1) {
              // The unique directed copy path orders these copies. Select its
              // last hadron-adjacent member, not an arbitrary collection row.
              v2PreHad[ih] = hadronAdjacent.back();
            }
          }
          if (v2Quality[ih] == kV2Unique && !path.ordered.empty()) {
            const auto row = prunedToQFrag.find(path.ordered.back());
            if (row != prunedToQFrag.end()) v2QFrag[ih] = row->second;
          }
          if (v2PreHad[ih] >= 0) {
            const auto& q = (*pruned)[v2PreHad[ih]];
            const auto& source = (*pruned)[v2SourcePruned];
            const double px = q.px(), py = q.py(), pz = q.pz();
            const double p2 = px*px + py*py + pz*pz;
            const double p = std::sqrt(p2), pt = q.pt();
            const double sourceEta = source.eta(), sourcePhi = source.phi();
            if (std::isfinite(px) && std::isfinite(py) && std::isfinite(pz) && std::isfinite(p2) && std::isfinite(pt) &&
                std::isfinite(source.px()) && std::isfinite(source.py()) && std::isfinite(source.pz()) &&
                p > kMinMomentumGeV && pt / p > kMinTransverseFraction &&
                std::isfinite(q.eta()) && std::isfinite(q.phi()) && std::isfinite(sourceEta) && std::isfinite(sourcePhi)) {
              const double z = (source.px()*px + source.py()*py + source.pz()*pz) / p2;
              const double dr = reco::deltaR(sourceEta, sourcePhi, q.eta(), q.phi());
              if (std::isfinite(z) && std::isfinite(dr)) {
                v2KinValid[ih] = 1;
                v2PreHadPt[ih] = q.pt(); v2PreHadEta[ih] = q.eta(); v2PreHadPhi[ih] = q.phi();
              }
            }
          }

        }
      }
    }
  }

  // Persist the same unclipped momentum observables used by the Stage 3
  // prompt-D study. The values are available for uniquely associated hadrons
  // with valid proxy kinematics; the two target masks below separately select
  // confirmed prompt-D and B fragmentation targets.
  for (size_t ih = 0; ih < hPruned.size(); ++ih) {
    const int qidx = v2QFrag[ih];
    if (v2Quality[ih] != kV2Unique || qidx < 0 || static_cast<size_t>(qidx) >= qRep.size() || !v2KinValid[ih]) continue;
    const auto& hadron = (*pruned)[hPruned[ih]];
    const auto& quark = (*pruned)[v2PreHad[ih]];
    const double q2 = quark.px() * quark.px() + quark.py() * quark.py() + quark.pz() * quark.pz();
    const double qnorm = std::sqrt(q2);
    const double qpt = quark.pt();
    const double hpt = hadron.pt();
    if (!(std::isfinite(hadron.px()) && std::isfinite(hadron.py()) && std::isfinite(hadron.pz()) &&
          std::isfinite(hadron.eta()) && std::isfinite(hadron.phi()) && std::isfinite(hpt) && hpt > 0. &&
          std::isfinite(quark.px()) && std::isfinite(quark.py()) && std::isfinite(quark.pz()) &&
          std::isfinite(quark.eta()) && std::isfinite(quark.phi()) && std::isfinite(q2) && q2 > 0. &&
          std::isfinite(qnorm) && qnorm > 0. && std::isfinite(qpt) && qpt > 0.)) continue;
    const double dot = hadron.px() * quark.px() + hadron.py() * quark.py() + hadron.pz() * quark.pz();
    const double crossX = hadron.py() * quark.pz() - hadron.pz() * quark.py();
    const double crossY = hadron.pz() * quark.px() - hadron.px() * quark.pz();
    const double crossZ = hadron.px() * quark.py() - hadron.py() * quark.px();
    const double crossNorm = std::sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
    const double zProj = dot / q2;
    const double zT = hadron.pt() / qpt;
    const double deltaR = reco::deltaR(hadron.eta(), hadron.phi(), quark.eta(), quark.phi());
    const double kTrel = crossNorm / qnorm;
    if (!(std::isfinite(zProj) && std::isfinite(zT) && std::isfinite(deltaR) && std::isfinite(kTrel))) continue;
    hObservableInputsValid[ih] = 1;
    hZProj[ih] = zProj;
    hZT[ih] = zT;
    hDeltaRQuark[ih] = deltaR;
    hKTrel[ih] = kTrel;

    const int quarkPdgId = (*pruned)[qRep[qidx]].pdgId();
    const bool validPromptD = hadCharmOrigin[ih] == kPromptCharm && std::abs(quarkPdgId) == 4 &&
                              v2SourceCertain[ih] == 1;
    const bool validB = isBHadron(hadron.pdgId()) && std::abs(quarkPdgId) == 5 && v2SourceCertain[ih] == 1;
    hAuxiliaryTargetValid[ih] = validPromptD;
    hBFragmentationTargetValid[ih] = validB;
  }

  std::vector<int> cPacked, cPdg, cCharge, cLastPruned, cLastSelf, cFromB, cFromPromptD, cFromBtoD;
  std::vector<int> cV2Hadron, cV2Quark, cV2Quality, cV2Membership, cV2SourceCertain;
  std::vector<float> cPt, cEta, cPhi, cMass;
  for (size_t ip = 0; ip < packed->size(); ++ip) {
    const auto& p = (*packed)[ip];
    if (p.status() != 1 || p.charge() == 0) continue;
    int last = -1;
    const auto& ref = p.lastPrunedRef();
    if (ref.isNonnull() && ref.isAvailable() && ref.id() == pruned.id()) last = static_cast<int>(ref.key());
    Ancestors ancestry;
    if (last >= 0) ancestry = ancestors(last, *pruned, prunedIndex); else ancestry.incomplete = true;
    bool hasB = false, hasD = false, hasPromptD = false;
    for (const int ancestor : ancestry.indices) {
      if (isBHadron((*pruned)[ancestor].pdgId())) hasB = true;
      if (!isDHadron((*pruned)[ancestor].pdgId())) continue;
      hasD = true;
      const auto d = prunedToHadron.find(ancestor);
      if (d != prunedToHadron.end() && hadCharmOrigin[d->second] == kPromptCharm) hasPromptD = true;
    }
    const bool hasBtoD = hasB && hasD;
    // PackedGenParticle::motherRef treats a stored status-1 pruned reference
    // as the packed particle's own pruned representation and steps past it.
    const bool lastSelf = last >= 0 && (*pruned)[last].status() == 1;

    // V2 charged assignment: distances are measured on the retained mother
    // graph from lastPrunedRef (which may be the packed particle's own
    // status-1 pruned record). The nearest final B/D ancestor owns the decay
    // association. Equal-distance distinct hadrons stay ambiguous.
    std::unordered_map<int, unsigned int> ancestorDistance;
    std::queue<std::pair<int, unsigned int>> distanceQueue;
    if (last >= 0) { ancestorDistance[last] = 0; distanceQueue.push({last, 0}); }
    while (!distanceQueue.empty()) {
      const int node = distanceQueue.front().first;
      const unsigned int depth = distanceQueue.front().second;
      distanceQueue.pop();
      const auto& nodeParticle = (*pruned)[node];
      for (size_t im = 0; im < nodeParticle.numberOfMothers(); ++im) {
        const auto found = prunedIndex.find(nodeParticle.mother(im));
        if (nodeParticle.mother(im) && found != prunedIndex.end() && !ancestorDistance.count(found->second)) {
          ancestorDistance[found->second] = depth + 1;
          distanceQueue.push({found->second, depth + 1});
        }
      }
    }
    unsigned int nearestHeavyDistance = std::numeric_limits<unsigned int>::max();
    std::unordered_set<int> nearestHeavyPruned;
    for (const auto& entry : ancestorDistance) {
      const int index = entry.first;
      if (!isHeavyHadron((*pruned)[index].pdgId())) continue;
      if (entry.second < nearestHeavyDistance) {
        nearestHeavyDistance = entry.second;
        nearestHeavyPruned.clear();
      }
      if (entry.second == nearestHeavyDistance) {
        nearestHeavyPruned.insert(index);
      }
    }
    int v2Hadron = -1;
    if (nearestHeavyPruned.size() == 1) {
      const auto found = prunedToHadron.find(*nearestHeavyPruned.begin());
      if (found != prunedToHadron.end()) v2Hadron = found->second;
    }
    int v2Quark = -1, v2AssociationQuality = kMissing;
    int v2Membership = kUnassigned, v2SourceIsCertain = 0;
    std::unordered_set<int> candidateLineages, candidateV2Systems;
    bool candidateAncestryIncomplete = ancestry.incomplete || ancestry.truncated;
    if (nearestHeavyPruned.empty()) {
      const QuarkFrontier frontier = firstQuarkFrontier(last, *pruned, prunedIndex, prunedToLineage, prunedToQFrag);
      candidateLineages = frontier.lineages;
      candidateV2Systems = frontier.qfragSystems;
      candidateAncestryIncomplete = candidateAncestryIncomplete || frontier.incomplete;
    } else {
      // Heavy-hadron decay rows use the validated QHadron association.
      std::unordered_set<int> ancestryQuarks(ancestry.indices.begin(), ancestry.indices.end());
      if (last >= 0 && isQuark((*pruned)[last].pdgId())) ancestryQuarks.insert(last);
      for (const int ancestor : ancestryQuarks) {
        const auto lineage = prunedToLineage.find(ancestor);
        if (lineage == prunedToLineage.end()) continue;
        candidateLineages.insert(lineage->second);
        const auto row = prunedToQFrag.find(ancestor);
        if (row != prunedToQFrag.end()) candidateV2Systems.insert(row->second);
      }
    }
    if (v2Hadron >= 0 && nearestHeavyPruned.size() == 1) {
      v2Membership = kHeavyHadronDecay;
      v2AssociationQuality = v2Quality[v2Hadron];
      v2Quark = v2QFrag[v2Hadron];
      v2SourceIsCertain = v2SourceCertain[v2Hadron];
    } else if (nearestHeavyPruned.size() > 1) {
      v2Membership = kAmbiguousMembership;
      v2AssociationQuality = kAmbiguous;
    } else if (nearestHeavyPruned.size() == 1) {
      // A retained heavy-hadron ancestor exists but is not represented by a
      // selected QHadron row. Keep decay membership while withholding indices.
      v2Membership = kHeavyHadronDecay;
      v2AssociationQuality = kV2Unknown;
    } else if (nearestHeavyPruned.empty()) {
      if (candidateLineages.size() > 1 || candidateV2Systems.size() > 1) {
        v2AssociationQuality = kAmbiguous;
        v2Membership = kAmbiguousMembership;
      } else if (candidateV2Systems.size() == 1 && candidateLineages.size() == 1 && !candidateAncestryIncomplete) {
        v2Quark = *candidateV2Systems.begin();
        v2AssociationQuality = kUnique;
        v2Membership = kRetainedFragmentation;
        // A unique retained ancestry does not establish the physical source;
        // additional-particle certainty remains false.
      } else if (candidateLineages.size() == 1 || candidateAncestryIncomplete || last < 0) {
        v2AssociationQuality = kV2Unknown;
      }
    }
    cPacked.push_back(ip); cPdg.push_back(p.pdgId()); cCharge.push_back(p.charge());
    cPt.push_back(p.pt()); cEta.push_back(p.eta()); cPhi.push_back(p.phi()); cMass.push_back(p.mass());
    cLastPruned.push_back(last); cLastSelf.push_back(lastSelf);
    cFromB.push_back(hasB); cFromPromptD.push_back(hasPromptD); cFromBtoD.push_back(hasBtoD);
    cV2Hadron.push_back(v2Hadron); cV2Quark.push_back(v2Quark); cV2Quality.push_back(v2AssociationQuality);
    cV2Membership.push_back(v2Membership); cV2SourceCertain.push_back(v2SourceIsCertain);
  }

  auto qTable = std::make_unique<nanoaod::FlatTable>(qRep.size(), "QFrag", false, false);
  std::vector<int> qPdg, qOriginCategory, qOriginValid;
  std::vector<float> qPt, qEta, qPhi, qMass, qPx, qPy, qPz, qEnergy, qHardPt, qHardEta, qHardPhi;
  for (size_t iq = 0; iq < qRep.size(); ++iq) {
    const auto& q = (*pruned)[qRep[iq]];
    qPdg.push_back(q.pdgId()); qPt.push_back(q.pt()); qEta.push_back(q.eta()); qPhi.push_back(q.phi()); qMass.push_back(q.mass());
    qPx.push_back(q.px()); qPy.push_back(q.py()); qPz.push_back(q.pz()); qEnergy.push_back(q.energy());
    qOriginCategory.push_back(qHard[iq] >= 0 ? 1 : 0); qOriginValid.push_back(qHard[iq] >= 0 ? 1 : 0);
    if (qHard[iq] >= 0) {
      const auto& hard = (*pruned)[qHard[iq]];
      qHardPt.push_back(hard.pt()); qHardEta.push_back(hard.eta()); qHardPhi.push_back(hard.phi());
    } else {
      qHardPt.push_back(kInvalid); qHardEta.push_back(kInvalid); qHardPhi.push_back(kInvalid);
    }
  }
  qTable->addColumn<int>("pdgId", qPdg, "signed PDG identifier");
  qTable->addColumn<float>("px", qPx, "representative quark px [GeV]");
  qTable->addColumn<float>("py", qPy, "representative quark py [GeV]");
  qTable->addColumn<float>("pz", qPz, "representative quark pz [GeV]");
  qTable->addColumn<float>("energy", qEnergy, "representative quark energy [GeV]");
  qTable->addColumn<int>("originCategory", qOriginCategory, "0 unavailable/ambiguous hard-copy ancestry, 1 unique hard-process copy retained");
  qTable->addColumn<int>("originCategoryValid", qOriginValid, "1 when exactly one hard-process copy is retained in the connected quark-copy component");
  if (!minimalQuarkTruth_) {
    qTable->addColumn<int>("prunedGenIdx", qRep, "representative terminal index in prunedGenParticles");
    qTable->addColumn<float>("pt", qPt, "representative transverse momentum [GeV]");
    qTable->addColumn<float>("eta", qEta, "representative pseudorapidity");
    qTable->addColumn<float>("phi", qPhi, "representative azimuth [rad]");
    qTable->addColumn<float>("mass", qMass, "representative mass [GeV]");
    qTable->addColumn<int>("hardPrunedIdx", qHard, "unique hard-process copy index in prunedGenParticles, or -1");
    qTable->addColumn<float>("hardPt", qHardPt, "hard-copy transverse momentum [GeV], or -999");
    qTable->addColumn<float>("hardEta", qHardEta, "hard-copy pseudorapidity, or -999");
    qTable->addColumn<float>("hardPhi", qHardPhi, "hard-copy azimuth [rad], or -999");
  }

  auto hTable = std::make_unique<nanoaod::FlatTable>(hPruned.size(), "QHadron", false, false);
  std::vector<int> hIndex,hPdg,hStatus,hIsPrompt,hIsHardProcess,hStatusLastCopy;
  std::vector<float> hPt, hEta, hPhi, hMass,hP,hPx,hPy,hPz,hEnergy;
  std::vector<float> prodX, prodY, prodZ, decayX, decayY, decayZ, trueLxy, trueL3D, trueCt;
  std::vector<int> prodValid, decayValid, decaySource, decayResolved, trueLengthValid, trueCtValid;
  for (const int index : hPruned) {
    const auto& h = (*pruned)[index];
    hPdg.push_back(h.pdgId()); hPt.push_back(h.pt()); hEta.push_back(h.eta()); hPhi.push_back(h.phi()); hMass.push_back(h.mass());
    if (!saveDecayVertexTruth_) continue;
    hIndex.push_back(hIndex.size());
    hStatus.push_back(h.status());hIsPrompt.push_back(h.statusFlags().isPrompt());hIsHardProcess.push_back(h.statusFlags().isHardProcess());hStatusLastCopy.push_back(h.statusFlags().isLastCopy());
    hP.push_back(h.p());hPx.push_back(h.px());hPy.push_back(h.py());hPz.push_back(h.pz());hEnergy.push_back(h.energy());
    const double x = h.vx(), y = h.vy(), z = h.vz();
    const bool validProd = std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    prodX.push_back(validProd ? x : kInvalid); prodY.push_back(validProd ? y : kInvalid); prodZ.push_back(validProd ? z : kInvalid);
    prodValid.push_back(validProd ? 1 : 0);
    std::vector<std::array<double, 3>> daughterVertices;
    for (size_t id = 0; id < h.numberOfDaughters(); ++id) {
      const reco::Candidate* daughter = h.daughter(id);
      if (!daughter) continue;
      const double dx = daughter->vx(), dy = daughter->vy(), dz = daughter->vz();
      if (std::isfinite(dx) && std::isfinite(dy) && std::isfinite(dz))
        daughterVertices.push_back(std::array<double, 3>{{dx, dy, dz}});
    }
    bool resolved = !daughterVertices.empty();
    if (resolved) {
      constexpr double kVertexToleranceCm = 1.e-4;
      for (const auto& v : daughterVertices) {
        const double distance = std::sqrt(std::pow(v[0] - daughterVertices.front()[0], 2) +
                                          std::pow(v[1] - daughterVertices.front()[1], 2) +
                                          std::pow(v[2] - daughterVertices.front()[2], 2));
        if (!std::isfinite(distance) || distance > kVertexToleranceCm) resolved = false;
      }
    }
    const bool validDecay = resolved;
    const double dx = validDecay ? daughterVertices.front()[0] : 0.;
    const double dy = validDecay ? daughterVertices.front()[1] : 0.;
    const double dz = validDecay ? daughterVertices.front()[2] : 0.;
    decayX.push_back(validDecay ? dx : kInvalid); decayY.push_back(validDecay ? dy : kInvalid); decayZ.push_back(validDecay ? dz : kInvalid);
    decayValid.push_back(validDecay ? 1 : 0); decayResolved.push_back(resolved ? 1 : 0);
    decaySource.push_back(validDecay ? 1 : 0);  // immediate daughter production vertex
    const double lx = validDecay && validProd ? dx - x : 0.;
    const double ly = validDecay && validProd ? dy - y : 0.;
    const double lz = validDecay && validProd ? dz - z : 0.;
    const bool lengthValid = validDecay && validProd && std::isfinite(lx + ly + lz);
    const double lxy = lengthValid ? std::hypot(lx, ly) : 0.;
    const double l3d = lengthValid ? std::sqrt(lx * lx + ly * ly + lz * lz) : 0.;
    trueLxy.push_back(lengthValid ? lxy : kInvalid); trueL3D.push_back(lengthValid ? l3d : kInvalid);
    trueLengthValid.push_back(lengthValid ? 1 : 0);
    const double p2 = h.px() * h.px() + h.py() * h.py() + h.pz() * h.pz();
    const double ct = lengthValid && p2 > 0. && std::isfinite(p2) ? h.mass() * (lx * h.px() + ly * h.py() + lz * h.pz()) / p2 : 0.;
    const bool ctValid = lengthValid && std::isfinite(ct);
    trueCt.push_back(ctValid ? ct : kInvalid); trueCtValid.push_back(ctValid ? 1 : 0);
  }
  hTable->addColumn<int>("prunedGenIdx", hPruned, "index in prunedGenParticles");
  hTable->addColumn<int>("pdgId", hPdg, "signed PDG identifier");
  hTable->addColumn<float>("pt", hPt, "transverse momentum [GeV]");
  hTable->addColumn<float>("eta", hEta, "pseudorapidity");
  hTable->addColumn<float>("phi", hPhi, "azimuth [rad]");
  hTable->addColumn<float>("mass", hMass, "invariant mass [GeV]");
  hTable->addColumn<int>("charmOrigin", hadCharmOrigin, "0 not D, 1 confirmed prompt D, 2 B-to-D, 3 unknown D origin");
  hTable->addColumn<int>("parentBHadronIdx", hadParentB, "event-local QHadron row for unique retained B ancestor, or -1");
  hTable->addColumn<int>("quarkIdx", v2QFrag, "validated QFrag row, or -1 when unresolved/unmapped");
  hTable->addColumn<int>("associationQuality", v2Quality, "0 missing, 1 unique ancestry, 2 ambiguous, 3 unknown/incomplete");
  hTable->addColumn<int>("fragmentationSourceCertain", v2SourceCertain, "1 only when a physical fragmentation source is established");
  if (saveFragmentationObservables_) {
    hTable->addColumn<int>("preHadPrunedIdx", v2PreHad, "validated pre-hadronization copy index in prunedGenParticles, or -1");
    hTable->addColumn<float>("preHadPt", v2PreHadPt, "pre-hadronization proxy transverse momentum [GeV], or -999");
    hTable->addColumn<float>("preHadEta", v2PreHadEta, "pre-hadronization proxy pseudorapidity, or -999");
    hTable->addColumn<float>("preHadPhi", v2PreHadPhi, "pre-hadronization proxy azimuth [rad], or -999");
    hTable->addColumn<int>("preHadKinematicsValid", v2KinValid, "1 valid pre-hadronization kinematics, 0 otherwise");
    hTable->addColumn<float>("zProj", hZProj, "dot(p_hadron,p_quark)/|p_quark|^2; unclipped, or -999 when unavailable");
    hTable->addColumn<float>("zT", hZT, "hadron pT / associated quark-proxy pT; unclipped, or -999 when unavailable");
    hTable->addColumn<float>("deltaR_DQuark", hDeltaRQuark, "deltaR between hadron and associated quark proxy (D or B), or -999");
    hTable->addColumn<float>("kTrel", hKTrel, "|p_hadron x p_quark|/|p_quark| [GeV], or -999 when unavailable");
    hTable->addColumn<int>("observableInputsValid", hObservableInputsValid, "1 for unique valid QFrag association, valid proxy kinematics, and finite nonzero observable inputs");
    hTable->addColumn<int>("auxiliaryTargetValid", hAuxiliaryTargetValid, "prompt-D fragmentation target mask; preserves prompt-D meaning");
    hTable->addColumn<int>("bFragmentationTargetValid", hBFragmentationTargetValid, "B fragmentation target mask; 1 only for uniquely associated, source-certain b-quark proxy and valid inputs");
  }
  if (saveDecayVertexTruth_) {
    hTable->addColumn<int>("hadronIdx",hIndex,"event-local QHadron row identifier");
    hTable->addColumn<int>("status",hStatus,"pruned generator status of selected final-copy hadron");
    hTable->addColumn<int>("isPrompt",hIsPrompt,"MiniAOD GenStatusFlags isPrompt bit");
    hTable->addColumn<int>("isHardProcess",hIsHardProcess,"MiniAOD GenStatusFlags isHardProcess bit");
    hTable->addColumn<int>("statusFlagLastCopy",hStatusLastCopy,"MiniAOD GenStatusFlags isLastCopy bit; selected final copy also follows the producer's daughter-copy test");
    hTable->addColumn<float>("p",hP,"total generator momentum [GeV]");hTable->addColumn<float>("px",hPx,"generator px [GeV]");hTable->addColumn<float>("py",hPy,"generator py [GeV]");hTable->addColumn<float>("pz",hPz,"generator pz [GeV]");hTable->addColumn<float>("energy",hEnergy,"generator energy [GeV]");
    hTable->addColumn<float>("prodX", prodX, "generator hadron production x [cm]");
    hTable->addColumn<float>("prodY", prodY, "generator hadron production y [cm]");
    hTable->addColumn<float>("prodZ", prodZ, "generator hadron production z [cm]");
    hTable->addColumn<int>("prodVertexValid", prodValid, "1 when generator hadron production position is finite");
    hTable->addColumn<float>("decayX", decayX, "generator decay x from consistent immediate-daughter production vertices [cm], or -999");
    hTable->addColumn<float>("decayY", decayY, "generator decay y from consistent immediate-daughter production vertices [cm], or -999");
    hTable->addColumn<float>("decayZ", decayZ, "generator decay z from consistent immediate-daughter production vertices [cm], or -999");
    hTable->addColumn<int>("decayVertexValid", decayValid, "1 only when immediate-daughter vertices are present and mutually consistent");
    hTable->addColumn<int>("decayVertexResolved", decayResolved, "1 when the decay vertex was unambiguous within 1e-4 cm");
    hTable->addColumn<int>("decayVertexSource", decaySource, "1 immediate daughter production vertex; 0 unavailable/ambiguous");
    if (saveDerivedDecayTruth_) {
      hTable->addColumn<float>("trueLxy", trueLxy, "transverse generator decay length [cm], or -999");
      hTable->addColumn<float>("trueL3D", trueL3D, "3D generator decay length [cm], or -999");
      hTable->addColumn<int>("trueDecayLengthValid", trueLengthValid, "1 when production and decay positions are valid");
      hTable->addColumn<float>("trueCt", trueCt, "projected proper decay length m*(L dot p)/|p|^2 [cm], or -999");
      hTable->addColumn<int>("trueCtValid", trueCtValid, "1 when true proper decay length inputs and result are valid");
    }
  }

  auto cTable = std::make_unique<nanoaod::FlatTable>(cPacked.size(), "QGenCharged", false, false);
  cTable->addColumn<int>("packedGenIdx", cPacked, "index in packedGenParticles");
  cTable->addColumn<int>("pdgId", cPdg, "signed PDG identifier");
  cTable->addColumn<int>("charge", cCharge, "electric charge in units of e");
  cTable->addColumn<float>("pt", cPt, "transverse momentum [GeV]");
  cTable->addColumn<float>("eta", cEta, "pseudorapidity");
  cTable->addColumn<float>("phi", cPhi, "azimuth [rad]");
  cTable->addColumn<float>("mass", cMass, "invariant mass [GeV]");
  cTable->addColumn<int>("lastPrunedIdx", cLastPruned, "lastPrunedRef key, or -1 if unavailable");
  cTable->addColumn<int>("lastPrunedIsSelf", cLastSelf, "1 when lastPrunedRef is the packed particle's status-1 self record");
  cTable->addColumn<int>("hadronIdx", cV2Hadron, "nearest uniquely selected heavy ancestor, event-local QHadron row; or -1");
  cTable->addColumn<int>("quarkIdx", cV2Quark, "validated QFrag row; or -1");
  cTable->addColumn<int>("associationQuality", cV2Quality, "0 missing, 1 unique ancestry, 2 ambiguous, 3 unknown/incomplete");
  cTable->addColumn<int>("membershipType", cV2Membership, "0 unassigned, 1 heavy-hadron decay, 2 additional candidate, 3 ambiguous");
  cTable->addColumn<int>("fragmentationSourceCertain", cV2SourceCertain, "physical source certainty; additional candidates remain 0");
  cTable->addColumn<int>("isFromB", cFromB, "retained B-hadron ancestry");
  cTable->addColumn<int>("isFromPromptD", cFromPromptD, "retained confirmed-prompt D ancestry");
  cTable->addColumn<int>("isFromBtoD", cFromBtoD, "both retained B and D ancestry");
  event.put(std::move(qTable), "QFragTable");
  event.put(std::move(hTable), "QHadronTable");
  event.put(std::move(cTable), "QGenChargedTable");
}

DEFINE_FWK_MODULE(QuarkFragmentationProducer);
