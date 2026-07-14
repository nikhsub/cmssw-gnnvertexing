#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>
#include <optional>
#include <algorithm>

typedef reco::Vertex::CovarianceMatrix CovMatrix;

class GenVertexProducer : public edm::stream::EDProducer<> {
public:
  explicit GenVertexProducer(const edm::ParameterSet&);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  int checkPDG(int abs_pdg) const;
  bool hasBHadronAncestor(const reco::Candidate* cand) const;

  int getDaughterOriginLabelNoPU(const reco::Candidate* dau) const;

  std::optional<std::tuple<float, float, float>> isAncestor(
      const reco::Candidate* mother,
      const reco::Candidate* daughter) const;

  std::vector<std::vector<float>> computeDistanceMatrix(
      const std::vector<float>& SV_x,
      const std::vector<float>& SV_y,
      const std::vector<float>& SV_z,
      std::vector<CovMatrix> SV_cov,
      const std::vector<float>& Hadron_GVx,
      const std::vector<float>& Hadron_GVy,
      const std::vector<float>& Hadron_GVz);

  std::tuple<std::vector<int>, std::vector<float>, std::vector<float>> matchHadronsToSV(
      std::vector<std::vector<float>> distances,
      const std::vector<float>& SVtrk_pt,
      const std::vector<float>& SVtrk_eta,
      const std::vector<float>& SVtrk_phi,
      const std::vector<int>& SVtrk_SVidx,
      const std::vector<float>& Daughters_pt,
      const std::vector<float>& Daughters_eta,
      const std::vector<float>& Daughters_phi,
      const std::vector<int>& Daughters_GVidx,
      int n_Hadrons,
      int nRequiredCommonTracks,
      double dR_max,
      double relPt_max,
      bool doubleMatching,
      double doubleMatching_maxSignificance,
      double doubleMatching_dR_max,
      double doubleMatching_relPt_max);

  const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;

  int nRequiredCommonTracks_;
  double dlenSigMin_;
  double dR_max_;
  double relPt_max_;
  bool doubleMatching_;
  double doubleMatching_maxSignificance_;
  double doubleMatching_dR_max_;
  double doubleMatching_relPt_max_;
};

GenVertexProducer::GenVertexProducer(const edm::ParameterSet& iConfig)
    : pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      genToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
      svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
      nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
      dR_max_(iConfig.getParameter<double>("dR_max")),
      relPt_max_(iConfig.getParameter<double>("relPt_max")),
      doubleMatching_(iConfig.getParameter<bool>("doubleMatching")),
      doubleMatching_maxSignificance_(iConfig.getParameter<double>("doubleMatching_maxSignificance")),
      doubleMatching_dR_max_(iConfig.getParameter<double>("doubleMatching_dR_max")),
      doubleMatching_relPt_max_(iConfig.getParameter<double>("doubleMatching_relPt_max")) {
  produces<nanoaod::FlatTable>("GVTable");
  produces<nanoaod::FlatTable>("GVDaughtersTable");
}

void GenVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::View<reco::Candidate>> genHandle;
  iEvent.getByToken(genToken_, genHandle);

  edm::Handle<std::vector<reco::Vertex>> svHandle;
  iEvent.getByToken(svToken_, svHandle);

  auto pvsIn = iEvent.getHandle(pvs_);

  if (!genHandle.isValid() || !svHandle.isValid() || !pvsIn.isValid() || pvsIn->empty()) {
    auto gvTable = std::make_unique<nanoaod::FlatTable>(0, "GV", false);
    auto dauTable = std::make_unique<nanoaod::FlatTable>(0, "GVDaughters", false);
    iEvent.put(std::move(gvTable), "GVTable");
    iEvent.put(std::move(dauTable), "GVDaughtersTable");
    return;
  }

  const auto& genParticles = genHandle;
  const auto& secondaryVertices = svHandle;
  const auto& PV0 = pvsIn->front();

  std::vector<float> Hadron_pt, Hadron_eta, Hadron_phi;
  std::vector<float> Hadron_GVx, Hadron_GVy, Hadron_GVz;
  std::vector<float> Hadron_GVx_i, Hadron_GVy_i, Hadron_GVz_i;

  std::vector<int> Hadron_pdgId;
  std::vector<int> Hadron_pdgClass;
  std::vector<int> Hadron_isB;
  std::vector<int> Hadron_isD;
  std::vector<int> Hadron_isBtoD;

  std::vector<int> Hadron_nDauNoRecognizedSecondaryAncestor;
  std::vector<int> Hadron_nDauFromB;
  std::vector<int> Hadron_nDauFromBC;
  std::vector<int> Hadron_nDauFromC;
  std::vector<int> Hadron_nDauOtherSecondary;
  std::vector<int> Hadron_nDauOriginUnknown;

  std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
  std::vector<int> Daughters_charge;
  std::vector<int> Daughters_GVidx;
  std::vector<int> Daughters_originLabel;

  std::vector<float> SV_x, SV_y, SV_z;
  std::vector<CovMatrix> SV_cov;

  VertexDistance3D vdist;

  for (const auto& sv : *secondaryVertices) {
    Measurement1D dl = vdist.distance(
        PV0,
        VertexState(
            RecoVertex::convertPos(sv.position()),
            RecoVertex::convertError(sv.error())));

    if (dl.value() > 0. && dl.significance() > dlenSigMin_) {
      SV_x.push_back(sv.x());
      SV_y.push_back(sv.y());
      SV_z.push_back(sv.z());
      SV_cov.push_back(sv.covariance());
    }
  }

  int ngv = 0;

  for (size_t i = 0; i < genParticles->size(); ++i) {
    const reco::Candidate* hadron = &(*genParticles)[i];

    if (!(hadron->pt() > 10. && std::abs(hadron->eta()) < 2.5))
      continue;

    const int hadPDG = checkPDG(std::abs(hadron->pdgId()));
    if (hadPDG == 0)
      continue;

    std::vector<float> temp_pt, temp_eta, temp_phi;
    std::vector<int> temp_charge, temp_GVidx, temp_originLabel;

    int nPack = 0;

    float vx = std::numeric_limits<float>::quiet_NaN();
    float vy = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();

    for (size_t j = 0; j < genParticles->size(); ++j) {
      const reco::Candidate* dau = &(*genParticles)[j];

      if (dau == hadron)
        continue;

      if (!(dau->status() == 1 &&
            dau->charge() != 0 &&
            dau->pt() > 0.8 &&
            std::abs(dau->eta()) < 2.5))
        continue;

      auto GV = isAncestor(hadron, dau);

      if (!GV.has_value())
        continue;

      std::tie(vx, vy, vz) = *GV;

      if (std::isnan(vx))
        continue;

      ++nPack;

      temp_pt.push_back(dau->pt());
      temp_eta.push_back(dau->eta());
      temp_phi.push_back(dau->phi());
      temp_charge.push_back(dau->charge());
      temp_GVidx.push_back(ngv);
      temp_originLabel.push_back(getDaughterOriginLabelNoPU(dau));
    }

    if (nPack < 2)
      continue;

    Hadron_pt.push_back(hadron->pt());
    Hadron_eta.push_back(hadron->eta());
    Hadron_phi.push_back(hadron->phi());
    Hadron_pdgId.push_back(hadron->pdgId());
    Hadron_pdgClass.push_back(hadPDG);

    Hadron_isB.push_back(hadPDG == 1 ? 1 : 0);
    Hadron_isD.push_back(hadPDG == 2 ? 1 : 0);

    if (hadPDG == 2)
      Hadron_isBtoD.push_back(hasBHadronAncestor(hadron) ? 1 : 0);
    else
      Hadron_isBtoD.push_back(0);

    Hadron_GVx.push_back(vx);
    Hadron_GVy.push_back(vy);
    Hadron_GVz.push_back(vz);

    Hadron_GVx_i.push_back(hadron->vx());
    Hadron_GVy_i.push_back(hadron->vy());
    Hadron_GVz_i.push_back(hadron->vz());

    int nNoRecognizedSecondaryAncestor = 0;
    int nFromB = 0;
    int nFromBC = 0;
    int nFromC = 0;
    int nOtherSecondary = 0;
    int nOriginUnknown = 0;

    for (const int lab : temp_originLabel) {
      if (lab == 0)
        ++nNoRecognizedSecondaryAncestor;
      else if (lab == 2)
        ++nFromB;
      else if (lab == 3)
        ++nFromBC;
      else if (lab == 4)
        ++nFromC;
      else if (lab == 5)
        ++nOtherSecondary;
      else
        ++nOriginUnknown;
    }

    Hadron_nDauNoRecognizedSecondaryAncestor.push_back(nNoRecognizedSecondaryAncestor);
    Hadron_nDauFromB.push_back(nFromB);
    Hadron_nDauFromBC.push_back(nFromBC);
    Hadron_nDauFromC.push_back(nFromC);
    Hadron_nDauOtherSecondary.push_back(nOtherSecondary);
    Hadron_nDauOriginUnknown.push_back(nOriginUnknown);

    Daughters_pt.insert(Daughters_pt.end(), temp_pt.begin(), temp_pt.end());
    Daughters_eta.insert(Daughters_eta.end(), temp_eta.begin(), temp_eta.end());
    Daughters_phi.insert(Daughters_phi.end(), temp_phi.begin(), temp_phi.end());
    Daughters_charge.insert(Daughters_charge.end(), temp_charge.begin(), temp_charge.end());
    Daughters_GVidx.insert(Daughters_GVidx.end(), temp_GVidx.begin(), temp_GVidx.end());
    Daughters_originLabel.insert(Daughters_originLabel.end(), temp_originLabel.begin(), temp_originLabel.end());

    ++ngv;
  }

  std::vector<float> SVtrk_pt, SVtrk_eta, SVtrk_phi;
  std::vector<int> SVtrk_SVidx;

  int SV_index = 0;

  for (const auto& sv : *secondaryVertices) {
    Measurement1D dl = vdist.distance(
        PV0,
        VertexState(
            RecoVertex::convertPos(sv.position()),
            RecoVertex::convertError(sv.error())));

    if (dl.value() > 0. && dl.significance() > dlenSigMin_) {
      for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
        const edm::RefToBase<reco::Track>& trkRef = *it;

        if (trkRef.isNull())
          continue;

        SVtrk_pt.push_back(trkRef->pt());
        SVtrk_eta.push_back(trkRef->eta());
        SVtrk_phi.push_back(trkRef->phi());
        SVtrk_SVidx.push_back(SV_index);
      }

      ++SV_index;
    }
  }

  auto distances = computeDistanceMatrix(
      SV_x,
      SV_y,
      SV_z,
      SV_cov,
      Hadron_GVx,
      Hadron_GVy,
      Hadron_GVz);

  std::vector<int> Hadron_SVIdx(ngv, -1);
  std::vector<float> Hadron_SVDistance(ngv, -1.f);
  std::vector<float> Hadron_minDistNotMatched(ngv, 999.f);

  auto result = matchHadronsToSV(
      distances,
      SVtrk_pt,
      SVtrk_eta,
      SVtrk_phi,
      SVtrk_SVidx,
      Daughters_pt,
      Daughters_eta,
      Daughters_phi,
      Daughters_GVidx,
      ngv,
      nRequiredCommonTracks_,
      dR_max_,
      relPt_max_,
      doubleMatching_,
      doubleMatching_maxSignificance_,
      doubleMatching_dR_max_,
      doubleMatching_relPt_max_);

  Hadron_SVIdx = std::get<0>(result);
  Hadron_SVDistance = std::get<1>(result);
  Hadron_minDistNotMatched = std::get<2>(result);

  auto gvTable = std::make_unique<nanoaod::FlatTable>(ngv, "GV", false);

  gvTable->addColumn<float>("pt", Hadron_pt, "Hadron pt");
  gvTable->addColumn<float>("eta", Hadron_eta, "Hadron eta");
  gvTable->addColumn<float>("phi", Hadron_phi, "Hadron phi");

  gvTable->addColumn<float>("x", Hadron_GVx, "GV x");
  gvTable->addColumn<float>("y", Hadron_GVy, "GV y");
  gvTable->addColumn<float>("z", Hadron_GVz, "GV z");

  gvTable->addColumn<float>("x_i", Hadron_GVx_i, "Hadron production x coordinate");
  gvTable->addColumn<float>("y_i", Hadron_GVy_i, "Hadron production y coordinate");
  gvTable->addColumn<float>("z_i", Hadron_GVz_i, "Hadron production z coordinate");

  gvTable->addColumn<int>("Hadron_SVIdx", Hadron_SVIdx, "Matched SV index");
  gvTable->addColumn<int>("Hadron_pdgId", Hadron_pdgId, "Hadron PDG id");
  gvTable->addColumn<float>("SV_distanceSig", Hadron_SVDistance, "Matched SV-GV distance significance");

  gvTable->addColumn<int>("isB", Hadron_isB, "Hadron is B");
  gvTable->addColumn<int>("isD", Hadron_isD, "Hadron is D");
  gvTable->addColumn<int>("isBtoD", Hadron_isBtoD, "D hadron has B-hadron ancestor");
  gvTable->addColumn<int>("pdgClass", Hadron_pdgClass, "Hadron class: 1 B, 2 D, 3 strange, 4 tau");
  gvTable->addColumn<float>("minDistNotMatched", Hadron_minDistNotMatched, "Minimum distance to SV among unmatched hadrons");

  gvTable->addColumn<int>(
      "nDauNoRecognizedSecondaryAncestor",
      Hadron_nDauNoRecognizedSecondaryAncestor,
      "Selected daughters with no recognized B/C/strange/tau/conversion-like ancestry; may include hard-scatter primary-like or pileup-like particles");

  gvTable->addColumn<int>("nDauFromB", Hadron_nDauFromB, "Selected daughters with B ancestry");
  gvTable->addColumn<int>("nDauFromBC", Hadron_nDauFromBC, "Selected daughters with B and C ancestry");
  gvTable->addColumn<int>("nDauFromC", Hadron_nDauFromC, "Selected daughters with C ancestry");
  gvTable->addColumn<int>("nDauOtherSecondary", Hadron_nDauOtherSecondary, "Selected daughters with strange/tau/conversion-like ancestry");
  gvTable->addColumn<int>("nDauOriginUnknown", Hadron_nDauOriginUnknown, "Selected daughters with unknown origin label");

  auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(), "GVDaughters", false);

  dauTable->addColumn<float>("pt", Daughters_pt, "Daughter pt");
  dauTable->addColumn<float>("eta", Daughters_eta, "Daughter eta");
  dauTable->addColumn<float>("phi", Daughters_phi, "Daughter phi");
  dauTable->addColumn<int>("charge", Daughters_charge, "Daughter charge");
  dauTable->addColumn<int>("hadronIndex", Daughters_GVidx, "Hadron index");

  dauTable->addColumn<int>(
      "originLabel",
      Daughters_originLabel,
      "Daughter origin label without collisionId: 0 noRecognizedSecondaryAncestor/primary-or-pileup-like, 2 fromB, 3 fromBC, 4 fromC, 5 otherSecondary, 9 unknown");

  iEvent.put(std::move(gvTable), "GVTable");
  iEvent.put(std::move(dauTable), "GVDaughtersTable");
}

int GenVertexProducer::checkPDG(int abs_pdg) const {
  std::vector<int> pdgList_B = {
      521, 511, 531, 541, 5122, 5132, 5232, 5332,
      5142, 5242, 5342, 5512, 5532, 5542, 5554};

  std::vector<int> pdgList_D = {
      411, 421, 431, 4122, 4232, 4132, 4332,
      4412, 4422, 4432, 4444};

  std::vector<int> pdgList_S = {
      3122, 3222, 3212, 3312, 3322, 3334};

  std::vector<int> pdgList_Tau = {15};

  if (std::find(pdgList_B.begin(), pdgList_B.end(), abs_pdg) != pdgList_B.end())
    return 1;

  if (std::find(pdgList_D.begin(), pdgList_D.end(), abs_pdg) != pdgList_D.end())
    return 2;

  if (std::find(pdgList_S.begin(), pdgList_S.end(), abs_pdg) != pdgList_S.end())
    return 3;

  if (std::find(pdgList_Tau.begin(), pdgList_Tau.end(), abs_pdg) != pdgList_Tau.end())
    return 4;

  return 0;
}

bool GenVertexProducer::hasBHadronAncestor(const reco::Candidate* cand) const {
  if (cand == nullptr)
    return false;

  const reco::Candidate* current = cand;

  while (current != nullptr && current->numberOfMothers() > 0) {
    const reco::Candidate* mother = current->mother(0);

    if (mother == nullptr || mother == current)
      break;

    if (checkPDG(std::abs(mother->pdgId())) == 1)
      return true;

    current = mother;
  }

  return false;
}

int GenVertexProducer::getDaughterOriginLabelNoPU(const reco::Candidate* dau) const {
  if (dau == nullptr)
    return 9;

  static const std::unordered_set<int> pdgSet_B = {
      521, 511, 531, 541, 5122, 5132, 5232, 5332,
      5142, 5242, 5342, 5512, 5532, 5542, 5554};

  static const std::unordered_set<int> pdgSet_C = {
      411, 421, 431, 4122, 4232, 4132, 4332,
      4412, 4422, 4432, 4444};

  static const std::unordered_set<int> pdgSet_S = {
      310, 130,
      3122, 3222, 3212, 3312, 3322, 3334};

  bool foundB = false;
  bool foundC = false;
  bool foundOtherSecondary = false;

  const reco::Candidate* cur = dau;
  int guard = 0;

  while (cur != nullptr && cur->numberOfMothers() > 0 && guard++ < 100) {
    const reco::Candidate* mom = cur->mother(0);

    if (mom == nullptr || mom == cur)
      break;

    const int apdg = std::abs(mom->pdgId());

    if (pdgSet_B.count(apdg))
      foundB = true;

    if (pdgSet_C.count(apdg))
      foundC = true;

    if (pdgSet_S.count(apdg) || apdg == 15 || apdg == 22)
      foundOtherSecondary = true;

    cur = mom;
  }

  if (foundB && foundC)
    return 3;

  if (foundB)
    return 2;

  if (foundC)
    return 4;

  if (foundOtherSecondary)
    return 5;

  return 0;
}

std::optional<std::tuple<float, float, float>> GenVertexProducer::isAncestor(
    const reco::Candidate* ancestor,
    const reco::Candidate* particle) const {
  std::unordered_set<int> pdgSet_B = {
      521, 511, 531, 541, 5122, 5132, 5232, 5332,
      5142, 5242, 5342, 5512, 5532, 5542, 5554};

  std::unordered_set<int> pdgSet_D = {
      411, 421, 431, 4122, 4232, 4132, 4332,
      4412, 4422, 4432, 4444};

  std::unordered_set<int> pdgSet_S = {
      3122, 3222, 3212, 3312, 3322, 3334};

  std::unordered_set<int> pdgSet_Tau = {15};

  const reco::Candidate* current = particle;

  while (current != nullptr && current->numberOfMothers() > 0) {
    const reco::Candidate* mother = current->mother(0);

    if (mother == ancestor) {
      return std::make_optional(std::make_tuple(
          current->vx(),
          current->vy(),
          current->vz()));
    }

    if (mother == nullptr || mother == current)
      break;

    const int mother_pdg = std::abs(mother->pdgId());

    if (pdgSet_B.count(mother_pdg) ||
        pdgSet_D.count(mother_pdg) ||
        pdgSet_S.count(mother_pdg) ||
        pdgSet_Tau.count(mother_pdg))
      break;

    current = mother;
  }

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
  const size_t nSV = SV_x.size();
  const size_t nHadron = Hadron_GVx.size();

  std::vector<std::vector<float>> distances(
      nSV,
      std::vector<float>(nHadron, 999.0f));

  for (size_t i = 0; i < nSV; ++i) {
    CovMatrix covInv = SV_cov[i];
    covInv.Invert();

    for (size_t j = 0; j < nHadron; ++j) {
      const float dx = SV_x[i] - Hadron_GVx[j];
      const float dy = SV_y[i] - Hadron_GVy[j];
      const float dz = SV_z[i] - Hadron_GVz[j];

      const float chi2 =
          dx * (covInv(0, 0) * dx +
                covInv(0, 1) * dy +
                covInv(0, 2) * dz) +
          dy * (covInv(1, 0) * dx +
                covInv(1, 1) * dy +
                covInv(1, 2) * dz) +
          dz * (covInv(2, 0) * dx +
                covInv(2, 1) * dy +
                covInv(2, 2) * dz);

      float dist = 999.0f;

      if (std::isfinite(chi2) && chi2 >= 0.0f)
        dist = std::sqrt(chi2);

      distances[i][j] = dist;
    }
  }

  return distances;
}

std::tuple<std::vector<int>, std::vector<float>, std::vector<float>>
GenVertexProducer::matchHadronsToSV(
    std::vector<std::vector<float>> distances,
    const std::vector<float>& SVtrk_pt,
    const std::vector<float>& SVtrk_eta,
    const std::vector<float>& SVtrk_phi,
    const std::vector<int>& SVtrk_SVidx,
    const std::vector<float>& Daughters_pt,
    const std::vector<float>& Daughters_eta,
    const std::vector<float>& Daughters_phi,
    const std::vector<int>& Daughters_GVidx,
    int n_Hadrons,
    int nRequiredCommonTracks,
    double dR_max,
    double relPt_max,
    bool doubleMatching,
    double doubleMatching_maxSignificance,
    double doubleMatching_dR_max,
    double doubleMatching_relPt_max) {
  const size_t nSV = distances.size();

  std::vector<int> Hadron_SVIdx(n_Hadrons, -1);
  std::vector<float> Hadron_SVDistance(n_Hadrons, -1.f);
  std::vector<size_t> svTrackIdxs_fromBestSV;

  std::vector<std::vector<float>> distancesOriginal = distances;

  while (true) {
    float minDist = 999.0f;
    int bestSV = -1;
    int bestHad = -1;

    for (size_t sv = 0; sv < nSV; ++sv) {
      for (int had = 0; had < n_Hadrons; ++had) {
        if (distances[sv][had] < minDist) {
          minDist = distances[sv][had];
          bestSV = static_cast<int>(sv);
          bestHad = had;
        }
      }
    }

    if (minDist >= 997.0f)
      break;

    svTrackIdxs_fromBestSV.clear();

    for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
      if (SVtrk_SVidx[i] == bestSV &&
          SVtrk_pt[i] > 0.8 &&
          std::fabs(SVtrk_eta[i]) < 2.5) {
        svTrackIdxs_fromBestSV.push_back(i);
      }
    }

    std::vector<size_t> GenDaughtersIdxs_fromBestHad;

    for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
      if (Daughters_GVidx[i] == bestHad)
        GenDaughtersIdxs_fromBestHad.push_back(i);
    }

    int common = 0;

    for (size_t iSV : svTrackIdxs_fromBestSV) {
      for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
        const float dR = deltaR(
            SVtrk_eta[iSV],
            SVtrk_phi[iSV],
            Daughters_eta[iHad],
            Daughters_phi[iHad]);

        const float relPt =
            std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) /
            std::max(Daughters_pt[iHad], 1e-6f);

        if (dR < dR_max && relPt < relPt_max) {
          ++common;

          if (common >= nRequiredCommonTracks)
            break;
        }
      }

      if (common >= nRequiredCommonTracks)
        break;
    }

    if (common >= nRequiredCommonTracks) {
      Hadron_SVIdx[bestHad] = bestSV;
      Hadron_SVDistance[bestHad] = minDist;

      for (int h = 0; h < n_Hadrons; ++h)
        distances[bestSV][h] = 1000.0f;

      for (size_t s = 0; s < nSV; ++s)
        distances[s][bestHad] = 1000.0f;

    } else {
      distances[bestSV][bestHad] = 998.0f;
    }
  }

  if (doubleMatching) {
    while (true) {
      float minDist = 999.0f;
      int bestSV = -1;
      int bestHad = -1;

      for (size_t sv = 0; sv < nSV; ++sv) {
        for (int had = 0; had < n_Hadrons; ++had) {
          if (distances[sv][had] < minDist) {
            minDist = distances[sv][had];
            bestSV = static_cast<int>(sv);
            bestHad = had;
          }
        }
      }

      if (minDist >= 998.5f)
        break;

      svTrackIdxs_fromBestSV.clear();

      for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
        if (SVtrk_SVidx[i] == bestSV &&
            SVtrk_pt[i] > 0.8 &&
            std::fabs(SVtrk_eta[i]) < 2.5) {
          svTrackIdxs_fromBestSV.push_back(i);
        }
      }

      std::vector<size_t> GenDaughtersIdxs_fromBestHad;

      for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
        if (Daughters_GVidx[i] == bestHad)
          GenDaughtersIdxs_fromBestHad.push_back(i);
      }

      int common = 0;

      for (size_t iSV : svTrackIdxs_fromBestSV) {
        for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
          const float dR = deltaR(
              SVtrk_eta[iSV],
              SVtrk_phi[iSV],
              Daughters_eta[iHad],
              Daughters_phi[iHad]);

          const float relPt =
              std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) /
              std::max(Daughters_pt[iHad], 1e-6f);

          if (dR < doubleMatching_dR_max &&
              relPt < doubleMatching_relPt_max) {
            ++common;

            if (common >= 1)
              break;
          }
        }

        if (common >= 1)
          break;
      }

      if (common >= 1 &&
          distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance) {
        Hadron_SVIdx[bestHad] = bestSV;
        Hadron_SVDistance[bestHad] = -distancesOriginal[bestSV][bestHad];

        for (int h = 0; h < n_Hadrons; ++h)
          distances[bestSV][h] = 1000.0f;

        for (size_t s = 0; s < nSV; ++s)
          distances[s][bestHad] = 1000.0f;

      } else {
        distances[bestSV][bestHad] = 999.0f;
      }
    }
  }

  std::vector<float> minDistNotMatched(n_Hadrons, 999.f);

  for (int had = 0; had < n_Hadrons; ++had) {
    for (size_t sv = 0; sv < nSV; ++sv) {
      if (Hadron_SVIdx[had] == -1 &&
          distancesOriginal[sv][had] < minDistNotMatched[had]) {
        minDistNotMatched[had] = distancesOriginal[sv][had];
      }
    }
  }

  return std::make_tuple(
      Hadron_SVIdx,
      Hadron_SVDistance,
      minDistNotMatched);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenVertexProducer);
