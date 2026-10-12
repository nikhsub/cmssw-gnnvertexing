// Kinematic one-to-one association between the package's unpacked reco tracks
// and the exact event-local QGenCharged rows written by Stage 1. This is not a
// simulated-hit or independently verified detector truth association.

#include "PhysicsTools/SimpleVertexTable/interface/QuarkRecoTrackAssignment.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

namespace {
  template <typename T>
  std::span<const T> requiredColumn(const nanoaod::FlatTable& table, const std::string& name) {
    const int index = table.columnIndex(name);
    if (index < 0)
      throw cms::Exception("ProductNotFound") << "QGenCharged FlatTable is missing required column '" << name << "'";
    return table.columnData<T>(static_cast<unsigned int>(index));
  }

  bool validKinematics(double pt, double eta, double phi) {
    return std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi) && pt > 0.;
  }
}  // namespace

class QuarkRecoTrackMatchProducer : public edm::stream::EDProducer<> {
public:
  explicit QuarkRecoTrackMatchProducer(const edm::ParameterSet& config)
      : tracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("tracks"))),
        genChargedToken_(consumes<nanoaod::FlatTable>(config.getParameter<edm::InputTag>("genChargedTable"))),
        maxDeltaR_(config.getParameter<double>("maxDeltaR")),
        maxRelPt_(config.getParameter<double>("maxRelPt")) {
    if (!std::isfinite(maxDeltaR_) || maxDeltaR_ <= 0. || !std::isfinite(maxRelPt_) || maxRelPt_ <= 0.)
      throw cms::Exception("Configuration") << "maxDeltaR and maxRelPt must be finite and positive";
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& tracks = event.get(tracksToken_);
    const auto& genTable = event.get(genChargedToken_);
    const auto genPt = requiredColumn<float>(genTable, "pt");
    const auto genEta = requiredColumn<float>(genTable, "eta");
    const auto genPhi = requiredColumn<float>(genTable, "phi");
    const auto genCharge = requiredColumn<int>(genTable, "charge");
    const unsigned int nTracks = tracks.size();
    const unsigned int nGen = genTable.nRows();
    if (genPt.size() != nGen || genEta.size() != nGen || genPhi.size() != nGen || genCharge.size() != nGen)
      throw cms::Exception("LogicError") << "QGenCharged FlatTable columns have inconsistent row counts";

    std::vector<int> nCandidates(nTracks, 0);
    std::vector<int> nTracksPerGen(nGen, 0);
    std::vector<qfrag::RecoGenCandidate> candidates;
    for (unsigned int itrack = 0; itrack < nTracks; ++itrack) {
      const auto& track = tracks[itrack];
      if (track.charge() == 0 || !validKinematics(track.pt(), track.eta(), track.phi())) continue;
      for (unsigned int igen = 0; igen < nGen; ++igen) {
        if (genCharge[igen] == 0 || track.charge() != genCharge[igen] ||
            !validKinematics(genPt[igen], genEta[igen], genPhi[igen]))
          continue;
        const double dr = reco::deltaR(track.eta(), track.phi(), genEta[igen], genPhi[igen]);
        const double relPt = std::abs(track.pt() - genPt[igen]) / genPt[igen];
        if (!std::isfinite(dr) || !std::isfinite(relPt) || dr > maxDeltaR_ || relPt > maxRelPt_) continue;
        const double normalizedCost = (dr / maxDeltaR_) * (dr / maxDeltaR_) +
                                      (relPt / maxRelPt_) * (relPt / maxRelPt_);
        candidates.push_back({static_cast<int>(itrack), static_cast<int>(igen),
                              static_cast<std::int64_t>(std::llround(normalizedCost * 1.e12))});
        ++nCandidates[itrack];
        ++nTracksPerGen[igen];
      }
    }

    const auto assignment = qfrag::maximumCardinalityMinimumCostAssignment(
        static_cast<int>(nTracks), static_cast<int>(nGen), std::move(candidates));
    std::vector<float> pt, eta, phi, matchDR, matchRelPt;
    std::vector<int> charge, isHighPurity, trackIdx, genChargedIdx, matchStatus;
    pt.reserve(nTracks); eta.reserve(nTracks); phi.reserve(nTracks); charge.reserve(nTracks);
    isHighPurity.reserve(nTracks); trackIdx.reserve(nTracks); genChargedIdx.reserve(nTracks); matchDR.reserve(nTracks);
    matchRelPt.reserve(nTracks); matchStatus.reserve(nTracks);
    for (unsigned int itrack = 0; itrack < nTracks; ++itrack) {
      const auto& track = tracks[itrack];
      const int igen = assignment[itrack];
      trackIdx.push_back(static_cast<int>(itrack));
      pt.push_back(track.pt()); eta.push_back(track.eta()); phi.push_back(track.phi()); charge.push_back(track.charge());
      isHighPurity.push_back(track.quality(reco::TrackBase::highPurity) ? 1 : 0);
      genChargedIdx.push_back(igen);
      if (nCandidates[itrack] == 0) {
        matchDR.push_back(-1.f);
        matchRelPt.push_back(-1.f);
        matchStatus.push_back(0);  // no admissible generator candidate
      } else if (igen < 0) {
        matchDR.push_back(-1.f);
        matchRelPt.push_back(-1.f);
        matchStatus.push_back(3);  // eligible, but lost a one-to-one assignment conflict
      } else {
        const double dr = reco::deltaR(track.eta(), track.phi(), genEta[igen], genPhi[igen]);
        const double relPt = std::abs(track.pt() - genPt[igen]) / genPt[igen];
        matchDR.push_back(static_cast<float>(dr));
        matchRelPt.push_back(static_cast<float>(relPt));
        // Code 2 records any competing admissible edge in the chosen pair:
        // either this track had other generator candidates or this generator
        // particle was also admissible for another track.
        matchStatus.push_back(nCandidates[itrack] > 1 || nTracksPerGen[igen] > 1 ? 2 : 1);
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nTracks, "QRecoTrack", false, false);
    table->addColumn<int>("trackIdx", trackIdx, "index in the unpackedTracksAndVertices reco::TrackCollection");
    table->addColumn<float>("pt", pt, "reconstructed track transverse momentum [GeV]");
    table->addColumn<float>("eta", eta, "reconstructed track pseudorapidity");
    table->addColumn<float>("phi", phi, "reconstructed track azimuth [rad]");
    table->addColumn<int>("charge", charge, "reconstructed track charge in units of e");
    table->addColumn<int>("isHighPurity", isHighPurity, "reco::Track highPurity quality bit; not required for matching");
    table->addColumn<int>("genChargedIdx", genChargedIdx, "event-local QGenCharged row, or -1 when unmatched");
    table->addColumn<float>("matchDR", matchDR, "reco-gen deltaR, or -1 when unmatched");
    table->addColumn<float>("matchRelPt", matchRelPt, "abs(recoPt-genPt)/genPt, or -1 when unmatched");
    table->addColumn<int>("matchStatus", matchStatus,
        "0 no candidate; 1 unopposed match; 2 match with competing candidates; 3 unmatched one-to-one conflict");
    table->addColumn<int>("nGenCandidates", nCandidates, "admissible QGenCharged rows before global assignment");
    event.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
    desc.add<edm::InputTag>("genChargedTable", edm::InputTag("quarkFragmentationProducer", "QGenChargedTable"));
    desc.add<double>("maxDeltaR", 0.03);
    desc.add<double>("maxRelPt", 0.5);
    descriptions.add("quarkRecoTrackMatchProducer", desc);
  }

private:
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<nanoaod::FlatTable> genChargedToken_;
  const double maxDeltaR_;
  const double maxRelPt_;
};

DEFINE_FWK_MODULE(QuarkRecoTrackMatchProducer);
