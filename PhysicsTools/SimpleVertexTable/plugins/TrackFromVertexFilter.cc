#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace pat {

  class PATTrackFromVertexFilter : public edm::global::EDProducer<> {

  public:
    explicit PATTrackFromVertexFilter(const edm::ParameterSet& iConfig);
    ~PATTrackFromVertexFilter() override = default;

    void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  private:
    const edm::EDGetTokenT<reco::TrackCollection>  tracks_;
    const edm::EDGetTokenT<reco::VertexCollection> pvs_;
    const double dzCut_;
  };

} // namespace pat

using pat::PATTrackFromVertexFilter;

PATTrackFromVertexFilter::PATTrackFromVertexFilter(const edm::ParameterSet& iConfig)
  : tracks_(consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("tracks"))),
    pvs_   (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    dzCut_ (iConfig.getParameter<double>("dzCut"))
{
  produces<reco::TrackCollection>();
  produces<reco::VertexCollection>();
  produces<edm::ValueMap<int>>("filteredIndex");
}

void PATTrackFromVertexFilter::produce(edm::StreamID,
                                       edm::Event& iEvent,
                                       const edm::EventSetup& iSetup) const
{
  using namespace edm;
  using namespace reco;

  Handle<TrackCollection> tracks;
  iEvent.getByToken(tracks_, tracks);

  Handle<VertexCollection> pvs;
  iEvent.getByToken(pvs_, pvs);

  auto outTracks = std::make_unique<TrackCollection>();
  auto outPVs    = std::make_unique<VertexCollection>();
  std::vector<int> inputToFilteredIndex(tracks->size(), -1);

  if (pvs->empty()) {
    iEvent.put(std::move(outTracks));
    iEvent.put(std::move(outPVs));
    auto indexMap = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*indexMap);
    filler.insert(tracks, inputToFilteredIndex.begin(), inputToFilteredIndex.end());
    filler.fill();
    iEvent.put(std::move(indexMap), "filteredIndex");
    return;
  }

  // Build key set from ALL pvs
  std::unordered_set<unsigned int> pvTrackKeys;
  for (const auto& vertex : *pvs)
    for (auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it)
      pvTrackKeys.insert(it->key());

  // dz reference point is always PV[0]
  const auto& pvPos = (*pvs)[0].position();

  // Build filtered track collection and old->new index map
  std::unordered_map<unsigned int, unsigned int> oldToNew;
  oldToNew.reserve(tracks->size());

  unsigned int newIdx = 0;
  for (unsigned int i = 0; i < tracks->size(); ++i) {
    const reco::Track& tk = (*tracks)[i];
    // method 1 keep if from PV or satisfiies dzcut
    //const bool keep = pvTrackKeys.count(i)
    //                  || std::abs(tk.dz(pvPos)) < dzCut_;
    const bool keep =  std::abs(tk.dz(pvPos)) < dzCut_;
    if (keep) {
      oldToNew[i] = newIdx++;
      inputToFilteredIndex[i] = static_cast<int>(oldToNew[i]);
      outTracks->push_back(tk);
    }
  }

  // Put tracks first to get OrphanHandle for TrackRef building
  OrphanHandle<TrackCollection> oh = iEvent.put(std::move(outTracks));

  // Map each input-track key to its row in the filtered output. Tracks rejected
  // by the dz selection retain the -1 sentinel.
  auto indexMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler filler(*indexMap);
  filler.insert(tracks, inputToFilteredIndex.begin(), inputToFilteredIndex.end());
  filler.fill();
  iEvent.put(std::move(indexMap), "filteredIndex");

  // Rebuild PV collection with remapped track refs
  for (const auto& oldPV : *pvs) {
    reco::Vertex newPV(oldPV.position(),
                       oldPV.covariance(),
                       oldPV.chi2(),
                       oldPV.ndof(),
                       0);
    for (auto it = oldPV.tracks_begin(); it != oldPV.tracks_end(); ++it) {
      auto found = oldToNew.find(it->key());
      if (found != oldToNew.end()) {
        TrackRef newRef(oh, found->second);
        newPV.add(TrackBaseRef(newRef), oldPV.trackWeight(*it));
      }
    }
    outPVs->push_back(newPV);
  }

  iEvent.put(std::move(outPVs));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATTrackFromVertexFilter);
