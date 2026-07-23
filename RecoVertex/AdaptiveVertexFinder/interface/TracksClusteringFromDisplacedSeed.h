#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//#define VTXDEBUG

class TracksClusteringFromDisplacedSeed {
public:
  struct Cluster {
    GlobalPoint seedPoint;
    reco::TransientTrack seedingTrack;
    std::vector<reco::TransientTrack> tracks;
  };
  TracksClusteringFromDisplacedSeed(const edm::ParameterSet &params);

  std::vector<Cluster> clusters(const reco::Vertex &pv,
    const std::vector<reco::TransientTrack> &selectedTracks,
    const edm::ValueMap<float> *svScores = nullptr,
    const edm::ValueMap<std::vector<float>> *edgeScores = nullptr,
    const edm::ValueMap<std::vector<int>> *edgeIndices = nullptr,
    float edgeScoreThreshold = 0.,
    float seedScoreThreshold = 0.);

private:
  bool trackFilter(const reco::TrackRef &track) const;
  std::pair<std::vector<reco::TransientTrack>, GlobalPoint> nearTracks(
        const reco::TransientTrack &seed,
        const std::vector<reco::TransientTrack> &tracks,
        const reco::Vertex &primaryVertex,
        const edm::ValueMap<float> *svScores,
        const edm::ValueMap<std::vector<float>> *edgeScores,
        const edm::ValueMap<std::vector<int>> *edgeIndices,
        float edgeScoreThreshold,
        float seedScoreThreshold) const;

  //	unsigned int				maxNTracks;
  double max3DIPSignificance;
  double max3DIPValue;
  double min3DIPSignificance;
  double min3DIPValue;
  double clusterMaxDistance;
  double clusterMaxSignificance;
  double distanceRatio;
  double clusterMinAngleCosine;
  double maxTimeSignificance;
};
