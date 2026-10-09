#ifndef TrackGVFeaturesProducer_h
#define TrackGVFeaturesProducer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <vector>

class TrackGVFeaturesProducer : public edm::stream::EDProducer<> {
public:
    explicit TrackGVFeaturesProducer(const edm::ParameterSet&);
    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    std::vector<std::vector<float>> computeDistanceMatrix(
        const std::vector<float>& SV_x, const std::vector<float>& SV_y, const std::vector<float>& SV_z,
        const std::vector<reco::Vertex::CovarianceMatrix>& SV_cov,
        const std::vector<float>& GV_x, const std::vector<float>& GV_y, const std::vector<float>& GV_z) const;

    edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
    edm::EDGetTokenT<nanoaod::FlatTable> gvTableToken_;
    edm::EDGetTokenT<nanoaod::FlatTable> gvDauTableToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

    double trkPtCut_;
    double dRMatchMax_;
    double relPtMatchMax_;
    double dlenSigMin_;
};
#endif
