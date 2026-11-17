/**
  \class 
  \brief
  \author   Andrea Rizzi
  \version  $Id: TrackAndVertexUnpacker.cc,v 1.2 2010/02/20 21:00:29 wmtan Exp $
*/


#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

//#include "DataFormats/Common/interface/View.h"
//#include "DataFormats/PatCandidates/interface/Vertexing.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include <vector>

//#include "PhysicsTools/PatAlgos/interface/VertexingHelper.h"


namespace pat {
  
  class PATTrackAndVertexUnpacker : public edm::global::EDProducer<> {
    
    
  public:
    
    explicit PATTrackAndVertexUnpacker(const edm::ParameterSet & iConfig);
    ~PATTrackAndVertexUnpacker() override;
    
    void produce(edm::StreamID, edm::Event & iEvent, const edm::EventSetup& iSetup) const override;
    
  private:
    typedef std::vector<edm::InputTag> VInputTag;
    // configurables
    const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >    Cands_;
    const edm::EDGetTokenT<reco::VertexCollection>         PVs_;
    const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection>         SVs_;
    const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >         AdditionalTracks_;
	const double minTrackPt_;   // <-- new member
	const double maxTrackEta_; 
    //////    std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > > particlesTokens_;
    
  };

}

using pat::PATTrackAndVertexUnpacker;

PATTrackAndVertexUnpacker::PATTrackAndVertexUnpacker(const edm::ParameterSet& iConfig) :
  Cands_(consumes< std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("packedCandidates"))),
  PVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("slimmedVertices"))),
  SVs_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("slimmedSecondaryVertices"))),
  AdditionalTracks_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("additionalTracks"))),
  minTrackPt_(iConfig.getParameter<double>("minTrackPt")),   // <-- new
  maxTrackEta_(iConfig.getParameter<double>("maxTrackEta"))   
{
    produces<reco::TrackCollection>();
    produces<reco::VertexCollection>();
    produces<reco::VertexCollection>("secondary");
    produces<edm::ValueMap<int>>("packedTrackIndex");
    produces<edm::ValueMap<int>>("lostTrackIndex");
}


PATTrackAndVertexUnpacker::~PATTrackAndVertexUnpacker() {
}


void PATTrackAndVertexUnpacker::produce(edm::StreamID, edm::Event & iEvent, const edm::EventSetup & iSetup) const {
	using namespace edm; using namespace std; using namespace reco;
	Handle<std::vector<pat::PackedCandidate> > cands;
	iEvent.getByToken(Cands_, cands);
	Handle<VertexCollection> pvs;
	iEvent.getByToken(PVs_, pvs);
	Handle<VertexCompositePtrCandidateCollection> svs;
	iEvent.getByToken(SVs_, svs);
	Handle<std::vector<pat::PackedCandidate> > addTracks;
	bool hasAdditionalTracks = iEvent.getByToken(AdditionalTracks_, addTracks);
	if (hasAdditionalTracks) {
	iEvent.getByToken(AdditionalTracks_, addTracks);
	}
	auto outTks = std::make_unique<std::vector<reco::Track>>();
	std::map<unsigned int, std::vector<unsigned int> > asso;
	std::map<unsigned int, unsigned int> trackKeys;
	std::map<unsigned int, unsigned int> addTrackKeys; // new index to map additional tracks index to the output track collection index
	std::vector<int> packedToUnpacked(cands->size(), -1);
	std::vector<int> lostToUnpacked(hasAdditionalTracks ? addTracks->size() : 0, -1);
	unsigned int j=0;
	
	// Packed PF Candidates
	for(unsigned int i=0;i<cands->size();i++)	{
		const pat::PackedCandidate & c = (*cands)[i];
		if(c.hasTrackDetails() && c.charge() != 0 && c.numberOfHits()> 0){
			if(c.pt() <= minTrackPt_)
        		continue;
			if(std::abs(c.eta()) >= maxTrackEta_)
				continue;
			outTks->push_back(c.pseudoTrack());
			for(size_t ipv=0;ipv< pvs->size(); ++ipv) {
				if(c.fromPV(ipv)==pat::PackedCandidate::PVUsedInFit)
					asso[ipv].push_back(j);
			}
			trackKeys[i]=j;
			packedToUnpacked[i] = static_cast<int>(j);
			j++;
		}	
	}

	int offsetAdd=j;

	// Additional Tracks
	if (hasAdditionalTracks){
		for(unsigned int i = 0; i < addTracks->size(); i++) {
			if(!(*addTracks)[i].hasTrackDetails())
				continue;
			if((*addTracks)[i].pseudoTrack().pt() <= minTrackPt_)
				continue;   // <-- skip, don't advance j
			if(std::abs((*addTracks)[i].pseudoTrack().eta()) >= maxTrackEta_)
				continue;

				outTks->push_back((*addTracks)[i].pseudoTrack());
					for(size_t ipv=0;ipv< pvs->size(); ++ipv) {
						if((*addTracks)[i].fromPV(ipv)==pat::PackedCandidate::PVUsedInFit)
						asso[ipv].push_back(j);
					}
				addTrackKeys[i] = j;
				lostToUnpacked[i] = static_cast<int>(j);
				j++;
		}
	}
	edm::OrphanHandle< std::vector<reco::Track>  > oh = iEvent.put(std::move(outTks));

	// Preserve source collection keys so SV constituents can be mapped through
	// the later filtered-track selection without guessing from kinematics.
	auto packedIndexMap = std::make_unique<edm::ValueMap<int>>();
	edm::ValueMap<int>::Filler packedFiller(*packedIndexMap);
	packedFiller.insert(cands, packedToUnpacked.begin(), packedToUnpacked.end());
	packedFiller.fill();
	iEvent.put(std::move(packedIndexMap), "packedTrackIndex");

	auto lostIndexMap = std::make_unique<edm::ValueMap<int>>();
	if (hasAdditionalTracks) {
		edm::ValueMap<int>::Filler lostFiller(*lostIndexMap);
		lostFiller.insert(addTracks, lostToUnpacked.begin(), lostToUnpacked.end());
		lostFiller.fill();
	}
	iEvent.put(std::move(lostIndexMap), "lostTrackIndex");
	
	auto outPv = std::make_unique<std::vector<reco::Vertex>>();
	
	for(size_t ipv=0;ipv< pvs->size(); ++ipv) {
		reco::Vertex  pv = (*pvs)[ipv];
		for(unsigned int i=0;i<asso[ipv].size();i++)
		{
			TrackRef r(oh,asso[ipv][i]); //asso[ipv][i] = index wrt oh collection of the track number i-th associated to the vertex number ipv
			TrackBaseRef rr(r);
			pv.add(rr);
		}
		outPv->push_back(pv);
	}
	iEvent.put(std::move(outPv));

	auto outSv = std::make_unique<std::vector<reco::Vertex>>();
	for(size_t i=0;i< svs->size(); i++) {
	const reco::VertexCompositePtrCandidate &sv = (*svs)[i];	
	outSv->push_back(reco::Vertex(sv.vertex(),sv.vertexCovariance(),sv.vertexChi2(),sv.vertexNdof(),0));
	// last 0 is for pre-reserving any capacity for tracks, but we will add them later
	for(size_t j=0;j<sv.numberOfDaughters();j++){
		TrackRef r;
		if(sv.daughterPtr(j).id() == cands.id()) {
			// Daughter from PackedCandidate collection
			// sv.daughterPtr(j).key() returns the index of the PackedCollection
			auto found = trackKeys.find(sv.daughterPtr(j).key());
			if (found != trackKeys.end()) r = TrackRef(oh, found->second);
		} else if (hasAdditionalTracks && sv.daughterPtr(j).id() == addTracks.id()) {
			// Daughter from additionalTracks collection
//				std::cout << "vertex " << i << " using lost Track " << sv.daughterPtr(j).key()  << "  " << offsetAdd+sv.daughterPtr(j).key() << std::endl;  
			auto found = addTrackKeys.find(sv.daughterPtr(j).key());
			if (found != addTrackKeys.end()) r = TrackRef(oh, found->second);
		}
		// Preserve the original candidate's full daughter list in downstream
		// tables; only add references here for daughters present in outTks.
		if (r.isNonnull()) outSv->back().add(TrackBaseRef(r));

		}	
	}   

       iEvent.put(std::move(outSv),"secondary");

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTrackAndVertexUnpacker);
