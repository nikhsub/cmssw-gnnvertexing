// MC-only supplementary fits for reco tracks assigned to the same QHadron.
// Ordinary trackPair rows are read-only inputs; this producer never changes them.
#include "PhysicsTools/SimpleVertexTable/interface/TwoTrackVertexFitter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
  template <typename T>
  std::span<const T> col(const nanoaod::FlatTable& t, const std::string& n) {
    const int i = t.columnIndex(n);
    if (i < 0) throw cms::Exception("ProductNotFound") << "Missing required flat-table column '" << n << "'";
    return t.columnData<T>(static_cast<unsigned int>(i));
  }
  std::uint64_t pairKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
  }
  bool isB(int pdg) {
    switch (std::abs(pdg)) {
      case 511: case 521: case 531: case 541: case 5122: case 5132: case 5232:
      case 5332: case 5142: case 5242: case 5342: case 5512: case 5532: case 5542: case 5554: return true;
      default: return false;
    }
  }
  bool isD(int pdg) {
    switch (std::abs(pdg)) {
      case 411: case 421: case 431: case 4122: case 4232: case 4132: case 4332:
      case 4412: case 4422: case 4432: case 4444: return true;
      default: return false;
    }
  }
  int hadronCategory(int pdg, int origin) { return isB(pdg) ? 1 : isD(pdg) ? (origin == 1 ? 2 : origin == 2 ? 3 : 4) : 0; }
  void readFit(const nanoaod::FlatTable& p, size_t row, TwoTrackVertexFitter::Result& r) {
    const auto i = [&](const char* n) { return col<int>(p, n)[row]; };
    const auto f = [&](const char* n) { return col<float>(p, n)[row]; };
    r.fitValid=i("vtxFitValid"); r.status=i("vtxFitStatus"); r.covarianceValid=i("vtxCovarianceValid");
    r.pvValid=i("vtxPVValid"); r.lxyValid=i("vtxLxyValid"); r.signedLxyValid=i("vtxSignedLxyValid");
    r.l3dValid=i("vtxL3DValid"); r.pointing2DValid=i("vtxPointing2DValid"); r.pointing3DValid=i("vtxPointing3DValid");
    r.perpendicularValid=i("vtxPerpendicularValid"); r.chi2=f("vtxChi2"); r.ndof=f("vtxNdof");
    r.chi2Ndof=f("vtxChi2Ndof"); r.probability=f("vtxProb"); r.x=f("vtxX"); r.y=f("vtxY"); r.z=f("vtxZ");
    r.covariance={{f("vtxCovXX"),f("vtxCovXY"),f("vtxCovXZ"),f("vtxCovYY"),f("vtxCovYZ"),f("vtxCovZZ")}};
    r.sigmaX=f("vtxSigmaX"); r.sigmaY=f("vtxSigmaY"); r.sigmaZ=f("vtxSigmaZ");
    r.lxy=f("vtxLxy"); r.lxyError=f("vtxLxyErr"); r.lxySignificance=f("vtxLxySig");
    r.signedLxy=f("vtxSignedLxy"); r.signedLxySignificance=f("vtxSignedLxySig");
    r.l3d=f("vtxL3D"); r.l3dError=f("vtxL3DErr"); r.l3dSignificance=f("vtxL3DSig");
    r.cosPointing2D=f("vtxCosPointing2D"); r.cosPointing3D=f("vtxCosPointing3D");
    r.perpendicular=f("vtxPerpDisplacement"); r.perpendicularSignificance=f("vtxPerpDisplacementSig");
  }
}

class TruthMatchedPairFitProducer : public edm::stream::EDProducer<> {
public:
  explicit TruthMatchedPairFitProducer(const edm::ParameterSet& c)
      : tracksToken_(consumes<reco::TrackCollection>(c.getParameter<edm::InputTag>("tracks"))),
        verticesToken_(consumes<reco::VertexCollection>(c.getParameter<edm::InputTag>("primaryVertices"))),
        genToken_(consumes<reco::GenParticleCollection>(c.getParameter<edm::InputTag>("prunedGenParticles"))),
        hadToken_(consumes<nanoaod::FlatTable>(c.getParameter<edm::InputTag>("hadronTable"))),
        chargedToken_(consumes<nanoaod::FlatTable>(c.getParameter<edm::InputTag>("genChargedTable"))),
        recoToken_(consumes<nanoaod::FlatTable>(c.getParameter<edm::InputTag>("recoTrackTable"))),
        pairToken_(consumes<nanoaod::FlatTable>(c.getParameter<edm::InputTag>("pairTable"))),
        transientToken_(consumes<std::vector<reco::TransientTrack>>(c.getParameter<edm::InputTag>("transientTracks"))) {
    produces<nanoaod::FlatTable>("QTruthMatchedTrack");
    produces<nanoaod::FlatTable>("QTruthMatchedHadron");
    produces<nanoaod::FlatTable>("QTruthMatchedPairFit");
    produces<nanoaod::FlatTable>("QTruthPairFitSummary");
  }

  void produce(edm::Event& ev, const edm::EventSetup&) override {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const auto& tracks=ev.get(tracksToken_); const auto& vertices=ev.get(verticesToken_); const auto& gen=ev.get(genToken_);
    const auto& transient=ev.get(transientToken_);
    const auto& had=ev.get(hadToken_); const auto& charged=ev.get(chargedToken_); const auto& reco=ev.get(recoToken_); const auto& pairs=ev.get(pairToken_);
    const auto hPruned=col<int>(had,"prunedGenIdx"), hPdg=col<int>(had,"pdgId"), hOrigin=col<int>(had,"charmOrigin");
    const bool hadDecayVertexAvailable=had.columnIndex("decayVertexValid")>=0;
    std::span<const int> hDecayValid;
    std::span<const float> hDecayX,hDecayY,hDecayZ;
    if(hadDecayVertexAvailable){hDecayValid=col<int>(had,"decayVertexValid");hDecayX=col<float>(had,"decayX");hDecayY=col<float>(had,"decayY");hDecayZ=col<float>(had,"decayZ");}
    const auto gLast=col<int>(charged,"lastPrunedIdx"), gHad=col<int>(charged,"hadronIdx"), gMember=col<int>(charged,"membershipType"), gPdg=col<int>(charged,"pdgId"), gCharge=col<int>(charged,"charge"), gQuark=col<int>(charged,"quarkIdx"), gQuality=col<int>(charged,"associationQuality"), gSource=col<int>(charged,"fragmentationSourceCertain");
    const auto rTrack=col<int>(reco,"trackIdx"), rGen=col<int>(reco,"genChargedIdx"), rStatus=col<int>(reco,"matchStatus"), rN=col<int>(reco,"nGenCandidates");
    const auto rMatchDR=col<float>(reco,"matchDR"), rMatchRelPt=col<float>(reco,"matchRelPt");
    const auto pA=col<int>(pairs,"trk1Idx"), pB=col<int>(pairs,"trk2Idx");
    const size_t nh=had.nRows(), ng=charged.nRows(), nr=reco.nRows();
    std::unordered_map<const reco::Candidate*, int> genIndex;
    genIndex.reserve(gen.size());
    for (size_t ig=0; ig<gen.size(); ++ig) genIndex.emplace(&gen[ig], static_cast<int>(ig));
    if (hPruned.size()!=nh || hPdg.size()!=nh || hOrigin.size()!=nh || gLast.size()!=ng || gHad.size()!=ng || gMember.size()!=ng ||
        gPdg.size()!=ng || gCharge.size()!=ng || gQuark.size()!=ng || gQuality.size()!=ng || gSource.size()!=ng ||
        rTrack.size()!=nr || rGen.size()!=nr || rStatus.size()!=nr || rN.size()!=nr || rMatchDR.size()!=nr || rMatchRelPt.size()!=nr)
      throw cms::Exception("LogicError") << "Truth matched table row mismatch";
    if(hadDecayVertexAvailable&&(hDecayValid.size()!=nh||hDecayX.size()!=nh||hDecayY.size()!=nh||hDecayZ.size()!=nh))
      throw cms::Exception("LogicError") << "QHadron decay-vertex columns have inconsistent row counts";
    std::vector<int> trackIdx, genIdx, hadIdx, quarkIdx, genPdg, genCharge, matchStatus, nCandidates, daughterRelation, immediateParent, immediateParentPdg, assocQuality, membership, sourceCertain;
    std::vector<float> genX,genY,genZ,matchDR,matchRelPt; std::vector<int> genVertexValid;
    std::vector<std::vector<int>> members(nh);
    trackIdx.reserve(nr); genIdx.reserve(nr); hadIdx.reserve(nr);
    for (size_t it=0; it<nr; ++it) {
      const int ig=rGen[it]; const bool vg=ig>=0 && static_cast<size_t>(ig)<ng;
      const int ih=vg?gHad[ig]:-1; const bool vh=vg && gMember[ig]==1 && ih>=0 && static_cast<size_t>(ih)<nh;
      const int last=vg?gLast[ig]:-1; const bool vlast=last>=0 && static_cast<size_t>(last)<gen.size();
      int rel=vh?3:0, parent=-1; bool hasParent=false, multipleParents=false, direct=false, reaches=false;
      if (vlast) {
        const auto& gp=gen[last];
        for (size_t m=0;m<gp.numberOfMothers();++m) {
          auto* mom=gp.mother(m); if (!mom) continue;
          const auto found=genIndex.find(mom); if(found==genIndex.end()) continue;
          const int k=found->second;
          if (!hasParent) {parent=k;hasParent=true;} else if (parent!=k) multipleParents=true;
          if (vh && k==hPruned[ih]) direct=true;
        }
        if (vh) {
          std::vector<int> todo{last}; std::vector<int> depth{0}; std::unordered_map<int,int> seen{{last,0}};
          for (size_t q=0;q<todo.size();++q) {
            const int node=todo[q]; if (node==hPruned[ih]) { reaches=true; break; }
            for (size_t m=0;m<gen[node].numberOfMothers();++m) {
              auto* mom=gen[node].mother(m); if (!mom) continue;
              const auto found=genIndex.find(mom); if(found!=genIndex.end() && !seen.count(found->second)) {seen[found->second]=depth[q]+1;todo.push_back(found->second);depth.push_back(depth[q]+1);}
            }
          }
          rel=direct?1:reaches?2:3;
        }
        genX.push_back(gp.vx()); genY.push_back(gp.vy()); genZ.push_back(gp.vz());
        genVertexValid.push_back(std::isfinite(gp.vx())&&std::isfinite(gp.vy())&&std::isfinite(gp.vz()));
      } else { genX.push_back(-999.f);genY.push_back(-999.f);genZ.push_back(-999.f);genVertexValid.push_back(0); }
      trackIdx.push_back(rTrack[it]); genIdx.push_back(vg?ig:-1); hadIdx.push_back(vh?ih:-1);quarkIdx.push_back(vg?gQuark[ig]:-1);
      genPdg.push_back(vg?gPdg[ig]:0); genCharge.push_back(vg?gCharge[ig]:0); matchStatus.push_back(rStatus[it]); nCandidates.push_back(rN[it]);
      daughterRelation.push_back(rel); immediateParent.push_back(hasParent&&!multipleParents?parent:-1);immediateParentPdg.push_back(hasParent&&!multipleParents?gen[parent].pdgId():0);
      assocQuality.push_back(vg?gQuality[ig]:0);membership.push_back(vg?gMember[ig]:0);sourceCertain.push_back(vg?gSource[ig]:0);matchDR.push_back(rMatchDR[it]);matchRelPt.push_back(rMatchRelPt[it]);
      if (vh) members[ih].push_back(static_cast<int>(it));
    }
    auto trackOut=std::make_unique<nanoaod::FlatTable>(nr,"QTruthMatchedTrack",false,false);
    trackOut->addColumn<int>("trackIdx",trackIdx,"row in the shared unpacked reco track collection / track table");
    trackOut->addColumn<int>("genChargedIdx",genIdx,"QGenCharged row assigned by Stage 2, or -1");
    trackOut->addColumn<int>("hadronIdx",hadIdx,"QHadron row for a selected heavy decay ancestor, or -1");
    trackOut->addColumn<int>("quarkIdx",quarkIdx,"QFrag row from QGenCharged, or -1");
    trackOut->addColumn<int>("genPdgId",genPdg,"matched generator charged-particle PDG ID");
    trackOut->addColumn<int>("genCharge",genCharge,"matched generator charge in units of e");
    trackOut->addColumn<int>("matchStatus",matchStatus,"unchanged Stage-2 matching status");
    trackOut->addColumn<int>("nGenCandidates",nCandidates,"Stage-2 candidate multiplicity");
    trackOut->addColumn<float>("matchDR",matchDR,"Stage-2 reco-to-generator deltaR");trackOut->addColumn<float>("matchRelPt",matchRelPt,"Stage-2 relative-pT residual");
    trackOut->addColumn<int>("associationQuality",assocQuality,"Stage-1 retained ancestry quality");trackOut->addColumn<int>("membershipType",membership,"Stage-1 charged-particle membership category");trackOut->addColumn<int>("fragmentationSourceCertain",sourceCertain,"Stage-1 physical-source certainty; additional candidates remain uncertain");
    trackOut->addColumn<float>("genProdX",genX,"matched generator particle production x [cm]");
    trackOut->addColumn<float>("genProdY",genY,"matched generator particle production y [cm]");
    trackOut->addColumn<float>("genProdZ",genZ,"matched generator particle production z [cm]");
    trackOut->addColumn<int>("genProdVertexValid",genVertexValid,"1 when matched generator production position is finite");
    trackOut->addColumn<int>("daughterRelation",daughterRelation,"0 unassigned, 1 immediate daughter, 2 indirect descendant, 3 unresolved ancestry");
    trackOut->addColumn<int>("immediateParentPrunedIdx",immediateParent,"unique immediate pruned mother index, or -1");
    trackOut->addColumn<int>("immediateParentPdgId",immediateParentPdg,"PDG ID of the unique immediate mother, or 0 when unresolved");
    ev.put(std::move(trackOut),"QTruthMatchedTrack");

    std::vector<int> truthHadronIdx(nh),truthTrainHadronIdx(nh),truthHadronPdg(nh),truthHadronCategory(nh),truthNReco(nh),truthNPairs(nh),truthNCategory(nh);
    int groups=0;
    for(size_t ih=0;ih<nh;++ih){truthHadronIdx[ih]=ih;truthTrainHadronIdx[ih]=ih;truthHadronPdg[ih]=hPdg[ih];truthHadronCategory[ih]=hadronCategory(hPdg[ih],hOrigin[ih]);truthNReco[ih]=members[ih].size();truthNPairs[ih]=members[ih].size()>=2?members[ih].size()*(members[ih].size()-1)/2:0;truthNCategory[ih]=members[ih].empty()?0:members[ih].size()==1?1:members[ih].size()==2?2:3;if(members[ih].size()>=2)++groups;}
    auto hadOut=std::make_unique<nanoaod::FlatTable>(nh,"QTruthMatchedHadron",false,false);
    hadOut->addColumn<int>("hadronIdx",truthHadronIdx,"row in QHadron");hadOut->addColumn<int>("qTrainHadronIdx",truthTrainHadronIdx,"corresponding row in one-to-one QTrainHadron");
    hadOut->addColumn<int>("pdgId",truthHadronPdg,"signed selected generator hadron PDG ID");hadOut->addColumn<int>("hadronCategory",truthHadronCategory,"0 other, 1 B, 2 prompt D, 3 D-from-B, 4 D-origin-unknown");
    hadOut->addColumn<int>("nMatchedRecoTracks",truthNReco,"number of Stage-2 assigned reco tracks mapped to this QHadron; includes competing status=2 matches");hadOut->addColumn<int>("nCategory",truthNCategory,"0=N0, 1=N1, 2=N2, 3=N3+ for matched reco tracks");
    hadOut->addColumn<int>("nPossibleTrackPairs",truthNPairs,"nMatchedRecoTracks choose 2");ev.put(std::move(hadOut),"QTruthMatchedHadron");

    std::unordered_map<std::uint64_t,int> selected; selected.reserve(pairs.nRows()*2);
    for(size_t ip=0;ip<pairs.nRows();++ip) selected.emplace(pairKey(pA[ip],pB[ip]),static_cast<int>(ip));
    const bool selectedHasFits=pairs.columnIndex("vtxFitValid")>=0;
    if(transient.size()!=tracks.size()) throw cms::Exception("LogicError") << "Shared transient-track vector size mismatch";
    std::vector<int> t1,t2,hidx,qhidx,cat,nGroup,nGroupCategory,gen1,gen2,status1,status2,rel1,rel2,commonDirect,selectedPass,selectedRow,reused,reusedTruth,selectionReason;
    std::vector<float> pairDR,pairMass,pairDca;
    std::vector<TwoTrackVertexFitter::Result> results; std::vector<reco::TransientTrack> scratch; scratch.reserve(2);
    TwoTrackVertexFitter fitter; VertexDistance3D distance3D;
    int already=0,fitCalls=0,success=0,failed=0,covBad=0;
    std::unordered_map<std::uint64_t,TwoTrackVertexFitter::Result> truthFitCache;
    const reco::Vertex* pv=!vertices.empty()&&!vertices.front().isFake()?&vertices.front():nullptr;
    for(size_t ih=0;ih<nh;++ih) {
      auto& group=members[ih]; if(group.size()<2) continue;
      const int category=hadronCategory(hPdg[ih],hOrigin[ih]);
      for(size_t a=0;a<group.size();++a) for(size_t b=a+1;b<group.size();++b) {
        const int ia=group[a],ib=group[b],ta=rTrack[ia],tb=rTrack[ib],ga=rGen[ia],gb=rGen[ib];
        if(ta<0||tb<0||static_cast<size_t>(ta)>=tracks.size()||static_cast<size_t>(tb)>=tracks.size()) continue;
        const auto found=selected.find(pairKey(ta,tb)); const bool inSelected=found!=selected.end();
        const int relA=daughterRelation[ia],relB=daughterRelation[ib];
        TwoTrackVertexFitter::Result fit; bool reusedFit=false,reusedTruthFit=false;
        if(inSelected&&selectedHasFits) { readFit(pairs,found->second,fit); reusedFit=true; ++already; }
        else {
          if(inSelected) ++already;
          const auto cached=truthFitCache.find(pairKey(ta,tb));
          if(cached!=truthFitCache.end()){fit=cached->second;reusedTruthFit=true;}
          else {
            ++fitCalls;
            fit=fitter.fit(transient[ta],transient[tb],tracks[ta],tracks[tb],pv,scratch);
            truthFitCache.emplace(pairKey(ta,tb),fit);
            if(fit.fitValid) ++success; else ++failed; if(fit.fitValid&&!fit.covarianceValid) ++covBad;
          }
        }
        const double dr=reco::deltaR(tracks[ta].eta(),tracks[ta].phi(),tracks[tb].eta(),tracks[tb].phi());
        const double massHyp=0.13957018;
        const double e1=std::sqrt(tracks[ta].p()*tracks[ta].p()+massHyp*massHyp),e2=std::sqrt(tracks[tb].p()*tracks[tb].p()+massHyp*massHyp);
        const double sx=tracks[ta].px()+tracks[tb].px(),sy=tracks[ta].py()+tracks[tb].py(),sz=tracks[ta].pz()+tracks[tb].pz();
        const double mass=std::sqrt(std::max(0.,(e1+e2)*(e1+e2)-sx*sx-sy*sy-sz*sz));
        double dca=-1.; int reason=0;
        if(inSelected) reason=0;
        else if(!transient[ta].isValid()||!transient[tb].isValid()) reason=1;
        else if(dr>1.5) reason=2;
        else if(mass>5.) reason=3;
        else {
          TwoTrackMinimumDistance md;
          if(md.calculate(transient[ta].impactPointState(),transient[tb].impactPointState())) {
            auto dist=distance3D.distance(VertexState(md.points().second,transient[ta].impactPointState().cartesianError().position()),
                                           VertexState(md.points().first,transient[tb].impactPointState().cartesianError().position()));
            dca=dist.value(); if(dca>0.1) reason=4;
          } else reason=5; // historical selection stores invalid DCA (-1) rather than rejecting it
        }
        t1.push_back(ta);t2.push_back(tb);hidx.push_back(ih);qhidx.push_back(ih);cat.push_back(category);nGroup.push_back(group.size());nGroupCategory.push_back(group.size()==2?2:group.size()>2?3:1);
        gen1.push_back(ga);gen2.push_back(gb);status1.push_back(rStatus[ia]);status2.push_back(rStatus[ib]);rel1.push_back(relA);rel2.push_back(relB);
        const double daughterVertexSeparation = std::sqrt(std::pow(genX[ia]-genX[ib],2)+std::pow(genY[ia]-genY[ib],2)+std::pow(genZ[ia]-genZ[ib],2));
        bool matchesHadronDecay=true;
        if(hadDecayVertexAvailable){
          const double distanceA=std::sqrt(std::pow(genX[ia]-hDecayX[ih],2)+std::pow(genY[ia]-hDecayY[ih],2)+std::pow(genZ[ia]-hDecayZ[ih],2));
          const double distanceB=std::sqrt(std::pow(genX[ib]-hDecayX[ih],2)+std::pow(genY[ib]-hDecayY[ih],2)+std::pow(genZ[ib]-hDecayZ[ih],2));
          matchesHadronDecay=hDecayValid[ih]&&std::isfinite(distanceA)&&std::isfinite(distanceB)&&distanceA<=1.e-4&&distanceB<=1.e-4;
        }
        const bool commonVertex = relA==1 && relB==1 && genVertexValid[ia] && genVertexValid[ib] &&
                                  std::isfinite(daughterVertexSeparation) && daughterVertexSeparation<=1.e-4 && matchesHadronDecay;
        commonDirect.push_back(commonVertex);selectedPass.push_back(inSelected);selectedRow.push_back(inSelected?found->second:-1);reused.push_back(reusedFit);reusedTruth.push_back(reusedTruthFit);
        selectionReason.push_back(reason);pairDR.push_back(dr);pairMass.push_back(mass);pairDca.push_back(dca);results.push_back(fit);
      }
    }
    auto out=std::make_unique<nanoaod::FlatTable>(t1.size(),"QTruthMatchedPairFit",false,false);
#define ADD_I(n,v,d) out->addColumn<int>(n,v,d)
#define ADD_F(n,v,d) out->addColumn<float>(n,v,d)
    ADD_I("trackIdx1",t1,"shared event-local reconstructed track index"); ADD_I("trackIdx2",t2,"shared event-local reconstructed track index");
    ADD_I("hadronIdx",hidx,"QHadron row"); ADD_I("qTrainHadronIdx",qhidx,"same numeric row in one-to-one QTrainHadron table; explicit mapping");
    ADD_I("hadronCategory",cat,"0 other, 1 B, 2 prompt D, 3 D-from-B, 4 D-origin-unknown"); ADD_I("nMatchedTracks",nGroup,"assigned reconstructed tracks in this hadron group");ADD_I("nCategory",nGroupCategory,"1=N1, 2=N2, 3=N3+ group");
    ADD_I("genChargedIdx1",gen1,"QGenCharged endpoint row"); ADD_I("genChargedIdx2",gen2,"QGenCharged endpoint row");
    ADD_I("matchStatus1",status1,"Stage-2 endpoint status"); ADD_I("matchStatus2",status2,"Stage-2 endpoint status");
    ADD_I("daughterRelation1",rel1,"0 unassigned, 1 direct, 2 indirect, 3 unresolved"); ADD_I("daughterRelation2",rel2,"0 unassigned, 1 direct, 2 indirect, 3 unresolved");
    ADD_I("commonImmediateDecayVertex",commonDirect,"1 only when both have selected hadron as immediate mother and finite production positions agree with each other and, when available, the resolved hadron decay vertex within 1e-4 cm");
    ADD_I("passedOrdinaryPairSelection",selectedPass,"1 when endpoints appear in unchanged trackPair table"); ADD_I("trackPairIdx",selectedRow,"trackPair row, or -1");
    ADD_I("fitReusedFromTrackPair",reused,"1 when fit observables copied from selected-pair fit branches"); ADD_I("fitReusedFromEarlierTruthPair",reusedTruth,"1 when a shared physical pair fit was computed once for another hadron relationship"); ADD_I("ordinarySelectionFailure",selectionReason,"0 passed, 1 invalid transient track, 2 deltaR, 3 mass, 4 DCA, 5 DCA computation invalid but historical selection retains it");
    ADD_F("pairDeltaR",pairDR,"track deltaR"); ADD_F("pairInvMass",pairMass,"pion-mass-hypothesis pair mass [GeV]"); ADD_F("pairDca",pairDca,"minimum-distance DCA [cm], -1 if unavailable");
    std::vector<int> fitValid,fitStatus,covValid,pvValid,lxyValid,signedValid,l3dValid,p2Valid,p3Valid,perpValid;
    std::vector<float> chi2,ndof,chi2ndof,prob,x,y,z,cxx,cxy,cxz,cyy,cyz,czz,sxv,syv,szv,lxy,lxyerr,lxysig,slxy,slxysig,l3d,l3derr,l3dsig,cp2,cp3,perp,perpsig;
    for(const auto& r:results){fitValid.push_back(r.fitValid);fitStatus.push_back(r.status);covValid.push_back(r.covarianceValid);pvValid.push_back(r.pvValid);lxyValid.push_back(r.lxyValid);signedValid.push_back(r.signedLxyValid);l3dValid.push_back(r.l3dValid);p2Valid.push_back(r.pointing2DValid);p3Valid.push_back(r.pointing3DValid);perpValid.push_back(r.perpendicularValid);chi2.push_back(r.chi2);ndof.push_back(r.ndof);chi2ndof.push_back(r.chi2Ndof);prob.push_back(r.probability);x.push_back(r.x);y.push_back(r.y);z.push_back(r.z);cxx.push_back(r.covariance[0]);cxy.push_back(r.covariance[1]);cxz.push_back(r.covariance[2]);cyy.push_back(r.covariance[3]);cyz.push_back(r.covariance[4]);czz.push_back(r.covariance[5]);sxv.push_back(r.sigmaX);syv.push_back(r.sigmaY);szv.push_back(r.sigmaZ);lxy.push_back(r.lxy);lxyerr.push_back(r.lxyError);lxysig.push_back(r.lxySignificance);slxy.push_back(r.signedLxy);slxysig.push_back(r.signedLxySignificance);l3d.push_back(r.l3d);l3derr.push_back(r.l3dError);l3dsig.push_back(r.l3dSignificance);cp2.push_back(r.cosPointing2D);cp3.push_back(r.cosPointing3D);perp.push_back(r.perpendicular);perpsig.push_back(r.perpendicularSignificance);}
    ADD_I("vtxFitValid",fitValid,"Kalman two-track fit validity"); ADD_I("vtxFitStatus",fitStatus,"same status codes as trackPair"); ADD_I("vtxCovarianceValid",covValid,"positive-semidefinite covariance mask");
    ADD_F("vtxChi2",chi2,"fit chi2");ADD_F("vtxNdof",ndof,"fit ndof");ADD_F("vtxChi2Ndof",chi2ndof,"fit chi2/ndof");ADD_F("vtxProb",prob,"fit probability");ADD_F("vtxX",x,"fitted x [cm]");ADD_F("vtxY",y,"fitted y [cm]");ADD_F("vtxZ",z,"fitted z [cm]");
    ADD_F("vtxCovXX",cxx,"covariance xx [cm^2]");ADD_F("vtxCovXY",cxy,"covariance xy [cm^2]");ADD_F("vtxCovXZ",cxz,"covariance xz [cm^2]");ADD_F("vtxCovYY",cyy,"covariance yy [cm^2]");ADD_F("vtxCovYZ",cyz,"covariance yz [cm^2]");ADD_F("vtxCovZZ",czz,"covariance zz [cm^2]");
    ADD_F("vtxSigmaX",sxv,"fitted x uncertainty [cm]");ADD_F("vtxSigmaY",syv,"fitted y uncertainty [cm]");ADD_F("vtxSigmaZ",szv,"fitted z uncertainty [cm]");
    ADD_I("vtxPVValid",pvValid,"leading reconstructed PV validity");ADD_I("vtxLxyValid",lxyValid,"Lxy validity");ADD_F("vtxLxy",lxy,"PV-to-fit transverse displacement [cm]");ADD_F("vtxLxyErr",lxyerr,"Lxy uncertainty [cm]");ADD_F("vtxLxySig",lxysig,"Lxy significance");
    ADD_I("vtxSignedLxyValid",signedValid,"signed Lxy validity");ADD_F("vtxSignedLxy",slxy,"signed projected Lxy [cm]");ADD_F("vtxSignedLxySig",slxysig,"signed Lxy significance");ADD_I("vtxL3DValid",l3dValid,"L3D validity");ADD_F("vtxL3D",l3d,"PV-to-fit 3D displacement [cm]");ADD_F("vtxL3DErr",l3derr,"L3D uncertainty [cm]");ADD_F("vtxL3DSig",l3dsig,"L3D significance");
    ADD_I("vtxPointing2DValid",p2Valid,"2D pointing validity");ADD_F("vtxCosPointing2D",cp2,"2D pointing cosine");ADD_I("vtxPointing3DValid",p3Valid,"3D pointing validity");ADD_F("vtxCosPointing3D",cp3,"3D pointing cosine");ADD_I("vtxPerpendicularValid",perpValid,"perpendicular displacement validity");ADD_F("vtxPerpDisplacement",perp,"perpendicular displacement [cm]");ADD_F("vtxPerpDisplacementSig",perpsig,"perpendicular displacement significance");
#undef ADD_I
#undef ADD_F
    ev.put(std::move(out),"QTruthMatchedPairFit");
    auto summary=std::make_unique<nanoaod::FlatTable>(1,"QTruthPairFitSummary",true,false);
    std::vector<uint32_t> run{ev.id().run()}, luminosityBlock{ev.luminosityBlock()};
    std::vector<uint64_t> eventNumber{ev.id().event()};
    std::vector<int> nHad{static_cast<int>(nh)},nGroups{groups},nPairs{static_cast<int>(t1.size())},nSelected{already},nExtra{fitCalls},nSuccess{success},nFailed{failed},nCovBad{covBad};
    std::vector<float> elapsed{std::chrono::duration<float,std::milli>(Clock::now()-start).count()};
    std::vector<int> pvIsValid{pv?1:0};
    std::vector<float> pvx{pv?static_cast<float>(pv->x()):-999.f},pvy{pv?static_cast<float>(pv->y()):-999.f},pvz{pv?static_cast<float>(pv->z()):-999.f};
    std::vector<float> pvcxx{pv?static_cast<float>(pv->covariance(0,0)):-999.f},pvcxy{pv?static_cast<float>(pv->covariance(0,1)):-999.f},pvcxz{pv?static_cast<float>(pv->covariance(0,2)):-999.f};
    std::vector<float> pvcyy{pv?static_cast<float>(pv->covariance(1,1)):-999.f},pvcyz{pv?static_cast<float>(pv->covariance(1,2)):-999.f},pvczz{pv?static_cast<float>(pv->covariance(2,2)):-999.f};
    summary->addColumn<int>("nSelectedHeavyHadrons",nHad,"QHadron rows (all supported B/D categories)");summary->addColumn<int>("nHadronGroupsAtLeast2Tracks",nGroups,"hadron groups with at least two assigned reco tracks");summary->addColumn<int>("nTruthAssociatedPairs",nPairs,"hadron-specific unordered pair relationships");summary->addColumn<int>("nPairsAlreadyInTrackPair",nSelected,"truth-pair rows with an ordinary selected pair identity");summary->addColumn<int>("nPairsRequiringFit",nExtra,"truth-pair fits run by this producer");summary->addColumn<int>("nSupplementaryFitSuccess",nSuccess,"successful supplementary fits");summary->addColumn<int>("nSupplementaryFitFailure",nFailed,"failed supplementary fits");summary->addColumn<int>("nSupplementaryCovarianceInvalid",nCovBad,"successful supplementary fits with invalid covariance");summary->addColumn<float>("supplementaryProducerMs",elapsed,"per-event producer time [ms]");
    summary->addColumn<uint32_t>("run",run,"CMS run identifier");summary->addColumn<uint32_t>("luminosityBlock",luminosityBlock,"CMS luminosity-block identifier");summary->addColumn<uint64_t>("event",eventNumber,"CMS event identifier");
    summary->addColumn<int>("recoPVValid",pvIsValid,"1 when the same leading non-fake reconstructed PV used by TrackGVFeaturesProducer is available");
    summary->addColumn<float>("recoPVX",pvx,"reconstructed PV x [cm]");summary->addColumn<float>("recoPVY",pvy,"reconstructed PV y [cm]");summary->addColumn<float>("recoPVZ",pvz,"reconstructed PV z [cm]");
    summary->addColumn<float>("recoPVCovXX",pvcxx,"reconstructed PV covariance xx [cm^2]");summary->addColumn<float>("recoPVCovXY",pvcxy,"reconstructed PV covariance xy [cm^2]");summary->addColumn<float>("recoPVCovXZ",pvcxz,"reconstructed PV covariance xz [cm^2]");
    summary->addColumn<float>("recoPVCovYY",pvcyy,"reconstructed PV covariance yy [cm^2]");summary->addColumn<float>("recoPVCovYZ",pvcyz,"reconstructed PV covariance yz [cm^2]");summary->addColumn<float>("recoPVCovZZ",pvczz,"reconstructed PV covariance zz [cm^2]");
    ev.put(std::move(summary),"QTruthPairFitSummary");
  }
private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_; edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_; edm::EDGetTokenT<nanoaod::FlatTable> hadToken_,chargedToken_,recoToken_,pairToken_;
  edm::EDGetTokenT<std::vector<reco::TransientTrack>> transientToken_;
};

DEFINE_FWK_MODULE(TruthMatchedPairFitProducer);
