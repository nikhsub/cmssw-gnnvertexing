#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TLorentzVector.h"

class SVGraphFeatureProducer : public edm::global::EDProducer<> {
public:
  explicit SVGraphFeatureProducer(const edm::ParameterSet& iConfig);
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

private:
  edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> globalTrackIdxToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  double dlenSigMin_;
};

SVGraphFeatureProducer::SVGraphFeatureProducer(const edm::ParameterSet& iConfig)
    : svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
      globalTrackIdxToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("globalTrackIdxMap"))),
      ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")) {
  produces<nanoaod::FlatTable>("SVGraphVertex");
  produces<nanoaod::FlatTable>("SVGraphTrack");
  produces<nanoaod::FlatTable>("SVGraphEdge");
}

void SVGraphFeatureProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  auto svs = iEvent.getHandle(svToken_);
  auto pvs = iEvent.getHandle(pvToken_);
  auto tracks = iEvent.getHandle(tracksToken_);
  auto globalIdxMap = iEvent.getHandle(globalTrackIdxToken_);

  const auto& pv0 = pvs->front();
  const auto& ttBuilder = iSetup.getData(ttbToken_);
  VertexDistance3D vdist;

  std::vector<int> v_svIdx, v_nTracks;
  std::vector<float> v_x, v_y, v_z, v_dlen, v_dlenSig, v_chi2, v_ndof, v_pt, v_eta, v_phi, v_mass;

  std::vector<int> t_svIdx, t_globalIdx, t_localIdx, t_charge, t_validHits;
  std::vector<float> t_pt, t_eta, t_phi, t_p, t_dz, t_dzSig;

  std::vector<int> e_svIdx, e_srcGlobalIdx, e_dstGlobalIdx, e_srcLocalIdx, e_dstLocalIdx;
  std::vector<float> e_deltaR, e_invMass, e_dca, e_dcaSig, e_cpToPv, e_pvToPcaSrc, e_pvToPcaDst, e_dotprodSrc, e_dotprodDst, e_pairMom;

  int keptSvIdx = 0;
  for (const auto& sv : *svs) {
    Measurement1D dl = vdist.distance(
        pv0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
    if (!(dl.value() > 0. && dl.significance() > dlenSigMin_))
      continue;

    std::vector<std::pair<int, const reco::Track*>> svTracks;
    std::vector<int> sv_t_globalIdx, sv_t_localIdx, sv_t_charge, sv_t_validHits;
    std::vector<float> sv_t_pt, sv_t_eta, sv_t_phi, sv_t_p, sv_t_dz, sv_t_dzSig;
    int localIdx = 0;
    TLorentzVector svP4;
    for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
      const edm::RefToBase<reco::Track>& trkRefBase = *it;
      if (trkRefBase.isNull())
        continue;
      reco::TrackRef trkRef = trkRefBase.castTo<reco::TrackRef>();
      if (trkRef.isNull())
        continue;
      if (trkRef.id() != tracks.id())
        continue;

      const int globalIdx = (*globalIdxMap)[trkRef];
      const reco::Track* trk = trkRef.get();
      svTracks.emplace_back(globalIdx, trk);

      TLorentzVector trkP4;
      trkP4.SetPtEtaPhiM(trk->pt(), trk->eta(), trk->phi(), 0.13957039);
      svP4 += trkP4;

      sv_t_globalIdx.push_back(globalIdx);
      sv_t_localIdx.push_back(localIdx++);
      sv_t_pt.push_back(trk->pt());
      sv_t_eta.push_back(trk->eta());
      sv_t_phi.push_back(trk->phi());
      sv_t_p.push_back(trk->p());
      sv_t_charge.push_back(trk->charge());
      sv_t_validHits.push_back(trk->numberOfValidHits());
      const float dz = trk->dz(pv0.position());
      sv_t_dz.push_back(dz);
      sv_t_dzSig.push_back(dz / trk->dzError());
    }

    if (svTracks.size() < 2)
      continue;

    for (size_t i = 0; i < svTracks.size(); ++i) {
      t_svIdx.push_back(keptSvIdx);
      t_globalIdx.push_back(sv_t_globalIdx[i]);
      t_localIdx.push_back(sv_t_localIdx[i]);
      t_pt.push_back(sv_t_pt[i]);
      t_eta.push_back(sv_t_eta[i]);
      t_phi.push_back(sv_t_phi[i]);
      t_p.push_back(sv_t_p[i]);
      t_charge.push_back(sv_t_charge[i]);
      t_validHits.push_back(sv_t_validHits[i]);
      t_dz.push_back(sv_t_dz[i]);
      t_dzSig.push_back(sv_t_dzSig[i]);
    }

    v_svIdx.push_back(keptSvIdx);
    v_x.push_back(sv.x());
    v_y.push_back(sv.y());
    v_z.push_back(sv.z());
    v_dlen.push_back(dl.value());
    v_dlenSig.push_back(dl.significance());
    v_chi2.push_back(sv.chi2());
    v_ndof.push_back(sv.ndof());
    v_nTracks.push_back(static_cast<int>(svTracks.size()));
    v_pt.push_back(svP4.Pt());
    v_eta.push_back(svP4.Eta());
    v_phi.push_back(svP4.Phi());
    v_mass.push_back(svP4.M());

    for (size_t i = 0; i < svTracks.size(); ++i) {
      TLorentzVector p4i;
      p4i.SetPtEtaPhiM(svTracks[i].second->pt(), svTracks[i].second->eta(), svTracks[i].second->phi(), 0.13957039);
      for (size_t j = i + 1; j < svTracks.size(); ++j) {
        TLorentzVector p4j;
        p4j.SetPtEtaPhiM(svTracks[j].second->pt(), svTracks[j].second->eta(), svTracks[j].second->phi(), 0.13957039);
        e_svIdx.push_back(keptSvIdx);
        e_srcGlobalIdx.push_back(svTracks[i].first);
        e_dstGlobalIdx.push_back(svTracks[j].first);
        e_srcLocalIdx.push_back(static_cast<int>(i));
        e_dstLocalIdx.push_back(static_cast<int>(j));
        e_deltaR.push_back(p4i.DeltaR(p4j));
        e_invMass.push_back((p4i + p4j).M());

        const auto ttrkI = ttBuilder.build(*svTracks[i].second);
        const auto ttrkJ = ttBuilder.build(*svTracks[j].second);
        float dca_val = -1.f, dcaSig_val = -1.f, cpToPv_val = -1.f, pvToPcaSrc_val = -1.f, pvToPcaDst_val = -1.f;
        float dotprodSrc_val = -999.f, dotprodDst_val = -999.f, pairMom_val = -1.f;
        TwoTrackMinimumDistance minDist;
        if (ttrkI.isValid() && ttrkJ.isValid() &&
            minDist.calculate(ttrkI.impactPointState(), ttrkJ.impactPointState())) {
          VertexDistance3D distanceComputer;
          auto m = distanceComputer.distance(
              VertexState(minDist.points().second, ttrkI.impactPointState().cartesianError().position()),
              VertexState(minDist.points().first, ttrkJ.impactPointState().cartesianError().position()));
          dca_val = m.value();
          if (m.error() > 0.f) {
            dcaSig_val = m.value() / m.error();
          }
          GlobalPoint pvp(pv0.position().x(), pv0.position().y(), pv0.position().z());
          GlobalPoint cp(minDist.crossingPoint());
          GlobalPoint srcPCA = minDist.points().second;
          GlobalPoint dstPCA = minDist.points().first;
          cpToPv_val = (cp - pvp).mag();
          pvToPcaSrc_val = (srcPCA - pvp).mag();
          pvToPcaDst_val = (dstPCA - pvp).mag();
          dotprodSrc_val = (srcPCA - pvp).unit().dot(ttrkI.impactPointState().globalDirection().unit());
          dotprodDst_val = (dstPCA - pvp).unit().dot(ttrkJ.impactPointState().globalDirection().unit());
          GlobalVector pairMomentum((Basic3DVector<float>)(ttrkI.track().momentum() + ttrkJ.track().momentum()));
          pairMom_val = pairMomentum.mag();
        }
        e_dca.push_back(dca_val);
        e_dcaSig.push_back(dcaSig_val);
        e_cpToPv.push_back(cpToPv_val);
        e_pvToPcaSrc.push_back(pvToPcaSrc_val);
        e_pvToPcaDst.push_back(pvToPcaDst_val);
        e_dotprodSrc.push_back(dotprodSrc_val);
        e_dotprodDst.push_back(dotprodDst_val);
        e_pairMom.push_back(pairMom_val);
      }
    }

    ++keptSvIdx;
  }

  auto vtxTable = std::make_unique<nanoaod::FlatTable>(v_svIdx.size(), "SVGraphVertex", false);
  vtxTable->addColumn<int>("svIdx", v_svIdx, "SV index aligned to filtered SV iteration order");
  vtxTable->addColumn<float>("x", v_x, "SV x");
  vtxTable->addColumn<float>("y", v_y, "SV y");
  vtxTable->addColumn<float>("z", v_z, "SV z");
  vtxTable->addColumn<float>("dlen", v_dlen, "SV 3D flight distance");
  vtxTable->addColumn<float>("dlenSig", v_dlenSig, "SV 3D flight distance significance");
  vtxTable->addColumn<float>("chi2", v_chi2, "SV fit chi2");
  vtxTable->addColumn<float>("ndof", v_ndof, "SV fit ndof");
  vtxTable->addColumn<int>("nTracks", v_nTracks, "SV associated tracks");
  vtxTable->addColumn<float>("pt", v_pt, "SV pT");
  vtxTable->addColumn<float>("eta", v_eta, "SV eta");
  vtxTable->addColumn<float>("phi", v_phi, "SV phi");
  vtxTable->addColumn<float>("mass", v_mass, "SV mass");

  auto trkTable = std::make_unique<nanoaod::FlatTable>(t_svIdx.size(), "SVGraphTrack", false);
  trkTable->addColumn<int>("svIdx", t_svIdx, "Parent SV index");
  trkTable->addColumn<int>("trk_globalIdx", t_globalIdx, "Canonical global track index");
  trkTable->addColumn<int>("trk_localIdxInSV", t_localIdx, "Track index within SV");
  trkTable->addColumn<float>("trk_pt", t_pt, "Track pt");
  trkTable->addColumn<float>("trk_eta", t_eta, "Track eta");
  trkTable->addColumn<float>("trk_phi", t_phi, "Track phi");
  trkTable->addColumn<float>("trk_p", t_p, "Track p");
  trkTable->addColumn<int>("trk_charge", t_charge, "Track charge");
  trkTable->addColumn<int>("trk_nValidHits", t_validHits, "Track valid hits");
  trkTable->addColumn<float>("trk_dz", t_dz, "Track dz wrt PV0");
  trkTable->addColumn<float>("trk_dzSig", t_dzSig, "Track dz significance wrt PV0");

  auto edgeTable = std::make_unique<nanoaod::FlatTable>(e_svIdx.size(), "SVGraphEdge", false);
  edgeTable->addColumn<int>("svIdx", e_svIdx, "Parent SV index");
  edgeTable->addColumn<int>("edge_src_globalIdx", e_srcGlobalIdx, "Edge src global track index");
  edgeTable->addColumn<int>("edge_dst_globalIdx", e_dstGlobalIdx, "Edge dst global track index");
  edgeTable->addColumn<int>("edge_src_localIdxInSV", e_srcLocalIdx, "Edge src local index");
  edgeTable->addColumn<int>("edge_dst_localIdxInSV", e_dstLocalIdx, "Edge dst local index");
  edgeTable->addColumn<float>("edge_deltaR", e_deltaR, "Track pair deltaR");
  edgeTable->addColumn<float>("edge_invMass", e_invMass, "Track pair invariant mass");
  edgeTable->addColumn<float>("edge_dca", e_dca, "Track pair DCA");
  edgeTable->addColumn<float>("edge_dcaSig", e_dcaSig, "Track pair DCA significance");
  edgeTable->addColumn<float>("edge_cpToPv", e_cpToPv, "Distance PV to pair crossing point");
  edgeTable->addColumn<float>("edge_pvToPcaSrc", e_pvToPcaSrc, "Distance PV to src track PCA");
  edgeTable->addColumn<float>("edge_pvToPcaDst", e_pvToPcaDst, "Distance PV to dst track PCA");
  edgeTable->addColumn<float>("edge_dotprodSrc", e_dotprodSrc, "Directionality dot product for src track");
  edgeTable->addColumn<float>("edge_dotprodDst", e_dotprodDst, "Directionality dot product for dst track");
  edgeTable->addColumn<float>("edge_pairMom", e_pairMom, "Track pair momentum magnitude");

  iEvent.put(std::move(vtxTable), "SVGraphVertex");
  iEvent.put(std::move(trkTable), "SVGraphTrack");
  iEvent.put(std::move(edgeTable), "SVGraphEdge");
}

DEFINE_FWK_MODULE(SVGraphFeatureProducer);
