#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <limits>

class SVGraphFeatureProducer : public edm::global::EDProducer<> {
public:
  explicit SVGraphFeatureProducer(const edm::ParameterSet& iConfig);
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

private:
  struct GraphTrackInfo {
    int globalIdx = -1;
    const reco::Track* trk = nullptr;

    int isFitConstituent = 0;
    int isConeConstituent = 0;

    float dRToSVAxis = -1.f;

    float dzPV0 = 0.f;
    float dzSigPV0 = 0.f;
    float absDzPV0 = 0.f;

    float dxyPV0 = 0.f;
    float dxySigPV0 = 0.f;

    int closestPVIdx = -1;
    float dzClosestPV = 999.f;
    int isClosestPV0 = 0;
  };

  GraphTrackInfo makeGraphTrackInfo(const reco::Track* trk,
                                    int globalIdx,
                                    int isFitConstituent,
                                    int isConeConstituent,
                                    float dRToSVAxis,
                                    const std::vector<reco::Vertex>& pvs) const;

  edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> globalTrackIdxToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

  double dlenSigMin_;

  bool includeNearbyTracks_;
  double nearbyTrackDR_;
  double nearbyTrackPtMin_;
  double nearbyTrackEtaMax_;
  int maxExtraTracks_;
  bool requireExtraHighPurity_;
};

SVGraphFeatureProducer::SVGraphFeatureProducer(const edm::ParameterSet& iConfig)
    : svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
      globalTrackIdxToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("globalTrackIdxMap"))),
      ttbToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
      includeNearbyTracks_(iConfig.existsAs<bool>("includeNearbyTracks")
                               ? iConfig.getParameter<bool>("includeNearbyTracks")
                               : false),
      nearbyTrackDR_(iConfig.existsAs<double>("nearbyTrackDR")
                         ? iConfig.getParameter<double>("nearbyTrackDR")
                         : 0.3),
      nearbyTrackPtMin_(iConfig.existsAs<double>("nearbyTrackPtMin")
                            ? iConfig.getParameter<double>("nearbyTrackPtMin")
                            : 0.8),
      nearbyTrackEtaMax_(iConfig.existsAs<double>("nearbyTrackEtaMax")
                             ? iConfig.getParameter<double>("nearbyTrackEtaMax")
                             : 2.5),
      maxExtraTracks_(iConfig.existsAs<int>("maxExtraTracks")
                          ? iConfig.getParameter<int>("maxExtraTracks")
                          : 4),
      requireExtraHighPurity_(iConfig.existsAs<bool>("requireExtraHighPurity")
                                  ? iConfig.getParameter<bool>("requireExtraHighPurity")
                                  : false) {
  produces<nanoaod::FlatTable>("SVGraphVertex");
  produces<nanoaod::FlatTable>("SVGraphTrack");
  produces<nanoaod::FlatTable>("SVGraphEdge");
}

SVGraphFeatureProducer::GraphTrackInfo SVGraphFeatureProducer::makeGraphTrackInfo(
    const reco::Track* trk,
    int globalIdx,
    int isFitConstituent,
    int isConeConstituent,
    float dRToSVAxis,
    const std::vector<reco::Vertex>& pvs) const {
  GraphTrackInfo info;

  info.globalIdx = globalIdx;
  info.trk = trk;
  info.isFitConstituent = isFitConstituent;
  info.isConeConstituent = isConeConstituent;
  info.dRToSVAxis = dRToSVAxis;

  if (!trk || pvs.empty())
    return info;

  const auto& pv0 = pvs.front();

  const float dzPV0 = trk->dz(pv0.position());
  const float dzErr = trk->dzError();

  info.dzPV0 = dzPV0;
  info.absDzPV0 = std::abs(dzPV0);
  info.dzSigPV0 = (dzErr > 0.f ? dzPV0 / dzErr : 0.f);

  const float dxyPV0 = trk->dxy(pv0.position());
  const float dxyErr = trk->dxyError();

  info.dxyPV0 = dxyPV0;
  info.dxySigPV0 = (dxyErr > 0.f ? dxyPV0 / dxyErr : 0.f);

  int bestPV = -1;
  float bestAbsDz = 999.f;

  for (size_t ipv = 0; ipv < pvs.size(); ++ipv) {
    const float dz = trk->dz(pvs[ipv].position());
    const float absDz = std::abs(dz);

    if (absDz < bestAbsDz) {
      bestAbsDz = absDz;
      bestPV = static_cast<int>(ipv);
    }
  }

  info.closestPVIdx = bestPV;
  info.dzClosestPV = bestAbsDz;
  info.isClosestPV0 = (bestPV == 0 ? 1 : 0);

  return info;
}

void SVGraphFeatureProducer::produce(edm::StreamID,
                                     edm::Event& iEvent,
                                     const edm::EventSetup& iSetup) const {
  auto svs = iEvent.getHandle(svToken_);
  auto pvs = iEvent.getHandle(pvToken_);
  auto tracks = iEvent.getHandle(tracksToken_);
  auto globalIdxMap = iEvent.getHandle(globalTrackIdxToken_);

  if (!svs.isValid() || !pvs.isValid() || pvs->empty() ||
      !tracks.isValid() || !globalIdxMap.isValid()) {
    return;
  }

  const auto& pv0 = pvs->front();
  const auto& ttBuilder = iSetup.getData(ttbToken_);
  VertexDistance3D vdist;

  // Vertex table.
  // Physics and PV aggregate variables are fitted-SV / fit-track based.
  // nTracks is total graph tracks = fit + cone.
  std::vector<int> v_svIdx;
  std::vector<int> v_nTracks;
  std::vector<int> v_nGraphTracks;
  std::vector<int> v_nFitTracks;
  std::vector<int> v_nExtraTracks;

  std::vector<int> v_closestPVIdx;
  std::vector<int> v_isClosestPV0;
  std::vector<int> v_nDistinctClosestPVs;
  std::vector<int> v_majorityTrackPVIdx;

  std::vector<float> v_x, v_y, v_z;
  std::vector<float> v_dlen, v_dlenSig, v_chi2, v_ndof;
  std::vector<float> v_pt, v_eta, v_phi, v_mass;

  std::vector<float> v_dzPV0, v_absDzPV0, v_dxyPV0;
  std::vector<float> v_dzClosestPV;
  std::vector<float> v_fracTracksClosestPV0;
  std::vector<float> v_fracTracksClosestNonPV0;
  std::vector<float> v_fracTracksMajorityPV;

  // Track table: graph tracks = fit tracks + optional nearby cone tracks.
  std::vector<int> t_svIdx, t_globalIdx, t_localIdx, t_charge, t_validHits;
  std::vector<int> t_isFitConstituent, t_isConeConstituent;
  std::vector<int> t_closestPVIdx, t_isClosestPV0;

  std::vector<float> t_pt, t_eta, t_phi, t_p;
  std::vector<float> t_dz, t_dzSig;
  std::vector<float> t_absDzPV0, t_dxyPV0, t_dxySigPV0;
  std::vector<float> t_dzClosestPV, t_dRToSVAxis;

  // Edge table: complete graph over graph tracks.
  std::vector<int> e_svIdx, e_srcGlobalIdx, e_dstGlobalIdx, e_srcLocalIdx, e_dstLocalIdx;
  std::vector<float> e_deltaR, e_invMass, e_dca, e_dcaSig, e_cpToPv;
  std::vector<float> e_pvToPcaSrc, e_pvToPcaDst, e_dotprodSrc, e_dotprodDst, e_pairMom;

  int keptSvIdx = 0;

  for (const auto& sv : *svs) {
    Measurement1D dl = vdist.distance(
        pv0,
        VertexState(
            RecoVertex::convertPos(sv.position()),
            RecoVertex::convertError(sv.error())));

    if (!(dl.value() > 0. && dl.significance() > dlenSigMin_))
      continue;

    // Axis used only for selecting nearby tracks and storing per-track dR.
    TVector3 flightVec(
        sv.x() - pv0.x(),
        sv.y() - pv0.y(),
        sv.z() - pv0.z());

    float svAxisEta = 0.f;
    float svAxisPhi = 0.f;

    if (flightVec.Mag() > 0.f) {
      svAxisEta = flightVec.Eta();
      svAxisPhi = flightVec.Phi();
    }

    std::vector<GraphTrackInfo> fitTracks;
    std::vector<GraphTrackInfo> graphTracks;
    std::unordered_set<int> usedGlobalIdx;

    TLorentzVector fitP4;

    // -------------------------
    // 1. Add original SV-fit tracks.
    // -------------------------
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
      if (globalIdx < 0)
        continue;

      const reco::Track* trk = trkRef.get();

      const float dRToAxis = reco::deltaR(
          trk->eta(),
          trk->phi(),
          static_cast<double>(svAxisEta),
          static_cast<double>(svAxisPhi));

      GraphTrackInfo info = makeGraphTrackInfo(
          trk,
          globalIdx,
          1,
          0,
          dRToAxis,
          *pvs);

      fitTracks.push_back(info);
      graphTracks.push_back(info);
      usedGlobalIdx.insert(globalIdx);

      TLorentzVector trkP4;
      trkP4.SetPtEtaPhiM(trk->pt(), trk->eta(), trk->phi(), 0.13957039);
      fitP4 += trkP4;
    }

    if (fitTracks.size() < 2)
      continue;

    // Fallback axis for pathological flight vector.
    if (!(std::isfinite(svAxisEta) && std::isfinite(svAxisPhi)) && fitP4.Pt() > 0.f) {
      svAxisEta = fitP4.Eta();
      svAxisPhi = fitP4.Phi();
    }

    // -------------------------
    // 2. Add optional nearby cone tracks.
    // These only affect track/edge graph constituents and nTracks/nGraphTracks.
    // They do not modify vertex physics features or vertex PV aggregate variables.
    // -------------------------
    if (includeNearbyTracks_) {
      std::vector<std::pair<float, GraphTrackInfo>> extraCandidates;

      for (size_t itrk = 0; itrk < tracks->size(); ++itrk) {
        reco::TrackRef trkRef(tracks, itrk);
        if (trkRef.isNull())
          continue;

        const int globalIdx = (*globalIdxMap)[trkRef];
        if (globalIdx < 0)
          continue;

        if (usedGlobalIdx.count(globalIdx))
          continue;

        const reco::Track* trk = trkRef.get();

        if (trk->pt() < nearbyTrackPtMin_)
          continue;

        if (std::abs(trk->eta()) > nearbyTrackEtaMax_)
          continue;

        if (requireExtraHighPurity_ &&
            !trk->quality(reco::TrackBase::highPurity))
          continue;

        const float dRToAxis = reco::deltaR(
            trk->eta(),
            trk->phi(),
            static_cast<double>(svAxisEta),
            static_cast<double>(svAxisPhi));

        if (dRToAxis > nearbyTrackDR_)
          continue;

        GraphTrackInfo info = makeGraphTrackInfo(
            trk,
            globalIdx,
            0,
            1,
            dRToAxis,
            *pvs);

        extraCandidates.emplace_back(dRToAxis, info);
      }

      std::sort(
          extraCandidates.begin(),
          extraCandidates.end(),
          [](const auto& a, const auto& b) {
            return a.first < b.first;
          });

      int nAdded = 0;
      for (const auto& cand : extraCandidates) {
        if (maxExtraTracks_ >= 0 && nAdded >= maxExtraTracks_)
          break;

        graphTracks.push_back(cand.second);
        usedGlobalIdx.insert(cand.second.globalIdx);
        ++nAdded;
      }
    }

    if (graphTracks.size() < 2)
      continue;

    const int nFitTracks = static_cast<int>(fitTracks.size());
    const int nGraphTracks = static_cast<int>(graphTracks.size());
    const int nExtraTracks = nGraphTracks - nFitTracks;

    // -------------------------
    // 3. Vertex-level PV/PU aggregate variables.
    // These are computed from fitTracks only.
    // -------------------------
    int nFitClosestPV0 = 0;
    std::unordered_set<int> fitDistinctPVs;
    std::unordered_map<int, int> fitPVCounts;

    for (const auto& info : fitTracks) {
      if (info.closestPVIdx >= 0) {
        fitDistinctPVs.insert(info.closestPVIdx);
        ++fitPVCounts[info.closestPVIdx];
      }

      if (info.isClosestPV0)
        ++nFitClosestPV0;
    }

    int fitMajorityPVIdx = -1;
    int fitMajorityPVCount = 0;

    for (const auto& kv : fitPVCounts) {
      if (kv.second > fitMajorityPVCount ||
          (kv.second == fitMajorityPVCount &&
           (fitMajorityPVIdx < 0 || kv.first < fitMajorityPVIdx))) {
        fitMajorityPVIdx = kv.first;
        fitMajorityPVCount = kv.second;
      }
    }

    const float fracFitClosestPV0 =
        nFitTracks > 0
            ? static_cast<float>(nFitClosestPV0) / static_cast<float>(nFitTracks)
            : -1.f;

    const float fracFitClosestNonPV0 =
        nFitTracks > 0 ? 1.f - fracFitClosestPV0 : -1.f;

    const float fracFitMajorityPV =
        nFitTracks > 0
            ? static_cast<float>(fitMajorityPVCount) / static_cast<float>(nFitTracks)
            : -1.f;

    // SV-position PV compatibility.
    const float svDzPV0 = sv.z() - pv0.z();
    const float svAbsDzPV0 = std::abs(svDzPV0);
    const float svDxyPV0 = std::hypot(sv.x() - pv0.x(), sv.y() - pv0.y());

    int svClosestPVIdx = -1;
    float svDzClosestPV = 999.f;

    for (size_t ipv = 0; ipv < pvs->size(); ++ipv) {
      const float dz = sv.z() - (*pvs)[ipv].z();
      const float absDz = std::abs(dz);

      if (absDz < svDzClosestPV) {
        svDzClosestPV = absDz;
        svClosestPVIdx = static_cast<int>(ipv);
      }
    }

    const int svIsClosestPV0 = (svClosestPVIdx == 0 ? 1 : 0);

    // -------------------------
    // 4. Fill vertex table.
    // Physics features = fit SV only.
    // nTracks = total graph tracks.
    // PV aggregate features = fit tracks only.
    // -------------------------
    v_svIdx.push_back(keptSvIdx);

    v_x.push_back(sv.x());
    v_y.push_back(sv.y());
    v_z.push_back(sv.z());

    v_dlen.push_back(dl.value());
    v_dlenSig.push_back(dl.significance());
    v_chi2.push_back(sv.chi2());
    v_ndof.push_back(sv.ndof());

    v_nTracks.push_back(nGraphTracks);
    v_nGraphTracks.push_back(nGraphTracks);
    v_nFitTracks.push_back(nFitTracks);
    v_nExtraTracks.push_back(nExtraTracks);

    v_pt.push_back(fitP4.Pt());
    v_eta.push_back(fitP4.Eta());
    v_phi.push_back(fitP4.Phi());
    v_mass.push_back(fitP4.M());

    v_dzPV0.push_back(svDzPV0);
    v_absDzPV0.push_back(svAbsDzPV0);
    v_dxyPV0.push_back(svDxyPV0);

    v_closestPVIdx.push_back(svClosestPVIdx);
    v_dzClosestPV.push_back(svDzClosestPV);
    v_isClosestPV0.push_back(svIsClosestPV0);

    v_fracTracksClosestPV0.push_back(fracFitClosestPV0);
    v_fracTracksClosestNonPV0.push_back(fracFitClosestNonPV0);
    v_nDistinctClosestPVs.push_back(static_cast<int>(fitDistinctPVs.size()));
    v_majorityTrackPVIdx.push_back(fitMajorityPVIdx);
    v_fracTracksMajorityPV.push_back(fracFitMajorityPV);

    // -------------------------
    // 5. Fill track table: fit + cone tracks.
    // -------------------------
    for (size_t i = 0; i < graphTracks.size(); ++i) {
      const auto& info = graphTracks[i];
      const reco::Track* trk = info.trk;

      t_svIdx.push_back(keptSvIdx);
      t_globalIdx.push_back(info.globalIdx);
      t_localIdx.push_back(static_cast<int>(i));

      t_isFitConstituent.push_back(info.isFitConstituent);
      t_isConeConstituent.push_back(info.isConeConstituent);
      t_dRToSVAxis.push_back(info.dRToSVAxis);

      t_pt.push_back(trk->pt());
      t_eta.push_back(trk->eta());
      t_phi.push_back(trk->phi());
      t_p.push_back(trk->p());
      t_charge.push_back(trk->charge());
      t_validHits.push_back(trk->numberOfValidHits());

      t_dz.push_back(info.dzPV0);
      t_dzSig.push_back(info.dzSigPV0);
      t_absDzPV0.push_back(info.absDzPV0);

      t_dxyPV0.push_back(info.dxyPV0);
      t_dxySigPV0.push_back(info.dxySigPV0);

      t_closestPVIdx.push_back(info.closestPVIdx);
      t_dzClosestPV.push_back(info.dzClosestPV);
      t_isClosestPV0.push_back(info.isClosestPV0);
    }

    // -------------------------
    // 6. Build edge table: complete graph over fit + cone tracks.
    // -------------------------
    for (size_t i = 0; i < graphTracks.size(); ++i) {
      TLorentzVector p4i;
      p4i.SetPtEtaPhiM(
          graphTracks[i].trk->pt(),
          graphTracks[i].trk->eta(),
          graphTracks[i].trk->phi(),
          0.13957039);

      for (size_t j = i + 1; j < graphTracks.size(); ++j) {
        TLorentzVector p4j;
        p4j.SetPtEtaPhiM(
            graphTracks[j].trk->pt(),
            graphTracks[j].trk->eta(),
            graphTracks[j].trk->phi(),
            0.13957039);

        e_svIdx.push_back(keptSvIdx);
        e_srcGlobalIdx.push_back(graphTracks[i].globalIdx);
        e_dstGlobalIdx.push_back(graphTracks[j].globalIdx);
        e_srcLocalIdx.push_back(static_cast<int>(i));
        e_dstLocalIdx.push_back(static_cast<int>(j));

        e_deltaR.push_back(p4i.DeltaR(p4j));
        e_invMass.push_back((p4i + p4j).M());

        const auto ttrkI = ttBuilder.build(*graphTracks[i].trk);
        const auto ttrkJ = ttBuilder.build(*graphTracks[j].trk);

        float dca_val = -1.f;
        float dcaSig_val = -1.f;
        float cpToPv_val = -1.f;
        float pvToPcaSrc_val = -1.f;
        float pvToPcaDst_val = -1.f;
        float dotprodSrc_val = -999.f;
        float dotprodDst_val = -999.f;
        float pairMom_val = -1.f;

        TwoTrackMinimumDistance minDist;

        if (ttrkI.isValid() &&
            ttrkJ.isValid() &&
            minDist.calculate(ttrkI.impactPointState(), ttrkJ.impactPointState())) {
          VertexDistance3D distanceComputer;

          auto m = distanceComputer.distance(
              VertexState(
                  minDist.points().second,
                  ttrkI.impactPointState().cartesianError().position()),
              VertexState(
                  minDist.points().first,
                  ttrkJ.impactPointState().cartesianError().position()));

          dca_val = m.value();

          if (m.error() > 0.f)
            dcaSig_val = m.value() / m.error();

          GlobalPoint pvp(pv0.position().x(), pv0.position().y(), pv0.position().z());
          GlobalPoint cp(minDist.crossingPoint());
          GlobalPoint srcPCA = minDist.points().second;
          GlobalPoint dstPCA = minDist.points().first;

          cpToPv_val = (cp - pvp).mag();
          pvToPcaSrc_val = (srcPCA - pvp).mag();
          pvToPcaDst_val = (dstPCA - pvp).mag();

          dotprodSrc_val =
              (srcPCA - pvp).unit().dot(ttrkI.impactPointState().globalDirection().unit());

          dotprodDst_val =
              (dstPCA - pvp).unit().dot(ttrkJ.impactPointState().globalDirection().unit());

          GlobalVector pairMomentum(
              (Basic3DVector<float>)(ttrkI.track().momentum() + ttrkJ.track().momentum()));

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

  vtxTable->addColumn<int>("nTracks", v_nTracks, "Total graph tracks: fit tracks plus nearby cone tracks");
  vtxTable->addColumn<int>("nGraphTracks", v_nGraphTracks, "Total graph tracks: fit tracks plus nearby cone tracks");
  vtxTable->addColumn<int>("nFitTracks", v_nFitTracks, "Original SV-fit constituent tracks");
  vtxTable->addColumn<int>("nExtraTracks", v_nExtraTracks, "Added nearby cone tracks");

  vtxTable->addColumn<float>("pt", v_pt, "SV pT from original fit tracks only");
  vtxTable->addColumn<float>("eta", v_eta, "SV eta from original fit tracks only");
  vtxTable->addColumn<float>("phi", v_phi, "SV phi from original fit tracks only");
  vtxTable->addColumn<float>("mass", v_mass, "SV mass from original fit tracks only");

  vtxTable->addColumn<float>("dzPV0", v_dzPV0, "SV z - PV0 z");
  vtxTable->addColumn<float>("absDzPV0", v_absDzPV0, "abs(SV z - PV0 z)");
  vtxTable->addColumn<float>("dxyPV0", v_dxyPV0, "SV transverse distance from PV0");

  vtxTable->addColumn<int>("closestPVIdx", v_closestPVIdx, "Reco PV closest in z to fitted SV");
  vtxTable->addColumn<float>("dzClosestPV", v_dzClosestPV, "Minimum abs(SV z - PV z)");
  vtxTable->addColumn<int>("isClosestPV0", v_isClosestPV0, "Fitted SV is closest in z to PV0");

  vtxTable->addColumn<float>("fracTracksClosestPV0", v_fracTracksClosestPV0,
                             "Fraction of original fit tracks closest in z to PV0");
  vtxTable->addColumn<float>("fracTracksClosestNonPV0", v_fracTracksClosestNonPV0,
                             "Fraction of original fit tracks closest in z to non-PV0");
  vtxTable->addColumn<int>("nDistinctClosestPVs", v_nDistinctClosestPVs,
                           "Number of distinct closest PV indices among original fit tracks");
  vtxTable->addColumn<int>("majorityTrackPVIdx", v_majorityTrackPVIdx,
                           "Most common closest PV index among original fit tracks");
  vtxTable->addColumn<float>("fracTracksMajorityPV", v_fracTracksMajorityPV,
                             "Fraction of original fit tracks assigned to majority closest PV");

  auto trkTable = std::make_unique<nanoaod::FlatTable>(t_svIdx.size(), "SVGraphTrack", false);
  trkTable->addColumn<int>("svIdx", t_svIdx, "Parent SV index");
  trkTable->addColumn<int>("trk_globalIdx", t_globalIdx, "Canonical global track index");
  trkTable->addColumn<int>("trk_localIdxInSV", t_localIdx, "Track index within graph constituents");

  trkTable->addColumn<int>("trk_isFitConstituent", t_isFitConstituent, "Track was used in SV fit");
  trkTable->addColumn<int>("trk_isConeConstituent", t_isConeConstituent, "Track added as nearby cone constituent");
  trkTable->addColumn<float>("trk_dRToSVAxis", t_dRToSVAxis, "Track deltaR to SV flight axis");

  trkTable->addColumn<float>("trk_pt", t_pt, "Track pt");
  trkTable->addColumn<float>("trk_eta", t_eta, "Track eta");
  trkTable->addColumn<float>("trk_phi", t_phi, "Track phi");
  trkTable->addColumn<float>("trk_p", t_p, "Track p");
  trkTable->addColumn<int>("trk_charge", t_charge, "Track charge");
  trkTable->addColumn<int>("trk_nValidHits", t_validHits, "Track valid hits");

  trkTable->addColumn<float>("trk_dz", t_dz, "Track dz wrt PV0");
  trkTable->addColumn<float>("trk_dzSig", t_dzSig, "Track dz significance wrt PV0");
  trkTable->addColumn<float>("trk_absDzPV0", t_absDzPV0, "abs(track dz wrt PV0)");

  trkTable->addColumn<float>("trk_dxyPV0", t_dxyPV0, "Track dxy wrt PV0");
  trkTable->addColumn<float>("trk_dxySigPV0", t_dxySigPV0, "Track dxy significance wrt PV0");

  trkTable->addColumn<int>("trk_closestPVIdx", t_closestPVIdx, "Closest reco PV index in track dz");
  trkTable->addColumn<float>("trk_dzClosestPV", t_dzClosestPV, "Minimum abs track dz to any reco PV");
  trkTable->addColumn<int>("trk_isClosestPV0", t_isClosestPV0, "Track closest in dz to PV0");

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
