#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include <limits>
#include <cmath>

class SVTruthTableProducer : public edm::global::EDProducer<> {
public:
  explicit SVTruthTableProducer(const edm::ParameterSet& iConfig)
      : svTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("svTable"))),
        svTrkTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("svTrkTable"))),
        gvTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvTable"))),
        gvDauTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvDaughtersTable"))),
        nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
        dRMax_(iConfig.getParameter<double>("dR_max")),
        relPtMax_(iConfig.getParameter<double>("relPt_max")) {
    produces<nanoaod::FlatTable>("SVTruthTable");
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    edm::Handle<nanoaod::FlatTable> svTable, svTrkTable, gvTable, gvDauTable;
    iEvent.getByToken(svTableToken_, svTable);
    iEvent.getByToken(svTrkTableToken_, svTrkTable);
    iEvent.getByToken(gvTableToken_, gvTable);
    iEvent.getByToken(gvDauTableToken_, gvDauTable);

    const unsigned int nSV = svTable->size();
    const unsigned int nGV = gvTable->size();

    const int idxPdgClass = gvTable->columnIndex("pdgClass");
    const int idxPdgId = gvTable->columnIndex("Hadron_pdgId");
    const int idxIsB = gvTable->columnIndex("isB");
    const int idxIsD = gvTable->columnIndex("isD");

    const int idxSvTrkPt = svTrkTable->columnIndex("trk_pt");
    const int idxSvTrkEta = svTrkTable->columnIndex("trk_eta");
    const int idxSvTrkPhi = svTrkTable->columnIndex("trk_phi");
    const int idxSvTrkSVIdx = svTrkTable->columnIndex("trk_SVidx");

    const int idxDauPt = gvDauTable->columnIndex("pt");
    const int idxDauEta = gvDauTable->columnIndex("eta");
    const int idxDauPhi = gvDauTable->columnIndex("phi");
    const int idxDauGVIdx = gvDauTable->columnIndex("hadronIndex");

    if (idxPdgClass < 0 || idxPdgId < 0 || idxSvTrkPt < 0 || idxSvTrkEta < 0 || idxSvTrkPhi < 0 || idxSvTrkSVIdx < 0 ||
        idxDauPt < 0 || idxDauEta < 0 || idxDauPhi < 0 || idxDauGVIdx < 0) {
      throw cms::Exception("MissingColumn") << "SVTruthTableProducer is missing required columns from input FlatTables";
    }

    std::vector<int> svIdx(nSV, -1), truthPdgClass(nSV, 0), truthPdgId(nSV, 0), truthIsB(nSV, 0), truthIsD(nSV, 0), nMatchedGV(nSV, 0);
    std::vector<float> bestMatchScore(nSV, -1.f);

    // Build per-SV track lists
    std::vector<std::vector<float>> svTrkPt(nSV), svTrkEta(nSV), svTrkPhi(nSV);
    for (unsigned int i = 0; i < svTrkTable->size(); ++i) {
      const int isv = static_cast<int>(svTrkTable->getAnyValue(i, idxSvTrkSVIdx));
      if (isv < 0 || isv >= static_cast<int>(nSV))
        continue;
      svTrkPt[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkPt)));
      svTrkEta[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkEta)));
      svTrkPhi[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkPhi)));
    }

    // Build per-GV daughter lists
    std::vector<std::vector<float>> dauPt(nGV), dauEta(nGV), dauPhi(nGV);
    for (unsigned int i = 0; i < gvDauTable->size(); ++i) {
      const int igv = static_cast<int>(gvDauTable->getAnyValue(i, idxDauGVIdx));
      if (igv < 0 || igv >= static_cast<int>(nGV))
        continue;
      dauPt[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauPt)));
      dauEta[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauEta)));
      dauPhi[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauPhi)));
    }

    for (unsigned int isv = 0; isv < nSV; ++isv) {
      svIdx[isv] = static_cast<int>(isv);
      int bestGV = -1;
      int bestCommon = -1;
      float bestScore = -1.f;

      for (unsigned int igv = 0; igv < nGV; ++igv) {
        int nCommon = 0;
        float score = 0.f;
        for (size_t it = 0; it < svTrkPt[isv].size(); ++it) {
          for (size_t id = 0; id < dauPt[igv].size(); ++id) {
            const float dR = reco::deltaR(svTrkEta[isv][it], svTrkPhi[isv][it], dauEta[igv][id], dauPhi[igv][id]);
            const float relPt = std::abs(svTrkPt[isv][it] - dauPt[igv][id]) / std::max(dauPt[igv][id], 1e-6f);
            if (dR < dRMax_ && relPt < relPtMax_) {
              ++nCommon;
              score += (1.f - dR / dRMax_) + (1.f - relPt / relPtMax_);
              break;
            }
          }
        }

        if (nCommon >= nRequiredCommonTracks_) {
          ++nMatchedGV[isv];
          if (nCommon > bestCommon || (nCommon == bestCommon && (score > bestScore || (score == bestScore && static_cast<int>(igv) < bestGV)))) {
            bestCommon = nCommon;
            bestScore = score;
            bestGV = static_cast<int>(igv);
          }
        }
      }

      if (bestGV >= 0) {
        truthPdgClass[isv] = static_cast<int>(gvTable->getAnyValue(bestGV, idxPdgClass));
        truthPdgId[isv] = static_cast<int>(gvTable->getAnyValue(bestGV, idxPdgId));
        if (idxIsB >= 0)
          truthIsB[isv] = static_cast<int>(gvTable->getAnyValue(bestGV, idxIsB));
        if (idxIsD >= 0)
          truthIsD[isv] = static_cast<int>(gvTable->getAnyValue(bestGV, idxIsD));
        bestMatchScore[isv] = bestScore;
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nSV, "SVTruth", false, false);
    table->addColumn<int>("svIdx", svIdx, "SV row index");
    table->addColumn<int>("truth_pdgClass", truthPdgClass, "Primary truth class label from matched GV");
    table->addColumn<int>("truth_pdgId", truthPdgId, "PDG id of best matched truth hadron");
    table->addColumn<int>("truth_isB", truthIsB, "Best matched GV is B hadron");
    table->addColumn<int>("truth_isD", truthIsD, "Best matched GV is D hadron");
    table->addColumn<int>("nMatchedGV", nMatchedGV, "Number of GV rows matched to this SV");
    table->addColumn<float>("bestMatchScore", bestMatchScore, "Best match score (track-sharing quality)");
    iEvent.put(std::move(table), "SVTruthTable");
  }

private:
  edm::EDGetTokenT<nanoaod::FlatTable> svTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> svTrkTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> gvTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> gvDauTableToken_;
  int nRequiredCommonTracks_;
  double dRMax_;
  double relPtMax_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVTruthTableProducer);
