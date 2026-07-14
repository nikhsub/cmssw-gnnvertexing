#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <algorithm>

typedef reco::Vertex::CovarianceMatrix CovMatrix;

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

    // -------------------------
    // SV position + covariance columns
    // -------------------------
    const int idxSVx = svTable->columnIndex("x");
    const int idxSVy = svTable->columnIndex("y");
    const int idxSVz = svTable->columnIndex("z");

    const int idxSVCovXX = svTable->columnIndex("covXX");
    const int idxSVCovXY = svTable->columnIndex("covXY");
    const int idxSVCovXZ = svTable->columnIndex("covXZ");
    const int idxSVCovYY = svTable->columnIndex("covYY");
    const int idxSVCovYZ = svTable->columnIndex("covYZ");
    const int idxSVCovZZ = svTable->columnIndex("covZZ");

    // -------------------------
    // SV track columns
    // -------------------------
    const int idxSvTrkPt = svTrkTable->columnIndex("trk_pt");
    const int idxSvTrkEta = svTrkTable->columnIndex("trk_eta");
    const int idxSvTrkPhi = svTrkTable->columnIndex("trk_phi");
    const int idxSvTrkSVIdx = svTrkTable->columnIndex("trk_SVidx");

    // -------------------------
    // GV columns
    // -------------------------
    const int idxGVx = gvTable->columnIndex("x");
    const int idxGVy = gvTable->columnIndex("y");
    const int idxGVz = gvTable->columnIndex("z");

    const int idxPdgClass = gvTable->columnIndex("pdgClass");
    const int idxPdgId = gvTable->columnIndex("Hadron_pdgId");
    const int idxIsB = gvTable->columnIndex("isB");
    const int idxIsD = gvTable->columnIndex("isD");
    const int idxIsBtoD = gvTable->columnIndex("isBtoD");

    const int idxGVNDauNoRecognizedSecondaryAncestor =
        gvTable->columnIndex("nDauNoRecognizedSecondaryAncestor");
    const int idxGVNDauFromB = gvTable->columnIndex("nDauFromB");
    const int idxGVNDauFromBC = gvTable->columnIndex("nDauFromBC");
    const int idxGVNDauFromC = gvTable->columnIndex("nDauFromC");
    const int idxGVNDauOtherSecondary = gvTable->columnIndex("nDauOtherSecondary");
    const int idxGVNDauOriginUnknown = gvTable->columnIndex("nDauOriginUnknown");

    // -------------------------
    // GV daughter columns
    // -------------------------
    const int idxDauPt = gvDauTable->columnIndex("pt");
    const int idxDauEta = gvDauTable->columnIndex("eta");
    const int idxDauPhi = gvDauTable->columnIndex("phi");
    const int idxDauGVIdx = gvDauTable->columnIndex("hadronIndex");
    const int idxDauOriginLabel = gvDauTable->columnIndex("originLabel");

    if (idxSVx < 0 || idxSVy < 0 || idxSVz < 0 ||
        idxSVCovXX < 0 || idxSVCovXY < 0 || idxSVCovXZ < 0 ||
        idxSVCovYY < 0 || idxSVCovYZ < 0 || idxSVCovZZ < 0 ||
        idxGVx < 0 || idxGVy < 0 || idxGVz < 0 ||
        idxPdgClass < 0 || idxPdgId < 0 ||
        idxSvTrkPt < 0 || idxSvTrkEta < 0 || idxSvTrkPhi < 0 || idxSvTrkSVIdx < 0 ||
        idxDauPt < 0 || idxDauEta < 0 || idxDauPhi < 0 || idxDauGVIdx < 0) {
      throw cms::Exception("MissingColumn")
          << "SVTruthTableProducer is missing required columns from input FlatTables. "
          << "Required: SV x/y/z/covXX/covXY/covXZ/covYY/covYZ/covZZ, "
          << "SV track trk_pt/trk_eta/trk_phi/trk_SVidx, "
          << "GV x/y/z/pdgClass/Hadron_pdgId, "
          << "and GVDaughters pt/eta/phi/hadronIndex.";
    }

    const bool hasIsB = idxIsB >= 0;
    const bool hasIsD = idxIsD >= 0;
    const bool hasIsBtoD = idxIsBtoD >= 0;
    const bool hasDaughterOrigin = idxDauOriginLabel >= 0;

    // -------------------------
    // Output vectors
    // -------------------------
    std::vector<int> svIdx(nSV, -1);

    std::vector<int> truthPdgClass(nSV, 0);
    std::vector<int> truthPdgId(nSV, 0);
    std::vector<int> truthIsB(nSV, 0);
    std::vector<int> truthIsD(nSV, 0);
    std::vector<int> truthIsBtoD(nSV, 0);

    std::vector<int> matchedGVIdx(nSV, -1);
    std::vector<int> nMatchedGV(nSV, 0);
    std::vector<int> nCommonTracks(nSV, 0);

    std::vector<float> bestMatchScore(nSV, -1.f);
    std::vector<float> bestMatchDistanceSig(nSV, -1.f);

    std::vector<int> nMatchedDaughters(nSV, 0);

    // originLabel semantics:
    // 0 = noRecognizedSecondaryAncestor / primary-or-pileup-like
    // 2 = fromB
    // 3 = fromBC
    // 4 = fromC
    // 5 = otherSecondary
    // 9 = unknown
    std::vector<int> nOriginNoRecognizedSecondaryAncestor(nSV, 0);
    std::vector<int> nOriginFromB(nSV, 0);
    std::vector<int> nOriginFromBC(nSV, 0);
    std::vector<int> nOriginFromC(nSV, 0);
    std::vector<int> nOriginOtherSecondary(nSV, 0);
    std::vector<int> nOriginUnknown(nSV, 0);

    std::vector<float> fracOriginNoRecognizedSecondaryAncestor(nSV, -1.f);
    std::vector<float> fracOriginHF(nSV, -1.f);
    std::vector<int> dominantOriginLabel(nSV, 9);

    // Matched-GV copied daughter-origin summaries.
    std::vector<int> matchedGV_nDauNoRecognizedSecondaryAncestor(nSV, -1);
    std::vector<int> matchedGV_nDauFromB(nSV, -1);
    std::vector<int> matchedGV_nDauFromBC(nSV, -1);
    std::vector<int> matchedGV_nDauFromC(nSV, -1);
    std::vector<int> matchedGV_nDauOtherSecondary(nSV, -1);
    std::vector<int> matchedGV_nDauOriginUnknown(nSV, -1);

    // -------------------------
    // Read SV positions and covariance
    // -------------------------
    std::vector<float> svX(nSV), svY(nSV), svZ(nSV);
    std::vector<CovMatrix> svCov(nSV);

    for (unsigned int isv = 0; isv < nSV; ++isv) {
      svIdx[isv] = static_cast<int>(isv);

      svX[isv] = static_cast<float>(svTable->getAnyValue(isv, idxSVx));
      svY[isv] = static_cast<float>(svTable->getAnyValue(isv, idxSVy));
      svZ[isv] = static_cast<float>(svTable->getAnyValue(isv, idxSVz));

      CovMatrix cov;
      cov(0, 0) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovXX));
      cov(0, 1) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovXY));
      cov(1, 0) = cov(0, 1);

      cov(0, 2) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovXZ));
      cov(2, 0) = cov(0, 2);

      cov(1, 1) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovYY));

      cov(1, 2) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovYZ));
      cov(2, 1) = cov(1, 2);

      cov(2, 2) = static_cast<float>(svTable->getAnyValue(isv, idxSVCovZZ));

      svCov[isv] = cov;
    }

    // -------------------------
    // Read GV positions
    // -------------------------
    std::vector<float> gvX(nGV), gvY(nGV), gvZ(nGV);

    for (unsigned int igv = 0; igv < nGV; ++igv) {
      gvX[igv] = static_cast<float>(gvTable->getAnyValue(igv, idxGVx));
      gvY[igv] = static_cast<float>(gvTable->getAnyValue(igv, idxGVy));
      gvZ[igv] = static_cast<float>(gvTable->getAnyValue(igv, idxGVz));
    }

    // -------------------------
    // Build per-SV track lists
    // -------------------------
    std::vector<std::vector<float>> svTrkPt(nSV), svTrkEta(nSV), svTrkPhi(nSV);

    for (unsigned int i = 0; i < svTrkTable->size(); ++i) {
      const int isv = static_cast<int>(svTrkTable->getAnyValue(i, idxSvTrkSVIdx));
      if (isv < 0 || isv >= static_cast<int>(nSV))
        continue;

      svTrkPt[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkPt)));
      svTrkEta[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkEta)));
      svTrkPhi[isv].push_back(static_cast<float>(svTrkTable->getAnyValue(i, idxSvTrkPhi)));
    }

    // -------------------------
    // Build per-GV daughter lists
    // -------------------------
    std::vector<std::vector<float>> dauPt(nGV), dauEta(nGV), dauPhi(nGV);
    std::vector<std::vector<int>> dauOriginLabel(nGV);

    for (unsigned int i = 0; i < gvDauTable->size(); ++i) {
      const int igv = static_cast<int>(gvDauTable->getAnyValue(i, idxDauGVIdx));
      if (igv < 0 || igv >= static_cast<int>(nGV))
        continue;

      dauPt[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauPt)));
      dauEta[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauEta)));
      dauPhi[igv].push_back(static_cast<float>(gvDauTable->getAnyValue(i, idxDauPhi)));

      if (hasDaughterOrigin)
        dauOriginLabel[igv].push_back(static_cast<int>(gvDauTable->getAnyValue(i, idxDauOriginLabel)));
      else
        dauOriginLabel[igv].push_back(9);
    }

    // -------------------------
    // Compute SV-GV distance matrix:
    // distanceSig = sqrt(dx^T CovSV^{-1} dx)
    // same convention as GenVertexProducer.
    // -------------------------
    std::vector<std::vector<float>> distances(nSV, std::vector<float>(nGV, 999.0f));
    std::vector<std::vector<float>> distancesOriginal(nSV, std::vector<float>(nGV, 999.0f));

    for (unsigned int isv = 0; isv < nSV; ++isv) {
      CovMatrix covInv = svCov[isv];
      covInv.Invert();

      for (unsigned int igv = 0; igv < nGV; ++igv) {
        const float dx = svX[isv] - gvX[igv];
        const float dy = svY[isv] - gvY[igv];
        const float dz = svZ[isv] - gvZ[igv];

        const float chi2 =
            dx * (covInv(0, 0) * dx + covInv(0, 1) * dy + covInv(0, 2) * dz) +
            dy * (covInv(1, 0) * dx + covInv(1, 1) * dy + covInv(1, 2) * dz) +
            dz * (covInv(2, 0) * dx + covInv(2, 1) * dy + covInv(2, 2) * dz);

        float dist = 999.0f;
        if (std::isfinite(chi2) && chi2 >= 0.f)
          dist = std::sqrt(chi2);

        distances[isv][igv] = dist;
        distancesOriginal[isv][igv] = dist;
      }
    }

    // -------------------------
    // Greedy matching, same basic scheme as GenVertexProducer.
    // No double matching in SVTruthTableProducer.
    // -------------------------
    while (true) {
      float minDist = 999.0f;
      int bestSV = -1;
      int bestGV = -1;

      for (unsigned int isv = 0; isv < nSV; ++isv) {
        for (unsigned int igv = 0; igv < nGV; ++igv) {
          if (distances[isv][igv] < minDist) {
            minDist = distances[isv][igv];
            bestSV = static_cast<int>(isv);
            bestGV = static_cast<int>(igv);
          }
        }
      }

      if (minDist >= 997.0f || bestSV < 0 || bestGV < 0)
        break;

      int common = 0;
      float trackScore = 0.f;
      std::vector<int> matchedOriginLabels;

      for (size_t iSV = 0; iSV < svTrkPt[bestSV].size(); ++iSV) {
        for (size_t iHad = 0; iHad < dauPt[bestGV].size(); ++iHad) {
          const float dR = reco::deltaR(
              svTrkEta[bestSV][iSV],
              svTrkPhi[bestSV][iSV],
              dauEta[bestGV][iHad],
              dauPhi[bestGV][iHad]);

          const float relPt =
              std::abs(svTrkPt[bestSV][iSV] - dauPt[bestGV][iHad]) /
              std::max(dauPt[bestGV][iHad], 1e-6f);

          if (dR < dRMax_ && relPt < relPtMax_) {
            ++common;
            trackScore += (1.f - dR / dRMax_) + (1.f - relPt / relPtMax_);

            if (iHad < dauOriginLabel[bestGV].size())
              matchedOriginLabels.push_back(dauOriginLabel[bestGV][iHad]);
            else
              matchedOriginLabels.push_back(9);

            if (common >= nRequiredCommonTracks_)
              break;
          }
        }

        if (common >= nRequiredCommonTracks_)
          break;
      }

      if (common >= nRequiredCommonTracks_) {
        matchedGVIdx[bestSV] = bestGV;
        nMatchedGV[bestSV] = 1;
        nCommonTracks[bestSV] = common;
        bestMatchScore[bestSV] = trackScore;
        bestMatchDistanceSig[bestSV] = minDist;

        truthPdgClass[bestSV] =
            static_cast<int>(gvTable->getAnyValue(bestGV, idxPdgClass));

        truthPdgId[bestSV] =
            static_cast<int>(gvTable->getAnyValue(bestGV, idxPdgId));

        if (hasIsB)
          truthIsB[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxIsB));

        if (hasIsD)
          truthIsD[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxIsD));

        if (hasIsBtoD)
          truthIsBtoD[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxIsBtoD));

        if (idxGVNDauNoRecognizedSecondaryAncestor >= 0)
          matchedGV_nDauNoRecognizedSecondaryAncestor[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauNoRecognizedSecondaryAncestor));

        if (idxGVNDauFromB >= 0)
          matchedGV_nDauFromB[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauFromB));

        if (idxGVNDauFromBC >= 0)
          matchedGV_nDauFromBC[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauFromBC));

        if (idxGVNDauFromC >= 0)
          matchedGV_nDauFromC[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauFromC));

        if (idxGVNDauOtherSecondary >= 0)
          matchedGV_nDauOtherSecondary[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauOtherSecondary));

        if (idxGVNDauOriginUnknown >= 0)
          matchedGV_nDauOriginUnknown[bestSV] =
              static_cast<int>(gvTable->getAnyValue(bestGV, idxGVNDauOriginUnknown));

        nMatchedDaughters[bestSV] = static_cast<int>(matchedOriginLabels.size());

        for (const int lab : matchedOriginLabels) {
          if (lab == 0)
            ++nOriginNoRecognizedSecondaryAncestor[bestSV];
          else if (lab == 2)
            ++nOriginFromB[bestSV];
          else if (lab == 3)
            ++nOriginFromBC[bestSV];
          else if (lab == 4)
            ++nOriginFromC[bestSV];
          else if (lab == 5)
            ++nOriginOtherSecondary[bestSV];
          else
            ++nOriginUnknown[bestSV];
        }

        const int nLab = nMatchedDaughters[bestSV];
        if (nLab > 0) {
          fracOriginNoRecognizedSecondaryAncestor[bestSV] =
              static_cast<float>(nOriginNoRecognizedSecondaryAncestor[bestSV]) /
              static_cast<float>(nLab);

          const int nHF =
              nOriginFromB[bestSV] +
              nOriginFromBC[bestSV] +
              nOriginFromC[bestSV];

          fracOriginHF[bestSV] =
              static_cast<float>(nHF) / static_cast<float>(nLab);

          int bestLab = 9;
          int bestCnt = -1;

          const std::vector<std::pair<int, int>> labCounts = {
              {0, nOriginNoRecognizedSecondaryAncestor[bestSV]},
              {2, nOriginFromB[bestSV]},
              {3, nOriginFromBC[bestSV]},
              {4, nOriginFromC[bestSV]},
              {5, nOriginOtherSecondary[bestSV]},
              {9, nOriginUnknown[bestSV]}};

          for (const auto& lc : labCounts) {
            if (lc.second > bestCnt) {
              bestCnt = lc.second;
              bestLab = lc.first;
            }
          }

          dominantOriginLabel[bestSV] = bestLab;
        }

        // Remove matched SV row and matched GV column.
        for (unsigned int igv = 0; igv < nGV; ++igv)
          distances[bestSV][igv] = 1000.0f;

        for (unsigned int isv = 0; isv < nSV; ++isv)
          distances[isv][bestGV] = 1000.0f;

      } else {
        // Pair spatially close but failed track-daughter matching.
        distances[bestSV][bestGV] = 998.0f;
      }
    }

    // -------------------------
    // Output table
    // -------------------------
    auto table = std::make_unique<nanoaod::FlatTable>(nSV, "SVTruth", false, false);

    table->addColumn<int>("svIdx", svIdx, "SV row index");
    table->addColumn<int>("truth_pdgClass", truthPdgClass, "Primary truth class label from matched GV");
    table->addColumn<int>("truth_pdgId", truthPdgId, "PDG id of best matched truth hadron");
    table->addColumn<int>("truth_isB", truthIsB, "Best matched GV is B hadron");
    table->addColumn<int>("truth_isD", truthIsD, "Best matched GV is D hadron");
    table->addColumn<int>("truth_isBtoD", truthIsBtoD, "Best matched GV is a D hadron with B-hadron ancestry");

    table->addColumn<int>("matchedGVIdx", matchedGVIdx, "GV row index matched to this SV");
    table->addColumn<int>("nMatchedGV", nMatchedGV, "One if this SV was greedily matched to a GV, else zero");
    table->addColumn<int>("nCommonTracks", nCommonTracks, "Number of common tracks used to accept the SV-GV match");

    table->addColumn<float>("bestMatchScore", bestMatchScore, "Track-sharing score for accepted SV-GV match");
    table->addColumn<float>("bestMatchDistanceSig", bestMatchDistanceSig, "SV-GV covariance-weighted distance significance for accepted match");

    table->addColumn<int>("nMatchedDaughters", nMatchedDaughters, "Number of matched GV daughters contributing to accepted SV-GV match");

    table->addColumn<int>(
        "nOriginNoRecognizedSecondaryAncestor",
        nOriginNoRecognizedSecondaryAncestor,
        "Matched daughters with originLabel 0; no recognized B/C/strange/tau/conversion ancestry; may include primary-like or pileup-like particles");

    table->addColumn<int>("nOriginFromB", nOriginFromB, "Matched daughters from B ancestry");
    table->addColumn<int>("nOriginFromBC", nOriginFromBC, "Matched daughters from charm with B ancestry");
    table->addColumn<int>("nOriginFromC", nOriginFromC, "Matched daughters from charm ancestry");
    table->addColumn<int>("nOriginOtherSecondary", nOriginOtherSecondary, "Matched daughters from strange/tau/conversion-like ancestry");
    table->addColumn<int>("nOriginUnknown", nOriginUnknown, "Matched daughters with unavailable or unknown origin label");

    table->addColumn<float>(
        "fracOriginNoRecognizedSecondaryAncestor",
        fracOriginNoRecognizedSecondaryAncestor,
        "Fraction of matched daughters with no recognized secondary ancestry; not a PV/PU discriminator");

    table->addColumn<float>("fracOriginHF", fracOriginHF, "Fraction of matched daughters with B/BC/C origin labels");
    table->addColumn<int>("dominantOriginLabel", dominantOriginLabel, "Most frequent matched-daughter origin label");

    table->addColumn<int>(
        "matchedGV_nDauNoRecognizedSecondaryAncestor",
        matchedGV_nDauNoRecognizedSecondaryAncestor,
        "nDauNoRecognizedSecondaryAncestor copied from matched GV");

    table->addColumn<int>("matchedGV_nDauFromB", matchedGV_nDauFromB, "nDauFromB copied from matched GV");
    table->addColumn<int>("matchedGV_nDauFromBC", matchedGV_nDauFromBC, "nDauFromBC copied from matched GV");
    table->addColumn<int>("matchedGV_nDauFromC", matchedGV_nDauFromC, "nDauFromC copied from matched GV");
    table->addColumn<int>("matchedGV_nDauOtherSecondary", matchedGV_nDauOtherSecondary, "nDauOtherSecondary copied from matched GV");
    table->addColumn<int>("matchedGV_nDauOriginUnknown", matchedGV_nDauOriginUnknown, "nDauOriginUnknown copied from matched GV");

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
