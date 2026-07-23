#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <memory>
#include <string>
#include <vector>

class SVTruthTableProducer : public edm::global::EDProducer<> {
public:
  explicit SVTruthTableProducer(const edm::ParameterSet& iConfig)
      : svTableToken_(consumes<nanoaod::FlatTable>(
            iConfig.getParameter<edm::InputTag>("svTable"))),
        svGVMatchTableToken_(consumes<nanoaod::FlatTable>(
            iConfig.getParameter<edm::InputTag>("svGVMatchTable"))) {
    produces<nanoaod::FlatTable>("SVTruthTable");
  }

  void produce(edm::StreamID,
               edm::Event& iEvent,
               const edm::EventSetup&) const override {
    edm::Handle<nanoaod::FlatTable> svTable;
    edm::Handle<nanoaod::FlatTable> svGVMatchTable;

    iEvent.getByToken(svTableToken_, svTable);
    iEvent.getByToken(svGVMatchTableToken_, svGVMatchTable);

    if (!svTable.isValid() || !svGVMatchTable.isValid()) {
      throw cms::Exception("MissingInput")
          << "SVTruthTableProducer requires valid svTable and "
          << "svGVMatchTable inputs.";
    }

    const unsigned int nSV = svTable->size();
    if (svGVMatchTable->size() != nSV) {
      throw cms::Exception("TableSizeMismatch")
          << "SVTruthTableProducer received " << nSV
          << " rows in svTable but " << svGVMatchTable->size()
          << " rows in svGVMatchTable. The GenVertexProducer SV selection "
          << "must use the same ordering and dlenSig selection as the SV table.";
    }

    const std::vector<std::string> requiredIntColumns = {
        "GVIdx",
        "isMatched",
        "nMatchedGV",
        "nCommonTracks",
        "nMatchedDaughters",
        "truth_pdgClass",
        "truth_pdgId",
        "truth_isB",
        "truth_isD",
        "truth_isBtoD",
        "nOriginNoRecognizedSecondaryAncestor",
        "nOriginFromB",
        "nOriginFromBC",
        "nOriginFromC",
        "nOriginOtherSecondary",
        "nOriginUnknown",
        "dominantOriginLabel",
        "matchedGV_nDauNoRecognizedSecondaryAncestor",
        "matchedGV_nDauFromB",
        "matchedGV_nDauFromBC",
        "matchedGV_nDauFromC",
        "matchedGV_nDauOtherSecondary",
        "matchedGV_nDauOriginUnknown"};

    const std::vector<std::string> requiredFloatColumns = {
        "bestMatchScore",
        "bestMatchDistanceSig",
        "fracOriginNoRecognizedSecondaryAncestor",
        "fracOriginHF"};

    for (const auto& name : requiredIntColumns) {
      requireColumn(*svGVMatchTable, name);
    }
    for (const auto& name : requiredFloatColumns) {
      requireColumn(*svGVMatchTable, name);
    }

    std::vector<int> svIdx(nSV, -1);
    for (unsigned int i = 0; i < nSV; ++i) {
      svIdx[i] = static_cast<int>(i);
    }

    const auto matchedGVIdx = readIntColumn(*svGVMatchTable, "GVIdx");
    const auto nMatchedGV = readIntColumn(*svGVMatchTable, "nMatchedGV");
    const auto nCommonTracks = readIntColumn(*svGVMatchTable, "nCommonTracks");
    const auto nMatchedDaughters = readIntColumn(*svGVMatchTable, "nMatchedDaughters");

    const auto truthPdgClass = readIntColumn(*svGVMatchTable, "truth_pdgClass");
    const auto truthPdgId = readIntColumn(*svGVMatchTable, "truth_pdgId");
    const auto truthIsB = readIntColumn(*svGVMatchTable, "truth_isB");
    const auto truthIsD = readIntColumn(*svGVMatchTable, "truth_isD");
    const auto truthIsBtoD = readIntColumn(*svGVMatchTable, "truth_isBtoD");

    const auto bestMatchScore = readFloatColumn(*svGVMatchTable, "bestMatchScore");
    const auto bestMatchDistanceSig =
        readFloatColumn(*svGVMatchTable, "bestMatchDistanceSig");

    const auto nOriginNoRecognizedSecondaryAncestor =
        readIntColumn(*svGVMatchTable, "nOriginNoRecognizedSecondaryAncestor");
    const auto nOriginFromB = readIntColumn(*svGVMatchTable, "nOriginFromB");
    const auto nOriginFromBC = readIntColumn(*svGVMatchTable, "nOriginFromBC");
    const auto nOriginFromC = readIntColumn(*svGVMatchTable, "nOriginFromC");
    const auto nOriginOtherSecondary =
        readIntColumn(*svGVMatchTable, "nOriginOtherSecondary");
    const auto nOriginUnknown = readIntColumn(*svGVMatchTable, "nOriginUnknown");

    const auto fracOriginNoRecognizedSecondaryAncestor =
        readFloatColumn(*svGVMatchTable,
                        "fracOriginNoRecognizedSecondaryAncestor");
    const auto fracOriginHF = readFloatColumn(*svGVMatchTable, "fracOriginHF");
    const auto dominantOriginLabel =
        readIntColumn(*svGVMatchTable, "dominantOriginLabel");

    const auto matchedGVNDauNoRecognizedSecondaryAncestor =
        readIntColumn(*svGVMatchTable,
                      "matchedGV_nDauNoRecognizedSecondaryAncestor");
    const auto matchedGVNDauFromB =
        readIntColumn(*svGVMatchTable, "matchedGV_nDauFromB");
    const auto matchedGVNDauFromBC =
        readIntColumn(*svGVMatchTable, "matchedGV_nDauFromBC");
    const auto matchedGVNDauFromC =
        readIntColumn(*svGVMatchTable, "matchedGV_nDauFromC");
    const auto matchedGVNDauOtherSecondary =
        readIntColumn(*svGVMatchTable, "matchedGV_nDauOtherSecondary");
    const auto matchedGVNDauOriginUnknown =
        readIntColumn(*svGVMatchTable, "matchedGV_nDauOriginUnknown");

    auto table = std::make_unique<nanoaod::FlatTable>(nSV, "SVTruth", false, false);

    table->addColumn<int>("svIdx", svIdx, "SV row index");
    table->addColumn<int>("truth_pdgClass", truthPdgClass,
                          "Primary truth class label from matched GV");
    table->addColumn<int>("truth_pdgId", truthPdgId,
                          "PDG id of matched truth hadron");
    table->addColumn<int>("truth_isB", truthIsB,
                          "Matched GV is a B hadron");
    table->addColumn<int>("truth_isD", truthIsD,
                          "Matched GV is a D hadron");
    table->addColumn<int>("truth_isBtoD", truthIsBtoD,
                          "Matched GV is a D hadron with B ancestry");

    table->addColumn<int>("matchedGVIdx", matchedGVIdx,
                          "GV row index matched to this SV");
    table->addColumn<int>("nMatchedGV", nMatchedGV,
                          "One if this SV was matched to a GV, else zero");
    table->addColumn<int>("nCommonTracks", nCommonTracks,
                          "Number of matched SV tracks from GenVertexProducer");
    table->addColumn<int>("nMatchedDaughters", nMatchedDaughters,
                          "Number of matched GV daughters from GenVertexProducer");

    table->addColumn<float>("bestMatchScore", bestMatchScore,
                            "Track-sharing score from GenVertexProducer");
    table->addColumn<float>("bestMatchDistanceSig", bestMatchDistanceSig,
                            "Matched SV-GV covariance-weighted distance significance");

    table->addColumn<int>("nOriginNoRecognizedSecondaryAncestor",
                          nOriginNoRecognizedSecondaryAncestor,
                          "Matched daughters with originLabel 0");
    table->addColumn<int>("nOriginFromB", nOriginFromB,
                          "Matched daughters with B ancestry");
    table->addColumn<int>("nOriginFromBC", nOriginFromBC,
                          "Matched daughters with B and C ancestry");
    table->addColumn<int>("nOriginFromC", nOriginFromC,
                          "Matched daughters with C ancestry");
    table->addColumn<int>("nOriginOtherSecondary", nOriginOtherSecondary,
                          "Matched daughters from other secondary sources");
    table->addColumn<int>("nOriginUnknown", nOriginUnknown,
                          "Matched daughters with unknown origin");

    table->addColumn<float>("fracOriginNoRecognizedSecondaryAncestor",
                            fracOriginNoRecognizedSecondaryAncestor,
                            "Fraction of matched daughters with originLabel 0");
    table->addColumn<float>("fracOriginHF", fracOriginHF,
                            "Fraction of matched daughters with B/BC/C origin labels");
    table->addColumn<int>("dominantOriginLabel", dominantOriginLabel,
                          "Most frequent matched-daughter origin label");

    table->addColumn<int>("matchedGV_nDauNoRecognizedSecondaryAncestor",
                          matchedGVNDauNoRecognizedSecondaryAncestor,
                          "Matched GV daughter count with originLabel 0");
    table->addColumn<int>("matchedGV_nDauFromB", matchedGVNDauFromB,
                          "Matched GV daughter count from B");
    table->addColumn<int>("matchedGV_nDauFromBC", matchedGVNDauFromBC,
                          "Matched GV daughter count from BC");
    table->addColumn<int>("matchedGV_nDauFromC", matchedGVNDauFromC,
                          "Matched GV daughter count from C");
    table->addColumn<int>("matchedGV_nDauOtherSecondary",
                          matchedGVNDauOtherSecondary,
                          "Matched GV other-secondary daughter count");
    table->addColumn<int>("matchedGV_nDauOriginUnknown",
                          matchedGVNDauOriginUnknown,
                          "Matched GV unknown-origin daughter count");

    iEvent.put(std::move(table), "SVTruthTable");
  }

private:
  static void requireColumn(const nanoaod::FlatTable& table,
                            const std::string& name) {
    if (table.columnIndex(name) < 0) {
      throw cms::Exception("MissingColumn")
          << "SVTruthTableProducer requires column '" << name
          << "' in the GenVertexProducer SVGVMatchTable.";
    }
  }

  static std::vector<int> readIntColumn(const nanoaod::FlatTable& table,
                                        const std::string& name) {
    const int column = table.columnIndex(name);
    std::vector<int> values(table.size(), 0);
    for (unsigned int row = 0; row < table.size(); ++row) {
      values[row] = static_cast<int>(table.getAnyValue(row, column));
    }
    return values;
  }

  static std::vector<float> readFloatColumn(const nanoaod::FlatTable& table,
                                            const std::string& name) {
    const int column = table.columnIndex(name);
    std::vector<float> values(table.size(), 0.f);
    for (unsigned int row = 0; row < table.size(); ++row) {
      values[row] = static_cast<float>(table.getAnyValue(row, column));
    }
    return values;
  }

  edm::EDGetTokenT<nanoaod::FlatTable> svTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> svGVMatchTableToken_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVTruthTableProducer);
