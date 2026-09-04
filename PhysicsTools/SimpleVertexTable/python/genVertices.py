import FWCore.ParameterSet.Config as cms

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import (
    selectedHadronsAndPartons,
)
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import (
    ak4JetFlavourInfos,
)


selectedHadronsAndPartonsForGVStudy = (
    selectedHadronsAndPartons.clone(
        particles = cms.InputTag("prunedGenParticles")
    )
)


genJetFlavourInfosForGVStudy = ak4JetFlavourInfos.clone(
    jets = cms.InputTag("slimmedGenJets"),

    bHadrons = cms.InputTag(
        "selectedHadronsAndPartonsForGVStudy",
        "bHadrons",
    ),
    cHadrons = cms.InputTag(
        "selectedHadronsAndPartonsForGVStudy",
        "cHadrons",
    ),
    partons = cms.InputTag(
        "selectedHadronsAndPartonsForGVStudy",
        "physicsPartons",
    ),
    leptons = cms.InputTag(
        "selectedHadronsAndPartonsForGVStudy",
        "leptons",
    ),
)

genCandidateVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(0.),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.5)                                 # dPt/pt between tracks and daughters to be considered matched
)

genCentralVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(1),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(3.0),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.2)                                 # dPt/pt between tracks and daughters to be considered matched
)

genVertexProducer = cms.EDProducer("GenVertexProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    jetFlavourInfos = cms.InputTag("genJetFlavourInfosForGVStudy"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(1),        # number of tracks required to match the genDaughters
    hadPt_min = cms.double(0.0), #Minimum hadron pt for saving gen info
    dlenSigMin = cms.double(0.),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.2),
    doubleMatching = cms.bool(False),
    doubleMatching_nRequiredCommonTracks = cms.int32(3),        # number of tracks required to match the genDaughters
    doubleMatching_maxSignificance = cms.double(999.),
    doubleMatching_dR_max = cms.double(0.05),
    doubleMatching_relPt_max = cms.double(0.4),
    tracks                 = cms.InputTag("unpackedTracksAndVertices"),  # or your favorite general track collection
    trkMaxDeltaR           = cms.double(0.03),
    trkMaxDPtRel           = cms.double(0.2),
    trkCheckCharge         = cms.bool(False),
    trkResolveAmbiguities  = cms.bool(True),
)


svTruthTableProducer = cms.EDProducer(
    "SVTruthTableProducer",
    svTable = cms.InputTag("svTable", "SVTable"),
    svGVMatchTable = cms.InputTag("gvProducer", "SVGVMatchTable")
)


def custom_GV_producer(process, collection="candidate"):
    if collection=="candidate":
        print("Candidate collection is running")
        process.genCandidateVertexProducer = genCandidateVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.genCandidateVertexProducer)
    elif collection=="track":
        print("Track collection is running")
        process.selectedHadronsAndPartonsForGVStudy = (
            selectedHadronsAndPartonsForGVStudy
        )
        process.genJetFlavourInfosForGVStudy = (
            genJetFlavourInfosForGVStudy
        )
        process.gvProducer = genVertexProducer
        process.svTruthTable = svTruthTableProducer
        print(genVertexProducer)
        process.genVertexProducer_sequence = cms.Sequence(
            process.selectedHadronsAndPartonsForGVStudy
            * process.genJetFlavourInfosForGVStudy
            * process.gvProducer
            * process.svTruthTable
        )
    elif collection=="central":
        print("Central collection is running")
        process.gvCentralProducer = genCentralVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.gvCentralProducer)
    return process


#def add_sv_truth_output_commands(process):
#    keep_cmd = "keep *_svTruthTable_SVTruthTable_*"
#    if hasattr(process, "NANOAODEventContent"):
#        process.NANOAODEventContent.outputCommands.append(keep_cmd)
#    if hasattr(process, "NANOEDMAODEventContent"):
#        process.NANOEDMAODEventContent.outputCommands.append(keep_cmd)
#    return process
