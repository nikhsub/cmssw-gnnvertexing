import FWCore.ParameterSet.Config as cms

genCandidateVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(3.0),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.5)                                 # dPt/pt between tracks and daughters to be considered matched
)

genCentralVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(3.0),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.5)                                 # dPt/pt between tracks and daughters to be considered matched
)

genVertexProducer = cms.EDProducer("GenVertexProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(1),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(3.0),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.2),
    doubleMatching = cms.bool(True),
    doubleMatching_maxSignificance = cms.double(999.),
    doubleMatching_dR_max = cms.double(0.01),
    doubleMatching_relPt_max = cms.double(0.1)
)

svTruthTableProducer = cms.EDProducer("SVTruthTableProducer",
    svTable = cms.InputTag("svTable", "SVTable"),
    svTrkTable = cms.InputTag("svTable", "SVtrksTable"),
    gvTable = cms.InputTag("gvProducer", "GVTable"),
    gvDaughtersTable = cms.InputTag("gvProducer", "GVDaughtersTable"),
    nRequiredCommonTracks = cms.int32(1),
    dR_max = cms.double(0.03),
    relPt_max = cms.double(0.2)
)


def custom_GV_producer(process, collection="candidate"):
    if collection=="candidate":
        print("Candidate collection is running")
        process.genCandidateVertexProducer = genCandidateVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.genCandidateVertexProducer)
    elif collection=="track":
        print("Track collection is running")
        process.gvProducer = genVertexProducer
        process.svTruthTable = svTruthTableProducer
        print(genVertexProducer)
        process.genVertexProducer_sequence = cms.Sequence(process.gvProducer*process.svTruthTable)
    elif collection=="central":
        print("Central collection is running")
        process.gvCentralProducer = genCentralVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.gvCentralProducer)
    return process


def add_sv_truth_output_commands(process):
    keep_cmd = "keep *_svTruthTable_SVTruthTable_*"
    if hasattr(process, "NANOAODEventContent"):
        process.NANOAODEventContent.outputCommands.append(keep_cmd)
    if hasattr(process, "NANOEDMAODEventContent"):
        process.NANOEDMAODEventContent.outputCommands.append(keep_cmd)
    return process
