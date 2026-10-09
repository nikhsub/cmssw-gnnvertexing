import FWCore.ParameterSet.Config as cms
# For candidate collection only
genCandidateVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(0.),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.5)                                 # dPt/pt between tracks and daughters to be considered matched
)
# For central collection only
genCentralVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(3.0),
    dR_max = cms.double(0.03),                                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.2)                                 # dPt/pt between tracks and daughters to be considered matched
)

# For track collection only
genVertexTable = cms.EDProducer("GVProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    minHadronPt = cms.double(1.0),  
    minDaughterPt = cms.double(0.4),  
    maxHadronEta = cms.double(2.5),   
    maxDaughterEta = cms.double(2.5)   
)

gvMatchTable = cms.EDProducer("GVMatchProducer",
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    tracks = cms.InputTag("unpackedTracksAndVertices"),  # or your favorite general track collection

    # outputs from genVertexTable that GVMatchProducer consumes
    hadronGVx      = cms.InputTag("genVertexTable", "hadronGVx"),
    hadronGVy      = cms.InputTag("genVertexTable", "hadronGVy"),
    hadronGVz      = cms.InputTag("genVertexTable", "hadronGVz"),
    hadronEta      = cms.InputTag("genVertexTable", "hadronEta"),
    hadronPhi      = cms.InputTag("genVertexTable", "hadronPhi"),
    daughterPt     = cms.InputTag("genVertexTable", "daughterPt"),
    daughterEta    = cms.InputTag("genVertexTable", "daughterEta"),
    daughterPhi    = cms.InputTag("genVertexTable", "daughterPhi"),
    daughterCharge = cms.InputTag("genVertexTable", "daughterCharge"),
    daughterGVidx  = cms.InputTag("genVertexTable", "daughterGVidx"),
    nGV            = cms.InputTag("genVertexTable", "nGV"),

    # matching parameters, unchanged from the original module
    nRequiredCommonTracks = cms.int32(2),
    dlenSigMin = cms.double(0.3),
    dR_max = cms.double(0.03),
    relPt_max = cms.double(0.2),
    doubleMatching = cms.bool(False),
    doubleMatching_nRequiredCommonTracks = cms.int32(3),
    doubleMatching_maxSignificance = cms.double(999.),
    doubleMatching_dR_max = cms.double(0.05),
    doubleMatching_relPt_max = cms.double(0.4),
    trkMaxDeltaR = cms.double(0.03),
    trkMaxDPtRel = cms.double(0.2),
    trkCheckCharge = cms.bool(False),
    trkResolveAmbiguities = cms.bool(True),
)



def custom_GV_producer(process, collection="candidate"):
    if collection=="candidate":
        print("Candidate collection is running")
        process.genCandidateVertexProducer = genCandidateVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.genCandidateVertexProducer)
    elif collection=="track":
        print("Track collection is running")
        process.gvProducer = genVertexTable
        #process.gvMatchProducer = gvMatchTable
        process.genVertexProducer_sequence = cms.Sequence(process.gvProducer)
        #process.genVertexProducer_sequence = cms.Sequence(
        #    process.gvProducer *
        #    process.gvMatchProducer
        #)
    elif collection=="central":
        print("Central collection is running")
        process.gvCentralProducer = genCentralVertexProducer
        process.genVertexProducer_sequence = cms.Sequence(process.gvCentralProducer)
    return process
