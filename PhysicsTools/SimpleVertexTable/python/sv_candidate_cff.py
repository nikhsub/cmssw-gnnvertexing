import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import P4Vars, Var
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import inclusiveVertexFinder

# Candidate-based modules (without process)
inclusiveCandidateVertexFinder = cms.EDProducer('InclusiveCandidateVertexFinder',
  beamSpot = cms.InputTag('offlineBeamSpot'),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20),
        maxTimeSignificance = cms.double(3.5),
        seedMax3DIPSignificance = cms.double(9999),
        seedMax3DIPValue = cms.double(9999),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3),
    fitterTini = cms.double(256),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3),
    minHits = cms.uint32(0),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag("packedPFCandidates"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1),
        seccut = cms.double(3),
        smoothing = cms.bool(True)
    )
  #mightGet = cms.optional.untracked.vstring
)

candidateVertexMerger = cms.EDProducer( "CandidateVertexMerger",
    secondaryVertices = cms.InputTag("inclusiveCandidateVertexFinder"),  
    maxFraction = cms.double(0.7), 
    minSignificance = cms.double(2.0)
)

CandidateVertexArbitrator = cms.EDProducer("CandidateVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag("packedPFCandidates"),
    secondaryVertices = cms.InputTag("candidateVertexMerger"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    sigCut = cms.double(5),
    fitterSigmacut =  cms.double(3),
    fitterTini = cms.double(256),
    fitterRatio = cms.double(0.25),
    trackMinLayers = cms.int32(4),
    trackMinPt = cms.double(0.4),
    trackMinPixels = cms.int32(1),
    maxTimeSignificance = cms.double(3.5)
)

myFinalInclusiveSecondaryVertices = candidateVertexMerger.clone(
    secondaryVertices = "CandidateVertexArbitrator",
    maxFraction = cms.double(0.2), 
    minSignificance = cms.double(10) )
    
slimmedSecondaryVertices = cms.EDProducer("PATSecondaryVertexSlimmer",
    packedPFCandidates = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag("myFinalInclusiveSecondaryVertices")
)
#mySlimmedSecondaryVertices = cms.EDProducer("PATSecondaryVertexSlimmer",
#    src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
#    packedPFCandidates = cms.InputTag("packedPFCandidates")
#)
myvertexTable = cms.EDProducer("VertexTableProducer",
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), 
    svSrc = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    svCut = cms.string(""),
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(0),
    #storeCharge=cms.bool(False), not present in cmssw15
    pvName = cms.string("myPV"),
    svName = cms.string("mySV"),
    pfcSrc = cms.InputTag("packedPFCandidates"),
    svDoc  = cms.string("secondary vertices from custom IVF algorithm"),
)

svCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("myvertexTable"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("mySV"),
    #svDoc = cms.string("candidate SVs"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        ndof    = Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
        chi2    = Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
        ntracks = Var("numberOfDaughters()", "uint8", doc = "number of tracks"),
    ),
)



def custom_sv_candidate(process):
    process.inclusiveCandidateVertexFinder = inclusiveCandidateVertexFinder
    process.candidateVertexMerger = candidateVertexMerger
    process.CandidateVertexArbitrator = CandidateVertexArbitrator
    process.myFinalInclusiveSecondaryVertices = myFinalInclusiveSecondaryVertices
    process.myvertexTable = myvertexTable
    process.svCandidateTable = svCandidateTable
    process.sv_candidate = cms.Sequence(        process.inclusiveCandidateVertexFinder *
                                            process.candidateVertexMerger *
                                            process.CandidateVertexArbitrator *
                                            process.myFinalInclusiveSecondaryVertices *
                                            process.myvertexTable *
                                            process.svCandidateTable)
    return process
