import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar


#packedPFcandidates are filtered (if(c.hasTrackDetails() && c.charge() != 0 && c.numberOfHits()> 0))
unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
    slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    additionalTracks = cms.InputTag("lostTracks"),
    packedCandidates = cms.InputTag("packedPFCandidates")
)

dummyValueMap = cms.EDProducer("DummyTrackValueMap",
    src = cms.InputTag("unpackedTracksAndVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    model_path = cms.FileInPath("PhysicsTools/data/submod_out128_hyper_1802.onnx"),
    #threshold = cms.double(0.)
)


# IVF parameter : https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexFinder_cfi.py#L15C1-L16C49
inclusiveVertexFinder = cms.EDProducer('InclusiveVertexFinder',
svScores = cms.InputTag("dummyValueMap", "SVscore"),
  edgeScores = cms.InputTag("dummyValueMap", "edgeScores"),
  edgeIndices = cms.InputTag("dummyValueMap", "edgeIndices"),
  svScoreThreshold = cms.double(0.1),   # minimal SV score to enter in IVF
  seedScoreThreshold = cms.double(0.8), # minimal SV score to be a seed for IVF
  edgeScoreThreshold = cms.double(0.7),  # minimal edgeScore_ij of track j to be included in seed of track i
  beamSpot = cms.InputTag('offlineBeamSpot'),
  clusterizer = cms.PSet(
   clusterMaxDistance = cms.double(0.1),   #default 0.05
    clusterMaxSignificance = cms.double(9.0), #default 4.5
    clusterMinAngleCosine = cms.double(0.25), #default 0.5
    distanceRatio = cms.double(10.),    # default 20
    maxTimeSignificance = cms.double(3.5), #default 3.5
    seedMax3DIPSignificance = cms.double(9999), #disabled
    seedMax3DIPValue = cms.double(9999), #disabled
    seedMin3DIPSignificance = cms.double(0.6), #defaul 1.2
    seedMin3DIPValue = cms.double(0.0025), # defailt 0.005

  ),
  fitterRatio = cms.double(0.25),
  fitterSigmacut = cms.double(3),
  fitterTini = cms.double(256),
  maxNTracks = cms.uint32(30),
  maximumLongitudinalImpactParameter = cms.double(0.3), # 0.3
  maximumTimeSignificance = cms.double(3), # default 3.0
  minHits = cms.uint32(8), #8
  minPt = cms.double(0.8), #0.8
  primaryVertices = cms.InputTag('unpackedTracksAndVertices'),
  #tracks = cms.InputTag('dummyValueMap', 'selectedTracks'),
  tracks = cms.InputTag('unpackedTracksAndVertices'),
  useDirectVertexFitter = cms.bool(True),
  useVertexReco = cms.bool(True),
  vertexMinAngleCosine = cms.double(0.25), #0.95 default
  vertexMinDLen2DSig = cms.double(0.125), #2.5 default
  vertexMinDLenSig = cms.double(0.125), #0.5 default
  vertexReco = cms.PSet(
    finder = cms.string('avr'),
    primcut = cms.double(1),
    seccut = cms.double(3),
    smoothing = cms.bool(True)
  ),
)

#Vertex Merger step1 https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/RecoVertex/AdaptiveVertexFinder/python/vertexMerger_cfi.py
vertexMerger = cms.EDProducer( "VertexMerger",
    secondaryVertices = cms.InputTag("inclusiveVertexFinder"),  
    maxFraction = cms.double(0.7), 
    minSignificance = cms.double(2.0),
    edgeScores = cms.InputTag("dummyValueMap", "edgeScores"),
    edgeIndices = cms.InputTag("dummyValueMap", "edgeIndices"),
    useEdgeScore = cms.bool(False)

)

#Arbitrator step
# https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/trackVertexArbitrator_cfi.py#L16
trackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("unpackedTracksAndVertices"),
    tracks = cms.InputTag("unpackedTracksAndVertices"),
    secondaryVertices = cms.InputTag("vertexMerger"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    sigCut = cms.double(5),
    fitterSigmacut =  cms.double(3),
    fitterTini = cms.double(256),
    fitterRatio = cms.double(0.25),
    trackMinLayers = cms.int32(4),
    trackMinPt = cms.double(0.4),
    trackMinPixels = cms.int32(1)
    # plus any additional parameters it requires
)

#Vertex Merger step2 https://github.com/cms-sw/cmssw/blob/557f39bce1d5cba35316c2358a89e888901a07e5/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexing_cff.py#L7
myFinalInclusiveSecondaryVertices = vertexMerger.clone(
    secondaryVertices = "trackVertexArbitrator",
    maxFraction = cms.double(1.0), #0.2 default
    minSignificance = cms.double(0.),  #10 default
    edgeScores = cms.InputTag("dummyValueMap", "edgeScores"),
    edgeIndices = cms.InputTag("dummyValueMap", "edgeIndices"),
    useEdgeScore = cms.bool(True)
    ) 



svTable = cms.EDProducer("SVTableProducer", 
                        pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                        src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
                        dlenSigMin = cms.double(0.0))


svGraphGNNInference = cms.EDProducer("SVGraphGNNInferenceProducer",
    src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    trackSrc = cms.InputTag("unpackedTracksAndVertices"),
    model_path = cms.FileInPath("PhysicsTools/data/vertex_gnn_0106.onnx"),
    maxTracks = cms.uint32(16),
    maxEdges = cms.uint32(128),
    dlenSigMin = cms.double(0.0),
    debug = cms.untracked.bool(False)
)

svGraphVertexGNNInference = cms.EDProducer("SVGraphVertexGNNInferenceProducer",
    src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    trackSrc = cms.InputTag("unpackedTracksAndVertices"),
    globalTrackIdxMap = cms.InputTag("dummyValueMap", "globalTrackIdxMap"),
    model_path = cms.FileInPath("PhysicsTools/data/vertex_gnn_newother_1407.onnx"),
    maxTracks = cms.uint32(32),
    maxEdges = cms.uint32(128),
    dlenSigMin = cms.double(0.0),
    includeNearbyTracks = cms.bool(True),
    nearbyTrackDR = cms.double(0.4),
    nearbyTrackPtMin = cms.double(0.8),
    nearbyTrackEtaMax = cms.double(2.5),
    maxExtraTracks = cms.int32(5),
    requireExtraHighPurity = cms.bool(True),
    debug = cms.untracked.bool(False)
)


# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig

trackVertexVars = cms.EDProducer("TrackVertexVars",
    tracks = cms.InputTag("unpackedTracksAndVertices"),
    primaryVertices = cms.InputTag("unpackedTracksAndVertices"),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    #Arguments to know if a track was a seed
    seedMax3DIPSignificance = cms.double(9999),   # disabled
    seedMax3DIPValue        = cms.double(9999),   # disabled
    seedMin3DIPSignificance = cms.double(0.6),    # default 1.2
    seedMin3DIPValue        = cms.double(0.0025), # default 0.005
)

trackGenMatch = cms.EDProducer('TrackGenMatcher',
    tracks       = cms.InputTag("unpackedTracksAndVertices"),
    genParticles = cms.InputTag("mergedGenParticles"),  # or your favorite gen collection
    mcPdgId      = cms.vint32(),                # abs(pdgId) to match against; empty vint32() = any
    mcStatus     = cms.vint32(1),                        # status==1 (stable); empty vint32() = any
    checkCharge  = cms.bool(False),
    maxDeltaR    = cms.double(0.03),
    maxDPtRel    = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),                 # one-to-one matching
)
trackTable = cms.EDProducer(
    "SimpleTrackFlatTableProducer",
    src  = cms.InputTag("unpackedTracksAndVertices"),
    name = cms.string("track"),
    doc  = cms.string("reconstructed tracks (PF candidates + lost tracks, unpacked)"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt     = Var("pt", "float", doc="track transverse momentum"),
        eta    = Var("eta", "float", doc="track pseudorapidity"),
        phi    = Var("phi", "float", doc="track azimuthal angle"),
        charge = Var("charge", "int", doc="track charge"),
        nHits = Var("hitPattern().numberOfValidHits()", "int", doc="number of valid track hits")
    ),
    externalVariables = cms.PSet(
        # ValueMap<int> produced by TrackGenMatcher, keyed to the SAME
        # "unpackedTracksAndVertices" collection used as src above.
        genPartIdx = ExtVar(cms.InputTag("trackGenMatch", "genPartIdx"),"int",doc="index of the matched gen particle in the genMatch source collection, -1 if unmatched"),
        SVscore = ExtVar(cms.InputTag("dummyValueMap", "SVscore"), "float",doc="GNN-based SV score for the track",),
        dz         = ExtVar(cms.InputTag("trackVertexVars", "dz"), "float",doc="track dz w.r.t. the primary vertex"),
        timeSig    = ExtVar(cms.InputTag("trackVertexVars", "timeSig"), "float",doc="track time significance w.r.t. the primary vertex"),
        ip3dValue        = ExtVar(cms.InputTag("trackVertexVars", "ip3dValue"), "float",doc="absolute 3D impact parameter value w.r.t. leading PV"),
        ip3dSignificance = ExtVar(cms.InputTag("trackVertexVars", "ip3dSignificance"), "float",doc="absolute 3D impact parameter significance w.r.t. leading PV"),
        passSeed         = ExtVar(cms.InputTag("trackVertexVars", "passSeed"), "int",doc="1 if track passes TracksClusteringFromDisplacedSeed's seed window cut"),
    ),
)


def custom_sv_tracks(process, threshold_values=(0.0, 0., 0.), f1=0.7, minSig1=2, f2=1., minSig2=0.):
  process.unpackedTracksAndVertices = unpackedTracksAndVertices
  process.inclusiveVertexFinder = inclusiveVertexFinder.clone(
    svScoreThreshold=cms.double(threshold_values[0]),
    seedScoreThreshold=cms.double(threshold_values[1]),
    edgeScoreThreshold=cms.double(threshold_values[2]))
  print(f"Using custom thresholds for SV reconstruction: svScoreThreshold={threshold_values[0]}, seedScoreThreshold={threshold_values[1]}, edgeScoreThreshold={threshold_values[2]}")
  print(f"Using custom values for f1: {f1}, minS1: {minSig1}, f2: {f2}, minS2: {minSig2} ")
  process.vertexMerger = vertexMerger.clone(
    maxFraction = cms.double(f1),
    minSignificance = cms.double(minSig1)
  )
  process.trackVertexArbitrator = trackVertexArbitrator
  process.myFinalInclusiveSecondaryVertices = myFinalInclusiveSecondaryVertices.clone(
    maxFraction = cms.double(f2),
    minSignificance = cms.double(minSig2)
  )
  process.svTable = svTable
  process.svGNN = svGraphGNNInference.clone()
  process.svVertexGNN = svGraphVertexGNNInference.clone()
  process.dummyValueMap = dummyValueMap
  process.trackGenMatch = trackGenMatch
  process.trackTable = trackTable
  process.trackVertexVars = trackVertexVars
  process.sv_track = cms.Sequence(    process.unpackedTracksAndVertices*
                                      process.dummyValueMap*
                                      process.inclusiveVertexFinder*
                                      process.vertexMerger*
                                      process.trackVertexArbitrator*
                                      process.myFinalInclusiveSecondaryVertices*
                                      process.svTable*
                                      process.trackVertexVars*
                                      process.trackGenMatch*
                                      process.trackTable
                                      )
  process.sv_track += process.svGNN
  process.sv_track += process.svVertexGNN
  return process
