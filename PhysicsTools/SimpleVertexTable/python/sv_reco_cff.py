import FWCore.ParameterSet.Config as cms

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
    threshold = cms.double(0.)
)


# IVF parameter : https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexFinder_cfi.py#L15C1-L16C49
inclusiveVertexFinder = cms.EDProducer('InclusiveVertexFinder',
  beamSpot = cms.InputTag('offlineBeamSpot'),
  clusterizer = cms.PSet(
    clusterMaxDistance = cms.double(0.05), # requirement on adding a track to the cluster
    clusterMaxSignificance = cms.double(4.5), # requirement on adding a track to the cluster
    clusterMinAngleCosine = cms.double(0.5), # requirement on adding a track to the cluster
    distanceRatio = cms.double(20),
    maxTimeSignificance = cms.double(3.5),
    seedMax3DIPSignificance = cms.double(9999), #disabled
    seedMax3DIPValue = cms.double(9999), #disabled
    seedMin3DIPSignificance = cms.double(1.2), # which tracks can start a cluster
    seedMin3DIPValue = cms.double(0.005),  # which tracks can start a cluster
  ),
  fitterRatio = cms.double(0.25),
  fitterSigmacut = cms.double(3),
  fitterTini = cms.double(256),
  maxNTracks = cms.uint32(30),
  maximumLongitudinalImpactParameter = cms.double(0.3),
  maximumTimeSignificance = cms.double(3), # new?
  minHits = cms.uint32(4), #8
  minPt = cms.double(0.8),
  primaryVertices = cms.InputTag('unpackedTracksAndVertices'),
  tracks = cms.InputTag('dummyValueMap', 'selectedTracks'),
  #tracks = cms.InputTag('unpackedTracksAndVertices'),
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
  ),
)

#Vertex Merger step1 https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/RecoVertex/AdaptiveVertexFinder/python/vertexMerger_cfi.py
vertexMerger = cms.EDProducer( "VertexMerger",
    secondaryVertices = cms.InputTag("inclusiveVertexFinder"),  
    maxFraction = cms.double(0.7), 
    minSignificance = cms.double(2.0)
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
    maxFraction = cms.double(0.2), 
    minSignificance = cms.double(10) )


svTable = cms.EDProducer("SVTableProducer", 
                        pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                        src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
                        dlenSigMin = cms.double(3.0))


svGraphGNNInference = cms.EDProducer("SVGraphGNNInferenceProducer",
    src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    trackSrc = cms.InputTag("unpackedTracksAndVertices"),
    model_path = cms.FileInPath("PhysicsTools/data/vertex_gnn_0106.onnx"),
    maxTracks = cms.uint32(16),
    maxEdges = cms.uint32(128),
    dlenSigMin = cms.double(3.0),
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
    dlenSigMin = cms.double(3.0),
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



def custom_sv_tracks(process, threshold_value=0.0):
    process.unpackedTracksAndVertices = unpackedTracksAndVertices
    process.inclusiveVertexFinder = inclusiveVertexFinder
    process.vertexMerger = vertexMerger
    process.trackVertexArbitrator = trackVertexArbitrator
    process.myFinalInclusiveSecondaryVertices = myFinalInclusiveSecondaryVertices
    process.svTable = svTable
    process.svGNN = svGraphGNNInference.clone()
    process.svVertexGNN = svGraphVertexGNNInference.clone()
    process.dummyValueMap = dummyValueMap.clone(
          threshold = cms.double(threshold_value)
      )
    process.sv_track = cms.Sequence(   process.unpackedTracksAndVertices*
                                        process.dummyValueMap*
                                        process.inclusiveVertexFinder*
                                        process.vertexMerger*
                                        process.trackVertexArbitrator*
                                        process.myFinalInclusiveSecondaryVertices*
                                        process.svTable)

    process.sv_track += process.svGNN
    process.sv_track += process.svVertexGNN
    return process
