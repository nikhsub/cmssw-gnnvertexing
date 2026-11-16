import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar

minTrackPt = 0.4


#packedPFcandidates are filtered (if(c.hasTrackDetails() && c.charge() != 0 && c.numberOfHits()> 0))
unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
    slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    additionalTracks = cms.InputTag("lostTracks"),
    #additionalTracks = cms.InputTag(""),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    minTrackPt = cms.double(minTrackPt),
    maxTrackEta = cms.double(2.5),
)

filteredTracks = cms.EDProducer(
    "PATTrackFromVertexFilter",
    tracks          = cms.InputTag("unpackedTracksAndVertices"),   # your unpacker output
    primaryVertices = cms.InputTag("unpackedTracksAndVertices"),   # same unpacker, default label
    dzCut           = cms.double(1.5),   # cm #1.5
)

modelInference = cms.EDProducer("modelInference",
    src = cms.InputTag("filteredTracks"),
    pvSrc = cms.InputTag("filteredTracks"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    dRMatchMax = cms.double(0.001),
    relPtMatchMax = cms.double(0.2),
    trackPairSrc = cms.InputTag("trackPairFeatures", "trackPair"),
    model_path = cms.FileInPath("PhysicsTools/data/best_20261001_142559.onnx"),
)
trackPairFeatures = cms.EDProducer("TrackPairFeaturesProducer",
    tracks = cms.InputTag("filteredTracks"),
    primaryVertices = cms.InputTag("filteredTracks"),
    maxDR = cms.double(0.8),
    maxDEta = cms.double(0.4),
    maxDPhi = cms.double(0.5),
    maxPairEta = cms.double(2.5),
    maxInvMass = cms.double(5.0),
    maxDca = cms.double(0.4),
    # Optional two-track Kalman fit for each selected trackPair row.
    # Off by default so existing reconstruction and label production are unchanged.
    enableTwoTrackVertexFit = cms.bool(False),
    enableVertexFitTiming = cms.bool(False),
)
# Truth-dependent labeling (per-track hadidx/label, per-pair edgeLabel). Needs
# gvTable/gvDaughtersTable/genPartIdx (generator truth) -- MC-only, training/ntuple
# production sequence only, never in the general reconstruction sequence below.
trackPairLabels = cms.EDProducer("TrackPairLabelProducer",
    tracks = cms.InputTag("filteredTracks"),
    gvTable = cms.InputTag("genVertexTable", "GVTable"),
    gvDaughtersTable = cms.InputTag("genVertexTable", "GVDaughtersTable"),
    genPartIdx = cms.InputTag("trackGenMatch", "genPartIdx"),
    trackPairSrc = cms.InputTag("trackPairFeatures", "trackPair"),
    relPtMatchMax = cms.double(0.2),
    dRMatchMax = cms.double(0.02),
    trkPtCut = cms.double(minTrackPt),
)
 

# IVF parameter : https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexFinder_cfi.py#L15C1-L16C49
inclusiveVertexFinder = cms.EDProducer('InclusiveVertexFinder',
  svScores = cms.InputTag("modelInference", "nodeScoreSVCombined"),
  edgeScores = cms.InputTag("modelInference", "edgeScores"),
  edgeIndices = cms.InputTag("modelInference", "edgeIndices"),
  svScoreThreshold = cms.double(0.),   # minimal SV score to enter in IVF
  seedScoreThreshold = cms.double(0.), # minimal SV score to be a seed for IVF
  edgeScoreThreshold = cms.double(0.),  # minimal edgeScore_ij of track j to be included in seed of track i
  beamSpot = cms.InputTag('offlineBeamSpot'),
  clusterizer = cms.PSet(
    #clusterMaxDistance = cms.double(0.1), # 0.05 default | requirement on adding a track to the cluster
    #clusterMaxSignificance = cms.double(9.), # 4.5 default |requirement on adding a track to the cluster
    #clusterMinAngleCosine = cms.double(0.25), # -2 to disable. | 0.5 default requirement on adding a track to the cluster
    #distanceRatio = cms.double(10.),
    clusterMaxDistance = cms.double(0.1),   #default 0.05
    clusterMaxSignificance = cms.double(4.5), #default 4.5 / 9
    clusterMinAngleCosine = cms.double(0.5), #default 0.5 / 0.25
    distanceRatio = cms.double(1.25),    # default 20 / 10
    maxTimeSignificance = cms.double(3.5), #default 3.5
    seedMax3DIPSignificance = cms.double(9999), #disabled
    seedMax3DIPValue = cms.double(9999), #disabled
    #seedMin3DIPSignificance = cms.double(0.6), # which tracks can start a cluster
    #seedMin3DIPValue = cms.double(0.0025),  # which tracks can start a cluster
    seedMin3DIPSignificance = cms.double(1.2), #defaul 1.2 / 0.6
    seedMin3DIPValue = cms.double(0.00125), # defailt 0.005 / 0.0025
  ),
  fitterRatio = cms.double(0.25),
  fitterSigmacut = cms.double(3),
  fitterTini = cms.double(256),
  maxNTracks = cms.uint32(30),
  maximumLongitudinalImpactParameter = cms.double(0.3), # 0.3
  maximumTimeSignificance = cms.double(3), # default 3.0
  minHits = cms.uint32(4), #8 / 4
  minPt = cms.double(minTrackPt), #0.8
  primaryVertices = cms.InputTag('filteredTracks'),
  tracks = cms.InputTag('filteredTracks'),
  useDirectVertexFitter = cms.bool(True),
  useVertexReco = cms.bool(True),
  vertexMinAngleCosine = cms.double(0.5), #0.95 default  / 0.25
  vertexMinDLen2DSig = cms.double(1.25), #2.5 default //2 / 0.125
  vertexMinDLenSig = cms.double(1.5), #0.5 default / 0.125
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
    edgeScores = cms.InputTag("modelInference", "edgeScores"),
    edgeIndices = cms.InputTag("modelInference", "edgeIndices"),
    useEdgeScore = cms.bool(False)
)

#Arbitrator step
# https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/trackVertexArbitrator_cfi.py#L16
trackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("filteredTracks"),
    tracks = cms.InputTag("filteredTracks"),
    secondaryVertices = cms.InputTag("vertexMerger"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    maxTimeSignificance = cms.double(3.5),
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
    maxFraction = cms.double(0.2), #0.2 default/1
    minSignificance = cms.double(10.),  #10 default/0
    edgeScores = cms.InputTag("modelInference", "edgeScores"),
    edgeIndices = cms.InputTag("modelInference", "edgeIndices"),
    useEdgeScore = cms.bool(False)
    ) 


svTable = cms.EDProducer("SVTableProducer", 
                        pvSrc=cms.InputTag("filteredTracks"),
                        src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
                        dlenSigMin = cms.double(2)) #3/0


# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig

trackVertexVars = cms.EDProducer("TrackVertexVars",
    tracks = cms.InputTag("filteredTracks"),
    primaryVertices = cms.InputTag("filteredTracks"),
    muons = cms.InputTag("slimmedMuons"),          # NEW
    electrons = cms.InputTag("slimmedElectrons"),  # NEW
    dRMatchMax = cms.double(0.001),      # NEW — tight by default, tracks are matched 1:1 to a lepton, not a jet-level cone
    relPtMatchMax = cms.double(0.2),    # NEW
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    #Arguments to know if a track was a seed
    seedMax3DIPSignificance = cms.double(9999),   # disabled
    seedMax3DIPValue        = cms.double(9999),   # disabled
    seedMin3DIPSignificance = cms.double(1.2),    # default 1.2
    seedMin3DIPValue        = cms.double(0.005), # default 0.005
)

trackGenMatch = cms.EDProducer('TrackGenMatcher',
    tracks       = cms.InputTag("filteredTracks"),
    genParticles = cms.InputTag("customMergedGenParticles"),  # or your favorite gen collection
    mcPdgId      = cms.vint32(),                # abs(pdgId) to match against; empty vint32() = any
    mcStatus     = cms.vint32(1),                        # status==1 (stable); empty vint32() = any
    checkCharge  = cms.bool(False),
    maxDeltaR    = cms.double(0.03),
    maxDPtRel    = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),                 # one-to-one matching
)
trackTable = cms.EDProducer(
    "SimpleTrackFlatTableProducer",
    src  = cms.InputTag("filteredTracks"),
    name = cms.string("track"),
    doc  = cms.string("reconstructed tracks (PF candidates + lost tracks, unpacked)"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt     = Var("pt", "float", doc="track transverse momentum"),
        px     = Var("px", "float", doc="track x momentum"),
        py     = Var("py", "float", doc="track y momentum"),
        pz     = Var("pz", "float", doc="track z momentum"),
        energy = Var("sqrt(p()*p()+0.13957*0.13957)", "float", doc="track energy (assuming pion mass hypothesis)"),
        qoverp = Var("qoverp", "float", doc="track charge/momentum"),
        qdotp = Var("charge*p()", "float", doc="charge times momentum"),
        chi2   = Var("chi2", "float", doc="track chi2"),
        isHighPurity = Var("quality('highPurity')", "bool", doc="track is high purity"),
        numberOfPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int", doc="number of valid pixel hits"),
        numberOfStripHits = Var("hitPattern().numberOfValidStripHits()", "int", doc="number of valid strip hits"),
        #lostInnerHits 
        eta    = Var("eta", "float", doc="track pseudorapidity"),
        phi    = Var("phi", "float", doc="track azimuthal angle"),
        charge = Var("charge", "int", doc="track charge"),
        nHits = Var("hitPattern().numberOfValidHits()", "int", doc="number of valid track hits")
        
    ),
    externalVariables = cms.PSet(
        # ValueMap<int> produced by TrackGenMatcher, keyed to the SAME
        # "filteredTracks" collection used as src above.
        # MC-only: removed for data in trackTableFor(isMC=False)
        genPartIdx = ExtVar(cms.InputTag("trackGenMatch", "genPartIdx"),"int",doc="index of the matched gen particle in the genMatch source collection, -1 if unmatched"),
        dz         = ExtVar(cms.InputTag("trackVertexVars", "dz"), "float",doc="track dz w.r.t. the primary vertex"),
        isMuon = ExtVar(cms.InputTag("trackVertexVars", "isMuon"), "int",doc="track is a muon"),
        isElectron = ExtVar(cms.InputTag("trackVertexVars", "isElectron"), "int",doc="track is an electron"),
        dzSignificance = ExtVar(cms.InputTag("trackVertexVars", "dzSignificance"), "float",doc="track dz significance w.r.t. the primary vertex"),
        timeSig    = ExtVar(cms.InputTag("trackVertexVars", "timeSig"), "float",doc="track time significance w.r.t. the primary vertex"),
        ip3dValue        = ExtVar(cms.InputTag("trackVertexVars", "ip3dValue"), "float",doc="absolute 3D impact parameter value w.r.t. leading PV"),
        ip3dSignificance = ExtVar(cms.InputTag("trackVertexVars", "ip3dSignificance"), "float",doc="absolute 3D impact parameter significance w.r.t. leading PV"),
        ip2dValue        = ExtVar(cms.InputTag("trackVertexVars", "ip2dValue"), "float",doc="absolute 2D (transverse) impact parameter value w.r.t. leading PV"),
        ip2dSignificance = ExtVar(cms.InputTag("trackVertexVars", "ip2dSignificance"), "float",doc="absolute 2D (transverse) impact parameter significance w.r.t. leading PV"),
        passSeed         = ExtVar(cms.InputTag("trackVertexVars", "passSeed"), "int",doc="1 if track passes TracksClusteringFromDisplacedSeed's seed window cut"),

        #Uncomment for inference
        nodeScore_0 = ExtVar(cms.InputTag("modelInference", "nodeScore0"), "float", doc="GNN node class 0 score"),
        nodeScore_1 = ExtVar(cms.InputTag("modelInference", "nodeScore1"), "float", doc="GNN node class 1 score"),
        nodeScore_2 = ExtVar(cms.InputTag("modelInference", "nodeScore2"), "float", doc="GNN node class 2 score"),
        nodeScore_3 = ExtVar(cms.InputTag("modelInference", "nodeScore3"), "float", doc="GNN node class 3 score"),
        nodeScore_4 = ExtVar(cms.InputTag("modelInference", "nodeScore4"), "float", doc="GNN node class 4 score"),
        nodeScore_5 = ExtVar(cms.InputTag("modelInference", "nodeScore5"), "float", doc="GNN node class 5 score"),
        svScores = ExtVar(cms.InputTag("modelInference", "nodeScoreSVCombined"), "float",doc="GNN combined SV score: (score2+score3+score4+score5)/sum(score0..5)"),

    ),
)


# ---------------------------------------------------------------------------
# MC/data helpers
# ---------------------------------------------------------------------------

def trackTableFor(isMC, withInference=True):
  """Return a clone of trackTable suitable for MC or data.
  Data: drop the gen-matching column (trackGenMatch does not exist).
  withInference=False: drop the GNN score columns (training ntuples)."""
  t = trackTable.clone()
  if not isMC:
    delattr(t.externalVariables, "genPartIdx")
  if not withInference:
    for name in ["nodeScore_%d" % i for i in range(6)] + ["svScores"]:
      delattr(t.externalVariables, name)
  return t


def _chain(modules):
  """Build a cms.Sequence from a list of modules, preserving order."""
  seq = modules[0]
  for m in modules[1:]:
    seq = seq * m
  return cms.Sequence(seq)


def apply_overrides(process, overrides):
  # overrides: {parameter name: value}; see the options registered in main_cfg.py
  for name, val in overrides.items():
    if name in ('clusterMaxDistance', 'clusterMaxSignificance', 'clusterMinAngleCosine', 'distanceRatio',
                'seedMin3DIPSignificance', 'seedMin3DIPValue'):
      setattr(process.inclusiveVertexFinder.clusterizer, name, cms.double(val))
    elif name in ('vertexMinAngleCosine', 'vertexMinDLen2DSig', 'vertexMinDLenSig'):
      setattr(process.inclusiveVertexFinder, name, cms.double(val))
    elif name == 'dlenSigMin': # svTable and the GV table must use the same cut, otherwise the mySV tables mismatch
      process.svTable.dlenSigMin = cms.double(val)
      if hasattr(process, "gvMatchTable"):   # MC only
        process.gvMatchTable.dlenSigMin = cms.double(val)
    else:
      raise ValueError("Unknown override: "+name)
  print("sv_reco overrides:", overrides)


def custom_sv_tracks(process, threshold_values=(0.0, 0., 0.), f1=0.7, minSig1=2, f2=1., minSig2=0., overrides=None, isMC=True,
                     enableTwoTrackVertexFit=False, enableVertexFitTiming=False):
  process.unpackedTracksAndVertices = unpackedTracksAndVertices
  process.filteredTracks = filteredTracks
  
  process.inclusiveVertexFinder = inclusiveVertexFinder.clone(
        svScoreThreshold=cms.double(threshold_values[0]),
        seedScoreThreshold=cms.double(threshold_values[1]),
        edgeScoreThreshold=cms.double(threshold_values[2])
        )
  #print(f"Using custom thresholds for SV reconstruction: svScoreThreshold={threshold_values[0]}, seedScoreThreshold={threshold_values[1]}, edgeScoreThreshold={threshold_values[2]}")
  #print(f"Using custom values for f1: {f1}, minS1: {minSig1}, f2: {f2}, minS2: {minSig2} ")
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
  process.modelInference = modelInference
  process.trackTable = trackTableFor(isMC)
  process.trackPairFeatures = trackPairFeatures.clone(
      enableTwoTrackVertexFit=cms.bool(enableTwoTrackVertexFit),
      enableVertexFitTiming=cms.bool(enableVertexFitTiming),
  )
  process.trackVertexVars = trackVertexVars

  modules = [
      process.unpackedTracksAndVertices,
      process.filteredTracks,
      # Add variables as value map for tracks
      process.trackVertexVars,
  ]
  if isMC:
    # Match any track to any genparticle in the event (needs customMergedGenParticles)
    process.trackGenMatch = trackGenMatch
    modules.append(process.trackGenMatch)
  modules += [
      # Derive pairfeatures
      process.trackPairFeatures,
      # Run inference on the track graph
      process.modelInference,
      process.trackTable,
  ]
  if isMC:
    # Derive labels for pairs and tracks (needs genVertexTable + trackGenMatch)
    process.trackPairLabels = trackPairLabels
    modules.append(process.trackPairLabels)
  modules += [
      process.inclusiveVertexFinder,
      process.vertexMerger,
      process.trackVertexArbitrator,
      process.myFinalInclusiveSecondaryVertices,
      process.svTable,
  ]
  process.sv_track = _chain(modules)

  if overrides:
    apply_overrides(process, overrides)
  return process


def custom_sv_tracks_training(process, isMC=True):
  process.unpackedTracksAndVertices = unpackedTracksAndVertices
  process.filteredTracks = filteredTracks
  process.trackVertexVars = trackVertexVars
  # no modelInference in this sequence -> drop the GNN score columns too
  process.trackTable = trackTableFor(isMC, withInference=False)

  modules = [
      process.unpackedTracksAndVertices,
      process.filteredTracks,
      process.trackVertexVars,
  ]
  if isMC:
    process.trackGenMatch = trackGenMatch
    modules.append(process.trackGenMatch)
  modules.append(process.trackTable)
  process.sv_track = _chain(modules)
  return process


def custom_train_features(process, gvProducerLabel="genVertexProducer"):
    # MC-only: needs trackGenMatch and the gen-vertex producer
    process.trackGVFeatures = cms.EDProducer("TrackGVFeaturesProducer",
        tracks = cms.InputTag("filteredTracks"),
        primaryVertices = cms.InputTag("filteredTracks"),
        genPartIdx = cms.InputTag("trackGenMatch", "genPartIdx"),
        #secondaryVertices = cms.InputTag("myFinalInclusiveSecondaryVertices"),
        gvTable = cms.InputTag(gvProducerLabel, "GVTable"),
        gvDaughtersTable = cms.InputTag(gvProducerLabel, "GVDaughtersTable"),
        trkPtCut = cms.double(minTrackPt),
        dRMatchMax = cms.double(0.02),
        relPtMatchMax = cms.double(0.2),
        trkMinPt = cms.double(minTrackPt),
        trkMaxEta = cms.double(2.5),
    )
    return process
