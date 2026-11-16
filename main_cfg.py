# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --scenario pp --era Run3_2024 --customise Configuration/DataProcessing/Utils.addMonitoring --step NANO:@JME --conditions 150X_mcRun3_2024_realistic_v2 --datatier NANOAODSIM --eventcontent NANOEDMAODSIM1 --python_filename GEN-RunIII2024Summer24NanoAODv15-00385_1_cfg.py --fileout file:GEN-RunIII2024Summer24NanoAODv15-00385.root --filein dbs:/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM --number 1763 --number_out 1763 --no_exec --mc
#
# Usage:
#   MC   : cmsRun main_cfg.py collection=track isMC=True
#   data : cmsRun main_cfg.py collection=track isMC=False
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
from FWCore.ParameterSet.VarParsing import VarParsing
import sys



process = cms.Process('NANO',Run3_2024)
options = VarParsing('python')

#options.register('outputName', "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "output Name"
#)
options.register('isMC',        True,        VarParsing.multiplicity.singleton,  VarParsing.varType.bool,    "True for MC, False for data")
options.register('globalTag',   "",          VarParsing.multiplicity.singleton,  VarParsing.varType.string,  "Global tag override (empty = default for MC/data)")
options.register('collection', "candidate", VarParsing.multiplicity.singleton,VarParsing.varType.string,    "collection chosen")
options.register('SVLoose',     0.75,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for SV score to enter in IVF")
options.register('SVSeed',      0.85,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for SV score to be a seed for IVF")
options.register('SVEdge',      0.95,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for edgeScore_ij of track j to be included in seed of track i")

options.register('skimEMu',     False,       VarParsing.multiplicity.singleton,  VarParsing.varType.bool,    "Keep only events with an opposite-charge mu-e pair and >=2 jets (loose, ID-free pre-selection for the BTVNanoCommissioning ttdilep e-mu selection); track collection only")
options.register('enableTwoTrackVertexFit', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Fit every selected trackPair row with an unconstrained two-track Kalman vertex fit")
options.register('enableVertexFitTiming', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Save per-event transient-track build and vertex-fit timing counters")
options.register('f1',          0.6,        VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "f1")
options.register('minS1',       4.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "minS1")
options.register('f2',          0.8,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "f2") # 1 to disable, 0.2 default
options.register('minS2',       20.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "minS2") # 0 to disable, 10 default

# Overrides for sv_reco_cff.py parameters ("track" collection). -999 = keep the value written in sv_reco_cff.py
UNSET = -999.
for _name in ['clusterMaxDistance', 'clusterMaxSignificance', 'clusterMinAngleCosine', 'distanceRatio',
              'seedMin3DIPSignificance', 'seedMin3DIPValue',
              'vertexMinAngleCosine', 'vertexMinDLen2DSig', 'vertexMinDLenSig',
              'dlenSigMin']: # dlenSigMin is applied to both svTable and genVertexTable
    options.register(_name, UNSET, VarParsing.multiplicity.singleton, VarParsing.varType.float, _name+" (override of sv_reco_cff.py value)")

options.parseArguments()
isMC = bool(options.isMC)
print("isMC =", isMC, "| collection =", options.collection)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
if isMC:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( (
        #'/store/data/Run2024C/MuonEG/MINIAOD/MINIv6NANOv15-v1/2530000/0058e4ed-6b17-4a79-abe8-afdc9f6d337b.root',
        '/store/mc/RunIII2024Summer24MiniAODv6/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/110000/01af4637-202e-4830-abe8-b129a6fc6a0c.root',
        #'/store/mc/Run3Summer19MiniAOD/TTToHadronic_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2024Scenario_106X_mcRun3_2024_realistic_v4-v2/60000/0173267F-D040-B14E-B16B-0602E6CBC48C.root',
        #"/store/mc/RunIII2024Summer24MiniAODv6/QCD_Bin-PT-15to7000_Par-PT-Flat_TuneCH3_13p6TeV_herwig7/MINIAODSIM/NoPU_150X_mcRun3_2024_realistic_v2_ext2-v2/100000/003bdd8a-afd1-402d-ab83-c300d729b8d9.root"
        #'file:/work/gcelotto/btv_mini_rerun/CMSSW_15_0_6/src/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8.root'
     ) ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--scenario nevts:1763'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

if options.outputFile=="output.root":
    if options.collection=="":
        processName = "outputBTV.root" 
    else:
        #processName = options.collection+".root"
        processName = "BTV_SV_Data24CD_skimmed_splitting2.root"
else:
    processName = options.outputFile
process.NANOEDMAODSIM1output = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM' if isMC else 'NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:'+processName),
    outputCommands = cms.untracked.vstring(
    #'keep *',   
    #'drop *',  # Drop everything by default
    "keep *_svTable_*_*",  # Keep event-level FlatTables
    "keep *_svCandidateTable_*_*",  # Keep event-level FlatTables
    'keep *_genVertexProducer_*_*',
    "keep *_genParticleTable_*_*",  # Keep the merged gen particles
    "keep *_trackTable_*_*",  
    "keep *_trackPairFeatures_*_*",  
    'keep nanoaodFlatTable_modelInference_*_*',     # works
    'keep nanoaodFlatTable_trackPairLabels_*_*',     # works
    'keep nanoaodFlatTable_gvMatchProducer_*_*',     # works
    'keep nanoaodFlatTable_genCentralVertexProducer_*_*',     # works
    'keep nanoaodFlatTable_vertexTable_*_*',     # works
    'keep nanoaodFlatTable_genVertexTable_*_*',     # works
    'keep nanoaodFlatTable_gvMatchTable_*_*',     
    


    #"keep *_gvProducer_*_*",
    #"keep *_gvCentralProducer_*_*",
    #"keep *_genCandidateVertexProducer_*_*",
    #"keep nanoaodFlatTable_*Table*_*_*",  # Keep event-level FlatTables
    'keep nanoaodFlatTable_electronTable_*_*',     # Electron_*
    'keep nanoaodFlatTable_muonTable_*_*',         # Muon_*
    'keep nanoaodFlatTable_jetPuppiTable_*_*', 
    #'keep nanoaodFlatTable_customGenParticleTable_*_*',


    'keep edmTriggerResults_*_*_HLT', # Trigger? 

    'keep nanoaodFlatTable_puTable_*_*',


    "keep nanoaodUniqueString_nanoMetadata_*_*",  # Keep basic metadata
    "keep nanoaodMergeableCounterTable_*_*_*",
    'keep nanoaodFlatTable_genWeightsTable_*_*'


)
)


if isMC:
    process.NANOEDMAODSIM1output.outputCommands.append('keep nanoaodFlatTable_customGenParticleTable_*_*')

# Additional output definition

# Other statements
# Global tag: MC uses the GT of the Summer24 MiniAODv6 production.
# For data, replace 'auto:run3_data' with the GT read from DAS:
#   dasgoclient -query="config dataset=/MuonEG/Run2024C-MINIv6NANOv15-v1/MINIAOD" | grep -i globaltag
# or pass it on the command line with globalTag=<GT>.
MC_GT   = '150X_mcRun3_2024_realistic_v2'
DATA_GT = 'auto:run3_data'
if options.globalTag:
    gt = options.globalTag
else:
    gt = MC_GT if isMC else DATA_GT
print("Using global tag:", gt)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, gt, '')



sys.path.append("/work/gcelotto/btv_new/CMSSW_15_0_6/src/PhysicsTools/SimpleVertexTable/python")
from sv_reco_cff import custom_sv_tracks, custom_train_features
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
process.load('PhysicsTools.NanoAOD.nano_cff')
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon

process = nanoAOD_customizeCommon(process)


#. 
#.    MERGED GEN PARTICLES + GEN VERTEX PRODUCER  (MC only)
#. 

if isMC:
    from mergedGenParticles import customMergedGenParticles, customFinalGenParticles, customGenParticleTable
    process.customMergedGenParticles = customMergedGenParticles 
    process.customFinalGenParticles = customFinalGenParticles 
    process.customGenParticleTable = customGenParticleTable 

    from genVertices import custom_GV_producer
    process = custom_GV_producer(process, collection=options.collection)

    #print("\n===== GV PRODUCER DEBUG =====")
    #print(process.gvProducer.dumpPython())
    #print("================================\n")


#.
#.    e-mu SKIM (optional)
#.
# Looser than BTVNanoCommissioning ttdilep_valid_sf (lepton ID/iso, "exactly one" lepton, jet ID and
# lepton-jet cleaning are left to the offline selection). No trigger requirement: data is MuonEG already.

if options.skimEMu:
    if options.collection != "track":
        raise RuntimeError("skimEMu is only implemented for collection=track")
    process.emuSkimMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("slimmedMuons"),
        cut = cms.string("pt > 25 && abs(eta) < 2.4"),
        filter = cms.bool(True))
    process.emuSkimElectrons = cms.EDFilter("PATElectronSelector",
        src = cms.InputTag("slimmedElectrons"),
        cut = cms.string("pt > 25 && abs(eta) < 2.5"),
        filter = cms.bool(True))
    #process.emuSkimPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    #    decay = cms.string("emuSkimMuons emuSkimElectrons"),
    #    cut = cms.string("daughter(0).charge * daughter(1).charge < 0"))     # opposite charge
    #process.emuSkimPairFilter = cms.EDFilter("CandViewCountFilter",
    #    src = cms.InputTag("emuSkimPairs"),
    #    minNumber = cms.uint32(1))
    process.emuSkimJets = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("slimmedJetsPuppi"),
        cut = cms.string("pt > 20 && abs(eta) < 2.5"),
        filter = cms.bool(True))
    process.emuSkimJetFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag("emuSkimJets"),
        minNumber = cms.uint32(2))
    process.emuSkim = cms.Sequence(
        process.emuSkimMuons * process.emuSkimElectrons * #process.emuSkimPairs * process.emuSkimPairFilter *
        process.emuSkimJets * process.emuSkimJetFilter)
    # the muon/electron/jet selectors filter on >=1 object, the count filters on the pair and >=2 jets
    process.NANOEDMAODSIM1output.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('nanoAOD_step'))


#.
#.    PATHS
#.

if options.collection=="candidate":
    if not isMC:
        raise RuntimeError("collection=candidate is MC-only (needs gen particles / gen vertices)")
    from sv_candidate_cff import custom_sv_candidate #inclusiveCandidateVertexFinder, candidateVertexMerger, CandidateVertexArbitrator, myCandidateInclusiveSecondaryVertices, myvertexTable, svCandidateTable, slimmedSecondaryVertices
    process = custom_sv_candidate(process)
    process.nanoAOD_step = cms.Path(
                                process.customMergedGenParticles*
                                process.customFinalGenParticles*
                                process.customGenParticleTable+
                                process.sv_candidate*
                                process.genVertexProducer_sequence )

elif options.collection=="track":
    process = custom_sv_tracks(process, threshold_values=(options.SVLoose, options.SVSeed, options.SVEdge), f1=options.f1, minSig1=options.minS1, f2=options.f2, minSig2=options.minS2,
                                 overrides={k: getattr(options, k) for k in ['clusterMaxDistance', 'clusterMaxSignificance', 'clusterMinAngleCosine', 'distanceRatio', 'seedMin3DIPSignificance', 'seedMin3DIPValue', 'vertexMinAngleCosine', 'vertexMinDLen2DSig', 'vertexMinDLenSig', 'dlenSigMin'] if getattr(options, k) != UNSET},
                                 isMC=isMC, enableTwoTrackVertexFit=options.enableTwoTrackVertexFit,
                                 enableVertexFitTiming=options.enableVertexFitTiming)
    # genWeightsTable keeps the run/lumi-level sum of weights: it has to see every event, so it stays before the skim
    skim = process.emuSkim if options.skimEMu else cms.Sequence()
    if isMC:
        process = custom_train_features(process, gvProducerLabel="gvProducer")
        process.nanoAOD_step = cms.Path(
            process.genWeightsTable+
            skim +
            process.puTable +
            process.customMergedGenParticles *
            process.customFinalGenParticles *
            process.customGenParticleTable +
            process.genVertexTable *     # PURELY GV TABLE
            process.sv_track *
            process.gvMatchTable
        )
    else:
        process.nanoAOD_step = cms.Path(skim + process.sv_track)

elif options.collection=="central":
    #process.nanoTableTaskCommon = cms.Task(process.boostedTauTask, process.electronTask, process.jetForMETTask, process.jetPuppiForMETTask, process.jetPuppiTablesTask, process.jetPuppiTask, process.jetTablesTask, process.jetTask, process.linkedObjects, process.lowPtElectronTablesTask, process.lowPtElectronTask, process.metTablesTask, process.muonTablesTask, process.muonTask, process.nanoMetadata, process.photonTask, process.svCandidateTable, process.tauTablesTask, process.tauTask, process.vertexTable, process.vertexTask)
    #process.nanoSequenceCommon = cms.Sequence(process.nanoTableTaskCommon)
    #process.nanoTableTaskFS = cms.Task()
    if isMC:
        process.nanoAOD_step = cms.Path(
                                    process.customMergedGenParticles*
                                    process.customFinalGenParticles*
                                    process.customGenParticleTable+
                                    process.nanoSequenceMC*
                                    process.genCentralVertexProducer 
                                    )
    else:
        process.nanoAOD_step = cms.Path(process.nanoSequence)
else:
    raise RuntimeError("Unknown collection '%s' (use candidate, track or central)" % options.collection)

# Add electrons and muons and jets
process.nanoAOD_step.associate(process.nanoTableTaskCommon)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOEDMAODSIM1output_step = cms.EndPath(process.NANOEDMAODSIM1output)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_step,    # central nanoAOD step already included in nanoSequenceMC
                                process.endjob_step,
                                process.NANOEDMAODSIM1output_step)



from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
#from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon 

#call to customisation function nanoAOD_customizeCommon imported from PhysicsTools.NanoAOD.nano_cff
#process = nanoAOD_customizeCommon(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.custom_jme_cff
from PhysicsTools.NanoAOD.custom_jme_cff import PrepJMECustomNanoAOD 

#call to customisation function PrepJMECustomNanoAOD imported from PhysicsTools.NanoAOD.custom_jme_cff
#process = PrepJMECustomNanoAOD(process)

# End of customisation functions


# Customisation from command line

process.source.delayReadingEventProducts = cms.untracked.bool(False)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
