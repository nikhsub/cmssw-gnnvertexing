import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import *

mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                            inputPruned = cms.InputTag("prunedGenParticles"),
                                            inputPacked = cms.InputTag("packedGenParticles"),
                                            )




finalGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("mergedGenParticles"),
    select = cms.vstring(
	"drop *",
        "keep++ abs(pdgId) == 2212  ",#  keep full tau decay chain for some taus
   )
)

##################### Tables for final output and docs ##########################
genParticleTable = cms.EDProducer("SimpleGenParticleFlatTableProducer",
    src = cms.InputTag("finalGenParticles"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("GenPart"),
    doc = cms.string("interesting gen particles "),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the taus
    variables = cms.PSet(
         pt  = Var("pt",  float, precision=8),
         phi = Var("phi", float,precision=8),
         eta  = Var("eta",  float,precision=8),
         #mass = Var("?mass>10 || (pdgId==22 && mass > 1) || abs(pdgId)==24 || pdgId==23?mass:0", float,precision="?8",doc="Mass stored for all particles with mass > 10 GeV and photons with mass > 1 GeV. For other particles you can lookup from PDGID"),
         pdgId  = Var("pdgId", int, doc="PDG id"),
         status  = Var("status", int, doc="Particle status. 1=stable"),
         genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", "int16", doc="index of the mother particle"),
         genPartIdxMother_2 = Var("?numberOfMothers>1?motherRef(1).key():-1", "int16", doc="index of the mother particle"),
         vx  = Var("vx", float, precision=14, doc="x coordinate of the production vertex"),
         vy  = Var("vy", float, precision=14, doc="x coordinate of the production vertex"),
         vz  = Var("vz", float, precision=14, doc="x coordinate of the production vertex"),
    )
)


