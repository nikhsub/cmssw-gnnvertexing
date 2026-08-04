import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4V1_cff import (
    pfUnifiedParticleTransformerAK4V1JetTags,
    pfUnifiedParticleTransformerAK4V1DiscriminatorsJetTags,
    _pfUnifiedParticleTransformerAK4V1JetTagsAll,
)


def add_upart_v1_inference(
    process,
    jet_source="slimmedJetsPuppi",
    postfix="UParT",
    table_name="UParTJet",
    min_jet_pt=2.0,
    max_jet_eta=2.5,
):
    """
    Add UParT AK4 V1 inference to a MiniAOD/MiniAODSIM process.

    The configured model is:
      RecoBTag/Combined/data/UParTAK4/PUPPI/V00/UParTAK4.onnx

    Parameters
    ----------
    process:
        cms.Process being configured.

    jet_source:
        Existing MiniAOD AK4 PUPPI PAT jet collection.

    postfix:
        Postfix used by updateJetCollection when cloning modules.

    table_name:
        NanoAOD FlatTable name. Output branches are
        <table_name>_<variable>.

    min_jet_pt:
        Minimum updated jet pT stored in the table.

    max_jet_eta:
        Maximum absolute updated jet eta stored in the table.

    Returns
    -------
    process:
        Updated cms.Process.
    """

    # ---------------------------------------------------------
    # Rerun UParT through the standard PAT jet update machinery.
    # ---------------------------------------------------------
    updateJetCollection(
        process,
        jetSource=cms.InputTag(jet_source),
        pvSource=cms.InputTag("offlineSlimmedPrimaryVertices"),
        pfCandidates=cms.InputTag("packedPFCandidates"),
        svSource=cms.InputTag("slimmedSecondaryVertices"),
        muSource=cms.InputTag("slimmedMuons"),
        elSource=cms.InputTag("slimmedElectrons"),
        jetCorrections=(
            "AK4PFPuppi",
            cms.vstring(
                "L2Relative",
                "L3Absolute",
            ),
            "None",
        ),
        btagDiscriminators=
            _pfUnifiedParticleTransformerAK4V1JetTagsAll,
        postfix=postfix,
        printWarning=False,
    )

    updated_jets = "selectedUpdatedPatJets" + postfix
    raw_tagger = (
        "pfUnifiedParticleTransformerAK4V1JetTags"
    )
    
    discriminator_tagger = (
        "pfUnifiedParticleTransformerAK4V1"
        "DiscriminatorsJetTags"
    )
    

    # ---------------------------------------------------------
    # Basic jet variables and generator-truth flavor labels.
    # ---------------------------------------------------------
    variables = cms.PSet(
        pt=Var(
            "pt",
            float,
            precision=10,
            doc="Updated AK4 PUPPI jet transverse momentum",
        ),
        eta=Var(
            "eta",
            float,
            precision=10,
            doc="Updated AK4 PUPPI jet pseudorapidity",
        ),
        phi=Var(
            "phi",
            float,
            precision=10,
            doc="Updated AK4 PUPPI jet azimuth",
        ),
        mass=Var(
            "mass",
            float,
            precision=10,
            doc="Updated AK4 PUPPI jet mass",
        ),

        # Standard PAT generator-flavor labels.
        #
        # hadronFlavour:
        #   5 = b jet
        #   4 = c jet
        #   0 = no b/c hadron, usually light flavor
        hadronFlavour=Var(
            "hadronFlavour",
            int,
            doc=(
                "Generator-level hadron flavor: "
                "5=b, 4=c, 0=light/no heavy hadron"
            ),
        ),

        # partonFlavour normally stores the signed PDG ID of
        # the associated initiating parton.
        partonFlavour=Var(
            "partonFlavour",
            int,
            doc=(
                "Generator-level parton flavor PDG ID; "
                "21=gluon, abs(1-3)=uds, abs(4)=c, abs(5)=b"
            ),
        ),

        # Convenient boolean truth labels.
        isB=Var(
            "hadronFlavour == 5",
            bool,
            doc="Jet contains a b hadron",
        ),
        isC=Var(
            "hadronFlavour == 4",
            bool,
            doc="Jet contains a c hadron and no b hadron",
        ),
        isLight=Var(
            "hadronFlavour == 0",
            bool,
            doc="Jet has no matched b or c hadron",
        ),
        isGluon=Var(
            "hadronFlavour == 0 && abs(partonFlavour) == 21",
            bool,
            doc="Light-flavor jet matched to a gluon",
        ),
        isUDS=Var(
            (
                "hadronFlavour == 0 && "
                "abs(partonFlavour) >= 1 && "
                "abs(partonFlavour) <= 3"
            ),
            bool,
            doc="Light-flavor jet matched to a u, d, or s quark",
        ),
    )

    # ---------------------------------------------------------
    # Save every raw UParT output class.
    # ---------------------------------------------------------
    for flav_name in (
        pfUnifiedParticleTransformerAK4V1JetTags.flav_names
    ):
        setattr(
            variables,
            flav_name,
            Var(
                "bDiscriminator('{}:{}')".format(
                    raw_tagger,
                    flav_name,
                ),
                float,
                precision=10,
                doc="Raw UParT V1 probability {}".format(
                    flav_name
                ),
            ),
        )
    # ---------------------------------------------------------
    # Save every combined UParT discriminator defined by CMSSW.
    # These normally include BvsAll, CvB, CvL, QvsG, etc.,
    # depending on the exact CMSSW model configuration.
    # ---------------------------------------------------------
    for discriminator in (
        pfUnifiedParticleTransformerAK4V1DiscriminatorsJetTags
        .discriminators
    ):
        discriminator_name = discriminator.name.value()

        # FlatTable variable names should be simple identifiers.
        branch_name = (
            discriminator_name
            .replace(":", "_")
            .replace("-", "_")
            .replace(".", "_")
        )

        setattr(
            variables,
            branch_name,
            Var(
                "bDiscriminator('{}:{}')".format(
                    discriminator_tagger,
                    discriminator_name,
                ),
                float,
                precision=10,
                doc="Combined UParT discriminator {}".format(
                    discriminator_name
                ),
            ),
        )

    # ---------------------------------------------------------
    # NanoAOD FlatTable.
    # ---------------------------------------------------------
    process.upartJetTable = cms.EDProducer(
        "SimplePATJetFlatTableProducer",
        src=cms.InputTag(updated_jets),
        cut=cms.string(
            "pt > {} && abs(eta) < {}".format(
                min_jet_pt,
                max_jet_eta,
            )
        ),
        name=cms.string(table_name),
        doc=cms.string(
            "AK4 PUPPI jets with UParT V1 outputs and MC truth"
        ),
        singleton=cms.bool(False),
        extension=cms.bool(False),
        variables=variables,
    )

    _patAlgosToolsTask = getPatAlgosToolsTask(process)

    process.upartPath = cms.Path(
        process.upartJetTable
    )

    process.upartPath.associate(
        _patAlgosToolsTask
    )

    return process 

    return process
