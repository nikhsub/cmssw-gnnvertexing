#!/usr/bin/env python3
"""Synthetic tests for Stage 3 event-local truth association and counts."""

import copy
import unittest

from analyze_hzzc import _metric_rows, analyze_event
from neighborhood_analysis import direction_from_tracks, delta_r
from truth_objects import GEN_FIELDS, HADRON_FIELDS, QUARK_FIELDS, RECO_FIELDS, accepted, n_category, validate_tables


def columns(rows, fields):
    return {key: [row[key] for row in rows] for key in fields}


def synthetic_event():
    qfrag = [
        dict(pdgId=4, pt=12., eta=.1, phi=.2, mass=0., prunedGenIdx=40, hardPrunedIdx=3),
        dict(pdgId=5, pt=20., eta=-.2, phi=-.3, mass=0., prunedGenIdx=50, hardPrunedIdx=5),
    ]
    hadrons = [
        dict(prunedGenIdx=100, pdgId=421, pt=3., eta=.1, phi=.2, mass=1.86, charmOrigin=1, parentBHadronIdx=-1, quarkIdx=0, associationQuality=1, fragmentationSourceCertain=1),
        dict(prunedGenIdx=101, pdgId=411, pt=1.2, eta=.2, phi=.25, mass=1.87, charmOrigin=1, parentBHadronIdx=-1, quarkIdx=0, associationQuality=1, fragmentationSourceCertain=1),
        dict(prunedGenIdx=102, pdgId=431, pt=4., eta=.15, phi=.3, mass=1.97, charmOrigin=1, parentBHadronIdx=-1, quarkIdx=0, associationQuality=1, fragmentationSourceCertain=1),
        dict(prunedGenIdx=103, pdgId=421, pt=2., eta=.4, phi=.4, mass=1.86, charmOrigin=3, parentBHadronIdx=-1, quarkIdx=0, associationQuality=1, fragmentationSourceCertain=0),
        dict(prunedGenIdx=104, pdgId=511, pt=8., eta=-.2, phi=-.3, mass=5.28, charmOrigin=0, parentBHadronIdx=-1, quarkIdx=1, associationQuality=1, fragmentationSourceCertain=1),
        dict(prunedGenIdx=105, pdgId=421, pt=2.5, eta=-.1, phi=-.2, mass=1.86, charmOrigin=2, parentBHadronIdx=4, quarkIdx=1, associationQuality=1, fragmentationSourceCertain=1),
    ]
    for row, qidx in zip(hadrons, (0, 0, 0, 0, 1, 1)):
        q = qfrag[qidx]
        row.update(preHadPt=q["pt"], preHadEta=q["eta"], preHadPhi=q["phi"], preHadKinematicsValid=1)
    values = []
    def gen(pdg, charge, pt, eta, phi, hadron, quark, membership, quality=1, certain=0, from_b=0, prompt=0, from_b_d=0):
        values.append(dict(packedGenIdx=len(values), pdgId=pdg, charge=charge, pt=pt, eta=eta, phi=phi, mass=.14,
                           lastPrunedIdx=-1, lastPrunedIsSelf=0, hadronIdx=hadron, quarkIdx=quark,
                           associationQuality=quality, membershipType=membership, fragmentationSourceCertain=certain,
                           isFromB=from_b, isFromPromptD=prompt, isFromBtoD=from_b_d))
    # D0 -> K pi, plus two quark-associated candidates. Extra tracks must not
    # change the D from N2.
    gen(321, -1, .8, .1, .20, 0, 0, 1, certain=1, prompt=1)
    gen(211, 1, .9, .11, .21, 0, 0, 1, certain=1, prompt=1)
    gen(211, 1, .6, .12, .22, -1, 0, 2, certain=0)
    gen(321, -1, .7, .09, .19, -1, 0, 2, certain=0)
    # A boundary track at pt=0.4, |eta|=2.5 is accepted; the other daughter fails.
    gen(321, -1, .4, 2.5, .5, 1, 0, 1, certain=1, prompt=1)
    gen(211, 1, .39, 0., .5, 1, 0, 1, certain=1, prompt=1)
    # Prompt D with exactly three accepted matched decay tracks.
    gen(321, -1, .8, .15, .3, 2, 0, 1, certain=1, prompt=1)
    gen(211, 1, .7, .16, .31, 2, 0, 1, certain=1, prompt=1)
    gen(211, 1, .6, .17, .32, 2, 0, 1, certain=1, prompt=1)
    # Unknown-origin D remains separate despite a unique ancestry link.
    gen(211, 1, .8, .4, .4, 3, 0, 1, certain=0)
    gen(211, 1, 1., -.2, -.3, 4, 1, 1, certain=1, from_b=1)
    # B -> D cascade; these tracks belong to the D row, not to a separate N.
    gen(321, -1, .9, -.1, -.2, 5, 1, 1, certain=1, from_b=1, from_b_d=1)
    gen(211, 1, .8, -.12, -.22, 5, 1, 1, certain=1, from_b=1, from_b_d=1)
    # One unrelated generator track plus an unmatched reconstructed track.
    gen(211, 1, 1.1, 1.8, 2.1, -1, -1, 0, quality=0)
    gen(211, -1, .9, -1.8, -2.1, -1, -1, 3, quality=2)
    recos = []
    def reco(gen_idx, pt, eta, phi, charge, status=1, high_purity=1):
        recos.append(dict(trackIdx=len(recos), pt=pt, eta=eta, phi=phi, charge=charge,
                          isHighPurity=high_purity, genChargedIdx=gen_idx, matchDR=.001,
                          matchRelPt=.01, matchStatus=status, nGenCandidates=1 if status == 1 else 2))
    for gen_idx in (0, 1, 2, 3):
        g = values[gen_idx]; reco(gen_idx, g["pt"], g["eta"], g["phi"], g["charge"])
    for gen_idx in (4, 6, 7, 8, 9, 10, 11, 12, 14):
        g = values[gen_idx]; reco(gen_idx, g["pt"], g["eta"], g["phi"], g["charge"])
    reco(-1, 1.2, 1.8, 2.2, 1, status=0, high_purity=0)
    return {
        "QFrag": columns(qfrag, ("pdgId", "pt", "eta", "phi", "mass", "prunedGenIdx", "hardPrunedIdx")),
        "QHadron": columns(hadrons, HADRON_FIELDS),
        "QGenCharged": columns(values, ("packedGenIdx", "pdgId", "charge", "pt", "eta", "phi", "mass", "lastPrunedIdx", "lastPrunedIsSelf", "hadronIdx", "quarkIdx", "associationQuality", "membershipType", "fragmentationSourceCertain", "isFromB", "isFromPromptD", "isFromBtoD")),
        "QRecoTrack": columns(recos, ("trackIdx", "pt", "eta", "phi", "charge", "isHighPurity", "genChargedIdx", "matchDR", "matchRelPt", "matchStatus", "nGenCandidates")),
    }


class Stage3SyntheticTest(unittest.TestCase):
    def setUp(self):
        self.tables = synthetic_event()
        self.event = {"run": 1, "luminosityBlock": 2, "event": 3}

    def test_inclusive_acceptance_and_multiplicity_categories(self):
        self.assertTrue(accepted(.4, -2.5))
        self.assertFalse(accepted(.399999, 0.))
        self.assertFalse(accepted(.8, 2.50001))
        self.assertEqual([n_category(x) for x in (0, 1, 2, 3, 4)], ["N0", "N1", "N2", "N3", "N3+"])

    def test_n2_ignores_extra_fragmentation_tracks_and_shared_qfrag(self):
        hadrons, extras, neighborhoods, metrics, _ = analyze_event(self.tables, self.event, "synthetic.root")
        d0 = next(r for r in hadrons if r["QHadron_row"] == 0)
        self.assertEqual(d0["N_category"], "N2")
        self.assertEqual(metrics["promptD_N2_accepted_additional_reco_tracks"], 2)
        self.assertEqual(len([r for r in extras if r["QHadron_row"] == 0]), 2)
        self.assertEqual(hadrons[0]["QFrag_row"], hadrons[1]["QFrag_row"])
        self.assertEqual(metrics["shared_prompt_D_systems"], 1)
        self.assertEqual(metrics["prompt_Ds_in_shared_prompt_D_systems"], 3)
        self.assertTrue(all(r["fragmentation_source_certain"] == 0 for r in extras))
        self.assertTrue(any(r["neighborhood_radius"] == "all" for r in neighborhoods))

    def test_n1_n3_unknown_origin_and_b_to_d_are_separate(self):
        hadrons, _, _, metrics, _ = analyze_event(self.tables, self.event, "synthetic.root")
        by_row = {r["QHadron_row"]: r for r in hadrons}
        self.assertEqual(by_row[1]["N_category"], "N1")
        self.assertEqual(by_row[2]["N_category"], "N3")
        self.assertEqual(metrics["hadron_D_prompt_N3"], 1)
        self.assertEqual(metrics["hadron_D_prompt_N3plus"], 1)
        self.assertEqual(metrics["promptD_N3_objects"], 1)
        self.assertEqual(metrics["promptD_N3+_objects"], 1)
        summary = _metric_rows(metrics, [], 1)
        n3plus = next(r for r in summary if r["metric"] == "D_prompt_N3+_objects")
        self.assertEqual(n3plus["value"], 1)
        self.assertEqual(by_row[3]["origin_category"], "D_unknown")
        self.assertEqual(by_row[5]["origin_category"], "D_fromB")
        self.assertEqual(by_row[5]["parentBHadronIdx"], 4)

    def test_auxiliary_targets_require_prompt_charm_and_preserve_diagnostics(self):
        hadrons, _, _, _, _ = analyze_event(self.tables, self.event, "synthetic.root")
        by_row = {r["QHadron_row"]: r for r in hadrons}
        d0 = by_row[0]
        self.assertEqual(d0["auxiliary_target_eligible"], 1)
        self.assertAlmostEqual(d0["zT"], 0.25)
        self.assertAlmostEqual(d0["deltaR_DQuark"], 0.0)
        self.assertEqual(d0["accepted_additional_reco_candidate_count"], 2)
        self.assertEqual(by_row[5]["observable_interpretation"], "D_vs_B_parent_lineage_not_B_fragmentation")
        self.assertEqual(by_row[5]["auxiliary_target_eligible"], 0)
        self.assertEqual(by_row[3]["observable_interpretation"], "unknown_origin_diagnostic_only")
        self.assertEqual(by_row[3]["auxiliary_target_eligible"], 0)
        damaged = copy.deepcopy(self.tables)
        damaged["QHadron"]["preHadKinematicsValid"][0] = 0
        damaged_rows = analyze_event(damaged, self.event, "synthetic.root")[0]
        self.assertEqual(damaged_rows[0]["auxiliary_target_eligible"], 0)
        self.assertEqual(damaged_rows[0]["zT_valid"], 0)

    def test_event_local_indices_and_unique_track_counting(self):
        validate_tables(self.tables)
        damaged = copy.deepcopy(self.tables)
        damaged["QGenCharged"]["hadronIdx"][0] = 99
        with self.assertRaisesRegex(ValueError, "invalid event-local QHadron"):
            validate_tables(damaged)
        damaged = copy.deepcopy(self.tables)
        damaged["QRecoTrack"]["genChargedIdx"][1] = 0
        with self.assertRaisesRegex(ValueError, "assigns QGenCharged row .* more than once"):
            validate_tables(damaged)

    def test_empty_event_collections_are_safe(self):
        empty = {
            "QFrag": {field: [] for field in QUARK_FIELDS},
            "QHadron": {field: [] for field in HADRON_FIELDS},
            "QGenCharged": {field: [] for field in GEN_FIELDS},
            "QRecoTrack": {field: [] for field in RECO_FIELDS},
        }
        hadrons, extras, neighborhoods, metrics, graph = analyze_event(empty, self.event, "empty.root")
        self.assertEqual(hadrons, [])
        self.assertEqual(extras, [])
        self.assertEqual(neighborhoods, [])
        self.assertEqual(metrics["QHadron_rows"], 0)
        self.assertEqual(graph["nodes"], 0)
        self.assertEqual(graph["all_pairs"], 0)

    def test_unresolved_tracks_not_confident_negatives(self):
        _, _, _, metrics, _ = analyze_event(self.tables, self.event, "synthetic.root")
        # Unknown-origin D remains a distinct origin population but its retained
        # hadron-decay link is usable for direct hadron identity labels. The
        # separate membershipType=3 track and unmatched reco track are excluded.
        self.assertEqual(metrics["OptionA_positive_same_hadron_track_pairs"], 5)
        self.assertEqual(metrics["OptionA_confident_different_hadron_negative_pairs"], 40)
        self.assertEqual(metrics["OptionA_uncertain_additional_reco_tracks"], 2)
        self.assertEqual(metrics["OptionA_D_prompt_N1_objects_without_pair_label"], 1)

    def test_high_purity_is_optional_and_direction_helpers_are_safe(self):
        standard = analyze_event(self.tables, self.event, "synthetic.root")
        high_purity = analyze_event(self.tables, self.event, "synthetic.root", require_high_purity=True)
        self.assertEqual(standard[0][0]["N_category"], "N2")
        self.assertEqual(high_purity[0][0]["N_category"], "N2")
        eta, phi = direction_from_tracks([{"pt": 1., "eta": 0., "phi": 0.}])
        self.assertAlmostEqual(eta, 0.)
        self.assertAlmostEqual(phi, 0.)
        self.assertAlmostEqual(float(delta_r(0., 0., 0., 0.)), 0.)


if __name__ == "__main__":
    unittest.main()
