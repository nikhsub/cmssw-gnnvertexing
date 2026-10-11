#!/usr/bin/env python3
"""Compare approved pre-cleanup tables with the compact producer output.

Run in a CMSSW runtime where PyROOT is available:
  python3 validate_quark_fragmentation_cleanup.py approved.root cleaned.root
"""
import argparse
import math
import sys
from collections import Counter

import ROOT

TREE = "Events"
# (new field suffix, approved post-ancestry-boundary field suffix, numeric kind)
SCHEMA = {
    "QFrag": [
        ("pdgId", "pdgId", "i"), ("prunedGenIdx", "prunedGenIdx", "i"),
        ("pt", "pt", "f"), ("eta", "eta", "f"), ("phi", "phi", "f"), ("mass", "mass", "f"),
        ("hardPrunedIdx", "hardQuarkIdx", "i"), ("hardPt", "hardPt", "f"),
        ("hardEta", "hardEta", "f"), ("hardPhi", "hardPhi", "f"),
    ],
    "QHadron": [
        ("prunedGenIdx", "prunedGenIdx", "i"), ("pdgId", "pdgId", "i"),
        ("pt", "pt", "f"), ("eta", "eta", "f"), ("phi", "phi", "f"), ("mass", "mass", "f"),
        ("charmOrigin", "charmOrigin", "i"), ("parentBHadronIdx", "parentBHadronIdx", "i"),
        ("quarkIdx", "lineageV2QuarkIdx", "i"),
        ("associationQuality", "lineageV2AssociationQuality", "i"),
        ("fragmentationSourceCertain", "lineageV2FragmentationSourceCertain", "i"),
        ("preHadPrunedIdx", "lineageV2PreHadPrunedIdx", "i"),
        ("preHadPt", "lineageV2PreHadPt", "f"), ("preHadEta", "lineageV2PreHadEta", "f"),
        ("preHadPhi", "lineageV2PreHadPhi", "f"),
        ("preHadKinematicsValid", "lineageV2PreHadKinematicsValid", "i"),
    ],
    "QGenCharged": [
        ("packedGenIdx", "packedGenIdx", "i"), ("pdgId", "pdgId", "i"), ("charge", "charge", "i"),
        ("pt", "pt", "f"), ("eta", "eta", "f"), ("phi", "phi", "f"), ("mass", "mass", "f"),
        ("lastPrunedIdx", "lastPrunedIdx", "i"), ("lastPrunedIsSelf", "lastPrunedIsSelf", "i"),
        ("hadronIdx", "v2HadronIdx", "i"), ("quarkIdx", "v2QuarkIdx", "i"),
        ("associationQuality", "v2AssociationQuality", "i"),
        ("membershipType", "v2MembershipType", "i"),
        ("fragmentationSourceCertain", "v2FragmentationSourceCertain", "i"),
        ("isFromB", "isFromB", "i"), ("isFromPromptD", "isFromPromptD", "i"),
        ("isFromBtoD", "isFromBtoD", "i"),
    ],
}
KEYS = ("run", "luminosityBlock", "event")
FLOAT_ATOL = 2.0e-6
FLOAT_RTOL = 2.0e-6


def fail(message):
    raise AssertionError(message)


def open_tree(path):
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        fail("cannot open ROOT file: " + path)
    tree = root_file.Get(TREE)
    if not tree:
        fail("missing Events tree in " + path)
    return root_file, tree


def branch_names(tree):
    return {branch.GetName() for branch in tree.GetListOfBranches()}


def read_rows(tree, is_reference):
    required = set(KEYS)
    expected_cleaned_table_branches = set()
    for table, fields in SCHEMA.items():
        for new_name, old_name, _ in fields:
            required.add(table + "_" + (old_name if is_reference else new_name))
            if not is_reference:
                expected_cleaned_table_branches.add(table + "_" + new_name)
    missing = required - branch_names(tree)
    if missing:
        fail("missing branches: " + ", ".join(sorted(missing)))
    if not is_reference:
        prefixes = tuple(table + "_" for table in SCHEMA)
        actual = {name for name in branch_names(tree) if name.startswith(prefixes)}
        extra = actual - expected_cleaned_table_branches
        if extra:
            fail("cleaned output contains branches outside the compact table schema: " + ", ".join(sorted(extra)))
    rows = {}
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        key = tuple(int(getattr(tree, name)) for name in KEYS)
        if key in rows:
            fail("duplicate event key %r" % (key,))
        event = {}
        for table, fields in SCHEMA.items():
            event[table] = {}
            for new_name, old_name, _ in fields:
                branch = table + "_" + (old_name if is_reference else new_name)
                event[table][new_name] = list(getattr(tree, branch))
        rows[key] = event
    return rows


def compare(reference, cleaned):
    if reference.keys() != cleaned.keys():
        missing = sorted(reference.keys() - cleaned.keys())[:5]
        extra = sorted(cleaned.keys() - reference.keys())[:5]
        fail("event keys differ; missing=%r extra=%r" % (missing, extra))
    stats = Counter()
    for key in sorted(reference):
        before, after = reference[key], cleaned[key]
        for table, fields in SCHEMA.items():
            for new_name, _, kind in fields:
                a, b = before[table][new_name], after[table][new_name]
                if len(a) != len(b):
                    fail("%s.%s row count differs in event %r: %d != %d" % (table, new_name, key, len(a), len(b)))
                for row, (old, new) in enumerate(zip(a, b)):
                    if kind == "i":
                        if int(old) != int(new):
                            fail("%s.%s differs at event %r row %d: %r != %r" % (table, new_name, key, row, old, new))
                    elif not math.isclose(float(old), float(new), rel_tol=FLOAT_RTOL, abs_tol=FLOAT_ATOL):
                        fail("%s.%s differs at event %r row %d: %.9g != %.9g" % (table, new_name, key, row, old, new))
                # Check every event-local reference, including sentinel handling.
        nq, nh = len(after["QFrag"]["pdgId"]), len(after["QHadron"]["pdgId"])
        for field in ("quarkIdx",):
            for table in ("QHadron", "QGenCharged"):
                for value in after[table][field]:
                    if value < -1 or value >= nq:
                        fail("invalid %s.%s=%d in event %r" % (table, field, value, key))
        for value in after["QHadron"]["parentBHadronIdx"]:
            if value < -1 or value >= nh:
                fail("invalid QHadron.parentBHadronIdx=%d in event %r" % (value, key))
        for value in after["QGenCharged"]["hadronIdx"]:
            if value < -1 or value >= nh:
                fail("invalid QGenCharged.hadronIdx=%d in event %r" % (value, key))
        charged = after["QGenCharged"]
        for row, (quality, membership, quark, certain) in enumerate(zip(
                charged["associationQuality"], charged["membershipType"], charged["quarkIdx"],
                charged["fragmentationSourceCertain"])):
            if membership == 2 and (quality != 1 or quark < 0 or certain != 0):
                fail("additional candidate violates unique-ancestry/source-uncertain policy in event %r row %d" % (key, row))
        stats["events"] += 1
        stats["QFrag rows"] += nq
        stats["QHadron rows"] += nh
        stats["QHadron unique ancestry assignments"] += sum(x == 1 for x in after["QHadron"]["associationQuality"])
        stats["QGenCharged rows"] += len(charged["pdgId"])
        stats["heavy decay charged rows"] += sum(x == 1 for x in charged["membershipType"])
        stats["unique additional candidates"] += sum(q == 1 and m == 2 for q, m in zip(charged["associationQuality"], charged["membershipType"]))
        # Accepted prompt-D and extra-candidate checks match the stated validation cuts.
        accepted_decay_by_hadron = Counter()
        for hadron, membership, pt, eta in zip(charged["hadronIdx"], charged["membershipType"], charged["pt"], charged["eta"]):
            if hadron >= 0 and membership == 1 and pt >= 0.4 and abs(eta) <= 2.5:
                accepted_decay_by_hadron[hadron] += 1
        prompt_hadrons_by_quark = {}
        for hrow, (origin, qidx) in enumerate(zip(after["QHadron"]["charmOrigin"], after["QHadron"]["quarkIdx"])):
            if origin == 1 and qidx >= 0:
                prompt_hadrons_by_quark.setdefault(qidx, []).append(hrow)
        for qidx, hadron_rows in prompt_hadrons_by_quark.items():
            stats["confirmed prompt-charm systems"] += 1
            n_decay = sum(accepted_decay_by_hadron[hrow] for hrow in hadron_rows)
            if n_decay == 2:
                stats["prompt-charm systems with exactly two accepted D tracks"] += 1
            n_unique = sum(1 for quality, membership, quark, pt, eta in zip(
                charged["associationQuality"], charged["membershipType"], charged["quarkIdx"], charged["pt"], charged["eta"])
                if quality == 1 and membership == 2 and quark == qidx and pt >= 0.4 and abs(eta) <= 2.5)
            n_certain = sum(1 for membership, quark, certain, pt, eta in zip(
                charged["membershipType"], charged["quarkIdx"], charged["fragmentationSourceCertain"], charged["pt"], charged["eta"])
                if membership == 2 and quark == qidx and certain == 1 and pt >= 0.4 and abs(eta) <= 2.5)
            stats["accepted unique additional particles in confirmed prompt-charm systems"] += n_unique
            stats["accepted source-certain additional particles in confirmed prompt-charm systems"] += n_certain
            if n_decay == 2 and n_unique:
                stats["two-D-track prompt systems with unique additional candidates"] += 1
            if n_decay == 2 and n_certain:
                stats["two-D-track prompt systems with source-certain additional particles"] += 1
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("approved", help="approved post-ancestry-boundary ROOT output")
    parser.add_argument("cleaned", help="ROOT output from cleaned producer")
    args = parser.parse_args()
    old_file, old_tree = open_tree(args.approved)
    new_file, new_tree = open_tree(args.cleaned)
    try:
        stats = compare(read_rows(old_tree, True), read_rows(new_tree, False))
    finally:
        old_file.Close()
        new_file.Close()
    print("PASS: event-by-event and row-by-row saved values match")
    for label in ("events", "QFrag rows", "QHadron rows", "QHadron unique ancestry assignments", "QGenCharged rows",
                  "heavy decay charged rows", "unique additional candidates", "confirmed prompt-charm systems",
                  "prompt-charm systems with exactly two accepted D tracks",
                  "accepted unique additional particles in confirmed prompt-charm systems",
                  "accepted source-certain additional particles in confirmed prompt-charm systems",
                  "two-D-track prompt systems with unique additional candidates",
                  "two-D-track prompt systems with source-certain additional particles"):
        print("%s: %d" % (label, stats[label]))
    denominator = stats["prompt-charm systems with exactly two accepted D tracks"]
    numerator = stats["two-D-track prompt systems with unique additional candidates"]
    fraction = float(numerator) / denominator if denominator else 0.0
    print("Exactly-two-D-track prompt systems with >=1 unique accepted extra: %d/%d (%.4f)" %
          (numerator, denominator, fraction))
    print("Float tolerance: rtol=%g atol=%g" % (FLOAT_RTOL, FLOAT_ATOL))
    return 0


if __name__ == "__main__":
    sys.exit(main())
