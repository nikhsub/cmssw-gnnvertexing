#!/usr/bin/env python3
"""Check V2 table invariants and compare Stage-1, V1, and V2 assignments."""
import argparse
import math
import ROOT

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("file", help="NanoAOD ROOT file produced with the lineage V2 plugin")
args = parser.parse_args()

root = ROOT.TFile.Open(args.file)
if not root or root.IsZombie():
    raise SystemExit("Could not open ROOT file: " + args.file)
events = root.Get("Events")
if not events:
    raise SystemExit("ROOT file has no Events tree")

counts = {key: 0 for key in (
    "v2Unique", "v2Ambiguous", "v2Missing", "v2Unknown", "v2KinematicsValid",
    "v2KinematicsInvalid", "unknownOriginUnique", "prompt", "bToD", "originUnknown",
)}
transitions = {}
checked = 0
for event in events:
    n_qfrag = len(event.QFrag_flavour)
    branches = (
        event.QHadron_associationQuality,
        event.QHadron_lineageAssociationQuality,
        event.QHadron_lineageV2AssociationQuality,
        event.QHadron_lineageV2AssociationReason,
        event.QHadron_lineageV2QuarkIdx,
        event.QHadron_lineageV2OriginCategory,
        event.QHadron_lineageV2SourceCategory,
        event.QHadron_lineageV2FragmentationSourceCertain,
        event.QHadron_lineageV2AncestryIncomplete,
        event.QHadron_lineageV2AncestryTruncated,
        event.QHadron_lineageV2PreHadPrunedIdx,
        event.QHadron_lineageV2PreHadSelectionReason,
        event.QHadron_lineageV2PreHadKinematicsValid,
        event.QHadron_lineageV2PreHadKinematicsReason,
        event.QHadron_lineageV2ZHadron,
        event.QHadron_lineageV2DeltaRHadron,
        event.QHadron_lineageV2NCompatibleHFLineages,
        event.QHadron_lineageV2NRoots,
        event.QHadron_lineageV2NTerminals,
        event.QHadron_lineageV2NBranchingNodes,
        event.QHadron_lineageV2NMultiMotherNodes,
        event.QHadron_lineageV2PathUnique,
        event.QHadron_flavour,
    )
    if len({len(values) for values in branches}) != 1:
        raise AssertionError("QHadron V2 branch lengths differ in event %s" % event.event)
    for (old, v1, v2, reason, qidx, origin, source, source_certain, ancestry_incomplete, ancestry_truncated, prehad,
         prehad_selection, valid, kin_reason, z, dr, n_compatible, n_roots, n_terminals, n_branching,
         n_multi_mother, path_unique, flavour) in zip(*branches):
        checked += 1
        transitions[(old, v1, v2)] = transitions.get((old, v1, v2), 0) + 1
        counts[{0: "v2Missing", 1: "v2Unique", 2: "v2Ambiguous", 3: "v2Unknown"}.get(v2, "v2Unknown")] += 1
        if qidx >= 0 and qidx >= n_qfrag:
            raise AssertionError("lineageV2QuarkIdx out of QFrag bounds: %d >= %d" % (qidx, n_qfrag))
        if v2 == 1 and qidx < 0 and reason != 7:
            raise AssertionError("unique V2 lineage without QFrag row must report reason 7")
        if v2 == 1 and path_unique != 1:
            raise AssertionError("unique V2 association must have one retained ancestry path")
        if v2 == 1 and (ancestry_incomplete or ancestry_truncated):
            raise AssertionError("unique V2 association cannot use incomplete/truncated ancestry")
        if path_unique == 1 and (n_roots != 1 or n_terminals != 1 or n_branching or n_multi_mother):
            raise AssertionError("unique V2 path has inconsistent copy topology diagnostics")
        if n_compatible > 1 and v2 != 2:
            raise AssertionError("multiple compatible heavy lineages must remain ambiguous")
        if bool(valid) != (kin_reason == 0):
            raise AssertionError("V2 kinematics-valid flag and reason code disagree")
        if prehad_selection in (1, 2) and prehad < 0:
            raise AssertionError("resolved pre-hadron copy selection has no pruned index")
        if prehad_selection in (3, 4) and prehad >= 0:
            raise AssertionError("unresolved pre-hadron copy selection has an index")
        if valid:
            if prehad < 0 or not math.isfinite(z) or not math.isfinite(dr) or z == -999.0 or dr == -999.0:
                raise AssertionError("valid V2 kinematics lack a reference or finite observables")
            counts["v2KinematicsValid"] += 1
        else:
            if z != -999.0 or dr != -999.0:
                raise AssertionError("invalid V2 kinematics must use -999 sentinels")
            counts["v2KinematicsInvalid"] += 1
        if flavour == 4:
            counts[{1: "prompt", 2: "bToD", 3: "originUnknown"}.get(origin, "originUnknown")] += 1
            if origin == 3 and v2 == 1:
                counts["unknownOriginUnique"] += 1
                if source != 1 or source_certain:
                    raise AssertionError("unknown-origin D must retain D ancestry without claiming a certain source")
            if origin == 2 and source not in (2, 3):
                raise AssertionError("B-to-D classification must retain a B or unresolved-B source category")
        if source == 3 and qidx >= 0:
            raise AssertionError("unresolved positive B source must not map to a QFrag row")

if checked == 0:
    print("No QHadron rows found")
else:
    print("V2 association: unique=%d ambiguous=%d missing=%d unknown=%d" %
          (counts["v2Unique"], counts["v2Ambiguous"], counts["v2Missing"], counts["v2Unknown"]))
    print("V2 pre-hadronization kinematics: valid=%d invalid=%d" %
          (counts["v2KinematicsValid"], counts["v2KinematicsInvalid"]))
    print("D origin: prompt=%d B-to-D=%d unknown=%d; unknown-origin with unique ancestry=%d" %
          (counts["prompt"], counts["bToD"], counts["originUnknown"], counts["unknownOriginUnique"]))
    print("(old Stage1 quality, V1 quality, V2 quality) counts:")
    for key, value in sorted(transitions.items()):
        print("  %s: %d" % (key, value))
    print("Validated %d QHadron rows" % checked)
