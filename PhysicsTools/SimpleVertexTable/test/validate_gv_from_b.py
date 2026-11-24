#!/usr/bin/env python3
"""Validate the GV_fromB truth column and optionally compare legacy truth tables.

Requires uproot, awkward, and NumPy. The optional baseline is matched by
run/luminosityBlock/event; GV_fromB is excluded from legacy-branch comparison.
"""
import argparse
import sys

try:
    import awkward as ak
    import uproot
except ImportError as exc:
    raise SystemExit("Install uproot and awkward to validate ROOT output: %s" % exc)


EVENT_ID = ("run", "luminosityBlock", "event")
B_HADRONS = {521, 511, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242,
             5342, 5512, 5532, 5542, 5554}
D_HADRONS = {411, 421, 431, 4122, 4232, 4132, 4332, 4412, 4422, 4432, 4444}
S_HADRONS = {310, 3122, 3112, 3222, 3212, 3312, 3322, 3334}
TAUS = {15}
RECOGNIZED_HADRONS = B_HADRONS | D_HADRONS | S_HADRONS | TAUS
LEGACY_REQUIRED = (
    "nGV", "GV_pt", "GV_eta", "GV_phi", "GV_x", "GV_y", "GV_z", "GV_x_i", "GV_y_i", "GV_z_i",
    "GV_Hadron_pdgId", "GV_isB", "GV_isTau", "GV_isD", "GV_fromHF", "GV_toHF", "GV_pdgClass",
    "GV_nDaughters", "GV_maxDaughterPairDeltaR",
    "nRejectedGV", "RejectedGV_pt", "RejectedGV_eta", "RejectedGV_phi", "RejectedGV_x", "RejectedGV_y",
    "RejectedGV_z", "RejectedGV_x_i", "RejectedGV_y_i", "RejectedGV_z_i",
    "nGVDaughters", "GVDaughters_pt", "GVDaughters_eta", "GVDaughters_phi", "GVDaughters_vx",
    "GVDaughters_vy", "GVDaughters_vz", "GVDaughters_charge", "GVDaughters_pdgId",
    "GVDaughters_statusFlags", "GVDaughters_hadronIndex",
    "nGVDirectDaughters", "GVDirectDaughters_pt", "GVDirectDaughters_eta", "GVDirectDaughters_phi",
    "GVDirectDaughters_charge", "GVDirectDaughters_pdgId", "GVDirectDaughters_hadronIndex",
    "track_hadidx", "track_label", "trackPair_edgeLabel",
)


def value(arrays, branch, event):
    return ak.to_list(arrays[branch][event])


def scalar(arrays, branch, event):
    result = value(arrays, branch, event)
    if isinstance(result, list):
        if len(result) != 1:
            raise ValueError("Expected scalar %s in event row %d; got %r" % (branch, event, result))
        result = result[0]
    return result


def legacy_branches(branches):
    prefixes = ("GV_", "GVDaughters_", "GVDirectDaughters_")
    explicit = {"nGV", "nGVDaughters", "nGVDirectDaughters", "track_hadidx",
                "track_label", "trackPair_edgeLabel"}
    return sorted(name for name in branches if name != "GV_fromB" and
                  (name in explicit or name.startswith(prefixes)))


def compare_baseline(current_tree, baseline_path, tree_name, errors):
    with uproot.open(baseline_path) as baseline_file:
        if tree_name not in baseline_file:
            errors.append("baseline has no tree %s" % tree_name)
            return
        baseline_tree = baseline_file[tree_name]
        current_names = set(current_tree.keys())
        baseline_names = set(baseline_tree.keys())
        old_branches = legacy_branches(current_names) + ["run", "luminosityBlock", "event"]
        missing = sorted(set(old_branches) - baseline_names)
        if missing:
            errors.append("baseline is missing legacy branch(es): %s" % missing)
            return
        compare = [name for name in old_branches if name not in EVENT_ID]
        current = current_tree.arrays(old_branches, library="ak")
        previous = baseline_tree.arrays(old_branches, library="ak")
        old_index = {}
        for row in range(baseline_tree.num_entries):
            key = tuple(int(scalar(previous, name, row)) for name in EVENT_ID)
            if key in old_index:
                errors.append("baseline has duplicate event identity %s" % (key,))
                return
            old_index[key] = row
        checked = 0
        for row in range(current_tree.num_entries):
            key = tuple(int(scalar(current, name, row)) for name in EVENT_ID)
            if key not in old_index:
                errors.append("baseline has no event %s" % (key,))
                continue
            old_row = old_index[key]
            for branch in compare:
                if value(current, branch, row) != value(previous, branch, old_row):
                    errors.append("legacy branch %s changed in event %s" % (branch, key))
            checked += 1
        if checked != baseline_tree.num_entries:
            errors.append("baseline/current event counts differ (%d vs %d)" %
                          (baseline_tree.num_entries, current_tree.num_entries))
        print("Legacy generator/track/edge branches compared event-by-event:", checked)


def validate(root_path, tree_name="Events", baseline_path=None):
    errors = []
    required = ["run", "luminosityBlock", "event", "GV_fromB"] + list(LEGACY_REQUIRED)
    with uproot.open(root_path) as root_file:
        if tree_name not in root_file:
            raise SystemExit("Missing tree %r in %s" % (tree_name, root_path))
        tree = root_file[tree_name]
        branch_names = set(tree.keys())
        missing = [name for name in required if name not in branch_names]
        if missing:
            raise SystemExit("Missing required branch(es): %s" % missing)
        arrays = tree.arrays(required, library="ak")
        total_rows = 0
        total_from_b = 0
        total_charm = 0
        charm_from_b = 0
        for event in range(tree.num_entries):
            key = tuple(int(scalar(arrays, name, event)) for name in EVENT_ID)
            n_gv = int(scalar(arrays, "nGV", event))
            from_b = [int(x) for x in value(arrays, "GV_fromB", event)]
            from_hf = [int(x) for x in value(arrays, "GV_fromHF", event)]
            is_b = [int(x) for x in value(arrays, "GV_isB", event)]
            is_d = [int(x) for x in value(arrays, "GV_isD", event)]
            is_tau = [int(x) for x in value(arrays, "GV_isTau", event)]
            pdg = [abs(int(x)) for x in value(arrays, "GV_Hadron_pdgId", event)]
            lengths = {len(from_b), len(from_hf), len(is_b), len(is_d), len(is_tau), len(pdg)}
            if lengths != {n_gv}:
                errors.append("event %s: GV column lengths %s disagree with nGV=%d" % (key, sorted(lengths), n_gv))
                continue
            for row, label in enumerate(from_b):
                if label not in (0, 1):
                    errors.append("event %s GV row %d: GV_fromB=%d is not binary" % (key, row, label))
                    continue
                if label == 0:
                    continue
                if pdg[row] not in RECOGNIZED_HADRONS:
                    errors.append("event %s GV row %d: fromB set for unrecognized hadron PDG %d" %
                                  (key, row, pdg[row]))
                if from_hf[row] != 1:
                    errors.append("event %s GV row %d: GV_fromB=1 but GV_fromHF!=1" % (key, row))
                if pdg[row] in B_HADRONS and is_b[row] != 1:
                    errors.append("event %s GV row %d: B-hadron PDG inconsistent with GV_isB" % (key, row))
                if pdg[row] in D_HADRONS and is_d[row] != 1:
                    errors.append("event %s GV row %d: D-hadron PDG inconsistent with GV_isD" % (key, row))
                if pdg[row] in TAUS and is_tau[row] != 1:
                    errors.append("event %s GV row %d: tau PDG inconsistent with GV_isTau" % (key, row))
            total_rows += n_gv
            total_from_b += sum(label == 1 for label in from_b)
            total_charm += sum(flag == 1 for flag in is_d)
            charm_from_b += sum(flag == 1 and label == 1 for flag, label in zip(is_d, from_b))
        print("Events:", tree.num_entries)
        print("Accepted GV rows:", total_rows)
        print("Recognized bottom ancestry (GV_fromB=1):", total_from_b)
        print("Charm GV rows:", total_charm)
        print("Charm GV rows with bottom ancestry:", charm_from_b)
        if baseline_path:
            compare_baseline(tree, baseline_path, tree_name, errors)
    if errors:
        print("Validation failed with %d issue(s):" % len(errors))
        for error in errors[:100]:
            print(" -", error)
        return 1
    print("GV_fromB validation passed")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--baseline", help="Previous ROOT file to compare existing truth branches against")
    args = parser.parse_args()
    return validate(args.root_file, args.tree, args.baseline)


if __name__ == "__main__":
    sys.exit(main())
