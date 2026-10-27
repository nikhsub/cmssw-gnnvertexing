#!/usr/bin/env python3
"""Validate optional pair vertex-fit columns and compare pair identity to baseline."""
import argparse
import csv
import math
import sys
import re

try:
    import ROOT
except ImportError as exc:
    raise SystemExit("Run this script inside a CMSSW environment with PyROOT available") from exc


PREFIX = "trackPair_"  # NanoAOD writes the FlatTable name, not its producer label, in branch names.
FIT_FLOATS = [
    "vtxChi2", "vtxNdof", "vtxChi2Ndof", "vtxProb", "vtxX", "vtxY", "vtxZ",
    "vtxCovXX", "vtxCovXY", "vtxCovXZ", "vtxCovYY", "vtxCovYZ", "vtxCovZZ",
    "vtxSigmaX", "vtxSigmaY", "vtxSigmaZ", "vtxLxy", "vtxLxyErr", "vtxLxySig",
    "vtxSignedLxy", "vtxSignedLxySig", "vtxL3D", "vtxL3DErr", "vtxL3DSig",
    "vtxCosPointing2D", "vtxCosPointing3D", "vtxPerpDisplacement", "vtxPerpDisplacementSig",
]
VALID_FLAGS = ["vtxFitValid", "vtxCovarianceValid", "vtxPVValid", "vtxLxyValid", "vtxSignedLxyValid", "vtxL3DValid",
               "vtxPointing2DValid", "vtxPointing3DValid", "vtxPerpendicularValid"]
FIT_BRANCHES = ["vtxFitStatus"] + FIT_FLOATS + VALID_FLAGS
SENTINEL = -9999.0


def branch_values(tree, name):
    value = getattr(tree, name)
    try:
        return list(value)
    except TypeError:
        return [value]


def find_branch(tree, suffix):
    exact = PREFIX + suffix
    if tree.GetBranch(exact):
        return exact
    matches = [b.GetName() for b in tree.GetListOfBranches() if b.GetName().endswith("_" + suffix)]
    if len(matches) == 1:
        return matches[0]
    raise RuntimeError("Could not uniquely locate ROOT branch for " + suffix + ": " + repr(matches))


def open_tree(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)
    tree = f.Get("Events")
    if not tree:
        raise RuntimeError(path + " has no Events tree")
    return f, tree


def read_time_report(path):
    if not path:
        return {}
    with open(path) as handle:
        content = handle.read()
    def number(pattern):
        match = re.search(pattern, content, re.MULTILINE)
        return float(match.group(1)) if match else -1.0
    elapsed = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", content)
    wall_seconds = -1.0
    if elapsed:
        parts = [float(x) for x in elapsed.group(1).split(":")]
        wall_seconds = sum(value * (60 ** index) for index, value in enumerate(reversed(parts)))
    user = number(r"^\s*User time \(seconds\):\s*([0-9.]+)")
    system = number(r"^\s*System time \(seconds\):\s*([0-9.]+)")
    rss = number(r"^\s*Maximum resident set size \(kbytes\):\s*([0-9.]+)")
    return {"userSeconds": user, "systemSeconds": system, "cpuSeconds": user + system if user >= 0 and system >= 0 else -1.0,
            "wallSeconds": wall_seconds, "peakRssKB": rss}


def collect(path, require_fit):
    f, tree = open_tree(path)
    names = {b.GetName() for b in tree.GetListOfBranches()}
    trk1_name = find_branch(tree, "trk1Idx")
    trk2_name = find_branch(tree, "trk2Idx")
    pair_prefix = trk1_name[:-len("trk1Idx")]
    fit_names = {}
    if require_fit:
        missing = [name for name in FIT_BRANCHES if not any(branch.endswith("_" + name) for branch in names)]
        if missing:
            raise RuntimeError(path + " missing fit branches: " + ", ".join(missing))
        fit_names = {name: find_branch(tree, name) for name in FIT_BRANCHES}
    summary_names = {}
    if require_fit:
        for name in ("nTracks", "nSelectedPairs", "nFitAttempts", "nFitSuccess", "nCovarianceInvalid",
                     "nInvalidTrackState", "nFitFailure",
                     "transientTrackBuildMs", "vertexFitMs", "producerTotalMs"):
            matches = [b.GetName() for b in tree.GetListOfBranches() if b.GetName().endswith("vertexFitSummary_" + name)]
            if len(matches) != 1:
                raise RuntimeError(path + " missing or ambiguous vertexFitSummary_" + name + ": " + repr(matches))
            summary_names[name] = matches[0]

    pair_ids, original_pair_columns, checks, covariance_examples = [], {}, {}, []
    def bad_check(name):
        checks[name] = checks.get(name, 0) + 1
    stats = {"events": int(tree.GetEntries()), "pairs": 0, "attempts": 0, "success": 0,
                           "failed": 0, "invalid_tracks": 0, "covariance_invalid": 0,
                           "bad_rows": 0, "timing_fit_ms": 0.0,
                           "timing_build_ms": 0.0, "timing_total_ms": 0.0,
                           "minimumNormalizedCovarianceDiagonal": math.inf,
                           "minimumNormalized2x2PrincipalMinor": math.inf,
                           "minimumNormalizedCovarianceDeterminant": math.inf}
    for entry in tree:
        first, second = branch_values(tree, trk1_name), branch_values(tree, trk2_name)
        if len(first) != len(second):
            raise RuntimeError(path + ": endpoint branch lengths differ")
        pairs = list(zip(map(int, first), map(int, second)))
        pair_ids.append(pairs)
        for branch in names:
            if branch.startswith(pair_prefix):
                column = branch[len(pair_prefix):]
                if column not in FIT_BRANCHES:
                    original_pair_columns.setdefault(column, []).append(branch_values(tree, branch))
        stats["pairs"] += len(pairs)
        if not require_fit:
            continue

        columns = {key: branch_values(tree, branch) for key, branch in fit_names.items()}
        for key, values in columns.items():
            if len(values) != len(pairs):
                raise RuntimeError(path + ": " + key + " has " + str(len(values)) + " rows for " + str(len(pairs)) + " pairs")
        valid = list(map(int, columns["vtxFitValid"]))
        status = list(map(int, columns["vtxFitStatus"]))
        stats["attempts"] += len(pairs)
        stats["success"] += sum(valid)
        covariance_valid = list(map(int, columns["vtxCovarianceValid"]))
        event_covariance_invalid = sum(ok and not cov_ok for ok, cov_ok in zip(valid, covariance_valid))
        stats["covariance_invalid"] += event_covariance_invalid
        stats["failed"] += sum(s not in (1, 2) for s in status)
        stats["invalid_tracks"] += sum(s == 2 for s in status)
        for i, ok in enumerate(valid):
            if status[i] not in (1, 2, 3, 4) or ok not in (0, 1) or ok != int(status[i] == 1):
                bad_check("status_validity_mismatch")
            for flag in VALID_FLAGS[1:]:
                flag_value = int(columns[flag][i])
                if flag_value not in (0, 1):
                    bad_check("non_boolean_validity_flag")
            if not covariance_valid[i] and any(int(columns[flag][i]) for flag in
                                               ("vtxLxyValid", "vtxSignedLxyValid", "vtxL3DValid", "vtxPerpendicularValid")):
                bad_check("covariance_dependent_feature_marked_valid_with_invalid_covariance")
            for flag, value_names in {
                "vtxLxyValid": ("vtxLxy", "vtxLxyErr", "vtxLxySig"),
                "vtxSignedLxyValid": ("vtxSignedLxy", "vtxSignedLxySig"),
                "vtxL3DValid": ("vtxL3D", "vtxL3DErr", "vtxL3DSig"),
                "vtxPointing2DValid": ("vtxCosPointing2D",),
                "vtxPointing3DValid": ("vtxCosPointing3D",),
                "vtxPerpendicularValid": ("vtxPerpDisplacement", "vtxPerpDisplacementSig"),
            }.items():
                expected_sentinel = int(columns[flag][i]) == 0
                if any((float(columns[name][i]) == SENTINEL) != expected_sentinel for name in value_names):
                    bad_check("derived_validity_sentinel_mismatch")
            for name in FIT_FLOATS:
                x = float(columns[name][i])
                if ok and not math.isfinite(x):
                    bad_check("nonfinite_valid_fit_value")
                if not ok and x != SENTINEL:
                    bad_check("invalid_fit_value_not_sentinel")
            probability = float(columns["vtxProb"][i])
            if ok and not 0.0 <= probability <= 1.0:
                bad_check("probability_out_of_range")
            if ok:
                xx, xy, xz = (float(columns[n][i]) for n in ("vtxCovXX", "vtxCovXY", "vtxCovXZ"))
                yy, yz, zz = (float(columns[n][i]) for n in ("vtxCovYY", "vtxCovYZ", "vtxCovZZ"))
                scale = max(abs(xx), abs(yy), abs(zz), 1e-30)
                tol1, tol2, tol3 = 1e-10 * scale, 1e-10 * scale**2, 1e-10 * scale**3
                minor2 = (xx * yy - xy * xy, xx * zz - xz * xz, yy * zz - yz * yz)
                determinant = xx * yy * zz + 2 * xy * xz * yz - xx * yz**2 - yy * xz**2 - zz * xy**2
                bad_diagonal = min(xx, yy, zz) < -tol1
                bad_minor2 = min(minor2) < -tol2
                bad_determinant = determinant < -tol3
                stats["minimumNormalizedCovarianceDiagonal"] = min(
                    stats["minimumNormalizedCovarianceDiagonal"], xx / scale, yy / scale, zz / scale)
                stats["minimumNormalized2x2PrincipalMinor"] = min(
                    stats["minimumNormalized2x2PrincipalMinor"], *(value / scale**2 for value in minor2))
                stats["minimumNormalizedCovarianceDeterminant"] = min(
                    stats["minimumNormalizedCovarianceDeterminant"], determinant / scale**3)
                if covariance_valid[i] and bad_diagonal:
                    bad_check("covariance_negative_diagonal")
                if covariance_valid[i] and bad_minor2:
                    bad_check("covariance_negative_2x2_principal_minor")
                if covariance_valid[i] and bad_determinant:
                    bad_check("covariance_negative_determinant")
                if (bad_diagonal or bad_minor2 or bad_determinant) and len(covariance_examples) < 5:
                    covariance_examples.append({
                        "event": len(pair_ids) - 1, "pairRow": i,
                        "trk1Idx": pairs[i][0], "trk2Idx": pairs[i][1],
                        "covXX": xx, "covXY": xy, "covXZ": xz,
                        "covYY": yy, "covYZ": yz, "covZZ": zz,
                        "normalized2x2Min": min(minor2) / scale**2,
                        "normalizedDeterminant": determinant / scale**3,
                    })
        sums = {name: branch_values(tree, branch) for name, branch in summary_names.items()}
        if any(len(values) != 1 for values in sums.values()):
            raise RuntimeError(path + ": expected singleton event summary rows")
        if int(sums["nSelectedPairs"][0]) != len(pairs) or int(sums["nFitAttempts"][0]) != len(pairs):
            bad_check("summary_pair_count_mismatch")
        if int(sums["nFitSuccess"][0]) != sum(valid):
            bad_check("summary_success_count_mismatch")
        if int(sums["nCovarianceInvalid"][0]) != event_covariance_invalid:
            bad_check("summary_covariance_invalid_count_mismatch")
        stats["timing_build_ms"] += max(0.0, float(sums["transientTrackBuildMs"][0]))
        stats["timing_fit_ms"] += max(0.0, float(sums["vertexFitMs"][0]))
        stats["timing_total_ms"] += max(0.0, float(sums["producerTotalMs"][0]))
    f.Close()
    for key in ("minimumNormalizedCovarianceDiagonal", "minimumNormalized2x2PrincipalMinor",
                "minimumNormalizedCovarianceDeterminant"):
        if math.isinf(stats[key]):
            stats[key] = -1.0
    stats["bad_rows"] = sum(checks.values())
    for name, count in checks.items():
        stats["check_" + name] = count
    stats["covariance_examples"] = covariance_examples
    return pair_ids, original_pair_columns, stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, help="fit-disabled ROOT output")
    parser.add_argument("--enabled", required=True, help="fit-enabled ROOT output")
    parser.add_argument("--baseline-time", help="GNU /usr/bin/time -v report for baseline run")
    parser.add_argument("--enabled-time", help="GNU /usr/bin/time -v report for enabled run")
    parser.add_argument("--csv", default="two_track_vertex_fit_validation.csv")
    args = parser.parse_args()
    base_ids, base_pair_columns, base = collect(args.baseline, False)
    enabled_ids, enabled_pair_columns, enabled = collect(args.enabled, True)
    base.update(read_time_report(args.baseline_time))
    enabled.update(read_time_report(args.enabled_time))
    enabled["fitMillisecondsPerAttempt"] = enabled["timing_fit_ms"] / enabled["attempts"] if enabled["attempts"] else 0.0
    if base.get("cpuSeconds", -1) >= 0 and enabled.get("cpuSeconds", -1) >= 0 and enabled["attempts"]:
        enabled["incrementalCpuSecondsPerAttempt"] = (enabled["cpuSeconds"] - base["cpuSeconds"]) / enabled["attempts"]
    else:
        enabled["incrementalCpuSecondsPerAttempt"] = -1.0
    identity_ok = len(base_ids) == len(enabled_ids) and all(a == b for a, b in zip(base_ids, enabled_ids))
    errors = []
    if not identity_ok:
        errors.append("Pair row counts or (trk1Idx,trk2Idx) identities differ between baseline and enabled files")
    if enabled["bad_rows"]:
        errors.append("Detected " + str(enabled["bad_rows"]) + " invalid fit/table consistency checks")
    if base["events"] != enabled["events"]:
        errors.append("Event counts differ between baseline and enabled outputs")
    common = set(base_pair_columns) & set(enabled_pair_columns)
    if not common:
        errors.append("No common original pair columns found for baseline comparison")
        common = set()
    mismatches = []
    for branch in common:
        if base_pair_columns[branch] != enabled_pair_columns[branch]:
            mismatches.append(branch)
    enabled["original_pair_column_mismatches"] = len(mismatches)
    rows = []
    for label, values in (("baseline", base), ("enabled", enabled)):
        row = {"sample": label}
        row.update(values)
        rows.append(row)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with open(args.csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if mismatches:
        errors.extend("Original pair feature changed between runs: " + branch for branch in mismatches)
    for label, values in (("baseline", base), ("enabled", enabled)):
        for key, value in sorted(values.items()):
            if key.startswith("check_") and value:
                print(label + " " + key + ":", value)
    for example in enabled.get("covariance_examples", []):
        print("Representative non-PSD covariance:", example)
    if errors:
        print("Validation CSV written despite failures:", args.csv)
        raise SystemExit("; ".join(errors))
    print("Pair identities match for", base["events"], "events and", base["pairs"], "pairs")
    print("Enabled fit success:", enabled["success"], "/", enabled["attempts"])
    print("Enabled cumulative fit time (ms):", round(enabled["timing_fit_ms"], 3))
    if enabled.get("incrementalCpuSecondsPerAttempt", -1.) >= 0.:
        print("Incremental CPU seconds per attempted fit:", enabled["incrementalCpuSecondsPerAttempt"])
    print("Validation CSV:", args.csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
