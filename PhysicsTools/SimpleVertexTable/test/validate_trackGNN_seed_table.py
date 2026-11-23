#!/usr/bin/env python3
"""Independently reconstruct and validate top-K TrackGNN Seed tables.

Requires uproot, awkward, and NumPy. Input is a NanoAOD-style Events tree from
TrackGNNSeedProducer. An optional --baseline ROOT file can supply the previous
all-pair trackPair_vtx* output for numerical fit comparison.
"""
import argparse
import math
import sys

try:
    import awkward as ak
    import numpy as np
    import uproot
except ImportError as exc:
    raise SystemExit("Install/import uproot, awkward, and NumPy in the validation environment: %s" % exc)


def values(arrays, name, event):
    return ak.to_list(arrays[name][event])


def sigmoid(value):
    if value >= 0:
        return 1.0 / (1.0 + math.exp(-value))
    exp_value = math.exp(value)
    return exp_value / (1.0 + exp_value)


def require(condition, errors, message):
    if not condition:
        errors.append(message)


def covariance_psd(row, prefix="Seed_"):
    xx, xy, xz = (float(row[prefix + x]) for x in ("vtxCovXX", "vtxCovXY", "vtxCovXZ"))
    yy, yz, zz = (float(row[prefix + x]) for x in ("vtxCovYY", "vtxCovYZ", "vtxCovZZ"))
    matrix = np.asarray([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=np.float64)
    if not np.isfinite(matrix).all():
        return False
    scale = max(float(np.max(np.abs(matrix))), 1.e-30)
    return bool(np.linalg.eigvalsh(matrix).min() >= -1.e-8 * scale)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--max-seeds", type=int, default=128)
    parser.add_argument("--baseline", help="Optional all-pair-fit ROOT file for selected-pair fit comparison")
    args = parser.parse_args()
    errors, ties = [], 0
    tie_disagreements = []

    seed_branches = [
        "nSeed", "Seed_rank", "Seed_trk1Idx", "Seed_trk2Idx", "Seed_edgeIdx12", "Seed_edgeIdx21",
        "Seed_logit12", "Seed_logit21", "Seed_score12", "Seed_score21", "Seed_meanLogit", "Seed_score",
        "Seed_pt", "Seed_eta", "Seed_phi", "Seed_mass", "Seed_deltaR", "Seed_dca", "Seed_dcaSig", "Seed_cptopv",
        "Seed_vtxFitValid", "Seed_vtxFitStatus", "Seed_vtxCovarianceValid", "Seed_vtxPVValid",
        "Seed_vtxLxyValid", "Seed_vtxSignedLxyValid", "Seed_vtxL3DValid", "Seed_vtxPointing2DValid",
        "Seed_vtxPointing3DValid", "Seed_vtxPerpendicularValid", "Seed_vtxChi2", "Seed_vtxNdof",
        "Seed_vtxChi2Ndof", "Seed_vtxProb", "Seed_vtxX", "Seed_vtxY", "Seed_vtxZ", "Seed_vtxCovXX",
        "Seed_vtxCovXY", "Seed_vtxCovXZ", "Seed_vtxCovYY", "Seed_vtxCovYZ", "Seed_vtxCovZZ",
        "Seed_vtxSigmaX", "Seed_vtxSigmaY", "Seed_vtxSigmaZ", "Seed_vtxLxy", "Seed_vtxLxyErr",
        "Seed_vtxLxySig", "Seed_vtxSignedLxy", "Seed_vtxSignedLxySig", "Seed_vtxL3D", "Seed_vtxL3DErr",
        "Seed_vtxL3DSig", "Seed_vtxCosPointing2D", "Seed_vtxCosPointing3D", "Seed_vtxPerpDisplacement",
        "Seed_vtxPerpDisplacementSig",
    ]
    summary_branches = ["nSeed", "SeedSummary_nTracks", "SeedSummary_nDirectedEdges", "SeedSummary_nUniquePairs",
                        "SeedSummary_nValidScoredPairs", "SeedSummary_nSeeds", "SeedSummary_nFitAttempts",
                        "SeedSummary_nFitSuccess", "SeedSummary_nFitFailure", "SeedSummary_nCovarianceInvalid",
                        "SeedSummary_nMissingReciprocal", "SeedSummary_nInvalidPairIndices",
                        "SeedSummary_nInvalidScores", "SeedSummary_nDuplicateDirectedEdges"]
    required = ["ntrack", "ntrackPair", "trackPair_trk1Idx", "trackPair_trk2Idx", "trackPair_seederLogit",
                "trackPair_seederValid"] + seed_branches + summary_branches[1:]
    with uproot.open(args.root_file) as root:
        if args.tree not in root:
            raise SystemExit("Missing tree %r" % args.tree)
        tree = root[args.tree]
        available = set(tree.keys())
        missing = [name for name in required if name not in available]
        if missing:
            raise SystemExit("Missing required branches: %s" % missing)
        optional_mc_branches = ("track_hadidx", "trackPair_edgeLabel", "GV_isD", "GV_fromHF")
        optional_mc = all(x in available for x in optional_mc_branches)
        read_branches = required + [x for x in ("run", "luminosityBlock", "event") if x in available]
        if optional_mc:
            read_branches += list(optional_mc_branches)
        arrays = tree.arrays(read_branches, library="ak")
        event_count = tree.num_entries
        total_seeds = total_fits = total_success = 0
        compared_pairs = 0
        recovery = {k: {"promptD": 0, "promptD_recovered": 0, "n2": 0, "n2_recovered": 0,
                        "n3plus": 0, "n3plus_recovered": 0, "true_fit": 0, "true_fit_success": 0,
                        "fake_fit": 0, "fake_fit_success": 0, "true_seeds": 0, "fake_seeds": 0}
                    for k in (16, 32, 64, 128)}

        for ev in range(event_count):
            ntracks = int(arrays["ntrack"][ev])
            nedges = int(arrays["ntrackPair"][ev])
            edge_i = [int(x) for x in values(arrays, "trackPair_trk1Idx", ev)]
            edge_j = [int(x) for x in values(arrays, "trackPair_trk2Idx", ev)]
            edge_logit = [float(x) for x in values(arrays, "trackPair_seederLogit", ev)]
            edge_valid = [int(x) for x in values(arrays, "trackPair_seederValid", ev)]
            by_pair, first_order = {}, []
            invalid_index_rows = invalid_score_rows = duplicate_rows = 0
            for row, (a, b, logit, valid) in enumerate(zip(edge_i, edge_j, edge_logit, edge_valid)):
                if a < 0 or b < 0 or a >= ntracks or b >= ntracks or a == b:
                    invalid_index_rows += 1
                    continue
                row_is_valid = valid == 1 and math.isfinite(logit)
                if not row_is_valid:
                    invalid_score_rows += 1
                i, j = sorted((a, b))
                key = (i, j)
                if key not in by_pair:
                    by_pair[key] = {"12": None, "21": None}
                    first_order.append(key)
                direction = "12" if (a, b) == (i, j) else "21"
                if by_pair[key][direction] is not None:
                    duplicate_rows += 1
                    continue
                if row_is_valid:
                    by_pair[key][direction] = (row, logit)
                else:
                    by_pair[key][direction] = (row, None)
            missing_reciprocal = sum(1 for k in first_order if by_pair[k]["12"] is None or by_pair[k]["21"] is None)
            ranked = []
            for order, key in enumerate(first_order):
                one, two = by_pair[key]["12"], by_pair[key]["21"]
                if one is None or two is None or one[1] is None or two[1] is None:
                    continue
                # Match the producer's float arithmetic: 0.5f * (logit12 + logit21).
                mean = float(np.float32(0.5) * np.float32(np.float32(one[1]) + np.float32(two[1])))
                ranked.append({"key": key, "edge12": one[0], "edge21": two[0], "logit12": one[1],
                               "logit21": two[1], "mean": mean, "score": sigmoid(mean), "order": order})
            # Python's stable sort preserves first physical-pair occurrence for ties.
            ranked.sort(key=lambda item: -item["mean"])
            selected = ranked[:args.max_seeds]
            ties += sum(1 for a, b in zip(ranked, ranked[1:]) if a["mean"] == b["mean"])
            expected_by_pair = {item["key"]: item for item in ranked}

            count = int(arrays["nSeed"][ev])
            total_seeds += count
            total_fits += int(arrays["SeedSummary_nFitAttempts"][ev])
            total_success += int(arrays["SeedSummary_nFitSuccess"][ev])
            require(count <= args.max_seeds, errors, "event %d: nSeed exceeds maxSeeds" % ev)
            require(count == len(selected), errors, "event %d: saved seed count %d != expected %d" % (ev, count, len(selected)))
            require(int(arrays["SeedSummary_nTracks"][ev]) == ntracks, errors, "event %d: summary nTracks mismatch" % ev)
            require(int(arrays["SeedSummary_nDirectedEdges"][ev]) == nedges, errors, "event %d: summary directed-edge count mismatch" % ev)
            require(int(arrays["SeedSummary_nUniquePairs"][ev]) == len(first_order), errors, "event %d: unique-pair count mismatch" % ev)
            require(int(arrays["SeedSummary_nValidScoredPairs"][ev]) == len(ranked), errors, "event %d: valid-scored-pair count mismatch" % ev)
            require(int(arrays["SeedSummary_nMissingReciprocal"][ev]) == missing_reciprocal, errors, "event %d: missing-reciprocal count mismatch" % ev)
            require(int(arrays["SeedSummary_nInvalidPairIndices"][ev]) == invalid_index_rows, errors, "event %d: invalid-index count mismatch" % ev)
            require(int(arrays["SeedSummary_nInvalidScores"][ev]) == invalid_score_rows, errors, "event %d: invalid-score count mismatch" % ev)
            require(int(arrays["SeedSummary_nDuplicateDirectedEdges"][ev]) == duplicate_rows, errors, "event %d: duplicate-edge count mismatch" % ev)
            require(int(arrays["SeedSummary_nFitAttempts"][ev]) == count, errors, "event %d: fit attempts != nSeed" % ev)
            require(int(arrays["SeedSummary_nFitSuccess"][ev]) + int(arrays["SeedSummary_nFitFailure"][ev]) == count,
                    errors, "event %d: fit success + failure != nSeed" % ev)

            names = [x for x in seed_branches if x != "nSeed"]
            rows = {name: values(arrays, name, ev) for name in names}
            for name, branch_rows in rows.items():
                require(len(branch_rows) == count, errors, "event %d: %s has %d rows, nSeed=%d" %
                        (ev, name, len(branch_rows), count))
            seen_pairs = set()
            n_compare = min([count, len(selected)] + [len(branch_rows) for branch_rows in rows.values()])
            for row in range(n_compare):
                expected = selected[row]
                actual_pair = (int(rows["Seed_trk1Idx"][row]), int(rows["Seed_trk2Idx"][row]))
                require(int(rows["Seed_rank"][row]) == row + 1, errors, "event %d row %d: rank not sequential" % (ev, row))
                if actual_pair != expected["key"]:
                    actual_expected = expected_by_pair.get(actual_pair)
                    if actual_expected is not None and actual_expected["mean"] == expected["mean"]:
                        tie_disagreements.append("event %d rank %d: found %s, expected %s in exact meanLogit tie %.9g" %
                                                 (ev, row + 1, actual_pair, expected["key"], expected["mean"]))
                    errors.append("event %d row %d: pair %s != ranked pair %s" % (ev, row, actual_pair, expected["key"]))
                require(actual_pair[0] < actual_pair[1], errors, "event %d row %d: pair is not canonical" % (ev, row))
                require(0 <= actual_pair[0] < ntracks and 0 <= actual_pair[1] < ntracks, errors,
                        "event %d row %d: track index outside track collection" % (ev, row))
                require(actual_pair not in seen_pairs, errors, "event %d row %d: duplicate Seed physical pair" % (ev, row))
                seen_pairs.add(actual_pair)
                require(int(rows["Seed_edgeIdx12"][row]) == expected["edge12"] and int(rows["Seed_edgeIdx21"][row]) == expected["edge21"],
                        errors, "event %d row %d: directed edge row references disagree" % (ev, row))
                for branch, exp in (("Seed_logit12", expected["logit12"]), ("Seed_logit21", expected["logit21"]),
                                    ("Seed_meanLogit", expected["mean"]), ("Seed_score", expected["score"])):
                    require(math.isfinite(float(rows[branch][row])) and abs(float(rows[branch][row]) - exp) <= 1.e-6,
                            errors, "event %d row %d: %s mismatch" % (ev, row, branch))
                require(abs(float(rows["Seed_score12"][row]) - sigmoid(expected["logit12"])) <= 1.e-6 and
                        abs(float(rows["Seed_score21"][row]) - sigmoid(expected["logit21"])) <= 1.e-6,
                        errors, "event %d row %d: directed score mismatch" % (ev, row))
                if row and float(rows["Seed_meanLogit"][row]) > float(rows["Seed_meanLogit"][row - 1]):
                    errors.append("event %d: seeds are not in descending meanLogit order" % ev)
                fit_valid = int(rows["Seed_vtxFitValid"][row])
                if not int(rows["Seed_vtxFitStatus"][row]):
                    errors.append("event %d row %d: missing fit status" % (ev, row))
                if fit_valid:
                    finite_fit = ("Seed_vtxX", "Seed_vtxY", "Seed_vtxZ", "Seed_vtxChi2", "Seed_vtxNdof",
                                  "Seed_vtxChi2Ndof", "Seed_vtxProb")
                    if not all(math.isfinite(float(rows[x][row])) for x in finite_fit):
                        errors.append("event %d row %d: valid fit has nonfinite quality/position" % (ev, row))
                    if not 0.0 <= float(rows["Seed_vtxProb"][row]) <= 1.0:
                        errors.append("event %d row %d: valid fit probability outside [0,1]" % (ev, row))
                    if int(rows["Seed_vtxCovarianceValid"][row]):
                        covrow = {key: rows[key][row] for key in rows}
                        if not covariance_psd(covrow):
                            errors.append("event %d row %d: covarianceValid matrix is not PSD" % (ev, row))
                    require(not int(rows["Seed_vtxCovarianceValid"][row]) or fit_valid, errors,
                            "event %d row %d: covariance marked valid for invalid fit" % (ev, row))
                else:
                    fit_float_branches = [x for x in rows if x.startswith("Seed_vtx") and
                                          not x.endswith("Valid") and x != "Seed_vtxFitStatus"]
                    for branch in fit_float_branches:
                        if float(rows[branch][row]) != -9999.0:
                            errors.append("event %d row %d: failed fit %s lacks -9999 sentinel" % (ev, row, branch))
                for flag, observables in (("Seed_vtxLxyValid", ("Seed_vtxLxy", "Seed_vtxLxyErr", "Seed_vtxLxySig")),
                                          ("Seed_vtxSignedLxyValid", ("Seed_vtxSignedLxy", "Seed_vtxSignedLxySig")),
                                          ("Seed_vtxL3DValid", ("Seed_vtxL3D", "Seed_vtxL3DErr", "Seed_vtxL3DSig")),
                                          ("Seed_vtxPointing2DValid", ("Seed_vtxCosPointing2D",)),
                                          ("Seed_vtxPointing3DValid", ("Seed_vtxCosPointing3D",)),
                                          ("Seed_vtxPerpendicularValid", ("Seed_vtxPerpDisplacement", "Seed_vtxPerpDisplacementSig"))):
                    valid = int(rows[flag][row])
                    for observable in observables:
                        value = float(rows[observable][row])
                        if valid and not math.isfinite(value):
                            errors.append("event %d row %d: %s valid mask with nonfinite observable" % (ev, row, observable))
                        if not valid and value != -9999.0:
                            errors.append("event %d row %d: %s invalid value lacks fitter sentinel" % (ev, row, observable))
                for flag in ("Seed_vtxLxyValid", "Seed_vtxSignedLxyValid", "Seed_vtxL3DValid",
                             "Seed_vtxPointing2DValid", "Seed_vtxPointing3DValid", "Seed_vtxPerpendicularValid"):
                    if int(rows[flag][row]) and not int(rows["Seed_vtxPVValid"][row]):
                        errors.append("event %d row %d: %s set without a valid PV" % (ev, row, flag))
                for flag in ("Seed_vtxLxyValid", "Seed_vtxSignedLxyValid", "Seed_vtxL3DValid", "Seed_vtxPerpendicularValid"):
                    if int(rows[flag][row]) and not int(rows["Seed_vtxCovarianceValid"][row]):
                        errors.append("event %d row %d: %s set with invalid covariance" % (ev, row, flag))
                compared_pairs += 1

            if optional_mc:
                hadidx = [int(x) for x in values(arrays, "track_hadidx", ev)]
                edge_labels = [float(x) for x in values(arrays, "trackPair_edgeLabel", ev)]
                if len(hadidx) != ntracks or len(edge_labels) != nedges:
                    errors.append("event %d: optional MC track/edge truth rows are not aligned" % ev)
                    continue
                is_d = [int(x) for x in values(arrays, "GV_isD", ev)]
                from_hf = [int(x) for x in values(arrays, "GV_fromHF", ev)]
                prompt_ds = [g for g in range(min(len(is_d), len(from_hf))) if is_d[g] == 1 and from_hf[g] == 0]
                multiplicity = {g: sum(1 for h in hadidx if h == g) for g in prompt_ds}
                n2 = [g for g in prompt_ds if multiplicity[g] == 2]
                n3 = [g for g in prompt_ds if multiplicity[g] >= 3]
                for k in recovery:
                    top = selected[:k]
                    recovered = {g for g in prompt_ds if any(a < len(hadidx) and b < len(hadidx) and hadidx[a] == g and hadidx[b] == g for a, b in (x["key"] for x in top))}
                    recovery[k]["promptD"] += len(prompt_ds); recovery[k]["promptD_recovered"] += len(recovered)
                    recovery[k]["n2"] += len(n2); recovery[k]["n2_recovered"] += sum(g in recovered for g in n2)
                    recovery[k]["n3plus"] += len(n3); recovery[k]["n3plus_recovered"] += sum(g in recovered for g in n3)
                    labels_by_pair = {}
                    for idx, (a, b) in enumerate(zip(values(arrays, "trackPair_trk1Idx", ev), values(arrays, "trackPair_trk2Idx", ev))):
                        labels_by_pair[(int(a), int(b))] = edge_labels[idx]
                    for seed_row, item in enumerate(top):
                        label = labels_by_pair.get(item["key"], labels_by_pair.get((item["key"][1], item["key"][0]), 0.0))
                        fit_rows = rows.get("Seed_vtxFitValid", [])
                        fit_ok = int(fit_rows[seed_row]) if seed_row < len(fit_rows) else 0
                        key_name = "true" if label > 0.5 else "fake"
                        recovery[k][key_name + "_seeds"] += 1
                        recovery[k][key_name + "_fit"] += 1
                        recovery[k][key_name + "_fit_success"] += int(fit_ok)

    if args.baseline:
        with uproot.open(args.baseline) as base:
            bt = base[args.tree]
            bnames = set(bt.keys())
            fit_fields = ("vtxFitValid", "vtxFitStatus", "vtxCovarianceValid", "vtxPVValid", "vtxLxyValid",
                          "vtxSignedLxyValid", "vtxL3DValid", "vtxPointing2DValid", "vtxPointing3DValid",
                          "vtxPerpendicularValid", "vtxChi2", "vtxNdof", "vtxChi2Ndof", "vtxProb", "vtxX",
                          "vtxY", "vtxZ", "vtxCovXX", "vtxCovXY", "vtxCovXZ", "vtxCovYY", "vtxCovYZ",
                          "vtxCovZZ", "vtxSigmaX", "vtxSigmaY", "vtxSigmaZ", "vtxLxy", "vtxLxyErr",
                          "vtxLxySig", "vtxSignedLxy", "vtxSignedLxySig", "vtxL3D", "vtxL3DErr",
                          "vtxL3DSig", "vtxCosPointing2D", "vtxCosPointing3D", "vtxPerpDisplacement",
                          "vtxPerpDisplacementSig")
            needed = {"trackPair_trk1Idx", "trackPair_trk2Idx"} | {"trackPair_" + x for x in fit_fields}
            if not needed.issubset(bnames):
                errors.append("baseline file lacks required old trackPair fit columns: %s" % sorted(needed - bnames))
            else:
                event_id_branches = ("run", "luminosityBlock", "event")
                if not all(name in available and name in bnames for name in event_id_branches):
                    errors.append("baseline comparison requires run/luminosityBlock/event in both files; comparison skipped")
                    b = None
                else:
                    b = bt.arrays(list(needed) + list(event_id_branches), library="ak")
                if b is None:
                    pass
                elif bt.num_entries != event_count:
                    errors.append("baseline event count differs; baseline comparison skipped")
                else:
                    baseline_event_index = {}
                    for bev in range(bt.num_entries):
                        identity = tuple(int(values(b, name, bev)) for name in event_id_branches)
                        if identity in baseline_event_index:
                            errors.append("baseline contains duplicate event identity %s; comparison skipped" % (identity,))
                            baseline_event_index = {}
                            break
                        baseline_event_index[identity] = bev
                    for ev in range(event_count):
                        identity = tuple(int(values(arrays, name, ev)) for name in event_id_branches)
                        if identity not in baseline_event_index:
                            errors.append("baseline lacks event identity %s; selected fit comparison skipped for event row %d" % (identity, ev))
                            continue
                        bev = baseline_event_index[identity]
                        seed_n = int(arrays["nSeed"][ev])
                        bi = [int(x) for x in values(b, "trackPair_trk1Idx", bev)]
                        bj = [int(x) for x in values(b, "trackPair_trk2Idx", bev)]
                        lookup = {}
                        for row, (x, y) in enumerate(zip(bi, bj)):
                            lookup.setdefault((min(x, y), max(x, y)), row)
                    n_baseline_compare = min(seed_n, len(values(arrays, "Seed_trk1Idx", ev)),
                                             len(values(arrays, "Seed_trk2Idx", ev)))
                    for row in range(n_baseline_compare):
                            key = (int(values(arrays, "Seed_trk1Idx", ev)[row]), int(values(arrays, "Seed_trk2Idx", ev)[row]))
                            if key not in lookup:
                                errors.append("baseline event %s: selected pair %s absent" % (identity, key)); continue
                            old = lookup[key]
                            for branch in fit_fields:
                                new_value = values(arrays, "Seed_" + branch, ev)[row]
                                old_value = values(b, "trackPair_" + branch, bev)[old]
                                if branch.endswith("Valid") or branch in ("vtxFitStatus",):
                                    if int(new_value) != int(old_value):
                                        errors.append("baseline event %s pair %s %s differs" % (identity, key, branch))
                                elif not math.isclose(float(new_value), float(old_value), rel_tol=1.e-5, abs_tol=1.e-5):
                                    errors.append("baseline event %s pair %s %s differs" % (identity, key, branch))
    print("Events:", event_count)
    print("Seeds:", total_seeds, "Fit attempts:", total_fits, "Fit successes:", total_success)
    print("Selected rows checked:", compared_pairs)
    print("Exact reciprocal mean-logit ties among valid physical pairs:", ties)
    if tie_disagreements:
        print("Tie-related pair/rank disagreements:", len(tie_disagreements))
        for message in tie_disagreements[:50]:
            print(" -", message)
    elif ties:
        print("No event-by-event rank disagreements occurred within exact mean-logit ties.")
    if optional_mc:
        print("Optional MC seed recovery (same-hadron labels; GV_fromHF is not exact D-from-B ancestry):")
        for k, result in recovery.items():
            print(" K=%d promptD=%d/%d N2=%d/%d N3plus=%d/%d trueSeeds/event=%.3f fakeSeeds/event=%.3f trueFit=%d/%d fakeFit=%d/%d" %
                  (k, result["promptD_recovered"], result["promptD"], result["n2_recovered"], result["n2"],
                   result["n3plus_recovered"], result["n3plus"],
                   result["true_seeds"] / event_count if event_count else 0.0,
                   result["fake_seeds"] / event_count if event_count else 0.0,
                   result["true_fit_success"], result["true_fit"],
                   result["fake_fit_success"], result["fake_fit"]))
    if errors:
        print("Validation failed with %d issue(s):" % len(errors))
        for message in errors[:100]:
            print(" -", message)
        if len(errors) > 100:
            print(" - ... %d more" % (len(errors) - 100))
        return 1
    print("Seed table validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
