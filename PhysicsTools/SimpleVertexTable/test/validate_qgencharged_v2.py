#!/usr/bin/env python3
"""Validate versioned QGenCharged assignments against Stage-1 columns.

Requires uproot and awkward. Rows are matched within each Events entry by
packedGenIdx; no flattened cross-event row alignment is assumed.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

import awkward as ak
import uproot


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_reference(filename, tree_name):
    """Index a pre-update output by event identity, preserving each table's row order."""
    names = ["run", "luminosityBlock", "event", "QFrag_prunedGenIdx", "QHadron_prunedGenIdx",
             "QHadron_lineageV2QuarkIdx", "QHadron_lineageV2QuarkPdgId", "QHadron_lineageV2AssociationQuality",
             "QHadron_lineageV2AssociationReason", "QHadron_lineageV2OriginCategory",
             "QHadron_lineageV2SourceCategory", "QHadron_lineageV2SourceReason",
             "QHadron_lineageV2FragmentationSourceCertain",
             "QGenCharged_packedGenIdx"]
    indexed = {}
    with uproot.open(filename) as root_file:
        tree = root_file[tree_name]
        missing = sorted(set(names) - set(tree.keys()))
        if missing:
            raise SystemExit("reference %s is missing required branches: %s" % (filename, ", ".join(missing)))
        for chunk in tree.iterate(expressions=names, library="ak", step_size="100 MB"):
            columns = {name: ak.to_list(chunk[name]) for name in names}
            for i in range(len(columns["event"])):
                key = (columns["run"][i], columns["luminosityBlock"][i], columns["event"][i])
                if key in indexed:
                    raise AssertionError("duplicate event identity in reference file %s: %r" % (filename, key))
                indexed[key] = {name: columns[name][i] for name in names[3:]}
    return indexed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="NanoAOD ROOT input file(s)")
    parser.add_argument("--reference", nargs="+", help="Optional pre-update ROOT file(s), one per input, for invariant regression checks")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for per-file CSV reports")
    parser.add_argument("--max-examples", type=int, default=100, help="Maximum detailed issue rows per input")
    parser.add_argument("--tree", default="Events", help="ROOT tree path (default: Events)")
    args = parser.parse_args()
    if args.reference and len(args.reference) != len(args.inputs):
        parser.error("--reference must provide exactly one file for each input")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    branches = [
        "run", "luminosityBlock", "event",
        "QFrag_flavour", "QFrag_pdgId",
        "QHadron_flavour", "QHadron_lineageV2QuarkIdx", "QHadron_lineageV2AssociationQuality",
        "QHadron_lineageV2OriginCategory", "QHadron_lineageV2FragmentationSourceCertain",
        "QGenCharged_packedGenIdx", "QGenCharged_pt", "QGenCharged_eta",
        "QGenCharged_quarkIdx", "QGenCharged_hadronIdx", "QGenCharged_isHFDecay",
        "QGenCharged_isFromB", "QGenCharged_isFromPromptD", "QGenCharged_isFromBtoD",
        "QGenCharged_membershipType", "QGenCharged_isAmbiguous",
        "QGenCharged_v2HadronIdx", "QGenCharged_v2NearestHeavyPrunedIdx",
        "QGenCharged_v2HadronAssociationQuality", "QGenCharged_v2HadronAssociationReason",
        "QGenCharged_v2QuarkIdx", "QGenCharged_v2AssociationQuality",
        "QGenCharged_v2AssociationReason", "QGenCharged_v2MembershipType",
        "QGenCharged_v2FragmentationSourceCertain", "QGenCharged_v2IsFromUnknownD",
        "QGenCharged_v2NNearestHeavyHadrons", "QGenCharged_v2NCandidateQuarkLineages",
        "QGenCharged_v2NCandidateQFragSystems",
    ]
    required = set(branches)
    for file_number, filename in enumerate(args.inputs):
        path = Path(filename)
        stem = path.stem
        summary, transitions, examples, acceptance_rows = Counter(), Counter(), [], []
        seen_current_events = set()
        reference = load_reference(args.reference[file_number], args.tree) if args.reference else None
        with uproot.open(path) as root_file:
            tree = root_file[args.tree]
            available = set(tree.keys())
            current_branches = list(branches)
            if reference is not None:
                current_branches += ["QFrag_prunedGenIdx", "QHadron_prunedGenIdx", "QHadron_lineageV2QuarkPdgId",
                                     "QHadron_lineageV2AssociationReason", "QHadron_lineageV2SourceCategory",
                                     "QHadron_lineageV2SourceReason"]
            missing = sorted(set(current_branches) - available)
            if missing:
                raise SystemExit("%s is missing required branches: %s" % (path, ", ".join(missing)))
            for chunk in tree.iterate(expressions=current_branches, library="ak", step_size="100 MB"):
                columns = {name: ak.to_list(chunk[name]) for name in current_branches}
                for iev in range(len(columns["event"])):
                    get = lambda name: columns[name][iev]
                    event_key = (get("run"), get("luminosityBlock"), get("event"))
                    seen_current_events.add(event_key)
                    qflav = get("QFrag_flavour")
                    qpdg = get("QFrag_pdgId")
                    hflav = get("QHadron_flavour")
                    hq = get("QHadron_lineageV2QuarkIdx")
                    hquality = get("QHadron_lineageV2AssociationQuality")
                    horigin = get("QHadron_lineageV2OriginCategory")
                    hcertain = get("QHadron_lineageV2FragmentationSourceCertain")
                    packed = get("QGenCharged_packedGenIdx")
                    arrays = {name: get(name) for name in branches if name.startswith("QGenCharged_")}
                    lengths = {name: len(value) for name, value in arrays.items()}
                    if len(set(lengths.values())) != 1:
                        raise AssertionError("QGenCharged branch lengths differ in event %s: %r" % (get("event"), lengths))
                    if len(set(packed)) != len(packed):
                        raise AssertionError("duplicate packedGenIdx in event %s" % get("event"))
                    summary["events"] += 1
                    summary["stable_charged_particles"] += len(packed)
                    summary["unique_selected_QGenCharged_rows"] += len(packed)
                    if reference is not None:
                        baseline = reference.get(event_key)
                        mismatch = None
                        if baseline is None:
                            mismatch = "event missing from reference"
                        elif get("QFrag_prunedGenIdx") != baseline["QFrag_prunedGenIdx"]:
                            mismatch = "QFrag representative/order changed"
                        elif get("QHadron_prunedGenIdx") != baseline["QHadron_prunedGenIdx"]:
                            mismatch = "QHadron selection/order changed"
                        elif get("QHadron_lineageV2QuarkIdx") != baseline["QHadron_lineageV2QuarkIdx"]:
                            mismatch = "QHadron V2 quark association changed"
                        elif get("QHadron_lineageV2QuarkPdgId") != baseline["QHadron_lineageV2QuarkPdgId"]:
                            mismatch = "QHadron V2 quark flavour changed"
                        elif get("QHadron_lineageV2AssociationQuality") != baseline["QHadron_lineageV2AssociationQuality"]:
                            mismatch = "QHadron V2 association quality changed"
                        elif get("QHadron_lineageV2AssociationReason") != baseline["QHadron_lineageV2AssociationReason"]:
                            mismatch = "QHadron V2 association reason changed"
                        elif get("QHadron_lineageV2OriginCategory") != baseline["QHadron_lineageV2OriginCategory"]:
                            mismatch = "D origin classification changed"
                        elif get("QHadron_lineageV2SourceCategory") != baseline["QHadron_lineageV2SourceCategory"]:
                            mismatch = "QHadron V2 source category changed"
                        elif get("QHadron_lineageV2SourceReason") != baseline["QHadron_lineageV2SourceReason"]:
                            mismatch = "QHadron V2 source reason changed"
                        elif get("QHadron_lineageV2FragmentationSourceCertain") != baseline["QHadron_lineageV2FragmentationSourceCertain"]:
                            mismatch = "QHadron source certainty changed"
                        elif packed != baseline["QGenCharged_packedGenIdx"]:
                            mismatch = "QGenCharged selected particle set/order changed"
                        if mismatch:
                            summary["invariant_failures"] += 1
                            if len(examples) < args.max_examples:
                                examples.append({"file": str(path), "run": event_key[0], "luminosityBlock": event_key[1],
                                                 "event": event_key[2], "packedGenIdx": "", "issue": mismatch})
                        summary["reference_events_compared"] += baseline is not None
                    qidx_old = get("QGenCharged_quarkIdx")
                    old_amb = get("QGenCharged_isAmbiguous")
                    old_hf = get("QGenCharged_isHFDecay")
                    v2h = get("QGenCharged_v2HadronIdx")
                    v2nearest = get("QGenCharged_v2NearestHeavyPrunedIdx")
                    v2hq = get("QGenCharged_v2HadronAssociationQuality")
                    v2HadronReason = get("QGenCharged_v2HadronAssociationReason")
                    v2q = get("QGenCharged_v2QuarkIdx")
                    v2quality = get("QGenCharged_v2AssociationQuality")
                    v2mem = get("QGenCharged_v2MembershipType")
                    v2certain = get("QGenCharged_v2FragmentationSourceCertain")
                    v2reason = get("QGenCharged_v2AssociationReason")
                    pt = get("QGenCharged_pt")
                    eta = get("QGenCharged_eta")
                    from_b = get("QGenCharged_isFromB")
                    from_prompt_d = get("QGenCharged_isFromPromptD")
                    from_btod = get("QGenCharged_isFromBtoD")
                    from_unknown_d = get("QGenCharged_v2IsFromUnknownD")
                    for hrow, (quark, quality) in enumerate(zip(hq, hquality)):
                        if quark < -1 or quark >= len(qflav):
                            raise AssertionError("QHadron V2 QFrag index out of bounds in event %s, row %d" % (get("event"), hrow))
                        if quality not in (0, 1, 2, 3):
                            raise AssertionError("unknown QHadron V2 quality %r in event %s, row %d" % (quality, get("event"), hrow))
                    prompt_qsystems = set()
                    accepted_unique_extra_by_q = Counter()
                    accepted_certain_extra_by_q = Counter()
                    accepted_daughters_by_hadron = Counter()
                    prompt_d_rows = set()
                    for hrow, (flav, origin, certain, quark) in enumerate(zip(hflav, horigin, hcertain, hq)):
                        if flav == 4 and origin == 1 and certain and 0 <= quark < len(qflav) and abs(qflav[quark]) == 4:
                            prompt_qsystems.add(quark)
                            prompt_d_rows.add(hrow)
                    for row, packed_idx in enumerate(packed):
                        old_class = "unique" if qidx_old[row] >= 0 and not old_amb[row] else ("ambiguous" if old_amb[row] else "unresolved")
                        new_class = {1: "unique", 2: "ambiguous", 3: "unresolved", 0: "missing"}.get(v2quality[row], "unresolved")
                        transitions[(old_class, new_class)] += 1
                        summary["old_HF_decay"] += bool(old_hf[row])
                        summary["new_HF_decay"] += (v2mem[row] == 1)
                        summary["particles_with_selected_HF_ancestry"] += (v2h[row] >= 0)
                        summary["particles_with_retained_nearest_heavy_ancestry"] += (v2nearest[row] >= 0)
                        summary["nearest_hadron_unique"] += (v2hq[row] == 1)
                        summary["nearest_hadrons_tied"] += (v2hq[row] == 2)
                        summary["old_ambiguous_to_new_HF_decay"] += bool(old_amb[row] and v2mem[row] == 1)
                        summary["old_HF_decay_to_new_ambiguous"] += bool(old_hf[row] and v2mem[row] == 3)
                        summary["unique_ancestry_with_QFrag"] += (v2quality[row] == 1 and v2q[row] >= 0)
                        summary["unique_ancestry_without_QFrag"] += (v2quality[row] == 1 and v2q[row] < 0)
                        summary["ambiguous_quark_ancestry"] += (v2quality[row] == 2)
                        summary["missing_or_unresolved_quark_ancestry"] += (v2quality[row] in (0, 3))
                        summary["additional_fragmentation_candidate"] += (v2mem[row] == 2)
                        summary["ambiguous_fragmentation_membership"] += (v2mem[row] == 3)
                        summary["unassigned_membership"] += (v2mem[row] == 0)
                        if from_b[row]: summary["B_decay_ancestry"] += 1
                        if from_prompt_d[row]: summary["prompt_D_ancestry"] += 1
                        if from_btod[row]: summary["B_to_D_ancestry"] += 1
                        if from_unknown_d[row]: summary["unknown_origin_D_ancestry"] += 1
                        issue = None
                        if v2h[row] >= len(hflav) or v2h[row] < -1:
                            issue = "v2HadronIdx outside event QHadron table"
                        elif v2hq[row] == 1 and v2h[row] < 0:
                            issue = "unique nearest heavy hadron has no selected QHadron row"
                        elif v2h[row] >= 0:
                            if v2hq[row] != 1:
                                issue = "valid v2HadronIdx without unique nearest-hadron quality"
                            elif v2q[row] != hq[v2h[row]]:
                                issue = "inherited quark index differs from QHadron V2 quark index"
                            elif v2quality[row] != hquality[v2h[row]]:
                                issue = "inherited association quality differs from QHadron V2 quality"
                            elif v2certain[row] != hcertain[v2h[row]]:
                                issue = "inherited source certainty differs from QHadron"
                        if v2q[row] >= len(qflav) or v2q[row] < -1:
                            issue = "v2QuarkIdx outside event QFrag table"
                        if v2mem[row] == 1 and v2hq[row] != 1 and not (v2hq[row] == 3 and v2nearest[row] >= 0):
                            issue = "HF-decay membership lacks a unique nearest hadron"
                        if v2hq[row] == 2 and (v2h[row] != -1 or v2q[row] != -1):
                            issue = "tied nearest hadrons must keep hadron and quark indices unresolved"
                        if v2hq[row] == 3 and v2HadronReason[row] == 3 and (v2nearest[row] < 0 or v2mem[row] != 1 or v2h[row] != -1):
                            issue = "unselected nearest heavy hadron must retain decay membership and no selected-row index"
                        if from_btod[row] and from_prompt_d[row]:
                            issue = "B-to-D descendant also marked as prompt-D descendant"
                        if from_unknown_d[row] and v2certain[row]:
                            issue = "unknown-origin D descendant promoted to certain source"
                        if issue:
                            summary["invariant_failures"] += 1
                            if len(examples) < args.max_examples:
                                examples.append({"file": str(path), "run": get("run"), "luminosityBlock": get("luminosityBlock"),
                                                 "event": get("event"), "packedGenIdx": packed_idx, "issue": issue,
                                                 "old_quarkIdx": qidx_old[row], "new_quarkIdx": v2q[row],
                                                 "new_hadronIdx": v2h[row], "hadronQuality": v2hq[row], "membership": v2mem[row],
                                                 "associationQuality": v2quality[row], "associationReason": v2reason[row]})
                        if v2mem[row] == 2:
                            extra_old = "unique" if qidx_old[row] >= 0 and not old_amb[row] else ("ambiguous" if old_amb[row] else "unresolved")
                            extra_new = "unique" if v2quality[row] == 1 and v2q[row] >= 0 else ("ambiguous" if v2quality[row] == 2 else "ambiguous/unresolved")
                            transitions[("additional:" + extra_old, "additional:" + extra_new)] += 1
                            summary["additional_old_unique"] += extra_old == "unique"
                            summary["additional_new_unique"] += extra_new == "unique"
                            summary["additional_old_unique_to_new_unique"] += extra_old == "unique" and extra_new == "unique"
                            summary["additional_old_unique_to_new_ambiguous_or_unresolved"] += extra_old == "unique" and extra_new != "unique"
                            summary["additional_old_ambiguous_to_new_unique"] += extra_old == "ambiguous" and extra_new == "unique"
                            summary["additional_old_ambiguous_to_new_ambiguous_or_unresolved"] += extra_old == "ambiguous" and extra_new != "unique"
                        if v2mem[row] == 2 and v2quality[row] == 1 and v2q[row] >= 0:
                            flavour = int(qpdg[v2q[row]])
                            names = {4: "c", -4: "anti-c", 5: "b", -5: "anti-b", 1: "d", -1: "anti-d", 2: "u", -2: "anti-u", 3: "s", -3: "anti-s"}
                            summary["extra_unique_flavour_" + names.get(flavour, str(flavour))] += 1
                            summary["extra_unique_quark_association"] += 1
                            summary["extra_unique_source_certain"] += bool(v2certain[row])
                            summary["extra_unique_source_uncertain"] += not bool(v2certain[row])
                            if pt[row] >= 0.4 and abs(eta[row]) <= 2.5 and v2q[row] in prompt_qsystems:
                                accepted_unique_extra_by_q[v2q[row]] += 1
                                summary["prompt_charm_additional_accepted_unique"] += 1
                                if v2certain[row]:
                                    accepted_certain_extra_by_q[v2q[row]] += 1
                                    summary["prompt_charm_additional_accepted_source_certain"] += 1
                                else:
                                    summary["prompt_charm_additional_accepted_unique_source_uncertain"] += 1
                        if v2mem[row] == 1 and v2q[row] >= 0 and v2certain[row] and pt[row] >= 0.4 and abs(eta[row]) <= 2.5:
                            hrow = v2h[row]
                            if hrow in prompt_d_rows:
                                accepted_daughters_by_hadron[hrow] += 1
                    unique_system_bins = Counter()
                    certain_system_bins = Counter()
                    d_bins = Counter()
                    for quark in prompt_qsystems:
                        unique_system_bins[min(accepted_unique_extra_by_q[quark], 3)] += 1
                        certain_system_bins[min(accepted_certain_extra_by_q[quark], 3)] += 1
                    for hrow in prompt_d_rows:
                        d_bins[min(accepted_daughters_by_hadron[hrow], 3)] += 1
                    summary["confirmed_prompt_charm_systems"] += len(prompt_qsystems)
                    for count in range(4):
                        suffix = "3plus" if count == 3 else str(count)
                        summary["prompt_systems_with_%s_unique_additional_accepted" % suffix] += unique_system_bins[count]
                        summary["prompt_systems_with_%s_source_certain_additional_accepted" % suffix] += certain_system_bins[count]
                        summary["prompt_D_with_%s_accepted_decay_descendants" % ("3plus" if count == 3 else count)] += d_bins[count]
                    for quark in sorted(prompt_qsystems):
                        drows = [hrow for hrow in prompt_d_rows if hq[hrow] == quark]
                        acceptance_rows.append({"file": str(path), "run": event_key[0], "luminosityBlock": event_key[1],
                                                "event": event_key[2], "qfragIdx": quark,
                                                "additionalAcceptedChargedUnique": accepted_unique_extra_by_q[quark],
                                                "additionalAcceptedChargedSourceCertain": accepted_certain_extra_by_q[quark],
                                                "additionalAcceptedChargedUniqueSourceUncertain": accepted_unique_extra_by_q[quark] - accepted_certain_extra_by_q[quark],
                                                "promptDCount": len(drows),
                                                "promptDObjectsWithExactlyTwoAcceptedChargedDaughters": sum(accepted_daughters_by_hadron[hrow] == 2 for hrow in drows),
                                                "sourceCertaintyNote": "unique count includes uncertain candidates; source-certain is the nested subset"})
        if reference is not None and seen_current_events != set(reference):
            missing_current = set(reference) - seen_current_events
            if missing_current:
                summary["invariant_failures"] += len(missing_current)
                for run, lumi, event in sorted(missing_current)[:max(0, args.max_examples - len(examples))]:
                    examples.append({"file": str(path), "run": run, "luminosityBlock": lumi, "event": event,
                                     "packedGenIdx": "", "issue": "reference event missing from updated output"})
        prefix = args.output_dir / stem
        write_csv(prefix.with_name(prefix.name + "_summary.csv"),
                  [{"metric": key, "value": value} for key, value in sorted(summary.items())], ["metric", "value"])
        write_csv(prefix.with_name(prefix.name + "_transitions.csv"),
                  [{"old_class": old, "new_class": new, "count": count} for (old, new), count in sorted(transitions.items())],
                  ["old_class", "new_class", "count"])
        write_csv(prefix.with_name(prefix.name + "_problematic.csv"), examples,
                  ["file", "run", "luminosityBlock", "event", "packedGenIdx", "issue", "old_quarkIdx", "new_quarkIdx", "new_hadronIdx", "hadronQuality", "membership", "associationQuality", "associationReason"])
        write_csv(prefix.with_name(prefix.name + "_prompt_charm_acceptance.csv"), acceptance_rows,
                  ["file", "run", "luminosityBlock", "event", "qfragIdx", "additionalAcceptedChargedUnique",
                   "additionalAcceptedChargedSourceCertain", "additionalAcceptedChargedUniqueSourceUncertain", "promptDCount",
                   "promptDObjectsWithExactlyTwoAcceptedChargedDaughters", "sourceCertaintyNote"])
        # Unique ancestry acceptance and the stricter source-certain subset are
        # reported independently to expose the effect of the certainty policy.
        if summary["invariant_failures"]:
            raise SystemExit("%s has %d invariant failures; see %s" % (path, summary["invariant_failures"], prefix.with_name(prefix.name + "_problematic.csv")))
        print("%s: %d events, %d QGenCharged rows; CSVs written with prefix %s" %
              (path, summary["events"], summary["stable_charged_particles"], prefix))
        if reference is None:
            print("  Pass --reference with the same events from the pre-update producer to verify selections and unchanged V2 hadron truth.")
        else:
            print("  Compared QFrag, QHadron, QGenCharged selection/order and QHadron V2 truth by event identity.")
        print("  Prompt-charm additional accepted: unique=%d, source-certain=%d, unique/source-uncertain=%d" %
              (summary["prompt_charm_additional_accepted_unique"],
               summary["prompt_charm_additional_accepted_source_certain"],
               summary["prompt_charm_additional_accepted_unique_source_uncertain"]))


if __name__ == "__main__":
    main()
