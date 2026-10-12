#!/usr/bin/env python3
"""Validate QRecoTrack output and produce matching/feasibility CSVs.

Requires uproot 5 and awkward. Example:
  python3 validate_quark_reco_track_matching.py --input HWWc_stage2.root \
      --sample HWWc --output-prefix HWWc_stage2
"""
import argparse
import csv
import math
from collections import Counter, defaultdict

import awkward as ak
import uproot


ACCEPTANCE_PT = 0.4
ACCEPTANCE_ABSETA = 2.5
EVENT_KEYS = ("run", "luminosityBlock", "event")
RECO_FIELDS = ("trackIdx", "pt", "eta", "phi", "charge", "isHighPurity", "genChargedIdx", "matchDR",
               "matchRelPt", "matchStatus", "nGenCandidates")
GEN_FIELDS = ("packedGenIdx", "pdgId", "charge", "pt", "eta", "phi", "hadronIdx", "quarkIdx",
              "associationQuality", "membershipType", "fragmentationSourceCertain", "isFromB",
              "isFromPromptD", "isFromBtoD")
HADRON_FIELDS = ("pdgId", "charmOrigin", "quarkIdx")
QUARK_FIELDS = ("pdgId",)


def accepted(pt, eta):
    return math.isfinite(float(pt)) and math.isfinite(float(eta)) and pt >= ACCEPTANCE_PT and abs(eta) <= ACCEPTANCE_ABSETA


def delta_phi(a, b):
    return math.remainder(float(a) - float(b), 2. * math.pi)


def candidate(reco, gen, max_dr, max_rel_pt):
    rpt, reta, rphi, rcharge = reco
    gpt, geta, gphi, gcharge = gen
    if rcharge == 0 or gcharge == 0 or rcharge != gcharge:
        return None
    values = (rpt, reta, rphi, gpt, geta, gphi)
    if not all(math.isfinite(float(x)) for x in values) or rpt <= 0 or gpt <= 0:
        return None
    dr = math.hypot(float(reta) - float(geta), delta_phi(rphi, gphi))
    relpt = abs(float(rpt) - float(gpt)) / float(gpt)
    return (dr, relpt) if dr <= max_dr and relpt <= max_rel_pt else None


def branch_table(arrays, prefix, fields, event):
    result = {}
    for field in fields:
        key = prefix + "_" + field
        if key not in arrays:
            raise ValueError("missing required branch " + key)
        result[field] = list(arrays[key][event])
    lengths = {len(value) for value in result.values()}
    if len(lengths) != 1:
        raise ValueError("inconsistent columns in " + prefix)
    return result


def add(summary, section, category, metric, value):
    summary.append((section, category, metric, value))


def quark_flavour(pdg):
    a = abs(int(pdg))
    if a == 4:
        return "charm" if pdg > 0 else "anticharm"
    if a == 5:
        return "bottom" if pdg > 0 else "antibottom"
    if 1 <= a <= 3:
        return "light"
    return "other_or_invalid"


def particle_type(pdg):
    return {211: "pion", 321: "kaon", 2212: "proton", 11: "electron", 13: "muon"}.get(abs(int(pdg)), "other")


def decay_category(gen, row, hadrons):
    if gen["membershipType"][row] != 1:
        return None
    if gen["isFromBtoD"][row] == 1:
        return "B_to_D"
    if gen["isFromPromptD"][row] == 1:
        return "confirmed_prompt_D"
    hadron = gen["hadronIdx"][row]
    if 0 <= hadron < len(hadrons["charmOrigin"]) and hadrons["charmOrigin"][hadron] == 3:
        return "unknown_origin_D_decay"
    if gen["isFromB"][row] == 1:
        return "B_hadron_decay"
    return "other_heavy_hadron_decay"


def csv_write(path, header, rows):
    with open(path, "w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="NanoAOD ROOT file containing QFrag, QHadron, QGenCharged, and QRecoTrack")
    parser.add_argument("--sample", required=True, choices=("HWWc", "ttbarhad"), help="sample label for CSV rows")
    parser.add_argument("--output-prefix", required=True, help="path prefix for CSV outputs")
    parser.add_argument("--maxDeltaR", type=float, default=0.03, help="matching cone used by the producer")
    parser.add_argument("--maxRelPt", type=float, default=0.5, help="relative-pT threshold used by the producer")
    parser.add_argument("--max-examples", type=int, default=100, help="maximum representative cases per case type")
    parser.add_argument("--max-events", type=int, default=None,
                        help="process only the first N Events entries (also limits ROOT reading); default: all events")
    args = parser.parse_args()
    if args.maxDeltaR <= 0 or args.maxRelPt <= 0 or args.max_examples < 0 or (args.max_events is not None and args.max_events <= 0):
        parser.error("thresholds must be positive, max-examples nonnegative, and max-events positive when specified")

    with uproot.open(args.input) as root_file:
        tree = root_file["Events"]
        input_events = tree.num_entries
        events_to_process = input_events if args.max_events is None else min(input_events, args.max_events)
        required = list(EVENT_KEYS)
        for prefix, fields in (("QRecoTrack", RECO_FIELDS), ("QGenCharged", GEN_FIELDS),
                               ("QHadron", HADRON_FIELDS), ("QFrag", QUARK_FIELDS)):
            required.extend(prefix + "_" + field for field in fields)
        missing = [name for name in required if name not in tree.keys()]
        if missing:
            raise ValueError("missing Events branches: " + ", ".join(missing))
        arrays = tree.arrays(required, entry_stop=events_to_process, library="ak", how=dict)

    summary = []
    systems = []
    hadron_rows = []
    cases = defaultdict(list)
    residual_dr = []
    residual_relpt = []
    total = Counter()
    total["input_events_in_file"] = input_events
    total["event_limit_requested"] = args.max_events if args.max_events is not None else -1
    total["events_processed"] = events_to_process
    heavy = defaultdict(Counter)
    heavy_bins = defaultdict(Counter)
    candidate_groups = defaultdict(Counter)
    flavour_groups = defaultdict(Counter)
    source_groups = defaultdict(Counter)
    purity_counts = Counter()

    def save_case(kind, row):
        if len(cases[kind]) < args.max_examples:
            cases[kind].append(row)

    for iev in range(len(arrays["event"])):
        run, lumi, event_id = (int(arrays[key][iev]) for key in EVENT_KEYS)
        reco = branch_table(arrays, "QRecoTrack", RECO_FIELDS, iev)
        gen = branch_table(arrays, "QGenCharged", GEN_FIELDS, iev)
        hadrons = branch_table(arrays, "QHadron", HADRON_FIELDS, iev)
        quarks = branch_table(arrays, "QFrag", QUARK_FIELDS, iev)
        ntracks, ngen, nhadrons, nquarks = len(reco["trackIdx"]), len(gen["pdgId"]), len(hadrons["pdgId"]), len(quarks["pdgId"])
        if reco["trackIdx"] != list(range(ntracks)):
            raise ValueError("trackIdx is not the positional index in event %d" % event_id)
        matched_by_gen = {}
        for it, igen in enumerate(reco["genChargedIdx"]):
            status = int(reco["matchStatus"][it])
            ncand = int(reco["nGenCandidates"][it])
            if status not in (0, 1, 2, 3) or ncand < 0:
                raise ValueError("invalid matchStatus/nGenCandidates in event %d track %d" % (event_id, it))
            if reco["isHighPurity"][it] not in (0, 1):
                raise ValueError("isHighPurity must be encoded as 0 or 1")
            if igen < 0:
                if igen != -1 or reco["matchDR"][it] != -1 or reco["matchRelPt"][it] != -1:
                    raise ValueError("unmatched track does not use documented -1 sentinels")
                if ncand > 0 and status != 3:
                    raise ValueError("eligible unmatched track must use one-to-one-conflict status 3")
            else:
                if igen >= ngen or igen in matched_by_gen:
                    raise ValueError("invalid or duplicated event-local QGenCharged reference")
                matched_by_gen[int(igen)] = it
                if status not in (1, 2):
                    raise ValueError("matched track has unmatched status")
            if (ncand == 0) != (status == 0):
                raise ValueError("status 0 must be equivalent to zero admissible candidates")
            if status == 3 and (igen != -1 or ncand == 0):
                raise ValueError("status 3 must be an unmatched track with eligible candidates")

        reco_kin = [(reco["pt"][i], reco["eta"][i], reco["phi"][i], reco["charge"][i]) for i in range(ntracks)]
        gen_kin = [(gen["pt"][j], gen["eta"][j], gen["phi"][j], gen["charge"][j]) for j in range(ngen)]
        edges = defaultdict(list)
        gen_degree = Counter()
        tracks_per_gen = defaultdict(list)
        for it, rkin in enumerate(reco_kin):
            for ig, gkin in enumerate(gen_kin):
                residual = candidate(rkin, gkin, args.maxDeltaR, args.maxRelPt)
                if residual is not None:
                    edges[it].append((ig, residual[0], residual[1]))
                    gen_degree[ig] += 1
                    tracks_per_gen[ig].append(it)
        if any(int(reco["nGenCandidates"][i]) != len(edges[i]) for i in range(ntracks)):
            raise ValueError("nGenCandidates differs from validator candidate count in event %d" % event_id)
        for it, igen in enumerate(reco["genChargedIdx"]):
            if igen >= 0:
                possible = {ig for ig, _, _ in edges[it]}
                if igen not in possible:
                    raise ValueError("assigned pair fails candidate requirements")
                dr, relpt = candidate(reco_kin[it], gen_kin[igen], args.maxDeltaR, args.maxRelPt)
                if not math.isclose(dr, float(reco["matchDR"][it]), rel_tol=1e-5, abs_tol=1e-6):
                    raise ValueError("stored matchDR disagrees with kinematics")
                if not math.isclose(relpt, float(reco["matchRelPt"][it]), rel_tol=1e-5, abs_tol=1e-6):
                    raise ValueError("stored matchRelPt disagrees with kinematics")
                if status := int(reco["matchStatus"][it]):
                    expected_competition = len(edges[it]) > 1 or gen_degree[igen] > 1
                    if status != (2 if expected_competition else 1):
                        raise ValueError("matchStatus does not reflect candidate competition")
                residual_dr.append(dr)
                residual_relpt.append(relpt)

        for idx, (hadron, quark) in enumerate(zip(gen["hadronIdx"], gen["quarkIdx"])):
            if hadron < -1 or hadron >= nhadrons or quark < -1 or quark >= nquarks:
                raise ValueError("invalid event-local QGenCharged reference in event %d row %d" % (event_id, idx))
            if gen["membershipType"][idx] == 2 and gen["fragmentationSourceCertain"][idx] != 0:
                raise ValueError("Stage 1 additional candidates must retain source-certain=0")
        for idx, quark in enumerate(hadrons["quarkIdx"]):
            if quark < -1 or quark >= nquarks:
                raise ValueError("invalid event-local QHadron quarkIdx")

        accepted_tracks = [accepted(reco["pt"][i], reco["eta"][i]) for i in range(ntracks)]
        accepted_gens = [accepted(gen["pt"][i], gen["eta"][i]) for i in range(ngen)]
        matched_gen_ids = set(matched_by_gen)
        local_counts = Counter()
        local_counts["events"] += 1
        local_counts["reconstructed_tracks"] += ntracks
        local_counts["QGenCharged_particles"] += ngen
        local_counts["accepted_QGenCharged_particles"] += sum(accepted_gens)
        local_counts["eligible_reco_gen_pairs"] += sum(len(edges[i]) for i in range(ntracks))
        local_counts["matched_reco_tracks"] += len(matched_by_gen)
        local_counts["unmatched_reco_tracks"] += ntracks - len(matched_by_gen)
        local_counts["matched_generator_particles"] += len(matched_gen_ids)
        local_counts["unmatched_generator_particles"] += ngen - len(matched_gen_ids)
        local_counts["tracks_zero_candidates"] += sum(len(edges[i]) == 0 for i in range(ntracks))
        local_counts["tracks_one_candidate"] += sum(len(edges[i]) == 1 for i in range(ntracks))
        local_counts["tracks_multiple_candidates"] += sum(len(edges[i]) > 1 for i in range(ntracks))
        local_counts["tracks_lost_one_to_one_conflict"] += sum(int(x) == 3 for x in reco["matchStatus"])
        local_counts["accepted_reco_tracks"] += sum(accepted_tracks)
        local_counts["accepted_reco_matched_tracks"] += sum(accepted_tracks[i] and reco["genChargedIdx"][i] >= 0 for i in range(ntracks))
        local_counts["accepted_reco_unmatched_tracks"] += sum(accepted_tracks[i] and reco["genChargedIdx"][i] < 0 for i in range(ntracks))
        local_counts["accepted_reco_eligible_pairs"] += sum(len(edges[i]) for i in range(ntracks) if accepted_tracks[i])
        local_counts["accepted_reco_zero_candidate_tracks"] += sum(accepted_tracks[i] and len(edges[i]) == 0 for i in range(ntracks))
        local_counts["accepted_reco_one_candidate_tracks"] += sum(accepted_tracks[i] and len(edges[i]) == 1 for i in range(ntracks))
        local_counts["accepted_reco_multiple_candidate_tracks"] += sum(accepted_tracks[i] and len(edges[i]) > 1 for i in range(ntracks))
        local_counts["accepted_reco_conflict_tracks"] += sum(accepted_tracks[i] and int(reco["matchStatus"][i]) == 3 for i in range(ntracks))
        for metric, value in local_counts.items():
            total[metric] += value

        purity_counts["high_purity_tracks"] += sum(bool(x) for x in reco["isHighPurity"])
        purity_counts["high_purity_matched_tracks"] += sum(bool(reco["isHighPurity"][i]) and reco["genChargedIdx"][i] >= 0 for i in range(ntracks))
        purity_counts["high_purity_accepted_tracks"] += sum(bool(reco["isHighPurity"][i]) and accepted_tracks[i] for i in range(ntracks))
        purity_counts["high_purity_accepted_matched_tracks"] += sum(bool(reco["isHighPurity"][i]) and accepted_tracks[i] and reco["genChargedIdx"][i] >= 0 for i in range(ntracks))

        # Per-category candidate availability, matches, and the Stage 1
        # certainty bit remain independent dimensions.
        for ig in range(ngen):
            membership = int(gen["membershipType"][ig])
            quality = int(gen["associationQuality"][ig])
            if membership == 2 and quality == 1 and gen["quarkIdx"][ig] >= 0:
                category = "unique_additional"
                qidx = int(gen["quarkIdx"][ig])
                qflavour = quark_flavour(quarks["pdgId"][qidx])
                matched_track = matched_by_gen.get(ig)
                accepted_reco_match = matched_track is not None and accepted(
                    reco["pt"][matched_track], reco["eta"][matched_track])
                flavour_groups[qflavour]["gen_particles"] += 1
                flavour_groups[qflavour]["matched_gen_particles"] += ig in matched_by_gen
                flavour_groups[qflavour]["accepted_gen_particles"] += accepted_gens[ig]
                flavour_groups[qflavour]["accepted_gen_matched_reco_tracks"] += accepted_gens[ig] and accepted_reco_match
            elif membership == 3 or quality == 2:
                category = "ambiguous_charged"
            elif membership == 1:
                category = "heavy_hadron_decay"
            else:
                category = "unassigned_charged"
            candidate_groups[category]["gen_particles"] += 1
            candidate_groups[category]["matched_gen_particles"] += ig in matched_by_gen
            candidate_groups[category]["accepted_gen_particles"] += accepted_gens[ig]
            matched_track = matched_by_gen.get(ig)
            accepted_reco_match = matched_track is not None and accepted(
                reco["pt"][matched_track], reco["eta"][matched_track])
            candidate_groups[category]["accepted_gen_matched_reco_tracks"] += accepted_gens[ig] and accepted_reco_match
            source_groups[category]["source_certain_0"] += gen["fragmentationSourceCertain"][ig] == 0
            source_groups[category]["source_certain_1"] += gen["fragmentationSourceCertain"][ig] == 1

            decay = decay_category(gen, ig, hadrons)
            if decay is not None:
                matched_track = matched_by_gen.get(ig)
                reco_acc = matched_track is not None and accepted(reco["pt"][matched_track], reco["eta"][matched_track])
                heavy[decay]["generator_particles"] += 1
                heavy[decay]["matched_reco_tracks"] += matched_track is not None
                heavy[decay]["generator_accepted_particles"] += accepted_gens[ig]
                heavy[decay]["accepted_gen_matched_reco_tracks"] += accepted_gens[ig] and reco_acc
                ptbin = "<1" if gen["pt"][ig] < 1 else "1-2" if gen["pt"][ig] < 2 else "2-5" if gen["pt"][ig] < 5 else "5+"
                etabin = "0-0.8" if abs(gen["eta"][ig]) < 0.8 else "0.8-1.5" if abs(gen["eta"][ig]) < 1.5 else "1.5-2.5" if abs(gen["eta"][ig]) <= 2.5 else ">2.5"
                ptype = particle_type(gen["pdgId"][ig])
                for dimension, label in (("gen_pt", ptbin), ("gen_abs_eta", etabin), ("particle_type", ptype)):
                    counter = heavy_bins[(decay, dimension, label)]
                    counter["particles"] += 1
                    counter["matched_reco_tracks"] += matched_track is not None
                    counter["generator_accepted_particles"] += accepted_gens[ig]
                    counter["accepted_gen_matched_reco_tracks"] += accepted_gens[ig] and reco_acc

        # Reconstruct system and per-hadron views from the unchanged Stage 1
        # event-local rows; count each charged generator row at most once.
        prompt_by_q = defaultdict(list)
        for ih, (origin, qidx) in enumerate(zip(hadrons["charmOrigin"], hadrons["quarkIdx"])):
            if origin == 1 and qidx >= 0:
                prompt_by_q[int(qidx)].append(ih)
        for qidx, prompt_hadrons in prompt_by_q.items():
            prompt_set = set(prompt_hadrons)
            dgen = [ig for ig in range(ngen) if gen["membershipType"][ig] == 1 and
                    gen["isFromPromptD"][ig] == 1 and gen["hadronIdx"][ig] in prompt_set and accepted_gens[ig]]
            dmatched = [ig for ig in dgen if ig in matched_by_gen and
                        accepted(reco["pt"][matched_by_gen[ig]], reco["eta"][matched_by_gen[ig]])]
            extras = [ig for ig in range(ngen) if gen["membershipType"][ig] == 2 and
                      gen["associationQuality"][ig] == 1 and gen["quarkIdx"][ig] == qidx and accepted_gens[ig]]
            extra_matched = [ig for ig in extras if ig in matched_by_gen and
                             accepted(reco["pt"][matched_by_gen[ig]], reco["eta"][matched_by_gen[ig]])]
            qpdg = quarks["pdgId"][qidx]
            system_row = [args.sample, run, lumi, event_id, qidx, int(qpdg), ";".join(map(str, prompt_hadrons)),
                          len(prompt_hadrons), len(dgen), len(dmatched), len(extras), len(extra_matched),
                          sum(gen["fragmentationSourceCertain"][ig] == 0 for ig in extras),
                          sum(gen["fragmentationSourceCertain"][ig] == 1 for ig in extras),
                          ";".join(map(str, dgen)), ";".join(map(str, extras))]
            systems.append(system_row)
            if len(dgen) == 2:
                total["generator_defined_exactly_two_D_track_systems"] += 1
                if extras:
                    total["generator_defined_two_D_systems_with_additional_gen_candidate"] += 1
                if extra_matched:
                    total["generator_defined_two_D_systems_with_matched_additional_track"] += 1
            if len(dmatched) == 2:
                total["reconstructed_exactly_two_D_track_systems"] += 1
                if extra_matched:
                    total["reconstructed_two_D_systems_with_matched_additional_track"] += 1
            if len(dgen) == 2 and len(dmatched) < 2:
                save_case("prompt_D_loses_track", [args.sample, run, lumi, event_id, ",".join(map(str, prompt_hadrons)),
                    qidx, "", "", ";".join(map(str, dgen)), ";".join(str(matched_by_gen.get(x, -1)) for x in dgen),
                    "accepted D gen rows=%d; matched reco-accepted D tracks=%d" % (len(dgen), len(dmatched))])
            lost_extras = [ig for ig in extras if ig not in matched_by_gen or
                           not accepted(reco["pt"][matched_by_gen[ig]], reco["eta"][matched_by_gen[ig]])]
            if lost_extras:
                save_case("additional_candidate_not_reconstructed", [args.sample, run, lumi, event_id, "",
                    qidx, "", "", ";".join(map(str, lost_extras)),
                    ";".join(str(matched_by_gen.get(x, -1)) for x in lost_extras),
                    "unique additional gen candidates not matched to accepted reco tracks"])
            for ih in prompt_hadrons:
                hgen = [ig for ig in dgen if gen["hadronIdx"][ig] == ih]
                hmatched = [ig for ig in hgen if ig in matched_by_gen and
                            accepted(reco["pt"][matched_by_gen[ig]], reco["eta"][matched_by_gen[ig]])]
                hadron_rows.append([args.sample, run, lumi, event_id, qidx, int(qpdg), ih,
                                    int(hadrons["pdgId"][ih]), len(hgen), len(hmatched),
                                    ";".join(map(str, hgen)), ";".join(str(matched_by_gen.get(x, -1)) for x in hgen)])

        for igen, candidate_tracks in tracks_per_gen.items():
            if len(candidate_tracks) > 1:
                save_case("generator_with_competing_tracks", [args.sample, run, lumi, event_id,
                    gen["hadronIdx"][igen], gen["quarkIdx"][igen], matched_by_gen.get(igen, -1), igen,
                    str(igen), ";".join(map(str, candidate_tracks)), "eligible tracks=%d; assigned track=%s" %
                    (len(candidate_tracks), matched_by_gen.get(igen, -1))])
        for it, gen_candidates in edges.items():
            if len(gen_candidates) > 1:
                save_case("track_with_competing_generators", [args.sample, run, lumi, event_id,
                    "", "", it, reco["genChargedIdx"][it], ";".join(str(x[0]) for x in gen_candidates), "",
                    "eligible QGenCharged rows=%d; assigned=%s" % (len(gen_candidates), reco["genChargedIdx"][it])])
            if int(reco["matchStatus"][it]) == 3:
                save_case("track_lost_one_to_one_conflict", [args.sample, run, lumi, event_id, "", "", it, -1,
                    ";".join(str(x[0]) for x in gen_candidates), "",
                    "eligible QGenCharged rows=%d; no assignment" % len(gen_candidates)])

    # Compact count and fraction tables.
    for key, value in sorted(total.items()):
        add(summary, "overall", "all", key, value)
    for key, value in sorted(purity_counts.items()):
        add(summary, "high_purity", "all", key, value)
    add(summary, "high_purity", "all", "all_high_purity_track_match_fraction",
        purity_counts["high_purity_matched_tracks"] / purity_counts["high_purity_tracks"]
        if purity_counts["high_purity_tracks"] else 0.)
    add(summary, "high_purity", "all", "accepted_high_purity_match_fraction",
        purity_counts["high_purity_accepted_matched_tracks"] / purity_counts["high_purity_accepted_tracks"]
        if purity_counts["high_purity_accepted_tracks"] else 0.)
    for category, values in sorted(heavy.items()):
        for key, value in sorted(values.items()):
            add(summary, "heavy_decay", category, key, value)
        denominator = values["generator_accepted_particles"]
        add(summary, "heavy_decay", category, "all_generator_to_reco_match_fraction",
            values["matched_reco_tracks"] / values["generator_particles"] if values["generator_particles"] else 0.)
        add(summary, "heavy_decay", category, "accepted_gen_to_accepted_reco_match_fraction",
            values["accepted_gen_matched_reco_tracks"] / denominator if denominator else 0.)
    for (category, dimension, label), values in sorted(heavy_bins.items()):
        for key, value in sorted(values.items()):
            add(summary, "heavy_decay_bins_" + dimension, category + ":" + label, key, value)
        add(summary, "heavy_decay_bins_" + dimension, category + ":" + label, "all_generator_to_reco_match_fraction",
            values["matched_reco_tracks"] / values["particles"] if values["particles"] else 0.)
        add(summary, "heavy_decay_bins_" + dimension, category + ":" + label, "accepted_gen_to_accepted_reco_match_fraction",
            values["accepted_gen_matched_reco_tracks"] / values["generator_accepted_particles"]
            if values["generator_accepted_particles"] else 0.)
    for category, values in sorted(candidate_groups.items()):
        for key, value in sorted(values.items()):
            add(summary, "additional_candidate_class", category, key, value)
        denominator = values["accepted_gen_particles"]
        add(summary, "additional_candidate_class", category, "all_generator_to_reco_match_fraction",
            values["matched_gen_particles"] / values["gen_particles"] if values["gen_particles"] else 0.)
        add(summary, "additional_candidate_class", category, "accepted_gen_to_accepted_reco_match_fraction",
            values["accepted_gen_matched_reco_tracks"] / denominator if denominator else 0.)
    for flavour, values in sorted(flavour_groups.items()):
        for key, value in sorted(values.items()):
            add(summary, "unique_additional_quark_flavour", flavour, key, value)
        add(summary, "unique_additional_quark_flavour", flavour, "all_generator_to_reco_match_fraction",
            values["matched_gen_particles"] / values["gen_particles"] if values["gen_particles"] else 0.)
        add(summary, "unique_additional_quark_flavour", flavour, "accepted_gen_to_accepted_reco_match_fraction",
            values["accepted_gen_matched_reco_tracks"] / values["accepted_gen_particles"]
            if values["accepted_gen_particles"] else 0.)
    for category, values in sorted(source_groups.items()):
        for key, value in sorted(values.items()):
            add(summary, "source_certainty", category, key, value)
    gen_den = total["generator_defined_exactly_two_D_track_systems"]
    reco_den = total["reconstructed_exactly_two_D_track_systems"]
    add(summary, "two_track_prompt_D", "generator_defined", "systems_with_exactly_two_accepted_D_gen_particles", gen_den)
    add(summary, "two_track_prompt_D", "generator_defined", "with_at_least_one_accepted_additional_gen_candidate",
        total["generator_defined_two_D_systems_with_additional_gen_candidate"])
    add(summary, "two_track_prompt_D", "generator_defined", "with_at_least_one_matched_additional_reco_track",
        total["generator_defined_two_D_systems_with_matched_additional_track"])
    add(summary, "two_track_prompt_D", "generator_defined", "Stage1_baseline_fraction",
        total["generator_defined_two_D_systems_with_additional_gen_candidate"] / gen_den if gen_den else 0.)
    add(summary, "two_track_prompt_D", "generator_defined", "matched_additional_reco_track_fraction",
        total["generator_defined_two_D_systems_with_matched_additional_track"] / gen_den if gen_den else 0.)
    add(summary, "two_track_prompt_D", "reconstructed_defined", "systems_with_exactly_two_matched_D_reco_tracks", reco_den)
    add(summary, "two_track_prompt_D", "reconstructed_defined", "with_at_least_one_matched_additional_reco_track",
        total["reconstructed_two_D_systems_with_matched_additional_track"])
    add(summary, "two_track_prompt_D", "reconstructed_defined", "matched_additional_reco_track_fraction",
        total["reconstructed_two_D_systems_with_matched_additional_track"] / reco_den if reco_den else 0.)
    expected = (239, 66) if args.sample == "HWWc" else (227, 99)
    observed = (gen_den, total["generator_defined_two_D_systems_with_additional_gen_candidate"])
    add(summary, "two_track_prompt_D", "generator_defined", "previous_Stage1_reference_denominator", expected[0])
    add(summary, "two_track_prompt_D", "generator_defined", "previous_Stage1_reference_numerator", expected[1])
    add(summary, "two_track_prompt_D", "generator_defined", "previous_Stage1_reference_matches_current", int(observed == expected))

    prefix = args.output_prefix
    csv_write(prefix + "_summary.csv", ["sample", "section", "category", "metric", "value"],
              [[args.sample, *row] for row in summary])
    csv_write(prefix + "_prompt_charm_systems.csv",
              ["sample", "run", "luminosityBlock", "event", "qfragIdx", "quarkPdgId", "qhadronRows",
               "nPromptDHadrons", "nAcceptedPromptDDecayGen", "nMatchedPromptDRecoTracks",
               "nAcceptedUniqueAdditionalGen", "nMatchedUniqueAdditionalRecoTracks",
               "nSourceUncertainAdditionalGen", "nSourceCertainAdditionalGen", "promptDGenRows", "additionalGenRows"], systems)
    csv_write(prefix + "_prompt_d_hadrons.csv",
              ["sample", "run", "luminosityBlock", "event", "qfragIdx", "quarkPdgId", "qhadronIdx",
               "hadronPdgId", "nAcceptedDecayGen", "nMatchedDecayRecoTracks", "genChargedRows", "matchedTrackRows"], hadron_rows)

    case_header = ["sample", "run", "luminosityBlock", "event", "qhadronIdx", "qfragIdx", "trackIdx",
                   "genChargedIdx", "relatedGenRows", "relatedTrackRows", "details"]
    case_output = []
    for case_type, rows in sorted(cases.items()):
        for row in rows:
            case_output.append([*row[:4], case_type, *row[4:]])
    csv_write(prefix + "_representative_cases.csv", ["sample", "run", "luminosityBlock", "event", "caseType",
              *case_header[4:]], case_output)

    def write_hist(path, values, maximum, bins):
        width = maximum / bins
        counts = Counter(min(int(value / width), bins - 1) for value in values)
        rows = [[i * width, (i + 1) * width, counts[i]] for i in range(bins)]
        csv_write(path, ["low", "high", "count"], rows)
    write_hist(prefix + "_match_dr_hist.csv", residual_dr, args.maxDeltaR, 30)
    write_hist(prefix + "_match_relpt_hist.csv", residual_relpt, args.maxRelPt, 20)
    print("Processed %d of %d Events entries%s" % (
        events_to_process, input_events,
        " (first-entry prefix limited by --max-events)" if events_to_process < input_events else ""))
    print("Wrote %s_{summary,prompt_charm_systems,prompt_d_hadrons,representative_cases,match_dr_hist,match_relpt_hist}.csv" % prefix)
    for row in summary:
        if row[0] in ("overall", "two_track_prompt_D"):
            print("%s,%s,%s,%s" % row)
    if args.sample == "HWWc":
        print("Stage 1 baseline observed: %d/%d; previous reference: 66/239" %
              (total["generator_defined_two_D_systems_with_additional_gen_candidate"], gen_den))
    else:
        print("Stage 1 baseline observed: %d/%d; previous reference: 99/227" %
              (total["generator_defined_two_D_systems_with_additional_gen_candidate"], gen_den))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
