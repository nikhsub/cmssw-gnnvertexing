#!/usr/bin/env python3
"""Build HZZc Stage 3 truth, neighborhood, and training-feasibility CSVs.

Requires uproot 5, awkward, NumPy, and pandas. Tables are joined only by their
event-local row references; run/lumi/event are retained in every output row.
"""

from __future__ import annotations

import argparse
import math
import os
from collections import Counter, defaultdict

import numpy as np

from neighborhood_analysis import delta_r, direction_from_tracks, upper_triangle_distances, nearest_neighbor_edge_pairs
from track_associations import additional_gen_rows, heavy_decay_gen_rows, matched_reco_by_gen
from truth_objects import (
    ASSOCIATION_UNIQUE, GEN_FIELDS, HADRON_FIELDS, MEMBERSHIP_ADDITIONAL,
    MEMBERSHIP_HEAVY_DECAY,
    QUARK_FIELDS, RECO_FIELDS, accepted, n_category, origin_category, reco_confident,
    successful_match, table_rows, validate_tables,
)

RADII = (0.2, 0.4, 0.6, 0.8, 1.0)
NEIGHBORHOOD_RADII = RADII + ("all",)
TREE_NAME = "Events"
EVENT_FIELDS = ("run", "luminosityBlock", "event")
TABLE_FIELDS = {
    "QFrag": QUARK_FIELDS,
    "QHadron": HADRON_FIELDS,
    "QGenCharged": GEN_FIELDS,
    "QRecoTrack": RECO_FIELDS,
}


def _bin(value, edges, labels):
    value = float(value)
    for edge, label in zip(edges[1:], labels):
        if value < edge:
            return label
    return labels[-1]


def _dr(row_a, row_b):
    return float(delta_r(row_a["eta"], row_a["phi"], row_b["eta"], row_b["phi"]))


def _wrap_phi(value):
    return (float(value) + math.pi) % (2.0 * math.pi) - math.pi


def _fragmentation_observables(had, quark):
    """Return unclipped three-momentum observables and finite-input status."""
    vals = [float(had[k]) for k in ("pt", "eta", "phi")] + [float(quark[k]) for k in ("pt", "eta", "phi")]
    if not all(math.isfinite(v) for v in vals):
        return {}, False
    pt_d, eta_d, phi_d, pt_q, eta_q, phi_q = vals
    if pt_d <= 0 or pt_q <= 0:
        return {}, False
    try:
        pd = np.asarray((pt_d * math.cos(phi_d), pt_d * math.sin(phi_d), pt_d * math.sinh(eta_d)), dtype=float)
        pq = np.asarray((pt_q * math.cos(phi_q), pt_q * math.sin(phi_q), pt_q * math.sinh(eta_q)), dtype=float)
    except OverflowError:
        return {}, False
    q2 = float(np.dot(pq, pq))
    qnorm = math.sqrt(q2) if q2 > 0 else 0.0
    if not math.isfinite(qnorm) or qnorm <= 0 or not np.all(np.isfinite(pd)) or not np.all(np.isfinite(pq)):
        return {}, False
    dot = float(np.dot(pd, pq))
    values = {
        "zProj": dot / q2,
        "zT": pt_d / pt_q,
        "deltaR_DQuark": math.hypot(eta_d - eta_q, _wrap_phi(phi_d - phi_q)),
        "kTrel": float(np.linalg.norm(np.cross(pd, pq))) / qnorm,
    }
    return values, all(math.isfinite(v) for v in values.values())


def _rankdata(values):
    """Small dependency-free average ranks for Spearman correlations."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def _correlation(x, y, rank=False):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(good)) < 3:
        return math.nan
    x, y = x[good], y[good]
    if rank:
        x, y = _rankdata(x), _rankdata(y)
    if np.std(x) == 0 or np.std(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _auxiliary_summary(rows):
    """Summary statistics/correlations; values are retained without clipping."""
    output = []
    variables = ("zProj", "zT", "deltaR_DQuark", "kTrel")
    for origin in ("D_prompt", "D_fromB", "D_unknown"):
        population = [r for r in rows if r["origin_category"] == origin]
        eligible = [r for r in population if int(r["auxiliary_target_eligible"] if origin == "D_prompt" else r["observable_inputs_valid"])]
        for ncat in ("N1", "N2", "N3+"):
            group = [r for r in population if r["N_category"] == ncat or (ncat == "N3+" and r["N_category"] == "N3")]
            valid = [r for r in group if int(r["auxiliary_target_eligible"] if origin == "D_prompt" else r["observable_inputs_valid"])]
            for var in variables:
                values = np.asarray([float(r[var]) for r in valid], dtype=float)
                values = values[np.isfinite(values)]
                count = len(values)
                q25, q75 = (np.percentile(values, [25, 75]) if count else (math.nan, math.nan))
                outliers = sum(int(r[var+"_outlier_tukey"]) == 1 for r in valid)
                output.append({
                    "origin_category": origin, "N_category": ncat, "observable": var,
                    "summary_type": "distribution", "pt_bin": "all",
                    "objects": len(group), "valid_objects": count,
                    "valid_fraction": count / len(group) if group else math.nan,
                    "mean": float(np.mean(values)) if count else math.nan,
                    "median": float(np.median(values)) if count else math.nan,
                    "q05": float(np.percentile(values, 5)) if count else math.nan,
                    "q95": float(np.percentile(values, 95)) if count else math.nan,
                    "outlier_count": outliers, "outlier_fraction": outliers/count if count else math.nan,
                    "pearson": math.nan, "spearman": math.nan,
                    "note": "Tukey 1.5-IQR outlier flag; raw observable is not clipped. Prompt-D rows require unique valid charm QFrag; other origins are diagnostic only."})
                if var == "zProj":
                    physical = sum(int(r["zProj_outside_0_1"]) for r in valid)
                    output.append({"origin_category": origin, "N_category": ncat, "observable": var,
                        "summary_type": "physical_boundary_outliers_outside_0_1", "pt_bin": "all",
                        "objects": len(group), "valid_objects": count, "valid_fraction": count/len(group) if group else math.nan,
                        "mean": math.nan, "median": math.nan, "q05": math.nan, "q95": math.nan,
                        "outlier_count": physical, "outlier_fraction": physical/count if count else math.nan,
                        "pearson": math.nan, "spearman": math.nan,
                        "note": "Diagnostic physical-range flag only; zProj values are retained unmodified."})
                if var == "zT":
                    physical = sum(int(r["zT_above_1"]) for r in valid)
                    output.append({"origin_category": origin, "N_category": ncat, "observable": var,
                        "summary_type": "physical_boundary_outliers_above_1", "pt_bin": "all",
                        "objects": len(group), "valid_objects": count, "valid_fraction": count/len(group) if group else math.nan,
                        "mean": math.nan, "median": math.nan, "q05": math.nan, "q95": math.nan,
                        "outlier_count": physical, "outlier_fraction": physical/count if count else math.nan,
                        "pearson": math.nan, "spearman": math.nan,
                        "note": "Diagnostic zT>1 flag only; values are retained unmodified."})
                # Correlate each target with available reconstructed or angular observables.
                for label, predictor in (("D_pt", "hadron_pt"), ("additional_candidate_count", "accepted_additional_reco_candidate_count"),
                                         ("additional_candidate_mean_deltaR_R0p4", "additional_candidate_mean_deltaR_R0p4"),
                                         ("local_track_multiplicity_R0p4", "local_reco_track_multiplicity_R0p4")):
                    output.append({
                        "origin_category": origin, "N_category": ncat, "observable": var,
                        "summary_type": "correlation_vs_"+label, "pt_bin": "all",
                        "objects": len(group), "valid_objects": count, "valid_fraction": count/len(group) if group else math.nan,
                        "mean": math.nan, "median": math.nan, "q05": math.nan, "q95": math.nan,
                        "outlier_count": "", "outlier_fraction": math.nan,
                        "pearson": _correlation([r[var] for r in valid], [r[predictor] for r in valid]),
                        "spearman": _correlation([r[var] for r in valid], [r[predictor] for r in valid], rank=True),
                        "note": "Pairwise finite observations; correlations require at least 3 objects and nonconstant inputs."})
                for low, high, label in ((0., 2., "0-2"), (2., 4., "2-4"), (4., 8., "4-8"), (8., math.inf, "8+")):
                    binned = [r for r in valid if low <= float(r["hadron_pt"]) < high]
                    output.append({
                        "origin_category": origin, "N_category": ncat, "observable": var,
                        "summary_type": "within_pt_bin_correlation_vs_additional_candidate_count", "pt_bin": label,
                        "objects": len(binned), "valid_objects": len(binned), "valid_fraction": math.nan,
                        "mean": math.nan, "median": math.nan, "q05": math.nan, "q95": math.nan,
                        "outlier_count": "", "outlier_fraction": math.nan,
                        "pearson": _correlation([r[var] for r in binned], [r["accepted_additional_reco_candidate_count"] for r in binned]),
                        "spearman": _correlation([r[var] for r in binned], [r["accepted_additional_reco_candidate_count"] for r in binned], rank=True),
                        "note": "Tests persistence after coarse D-pT stratification; interpret only when n>=3."})
    return output


def _new_metrics():
    return Counter()


def analyze_event(tables: dict, event_id: dict, source_file: str, require_high_purity: bool = False):
    """Analyze a single event. Kept pure of ROOT I/O for synthetic testing."""
    validate_tables(tables)
    quarks = table_rows(tables["QFrag"])
    hadrons = table_rows(tables["QHadron"])
    gens = table_rows(tables["QGenCharged"])
    recos = table_rows(tables["QRecoTrack"])
    matched = matched_reco_by_gen(recos)
    reco_by_track = {int(r["trackIdx"]): r for r in recos}
    reco_accept = lambda r: accepted(r["pt"], r["eta"]) and (not require_high_purity or int(r["isHighPurity"]) == 1)
    accepted_recos = {int(r["trackIdx"]): r for r in recos if reco_accept(r)}
    event_key = {key: int(value) for key, value in event_id.items()}
    identity = "%d:%d:%d" % (event_key["run"], event_key["luminosityBlock"], event_key["event"])
    base = {**event_key, "event_id": identity, "source_file": source_file}
    metrics = _new_metrics()
    metrics["events"] = 1
    metrics["QFrag_rows"] = len(quarks)
    metrics["QHadron_rows"] = len(hadrons)
    metrics["QGenCharged_rows"] = len(gens)
    metrics["QRecoTrack_rows"] = len(recos)
    for quality in (0, 1, 2, 3):
        metrics["QGenCharged_associationQuality_%d" % quality] += sum(int(r["associationQuality"]) == quality for r in gens)
    metrics["QGenCharged_source_certain"] += sum(int(r["fragmentationSourceCertain"]) == 1 for r in gens)
    metrics["QGenCharged_source_uncertain"] += sum(int(r["fragmentationSourceCertain"]) == 0 for r in gens)
    metrics["QGenCharged_additional_membership"] += sum(int(r["membershipType"]) == MEMBERSHIP_ADDITIONAL for r in gens)
    metrics["QGenCharged_heavy_decay_membership"] += sum(int(r["membershipType"]) == MEMBERSHIP_HEAVY_DECAY for r in gens)
    metrics["accepted_reco_nodes"] = len(accepted_recos)
    metrics["accepted_reco_high_purity_nodes"] = sum(int(r["isHighPurity"]) == 1 for r in accepted_recos.values())
    metrics["Stage2_matches_unopposed"] += sum(successful_match(r) and int(r["matchStatus"]) == 1 for r in recos)
    metrics["Stage2_matches_competing"] += sum(successful_match(r) and int(r["matchStatus"]) == 2 for r in recos)
    metrics["Stage2_no_candidate"] += sum(int(r["matchStatus"]) == 0 for r in recos)
    metrics["Stage2_lost_assignment_conflict"] += sum(int(r["matchStatus"]) == 3 for r in recos)
    match_dr_values = [float(r["matchDR"]) for r in recos if successful_match(r)]
    match_rel_values = [float(r["matchRelPt"]) for r in recos if successful_match(r)]
    metrics["Stage2_matchDR_sum"] += sum(match_dr_values)
    metrics["Stage2_matchRelPt_sum"] += sum(match_rel_values)
    metrics["Stage2_successful_matches"] += len(match_dr_values)
    for value, edges, name in ((match_dr_values, (0., .005, .01, .02, .03, math.inf), "matchDR"),
                               (match_rel_values, (0., .05, .1, .2, .5, math.inf), "matchRelPt")):
        for x in value:
            metrics["Stage2_%s_%s" % (name, _bin(x, edges, ["low", "midlow", "mid", "high", "over"]))] += 1

    # Event-level reco labels are used consistently for the neighborhoods,
    # graph estimates, and conservative supervised-pair counts.
    track_labels = {}
    for track_idx, reco in accepted_recos.items():
        gen_idx = int(reco["genChargedIdx"])
        if not successful_match(reco):
            track_labels[track_idx] = {"kind": "ambiguous_or_unmatched", "hadron": -1, "quark": -1}
            continue
        gen = gens[gen_idx]
        hidx, qidx = int(gen["hadronIdx"]), int(gen["quarkIdx"])
        if int(gen["membershipType"]) == MEMBERSHIP_HEAVY_DECAY and hidx >= 0:
            had = hadrons[hidx]
            high_conf = (int(had["associationQuality"]) == ASSOCIATION_UNIQUE
                         and int(had["fragmentationSourceCertain"]) == 1
                         and int(gen["associationQuality"]) == ASSOCIATION_UNIQUE
                         and int(gen["quarkIdx"]) == int(had["quarkIdx"])
                         and int(gen["fragmentationSourceCertain"]) == 1
                         and reco_confident(reco))
            track_labels[track_idx] = {
                "kind": "heavy_decay", "hadron": hidx, "quark": qidx,
                "high_conf": high_conf, "gen": gen_idx,
            }
        elif (int(gen["membershipType"]) == MEMBERSHIP_ADDITIONAL
              and int(gen["associationQuality"]) == ASSOCIATION_UNIQUE and qidx >= 0):
            track_labels[track_idx] = {
                "kind": "additional_uncertain", "hadron": -1, "quark": qidx,
                "high_conf": False, "gen": gen_idx,
            }
        elif int(gen["associationQuality"]) in (2, 3) or int(gen["membershipType"]) == 3:
            track_labels[track_idx] = {"kind": "ambiguous_or_unmatched", "hadron": -1, "quark": -1, "gen": gen_idx}
        elif qidx >= 0:
            track_labels[track_idx] = {"kind": "other_quark_associated", "hadron": hidx, "quark": qidx, "gen": gen_idx}
        else:
            track_labels[track_idx] = {"kind": "ambiguous_or_unmatched", "hadron": hidx, "quark": -1, "gen": gen_idx}

    hadron_rows = []
    extra_rows = []
    extra_row_index = {}
    neighborhood_rows = []
    per_hadron = {}
    prompt_d_by_q = defaultdict(list)
    high_conf_q_tracks = defaultdict(set)
    uncertain_q_tracks = defaultdict(set)
    q_system_d_categories = defaultdict(Counter)

    for hidx, had in enumerate(hadrons):
        category = origin_category(had["pdgId"], had["charmOrigin"])
        decay_gen = heavy_decay_gen_rows(gens, hidx)
        accepted_gen = [ig for ig in decay_gen if accepted(gens[ig]["pt"], gens[ig]["eta"])]
        decay_reco = {}
        for igen in decay_gen:
            reco = matched.get(igen)
            if reco is not None and reco_accept(reco):
                decay_reco[int(reco["trackIdx"])] = reco
        decay_tracks = list(decay_reco.values())
        n_tracks = len(decay_tracks)
        ncat = n_category(n_tracks)
        qidx = int(had["quarkIdx"])
        quality = int(had["associationQuality"])
        qvalid = 0 <= qidx < len(quarks)
        q_reliable = qvalid and quality == ASSOCIATION_UNIQUE
        hrow = {
            **base, "QHadron_row": hidx, "QFrag_row": qidx,
            "hadron_pdgId": int(had["pdgId"]), "origin_category": category,
            "hadron_pt": float(had["pt"]), "hadron_eta": float(had["eta"]), "hadron_phi": float(had["phi"]),
            "generator_decay_charged_descendants": len(decay_gen),
            "accepted_generator_decay_particles": len(accepted_gen),
            "accepted_matched_reco_decay_tracks": n_tracks, "N_category": ncat,
            "N_ge3_inclusive_category": "N3+" if n_tracks >= 3 else ncat,
            "quark_association_quality": quality,
            "fragmentation_source_certain": int(had["fragmentationSourceCertain"]),
            "valid_QFrag_index": int(qvalid), "unique_QFrag_association": int(q_reliable),
            "parentBHadronIdx": int(had["parentBHadronIdx"]),
            "decay_reco_matchStatus2_count": sum(int(r["matchStatus"]) == 2 for r in decay_tracks),
            "preHadPt": float(had["preHadPt"]), "preHadEta": float(had["preHadEta"]),
            "preHadPhi": float(had["preHadPhi"]),
            "preHadKinematicsValid": int(had["preHadKinematicsValid"]),
        }
        hadron_rows.append(hrow)
        per_hadron[hidx] = {"row": hrow, "decay_gen": decay_gen, "accepted_gen": accepted_gen,
                             "decay_tracks": decay_tracks, "track_ids": set(decay_reco), "category": category,
                             "qidx": qidx, "q_reliable": q_reliable, "hadron": had}
        metrics["hadron_" + category] += 1
        if n_tracks < 3:
            metrics["hadron_%s_N%d" % (category, n_tracks)] += 1
        elif n_tracks == 3:
            metrics["hadron_%s_N3" % category] += 1
        metrics["hadron_%s_N3plus" % category] += int(n_tracks >= 3)
        metrics["hadron_%s_with_valid_unique_QFrag" % category] += int(q_reliable)
        metrics["hadron_%s_with_source_certain_QFrag" % category] += int(q_reliable and int(had["fragmentationSourceCertain"]) == 1)
        metrics["hadron_%s_associationQuality_%d" % (category, quality)] += 1
        metrics["hadron_%s_source_certain" % category] += int(had["fragmentationSourceCertain"] == 1)
        if category == "D_prompt":
            metrics["prompt_D_pdg_%d" % int(had["pdgId"])] += 1
            metrics["prompt_D_pt_%s" % _bin(had["pt"], [0, 1, 2, 4, 8, math.inf], ["0-1", "1-2", "2-4", "4-8", "8+"]) ] += 1
            metrics["prompt_D_absEta_%s" % _bin(abs(had["eta"]), [0, .8, 1.6, 2.5, math.inf], ["0-.8", ".8-1.6", "1.6-2.5", "2.5+"]) ] += 1
            prompt_d_by_q[(qidx if q_reliable else -1)].append(hidx)
        if category in ("D_prompt", "D_fromB", "D_unknown") and q_reliable:
            q_system_d_categories[qidx]["%s_%s" % (category, ncat)] += 1
            if category == "D_prompt" and n_tracks >= 3:
                q_system_d_categories[qidx]["D_prompt_N3plus"] += 1
    # Additional membership rows are counted per D explicitly, then system
    # totals are deduplicated by (event,QFrag,QGenCharged) below.
    extra_by_hadron = {}
    for hidx, info in per_hadron.items():
        hrow = info["row"]
        category = info["category"]
        qidx = info["qidx"]
        use_for_extra = category == "D_prompt" and info["q_reliable"]
        gen_candidates = additional_gen_rows(gens, qidx) if use_for_extra else []
        gen_candidates = [ig for ig in gen_candidates if int(gens[ig]["hadronIdx"]) != hidx
                          and int(gens[ig]["membershipType"]) == MEMBERSHIP_ADDITIONAL]
        accepted_gen_candidates = [ig for ig in gen_candidates if accepted(gens[ig]["pt"], gens[ig]["eta"])]
        candidate_recos = {}
        for igen in gen_candidates:
            reco = matched.get(igen)
            if reco is not None:
                candidate_recos[int(reco["trackIdx"])] = (igen, reco)
                if int(gens[igen]["fragmentationSourceCertain"]) == 0:
                    uncertain_q_tracks[qidx].add(int(reco["trackIdx"]))
                extra_rows.append({
                    **base, "QHadron_row": hidx, "QFrag_row": qidx,
                    "QGenCharged_row": igen, "QRecoTrack_row": int(reco["trackIdx"]),
                    "reco_pt": float(reco["pt"]), "reco_eta": float(reco["eta"]), "reco_phi": float(reco["phi"]),
                    "reco_charge": int(reco["charge"]), "reco_isHighPurity": int(reco["isHighPurity"]),
                    "generator_pt": float(gens[igen]["pt"]), "generator_eta": float(gens[igen]["eta"]), "generator_phi": float(gens[igen]["phi"]),
                    "generator_accepted": int(accepted(gens[igen]["pt"], gens[igen]["eta"])),
                    "reco_kinematic_accepted": int(accepted(reco["pt"], reco["eta"])),
                    "reco_accepted": int(reco_accept(reco)),
                    "deltaR_to_generator_D": _dr(gens[igen], had),
                    "min_deltaR_to_D_decay_tracks": min((_dr(reco, d) for d in info["decay_tracks"]), default=-1.0),
                    "deltaR_to_reconstructed_D_direction": -1.0,
                    "association_quality": int(gens[igen]["associationQuality"]),
                    "fragmentation_source_certain": int(gens[igen]["fragmentationSourceCertain"]),
                    "reco_gen_match_status": int(reco["matchStatus"]),
                    "reco_gen_deltaR": float(reco["matchDR"]), "reco_gen_relative_pt": float(reco["matchRelPt"]),
                })
                extra_row_index[(hidx, igen)] = len(extra_rows) - 1
        extra_by_hadron[hidx] = {"gen": gen_candidates, "accepted_gen": accepted_gen_candidates,
                                 "matched": candidate_recos,
                                 "accepted_reco": {tid: v for tid, v in candidate_recos.items() if reco_accept(v[1])} }
        if category == "D_prompt":
            metrics["promptD_unique_objects"] += 1
            metrics["promptD_unique_objects_with_extra_reco"] += int(bool(extra_by_hadron[hidx]["accepted_reco"]))
            metrics["promptD_unique_accepted_additional_reco_tracks"] += len(extra_by_hadron[hidx]["accepted_reco"])
            study_categories = [hrow["N_category"]]
            if len(info["track_ids"]) >= 3 and hrow["N_category"] != "N3+":
                study_categories.append("N3+")
            for study_ncat in study_categories:
                key = "promptD_%s_" % study_ncat
                metrics[key + "objects"] += 1
                metrics[key + "objects_with_valid_unique_QFrag"] += int(info["q_reliable"])
                metrics[key + "additional_gen_candidates"] += len(gen_candidates)
                metrics[key + "accepted_additional_gen_candidates"] += len(accepted_gen_candidates)
                metrics[key + "matched_additional_reco_tracks"] += len(candidate_recos)
                metrics[key + "accepted_additional_reco_tracks"] += len(extra_by_hadron[hidx]["accepted_reco"])
                metrics[key + "objects_with_extra_reco"] += int(bool(extra_by_hadron[hidx]["accepted_reco"]))
                metrics[key + "objects_with_source_certain_extra"] += int(any(
                    int(gens[ig]["fragmentationSourceCertain"]) == 1 for ig in gen_candidates))

    # Neighborhood scan: N1 is centered on its one accepted daughter track as
    # an oracle diagnostic; N2+ uses the accepted reco-daughter momentum sum.
    accepted_track_rows = list(accepted_recos.values())
    event_eta = np.asarray([r["eta"] for r in accepted_track_rows], dtype=float)
    event_phi = np.asarray([r["phi"] for r in accepted_track_rows], dtype=float)
    event_track_ids = [int(r["trackIdx"]) for r in accepted_track_rows]
    track_pos = {tid: i for i, tid in enumerate(event_track_ids)}
    d_centers = []
    for hidx, info in per_hadron.items():
        if info["category"] != "D_prompt" or info["row"]["N_category"] not in ("N1", "N2", "N3", "N3+"):
            continue
        ncat = info["row"]["N_category"]
        dtracks = info["decay_tracks"]
        if len(dtracks) == 1:
            center_eta, center_phi = float(dtracks[0]["eta"]), float(dtracks[0]["phi"])
            center_kind = "single_decay_track_oracle"
        else:
            center_eta, center_phi = direction_from_tracks(dtracks)
            center_kind = "reconstructed_D_momentum_from_decay_tracks"
        if center_eta is not None and center_phi is not None:
            d_centers.append((hidx, center_eta, center_phi))
        reco_d_eta, reco_d_phi = direction_from_tracks(dtracks) if len(dtracks) >= 2 else (None, None)
        extra = extra_by_hadron.get(hidx, {"accepted_reco": {}})["accepted_reco"]
        decay_ids = info["track_ids"]
        for tid, (igen, reco) in extra.items():
            if reco_d_eta is not None:
                extra_rows_index = extra_row_index.get((hidx, igen))
                if extra_rows_index is not None:
                    extra_rows[extra_rows_index]["deltaR_to_reconstructed_D_direction"] = float(
                        delta_r(reco["eta"], reco["phi"], reco_d_eta, reco_d_phi))
        if center_eta is None or center_phi is None:
            continue
        distances = delta_r(event_eta, event_phi, center_eta, center_phi) if len(event_eta) else np.empty(0)
        for radius in NEIGHBORHOOD_RADII:
            nearby = (set(event_track_ids) if radius == "all" else
                      {event_track_ids[i] for i, d in enumerate(distances) if float(d) < radius})
            true_decay = nearby & decay_ids
            extra_near = nearby & set(extra)
            other = nearby - decay_ids - set(extra)
            other_classes = Counter(track_labels.get(tid, {}).get("kind", "ambiguous_or_unmatched") for tid in other)
            candidate_track_rows = [accepted_recos[tid] for tid in sorted(extra_near)]
            other_track_rows = [accepted_recos[tid] for tid in sorted(other)]
            denominator = len(extra_near) + len(other)
            n_extra_near = len(extra_near)
            neighborhood_row = {
                **base, "QHadron_row": hidx, "QFrag_row": info["qidx"], "N_category": ncat,
                "neighborhood_radius": radius, "center_definition": center_kind,
                "number_of_D_objects": 1, "objects_with_ge1_extra_candidate": int(bool(extra_near)),
                "objects_with_0_extra_candidates": int(n_extra_near == 0),
                "objects_with_1_extra_candidate": int(n_extra_near == 1),
                "objects_with_2_extra_candidates": int(n_extra_near == 2),
                "objects_with_3plus_extra_candidates": int(n_extra_near >= 3),
                "additional_candidate_tracks": len(extra_near), "unrelated_reco_tracks": len(other),
                "other_heavy_decay_tracks": other_classes["heavy_decay"],
                "other_quark_associated_tracks": other_classes["other_quark_associated"],
                "ambiguous_or_unmatched_tracks": other_classes["ambiguous_or_unmatched"],
                "candidate_to_other_track_ratio": len(extra_near) / len(other) if other else (math.inf if extra_near else 0.0),
                "candidate_fraction_nondecay": len(extra_near) / denominator if denominator else 0.0,
                "candidate_mean_pt": float(np.mean([r["pt"] for r in candidate_track_rows])) if candidate_track_rows else math.nan,
                "candidate_median_pt": float(np.median([r["pt"] for r in candidate_track_rows])) if candidate_track_rows else math.nan,
                "candidate_mean_deltaR_to_center": float(np.mean([distances[track_pos[tid]] for tid in sorted(extra_near)])) if extra_near else math.nan,
                "candidate_positive_charge_tracks": sum(int(r["charge"]) > 0 for r in candidate_track_rows),
                "candidate_negative_charge_tracks": sum(int(r["charge"]) < 0 for r in candidate_track_rows),
                "candidate_highPurity_fraction": float(np.mean([int(r["isHighPurity"]) for r in candidate_track_rows])) if candidate_track_rows else math.nan,
                "other_mean_pt": float(np.mean([r["pt"] for r in other_track_rows])) if other_track_rows else math.nan,
                "other_median_pt": float(np.median([r["pt"] for r in other_track_rows])) if other_track_rows else math.nan,
                "other_positive_charge_tracks": sum(int(r["charge"]) > 0 for r in other_track_rows),
                "other_negative_charge_tracks": sum(int(r["charge"]) < 0 for r in other_track_rows),
                "other_highPurity_fraction": float(np.mean([int(r["isHighPurity"]) for r in other_track_rows])) if other_track_rows else math.nan,
                "N_true_D_decay_tracks_in_neighborhood": len(true_decay),
                "true_D_decay_track_coverage": len(true_decay) / len(decay_ids) if decay_ids else 0.0,
                "total_neighborhood_multiplicity": len(nearby),
            }
            neighborhood_rows.append(neighborhood_row)
            if ncat == "N3":
                neighborhood_inclusive = dict(neighborhood_row)
                neighborhood_inclusive["N_category"] = "N3+"
                neighborhood_rows.append(neighborhood_inclusive)

    for radius in RADII:
        for i, (ha, eta_a, phi_a) in enumerate(d_centers):
            for hb, eta_b, phi_b in d_centers[i+1:]:
                if float(delta_r(eta_a, phi_a, eta_b, phi_b)) < 2.0 * radius:
                    metrics["prompt_D_overlapping_neighborhood_pairs_R%s" % str(radius).replace(".", "p")] += 1

    # Per-D observables use the saved pre-hadronization proxy direction and
    # momentum. Only prompt D + unique valid charm QFrag rows are eligible
    # auxiliary labels; B->D and unknown-origin values are diagnostic only.
    neighborhood04 = {int(r["QHadron_row"]): r for r in neighborhood_rows
                      if r["neighborhood_radius"] == 0.4}
    for hidx, info in per_hadron.items():
        if not info["category"].startswith("D_"):
            continue
        had, hrow = info["hadron"], info["row"]
        qidx = info["qidx"]
        qpdg = int(quarks[qidx]["pdgId"]) if 0 <= qidx < len(quarks) else 0
        prehad_valid = int(had["preHadKinematicsValid"]) == 1
        unique_assoc = info["q_reliable"]
        charm_q = qpdg != 0 and abs(qpdg) == 4
        obs, finite_inputs = _fragmentation_observables(
            {"pt": had["pt"], "eta": had["eta"], "phi": had["phi"]},
            {"pt": had["preHadPt"], "eta": had["preHadEta"], "phi": had["preHadPhi"]})
        obs_valid = bool(prehad_valid and unique_assoc and finite_inputs)
        prompt_eligible = bool(info["category"] == "D_prompt" and unique_assoc and charm_q and prehad_valid and finite_inputs)
        extra = extra_by_hadron.get(hidx, {"accepted_reco": {}})["accepted_reco"]
        near = neighborhood04.get(hidx, {})
        hrow.update({
            "QFrag_pdgId": qpdg, "unique_valid_charm_QFrag": int(unique_assoc and charm_q),
            "preHad_valid_mask": int(prehad_valid), "finite_nonzero_momentum_mask": int(finite_inputs),
            "observable_inputs_valid": int(obs_valid), "auxiliary_target_eligible": int(prompt_eligible),
            "accepted_additional_reco_candidate_count": len(extra),
            "local_reco_track_multiplicity_R0p4": near.get("total_neighborhood_multiplicity", math.nan),
            "additional_candidate_mean_deltaR_R0p4": near.get("candidate_mean_deltaR_to_center", math.nan),
            "observable_interpretation": ("prompt_D_charm_lineage_candidate_target" if info["category"] == "D_prompt"
                else "D_vs_B_parent_lineage_not_B_fragmentation" if info["category"] == "D_fromB"
                else "unknown_origin_diagnostic_only"),
        })
        for name in ("zProj", "zT", "deltaR_DQuark", "kTrel"):
            hrow[name] = obs.get(name, math.nan) if obs_valid else math.nan
            hrow[name + "_valid"] = int(obs_valid and name in obs)
            hrow[name + "_outlier_tukey"] = 0
        hrow["zProj_outside_0_1"] = int(obs_valid and (obs["zProj"] < 0.0 or obs["zProj"] > 1.0))
        hrow["zT_above_1"] = int(obs_valid and obs["zT"] > 1.0)
        hrow["prompt_auxiliary_valid"] = int(prompt_eligible)
        hrow["preHad_pt_nonpositive"] = int(not math.isfinite(float(had["preHadPt"])) or float(had["preHadPt"]) <= 0)
        metrics["aux_%s_objects" % info["category"]] += 1
        metrics["aux_%s_valid" % info["category"]] += int(obs_valid)
        metrics["aux_%s_eligible" % info["category"]] += int(prompt_eligible)

    # Deduplicated quark-system membership accounting (one track/system/event).
    for hidx, info in per_hadron.items():
        if info["q_reliable"] and int(info["row"]["fragmentation_source_certain"]) == 1:
            qidx = info["qidx"]
            for tid in info["track_ids"]:
                if track_labels.get(tid, {}).get("high_conf", False):
                    high_conf_q_tracks[qidx].add(tid)
    for hidx, info in per_hadron.items():
        if info["category"] != "D_prompt" or not info["q_reliable"]:
            continue
        for tid, (igen, reco) in extra_by_hadron.get(hidx, {"accepted_reco": {}})["accepted_reco"].items():
            uncertain_q_tracks[info["qidx"]].add(tid)

    # Conservative direct-hadron pair labels: require unopposed Stage 2
    # matches and a known QHadron decay descendant. QFrag certainty is not
    # required for a same-hadron decay label; quark-system labels below are
    # evaluated with their separate stricter criteria.
    strict_labels = {tid: lab for tid, lab in track_labels.items()
                     if lab.get("kind") == "heavy_decay" and reco_confident(reco_by_track[tid])}
    strict_ids = sorted(strict_labels)
    positive_pairs = set()
    negative_pairs = set()
    for i, tid_a in enumerate(strict_ids):
        ha = strict_labels[tid_a]["hadron"]
        for tid_b in strict_ids[i+1:]:
            hb = strict_labels[tid_b]["hadron"]
            if ha == hb:
                positive_pairs.add((tid_a, tid_b))
            else:
                negative_pairs.add((tid_a, tid_b))
    metrics["OptionA_positive_same_hadron_track_pairs"] += len(positive_pairs)
    metrics["OptionA_confident_different_hadron_negative_pairs"] += len(negative_pairs)
    metrics["OptionA_positive_D_track_pairs"] += sum(
        per_hadron[strict_labels[a]["hadron"]]["category"].startswith("D_") for a, b in positive_pairs)
    metrics["OptionA_positive_B_track_pairs"] += sum(
        per_hadron[strict_labels[a]["hadron"]]["category"] == "B" for a, b in positive_pairs)
    metrics["OptionA_uncertain_additional_reco_tracks"] += sum(lab.get("kind") == "additional_uncertain" for lab in track_labels.values())
    metrics["OptionA_confident_pairs_from_N2_objects"] += sum(
        len(info["track_ids"]) * (len(info["track_ids"])-1) // 2 for info in per_hadron.values()
        if info["category"] == "D_prompt" and info["row"]["N_category"] == "N2"
        and all(tid in strict_labels for tid in info["track_ids"]))
    metrics["OptionA_confident_pairs_from_N3plus_objects"] += sum(
        len(info["track_ids"]) * (len(info["track_ids"])-1) // 2 for info in per_hadron.values()
        if info["category"] == "D_prompt" and len(info["track_ids"]) >= 3
        and all(tid in strict_labels for tid in info["track_ids"]))
    metrics["OptionA_D_prompt_N1_objects_without_pair_label"] += sum(
        info["category"] == "D_prompt" and info["row"]["N_category"] == "N1" for info in per_hadron.values())

    # Quark-system object availability; source-uncertain additional candidates
    # are retained as a separate context population, never promoted to labels.
    q_systems = set(high_conf_q_tracks) | set(uncertain_q_tracks)
    metrics["OptionB_quark_systems_with_high_conf_decay_tracks"] += sum(bool(high_conf_q_tracks[q]) for q in q_systems)
    metrics["OptionB_quark_systems_with_uncertain_additional_tracks"] += sum(bool(uncertain_q_tracks[q]) for q in q_systems)
    metrics["OptionB_high_conf_decay_reco_tracks"] += sum(len(v) for v in high_conf_q_tracks.values())
    metrics["OptionB_uncertain_additional_reco_tracks"] += sum(len(v) for v in uncertain_q_tracks.values())
    metrics["OptionB_quark_systems_with_ge2_high_conf_tracks"] += sum(len(v) >= 2 for v in high_conf_q_tracks.values())
    metrics["OptionB_quark_systems_with_ge2_total_tracks"] += sum(
        len(high_conf_q_tracks[q] | uncertain_q_tracks[q]) >= 2 for q in q_systems)
    for q in q_systems:
        n_high = len(high_conf_q_tracks[q])
        n_all = len(high_conf_q_tracks[q] | uncertain_q_tracks[q])
        metrics["OptionB_systems_highConfTracks_%s" % ("N3plus" if n_high >= 3 else "N%d" % n_high)] += 1
        metrics["OptionB_systems_totalTracks_%s" % ("N3plus" if n_all >= 3 else "N%d" % n_all)] += 1
        for ncat in ("D_prompt_N1", "D_prompt_N2", "D_prompt_N3", "D_prompt_N3+"):
            if q_system_d_categories[q][ncat]:
                metrics["OptionB_systems_with_%s" % ncat] += 1
    metrics["shared_prompt_D_systems"] += sum(len(ds) > 1 for q, ds in prompt_d_by_q.items() if q >= 0)
    metrics["prompt_Ds_in_shared_prompt_D_systems"] += sum(len(ds) for q, ds in prompt_d_by_q.items() if q >= 0 and len(ds) > 1)
    metrics["prompt_D_objects_without_valid_unique_QFrag"] += sum(
        info["category"] == "D_prompt" and not info["q_reliable"] for info in per_hadron.values())

    # Event graph-size and labeled-pair coverage estimates. Complete pairs are
    # calculated vectorially per event; no graph is materialized or trained.
    graph_track_rows = list(accepted_recos.values())
    left, right, pair_dr = upper_triangle_distances(graph_track_rows)
    n_nodes = len(graph_track_rows)
    graph_pt = np.asarray([float(r["pt"]) for r in graph_track_rows], dtype=float)
    graph_charge = np.asarray([int(r["charge"]) for r in graph_track_rows], dtype=int)
    radius04 = pair_dr < .4
    edge_log_pt_ratio = np.abs(np.log(graph_pt[left] / graph_pt[right])) if len(left) else np.empty(0)
    edge_charge_product = graph_charge[left] * graph_charge[right] if len(left) else np.empty(0, dtype=int)
    graph_additional_ids = {tid for tid, label in track_labels.items() if label.get("kind") == "additional_uncertain"}
    graph_additional_mask = np.asarray([int(r["trackIdx"]) in graph_additional_ids for r in graph_track_rows], dtype=bool)
    candidate_incident = (graph_additional_mask[left] | graph_additional_mask[right]) if len(left) else np.zeros(0, dtype=bool)
    knn_edge_pairs = nearest_neighbor_edge_pairs(left, right, pair_dr, n_nodes, 8)
    knn_additional_incident = sum(bool(graph_additional_mask[a] or graph_additional_mask[b]) for a, b in knn_edge_pairs)
    metrics["graph_all_nodes"] += n_nodes
    metrics["graph_all_pairs"] += n_nodes * (n_nodes - 1) // 2
    metrics["graph_additional_candidate_nodes"] += len(graph_additional_ids)
    metrics["graph_known_heavy_decay_nodes"] += sum(label.get("kind") == "heavy_decay" for label in track_labels.values())
    metrics["graph_ambiguous_or_unmatched_nodes"] += sum(label.get("kind") == "ambiguous_or_unmatched" for label in track_labels.values())
    metrics["graph_radius_edges_lt_0p4"] += int(np.sum(pair_dr < 0.4))
    metrics["graph_radius_edges_lt_0p6"] += int(np.sum(pair_dr < 0.6))
    metrics["graph_radius_edges_lt_1p0"] += int(np.sum(pair_dr < 1.0))
    metrics["graph_knn8_edges"] += len(knn_edge_pairs)
    metrics["graph_nodes_event_max"] = max(metrics["graph_nodes_event_max"], n_nodes)
    metrics["graph_all_pairs_event_max"] = max(metrics["graph_all_pairs_event_max"], n_nodes * (n_nodes - 1) // 2)
    metrics["graph_radius0p4_event_max"] = max(metrics["graph_radius0p4_event_max"], int(np.sum(pair_dr < 0.4)))
    metrics["graph_radius0p6_event_max"] = max(metrics["graph_radius0p6_event_max"], int(np.sum(pair_dr < 0.6)))
    graph_event = {**base, "nodes": n_nodes, "all_pairs": n_nodes*(n_nodes-1)//2,
                   "additional_candidate_nodes": len(graph_additional_ids),
                   "known_heavy_decay_nodes": sum(label.get("kind") == "heavy_decay" for label in track_labels.values()),
                   "ambiguous_or_unmatched_nodes": sum(label.get("kind") == "ambiguous_or_unmatched" for label in track_labels.values()),
                   "edges_dr_lt_0p2": int(np.sum(pair_dr < .2)), "edges_dr_lt_0p4": int(np.sum(pair_dr < .4)),
                   "edges_dr_lt_0p6": int(np.sum(pair_dr < .6)), "edges_dr_lt_0p8": int(np.sum(pair_dr < .8)),
                   "edges_dr_lt_1p0": int(np.sum(pair_dr < 1.0)), "edges_knn8_union": len(knn_edge_pairs),
                   "additional_incident_edges_knn8_union": knn_additional_incident,
                   "additional_incident_edges_dr_lt_0p2": int(np.sum((pair_dr < .2) & candidate_incident)),
                   "additional_incident_edges_dr_lt_0p4": int(np.sum((pair_dr < .4) & candidate_incident)),
                   "additional_incident_edges_dr_lt_0p6": int(np.sum((pair_dr < .6) & candidate_incident)),
                   "additional_incident_edges_dr_lt_0p8": int(np.sum((pair_dr < .8) & candidate_incident)),
                   "additional_incident_edges_dr_lt_1p0": int(np.sum((pair_dr < 1.0) & candidate_incident)),
                   "dr0p4_mean_deltaR": float(np.mean(pair_dr[radius04])) if np.any(radius04) else math.nan,
                   "dr0p4_mean_abs_log_pt_ratio": float(np.mean(edge_log_pt_ratio[radius04])) if np.any(radius04) else math.nan,
                   "dr0p4_same_charge_edges": int(np.sum((edge_charge_product > 0) & radius04)),
                   "dr0p4_opposite_charge_edges": int(np.sum((edge_charge_product < 0) & radius04))}
    pos = {int(row["trackIdx"]): i for i, row in enumerate(graph_track_rows)}
    for radius in RADII:
        suffix = str(radius).replace(".", "p")
        coverage = sum(1 for a, b in positive_pairs if a in pos and b in pos
                       and _dr(reco_by_track[a], reco_by_track[b]) < radius)
        prompt_coverage = sum(1 for a, b in positive_pairs
                              if per_hadron[strict_labels[a]["hadron"]]["category"] == "D_prompt"
                              and a in pos and b in pos and _dr(reco_by_track[a], reco_by_track[b]) < radius)
        negative_coverage = sum(1 for a, b in negative_pairs if a in pos and b in pos
                                and _dr(reco_by_track[a], reco_by_track[b]) < radius)
        metrics["OptionA_positive_pairs_retained_dr_lt_%s" % suffix] += coverage
        metrics["OptionA_positive_prompt_D_pairs_retained_dr_lt_%s" % suffix] += prompt_coverage
        metrics["OptionA_confident_negative_pairs_retained_dr_lt_%s" % suffix] += negative_coverage
    metrics["OptionA_positive_pairs_total_for_graph_coverage"] += len(positive_pairs)
    metrics["OptionA_positive_pairs_prompt_D"] += sum(
        per_hadron[strict_labels[a]["hadron"]]["category"] == "D_prompt" for a, b in positive_pairs)
    return hadron_rows, extra_rows, neighborhood_rows, metrics, graph_event


def _metric_rows(total: Counter, graph_events: list[dict], n_events: int):
    rows = []
    def add(section, metric, value, numerator="", denominator="", note=""):
        rows.append({"section": section, "metric": metric, "value": value,
                     "numerator": numerator, "denominator": denominator, "note": note})
    add("input", "events_processed", n_events)
    for category in ("B", "D_prompt", "D_fromB", "D_unknown"):
        denominator = int(total["hadron_" + category])
        for ncat in ("N0", "N1", "N2", "N3", "N3+"):
            metric_ncat = "N3plus" if ncat == "N3+" else ncat
            count = int(total["hadron_%s_%s" % (category, metric_ncat)])
            add("truth_population", "%s_%s_objects" % (category, ncat), count, count, denominator,
                "N3 is exactly three; N3+ is three or more and overlaps exact N3 by design." if ncat in ("N3", "N3+") else "")
        add("truth_population", category + "_objects_with_unique_QFrag", total[category+"_with_valid_unique_QFrag"],
            total[category+"_with_valid_unique_QFrag"], denominator)
    for ncat in ("N1", "N2", "N3", "N3+"):
        key = "promptD_%s_" % ncat
        denom = int(total[key+"objects"])
        add("prompt_D_additional_tracks", ncat+"_objects", denom, denom, "")
        for suffix, label in (("objects_with_extra_reco", "objects_with_ge1_accepted_additional_reco_track"),
                              ("accepted_additional_gen_candidates", "accepted_additional_generator_candidates"),
                              ("matched_additional_reco_tracks", "matched_additional_reconstructed_tracks"),
                              ("accepted_additional_reco_tracks", "accepted_additional_reconstructed_tracks")):
            value = int(total[key+suffix])
            add("prompt_D_additional_tracks", ncat+"_"+label, value, value, denom if suffix == "objects_with_extra_reco" else "")
    for metric in sorted(total):
        if metric.startswith(("prompt_D_pdg_", "prompt_D_pt_", "prompt_D_absEta_", "shared_prompt_D", "prompt_D_objects_without",
                              "prompt_D_overlapping_neighborhood_pairs_", "OptionB_systems_")):
            add("prompt_D_breakdown", metric, int(total[metric]))
        elif metric.startswith("hadron_") and ("_associationQuality_" in metric or "_source_certain" in metric):
            add("truth_association_quality", metric, int(total[metric]))
    add("label_feasibility", "OptionA_positive_same_hadron_track_pairs", int(total["OptionA_positive_same_hadron_track_pairs"]))
    add("label_feasibility", "OptionA_positive_D_track_pairs", int(total["OptionA_positive_D_track_pairs"]))
    add("label_feasibility", "OptionA_positive_B_track_pairs_not_D_targets", int(total["OptionA_positive_B_track_pairs"]))
    add("label_feasibility", "OptionA_confident_different_hadron_negative_pairs", int(total["OptionA_confident_different_hadron_negative_pairs"]))
    add("label_feasibility", "OptionA_D_prompt_N1_objects_without_pair_label", int(total["OptionA_D_prompt_N1_objects_without_pair_label"]))
    add("label_feasibility", "OptionA_uncertain_additional_reco_tracks", int(total["OptionA_uncertain_additional_reco_tracks"]))
    positive = int(total["OptionA_positive_D_track_pairs"])
    negative = int(total["OptionA_confident_different_hadron_negative_pairs"])
    add("label_feasibility", "OptionA_negative_to_positive_ratio", negative/positive if positive else "undefined", negative, positive)
    add("label_feasibility", "OptionA_positive_pairs_from_N2_objects", int(total["OptionA_confident_pairs_from_N2_objects"]))
    add("label_feasibility", "OptionA_positive_pairs_from_N3plus_objects", int(total["OptionA_confident_pairs_from_N3plus_objects"]))
    add("label_feasibility", "OptionB_quark_systems_with_high_conf_decay_tracks", int(total["OptionB_quark_systems_with_high_conf_decay_tracks"]))
    add("label_feasibility", "OptionB_quark_systems_with_uncertain_additional_tracks", int(total["OptionB_quark_systems_with_uncertain_additional_tracks"]))
    add("label_feasibility", "OptionB_high_conf_decay_reco_tracks", int(total["OptionB_high_conf_decay_reco_tracks"]))
    add("label_feasibility", "OptionB_uncertain_additional_reco_tracks", int(total["OptionB_uncertain_additional_reco_tracks"]))
    add("label_feasibility", "OptionB_systems_with_ge2_high_conf_tracks", int(total["OptionB_quark_systems_with_ge2_high_conf_tracks"]))
    add("label_feasibility", "OptionB_systems_with_ge2_total_tracks_including_uncertain_context", int(total["OptionB_quark_systems_with_ge2_total_tracks"]))
    add("label_feasibility", "prompt_D_objects_in_shared_QFrag_systems", int(total["prompt_Ds_in_shared_prompt_D_systems"]))
    add("label_feasibility", "shared_prompt_D_QFrag_systems", int(total["shared_prompt_D_systems"]))
    if graph_events:
        graph_metrics = ("nodes", "additional_candidate_nodes", "known_heavy_decay_nodes", "ambiguous_or_unmatched_nodes",
                         "all_pairs", "edges_dr_lt_0p2", "edges_dr_lt_0p4", "edges_dr_lt_0p6", "edges_dr_lt_0p8",
                         "edges_dr_lt_1p0", "edges_knn8_union", "additional_incident_edges_dr_lt_0p2",
                         "additional_incident_edges_knn8_union", "additional_incident_edges_dr_lt_0p4", "additional_incident_edges_dr_lt_0p6",
                         "additional_incident_edges_dr_lt_0p8", "additional_incident_edges_dr_lt_1p0",
                         "dr0p4_same_charge_edges", "dr0p4_opposite_charge_edges")
        for key in graph_metrics:
            values = np.asarray([e[key] for e in graph_events], dtype=np.int64)
            payload_bytes = 20 if key in ("nodes", "additional_candidate_nodes", "known_heavy_decay_nodes", "ambiguous_or_unmatched_nodes") else 12
            add("graph_feasibility", "%s_mean_per_event" % key, float(np.mean(values)), note="raw payload estimate: %d bytes per node or edge; excludes overhead." % payload_bytes)
            add("graph_feasibility", "%s_median_per_event" % key, int(np.median(values)))
            add("graph_feasibility", "%s_p95_per_event" % key, int(np.percentile(values, 95)))
            add("graph_feasibility", "%s_max_per_event" % key, int(np.max(values)))
        total_nodes = sum(e["nodes"] for e in graph_events)
        add("graph_feasibility", "approx_node_plus_dr0p4_edge_raw_bytes", total_nodes*20 + sum(e["edges_dr_lt_0p4"] for e in graph_events)*12,
            note="lower-bound NumPy payload estimate (20 B/node, 12 B/undirected edge); excludes Python/metadata overhead.")
        add("graph_feasibility", "approx_all_pairs_edge_raw_bytes", sum(e["all_pairs"] for e in graph_events)*12,
            note="lower-bound edge payload estimate; highlights quadratic complete-graph cost.")
    add("label_semantics", "high_confidence_quark_membership_rule", "associationQuality==1 AND quarkIdx>=0 AND fragmentationSourceCertain==1")
    add("label_semantics", "additional_candidate_rule", "membershipType==2 AND associationQuality==1 AND matching QFrag row; preserve source certainty (usually 0)")
    add("track_features", "available", "pt, eta, phi, charge, isHighPurity, match diagnostics")
    add("track_features", "unavailable_in_QRecoTrack", "dxy, dz, impact-parameter significances, reference-point coordinates, hit/fit details")
    add("limitations", "training", "No model is trained; all pair labels are feasibility counts and strict subsets.")
    add("limitations", "reco_gen", "Successful Stage 2 assignment is kinematic truth matching, not simulated-hit-validated tracking truth.")
    matched_count = int(total["Stage2_successful_matches"])
    add("stage2_matching", "successful_matchStatus_1", int(total["Stage2_matches_unopposed"]))
    add("stage2_matching", "successful_matchStatus_2_competing", int(total["Stage2_matches_competing"]))
    add("stage2_matching", "no_candidate_status_0", int(total["Stage2_no_candidate"]))
    add("stage2_matching", "lost_assignment_status_3", int(total["Stage2_lost_assignment_conflict"]))
    add("stage2_matching", "mean_deltaR_for_successful_matches", float(total["Stage2_matchDR_sum"])/matched_count if matched_count else "undefined", total["Stage2_matchDR_sum"], matched_count)
    add("stage2_matching", "mean_relative_pt_for_successful_matches", float(total["Stage2_matchRelPt_sum"])/matched_count if matched_count else "undefined", total["Stage2_matchRelPt_sum"], matched_count)
    for key in sorted(total):
        if ((key.startswith(("Stage2_matchDR_", "Stage2_matchRelPt_")) and not key.endswith("_sum"))
                or key.startswith("QGenCharged_associationQuality_")):
            add("stage2_matching" if key.startswith("Stage2") else "truth_association_quality", key, int(total[key]))
    return rows


def _write_outputs(outdir, prefix, hadrons, extras, neighborhoods, total, graph_events, n_events):
    import pandas as pd

    os.makedirs(outdir, exist_ok=True)
    hadron_columns = list(hadrons[0]) if hadrons else [
        "run", "luminosityBlock", "event", "event_id", "source_file", "QHadron_row", "QFrag_row",
        "hadron_pdgId", "origin_category", "hadron_pt", "hadron_eta", "hadron_phi",
        "generator_decay_charged_descendants", "accepted_generator_decay_particles",
        "accepted_matched_reco_decay_tracks", "N_category", "N_ge3_inclusive_category", "quark_association_quality",
        "fragmentation_source_certain", "valid_QFrag_index", "unique_QFrag_association", "parentBHadronIdx",
        "decay_reco_matchStatus2_count",
    ]
    extra_columns = list(extras[0]) if extras else [
        "run", "luminosityBlock", "event", "event_id", "source_file", "QHadron_row", "QFrag_row",
        "QGenCharged_row", "QRecoTrack_row", "reco_pt", "reco_eta", "reco_phi", "reco_charge", "reco_isHighPurity",
        "generator_pt", "generator_eta", "generator_phi", "generator_accepted", "reco_kinematic_accepted", "reco_accepted",
        "deltaR_to_generator_D",
        "min_deltaR_to_D_decay_tracks", "deltaR_to_reconstructed_D_direction", "association_quality",
        "fragmentation_source_certain", "reco_gen_match_status", "reco_gen_deltaR", "reco_gen_relative_pt",
    ]
    pd.DataFrame(hadrons, columns=hadron_columns).to_csv(os.path.join(outdir, prefix+"_truth_hadrons.csv"), index=False)
    pd.DataFrame(extras, columns=extra_columns).to_csv(os.path.join(outdir, prefix+"_fragmentation_tracks.csv"), index=False)
    pd.DataFrame(_metric_rows(total, graph_events, n_events)).to_csv(os.path.join(outdir, prefix+"_feasibility_summary.csv"), index=False)
    auxiliary = [r for r in hadrons if str(r.get("origin_category", "")).startswith("D_")]
    # Flag sample-defined Tukey outliers by origin and observable. Keep the
    # unmodified values in CSV; zProj/zT also have explicit physical-boundary flags.
    aux_variables = ("zProj", "zT", "deltaR_DQuark", "kTrel")
    for origin in ("D_prompt", "D_fromB", "D_unknown"):
        population = [r for r in auxiliary if r["origin_category"] == origin
                      and int(r["auxiliary_target_eligible"] if origin == "D_prompt" else r["observable_inputs_valid"])]
        for var in aux_variables:
            vals = np.asarray([float(r[var]) for r in population if math.isfinite(float(r[var]))], dtype=float)
            if not len(vals):
                continue
            q25, q75 = np.percentile(vals, [25, 75])
            lo, hi = q25 - 1.5*(q75-q25), q75 + 1.5*(q75-q25)
            for row in population:
                value = float(row[var])
                row[var+"_outlier_tukey"] = int(math.isfinite(value) and (value < lo or value > hi))
    aux_columns = ["run", "luminosityBlock", "event", "event_id", "source_file", "QHadron_row", "QFrag_row",
        "hadron_pdgId", "origin_category", "N_category", "N_ge3_inclusive_category", "hadron_pt", "hadron_eta", "hadron_phi",
        "preHadPt", "preHadEta", "preHadPhi", "preHadKinematicsValid", "QFrag_pdgId",
        "quark_association_quality", "valid_QFrag_index", "unique_QFrag_association", "unique_valid_charm_QFrag",
        "preHad_valid_mask", "finite_nonzero_momentum_mask", "observable_inputs_valid", "auxiliary_target_eligible",
        "prompt_auxiliary_valid", "accepted_additional_reco_candidate_count", "local_reco_track_multiplicity_R0p4",
        "additional_candidate_mean_deltaR_R0p4", "observable_interpretation"]
    for var in aux_variables:
        aux_columns.extend((var, var+"_valid", var+"_outlier_tukey"))
    aux_columns.extend(("zProj_outside_0_1", "zT_above_1", "preHad_pt_nonpositive"))
    pd.DataFrame(auxiliary, columns=aux_columns).to_csv(os.path.join(outdir, prefix+"_fragmentation_auxiliary.csv"), index=False)
    aux_summary = _auxiliary_summary(auxiliary)
    aux_summary_columns = ["origin_category", "N_category", "observable", "summary_type", "pt_bin", "objects", "valid_objects",
        "valid_fraction", "mean", "median", "q05", "q95", "outlier_count", "outlier_fraction", "pearson", "spearman", "note"]
    pd.DataFrame(aux_summary, columns=aux_summary_columns).to_csv(os.path.join(outdir, prefix+"_fragmentation_auxiliary_summary.csv"), index=False)
    nd = pd.DataFrame(neighborhoods)
    if len(nd):
        group_cols = ["N_category", "neighborhood_radius", "center_definition"]
        agg = nd.groupby(group_cols, dropna=False).agg(
            number_of_D_objects=("number_of_D_objects", "sum"),
            objects_with_ge1_extra_candidate=("objects_with_ge1_extra_candidate", "sum"),
            objects_with_0_extra_candidates=("objects_with_0_extra_candidates", "sum"),
            objects_with_1_extra_candidate=("objects_with_1_extra_candidate", "sum"),
            objects_with_2_extra_candidates=("objects_with_2_extra_candidates", "sum"),
            objects_with_3plus_extra_candidates=("objects_with_3plus_extra_candidates", "sum"),
            additional_candidate_tracks=("additional_candidate_tracks", "sum"),
            unrelated_reco_tracks=("unrelated_reco_tracks", "sum"),
            other_heavy_decay_tracks=("other_heavy_decay_tracks", "sum"),
            other_quark_associated_tracks=("other_quark_associated_tracks", "sum"),
            ambiguous_or_unmatched_tracks=("ambiguous_or_unmatched_tracks", "sum"),
            candidate_mean_pt=("candidate_mean_pt", "mean"),
            candidate_median_pt=("candidate_median_pt", "median"),
            candidate_mean_deltaR_to_center=("candidate_mean_deltaR_to_center", "mean"),
            candidate_positive_charge_tracks=("candidate_positive_charge_tracks", "sum"),
            candidate_negative_charge_tracks=("candidate_negative_charge_tracks", "sum"),
            candidate_highPurity_fraction_mean=("candidate_highPurity_fraction", "mean"),
            other_mean_pt=("other_mean_pt", "mean"),
            other_median_pt=("other_median_pt", "median"),
            other_positive_charge_tracks=("other_positive_charge_tracks", "sum"),
            other_negative_charge_tracks=("other_negative_charge_tracks", "sum"),
            other_highPurity_fraction_mean=("other_highPurity_fraction", "mean"),
            N_true_D_decay_tracks=("N_true_D_decay_tracks_in_neighborhood", "sum"),
            true_D_decay_track_coverage_mean=("true_D_decay_track_coverage", "mean"),
            total_neighborhood_multiplicity=("total_neighborhood_multiplicity", "sum"),
        ).reset_index()
        agg["candidate_to_other_track_ratio"] = agg.apply(lambda r: r.additional_candidate_tracks/r.unrelated_reco_tracks if r.unrelated_reco_tracks else (math.inf if r.additional_candidate_tracks else 0.0), axis=1)
        agg["candidate_fraction_nondecay"] = agg.apply(lambda r: r.additional_candidate_tracks/(r.additional_candidate_tracks+r.unrelated_reco_tracks) if r.additional_candidate_tracks+r.unrelated_reco_tracks else 0.0, axis=1)
    else:
        agg = pd.DataFrame(columns=[
            "N_category", "neighborhood_radius", "center_definition", "number_of_D_objects",
            "objects_with_ge1_extra_candidate", "objects_with_0_extra_candidates", "objects_with_1_extra_candidate",
            "objects_with_2_extra_candidates", "objects_with_3plus_extra_candidates", "additional_candidate_tracks",
            "unrelated_reco_tracks", "other_heavy_decay_tracks", "other_quark_associated_tracks",
            "ambiguous_or_unmatched_tracks", "candidate_mean_pt", "candidate_median_pt",
            "candidate_mean_deltaR_to_center", "candidate_positive_charge_tracks", "candidate_negative_charge_tracks",
            "candidate_highPurity_fraction_mean", "other_mean_pt", "other_median_pt", "other_positive_charge_tracks",
            "other_negative_charge_tracks", "other_highPurity_fraction_mean", "N_true_D_decay_tracks",
            "true_D_decay_track_coverage_mean", "total_neighborhood_multiplicity",
            "candidate_to_other_track_ratio", "candidate_fraction_nondecay",
        ])
    agg.to_csv(os.path.join(outdir, prefix+"_neighborhood_scan.csv"), index=False)
    graph_df = pd.DataFrame(graph_events, columns=list(graph_events[0]) if graph_events else [
        "run", "luminosityBlock", "event", "event_id", "source_file", "nodes", "all_pairs",
        "additional_candidate_nodes", "known_heavy_decay_nodes", "ambiguous_or_unmatched_nodes",
        "edges_dr_lt_0p2", "edges_dr_lt_0p4", "edges_dr_lt_0p6", "edges_dr_lt_0p8", "edges_dr_lt_1p0", "edges_knn8_union",
        "additional_incident_edges_dr_lt_0p2", "additional_incident_edges_dr_lt_0p4", "additional_incident_edges_dr_lt_0p6",
        "additional_incident_edges_dr_lt_0p8", "additional_incident_edges_dr_lt_1p0", "additional_incident_edges_knn8_union",
        "dr0p4_mean_deltaR", "dr0p4_mean_abs_log_pt_ratio", "dr0p4_same_charge_edges", "dr0p4_opposite_charge_edges",
    ])
    graph_df.to_csv(os.path.join(outdir, prefix+"_graph_event_sizes.csv"), index=False)
    objective_rows = []
    def objective(name, metric, value, numerator="", denominator="", limitation=""):
        objective_rows.append({"objective": name, "metric": metric, "value": value,
                               "numerator": numerator, "denominator": denominator, "limitation": limitation})
    for ncat in ("N1", "N2", "N3+"):
        group = [r for r in auxiliary if r["origin_category"] == "D_prompt" and r["N_category"] in (("N3", "N3+") if ncat == "N3+" else (ncat,))]
        good = [r for r in group if int(r["auxiliary_target_eligible"])]
        objective("D_Auxiliary_fragmentation_targets", ncat+"_prompt_D_target_valid_fraction",
                  len(good)/len(group) if group else "undefined", len(good), len(group),
                  "Generator-truth regression labels for training only; inference inputs must remain reconstructed-track quantities.")
        for var in aux_variables:
            vals = [float(r[var]) for r in good if math.isfinite(float(r[var]))]
            objective("D_Auxiliary_fragmentation_targets", ncat+"_"+var+"_outlier_fraction",
                      sum(int(r[var+"_outlier_tukey"]) for r in good)/len(vals) if vals else "undefined",
                      sum(int(r[var+"_outlier_tukey"]) for r in good), len(vals),
                      "Tukey 1.5-IQR diagnostic; raw labels are preserved and never clipped.")
    n_prompt = sum(int(total["hadron_D_prompt_"+n]) for n in ("N0", "N1", "N2", "N3+"))
    positive = int(total["OptionA_positive_same_hadron_track_pairs"])
    negative = int(total["OptionA_confident_different_hadron_negative_pairs"])
    objective("A_Decay_track_pair_reconstruction", "available_events", n_events)
    objective("A_Decay_track_pair_reconstruction", "target_prompt_D_objects", n_prompt)
    objective("A_Decay_track_pair_reconstruction", "usable_positive_same_D_pairs", positive,
              limitation="Includes known D objects from prompt, B-decay, and unknown-origin categories; origin populations remain separately reported.")
    objective("A_Decay_track_pair_reconstruction", "positive_same_B_pairs_not_D_targets", int(total["OptionA_positive_B_track_pairs"]))
    objective("A_Decay_track_pair_reconstruction", "confident_different_hadron_negative_pairs", negative,
              limitation="Only successful unopposed matches to known selected-hadron decay descendants count; unresolved tracks are excluded.")
    objective("A_Decay_track_pair_reconstruction", "negative_to_positive_imbalance", negative/positive if positive else "undefined", negative, positive)
    objective("A_Decay_track_pair_reconstruction", "N2_positive_pairs", int(total["OptionA_confident_pairs_from_N2_objects"]))
    objective("A_Decay_track_pair_reconstruction", "N3plus_positive_pairs", int(total["OptionA_confident_pairs_from_N3plus_objects"]))
    objective("A_Decay_track_pair_reconstruction", "N1_without_same_hadron_pair", int(total["OptionA_D_prompt_N1_objects_without_pair_label"]),
              limitation="Requires a single-track or object-level target.")
    objective("A_Decay_track_pair_reconstruction", "prompt_D_positive_pairs_retained_with_dr_lt_0p4_graph", int(total["OptionA_positive_prompt_D_pairs_retained_dr_lt_0p4"]))
    objective("A_Decay_track_pair_reconstruction", "confident_negative_pairs_retained_with_dr_lt_0p4_graph", int(total["OptionA_confident_negative_pairs_retained_dr_lt_0p4"]))
    graph_pos, graph_neg = int(total["OptionA_positive_prompt_D_pairs_retained_dr_lt_0p4"]), int(total["OptionA_confident_negative_pairs_retained_dr_lt_0p4"])
    objective("A_Decay_track_pair_reconstruction", "dr_lt_0p4_edge_negative_to_positive_ratio", graph_neg/graph_pos if graph_pos else "undefined", graph_neg, graph_pos)
    objective("A_Decay_track_pair_reconstruction", "uncertain_additional_candidate_tracks_not_used_as_pair_labels", int(total["OptionA_uncertain_additional_reco_tracks"]))
    objective("B_Heavy_quark_system_reconstruction", "available_events", n_events)
    objective("B_Heavy_quark_system_reconstruction", "systems_with_high_confidence_decay_tracks", int(total["OptionB_quark_systems_with_high_conf_decay_tracks"]))
    objective("B_Heavy_quark_system_reconstruction", "systems_with_ge2_high_confidence_tracks", int(total["OptionB_quark_systems_with_ge2_high_conf_tracks"]))
    objective("B_Heavy_quark_system_reconstruction", "systems_with_ge2_tracks_including_uncertain_context", int(total["OptionB_quark_systems_with_ge2_total_tracks"]))
    objective("B_Heavy_quark_system_reconstruction", "high_confidence_decay_tracks", int(total["OptionB_high_conf_decay_reco_tracks"]))
    objective("B_Heavy_quark_system_reconstruction", "uncertain_additional_tracks", int(total["OptionB_uncertain_additional_reco_tracks"]),
              limitation="Unique ancestry with source certainty 0; never promoted to a supervised positive/negative label.")
    objective("B_Heavy_quark_system_reconstruction", "class_imbalance", "undefined",
              limitation="No confidently labeled negative quark-membership population is defined; unresolved/unassigned tracks are not negatives.")
    objective("B_Heavy_quark_system_reconstruction", "QFrag_systems_shared_by_multiple_prompt_D_objects", int(total["shared_prompt_D_systems"]))
    for metric in sorted(total):
        if metric.startswith("OptionB_systems_"):
            objective("B_Heavy_quark_system_reconstruction", metric, int(total[metric]))
    objective("C_Additional_track_context", "prompt_D_objects_with_accepted_extra_track", int(total["promptD_unique_objects_with_extra_reco"]),
              denominator=int(total["promptD_unique_objects"]))
    objective("C_Additional_track_context", "accepted_additional_reco_tracks", int(total["promptD_unique_accepted_additional_reco_tracks"]))
    candidate_context = sum(r["additional_candidate_tracks"] for r in neighborhoods if r["neighborhood_radius"] == .4)
    background_context = sum(r["unrelated_reco_tracks"] for r in neighborhoods if r["neighborhood_radius"] == .4)
    objective("C_Additional_track_context", "candidate_tracks_in_dr_lt_0p4_neighborhoods", candidate_context)
    objective("C_Additional_track_context", "unrelated_tracks_in_dr_lt_0p4_neighborhoods", background_context)
    objective("C_Additional_track_context", "candidate_to_unrelated_track_ratio_dr_lt_0p4", candidate_context/background_context if background_context else "undefined", candidate_context, background_context)
    objective("C_Additional_track_context", "overlapping_prompt_D_neighborhood_pairs_R0p4", int(total["prompt_D_overlapping_neighborhood_pairs_R0p4"]),
              limitation="Centers are oracle diagnostics: one daughter for N1 and summed D-decay momentum for N2+.")
    objective("C_Additional_track_context", "all_event_nodes_mean", float(np.mean([g["nodes"] for g in graph_events])) if graph_events else 0)
    objective("C_Additional_track_context", "all_event_dr_lt_0p4_edges_mean", float(np.mean([g["edges_dr_lt_0p4"] for g in graph_events])) if graph_events else 0,
              limitation="Neighborhoods are truth-centered oracle diagnostics; inference graph edges must use reco quantities only.")
    objective("graph_feasibility", "additional_candidate_nodes_mean_per_event", float(np.mean([g["additional_candidate_nodes"] for g in graph_events])) if graph_events else 0)
    objective("graph_feasibility", "additional_candidate_incident_dr_lt_0p4_edges_mean_per_event", float(np.mean([g["additional_incident_edges_dr_lt_0p4"] for g in graph_events])) if graph_events else 0)
    objective("graph_feasibility", "additional_candidate_incident_knn8_edges_mean_per_event", float(np.mean([g["additional_incident_edges_knn8_union"] for g in graph_events])) if graph_events else 0)
    dr04_edge_count = sum(g["edges_dr_lt_0p4"] for g in graph_events)
    objective("graph_feasibility", "dr_lt_0p4_edges_mean_abs_log_pt_ratio", sum(g["dr0p4_mean_abs_log_pt_ratio"] * g["edges_dr_lt_0p4"] for g in graph_events if math.isfinite(g["dr0p4_mean_abs_log_pt_ratio"])) / dr04_edge_count if dr04_edge_count else "undefined",
              limitation="Reconstructed-only pair feature diagnostic; no pT-ratio cutoff is imposed.")
    objective("graph_feasibility", "dr_lt_0p4_same_charge_edges", sum(g["dr0p4_same_charge_edges"] for g in graph_events))
    objective("graph_feasibility", "dr_lt_0p4_opposite_charge_edges", sum(g["dr0p4_opposite_charge_edges"] for g in graph_events))
    objective("graph_feasibility", "all_pairs_edge_raw_payload_bytes", sum(g["all_pairs"] for g in graph_events)*12,
              limitation="Approximate 12 bytes per undirected edge; quadratic all-pairs payload excludes framework overhead.")
    objective("graph_feasibility", "dr_lt_0p4_edge_raw_payload_bytes", sum(g["edges_dr_lt_0p4"] for g in graph_events)*12,
              limitation="Sparse geometric graph estimate using reco-only deltaR; 12 bytes per edge before metadata.")
    pd.DataFrame(objective_rows, columns=["objective", "metric", "value", "numerator", "denominator", "limitation"]).to_csv(os.path.join(outdir, prefix+"_training_feasibility.csv"), index=False)
    report = _make_report(prefix, total, n_events, graph_events, auxiliary, aux_summary)
    with open(os.path.join(outdir, prefix+"_report.md"), "w", encoding="utf-8") as output:
        output.write(report)


def _make_report(prefix, total, n_events, graph_events, auxiliary=None, aux_summary=None):
    lines = ["# HZZc Stage 3 feasibility report", "", "Input label: `%s`; events analyzed: **%d**." % (prefix, n_events), "",
             "No model was trained. N is the number of accepted reconstructed tracks matched to charged QGenCharged descendants with `membershipType=1` and the target `hadronIdx`; additional candidates never contribute to N.", "",
             "## Prompt-D reconstruction population", "",
             "| Reco decay tracks | Objects | With ≥1 accepted additional reco candidate | Fraction |", "|---|---:|---:|---:|"]
    for ncat in ("N0", "N1", "N2", "N3", "N3+"):
        metric_ncat = "N3plus" if ncat == "N3+" else ncat
        cat_count = int(total["hadron_D_prompt_"+metric_ncat])
        extra_count = int(total["promptD_%s_objects_with_extra_reco" % ncat])
        fraction = extra_count / cat_count if cat_count else float("nan")
        lines.append("| %s | %d | %d | %s |" % (ncat, cat_count, extra_count, ("%.4f" % fraction) if cat_count else "undefined"))
    lines += ["", "N3 means exactly three; N3+ means at least three and therefore includes N3. Additional candidates require a unique ancestry association to the same QFrag row, but their source-certainty flag is preserved and they are not asserted to be exclusive physical fragmentation products.", "",
              "## Label and graph feasibility", "",
              "- Strict same-D positive pairs: **%d** (including **%d** confirmed prompt-D pairs); same-B pairs are reported separately (**%d**). Confident different-hadron negative pairs: **%d**." % (total["OptionA_positive_D_track_pairs"], total["OptionA_positive_pairs_prompt_D"], total["OptionA_positive_B_track_pairs"], total["OptionA_confident_different_hadron_negative_pairs"]),
              "- Quark systems with source-certain decay tracks: **%d**; systems with source-uncertain additional matched tracks: **%d**." % (total["OptionB_quark_systems_with_high_conf_decay_tracks"], total["OptionB_quark_systems_with_uncertain_additional_tracks"]),
              "- Prompt D objects sharing a QFrag row with another prompt D: **%d** objects across **%d** systems." % (total["prompt_Ds_in_shared_prompt_D_systems"], total["shared_prompt_D_systems"]),
              "- Maximum accepted reco-track nodes in one event: **%d**; maximum ΔR<0.4 sparse edges: **%d**." % (total["graph_nodes_event_max"], total["graph_radius0p4_event_max"]),
              "- QRecoTrack provides `pt`, `eta`, `phi`, `charge`, `isHighPurity`, and matching diagnostics. It does not store dxy/dz, impact-parameter significances, reference point, or hit/fit features.",
              "- Stage 2 matches are kinematic associations, not simulated-hit-validated truth. Competing match status 2 is retained in the broad counts; strict pair labels use status 1.", "",
              "## Interpretation", "",
              "Option A has direct same-hadron labels only when both reco matches and selected-hadron decay ancestry are unambiguous; N1 cannot provide a same-hadron pair. Option B can be evaluated using source-certain decay descendants while preserving `membershipType=2` candidates with `fragmentationSourceCertain=0` as uncertain. Option C can evaluate those candidates as context. These counts alone do not establish discrimination power or select a final architecture; compare all three after running on HZZc.", "",
              "Input event identifiers use run/luminosityBlock/event. Every cross-table reference is event-local and validated before use.", ""]
    auxiliary = auxiliary or []
    aux_summary = aux_summary or []
    lines += ["", "## Auxiliary fragmentation observables", "",
              "`zProj`, `zT`, `deltaR_DQuark`, and `kTrel` use the QHadron momentum and saved pre-hadronization proxy momentum. Prompt-D training-target eligibility requires confirmed prompt origin, unique association quality, a valid event-local QFrag row with |pdgId|=4, `preHadKinematicsValid`, and finite nonzero momenta. Values are not clipped; Tukey outliers and the physical flags zProj outside [0,1] and zT>1 remain visible.", "",
              "| Prompt-D N | Objects | Valid targets | Fraction |", "|---|---:|---:|---:|"]
    for ncat in ("N1", "N2", "N3+"):
        group = [r for r in auxiliary if r["origin_category"] == "D_prompt" and r["N_category"] in (("N3", "N3+") if ncat == "N3+" else (ncat,))]
        valid = sum(int(r["auxiliary_target_eligible"]) for r in group)
        lines.append("| %s | %d | %d | %s |" % (ncat, len(group), valid, ("%.3f" % (valid/len(group))) if group else "undefined"))
    lines += ["", "For each N category, `*_fragmentation_auxiliary_summary.csv` reports distributions, Tukey rates, Pearson/Spearman correlations with D pT, accepted additional matched reco-candidate count, candidate-to-D angular separation (R<0.4), and local R<0.4 track multiplicity. It also reports additional-candidate correlations in D-pT bins. The additional-candidate count and candidate separation are truth-assisted evaluation quantities because membership comes from QGenCharged; they must not be inference inputs. Small or constant samples are marked with undefined (NaN) correlations.", "",
              "B→D rows are reported separately as D-versus-parent-bottom-lineage diagnostics; their momentum ratios are not B fragmentation observables. Unknown-origin D rows are diagnostic only. Reconstructed tracks contain pT, eta, phi, charge, purity and match diagnostics, so angular/multiplicity structure is directly observable, while all four observables also depend on the unavailable quark reference and can only be statistically predicted. The truth-assisted extra-candidate count is not an allowed feature. These are plausible masked training-only auxiliary targets, never inference inputs. Test them in Stage 4 as a limited auxiliary-loss ablation only if validity and within-pT correlations provide adequate statistics; this analysis alone does not establish predictability or benefit.", "",
              "Input event identifiers use run/luminosityBlock/event. Every cross-table reference is event-local and validated before use.", ""]
    return "\n".join(lines)


def _branch_names(prefix, fields):
    return [prefix + "_" + field for field in fields]


def run_analysis(input_files, outdir, prefix, entry_start=0, max_events=None, step_size="50 MB", require_high_purity=False):
    import awkward as ak
    import uproot

    if entry_start < 0 or (max_events is not None and max_events <= 0):
        raise ValueError("entry-start must be nonnegative and max-events positive when specified")
    if isinstance(step_size, str) and step_size.strip().isdigit():
        step_size = int(step_size.strip())
    if isinstance(step_size, int) and step_size <= 0:
        raise ValueError("step-size entry count must be positive")
    file_entries = []
    for path in input_files:
        with uproot.open(path) as root_file:
            if TREE_NAME not in root_file:
                raise ValueError("%s has no %s tree" % (path, TREE_NAME))
            file_entries.append(int(root_file[TREE_NAME].num_entries))
    total_entries = sum(file_entries)
    start = min(entry_start, total_entries)
    stop = total_entries if max_events is None else min(total_entries, start + max_events)
    required = list(EVENT_FIELDS)
    for table, fields in TABLE_FIELDS.items():
        required.extend(_branch_names(table, fields))
    hadrons, extras, neighborhoods, graph_events = [], [], [], []
    total = Counter()
    global_pos = 0
    duplicate_ids = set()
    for path, entries in zip(input_files, file_entries):
        file_begin, file_end = global_pos, global_pos + entries
        global_pos = file_end
        local_start = max(0, start - file_begin)
        local_stop = min(entries, stop - file_begin)
        if local_stop <= local_start:
            continue
        with uproot.open(path) as root_file:
            tree = root_file[TREE_NAME]
            missing = [name for name in required if name not in tree.keys()]
            if missing:
                raise ValueError("%s Events tree is missing branches: %s" % (path, ", ".join(missing)))
            for chunk in tree.iterate(expressions=required, entry_start=local_start, entry_stop=local_stop,
                                      step_size=step_size, library="ak", how=dict):
                chunk_events = len(chunk["event"])
                for iev in range(chunk_events):
                    event_id = {field: int(ak.to_list(chunk[field][iev])) for field in EVENT_FIELDS}
                    key = (event_id["run"], event_id["luminosityBlock"], event_id["event"])
                    if key in duplicate_ids:
                        raise ValueError("duplicate run/luminosityBlock/event identifier across inputs: %s" % (key,))
                    duplicate_ids.add(key)
                    tables = {}
                    for table, fields in TABLE_FIELDS.items():
                        tables[table] = {field: ak.to_list(chunk[table+"_"+field][iev]) for field in fields}
                    result = analyze_event(tables, event_id, str(path), require_high_purity=require_high_purity)
                    hr, er, nr, m, gr = result
                    hadrons.extend(hr); extras.extend(er); neighborhoods.extend(nr); graph_events.append(gr)
                    total.update(m)
    processed = max(0, stop-start)
    _write_outputs(outdir, prefix, hadrons, extras, neighborhoods, total, graph_events, processed)
    print("Processed %d of %d total Events entries (global entry range [%d, %d))." % (processed, total_entries, start, stop))
    print("Wrote %s_{truth_hadrons,fragmentation_tracks,fragmentation_auxiliary,fragmentation_auxiliary_summary,feasibility_summary,neighborhood_scan,training_feasibility,graph_event_sizes}.csv and %s_report.md under %s" % (prefix, prefix, outdir))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, nargs="+", help="one or more Stage 1/2 NanoAOD ROOT files")
    parser.add_argument("--output-dir", default="stage3_hzzc", help="directory for uncompressed CSV and report outputs")
    parser.add_argument("--prefix", default="hzzc", help="output filename prefix")
    parser.add_argument("--entry-start", type=int, default=0, help="zero-based global entry to start at across the input list")
    parser.add_argument("--max-events", type=int, default=None, help="maximum total entries across all inputs; default: through end")
    parser.add_argument("--step-size", default="50 MB", help="uproot chunk size, e.g. '50 MB' or 100 events")
    parser.add_argument("--require-high-purity", action="store_true",
                        help="optional robustness variant requiring QRecoTrack.isHighPurity==1 for reco acceptance")
    args = parser.parse_args()
    if args.entry_start < 0 or (args.max_events is not None and args.max_events <= 0):
        parser.error("entry-start must be nonnegative and max-events positive when specified")
    run_analysis(args.input, args.output_dir, args.prefix, args.entry_start, args.max_events,
                 args.step_size, args.require_high_purity)


if __name__ == "__main__":
    main()
