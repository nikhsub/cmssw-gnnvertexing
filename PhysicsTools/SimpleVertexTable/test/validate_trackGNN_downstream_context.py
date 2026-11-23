#!/usr/bin/env python3
"""Validate downstream PF/PV/SV context tables and exact GNN index mappings.

Requires uproot, awkward, and NumPy. An optional --baseline file is matched by
run/luminosityBlock/event to check that adding context preserves frozen seeder
outputs event by event.
"""
import argparse
import math
import os
import sys

try:
    import awkward as ak
    import numpy as np
    import uproot
except ImportError as exc:
    raise SystemExit("Install/import uproot, awkward, and NumPy in the validation environment: %s" % exc)


EVENT_ID = ("run", "luminosityBlock", "event")


def event_values(arrays, name, event):
    return ak.to_list(arrays[name][event])


def scalar(arrays, name, event):
    value = event_values(arrays, name, event)
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError("Expected scalar/singleton %s in event row %d, got %s" % (name, event, value))
        value = value[0]
    return value


def require(condition, errors, message):
    if not condition:
        errors.append(message)


def report_multiplicity(name, values):
    data = np.asarray(values, dtype=np.float64)
    if not len(data):
        print("%-32s no events" % name)
        return
    print("%-32s mean=%.3f median=%.3f p95=%.3f max=%d" %
          (name, float(np.mean(data)), float(np.median(data)), float(np.percentile(data, 95)), int(np.max(data))))


def validate(root_path, tree_name, baseline_path=None):
    errors = []
    required = [
        "run", "luminosityBlock", "event", "ntrack", "ntrackPair", "nSeed",
        "trackPair_trk1Idx", "trackPair_trk2Idx", "trackPair_deltaR", "trackPair_seederLogit",
        "trackPair_seederScore", "trackPair_seederValid",
        "nSeederNode", "SeederNode_trackIdx", "SeederNode_score", "SeederNode_valid",
        "nCtxChargedPF", "CtxChargedPF_sourceCollectionId", "CtxChargedPF_sourceIdx",
        "CtxChargedPF_unpackedTrackIdx", "CtxChargedPF_filteredTrackIdx", "CtxChargedPF_embeddingValid",
        "CtxChargedPF_pt", "CtxChargedPF_eta", "CtxChargedPF_phi", "CtxChargedPF_mass",
        "CtxChargedPF_energy", "CtxChargedPF_px", "CtxChargedPF_py", "CtxChargedPF_pz",
        "CtxChargedPF_charge", "CtxChargedPF_pdgId", "CtxChargedPF_hasTrackDetails",
        "CtxChargedPF_trackHighPurity", "CtxChargedPF_lostInnerHits", "CtxChargedPF_numberOfHits",
        "CtxChargedPF_numberOfPixelHits", "CtxChargedPF_numberOfStripHits", "CtxChargedPF_fromPV",
        "CtxChargedPF_pvAssociationQuality", "CtxChargedPF_associatedPVIdx", "CtxChargedPF_associatedPVValid",
        "CtxChargedPF_dxy", "CtxChargedPF_dz", "CtxChargedPF_dxyError", "CtxChargedPF_dzError",
        "CtxChargedPF_ipValid", "CtxChargedPF_ipErrorValid", "CtxChargedPF_ipReference",
        "nCtxNeutralPF", "CtxNeutralPF_sourceIdx", "CtxNeutralPF_pt", "CtxNeutralPF_eta",
        "CtxNeutralPF_phi", "CtxNeutralPF_mass", "CtxNeutralPF_energy", "CtxNeutralPF_px",
        "CtxNeutralPF_py", "CtxNeutralPF_pz", "CtxNeutralPF_pdgId", "CtxNeutralPF_charge",
        "CtxPFSummary_nPackedInput", "CtxPFSummary_nLostInput", "CtxPFSummary_nUnpackedTracks",
        "CtxPFSummary_nFilteredTracks", "CtxPFSummary_nChargedPF", "CtxPFSummary_nNeutralPF",
        "CtxPFSummary_nChargedMappedToFilteredTracks", "CtxPFSummary_nChargedUnmapped",
        "CtxPFSummary_nInvalidReferences",
        "nCtxPV", "CtxPV_sourceIdx", "CtxPV_x", "CtxPV_y", "CtxPV_z", "CtxPV_chi2",
        "CtxPV_ndof", "CtxPV_chi2Ndof", "CtxPV_isFake", "CtxPV_nTracks", "CtxPV_covXX",
        "CtxPV_covXY", "CtxPV_covXZ", "CtxPV_covYY", "CtxPV_covYZ", "CtxPV_covZZ",
        "CtxPV_positionValid", "CtxPV_covarianceValid", "CtxPV_chi2NdofValid",
        "nCtxSV", "CtxSV_sourceIdx", "CtxSV_pt", "CtxSV_eta", "CtxSV_phi", "CtxSV_mass",
        "CtxSV_energy", "CtxSV_x", "CtxSV_y", "CtxSV_z", "CtxSV_chi2", "CtxSV_ndof",
        "CtxSV_chi2Ndof", "CtxSV_nDaughters", "CtxSV_nMappedTracks", "CtxSV_covXX",
        "CtxSV_covXY", "CtxSV_covXZ", "CtxSV_covYY", "CtxSV_covYZ", "CtxSV_covZZ",
        "CtxSV_covarianceValid", "CtxSV_chi2NdofValid", "CtxSV_pvGeometryValid", "CtxSV_lxy",
        "CtxSV_lxyErr", "CtxSV_lxySig", "CtxSV_lxyValid", "CtxSV_lxyErrValid", "CtxSV_lxySigValid",
        "CtxSV_l3d", "CtxSV_l3dErr", "CtxSV_l3dSig", "CtxSV_l3dValid", "CtxSV_l3dErrValid",
        "CtxSV_l3dSigValid", "CtxSV_cosPointing2D", "CtxSV_cosPointing3D",
        "CtxSV_cosPointing2DValid", "CtxSV_cosPointing3DValid",
        "nCtxSVConstituent", "CtxSVConstituent_svIdx", "CtxSVConstituent_sourceSVIdx",
        "CtxSVConstituent_daughterIdx", "CtxSVConstituent_sourceCollectionId",
        "CtxSVConstituent_sourceIdx", "CtxSVConstituent_unpackedTrackIdx",
        "CtxSVConstituent_chargedPFIdx", "CtxSVConstituent_filteredTrackIdx",
        "CtxSVConstituent_mappingValid", "CtxSVConstituent_chargedPFMappingValid",
        "CtxSVConstituent_filteredTrackMappingValid", "CtxSVConstituent_pt", "CtxSVConstituent_eta",
        "CtxSVConstituent_phi", "CtxSVConstituent_charge",
        "nSeed", "Seed_rank", "Seed_trk1Idx", "Seed_trk2Idx", "Seed_meanLogit", "Seed_score",
        "Seed_vtxFitValid", "Seed_vtxProb",
        "nmyPV", "myPV_x", "myPV_y", "myPV_z",
        "nmySV", "mySV_nTracks", "mySV_nMappedTracks",
        "nmySVtrks", "mySVtrks_SVidx", "mySVtrks_trk_trackIdx",
        "mySVtrks_trk_sourceCollectionId", "mySVtrks_trk_sourceIdx",
    ]
    embedding_components = ["SeederEmbedding_emb%d" % i for i in range(128)]
    with uproot.open(root_path) as root:
        if tree_name not in root:
            raise SystemExit("Missing tree %r" % tree_name)
        tree = root[tree_name]
        branch_names = set(tree.keys())
        missing = [name for name in required if name not in branch_names]
        if missing:
            raise SystemExit("Missing required branches: %s" % missing)
        has_embeddings = all(name in branch_names for name in
                             ["nSeederEmbedding", "SeederEmbedding_trackIdx", "SeederEmbedding_valid"] + embedding_components)
        partial_embeddings = any(name.startswith("SeederEmbedding") or name == "nSeederEmbedding" for name in branch_names) and not has_embeddings
        if partial_embeddings:
            errors.append("SeederEmbedding table is only partially present")
        if has_embeddings:
            require_branches = ["nSeederEmbedding", "SeederEmbedding_trackIdx", "SeederEmbedding_valid"] + embedding_components
        else:
            require_branches = []
        read_names = list(dict.fromkeys(required + require_branches))
        arrays = tree.arrays(read_names, library="ak")
        event_count = tree.num_entries

        multiplicities = {name: [] for name in ("Tracks/event", "Charged PF/event", "Lost tracks/event",
                                                  "Neutral PF/event", "PVs/event", "SVs/event",
                                                  "SV constituents/event")}
        numerator_pf_mapped = denominator_pf = 0
        numerator_sv_track = denominator_sv = 0
        numerator_sv_pf = 0
        total_embedding_valid = total_embedding_rows = 0
        total_sv_cov_valid = total_sv_rows = 0
        total_pf_invalid_refs = 0

        for ev in range(event_count):
            identity = tuple(int(scalar(arrays, x, ev)) for x in EVENT_ID)
            ntracks = int(scalar(arrays, "ntrack", ev))
            nedges = int(scalar(arrays, "ntrackPair", ev))
            nseed = int(scalar(arrays, "nSeed", ev))
            ncharged = int(scalar(arrays, "nCtxChargedPF", ev))
            nneutral = int(scalar(arrays, "nCtxNeutralPF", ev))
            npv = int(scalar(arrays, "nCtxPV", ev))
            nsv = int(scalar(arrays, "nCtxSV", ev))
            nconst = int(scalar(arrays, "nCtxSVConstituent", ev))
            for name, value in (("Tracks/event", ntracks), ("Charged PF/event", ncharged),
                                ("Neutral PF/event", nneutral), ("PVs/event", npv),
                                ("SVs/event", nsv), ("SV constituents/event", nconst)):
                multiplicities[name].append(value)

            charged_src = [int(x) for x in event_values(arrays, "CtxChargedPF_sourceCollectionId", ev)]
            charged_idx = [int(x) for x in event_values(arrays, "CtxChargedPF_sourceIdx", ev)]
            charged_unpacked = [int(x) for x in event_values(arrays, "CtxChargedPF_unpackedTrackIdx", ev)]
            charged_filtered = [int(x) for x in event_values(arrays, "CtxChargedPF_filteredTrackIdx", ev)]
            charged_embedding_valid = [int(x) for x in event_values(arrays, "CtxChargedPF_embeddingValid", ev)]
            charged_charge = [int(x) for x in event_values(arrays, "CtxChargedPF_charge", ev)]
            charged_pt = event_values(arrays, "CtxChargedPF_pt", ev)
            charged_eta = event_values(arrays, "CtxChargedPF_eta", ev)
            require(len(charged_src) == ncharged, errors, "event %s: charged PF table count mismatch" % (identity,))
            charged_keys = list(zip(charged_src, charged_idx))
            require(len(charged_keys) == len(set(charged_keys)), errors,
                    "event %s: duplicate charged PF source identities" % (identity,))
            npacked = int(scalar(arrays, "CtxPFSummary_nPackedInput", ev))
            nlost = int(scalar(arrays, "CtxPFSummary_nLostInput", ev))
            multiplicities["Lost tracks/event"].append(sum(source == 1 for source in charged_src))
            require(int(scalar(arrays, "CtxPFSummary_nChargedPF", ev)) == ncharged, errors,
                    "event %s: summary nChargedPF mismatch" % (identity,))
            require(int(scalar(arrays, "CtxPFSummary_nNeutralPF", ev)) == nneutral, errors,
                    "event %s: summary nNeutralPF mismatch" % (identity,))
            require(int(scalar(arrays, "CtxPFSummary_nUnpackedTracks", ev)) >= 0, errors,
                    "event %s: negative unpacked track count" % (identity,))
            require(int(scalar(arrays, "CtxPFSummary_nFilteredTracks", ev)) == ntracks, errors,
                    "event %s: summary filtered track count mismatch" % (identity,))
            require(int(scalar(arrays, "CtxPFSummary_nChargedUnmapped", ev)) == sum(i < 0 for i in charged_filtered), errors,
                    "event %s: nChargedUnmapped mismatch" % (identity,))
            require(int(scalar(arrays, "CtxPFSummary_nChargedMappedToFilteredTracks", ev)) == sum(i >= 0 for i in charged_filtered), errors,
                    "event %s: nChargedMappedToFilteredTracks mismatch" % (identity,))
            total_pf_invalid_refs += int(scalar(arrays, "CtxPFSummary_nInvalidReferences", ev))
            for row, (source, source_index, unpacked_idx, filtered_idx, charge) in enumerate(
                    zip(charged_src, charged_idx, charged_unpacked, charged_filtered, charged_charge)):
                require(source in (0, 1), errors, "event %s charged row %d: invalid source id %d" % (identity, row, source))
                require(source_index >= 0 and source_index < (npacked if source == 0 else nlost), errors,
                        "event %s charged row %d: source index out of range" % (identity, row))
                require(charge != 0, errors, "event %s charged row %d: neutral charge in charged table" % (identity, row))
                require(unpacked_idx == -1 or 0 <= unpacked_idx < int(scalar(arrays, "CtxPFSummary_nUnpackedTracks", ev)), errors,
                        "event %s charged row %d: unpackedTrackIdx out of range" % (identity, row))
                require(filtered_idx == -1 or 0 <= filtered_idx < ntracks, errors,
                        "event %s charged row %d: filteredTrackIdx out of range" % (identity, row))
                require(filtered_idx < 0 or unpacked_idx >= 0, errors,
                        "event %s charged row %d: filtered mapping has no unpacked source" % (identity, row))
                numerator_pf_mapped += int(filtered_idx >= 0)
                denominator_pf += 1
                if int(event_values(arrays, "CtxChargedPF_ipValid", ev)[row]):
                    require(int(event_values(arrays, "CtxChargedPF_ipReference", ev)[row]) == 0, errors,
                            "event %s charged row %d: impact parameters are not marked PV0-relative" % (identity, row))
                if filtered_idx >= 0 and has_embeddings:
                    emb_valid = [int(x) for x in event_values(arrays, "SeederEmbedding_valid", ev)]
                    require(charged_embedding_valid[row] == emb_valid[filtered_idx], errors,
                            "event %s charged row %d: embedding validity disagrees with mapped table row" % (identity, row))
                elif filtered_idx < 0:
                    require(charged_embedding_valid[row] == 0, errors,
                            "event %s charged row %d: unmapped candidate claims a valid embedding" % (identity, row))
            mapped_track_sources = {}
            for row, track_idx in enumerate(charged_filtered):
                if track_idx < 0:
                    continue
                identity_source = (charged_src[row], charged_idx[row])
                if track_idx in mapped_track_sources and mapped_track_sources[track_idx] != identity_source:
                    errors.append("event %s: filtered track %d has conflicting PF sources %s and %s" %
                                  (identity, track_idx, mapped_track_sources[track_idx], identity_source))
                mapped_track_sources[track_idx] = identity_source
            neutral_srcidx = [int(x) for x in event_values(arrays, "CtxNeutralPF_sourceIdx", ev)]
            neutral_charge = [int(x) for x in event_values(arrays, "CtxNeutralPF_charge", ev)]
            require(len(neutral_srcidx) == nneutral, errors, "event %s: neutral PF table count mismatch" % (identity,))
            require(len(neutral_srcidx) == len(set(neutral_srcidx)), errors, "event %s: duplicate neutral PF source indices" % (identity,))
            for idx, charge in zip(neutral_srcidx, neutral_charge):
                require(0 <= idx < npacked, errors, "event %s: neutral packed source index out of range" % (identity,))
                require(charge == 0, errors, "event %s: charged candidate in neutral PF table" % (identity,))

            # Source indices are the table row indices by construction; rows may be
            # filtered only in PF context, while all original PV/SV rows are retained.
            pv_source = [int(x) for x in event_values(arrays, "CtxPV_sourceIdx", ev)]
            require(pv_source == list(range(npv)), errors, "event %s: CtxPV source order is not original order" % (identity,))
            sv_source = [int(x) for x in event_values(arrays, "CtxSV_sourceIdx", ev)]
            require(sv_source == list(range(nsv)), errors,
                    "event %s: CtxSV source indices are not in original source order" % (identity,))
            for prefix, count, valid_branch, cov_branches in (
                    ("CtxPV", npv, "CtxPV_covarianceValid", ("covXX", "covXY", "covXZ", "covYY", "covYZ", "covZZ")),
                    ("CtxSV", nsv, "CtxSV_covarianceValid", ("covXX", "covXY", "covXZ", "covYY", "covYZ", "covZZ"))):
                valid_rows = [int(x) for x in event_values(arrays, valid_branch, ev)]
                cov_values = [event_values(arrays, prefix + "_" + branch, ev) for branch in cov_branches]
                for row in range(count):
                    if valid_rows[row]:
                        require(all(math.isfinite(float(values[row])) for values in cov_values), errors,
                                "event %s %s row %d: covariance marked valid but has nonfinite values" % (identity, prefix, row))

            pv_count = int(scalar(arrays, "nmyPV", ev))
            require(pv_count == (1 if npv > 0 else 0), errors,
                    "event %s: standard myPV singleton count disagrees with CtxPV availability" % (identity,))
            if npv and pv_count:
                for ctx_name, my_name in (("CtxPV_x", "myPV_x"), ("CtxPV_y", "myPV_y"), ("CtxPV_z", "myPV_z")):
                    ctx_value = float(event_values(arrays, ctx_name, ev)[0])
                    my_value = float(scalar(arrays, my_name, ev))
                    require(abs(ctx_value - my_value) < 1.e-5, errors,
                            "event %s: PV0 mismatch between %s and %s" % (identity, ctx_name, my_name))

            # Preserve the established, significance-selected slimmed-SV table
            # and its source daughter identities alongside the complete CtxSV.
            n_my_sv = int(scalar(arrays, "nmySV", ev))
            n_my_sv_tracks = int(scalar(arrays, "nmySVtrks", ev))
            my_sv_ntracks = [int(x) for x in event_values(arrays, "mySV_nTracks", ev)]
            my_sv_nmapped = [int(x) for x in event_values(arrays, "mySV_nMappedTracks", ev)]
            my_trk_sv = [int(x) for x in event_values(arrays, "mySVtrks_SVidx", ev)]
            my_trk_idx = [int(x) for x in event_values(arrays, "mySVtrks_trk_trackIdx", ev)]
            my_trk_source_id = [int(x) for x in event_values(arrays, "mySVtrks_trk_sourceCollectionId", ev)]
            my_trk_source_idx = [int(x) for x in event_values(arrays, "mySVtrks_trk_sourceIdx", ev)]
            require(len(my_sv_ntracks) == n_my_sv and len(my_sv_nmapped) == n_my_sv, errors,
                    "event %s: selected mySV table row count mismatch" % (identity,))
            require(all(0 <= idx < n_my_sv for idx in my_trk_sv), errors,
                    "event %s: mySVtrks has invalid selected-SV index" % (identity,))
            require(len(my_trk_idx) == n_my_sv_tracks and len(my_trk_sv) == n_my_sv_tracks and
                    len(my_trk_source_id) == n_my_sv_tracks and len(my_trk_source_idx) == n_my_sv_tracks,
                    errors, "event %s: mySVtrks branch lengths disagree with table count" % (identity,))
            require(all(idx == -1 or 0 <= idx < ntracks for idx in my_trk_idx), errors,
                    "event %s: mySVtrks contains an invalid filtered-track reference" % (identity,))
            require(all(value == -1 or value in (0, 1) for value in my_trk_source_id), errors,
                    "event %s: mySVtrks contains invalid source collection ID" % (identity,))
            require(all(value >= -1 for value in my_trk_source_idx), errors,
                    "event %s: mySVtrks contains invalid original source index" % (identity,))

            con_sv = [int(x) for x in event_values(arrays, "CtxSVConstituent_svIdx", ev)]
            con_source_sv = [int(x) for x in event_values(arrays, "CtxSVConstituent_sourceSVIdx", ev)]
            con_daughter = [int(x) for x in event_values(arrays, "CtxSVConstituent_daughterIdx", ev)]
            con_source_id = [int(x) for x in event_values(arrays, "CtxSVConstituent_sourceCollectionId", ev)]
            con_source_idx = [int(x) for x in event_values(arrays, "CtxSVConstituent_sourceIdx", ev)]
            con_unpacked = [int(x) for x in event_values(arrays, "CtxSVConstituent_unpackedTrackIdx", ev)]
            con_charged_pf = [int(x) for x in event_values(arrays, "CtxSVConstituent_chargedPFIdx", ev)]
            con_filtered = [int(x) for x in event_values(arrays, "CtxSVConstituent_filteredTrackIdx", ev)]
            con_mapping = [int(x) for x in event_values(arrays, "CtxSVConstituent_mappingValid", ev)]
            con_pf_valid = [int(x) for x in event_values(arrays, "CtxSVConstituent_chargedPFMappingValid", ev)]
            con_track_valid = [int(x) for x in event_values(arrays, "CtxSVConstituent_filteredTrackMappingValid", ev)]
            con_charge = [int(x) for x in event_values(arrays, "CtxSVConstituent_charge", ev)]
            require(len(con_sv) == nconst, errors, "event %s: SV constituent count mismatch" % (identity,))
            daughters_per_sv = [0] * nsv
            mapped_per_sv = [0] * nsv
            for row, sv_idx in enumerate(con_sv):
                require(0 <= sv_idx < nsv, errors, "event %s constituent row %d: svIdx out of range" % (identity, row))
                if not 0 <= sv_idx < nsv:
                    continue
                daughters_per_sv[sv_idx] += 1
                require(con_source_sv[row] == sv_source[sv_idx], errors,
                        "event %s constituent row %d: sourceSVIdx does not match CtxSV" % (identity, row))
                require(con_daughter[row] >= 0, errors, "event %s constituent row %d: invalid daughter ordinal" % (identity, row))
                require(con_mapping[row] == int(con_source_id[row] in (0, 1) and con_source_idx[row] >= 0), errors,
                        "event %s constituent row %d: source identity validity disagrees with sentinels" % (identity, row))
                require(con_filtered[row] == -1 or 0 <= con_filtered[row] < ntracks, errors,
                        "event %s constituent row %d: filteredTrackIdx out of range" % (identity, row))
                require(con_unpacked[row] == -1 or 0 <= con_unpacked[row] < int(scalar(arrays, "CtxPFSummary_nUnpackedTracks", ev)), errors,
                        "event %s constituent row %d: unpackedTrackIdx out of range" % (identity, row))
                require(con_track_valid[row] == int(con_filtered[row] >= 0), errors,
                        "event %s constituent row %d: filteredTrackMappingValid mismatch" % (identity, row))
                if con_filtered[row] >= 0:
                    mapped_per_sv[sv_idx] += 1
                    numerator_sv_track += 1
                denominator_sv += 1
                if con_charged_pf[row] >= 0:
                    require(con_charged_pf[row] < ncharged, errors,
                            "event %s constituent row %d: chargedPFIdx out of range" % (identity, row))
                    if con_charged_pf[row] < ncharged:
                        pf_row = con_charged_pf[row]
                        require(con_pf_valid[row] == 1, errors,
                                "event %s constituent row %d: chargedPFIdx lacks validity flag" % (identity, row))
                        require((charged_src[pf_row], charged_idx[pf_row]) == (con_source_id[row], con_source_idx[row]), errors,
                                "event %s constituent row %d: chargedPF source identity mismatch" % (identity, row))
                        require(charged_charge[pf_row] != 0, errors,
                                "event %s constituent row %d: chargedPFIdx points to neutral object" % (identity, row))
                        require(charged_filtered[pf_row] == con_filtered[row], errors,
                                "event %s constituent row %d: PF and SV filtered-track mapping disagree" % (identity, row))
                        numerator_sv_pf += 1
                else:
                    require(con_pf_valid[row] == 0, errors,
                            "event %s constituent row %d: invalid chargedPFIdx has valid flag" % (identity, row))
            sv_n_daughters = [int(x) for x in event_values(arrays, "CtxSV_nDaughters", ev)]
            sv_n_mapped = [int(x) for x in event_values(arrays, "CtxSV_nMappedTracks", ev)]
            for sv_idx in range(nsv):
                require(daughters_per_sv[sv_idx] == sv_n_daughters[sv_idx], errors,
                        "event %s SV row %d: saved daughters %d != constituent rows %d" %
                        (identity, sv_idx, sv_n_daughters[sv_idx], daughters_per_sv[sv_idx]))
                require(mapped_per_sv[sv_idx] == sv_n_mapped[sv_idx], errors,
                        "event %s SV row %d: nMappedTracks mismatch" % (identity, sv_idx))
            total_sv_rows += nsv
            total_sv_cov_valid += sum(int(x) for x in event_values(arrays, "CtxSV_covarianceValid", ev))

            if has_embeddings:
                nemb = int(scalar(arrays, "nSeederEmbedding", ev))
                emb_idx = [int(x) for x in event_values(arrays, "SeederEmbedding_trackIdx", ev)]
                emb_valid = [int(x) for x in event_values(arrays, "SeederEmbedding_valid", ev)]
                require(nemb == ntracks and emb_idx == list(range(ntracks)), errors,
                        "event %s: SeederEmbedding rows/indices do not align with tracks" % (identity,))
                total_embedding_rows += nemb
                total_embedding_valid += sum(emb_valid)
                for row, valid in enumerate(emb_valid):
                    if not valid:
                        continue
                    components = [event_values(arrays, branch, ev)[row] for branch in embedding_components]
                    require(all(math.isfinite(float(value)) for value in components), errors,
                            "event %s embedding row %d: valid vector contains nonfinite components" % (identity, row))
            else:
                # No vectors are duplicated in the PF table. If the optional table
                # is disabled, all mapped candidates correctly advertise invalid.
                require(all(value == 0 for value in charged_embedding_valid), errors,
                        "event %s: charged PF row claims embedding valid but SeederEmbedding is absent" % (identity,))

            # Frozen seeder row counts and table sizes.
            require(int(scalar(arrays, "nSeederNode", ev)) == ntracks, errors,
                    "event %s: SeederNode row count differs from ntrack" % (identity,))
            require(len(event_values(arrays, "trackPair_trk1Idx", ev)) == nedges, errors,
                    "event %s: directed pair table count mismatch" % (identity,))
            require(len(event_values(arrays, "Seed_rank", ev)) == nseed, errors,
                    "event %s: Seed table row count mismatch" % (identity,))
            require(nseed <= 128, errors, "event %s: Seed table exceeds top-128 maximum" % (identity,))

        print("Events processed:", event_count)
        for name, vals in multiplicities.items():
            report_multiplicity(name, vals)
        unique_graph_tracks = sum(
            len(set(int(x) for x in event_values(arrays, "CtxChargedPF_filteredTrackIdx", ev) if int(x) >= 0))
            for ev in range(event_count))
        all_graph_tracks = sum(int(scalar(arrays, "ntrack", ev)) for ev in range(event_count))
        total_charged = sum(int(scalar(arrays, "nCtxChargedPF", ev)) for ev in range(event_count))
        mapped_charged = sum(int(scalar(arrays, "CtxPFSummary_nChargedMappedToFilteredTracks", ev)) for ev in range(event_count))
        print("Fraction of graph tracks mapped to PF: %.6f (%d/%d)" %
              (unique_graph_tracks / all_graph_tracks if all_graph_tracks else 0.0, unique_graph_tracks, all_graph_tracks))
        print("Fraction of charged PF mapped to graph tracks: %.6f (%d/%d)" %
              (mapped_charged / total_charged if total_charged else 0.0, mapped_charged, total_charged))
        print("Fraction of SV constituents mapped to filtered tracks: %.6f (%d/%d)" %
              (numerator_sv_track / denominator_sv if denominator_sv else 0.0, numerator_sv_track, denominator_sv))
        print("Fraction of SV constituents mapped to charged PF: %.6f (%d/%d)" %
              (numerator_sv_pf / denominator_sv if denominator_sv else 0.0, numerator_sv_pf, denominator_sv))
        print("GNN embedding validity fraction: %.6f (%d/%d)" %
              (total_embedding_valid / total_embedding_rows if total_embedding_rows else 0.0,
               total_embedding_valid, total_embedding_rows))
        print("SV covariance validity fraction: %.6f (%d/%d)" %
              (total_sv_cov_valid / total_sv_rows if total_sv_rows else 0.0, total_sv_cov_valid, total_sv_rows))
        print("PF invalid-reference counter total:", total_pf_invalid_refs)
        print("Output file size [bytes]:", os.path.getsize(root_path))

        if baseline_path:
            compare_frozen_outputs(root, arrays, baseline_path, tree_name, errors)

    if errors:
        print("Validation failed with %d issue(s):" % len(errors))
        for message in errors[:100]:
            print(" -", message)
        if len(errors) > 100:
            print(" - ... %d additional issues" % (len(errors) - 100))
        return 1
    print("Downstream context validation passed")
    return 0


def compare_frozen_outputs(current_file, current_arrays, baseline_path, tree_name, errors):
    frozen = ["ntrack", "ntrackPair", "nSeed", "trackPair_seederLogit", "SeederNode_score",
              "Seed_trk1Idx", "Seed_trk2Idx", "Seed_rank", "Seed_meanLogit", "Seed_score",
              "Seed_vtxFitValid", "Seed_vtxProb"]
    with uproot.open(baseline_path) as base:
        if tree_name not in base:
            errors.append("baseline file missing tree %s" % tree_name)
            return
        tree = base[tree_name]
        branches = set(tree.keys())
        required = list(EVENT_ID) + frozen
        missing = [name for name in required if name not in branches]
        if missing:
            errors.append("baseline file missing frozen-output branches: %s" % missing)
            return
        base_arrays = tree.arrays(required, library="ak")
        baseline_index = {}
        for row in range(tree.num_entries):
            key = tuple(int(scalar(base_arrays, name, row)) for name in EVENT_ID)
            if key in baseline_index:
                errors.append("baseline contains duplicate event identity %s" % (key,))
                return
            baseline_index[key] = row
        current_tree = current_file[tree_name]
        current_keys = set()
        for row in range(current_tree.num_entries):
            key = tuple(int(scalar(current_arrays, name, row)) for name in EVENT_ID)
            if key in current_keys:
                errors.append("current file contains duplicate event identity %s" % (key,))
                continue
            current_keys.add(key)
            if key not in baseline_index:
                errors.append("baseline lacks event identity %s" % (key,))
                continue
            base_row = baseline_index[key]
            for name in frozen:
                current_value = event_values(current_arrays, name, row)
                baseline_value = event_values(base_arrays, name, base_row)
                if isinstance(current_value, list) or isinstance(baseline_value, list):
                    if len(current_value) != len(baseline_value):
                        errors.append("event %s: frozen branch %s length changed" % (key, name))
                        continue
                    for idx, (a, b) in enumerate(zip(current_value, baseline_value)):
                        if isinstance(a, (float, np.floating)) or isinstance(b, (float, np.floating)):
                            if not math.isclose(float(a), float(b), rel_tol=1.e-6, abs_tol=1.e-6):
                                errors.append("event %s: frozen branch %s row %d changed" % (key, name, idx))
                        elif a != b:
                            errors.append("event %s: frozen branch %s row %d changed" % (key, name, idx))
                elif current_value != baseline_value:
                    errors.append("event %s: frozen scalar branch %s changed" % (key, name))
        missing_current = set(baseline_index) - current_keys
        if missing_current:
            errors.append("current file is missing %d baseline event(s), first identity %s" %
                          (len(missing_current), next(iter(missing_current))))
        print("Frozen seeder baseline comparison: event-keyed rows checked", current_tree.num_entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--baseline", help="Optional earlier seeder ROOT file for event-keyed comparison")
    args = parser.parse_args()
    return validate(args.root_file, args.tree, args.baseline)


if __name__ == "__main__":
    sys.exit(main())
