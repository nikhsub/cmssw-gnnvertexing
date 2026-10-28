#!/usr/bin/env python3
"""Structural checks for generator decay geometry and truth-associated pair fits."""
import argparse
import math
import sys

import awkward as ak
import numpy as np
import uproot


def branches(tree):
    return set(tree.keys())


def find(tree, suffix, required=True, prefix=None):
    target = (prefix + "_" if prefix else "") + suffix
    hits = [name for name in tree.keys() if name == target or name.endswith("_" + target)]
    if len(hits) != 1:
        if not required and not hits:
            return None
        raise ValueError("Expected one branch ending in %r; found %s" % (suffix, hits))
    return hits[0]


def rows(tree, names, prefix):
    found = {n: find(tree, n, prefix=prefix) for n in names}
    arrays = tree.arrays(list(found.values()), library="ak")
    return {n: ak.to_list(arrays[b]) for n, b in found.items()}


def singleton_rows(tree, names, prefix):
    """Normalize FlatTable(singleton=True) scalar branches to one-item rows."""
    values = rows(tree, names, prefix)
    return {name: [row if isinstance(row, list) else [row] for row in events]
            for name, events in values.items()}


def require(condition, message, errors):
    if not condition:
        errors.append(message)


def validate(path, baseline=None):
    errors = []
    with uproot.open(path) as f:
        tree = f["Events"]
        required_tables = ["QTrainHadron", "QTruthMatchedTrack", "QTruthMatchedHadron", "QTruthMatchedPairFit"]
        for table in required_tables:
            require(any(k.startswith(table + "_") for k in tree.keys()), "missing %s table" % table, errors)
        had = rows(tree, ["hadronIdx", "pdgId", "pt", "eta", "phi", "p", "px", "py", "pz", "energy", "status", "isPrompt", "isHardProcess", "statusFlagLastCopy", "prodX", "prodY", "prodZ", "prodVertexValid",
                          "decayX", "decayY", "decayZ", "decayVertexValid", "decayVertexResolved", "decayVertexSource",
                          ], "QTrainHadron")
        # Derived flight-length columns are optional in the minimal study
        # configuration; raw production/decay coordinates are authoritative.
        derived_hadron_branches = {}
        for name in ("trueLxy", "trueL3D", "trueDecayLengthValid", "trueCt", "trueCtValid"):
            branch = find(tree, name, required=False, prefix="QTrainHadron")
            if branch:
                derived_hadron_branches[name] = ak.to_list(tree[branch].array(library="ak"))
        had.update(derived_hadron_branches)
        tracks = rows(tree, ["trackIdx", "genChargedIdx", "hadronIdx", "matchStatus", "daughterRelation", "genProdX", "genProdY", "genProdZ", "genProdVertexValid"], "QTruthMatchedTrack")
        gen_count_branch = find(tree, "packedGenIdx", prefix="QGenCharged")
        n_gen = ak.to_list(tree[gen_count_branch].array(library="ak"))
        hgroups = rows(tree, ["hadronIdx", "qTrainHadronIdx", "nMatchedRecoTracks", "nPossibleTrackPairs"], "QTruthMatchedHadron")
        pairs = rows(tree, ["trackIdx1", "trackIdx2", "hadronIdx", "qTrainHadronIdx", "genChargedIdx1", "genChargedIdx2", "matchStatus1", "matchStatus2",
                            "daughterRelation1", "daughterRelation2", "commonImmediateDecayVertex", "passedOrdinaryPairSelection", "trackPairIdx",
                            "fitReusedFromTrackPair", "fitReusedFromEarlierTruthPair", "vtxFitValid", "vtxFitStatus", "vtxCovarianceValid",
                            "vtxX", "vtxY", "vtxZ", "vtxChi2", "vtxProb", "vtxCovXX", "vtxCovXY", "vtxCovXZ", "vtxCovYY", "vtxCovYZ", "vtxCovZZ",
                            "vtxLxyValid", "vtxLxy", "vtxL3DValid", "vtxL3D", "vtxPointing2DValid", "vtxPointing3DValid"], "QTruthMatchedPairFit")
        summary = singleton_rows(tree, ["run", "luminosityBlock", "event", "recoPVValid", "recoPVX", "recoPVY", "recoPVZ", "recoPVCovXX", "recoPVCovXY", "recoPVCovXZ", "recoPVCovYY", "recoPVCovYZ", "recoPVCovZZ"], "QTruthPairFitSummary")
        original_hadrons = rows(tree, ["hadronIdx", "prunedGenIdx"], "QHadron")
        tp_a = find(tree, "trk1Idx", prefix="trackPair", required=False)
        tp_b = find(tree, "trk2Idx", prefix="trackPair", required=False)
        tp_fit = {}
        required_tp_fit = ("vtxX", "vtxY", "vtxZ", "vtxChi2", "vtxProb", "vtxCovXX", "vtxCovYY", "vtxCovZZ")
        for name in required_tp_fit:
            branch = find(tree, name, prefix="trackPair", required=False)
            if branch:
                tp_fit[name] = ak.to_list(tree[branch].array(library="ak"))
        require(set(tp_fit) == set(required_tp_fit),
                "trackPair is missing one or more expected two-track fit branches", errors)
        tp_endpoints = None
        if tp_a and tp_b:
            aa, bb = ak.to_list(tree[tp_a].array(library="ak")), ak.to_list(tree[tp_b].array(library="ak"))
            tp_endpoints = [list(zip(a, b)) for a, b in zip(aa, bb)]
        nevents = tree.num_entries
        def same_events(t):
            return len(next(iter(t.values()))) == nevents
        require(same_events(had) and same_events(tracks) and same_events(hgroups) and same_events(pairs) and same_events(summary) and same_events(original_hadrons), "table event counts differ", errors)
        counts = {"hadrons": 0, "valid_decay_vertices": 0, "truth_tracks": 0, "matched_tracks": 0, "truth_pairs": 0,
                  "selected_pairs": 0, "reused_selected_fits": 0, "fit_valid": 0, "covariance_valid": 0}
        for ie in range(nevents):
            nh = len(had["hadronIdx"][ie]); nt = len(tracks["trackIdx"][ie]); npair = len(pairs["trackIdx1"][ie])
            counts["hadrons"] += nh; counts["valid_decay_vertices"] += sum(had["decayVertexValid"][ie])
            counts["truth_tracks"] += nt; counts["truth_pairs"] += npair
            track_ids = set(tracks["trackIdx"][ie]); had_ids = set(hgroups["hadronIdx"][ie])
            track_pos = {tid: k for k, tid in enumerate(tracks["trackIdx"][ie])}
            require(len(summary["event"][ie]) == 1 and len(summary["run"][ie]) == 1, "event summary is not singleton in event %d" % ie, errors)
            if summary["recoPVValid"][ie][0]:
                require(all(math.isfinite(summary[k][ie][0]) for k in ("recoPVX","recoPVY","recoPVZ","recoPVCovXX","recoPVCovXY","recoPVCovXZ","recoPVCovYY","recoPVCovYZ","recoPVCovZZ")), "valid PV has nonfinite coordinates/covariance", errors)
            original_map = {idx: pruned for idx, pruned in zip(original_hadrons["hadronIdx"][ie], original_hadrons["prunedGenIdx"][ie])}
            require(len(track_ids) == len(tracks["trackIdx"][ie]), "duplicate shared track index in event %d" % ie, errors)
            for igen in tracks["genChargedIdx"][ie]:
                require(igen == -1 or 0 <= igen < len(n_gen[ie]), "QTruthMatchedTrack has out-of-range QGenCharged reference in event %d" % ie, errors)
            require(len(had_ids) == len(hgroups["hadronIdx"][ie]), "duplicate QTruthMatchedHadron row in event %d" % ie, errors)
            for j, ih in enumerate(hgroups["hadronIdx"][ie]):
                valid_hadron = 0 <= ih < nh
                require(valid_hadron and hgroups["qTrainHadronIdx"][ie][j] == ih, "invalid hadron/QTrainHadron map in event %d" % ie, errors)
                if not valid_hadron:
                    continue
                require(had["hadronIdx"][ie][ih] == ih and ih in original_map, "QHadron/QTrainHadron row mapping mismatch in event %d" % ie, errors)
                members = [k for k, hidx in enumerate(tracks["hadronIdx"][ie]) if hidx == ih]
                counts["matched_tracks"] += len(members)
                require(hgroups["nMatchedRecoTracks"][ie][j] == len(members), "hadron track count mismatch in event %d" % ie, errors)
            seen = set()
            for j in range(npair):
                a, b, ih = pairs["trackIdx1"][ie][j], pairs["trackIdx2"][ie][j], pairs["hadronIdx"][ie][j]
                key = (ih, min(a, b), max(a, b))
                require(key not in seen, "duplicate hadron-specific track pair %r event %d" % (key, ie), errors); seen.add(key)
                valid_hadron = 0 <= ih < nh
                valid_endpoints = a in track_ids and b in track_ids and a != b
                require(valid_endpoints, "pair endpoint out of range in event %d" % ie, errors)
                require(valid_hadron and ih in had_ids and pairs["qTrainHadronIdx"][ie][j] == ih,
                        "pair has invalid hadron row in event %d" % ie, errors)
                if not valid_hadron or not valid_endpoints:
                    continue
                ta, tb = track_pos.get(a), track_pos.get(b)
                if ta is not None and tb is not None:
                    require(tracks["hadronIdx"][ie][ta] == ih and tracks["hadronIdx"][ie][tb] == ih,
                            "pair endpoint hadron association mismatch in event %d" % ie, errors)
                    require(pairs["genChargedIdx1"][ie][j] == tracks["genChargedIdx"][ie][ta] and
                            pairs["genChargedIdx2"][ie][j] == tracks["genChargedIdx"][ie][tb],
                            "pair QGenCharged references disagree with track rows in event %d" % ie, errors)
                    separation = math.sqrt(sum((tracks[axis][ie][ta]-tracks[axis][ie][tb])**2 for axis in ("genProdX","genProdY","genProdZ")))
                    distance_a = math.sqrt(sum((tracks[axis][ie][ta]-had[axis.replace("genProd", "decay")][ie][ih])**2 for axis in ("genProdX","genProdY","genProdZ")))
                    distance_b = math.sqrt(sum((tracks[axis][ie][tb]-had[axis.replace("genProd", "decay")][ie][ih])**2 for axis in ("genProdX","genProdY","genProdZ")))
                    expected_direct = int(pairs["daughterRelation1"][ie][j] == 1 and pairs["daughterRelation2"][ie][j] == 1 and
                                          tracks["genProdVertexValid"][ie][ta] and tracks["genProdVertexValid"][ie][tb] and had["decayVertexValid"][ie][ih] and
                                          separation <= 1.e-4 and distance_a <= 1.e-4 and distance_b <= 1.e-4)
                    require(pairs["commonImmediateDecayVertex"][ie][j] == expected_direct,
                            "common direct-decay vertex flag inconsistent in event %d" % ie, errors)
                if pairs["passedOrdinaryPairSelection"][ie][j]:
                    counts["selected_pairs"] += 1
                    require(pairs["trackPairIdx"][ie][j] >= 0, "selected pair lacks trackPairIdx in event %d" % ie, errors)
                    row = pairs["trackPairIdx"][ie][j]
                    if tp_endpoints is not None and 0 <= row < len(tp_endpoints[ie]):
                        x, y = tp_endpoints[ie][row]
                        require((x, y) == (a, b) or (x, y) == (b, a), "selected pair endpoint identity mismatch in event %d" % ie, errors)
                if pairs["fitReusedFromTrackPair"][ie][j]:
                    counts["reused_selected_fits"] += 1
                    row = pairs["trackPairIdx"][ie][j]
                    require(pairs["passedOrdinaryPairSelection"][ie][j] and row >= 0,
                            "fit reuse flag without selected-pair reference in event %d" % ie, errors)
                    require(bool(tp_fit), "fit reuse flag set but trackPair fit branches are absent", errors)
                    valid_pair_row = tp_endpoints is not None and 0 <= row < len(tp_endpoints[ie])
                    valid_fit_row = valid_pair_row and all(row < len(tp_fit[name][ie]) for name in tp_fit)
                    if not valid_fit_row:
                        require(False, "fit reuse trackPairIdx out of range in event %d" % ie, errors)
                    else:
                        for name in tp_fit:
                            require(abs(pairs["vtx" + name[3:]][ie][j] - tp_fit[name][ie][row]) < 1.e-5,
                                    "reused %s differs from trackPair in event %d" % (name, ie), errors)
                if pairs["vtxFitValid"][ie][j]:
                    counts["fit_valid"] += 1
                    require(all(math.isfinite(pairs[k][ie][j]) for k in ("vtxX", "vtxY", "vtxZ")), "valid fit has nonfinite vertex in event %d" % ie, errors)
                if pairs["vtxCovarianceValid"][ie][j]:
                    counts["covariance_valid"] += 1
                    xx,xy,xz,yy,yz,zz=[pairs[k][ie][j] for k in ("vtxCovXX","vtxCovXY","vtxCovXZ","vtxCovYY","vtxCovYZ","vtxCovZZ")]
                    scale=max(abs(xx),abs(yy),abs(zz),1.e-30); tol2=1.e-10*scale*scale; tol3=tol2*scale
                    minors=(xx*yy-xy*xy,xx*zz-xz*xz,yy*zz-yz*yz)
                    det=xx*yy*zz+2*xy*xz*yz-xx*yz*yz-yy*xz*xz-zz*xy*xy
                    require(min(minors)>=-tol2 and det>=-tol3, "covariance-valid row fails PSD tolerance at event %d pair %d" % (ie,j), errors)
            # All valid positions and derived lengths must be finite; unavailable values keep their masks.
            for ih in range(nh):
                if had["decayVertexValid"][ie][ih]:
                    require(had["decayVertexResolved"][ie][ih] == 1 and had["decayVertexSource"][ie][ih] == 1,
                            "decay vertex provenance mismatch in event %d" % ie, errors)
                    require(all(math.isfinite(had[k][ie][ih]) for k in ("decayX","decayY","decayZ")), "valid decay vertex nonfinite", errors)
                if "trueDecayLengthValid" in had and had["trueDecayLengthValid"][ie][ih]:
                    require(math.isfinite(had["trueLxy"][ie][ih]) and had["trueLxy"][ie][ih] >= 0 and
                            math.isfinite(had["trueL3D"][ie][ih]) and had["trueL3D"][ie][ih] >= 0, "invalid true decay length", errors)
                    lx=had["decayX"][ie][ih]-had["prodX"][ie][ih];ly=had["decayY"][ie][ih]-had["prodY"][ie][ih];lz=had["decayZ"][ie][ih]-had["prodZ"][ie][ih]
                    require(abs(had["trueLxy"][ie][ih]-math.hypot(lx,ly))<1.e-4 and abs(had["trueL3D"][ie][ih]-math.sqrt(lx*lx+ly*ly+lz*lz))<1.e-4,
                            "true decay lengths disagree with saved coordinates in event %d"%ie,errors)
                if "trueCtValid" in had and had["trueCtValid"][ie][ih]:
                    require(math.isfinite(had["trueCt"][ie][ih]), "valid trueCt is nonfinite", errors)
            for j, ih in enumerate(hgroups["hadronIdx"][ie]):
                expected = hgroups["nPossibleTrackPairs"][ie][j]
                actual = sum(x == ih for x in pairs["hadronIdx"][ie])
                require(actual == expected, "truth pair enumeration incomplete for hadron %d event %d" % (ih, ie), errors)
        if baseline:
            with uproot.open(baseline) as bf:
                bt=bf["Events"]
                for n in ("trackPair_trk1Idx","trackPair_trk2Idx"):
                    a=find(tree,n,False); b=find(bt,n,False)
                    require(a is not None and b is not None, "baseline comparison missing %s" % n, errors)
                    if a and b:
                        av,bv=tree[a].array(library="ak"),bt[b].array(library="ak")
                        require(len(av)==len(bv) and ak.to_list(av)==ak.to_list(bv), "ordinary trackPair endpoints changed vs baseline", errors)
                base_cols=[k for k in bt.keys() if k.startswith("trackPair_") and not k.startswith("trackPair_vtx")]
                for b in base_cols:
                    a=find(tree,b,False)
                    require(a is not None, "original pair branch missing from extended output: %s"%b, errors)
                    if a:
                        require(ak.to_list(tree[a].array(library="ak"))==ak.to_list(bt[b].array(library="ak")), "original pair feature differs: %s"%b, errors)
        print("Events: %d" % nevents)
        for k,v in counts.items(): print("%s: %d" % (k,v))
    if errors:
        print("Validation failures (%d):" % len(errors), file=sys.stderr)
        for e in errors[:50]: print(" - " + e, file=sys.stderr)
        return 1
    print("Validation passed")
    return 0


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("root");p.add_argument("--baseline",help="previous output for trackPair invariance comparison")
    a=p.parse_args();return validate(a.root,a.baseline)


if __name__ == "__main__":
    raise SystemExit(main())
