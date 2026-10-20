#!/usr/bin/env python3
"""Compare full and compact ROOT productions made from identical events.

Requires uproot 5 and NumPy. Compares all AssociationGNN reconstructed
features/labels and all compact QTrain tables, then independently rebuilds
the fragmentation observables from the full QHadron/QFrag tables.
"""
import argparse
import math

import awkward as ak
import numpy as np
import uproot


NODE = (
    "isHighPurity", "charge", "nHits", "numberOfPixelHits", "numberOfStripHits", "passSeed",
    "chi2", "energy", "eta", "phi", "pt", "px", "py", "pz", "qdotp", "qoverp",
    "dz", "dzSignificance", "ip2dSignificance", "ip2dValue", "ip3dSignificance", "ip3dValue",
)
EDGE = (
    "trk1Idx", "trk2Idx", "deltaR", "deltaPhi", "deltaEta", "dca", "dcaSig", "cptopv",
    "pvtoPCA1", "pvtoPCA2", "dotprod1", "dotprod2", "pairMom", "pairInvMass", "pt", "eta", "edgeLabel",
)
QTRAIN = {
    "QTrainTrack": (
        "trackIdx", "hadronIdx", "quarkIdx", "matchStatus", "matchDR", "matchRelPt", "nGenCandidates",
        "associationQuality", "membershipType", "fragmentationSourceCertain", "quarkWeakValid", "quarkStrictValid",
        "hadronClass", "hadronTruthValid", "isB", "isPromptD", "isDFromB", "isDUnknown",
    ),
    "QTrainPair": (
        "hadronTruthValid", "sameHadron", "sameB", "samePromptD", "sameDFromB", "quarkTruthValid", "sameQuark",
        "quarkWeakTruthValid", "sameQuarkWeak",
    ),
    "QTrainHadron": (
        "hadronIdx", "pdgId", "charmOrigin", "pt", "eta", "phi", "parentBHadronIdx", "quarkIdx",
        "associationQuality", "fragmentationSourceCertain", "preHadPt", "preHadEta", "preHadPhi",
        "preHadKinematicsValid", "nAcceptedDecayTracks", "nCategory", "zProj", "zT", "deltaR_DQuark",
        "kTrel", "observableInputsValid", "auxiliaryTargetValid", "bFragmentationTargetValid",
    ),
    "QTrainQuark": ("pdgId",),
    "QTrainLegacyHadron": ("hadronIdx", "pdgId", "isB", "isD", "fromHF"),
}
D_HADRONS = {411, 421, 431, 4122, 4232, 4132, 4332, 4412, 4422, 4432, 4444}
B_HADRONS = {521, 511, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242, 5342, 5512, 5532, 5542, 5554}


def compare_jagged(a, b, branch):
    if len(a) != len(b):
        raise ValueError("event count differs for " + branch)
    if a.ndim == 1 and b.ndim == 1:
        va, vb = ak.to_numpy(a), ak.to_numpy(b)
        if va.dtype.kind in "fc" or vb.dtype.kind in "fc":
            equal = np.array_equal(va, vb, equal_nan=True)
        else:
            equal = np.array_equal(va, vb)
        if not equal:
            raise ValueError("values differ for " + branch)
        return
    counts_a, counts_b = ak.to_numpy(ak.num(a, axis=1)), ak.to_numpy(ak.num(b, axis=1))
    if not np.array_equal(counts_a, counts_b):
        raise ValueError("per-event row counts differ for " + branch)
    va = ak.to_numpy(ak.flatten(a, axis=None))
    vb = ak.to_numpy(ak.flatten(b, axis=None))
    if va.dtype.kind in "fc" or vb.dtype.kind in "fc":
        equal = np.array_equal(va, vb, equal_nan=True)
    else:
        equal = np.array_equal(va, vb)
    if not equal:
        raise ValueError("values differ for " + branch)


def wrap_phi(value):
    return (value + math.pi) % (2.0 * math.pi) - math.pi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", required=True)
    parser.add_argument("--compact", required=True)
    parser.add_argument("--rtol", type=float, default=2e-6)
    parser.add_argument("--atol", type=float, default=2e-6)
    args = parser.parse_args()

    with uproot.open(args.full) as rf, uproot.open(args.compact) as rc:
        tf, tc = rf["Events"], rc["Events"]
        if tf.num_entries != tc.num_entries:
            raise ValueError("full and compact files have different event counts")
        common = ["run", "luminosityBlock", "event"]
        reco = (["track_" + x for x in NODE] + ["track_hadidx", "track_label"] +
                ["trackPair_" + x for x in EDGE])
        derived = [table + "_" + name for table, names in QTRAIN.items() for name in names]
        check = common + reco + derived
        missing = [(tag, [x for x in check if x not in tree.keys()])
                   for tag, tree in (("full", tf), ("compact", tc))]
        missing = [(tag, xs) for tag, xs in missing if xs]
        if missing:
            raise ValueError("missing required branches: " + repr(missing))
        af, ac = tf.arrays(check, library="ak"), tc.arrays(check, library="ak")
        for name in common + reco + derived:
            compare_jagged(af[name], ac[name], name)

        # Confirm the explicitly named legacy and QFrag-derived label sets
        # stay independently available and aligned in the full file.
        legacy = ["track_hadidx", "track_label", "trackPair_edgeLabel"]
        if any(name not in tf.keys() for name in legacy):
            raise ValueError("full file lacks the original legacy truth labels")
        qgnn_check = [
            "QGNNTrack_trackIdx", "QGNNTrack_hadronIdx", "QGNNTrack_quarkIdx", "QGNNTrack_hadronClass",
            "QGNNTrack_hadronTruthValid", "QGNNTrack_associationQuality", "QGNNTrack_membershipType",
            "QGNNTrack_fragmentationSourceCertain", "QGNNPairTruth_pairIdx", "QGNNPairTruth_sameHadron",
            "QGNNPairTruth_sameB", "QGNNPairTruth_samePromptD", "QGNNPairTruth_sameDFromB",
            "QGNNPairTruth_quarkTruthValid", "QGNNPairTruth_sameQuark", "QGNNHadron_hadronIdx",
            "QGNNHadron_nAcceptedDecayTracks", "QGNNHadron_nCategory",
        ]
        absent = [name for name in qgnn_check if name not in tf.keys()]
        if absent:
            raise ValueError("full file lacks QGNN validation branches: " + ", ".join(absent))
        qf = tf.arrays(qgnn_check, library="ak")
        for i in range(tf.num_entries):
            ntracks = len(ak.to_list(af["QTrainTrack_trackIdx"][i]))
            n_pairs = len(ak.to_list(af["QTrainPair_sameHadron"][i]))
            train_track = {name: ak.to_list(af["QTrainTrack_" + name][i]) for name in QTRAIN["QTrainTrack"]}
            train_pair = {name: ak.to_list(af["QTrainPair_" + name][i]) for name in QTRAIN["QTrainPair"]}
            train_hadron = {name: ak.to_list(af["QTrainHadron_" + name][i]) for name in QTRAIN["QTrainHadron"]}
            n_hadrons = len(train_hadron["hadronIdx"])
            n_quarks = len(ak.to_list(af["QTrainQuark_pdgId"][i]))
            if train_track["trackIdx"] != list(range(ntracks)):
                raise ValueError("QTrainTrack does not preserve reconstructed track order in event %d" % i)
            legacy_hadron = {name: ak.to_list(af["QTrainLegacyHadron_" + name][i])
                             for name in QTRAIN["QTrainLegacyHadron"]}
            legacy_idx = ak.to_list(af["track_hadidx"][i])
            legacy_label = ak.to_list(af["track_label"][i])
            if legacy_hadron["hadronIdx"] != list(range(len(legacy_hadron["hadronIdx"]))):
                raise ValueError("QTrainLegacyHadron does not preserve GVTable row indices in event %d" % i)
            for idx, label in zip(legacy_idx, legacy_label):
                if idx >= 0:
                    if idx >= len(legacy_hadron["hadronIdx"]):
                        raise ValueError("legacy hadidx is out of range in event %d" % i)
                    expected = (2 if legacy_hadron["isB"][idx] == 1 else
                                3 if legacy_hadron["fromHF"][idx] == 1 else
                                4 if legacy_hadron["isD"][idx] == 1 else 5)
                    if label != expected:
                        raise ValueError("legacy track.label does not match hadidx in event %d" % i)
            if train_hadron["hadronIdx"] != list(range(n_hadrons)):
                raise ValueError("QTrainHadron does not preserve object-row indices in event %d" % i)
            if any(x < -1 or x >= n_hadrons for x in train_track["hadronIdx"]):
                raise ValueError("QTrainTrack has invalid hadron reference in event %d" % i)
            if any(x < -1 or x >= n_quarks for x in train_track["quarkIdx"]):
                raise ValueError("QTrainTrack has invalid quark reference in event %d" % i)
            if any(x < -1 or x >= n_quarks for x in train_hadron["quarkIdx"]):
                raise ValueError("QTrainHadron has invalid quark reference in event %d" % i)
            counts = [0] * n_hadrons
            track_pts = ak.to_list(af["track_pt"][i])
            track_etas = ak.to_list(af["track_eta"][i])
            for row, ih in enumerate(train_track["hadronIdx"]):
                if (train_track["hadronTruthValid"][row] == 1 and ih >= 0 and
                        train_track["membershipType"][row] == 1 and track_pts[row] >= 0.4 and abs(track_etas[row]) <= 2.5):
                    counts[ih] += 1
            expected_cat = [0 if n == 0 else 1 if n == 1 else 2 if n == 2 else 3 for n in counts]
            if counts != train_hadron["nAcceptedDecayTracks"] or expected_cat != train_hadron["nCategory"]:
                raise ValueError("QTrainHadron N-category/count mismatch in event %d" % i)
            if any((valid == 0 and label != -1) or (valid == 1 and label not in (0, 1))
                   for valid, label in zip(train_pair["hadronTruthValid"], train_pair["sameHadron"])):
                raise ValueError("invalid sameHadron mask/sentinel in event %d" % i)
            if any((valid == 0 and label != -1) or (valid == 1 and label not in (0, 1))
                   for valid, label in zip(train_pair["quarkTruthValid"], train_pair["sameQuark"])):
                raise ValueError("invalid strict same-quark mask/sentinel in event %d" % i)
            if any((valid == 0 and label != -1) or (valid == 1 and label not in (0, 1))
                   for valid, label in zip(train_pair["quarkWeakTruthValid"], train_pair["sameQuarkWeak"])):
                raise ValueError("invalid weak same-quark mask/sentinel in event %d" % i)
            src, dst = ak.to_list(af["trackPair_trk1Idx"][i]), ak.to_list(af["trackPair_trk2Idx"][i])
            for label, ti, tj in zip(train_pair["sameQuark"], src, dst):
                if label == 1 and not (train_track["quarkStrictValid"][ti] == train_track["quarkStrictValid"][tj] == 1 and
                                       train_track["quarkIdx"][ti] == train_track["quarkIdx"][tj] and
                                       train_track["fragmentationSourceCertain"][ti] == train_track["fragmentationSourceCertain"][tj] == 1):
                    raise ValueError("strict same-quark positive violates certainty in event %d" % i)
            for label, ti, tj in zip(train_pair["sameQuarkWeak"], src, dst):
                if label == 1 and not (train_track["quarkWeakValid"][ti] == train_track["quarkWeakValid"][tj] == 1 and
                                       train_track["quarkIdx"][ti] == train_track["quarkIdx"][tj]):
                    raise ValueError("weak same-quark positive violates weak ancestry in event %d" % i)
            if ak.to_list(qf["QGNNTrack_trackIdx"][i]) != list(range(ntracks)):
                raise ValueError("QGNNTrack row order mismatch in event entry %d" % i)
            for old, new in (("QGNNTrack_hadronIdx", "QTrainTrack_hadronIdx"),
                             ("QGNNTrack_quarkIdx", "QTrainTrack_quarkIdx"),
                             ("QGNNTrack_hadronClass", "QTrainTrack_hadronClass"),
                             ("QGNNTrack_hadronTruthValid", "QTrainTrack_hadronTruthValid"),
                             ("QGNNTrack_associationQuality", "QTrainTrack_associationQuality"),
                             ("QGNNTrack_membershipType", "QTrainTrack_membershipType"),
                             ("QGNNTrack_fragmentationSourceCertain", "QTrainTrack_fragmentationSourceCertain")):
                if ak.to_list(qf[old][i]) != ak.to_list(af[new][i]):
                    raise ValueError("full/compact truth-link mismatch for %s in event %d" % (old, i))
            if ak.to_list(qf["QGNNPairTruth_pairIdx"][i]) != list(range(n_pairs)):
                raise ValueError("QGNNPairTruth row order mismatch in event entry %d" % i)
            for old, new in (("QGNNPairTruth_sameHadron", "QTrainPair_sameHadron"),
                             ("QGNNPairTruth_sameB", "QTrainPair_sameB"),
                             ("QGNNPairTruth_samePromptD", "QTrainPair_samePromptD"),
                             ("QGNNPairTruth_sameDFromB", "QTrainPair_sameDFromB"),
                             ("QGNNPairTruth_quarkTruthValid", "QTrainPair_quarkTruthValid"),
                             ("QGNNPairTruth_sameQuark", "QTrainPair_sameQuark")):
                if ak.to_list(qf[old][i]) != ak.to_list(af[new][i]):
                    raise ValueError("full/compact pair-label mismatch for %s in event %d" % (old, i))
            if ak.to_list(qf["QGNNHadron_hadronIdx"][i]) != train_hadron["hadronIdx"]:
                raise ValueError("full/compact hadron object indices differ in event %d" % i)
            if (ak.to_list(qf["QGNNHadron_nAcceptedDecayTracks"][i]) != train_hadron["nAcceptedDecayTracks"] or
                    ak.to_list(qf["QGNNHadron_nCategory"][i]) != train_hadron["nCategory"]):
                raise ValueError("full/compact hadron N counts/categories differ in event %d" % i)

        # Independently recompute all four target values using full tables.
        target_fields = ["QHadron_pdgId", "QHadron_pt", "QHadron_eta", "QHadron_phi", "QHadron_preHadPt", "QHadron_preHadEta",
                         "QHadron_preHadPhi", "QHadron_preHadKinematicsValid", "QHadron_associationQuality",
                         "QHadron_quarkIdx", "QHadron_charmOrigin", "QHadron_fragmentationSourceCertain", "QFrag_pdgId"]
        absent = [name for name in target_fields if name not in tf.keys()]
        if absent:
            raise ValueError("full file lacks auxiliary-target source branches: " + ", ".join(absent))
        full_targets = tf.arrays(target_fields, library="ak")
        compared_targets = 0
        for i in range(tf.num_entries):
            cols = {k: ak.to_list(full_targets[k][i]) for k in target_fields}
            compact_cols = {k: ak.to_list(af["QTrainHadron_" + k][i]) for k in QTRAIN["QTrainHadron"]}
            for h in range(len(cols["QHadron_pt"])):
                iq = cols["QHadron_quarkIdx"][h]
                validq = 0 <= iq < len(cols["QFrag_pdgId"])
                is_d = abs(cols["QHadron_pdgId"][h]) in D_HADRONS
                is_b = abs(cols["QHadron_pdgId"][h]) in B_HADRONS
                valid_ancestry = (validq and cols["QHadron_preHadKinematicsValid"][h] == 1 and
                                  cols["QHadron_associationQuality"][h] == 1)
                d_inputs = valid_ancestry and is_d
                b_inputs = (valid_ancestry and is_b and cols["QHadron_fragmentationSourceCertain"][h] == 1 and
                            abs(cols["QFrag_pdgId"][iq]) == 5)
                vals = [cols[f][h] for f in ("QHadron_pt", "QHadron_eta", "QHadron_phi",
                                             "QHadron_preHadPt", "QHadron_preHadEta", "QHadron_preHadPhi")]
                finite_kinematics = (all(math.isfinite(float(x)) for x in vals) and
                                     vals[0] > 0 and vals[3] > 0)
                valid_d = False
                valid_b = False
                target = [math.nan] * 4
                if finite_kinematics and (d_inputs or b_inputs):
                    pt, eta, phi, qpt, qeta, qphi = map(float, vals)
                    pd = (pt*math.cos(phi), pt*math.sin(phi), pt*math.sinh(eta))
                    pq = (qpt*math.cos(qphi), qpt*math.sin(qphi), qpt*math.sinh(qeta))
                    q2 = sum(x*x for x in pq)
                    qnorm = math.sqrt(q2) if q2 > 0 else 0.0
                    dot = sum(x*y for x, y in zip(pd, pq))
                    cross = (pd[1]*pq[2]-pd[2]*pq[1], pd[2]*pq[0]-pd[0]*pq[2], pd[0]*pq[1]-pd[1]*pq[0])
                    target = [dot/q2, pt/qpt, math.hypot(eta-qeta, wrap_phi(phi-qphi)),
                              math.sqrt(sum(x*x for x in cross))/qnorm] if qnorm > 0 else target
                    valid_d = d_inputs and all(math.isfinite(x) for x in target)
                    valid_b = b_inputs and all(math.isfinite(x) for x in target)
                valid_values = valid_d or valid_b
                if int(compact_cols["observableInputsValid"][h]) != int(valid_d):
                    raise ValueError("auxiliary validity mismatch at event %d hadron %d" % (i, h))
                for name, expected in zip(("zProj", "zT", "deltaR_DQuark", "kTrel"), target):
                    got = compact_cols[name][h]
                    if valid_values and not math.isclose(float(got), expected, rel_tol=args.rtol, abs_tol=args.atol):
                        raise ValueError("%s mismatch at event %d hadron %d" % (name, i, h))
                    if not valid_values and not math.isnan(float(got)):
                        raise ValueError("%s should be NaN for invalid inputs at event %d hadron %d" % (name, i, h))
                prompt_valid = bool(valid_d and cols["QHadron_charmOrigin"][h] == 1 and
                                    abs(cols["QFrag_pdgId"][iq]) == 4)
                if int(compact_cols["auxiliaryTargetValid"][h]) != int(prompt_valid):
                    raise ValueError("prompt auxiliary-target mask mismatch at event %d hadron %d" % (i, h))
                if int(compact_cols["bFragmentationTargetValid"][h]) != int(valid_b):
                    raise ValueError("B fragmentation-target mask mismatch at event %d hadron %d" % (i, h))
                compared_targets += 1

        print("Full/compact validation passed: %d identical events, %d hadron target rows; all node/edge inputs, legacy labels, QTrain targets and strict QGNN labels agree." %
              (tf.num_entries, compared_targets))


if __name__ == "__main__":
    main()
