#!/usr/bin/env python3
"""Validate AssociationGNN and QFrag table alignment in a produced ROOT file.

Requires uproot 5 and awkward. See README_association_gnn_root.md for usage.
"""
import argparse

import awkward as ak
import uproot


NODE = (
    "isHighPurity", "charge", "nHits", "numberOfPixelHits", "numberOfStripHits", "passSeed",
    "chi2", "energy", "eta", "phi", "pt", "px", "py", "pz", "qdotp", "qoverp",
    "dz", "dzSignificance", "ip2dSignificance", "ip2dValue", "ip3dSignificance", "ip3dValue",
)
EDGE = (
    "deltaR", "deltaPhi", "deltaEta", "dca", "dcaSig", "cptopv", "pvtoPCA1", "pvtoPCA2",
    "dotprod1", "dotprod2", "pairMom", "pairInvMass", "pt", "eta",
)


def values(arrays, name, event):
    if name not in arrays.fields:
        raise ValueError("missing required branch " + name)
    return ak.to_list(arrays[name][event])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    args = parser.parse_args()
    node_fields = ["track_" + name for name in NODE]
    edge_fields = ["trackPair_" + name for name in EDGE]
    required = (
        ["track_" + name for name in NODE] + ["track_hadidx", "track_label"] + edge_fields +
        ["trackPair_trk1Idx", "trackPair_trk2Idx", "trackPair_edgeLabel"] +
        ["QRecoTrack_trackIdx", "QRecoTrack_genChargedIdx", "QRecoTrack_pt", "QRecoTrack_eta"] +
        ["QGenCharged_membershipType", "QGenCharged_associationQuality", "QGenCharged_fragmentationSourceCertain",
         "QGenCharged_hadronIdx", "QGenCharged_quarkIdx"] +
        ["QHadron_pdgId", "QHadron_charmOrigin"] +
        ["QGNNTrack_trackIdx", "QGNNTrack_genChargedIdx", "QGNNTrack_hadronIdx", "QGNNTrack_quarkIdx",
         "QGNNTrack_associationQuality", "QGNNTrack_membershipType", "QGNNTrack_fragmentationSourceCertain",
         "QGNNTrack_hadronClass", "QGNNTrack_hadronTruthValid"] +
        ["QGNNPairTruth_pairIdx", "QGNNPairTruth_hadronTruthValid", "QGNNPairTruth_sameHadron", "QGNNPairTruth_quarkTruthValid",
         "QGNNPairTruth_sameQuark"] +
        ["QGNNHadron_hadronIdx", "QGNNHadron_nAcceptedDecayTracks", "QGNNHadron_nCategory"]
    )
    with uproot.open(args.input) as root:
        tree = root["Events"]
        missing = [name for name in required if name not in tree.keys()]
        if missing:
            raise ValueError("missing required output branches: " + ", ".join(missing))
        arrays = tree.arrays(required, library="ak")
        totals = {"tracks": 0, "pairs": 0, "hadrons": 0}
        for event in range(tree.num_entries):
            reco_idx = values(arrays, "QRecoTrack_trackIdx", event)
            n = len(reco_idx)
            if reco_idx != list(range(n)):
                raise ValueError("QRecoTrack indices are not the unchanged unpacked-track row order in event %d" % event)
            for name in node_fields + ["track_hadidx", "track_label"]:
                if len(values(arrays, name, event)) != n:
                    raise ValueError("track feature/truth row count mismatch: %s event %d" % (name, event))
            qgnn_idx = values(arrays, "QGNNTrack_trackIdx", event)
            if qgnn_idx != reco_idx:
                raise ValueError("QGNNTrack and QRecoTrack row indices differ in event %d" % event)
            qgen_idx = values(arrays, "QRecoTrack_genChargedIdx", event)
            reco_pt = values(arrays, "QRecoTrack_pt", event)
            reco_eta = values(arrays, "QRecoTrack_eta", event)
            qgen_n = len(values(arrays, "QGenCharged_membershipType", event))
            if any(i < -1 or i >= qgen_n for i in qgen_idx):
                raise ValueError("QRecoTrack has out-of-range QGenCharged index in event %d" % event)
            for suffix in ("genChargedIdx", "hadronIdx", "quarkIdx", "associationQuality", "membershipType",
                           "fragmentationSourceCertain", "hadronClass", "hadronTruthValid"):
                if len(values(arrays, "QGNNTrack_" + suffix, event)) != n:
                    raise ValueError("QGNNTrack row count mismatch for %s in event %d" % (suffix, event))
            if values(arrays, "QGNNTrack_genChargedIdx", event) != qgen_idx:
                raise ValueError("QGNNTrack does not preserve QRecoTrack assignments in event %d" % event)
            qmembership = values(arrays, "QGenCharged_membershipType", event)
            qhadron = values(arrays, "QGenCharged_hadronIdx", event)
            qtrack_hadron = values(arrays, "QGNNTrack_hadronIdx", event)
            nhad = len(values(arrays, "QHadron_pdgId", event))
            expected_n = [0] * nhad
            for i, ig in enumerate(qgen_idx):
                expected_hadron = -1
                if ig >= 0 and 0 <= qhadron[ig] < nhad:
                    expected_hadron = qhadron[ig]
                    if qmembership[ig] == 1 and reco_pt[i] >= 0.4 and abs(reco_eta[i]) <= 2.5:
                        expected_n[expected_hadron] += 1
                if qtrack_hadron[i] != expected_hadron:
                    raise ValueError("QGNNTrack hadron join mismatch in event %d row %d" % (event, i))

            src = values(arrays, "trackPair_trk1Idx", event)
            dst = values(arrays, "trackPair_trk2Idx", event)
            npair = len(src)
            if len(dst) != npair or any(i < 0 or j >= n or i >= j for i, j in zip(src, dst)):
                raise ValueError("trackPair has invalid physical-pair endpoints in event %d" % event)
            for name in edge_fields + ["trackPair_edgeLabel"]:
                if len(values(arrays, name, event)) != npair:
                    raise ValueError("trackPair feature/label row count mismatch: %s event %d" % (name, event))
            if values(arrays, "QGNNPairTruth_pairIdx", event) != list(range(npair)):
                raise ValueError("QGNNPairTruth row order does not match trackPair in event %d" % event)
            had_valid = values(arrays, "QGNNPairTruth_hadronTruthValid", event)
            same_had = values(arrays, "QGNNPairTruth_sameHadron", event)
            q_valid = values(arrays, "QGNNPairTruth_quarkTruthValid", event)
            same_q = values(arrays, "QGNNPairTruth_sameQuark", event)
            if any((v == 0 and y != -1) or (v == 1 and y not in (0, 1)) for v, y in zip(had_valid, same_had)):
                raise ValueError("invalid sameHadron sentinel/mask in event %d" % event)
            if any((v == 0 and y != -1) or (v == 1 and y not in (0, 1)) for v, y in zip(q_valid, same_q)):
                raise ValueError("invalid sameQuark sentinel/mask in event %d" % event)
            gquark = values(arrays, "QGNNTrack_quarkIdx", event)
            gquality = values(arrays, "QGNNTrack_associationQuality", event)
            gcertain = values(arrays, "QGNNTrack_fragmentationSourceCertain", event)
            if any(y == 1 and not (gquality[i] == gquality[j] == 1 and gcertain[i] == gcertain[j] == 1 and
                                   gquark[i] >= 0 and gquark[i] == gquark[j])
                   for y, i, j in zip(same_q, src, dst)):
                raise ValueError("sameQuark positive violates source-certainty rule in event %d" % event)

            hidx = values(arrays, "QGNNHadron_hadronIdx", event)
            hn = values(arrays, "QGNNHadron_nAcceptedDecayTracks", event)
            hcat = values(arrays, "QGNNHadron_nCategory", event)
            if hidx != list(range(nhad)) or len(hn) != nhad or len(hcat) != nhad or hn != expected_n:
                raise ValueError("QGNNHadron is not aligned with QHadron in event %d" % event)
            if any(c != (0 if count == 0 else 1 if count == 1 else 2 if count == 2 else 3)
                   for count, c in zip(hn, hcat)):
                raise ValueError("invalid N-category in event %d" % event)
            totals["tracks"] += n
            totals["pairs"] += npair
            totals["hadrons"] += nhad
        print("Validated %d events: %d track rows, %d historical physical pairs, %d hadron rows" %
              (tree.num_entries, totals["tracks"], totals["pairs"], totals["hadrons"]))


if __name__ == "__main__":
    main()
