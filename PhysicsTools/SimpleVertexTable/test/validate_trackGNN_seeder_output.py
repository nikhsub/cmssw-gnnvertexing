#!/usr/bin/env python3
"""Integrity checks for TrackGNNSeederInferenceProducer NanoAOD output.

Requires uproot and awkward. This checks ROOT table alignment and value
validity only; it does not claim PyTorch/CMSSW numerical parity.
"""
import argparse
import math
import statistics
import sys

try:
    import awkward as ak
    import uproot
except ImportError as exc:
    raise SystemExit("Install/import uproot and awkward in the validation environment: %s" % exc)


def as_list(array, event):
    return ak.to_list(array[event])


def finite_values(values):
    return [float(x) for x in values if isinstance(x, (int, float)) and math.isfinite(float(x))]


def summarize(values):
    vals = finite_values(values)
    if not vals:
        return "n=0"
    return "n=%d min=%.6g mean=%.6g max=%.6g" % (len(vals), min(vals), statistics.fmean(vals), max(vals))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file")
    parser.add_argument("--tree", default="Events")
    args = parser.parse_args()

    errors = []
    with uproot.open(args.root_file) as root_file:
        if args.tree not in root_file:
            raise SystemExit("Missing tree %r; keys: %s" % (args.tree, root_file.keys()))
        tree = root_file[args.tree]
        names = set(tree.keys())
        required = [
            "ntrack", "track_pt", "ntrackPair", "trackPair_trk1Idx", "trackPair_trk2Idx",
            "trackPair_seederLogit", "trackPair_seederScore", "trackPair_seederValid",
            "nSeederNode", "SeederNode_trackIdx", "SeederNode_valid",
        ]
        required += ["SeederNode_logit%d" % i for i in range(6)]
        required += ["SeederNode_prob%d" % i for i in range(6)]
        missing = [name for name in required if name not in names]
        if missing:
            raise SystemExit("Missing expected NanoAOD branches: %s\nAvailable relevant branches: %s" %
                             (missing, sorted(n for n in names if "Seeder" in n or "trackPair" in n or n.startswith("ntrack"))))

        embedding_branches = ["nSeederEmbedding", "SeederEmbedding_trackIdx", "SeederEmbedding_valid"]
        embedding_branches += ["SeederEmbedding_emb%d" % i for i in range(128)]
        has_embeddings = all(name in names for name in embedding_branches)
        partial_embedding = any(name in names for name in embedding_branches) and not has_embeddings
        if partial_embedding:
            errors.append("SeederEmbedding table is incomplete")

        read = required + (embedding_branches if has_embeddings else [])
        arrays = tree.arrays(read, library="ak")
        n_events = tree.num_entries
        edge_scores, node_probs = [], [[] for _ in range(6)]
        node_counts, edge_counts = [], []
        invalid_nodes = invalid_edges = total_nodes = total_edges = no_edge_events = 0
        checked_edges = 0
        for ev in range(n_events):
            ntrack = int(arrays["ntrack"][ev])
            npair = int(arrays["ntrackPair"][ev])
            nnode = int(arrays["nSeederNode"][ev])
            tracks_pt = as_list(arrays["track_pt"], ev)
            edge_i = as_list(arrays["trackPair_trk1Idx"], ev)
            edge_j = as_list(arrays["trackPair_trk2Idx"], ev)
            logits = as_list(arrays["trackPair_seederLogit"], ev)
            scores = as_list(arrays["trackPair_seederScore"], ev)
            edge_valid = as_list(arrays["trackPair_seederValid"], ev)
            node_idx = as_list(arrays["SeederNode_trackIdx"], ev)
            node_valid = as_list(arrays["SeederNode_valid"], ev)
            node_logits = [as_list(arrays["SeederNode_logit%d" % k], ev) for k in range(6)]
            probs = [as_list(arrays["SeederNode_prob%d" % k], ev) for k in range(6)]

            node_counts.append(nnode)
            edge_counts.append(npair)
            total_nodes += nnode
            total_edges += npair
            if npair == 0:
                no_edge_events += 1
            if len(tracks_pt) != ntrack or nnode != ntrack:
                errors.append("event %d: ntrack=%d track rows=%d nSeederNode=%d" % (ev, ntrack, len(tracks_pt), nnode))
            if not (len(edge_i) == len(edge_j) == len(logits) == len(scores) == len(edge_valid) == npair):
                errors.append("event %d: trackPair row counts differ from ntrackPair=%d" % (ev, npair))
            if not (len(node_idx) == len(node_valid) == nnode):
                errors.append("event %d: SeederNode row count mismatch" % ev)
            if node_idx != list(range(ntrack)):
                errors.append("event %d: SeederNode_trackIdx is not 0..ntrack-1" % ev)

            directed = set()
            for row, (i, j, logit, score, valid) in enumerate(zip(edge_i, edge_j, logits, scores, edge_valid)):
                i, j, valid = int(i), int(j), int(valid)
                if i < 0 or i >= ntrack or j < 0 or j >= ntrack or i == j:
                    errors.append("event %d pair row %d: invalid directed endpoint (%d,%d)" % (ev, row, i, j))
                if (i, j) in directed:
                    errors.append("event %d: duplicate directed edge (%d,%d)" % (ev, i, j))
                directed.add((i, j))
                if valid:
                    if not math.isfinite(float(logit)) or not math.isfinite(float(score)):
                        errors.append("event %d pair row %d: valid edge has nonfinite output" % (ev, row))
                    elif not (0.0 <= float(score) <= 1.0):
                        errors.append("event %d pair row %d: score outside [0,1]" % (ev, row))
                    else:
                        edge_scores.append(float(score))
                elif float(logit) != -1.0 or float(score) != -1.0:
                    errors.append("event %d pair row %d: invalid edge does not use -1 sentinels" % (ev, row))
                    invalid_edges += 1
                else:
                    invalid_edges += 1
                checked_edges += 1
            missing_reciprocals = sum(1 for i, j in directed if (j, i) not in directed)
            if missing_reciprocals:
                errors.append("event %d: %d directed edges lack reciprocal rows" % (ev, missing_reciprocals))

            for row, valid in enumerate(node_valid):
                valid = int(valid)
                if valid:
                    values = [float(node_logits[k][row]) for k in range(6)]
                    p = [float(probs[k][row]) for k in range(6)]
                    if not all(math.isfinite(x) for x in values + p):
                        errors.append("event %d node row %d: valid output is nonfinite" % (ev, row))
                    elif abs(sum(p) - 1.0) > 1.e-4 or any(x < 0.0 or x > 1.0 for x in p):
                        errors.append("event %d node row %d: softmax probabilities invalid" % (ev, row))
                    else:
                        for k in range(6):
                            node_probs[k].append(p[k])
                else:
                    invalid_nodes += 1
                    if any(float(node_logits[k][row]) != -1.0 or float(probs[k][row]) != -1.0 for k in range(6)):
                        errors.append("event %d node row %d: invalid output does not use -1 sentinels" % (ev, row))

            if has_embeddings:
                nemb = int(arrays["nSeederEmbedding"][ev])
                emb_idx = as_list(arrays["SeederEmbedding_trackIdx"], ev)
                emb_valid = as_list(arrays["SeederEmbedding_valid"], ev)
                if nemb != ntrack or emb_idx != list(range(ntrack)) or len(emb_valid) != ntrack:
                    errors.append("event %d: embedding rows do not align with track rows" % ev)
                for row, valid in enumerate(emb_valid):
                    components = [float(as_list(arrays["SeederEmbedding_emb%d" % k], ev)[row]) for k in range(128)]
                    if int(valid):
                        if not all(math.isfinite(x) for x in components):
                            errors.append("event %d embedding row %d: valid embedding is nonfinite" % (ev, row))
                    elif any(x != -1.0 for x in components):
                        errors.append("event %d embedding row %d: invalid values do not use -1 sentinels" % (ev, row))

        print("Events:", n_events)
        print("Tracks/event:", summarize(node_counts))
        print("Directed edges/event:", summarize(edge_counts))
        print("Directed edge scores:", summarize(edge_scores))
        for k in range(6):
            print("Node class %d probabilities:" % k, summarize(node_probs[k]))
        print("Invalid node rows: %d / %d (%.6f)" %
              (invalid_nodes, total_nodes, invalid_nodes / total_nodes if total_nodes else 0.0))
        print("Invalid edge rows: %d / %d (%.6f)" %
              (invalid_edges, checked_edges, invalid_edges / checked_edges if checked_edges else 0.0))
        print("Events with no directed edges:", no_edge_events)
        print("Embedding table:", "present and aligned" if has_embeddings else "not present")
        if errors:
            print("Validation failed with %d issue(s):" % len(errors))
            for error in errors[:100]:
                print(" -", error)
            if len(errors) > 100:
                print(" - ... %d additional issue(s)" % (len(errors) - 100))
            return 1
        print("Output integrity validation passed")
        return 0


if __name__ == "__main__":
    sys.exit(main())
