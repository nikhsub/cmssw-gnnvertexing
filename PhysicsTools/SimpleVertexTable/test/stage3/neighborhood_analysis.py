"""Reconstructed-quantity-only neighborhood and sparse-graph helpers."""

from __future__ import annotations

import math

import numpy as np


def delta_phi(phi1, phi2):
    return (np.asarray(phi1) - np.asarray(phi2) + np.pi) % (2.0 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    return np.hypot(np.asarray(eta1) - np.asarray(eta2), delta_phi(phi1, phi2))


def direction_from_tracks(rows: list[dict]):
    """Return eta, phi from the vector sum; (None, None) when undefined."""
    if not rows:
        return None, None
    px = sum(float(r["pt"]) * math.cos(float(r["phi"])) for r in rows)
    py = sum(float(r["pt"]) * math.sin(float(r["phi"])) for r in rows)
    pz = sum(float(r["pt"]) * math.sinh(float(r["eta"])) for r in rows)
    pt = math.hypot(px, py)
    momentum = math.sqrt(pt * pt + pz * pz)
    if not all(math.isfinite(x) for x in (px, py, pz, pt, momentum)) or pt <= 0 or momentum <= abs(pz):
        return None, None
    return math.asinh(pz / pt), math.atan2(py, px)


def upper_triangle_distances(track_rows: list[dict]):
    """Vectorized event-local pairwise deltaR for accepted reco tracks."""
    n = len(track_rows)
    if n < 2:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64))
    eta = np.asarray([r["eta"] for r in track_rows], dtype=np.float64)
    phi = np.asarray([r["phi"] for r in track_rows], dtype=np.float64)
    left, right = np.triu_indices(n, k=1)
    dr = delta_r(eta[left], phi[left], eta[right], phi[right])
    return left, right, dr


def nearest_neighbor_edge_pairs(left, right, dr, n_nodes: int, k: int = 8):
    """Undirected union of each node's k closest neighbors, using reco only."""
    if n_nodes < 2:
        return set()
    n_neighbors = min(max(int(k), 1), n_nodes - 1)
    distances = np.full((n_nodes, n_nodes), np.inf, dtype=np.float64)
    distances[left, right] = dr
    distances[right, left] = dr
    # Stable full sort makes exact deltaR ties reproducible by event-local index.
    indices = np.argsort(distances, axis=1, kind="stable")[:, :n_neighbors]
    return {(min(node, int(other)), max(node, int(other)))
            for node, neighbors in enumerate(indices) for other in neighbors}


def nearest_neighbor_edges(left, right, dr, n_nodes: int, k: int = 8):
    return len(nearest_neighbor_edge_pairs(left, right, dr, n_nodes, k))

