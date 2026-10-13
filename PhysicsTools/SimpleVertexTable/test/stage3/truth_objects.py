"""Event-local truth-object and reco-track linking for the Stage 3 study."""

from __future__ import annotations

import math

PT_MIN = 0.4
ETA_MAX = 2.5

ORIGIN_NOT_APPLICABLE = 0
ORIGIN_PROMPT_D = 1
ORIGIN_FROM_B = 2
ORIGIN_UNKNOWN_D = 3

ASSOCIATION_MISSING = 0
ASSOCIATION_UNIQUE = 1
ASSOCIATION_AMBIGUOUS = 2
ASSOCIATION_UNKNOWN = 3

MEMBERSHIP_UNASSIGNED = 0
MEMBERSHIP_HEAVY_DECAY = 1
MEMBERSHIP_ADDITIONAL = 2
MEMBERSHIP_AMBIGUOUS = 3

MATCH_NO_CANDIDATE = 0
MATCH_UNOPPOSED = 1
MATCH_COMPETING = 2
MATCH_LOST_CONFLICT = 3

HADRON_FIELDS = (
    "prunedGenIdx", "pdgId", "pt", "eta", "phi", "mass", "charmOrigin",
    "parentBHadronIdx", "quarkIdx", "associationQuality", "fragmentationSourceCertain",
    "preHadPt", "preHadEta", "preHadPhi", "preHadKinematicsValid",
)
QUARK_FIELDS = ("pdgId", "pt", "eta", "phi", "mass", "prunedGenIdx", "hardPrunedIdx")
GEN_FIELDS = (
    "packedGenIdx", "pdgId", "charge", "pt", "eta", "phi", "mass", "lastPrunedIdx",
    "lastPrunedIsSelf", "hadronIdx", "quarkIdx", "associationQuality", "membershipType",
    "fragmentationSourceCertain", "isFromB", "isFromPromptD", "isFromBtoD",
)
RECO_FIELDS = (
    "trackIdx", "pt", "eta", "phi", "charge", "isHighPurity", "genChargedIdx", "matchDR",
    "matchRelPt", "matchStatus", "nGenCandidates",
)


def accepted(pt: float, eta: float) -> bool:
    """Analysis acceptance; boundaries are inclusive."""
    return math.isfinite(float(pt)) and math.isfinite(float(eta)) and float(pt) >= PT_MIN and abs(float(eta)) <= ETA_MAX


def origin_category(pdg_id: int, charm_origin: int) -> str:
    """Map the producer's QHadron codes to the analysis categories."""
    if int(charm_origin) == ORIGIN_NOT_APPLICABLE:
        return "B"
    return {
        ORIGIN_PROMPT_D: "D_prompt",
        ORIGIN_FROM_B: "D_fromB",
        ORIGIN_UNKNOWN_D: "D_unknown",
    }.get(int(charm_origin), "D_unknown")


def n_category(n: int) -> str:
    if n < 0:
        raise ValueError("track multiplicity cannot be negative")
    if n < 3:
        return "N%d" % n
    return "N3" if n == 3 else "N3+"


def validate_tables(tables: dict[str, dict[str, list]]) -> None:
    """Check schema lengths and all required event-local references."""
    for name, fields in (("QHadron", HADRON_FIELDS), ("QFrag", QUARK_FIELDS),
                         ("QGenCharged", GEN_FIELDS), ("QRecoTrack", RECO_FIELDS)):
        if name not in tables:
            raise ValueError("missing table " + name)
        columns = tables[name]
        missing = [field for field in fields if field not in columns]
        if missing:
            raise ValueError("%s missing columns: %s" % (name, ", ".join(missing)))
        lengths = {len(columns[field]) for field in fields}
        if len(lengths) != 1:
            raise ValueError("inconsistent column lengths in " + name)
    nh, nq, ng, nr = (len(tables[name][fields[0]]) for name, fields in (
        ("QHadron", HADRON_FIELDS), ("QFrag", QUARK_FIELDS),
        ("QGenCharged", GEN_FIELDS), ("QRecoTrack", RECO_FIELDS)))
    hadrons, gens, recos = tables["QHadron"], tables["QGenCharged"], tables["QRecoTrack"]
    for row, (origin, quality, certain) in enumerate(zip(hadrons["charmOrigin"], hadrons["associationQuality"], hadrons["fragmentationSourceCertain"])):
        if int(origin) not in (0, 1, 2, 3) or int(quality) not in (0, 1, 2, 3) or int(certain) not in (0, 1):
            raise ValueError("QHadron row %d has an unrecognized origin/quality/certainty code" % row)
    for row, (qidx, parent) in enumerate(zip(hadrons["quarkIdx"], hadrons["parentBHadronIdx"])):
        if int(qidx) < -1 or int(qidx) >= nq:
            raise ValueError("QHadron row %d has invalid event-local QFrag index %s" % (row, qidx))
        if int(parent) < -1 or int(parent) >= nh:
            raise ValueError("QHadron row %d has invalid event-local parentBHadronIdx %s" % (row, parent))
    for row, (hidx, qidx, quality, membership, certain) in enumerate(zip(
            gens["hadronIdx"], gens["quarkIdx"], gens["associationQuality"],
            gens["membershipType"], gens["fragmentationSourceCertain"])):
        if int(hidx) < -1 or int(hidx) >= nh:
            raise ValueError("QGenCharged row %d has invalid event-local QHadron index %s" % (row, hidx))
        if int(qidx) < -1 or int(qidx) >= nq:
            raise ValueError("QGenCharged row %d has invalid event-local QFrag index %s" % (row, qidx))
        if int(quality) not in (0, 1, 2, 3) or int(membership) not in (0, 1, 2, 3) or int(certain) not in (0, 1):
            raise ValueError("QGenCharged row %d has an unrecognized quality/membership/certainty code" % row)
    seen_gen = set()
    seen_track = set()
    for row, (track_idx, gen_idx, status) in enumerate(zip(recos["trackIdx"], recos["genChargedIdx"], recos["matchStatus"])):
        track_idx, gen_idx, status = int(track_idx), int(gen_idx), int(status)
        if track_idx < 0 or track_idx in seen_track:
            raise ValueError("QRecoTrack row %d has invalid or duplicate trackIdx %d" % (row, track_idx))
        seen_track.add(track_idx)
        if gen_idx < -1 or gen_idx >= ng:
            raise ValueError("QRecoTrack row %d has invalid event-local QGenCharged index %d" % (row, gen_idx))
        if status not in (0, 1, 2, 3):
            raise ValueError("QRecoTrack row %d has unknown matchStatus %d" % (row, status))
        if (status in (1, 2)) != (gen_idx >= 0):
            raise ValueError("QRecoTrack row %d has inconsistent matchStatus/genChargedIdx sentinels" % row)
        if gen_idx >= 0 and gen_idx in seen_gen:
            raise ValueError("QRecoTrack assigns QGenCharged row %d more than once" % gen_idx)
        if gen_idx >= 0:
            seen_gen.add(gen_idx)


def table_rows(columns: dict[str, list]) -> list[dict]:
    """Convert equal-length columns into row dictionaries."""
    if not columns:
        return []
    lengths = {len(values) for values in columns.values()}
    if len(lengths) > 1:
        raise ValueError("cannot make rows from unequal-length columns")
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def successful_match(reco_row: dict) -> bool:
    return int(reco_row["genChargedIdx"]) >= 0 and int(reco_row["matchStatus"]) in (MATCH_UNOPPOSED, MATCH_COMPETING)


def reco_confident(reco_row: dict) -> bool:
    """Strict Stage 2 subset: successful without any competing admissible edge."""
    return successful_match(reco_row) and int(reco_row["matchStatus"]) == MATCH_UNOPPOSED

