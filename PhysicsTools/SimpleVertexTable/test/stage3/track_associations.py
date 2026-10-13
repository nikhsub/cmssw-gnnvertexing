"""Association helpers for Stage 2 QRecoTrack to Stage 1 truth tables."""

from __future__ import annotations

from truth_objects import MEMBERSHIP_HEAVY_DECAY, MEMBERSHIP_ADDITIONAL, successful_match


def matched_reco_by_gen(reco_rows: list[dict]) -> dict[int, dict]:
    """Return the successful one-to-one reco match for each QGenCharged row."""
    result = {}
    for row in reco_rows:
        if not successful_match(row):
            continue
        gen_idx = int(row["genChargedIdx"])
        if gen_idx in result:
            raise ValueError("Stage 2 one-to-one assignment duplicated QGenCharged row %d" % gen_idx)
        result[gen_idx] = row
    return result


def heavy_decay_gen_rows(gen_rows: list[dict], hadron_idx: int) -> list[int]:
    """QGenCharged rows classified as decay descendants of one QHadron row."""
    return [i for i, row in enumerate(gen_rows)
            if int(row["membershipType"]) == MEMBERSHIP_HEAVY_DECAY
            and int(row["hadronIdx"]) == int(hadron_idx)]


def additional_gen_rows(gen_rows: list[dict], quark_idx: int) -> list[int]:
    """Unique ancestry candidates for a QFrag row, excluding all heavy decays.

    Stage 1's membershipType=2 is retained ancestry evidence, not guaranteed
    exclusive physical fragmentation ownership. The source-certainty branch
    is returned unchanged for every consumer to preserve that distinction.
    """
    if int(quark_idx) < 0:
        return []
    return [i for i, row in enumerate(gen_rows)
            if int(row["membershipType"]) == MEMBERSHIP_ADDITIONAL
            and int(row["associationQuality"]) == 1
            and int(row["quarkIdx"]) == int(quark_idx)]

