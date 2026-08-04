#!/usr/bin/env python3
"""Independent exhaustive verifier for the GF(2) <15,15,17> rank-2256 certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Iterable

N, M, P = 15, 15, 17
EXPECTED_RANK = 2256
AB = N * M
BC = M * P
AC = N * P


def set_bits(value: int) -> Iterable[int]:
    while value:
        least = value & -value
        yield least.bit_length() - 1
        value ^= least


def parse_certificate(path: Path) -> tuple[int, list[tuple[int, int, int]]]:
    declared_rank: int | None = None
    terms: list[tuple[int, int, int]] = []

    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if fields[0] == "R":
            fields = fields[1:]
        if declared_rank is None:
            if len(fields) != 1:
                raise ValueError(f"{path}:{line_number}: expected the rank header")
            declared_rank = int(fields[0], 10)
            continue
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: expected three decimal masks")
        try:
            u, v, w = (int(field, 10) for field in fields)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid decimal mask") from exc
        terms.append((u, v, w))

    if declared_rank is None:
        raise ValueError(f"{path}: missing rank header")
    return declared_rank, terms


def verify(path: Path) -> dict[str, object]:
    started = time.perf_counter()
    declared_rank, terms = parse_certificate(path)

    limit_u = 1 << AB
    limit_v = 1 << BC
    limit_w = 1 << AC
    bad_terms = [
        index
        for index, (u, v, w) in enumerate(terms)
        if not (0 < u < limit_u and 0 < v < limit_v and 0 < w < limit_w)
    ]
    duplicate_triples = len(terms) - len(set(terms))

    residual = [0] * (AB * BC)
    expanded_uv_pairs = 0
    for u, v, w in terms:
        u_bits = tuple(set_bits(u))
        v_bits = tuple(set_bits(v))
        expanded_uv_pairs += len(u_bits) * len(v_bits)
        for a in u_bits:
            base = a * BC
            for b in v_bits:
                residual[base + b] ^= w

    for i in range(N):
        for j in range(M):
            a = i * M + j
            base = a * BC
            for k in range(P):
                b = j * P + k
                residual[base + b] ^= 1 << (i * P + k)

    residual_pair_slices = sum(mask != 0 for mask in residual)
    residual_coefficients = sum(mask.bit_count() for mask in residual)
    verified = (
        declared_rank == EXPECTED_RANK
        and len(terms) == EXPECTED_RANK
        and not bad_terms
        and duplicate_triples == 0
        and residual_coefficients == 0
    )

    return {
        "certificate": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "tensor": [N, M, P],
        "expected_rank": EXPECTED_RANK,
        "declared_rank": declared_rank,
        "parsed_terms": len(terms),
        "nonzero_and_in_range": not bad_terms,
        "bad_term_indices": bad_terms[:20],
        "duplicate_triples": duplicate_triples,
        "expanded_uv_pairs": expanded_uv_pairs,
        "target_coefficients": N * M * P,
        "brent_coefficients": AB * BC * AC,
        "residual_pair_slices": residual_pair_slices,
        "residual_coefficients": residual_coefficients,
        "verified": verified,
        "seconds": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name("matmul_15x15x17_rank2256_gf2.txt"),
    )
    args = parser.parse_args()
    try:
        result = verify(args.certificate)
    except (OSError, ValueError) as exc:
        print(json.dumps({"verified": False, "error": str(exc)}, indent=2))
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    sys.exit(main())
