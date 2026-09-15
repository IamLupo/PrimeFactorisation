#!/usr/bin/env python3

import math
import time
from dataclasses import dataclass
from typing import Dict, Tuple, List

import numpy as np


# ============================================================
# START EXPERIMENT 198
#
# Bit-packed QR beam search for Fermat factorization.
#
# Fermat searches:
#
#     x = ceil(sqrt(n)), ..., (p+q)/2
#
# and tests:
#
#     x^2 - n = y^2
#
# Therefore for every odd prime l:
#
#     x^2 - n
#
# must be a quadratic residue modulo l.
#
# This experiment:
#
#   * works on the correct Fermat x variable
#   * precomputes QR masks
#   * packs masks into bits
#   * combines candidates with bitwise AND
#   * uses measured bitset runtime as the primary score
#   * uses survivor count as a secondary score
#   * uses beam search instead of greedy search
#
# ============================================================


# ------------------------------------------------------------
# TEST INSTANCES
# ------------------------------------------------------------

INSTANCES = {
    48: {
        "p": 8_390_069,
        "q": 33_547_589,
        "n": 281_466_586_493_641,
    },
    54: {
        "p": 124_517_461,
        "q": 144_517_463,
        "n": 17_994_947_562_921_443,
    },
    60: {
        "p": 1_058_841_403,
        "q": 1_088_841_421,
        "n": 1_152_910_377_856_153_663,
    },
    66: {
        "p": 8_569_934_017,
        "q": 8_609_934_041,
        "n": 73_786_566_622_092_172_697,
    },
}


# ------------------------------------------------------------
# Candidate primes
# ------------------------------------------------------------

PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
]


# ------------------------------------------------------------
# SEARCH PARAMETERS
# ------------------------------------------------------------

BEAM_WIDTH = 16

MAX_DEPTH = 12

# Calibration interval in Fermat x values.
CALIBRATION_SIZE = 16 * 1024 * 1024

# Repeated AND operations used to make timing measurable.
# This deliberately performs the same logical sieve several
# times, but does NOT rebuild the prime masks.
TIMING_ROUNDS = 20

# Number of final candidates retained from every depth.
FINALISTS_PER_DEPTH = 4

# Repeat actual final benchmarks.
FULL_REPEATS = 3

# Print this many beam members.
PRINT_TOP = 10


# ------------------------------------------------------------
# DATA STRUCTURES
# ------------------------------------------------------------

@dataclass(frozen=True)
class PrimeMask:
    prime: int

    # Packed bits.
    calibration: np.ndarray

    # Packed bits for the actual Fermat interval.
    full: np.ndarray


@dataclass
class Candidate:
    primes: Tuple[int, ...]

    # Measured calibration time.
    calibration_time: float

    # Number of surviving calibration x values.
    calibration_survivors: int


# ------------------------------------------------------------
# FER​​MAT HELPERS
# ------------------------------------------------------------

def fermat_x_start(n: int) -> int:
    """
    Smallest integer x such that x^2 > n.
    """
    return math.isqrt(n - 1) + 1


def planted_x(p: int, q: int) -> int:
    """
    True Fermat x = (p+q)/2.
    """
    if (p + q) & 1:
        raise ValueError("p and q do not produce integer Fermat x")

    return (p + q) // 2


def planted_y(p: int, q: int) -> int:
    """
    True Fermat y = (q-p)/2.
    """
    if (q - p) & 1:
        raise ValueError("p and q do not produce integer Fermat y")

    return (q - p) // 2


# ------------------------------------------------------------
# QUADRATIC-RESIDUE TABLE
# ------------------------------------------------------------

def build_qr_allowed_table(prime: int, n: int) -> np.ndarray:
    """
    allowed[r] == 1 iff

        r^2 - n

    is a quadratic residue modulo prime.
    """

    p = prime

    squares = np.zeros(p, dtype=np.uint8)

    r = np.arange(p, dtype=np.int64)

    squares[(r * r) % p] = 1

    residues = np.arange(p, dtype=np.int64)

    values = (
        (residues * residues - (n % p))
        % p
    )

    allowed = squares[values]

    return allowed


# ------------------------------------------------------------
# BUILD A PACKED MASK FOR AN INTERVAL
# ------------------------------------------------------------

def build_packed_mask(
    start_x: int,
    length: int,
    prime: int,
    qr_table: np.ndarray,
) -> np.ndarray:
    """
    Build bit-packed survival mask for:

        x = start_x ... start_x + length - 1
    """

    x = np.arange(
        start_x,
        start_x + length,
        dtype=np.int64,
    )

    residues = x % prime

    allowed = qr_table[residues]

    return np.packbits(
        allowed,
        bitorder="little",
    )


# ------------------------------------------------------------
# POPCOUNT LOOKUP
# ------------------------------------------------------------

POPCOUNT8 = np.array(
    [bin(i).count("1") for i in range(256)],
    dtype=np.uint8,
)


def popcount_bytes(bits: np.ndarray) -> int:
    """
    Count set bits in packed uint8 representation.
    """
    return int(
        POPCOUNT8[bits].sum(dtype=np.int64)
    )


# ------------------------------------------------------------
# MASK CACHE CONSTRUCTION
# ------------------------------------------------------------

def build_mask_cache(
    n: int,
    start_x: int,
    full_length: int,
) -> Dict[int, PrimeMask]:

    print()
    print("Building packed QR masks...")
    print(
        "  calibration values =",
        CALIBRATION_SIZE,
    )
    print(
        "  full values        =",
        full_length,
    )

    cache: Dict[int, PrimeMask] = {}

    for prime in PRIMES:

        qr = build_qr_allowed_table(
            prime,
            n,
        )

        calibration_mask = build_packed_mask(
            start_x,
            CALIBRATION_SIZE,
            prime,
            qr,
        )

        full_mask = build_packed_mask(
            start_x,
            full_length,
            prime,
            qr,
        )

        cache[prime] = PrimeMask(
            prime=prime,
            calibration=calibration_mask,
            full=full_mask,
        )

        density = (
            int(calibration_mask.sum())
        )

        print(
            f"  prime={prime:2d} "
            f"mask_bytes={len(full_mask):,}"
        )

    return cache


# ------------------------------------------------------------
# INITIAL ALL-ONES MASK
# ------------------------------------------------------------

def make_full_ones(length: int) -> np.ndarray:
    """
    Packed bitset containing `length` ones.
    """

    size = (length + 7) // 8

    bits = np.full(
        size,
        0xFF,
        dtype=np.uint8,
    )

    # Clear padding bits beyond length.
    extra = size * 8 - length

    if extra:
        bits[-1] &= (
            0xFF >> extra
        )

    return bits


# ------------------------------------------------------------
# APPLY PRIME SUBSET
# ------------------------------------------------------------

def apply_subset(
    primes: Tuple[int, ...],
    masks: Dict[int, np.ndarray],
    initial: np.ndarray,
    timing_rounds: int = 1,
) -> Tuple[float, np.ndarray]:
    """
    Apply:

        result &= mask_prime

    for every prime in the subset.

    Multiple timing rounds make very short bitwise operations
    measurable.
    """

    result = None
    elapsed = 0.0

    for _ in range(timing_rounds):

        current = initial.copy()

        t0 = time.perf_counter()

        for prime in primes:
            current &= masks[prime]

        elapsed += time.perf_counter() - t0

        result = current

    assert result is not None

    return elapsed / timing_rounds, result


# ------------------------------------------------------------
# CALIBRATION
# ------------------------------------------------------------

def benchmark_calibration(
    primes: Tuple[int, ...],
    mask_cache: Dict[int, PrimeMask],
    initial: np.ndarray,
) -> Tuple[float, int]:

    masks = {
        p: mask_cache[p].calibration
        for p in primes
    }

    dt, result = apply_subset(
        primes,
        masks,
        initial,
        timing_rounds=TIMING_ROUNDS,
    )

    survivors = popcount_bytes(result)

    return dt, survivors


# ------------------------------------------------------------
# FULL-RANGE BENCHMARK
# ------------------------------------------------------------

def benchmark_full(
    primes: Tuple[int, ...],
    mask_cache: Dict[int, PrimeMask],
    initial: np.ndarray,
    repeats: int,
) -> Tuple[float, int]:

    masks = {
        p: mask_cache[p].full
        for p in primes
    }

    best_time = float("inf")
    best_survivors = 0

    for _ in range(repeats):

        t0 = time.perf_counter()

        current = initial.copy()

        for prime in primes:
            current &= masks[prime]

        elapsed = (
            time.perf_counter()
            - t0
        )

        if elapsed < best_time:
            best_time = elapsed
            best_survivors = popcount_bytes(
                current
            )

    return best_time, best_survivors


# ------------------------------------------------------------
# TRUE X SURVIVAL
# ------------------------------------------------------------

def survives_at_x(
    x: int,
    primes: Tuple[int, ...],
    mask_cache: Dict[int, PrimeMask],
    start_x: int,
) -> bool:

    offset = x - start_x

    if offset < 0:
        return False

    byte_index = offset >> 3
    bit_index = offset & 7

    bit = 1 << bit_index

    for prime in primes:

        mask = mask_cache[prime].full

        if not (
            mask[byte_index]
            & bit
        ):
            return False

    return True


# ------------------------------------------------------------
# BEAM SEARCH
# ------------------------------------------------------------

def beam_search(
    n: int,
    mask_cache: Dict[int, PrimeMask],
    calibration_initial: np.ndarray,
) -> List[Candidate]:

    print()
    print("=" * 78)
    print("BEAM SEARCH")
    print("=" * 78)

    root = Candidate(
        primes=tuple(),
        calibration_time=0.0,
        calibration_survivors=CALIBRATION_SIZE,
    )

    beam = [root]

    visited = {tuple()}

    all_candidates: List[Candidate] = []

    for depth in range(
        1,
        MAX_DEPTH + 1,
    ):

        expanded: List[Candidate] = []

        for parent in beam:

            last = (
                parent.primes[-1]
                if parent.primes
                else 0
            )

            for prime in PRIMES:

                if prime <= last:
                    continue

                candidate_primes = (
                    parent.primes + (prime,)
                )

                if candidate_primes in visited:
                    continue

                visited.add(candidate_primes)

                dt, survivors = benchmark_calibration(
                    candidate_primes,
                    mask_cache,
                    calibration_initial,
                )

                expanded.append(
                    Candidate(
                        primes=candidate_primes,
                        calibration_time=dt,
                        calibration_survivors=survivors,
                    )
                )

        if not expanded:
            break

        # ----------------------------------------------------
        # IMPORTANT:
        #
        # Runtime is PRIMARY.
        #
        # Survivor count is SECONDARY.
        #
        # This avoids repeating the mistake where density was
        # treated as equivalent to runtime.
        # ----------------------------------------------------

        expanded.sort(
            key=lambda c: (
                c.calibration_time,
                c.calibration_survivors,
            )
        )

        beam = expanded[:BEAM_WIDTH]

        all_candidates.extend(beam)

        print()
        print("-" * 78)
        print(f"DEPTH {depth}")
        print("-" * 78)

        for rank, candidate in enumerate(
            beam[:PRINT_TOP],
            1,
        ):

            print(
                f"{rank:2d}. "
                f"{str(list(candidate.primes)):<48} "
                f"time={candidate.calibration_time:.9f}s "
                f"survivors={candidate.calibration_survivors:,}"
            )

    return all_candidates


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 198")
    print("Bit-packed QR beam search")
    print("=" * 78)

    for bits in sorted(INSTANCES):

        inst = INSTANCES[bits]

        p = inst["p"]
        q = inst["q"]
        n = inst["n"]

        start_x = fermat_x_start(n)
        true_x = planted_x(p, q)
        true_y = planted_y(p, q)

        full_length = (
            true_x
            - start_x
            + 1
        )

        print()
        print("#" * 78)
        print(f"INSTANCE {bits}-BIT")
        print("#" * 78)

        print("p       =", p)
        print("q       =", q)
        print("n       =", n)
        print("start_x =", start_x)
        print("true_x  =", true_x)
        print("true_y  =", true_y)
        print("gap     =", true_x - start_x)
        print("range   =", full_length)

        if true_x < start_x:
            raise RuntimeError(
                "FATAL: true_x < start_x"
            )

        # ----------------------------------------------------
        # Build masks.
        # ----------------------------------------------------

        mask_cache = build_mask_cache(
            n,
            start_x,
            full_length,
        )

        # ----------------------------------------------------
        # Initial masks.
        # ----------------------------------------------------

        calibration_initial = make_full_ones(
            CALIBRATION_SIZE
        )

        full_initial = make_full_ones(
            full_length
        )

        # ----------------------------------------------------
        # Verify true x against every prime.
        # ----------------------------------------------------

        rejected = []

        for prime in PRIMES:

            if not survives_at_x(
                true_x,
                (prime,),
                mask_cache,
                start_x,
            ):
                rejected.append(prime)

        if rejected:

            print()
            print(
                "ERROR: true_x rejected by:",
                rejected,
            )

            raise RuntimeError(
                "QR mask construction is inconsistent."
            )

        print()
        print(
            "True Fermat x survives all individual QR tests."
        )

        # ----------------------------------------------------
        # Baseline.
        # ----------------------------------------------------

        baseline_t0 = time.perf_counter()

        baseline_result = (
            full_initial.copy()
        )

        baseline_elapsed = (
            time.perf_counter()
            - baseline_t0
        )

        baseline_survivors = popcount_bytes(
            baseline_result
        )

        print()
        print("BASELINE")
        print(
            f"  no-sieve copy time = "
            f"{baseline_elapsed:.9f}s"
        )
        print(
            f"  candidates        = "
            f"{baseline_survivors:,}"
        )

        # ----------------------------------------------------
        # Beam search.
        # ----------------------------------------------------

        candidates = beam_search(
            n,
            mask_cache,
            calibration_initial,
        )

        # ----------------------------------------------------
        # Deduplicate.
        # ----------------------------------------------------

        unique = {}

        for candidate in candidates:

            old = unique.get(
                candidate.primes
            )

            if old is None:

                unique[
                    candidate.primes
                ] = candidate

            else:

                if (
                    candidate.calibration_time,
                    candidate.calibration_survivors,
                ) < (
                    old.calibration_time,
                    old.calibration_survivors,
                ):
                    unique[
                        candidate.primes
                    ] = candidate

        candidates = list(
            unique.values()
        )

        # ----------------------------------------------------
        # Select finalists.
        # ----------------------------------------------------

        by_depth: Dict[
            int,
            List[Candidate]
        ] = {}

        for candidate in candidates:

            by_depth.setdefault(
                len(candidate.primes),
                []
            ).append(candidate)

        finalists = {}

        for depth, values in by_depth.items():

            values.sort(
                key=lambda c: (
                    c.calibration_time,
                    c.calibration_survivors,
                )
            )

            for candidate in values[
                :FINALISTS_PER_DEPTH
            ]:
                finalists[
                    candidate.primes
                ] = candidate

        # ----------------------------------------------------
        # Full-range validation.
        # ----------------------------------------------------

        print()
        print("=" * 78)
        print("FULL-RANGE VALIDATION")
        print("=" * 78)

        results = []

        for candidate in finalists.values():

            elapsed, survivors = benchmark_full(
                candidate.primes,
                mask_cache,
                full_initial,
                FULL_REPEATS,
            )

            target_survives = survives_at_x(
                true_x,
                candidate.primes,
                mask_cache,
                start_x,
            )

            results.append(
                (
                    elapsed,
                    candidate,
                    survivors,
                    target_survives,
                )
            )

            print(
                f"{str(list(candidate.primes)):<50} "
                f"full={elapsed:.9f}s "
                f"survivors={survivors:,} "
                f"target={target_survives}"
            )

        # ----------------------------------------------------
        # Sort actual results.
        # ----------------------------------------------------

        results.sort(
            key=lambda x: x[0]
        )

        print()
        print("-" * 78)
        print("BEST FULL-RANGE RESULTS")
        print("-" * 78)

        for rank, (
            elapsed,
            candidate,
            survivors,
            target_survives,
        ) in enumerate(
            results[:10],
            1,
        ):

            print(
                f"{rank:2d}. "
                f"{str(list(candidate.primes)):<50} "
                f"{elapsed:.9f}s "
                f"survivors={survivors:,} "
                f"target={target_survives}"
            )

        # ----------------------------------------------------
        # Best valid candidate.
        # ----------------------------------------------------

        valid = [
            r
            for r in results
            if r[3]
        ]

        print()

        if valid:

            best = valid[0]

            elapsed = best[0]
            candidate = best[1]
            survivors = best[2]

            print(
                "BEST TARGET-PRESERVING SUBSET"
            )
            print(
                "  primes    =",
                list(candidate.primes),
            )
            print(
                "  runtime   =",
                elapsed,
            )
            print(
                "  survivors =",
                survivors,
            )

        else:

            print(
                "WARNING: no target-preserving "
                "finalist found."
            )

        # ----------------------------------------------------
        # Print beam winner for comparison.
        # ----------------------------------------------------

        if candidates:

            beam_best = min(
                candidates,
                key=lambda c: (
                    c.calibration_time,
                    c.calibration_survivors,
                )
            )

            print()
            print(
                "BEST CALIBRATION CANDIDATE"
            )
            print(
                "  primes    =",
                list(beam_best.primes),
            )
            print(
                "  time      =",
                beam_best.calibration_time,
            )
            print(
                "  survivors =",
                beam_best.calibration_survivors,
            )

        # ----------------------------------------------------
        # Direct verification of the mathematical factor.
        # ----------------------------------------------------

        check = (
            true_x * true_x
            - n
        )

        expected_y2 = (
            true_y * true_y
        )

        print()
        print("FERMAT SANITY CHECK")
        print(
            "  x^2 - n =", check
        )
        print(
            "  y^2     =", expected_y2
        )
        print(
            "  equal   =",
            check == expected_y2,
        )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 198")
    print("=" * 78)


if __name__ == "__main__":
    main()

