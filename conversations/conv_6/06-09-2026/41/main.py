#!/usr/bin/env python3

import math
import time
from dataclasses import dataclass
from typing import List, Tuple, Dict


# ============================================================
# START EXPERIMENT 197
# Corrected beam search over QR sieve prime subsets.
#
# IMPORTANT:
# Fermat factorization searches
#
#     x = ceil(sqrt(n)), ceil(sqrt(n))+1, ...
#
# and tests whether
#
#     x^2 - n = y^2
#
# is a square.
#
# Therefore the QR sieve is applied directly to x:
#
#     x^2 - n  must be a quadratic residue modulo l.
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
# Search parameters
# ------------------------------------------------------------

BEAM_WIDTH = 12
MAX_DEPTH = 12

# Number of x values used during calibration.
CALIBRATION_SIZE = 4_194_304

# Candidates retained for final full-range testing per depth.
FINALISTS_PER_DEPTH = 3

# Timing repetitions.
CALIBRATION_REPEATS = 1
FULL_REPEATS = 1

# Number printed from each beam.
PRINT_TOP = 8


# ------------------------------------------------------------
# Data structures
# ------------------------------------------------------------

@dataclass(frozen=True)
class PrimeMask:
    prime: int
    allowed: bytes


@dataclass
class Candidate:
    primes: Tuple[int, ...]
    calibration_time: float
    calibration_survivors: int


# ------------------------------------------------------------
# Fermat helpers
# ------------------------------------------------------------

def fermat_x_start(n: int) -> int:
    """
    Smallest integer x satisfying x^2 > n.

    This is ceil(sqrt(n)).
    """
    return math.isqrt(n - 1) + 1


def planted_x(p: int, q: int) -> int:
    """
    Fermat x for the known benchmark factors.
    """
    x = (p + q) // 2

    if (p + q) % 2 != 0:
        raise ValueError("p and q do not produce an integer Fermat x")

    return x


def planted_y(p: int, q: int) -> int:
    """
    Corresponding Fermat y.
    """
    y = (q - p) // 2

    if (q - p) % 2 != 0:
        raise ValueError("p and q do not produce an integer Fermat y")

    return y


# ------------------------------------------------------------
# QR mask construction
# ------------------------------------------------------------

def build_prime_mask(prime: int, n: int) -> PrimeMask:
    """
    For each residue r mod prime, determine whether

        r^2 - n

    is a quadratic residue modulo prime.

    A true Fermat x must satisfy

        x^2 - n = y^2,

    so every true x must survive this test.
    """

    p = prime

    # squares[t] = 1 iff t is a quadratic residue mod p
    squares = bytearray(p)

    for x in range(p):
        squares[(x * x) % p] = 1

    allowed = bytearray(p)

    nmod = n % p

    for r in range(p):
        value = (r * r - nmod) % p

        if squares[value]:
            allowed[r] = 1

    return PrimeMask(
        prime=p,
        allowed=bytes(allowed),
    )


# ------------------------------------------------------------
# Sieve
# ------------------------------------------------------------

def sieve_range(
    start_x: int,
    stop_x: int,
    masks: List[PrimeMask],
) -> Tuple[int, int]:
    """
    Segmented QR sieve on

        [start_x, stop_x)

    Returns:

        survivors
        total newly marked positions
    """

    length = stop_x - start_x

    if length <= 0:
        return 0, 0

    alive = bytearray(b"\x01") * length

    marked_total = 0

    for mask in masks:

        p = mask.prime
        allowed = mask.allowed

        # Residue corresponding to start_x.
        start_residue = start_x % p

        for i in range(length):

            if not alive[i]:
                continue

            r = (start_residue + i) % p

            if not allowed[r]:

                alive[i] = 0
                marked_total += 1

    survivors = sum(alive)

    return survivors, marked_total


# ------------------------------------------------------------
# True-x survival
# ------------------------------------------------------------

def survives_at_x(
    x: int,
    masks: List[PrimeMask],
) -> bool:

    for mask in masks:

        if not mask.allowed[x % mask.prime]:
            return False

    return True


# ------------------------------------------------------------
# Calibration benchmark
# ------------------------------------------------------------

def calibration_benchmark(
    n: int,
    primes: Tuple[int, ...],
    mask_cache: Dict[int, PrimeMask],
) -> Tuple[float, int]:

    start_x = fermat_x_start(n)
    stop_x = start_x + CALIBRATION_SIZE

    masks = [
        mask_cache[p]
        for p in primes
    ]

    best_time = float("inf")
    best_survivors = 0

    for _ in range(CALIBRATION_REPEATS):

        t0 = time.perf_counter()

        survivors, _ = sieve_range(
            start_x,
            stop_x,
            masks,
        )

        elapsed = time.perf_counter() - t0

        if elapsed < best_time:

            best_time = elapsed
            best_survivors = survivors

    return best_time, best_survivors


# ------------------------------------------------------------
# Beam search
# ------------------------------------------------------------

def beam_search(
    n: int,
    mask_cache: Dict[int, PrimeMask],
) -> List[Candidate]:

    print()
    print("Beam search starting point:")
    print("  x_start       =", fermat_x_start(n))
    print("  beam width    =", BEAM_WIDTH)
    print("  max depth     =", MAX_DEPTH)
    print("  calibration   =", CALIBRATION_SIZE)
    print()

    # Root candidate.
    root_time, root_survivors = calibration_benchmark(
        n,
        tuple(),
        mask_cache,
    )

    beam = [
        Candidate(
            primes=tuple(),
            calibration_time=root_time,
            calibration_survivors=root_survivors,
        )
    ]

    visited = {tuple()}

    all_candidates: List[Candidate] = []

    for depth in range(1, MAX_DEPTH + 1):

        expanded: List[Candidate] = []

        for parent in beam:

            last_prime = (
                parent.primes[-1]
                if parent.primes
                else 0
            )

            for prime in PRIMES:

                # Enforce canonical ordering.
                if prime <= last_prime:
                    continue

                candidate_primes = (
                    parent.primes + (prime,)
                )

                if candidate_primes in visited:
                    continue

                visited.add(candidate_primes)

                dt, survivors = calibration_benchmark(
                    n,
                    candidate_primes,
                    mask_cache,
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
        # Ranking
        #
        # First minimize survivors.
        # Then minimize measured calibration time.
        #
        # We intentionally do not use an analytic density model
        # as the primary ranking criterion.
        # ----------------------------------------------------

        expanded.sort(
            key=lambda c: (
                c.calibration_survivors,
                c.calibration_time,
                len(c.primes),
            )
        )

        beam = expanded[:BEAM_WIDTH]

        all_candidates.extend(beam)

        print()
        print("=" * 78)
        print(f"DEPTH {depth}")
        print("=" * 78)

        for rank, candidate in enumerate(
            beam[:PRINT_TOP],
            1,
        ):

            print(
                f"{rank:2d}. "
                f"{str(list(candidate.primes)):<45} "
                f"cal={candidate.calibration_time:.7f}s "
                f"survivors={candidate.calibration_survivors}"
            )

    return all_candidates


# ------------------------------------------------------------
# Full-range benchmark
# ------------------------------------------------------------

def full_benchmark(
    start_x: int,
    true_x: int,
    masks: List[PrimeMask],
) -> Tuple[float, int]:

    # Include the true Fermat x itself.
    stop_x = true_x + 1

    best_time = float("inf")
    best_survivors = 0

    for _ in range(FULL_REPEATS):

        t0 = time.perf_counter()

        survivors, _ = sieve_range(
            start_x,
            stop_x,
            masks,
        )

        elapsed = time.perf_counter() - t0

        if elapsed < best_time:

            best_time = elapsed
            best_survivors = survivors

    return best_time, best_survivors


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("START EXPERIMENT 197")
    print("Corrected Fermat-x QR beam search")
    print("=" * 78)

    for bits in sorted(INSTANCES):

        inst = INSTANCES[bits]

        p = inst["p"]
        q = inst["q"]
        n = inst["n"]

        start_x = fermat_x_start(n)
        true_x = planted_x(p, q)
        true_y = planted_y(p, q)

        gap = true_x - start_x

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
        print("gap     =", gap)

        if gap < 0:
            raise RuntimeError(
                "FATAL: Fermat target lies before search start."
            )

        # ----------------------------------------------------
        # Build masks.
        # ----------------------------------------------------

        print()
        print("Building QR masks...")

        mask_cache: Dict[int, PrimeMask] = {}

        for prime in PRIMES:

            mask_cache[prime] = build_prime_mask(
                prime,
                n,
            )

        print(
            "Masks built for",
            len(mask_cache),
            "primes."
        )

        # ----------------------------------------------------
        # Verify true x survives every individual test.
        # ----------------------------------------------------

        rejected = []

        for prime in PRIMES:

            mask = mask_cache[prime]

            if not survives_at_x(
                true_x,
                [mask],
            ):
                rejected.append(prime)

        if rejected:

            print()
            print(
                "ERROR: true Fermat x rejected by:",
                rejected,
            )

            print(
                "The QR sieve construction is inconsistent."
            )

            continue

        print()
        print("True Fermat x survives all individual QR tests.")

        # ----------------------------------------------------
        # Beam search.
        # ----------------------------------------------------

        all_candidates = beam_search(
            n,
            mask_cache,
        )

        # ----------------------------------------------------
        # Deduplicate.
        # ----------------------------------------------------

        unique: Dict[Tuple[int, ...], Candidate] = {}

        for candidate in all_candidates:

            old = unique.get(candidate.primes)

            if old is None:

                unique[candidate.primes] = candidate

            else:

                old_key = (
                    old.calibration_survivors,
                    old.calibration_time,
                )

                new_key = (
                    candidate.calibration_survivors,
                    candidate.calibration_time,
                )

                if new_key < old_key:
                    unique[candidate.primes] = candidate

        candidates = list(unique.values())

        # ----------------------------------------------------
        # Finalists:
        #
        # top few calibration candidates from every depth.
        # ----------------------------------------------------

        by_depth: Dict[int, List[Candidate]] = {}

        for candidate in candidates:

            depth = len(candidate.primes)

            by_depth.setdefault(
                depth,
                [],
            ).append(candidate)

        finalists: Dict[Tuple[int, ...], Candidate] = {}

        for depth, arr in by_depth.items():

            arr.sort(
                key=lambda c: (
                    c.calibration_survivors,
                    c.calibration_time,
                )
            )

            for candidate in arr[:FINALISTS_PER_DEPTH]:

                finalists[candidate.primes] = candidate

        # Also include the globally best calibration candidate.
        if candidates:

            best_calibration = min(
                candidates,
                key=lambda c: (
                    c.calibration_survivors,
                    c.calibration_time,
                )
            )

            finalists[
                best_calibration.primes
            ] = best_calibration

        # ----------------------------------------------------
        # Full-range validation.
        # ----------------------------------------------------

        print()
        print("-" * 78)
        print("FULL-RANGE VALIDATION")
        print("-" * 78)

        results = []

        for candidate in finalists.values():

            masks = [
                mask_cache[p]
                for p in candidate.primes
            ]

            elapsed, survivors = full_benchmark(
                start_x,
                true_x,
                masks,
            )

            target_survives = survives_at_x(
                true_x,
                masks,
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
                f"{str(list(candidate.primes)):<45} "
                f"full={elapsed:.7f}s "
                f"survivors={survivors:>9d} "
                f"target_survives={target_survives}"
            )

        # ----------------------------------------------------
        # Sort by actual runtime.
        # ----------------------------------------------------

        results.sort(
            key=lambda r: r[0]
        )

        print()
        print("-" * 78)
        print("BEST FULL-RANGE CANDIDATES")
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
                f"{str(list(candidate.primes)):<45} "
                f"{elapsed:.7f}s "
                f"survivors={survivors:>9d} "
                f"target={target_survives}"
            )

        # ----------------------------------------------------
        # Best valid target-preserving candidate.
        # ----------------------------------------------------

        valid = [
            r
            for r in results
            if r[3]
        ]

        print()

        if valid:

            best = valid[0]

            print("BEST TARGET-PRESERVING SUBSET")
            print("  primes    =", list(best[1].primes))
            print("  runtime   =", best[0])
            print("  survivors =", best[2])

        else:

            print(
                "WARNING: no finalist preserved the true Fermat x."
            )

        # ----------------------------------------------------
        # Baseline: no QR sieve.
        # ----------------------------------------------------

        baseline_masks: List[PrimeMask] = []

        baseline_time, baseline_survivors = full_benchmark(
            start_x,
            true_x,
            baseline_masks,
        )

        print()
        print("NO-SIEVE BASELINE")
        print("  runtime   =", baseline_time)
        print("  survivors =", baseline_survivors)

        if valid:

            speedup = (
                baseline_time / valid[0][0]
            )

            print(
                "  best speedup =",
                speedup,
                "x"
            )

    print()
    print("=" * 78)
    print("FINISHED EXPERIMENT 197")
    print("=" * 78)


if __name__ == "__main__":
    main()
