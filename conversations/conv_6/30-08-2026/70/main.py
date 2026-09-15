#!/usr/bin/env python3

import math
import random
import time
from collections import Counter

import sympy as sp


# =============================================================================
# START EXPERIMENT 70
# =============================================================================

print("=" * 100)
print("START EXPERIMENT 70")
print("RECURSIVE Kx -> K -> (k,l) END-TO-END FACTORING TEST")
print("=" * 100)


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 1_511_464_998

# Main scaling cases.
SCALES = (
    10**9,
    10**12,
    10**16,
)

# Keep this modest initially.
ANCHORS_PER_SCALE = 3

# Desired recursive object.
K_TARGET = 1000

# Try several R positions around n / K_TARGET.
R_OFFSETS = (
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
)

# Deterministic auxiliary constructions.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Search distance around sqrt(R) when finding balanced prime moduli.
PAIR_WINDOW = 3000

# Maximum number of auxiliary divisor pairs inspected.
MAX_AUX_PAIRS = 100_000

# Maximum number of K candidates generated from one auxiliary trajectory.
MAX_K_CANDIDATES_PER_AUX = 500

# Maximum residue states tested for one (k,l) candidate.
MAX_RESIDUE_TESTS = 5_000_000

# Number of output examples.
MAX_EXAMPLES = 8


# =============================================================================
# PRIME UTILITIES
# =============================================================================

def next_prime(x: int) -> int:
    return int(
        sp.nextprime(
            max(2, x - 1)
        )
    )


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchor(
    scale: int,
    rng: random.Random,
):
    """
    Generate:

        n = p*q

    with p and q balanced around sqrt(scale).

    p and q are retained ONLY for final validation.
    The solver itself never receives them.
    """

    root = math.isqrt(scale)

    lo = max(
        101,
        int(root * 0.75),
    )

    hi = max(
        lo + 100,
        int(root * 1.25),
    )

    p = int(
        sp.randprime(
            lo,
            hi,
        )
    )

    q = int(
        sp.randprime(
            lo,
            hi,
        )
    )

    while q == p:
        q = int(
            sp.randprime(
                lo,
                hi,
            )
        )

    return p, q, p * q


# =============================================================================
# BALANCED PRIME MODULUS PAIR
# =============================================================================

def choose_balanced_prime_pair(
    target_R: int,
):
    """
    Find prime r1,r2 such that:

        r1*r2 ~= target_R

    while keeping them reasonably balanced.
    """

    target_R = max(
        target_R,
        4,
    )

    center = math.isqrt(
        target_R
    )

    best = None

    for delta in range(
        PAIR_WINDOW + 1
    ):

        for c in (
            center - delta,
            center + delta,
        ):

            if c < 2:
                continue

            r1 = next_prime(c)

            desired = max(
                2,
                target_R // r1,
            )

            for probe in (
                desired - 2,
                desired - 1,
                desired,
                desired + 1,
                desired + 2,
            ):

                if probe < 2:
                    continue

                r2 = next_prime(
                    probe
                )

                lo = min(r1, r2)
                hi = max(r1, r2)

                if (
                    hi / lo
                    > 1.35
                ):
                    continue

                R = r1 * r2

                score = (
                    abs(R - target_R),
                    abs(r1 - r2),
                )

                if (
                    best is None
                    or score < best[0]
                ):
                    best = (
                        score,
                        r1,
                        r2,
                        R,
                    )

        if (
            best is not None
            and best[0][0] == 0
        ):
            break

    if best is None:
        raise RuntimeError(
            f"Could not find balanced prime pair "
            f"near R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# FACTOR / DIVISOR UTILITIES
# =============================================================================

def factorint_exact(
    n: int,
):
    return dict(
        sp.factorint(n)
    )


def divisor_pairs_from_factors(
    factors,
):
    """
    Generate all unordered divisor pairs d * e = N.
    """

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(
            divisors
        )

        powers = [1]
        cur = 1

        for _ in range(exponent):
            cur *= prime
            powers.append(cur)

        divisors = [
            d * p
            for d in old
            for p in powers
        ]

    divisors = sorted(
        set(divisors)
    )

    N = 1

    for p, e in factors.items():
        N *= p ** e

    root = math.isqrt(
        N
    )

    pairs = []

    for d in divisors:

        if d > root:
            break

        if N % d == 0:
            pairs.append(
                (
                    d,
                    N // d,
                )
            )

    return pairs


def factor_value(
    n: int,
):
    factors = factorint_exact(
        n
    )

    pairs = divisor_pairs_from_factors(
        factors
    )

    return factors, pairs


# =============================================================================
# DETERMINISTIC n+x
# =============================================================================

def deterministic_x(
    n: int,
    S: int,
):
    """
    x = (-n) mod S

    therefore:

        n+x == 0 mod S
    """

    return (-n) % S


# =============================================================================
# RESIDUE INTERVAL
# =============================================================================

def quotient_residue_interval(
    quotient: int,
    radix: int,
    carry: int,
):
    """
    From

        carry = floor(residue * quotient / radix)

    derive all residue values producing that carry.
    """

    # carry <= residue*q/radix < carry+1

    lo = (
        carry * radix
        + quotient
        - 1
    ) // quotient

    hi = (
        ((carry + 1) * radix) - 1
    ) // quotient

    lo = max(
        0,
        lo,
    )

    hi = min(
        radix - 1,
        hi,
    )

    return lo, hi


# =============================================================================
# TEST ONE (k,l) CANDIDATE
# =============================================================================

def test_quotient_candidate(
    n: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
):
    """
    Test one candidate:

        p = k*r1 + a
        q = l*r2 + b

    using the carry structure.

    Returns:

        (p,q)
        or None

    No factorization of n is performed here.
    """

    if k <= 0 or l <= 0:
        return None

    R = r1 * r2

    T = n // R

    K = k * l

    E = T - K

    # Fundamental sanity conditions.
    if E < 0:
        return None

    # From:

    #   E = c1 + c2 + c3

    # with:

    #   c1 <= l-1
    #   c2 <= k-1
    #   c3 <= 2

    # therefore:

    if E > k + l:
        return None

    residue_tests = 0

    # c3 is tiny.
    for c3 in (
        0,
        1,
        2,
    ):

        remaining = E - c3

        if remaining < 0:
            continue

        # c1 + c2 = remaining
        #
        # c1 < l
        # c2 < k

        c1_lo = max(
            0,
            remaining - (k - 1),
        )

        c1_hi = min(
            l - 1,
            remaining,
        )

        for c1 in range(
            c1_lo,
            c1_hi + 1,
        ):

            c2 = (
                remaining - c1
            )

            if not (
                0 <= c2 < k
            ):
                continue

            a_lo, a_hi = (
                quotient_residue_interval(
                    l,
                    r1,
                    c1,
                )
            )

            b_lo, b_hi = (
                quotient_residue_interval(
                    k,
                    r2,
                    c2,
                )
            )

            if (
                a_lo > a_hi
                or b_lo > b_hi
            ):
                continue

            a_width = (
                a_hi - a_lo + 1
            )

            b_width = (
                b_hi - b_lo + 1
            )

            # Scan the smaller dimension.
            if a_width <= b_width:

                if (
                    residue_tests
                    + a_width
                    > MAX_RESIDUE_TESTS
                ):
                    return None

                for a in range(
                    a_lo,
                    a_hi + 1,
                ):

                    residue_tests += 1

                    p = (
                        k * r1
                        + a
                    )

                    if p <= 1:
                        continue

                    if n % p:
                        continue

                    q = n // p

                    if not (
                        l * r2
                        <= q
                        < (l + 1) * r2
                    ):
                        continue

                    b = (
                        q
                        - l * r2
                    )

                    if not (
                        b_lo
                        <= b
                        <= b_hi
                    ):
                        continue

                    return p, q

            else:

                if (
                    residue_tests
                    + b_width
                    > MAX_RESIDUE_TESTS
                ):
                    return None

                for b in range(
                    b_lo,
                    b_hi + 1,
                ):

                    residue_tests += 1

                    q = (
                        l * r2
                        + b
                    )

                    if q <= 1:
                        continue

                    if n % q:
                        continue

                    p = n // q

                    if not (
                        k * r1
                        <= p
                        < (k + 1) * r1
                    ):
                        continue

                    a = (
                        p
                        - k * r1
                    )

                    if not (
                        a_lo
                        <= a
                        <= a_hi
                    ):
                        continue

                    return p, q

    return None


# =============================================================================
# EXTRACT K CANDIDATES FROM AUXILIARY FACTORIZATION
# =============================================================================

def derive_K_candidates(
    nx: int,
    r1: int,
    r2: int,
):
    """
    Factor n+x and transform every divisor pair into Kx.

    Then factor each small Kx and return all candidate K values
    encountered through its divisor structure.
    """

    nx_factors = factorint_exact(
        nx
    )

    divisor_pairs = (
        divisor_pairs_from_factors(
            nx_factors
        )
    )

    results = []

    seen_Kx = set()
    seen_K = set()

    for px, qx in divisor_pairs:

        if len(results) >= MAX_AUX_PAIRS:
            break

        for a, b in (
            (px, qx),
            (qx, px),
        ):

            kx = a // r1
            lx = b // r2

            if kx <= 0 or lx <= 0:
                continue

            Kx = kx * lx

            if Kx in seen_Kx:
                continue

            seen_Kx.add(
                Kx
            )

            if Kx <= 0:
                continue

            Kx_factors = factorint_exact(
                Kx
            )

            Kx_pairs = (
                divisor_pairs_from_factors(
                    Kx_factors
                )
            )

            for u, v in Kx_pairs:

                for K_candidate in (
                    u * v,
                ):

                    if K_candidate in seen_K:
                        continue

                    seen_K.add(
                        K_candidate
                    )

                    results.append(
                        {
                            "Kx": Kx,
                            "kx": kx,
                            "lx": lx,
                            "K_factors":
                                Kx_factors,
                            "candidate_K":
                                K_candidate,
                        }
                    )

                    if (
                        len(results)
                        >= MAX_K_CANDIDATES_PER_AUX
                    ):
                        break

                if (
                    len(results)
                    >= MAX_K_CANDIDATES_PER_AUX
                ):
                    break

            if (
                len(results)
                >= MAX_K_CANDIDATES_PER_AUX
            ):
                break

    return (
        nx_factors,
        divisor_pairs,
        results,
    )


# =============================================================================
# ONE R REGIME
# =============================================================================

def run_R_regime(
    n: int,
    r1: int,
    r2: int,
    invalid_K: set,
):
    """
    Actual recursive solver at one modulus product R.
    """

    R = r1 * r2

    T = n // R

    tested_K = set()

    aux_count = 0
    Kx_values = set()
    K_candidates = set()

    for S in S_VALUES:

        x = deterministic_x(
            n,
            S,
        )

        if x == 0:
            continue

        nx = n + x

        (
            nx_factors,
            aux_pairs,
            derived,
        ) = derive_K_candidates(
            nx,
            r1,
            r2,
        )

        aux_count += len(
            aux_pairs
        )

        for item in derived:

            Kx_values.add(
                item["Kx"]
            )

            K_candidate = item[
                "candidate_K"
            ]

            if K_candidate <= 0:
                continue

            if K_candidate > T:
                continue

            if K_candidate in invalid_K:
                continue

            if K_candidate in tested_K:
                continue

            tested_K.add(
                K_candidate
            )

            K_candidates.add(
                K_candidate
            )

            k_pairs = factor_value(
                K_candidate
            )[1]

            for k, l in k_pairs:

                for kk, ll in (
                    (k, l),
                    (l, k),
                ):

                    key = (
                        kk,
                        ll,
                    )

                    if key in tested_K:
                        pass

                    result = (
                        test_quotient_candidate(
                            n,
                            r1,
                            r2,
                            kk,
                            ll,
                        )
                    )

                    if result is not None:

                        return {
                            "success": True,
                            "p": result[0],
                            "q": result[1],
                            "R": R,
                            "r1": r1,
                            "r2": r2,
                            "T": T,
                            "K": K_candidate,
                            "k": kk,
                            "l": ll,
                            "aux_count":
                                aux_count,
                            "Kx_count":
                                len(Kx_values),
                            "K_candidates":
                                len(K_candidates),
                        }

            # If this K has produced no factor pair capable of
            # recovering n, permanently reject it for this n.
            invalid_K.add(
                K_candidate
            )

    return {
        "success": False,
        "R": R,
        "r1": r1,
        "r2": r2,
        "T": T,
        "aux_count": aux_count,
        "Kx_count": len(Kx_values),
        "K_candidates": len(K_candidates),
    }


# =============================================================================
# COMPLETE SOLVER
# =============================================================================

def solve_n(
    n: int,
):
    """
    Try to factor n without giving the solver the true p,q.

    The solver repeatedly changes R.

        n
         |
         v
        R
         |
         v
        n+x
         |
         v
        factor(n+x)
         |
         v
        Kx
         |
         v
        factor(Kx)
         |
         v
        candidate K
         |
         v
        candidate (k,l)
         |
         v
        carry cell
         |
         v
        exact pq=n
    """

    invalid_K = set()

    total_aux = 0
    total_K_tests = 0
    total_R = 0

    tried_R = set()

    # The target is R ~= n / 1000.
    base_R = (
        n // K_TARGET
    )

    for offset in R_OFFSETS:

        target_R = int(
            round(
                base_R * offset
            )
        )

        r1, r2, R = (
            choose_balanced_prime_pair(
                target_R
            )
        )

        if R in tried_R:
            continue

        tried_R.add(
            R
        )

        total_R += 1

        result = run_R_regime(
            n,
            r1,
            r2,
            invalid_K,
        )

        total_aux += result[
            "aux_count"
        ]

        total_K_tests += result[
            "K_candidates"
        ]

        if result[
            "success"
        ]:

            result[
                "total_R_tests"
            ] = total_R

            result[
                "total_aux"
            ] = total_aux

            result[
                "total_K_tests"
            ] = total_K_tests

            result[
                "invalid_K"
            ] = len(
                invalid_K
            )

            return result

    return {
        "success": False,
        "total_R_tests":
            total_R,
        "total_aux":
            total_aux,
        "total_K_tests":
            total_K_tests,
        "invalid_K":
            len(invalid_K),
    }


# =============================================================================
# ONE SCALE
# =============================================================================

def run_scale(
    scale: int,
    rng: random.Random,
):
    print()
    print("=" * 100)
    print(
        f"SCALE {scale:.0e}"
    )
    print("=" * 100)

    results = []

    for i in range(
        ANCHORS_PER_SCALE
    ):

        p_true, q_true, n = (
            generate_anchor(
                scale,
                rng,
            )
        )

        t0 = time.perf_counter()

        solution = solve_n(
            n
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        # Validation only.
        recovered = (
            solution.get(
                "success",
                False,
            )
            and (
                solution["p"]
                * solution["q"]
                == n
            )
        )

        same_orientation = (
            recovered
            and (
                (
                    solution["p"]
                    == p_true
                    and solution["q"]
                    == q_true
                )
                or
                (
                    solution["p"]
                    == q_true
                    and solution["q"]
                    == p_true
                )
            )
        )

        row = {
            "scale": scale,
            "n": n,
            "true_p": p_true,
            "true_q": q_true,
            "solution": solution,
            "time": elapsed,
            "recovered":
                recovered,
            "same_orientation":
                same_orientation,
        }

        results.append(
            row
        )

        if solution.get(
            "success",
            False,
        ):

            print(
                f"anchor {i+1:2d}/"
                f"{ANCHORS_PER_SCALE} "
                f"n={n:,}"
            )

            print(
                f"    SOLVED "
                f"R={solution['R']:,} "
                f"(r1,r2)=("
                f"{solution['r1']},"
                f"{solution['r2']})"
            )

            print(
                f"    K={solution['K']:,} "
                f"(k,l)=("
                f"{solution['k']},"
                f"{solution['l']}) "
                f"T={solution['T']:,}"
            )

            print(
                f"    recovered=("
                f"{solution['p']:,},"
                f"{solution['q']:,}) "
                f"valid={recovered}"
            )

            print(
                f"    R-tests="
                f"{solution['total_R_tests']} "
                f"aux-tests="
                f"{solution['total_aux']} "
                f"K-tests="
                f"{solution['total_K_tests']} "
                f"time="
                f"{elapsed:.4f}s"
            )

        else:

            print(
                f"anchor {i+1:2d}/"
                f"{ANCHORS_PER_SCALE} "
                f"n={n:,} "
                f"NOT SOLVED "
                f"R-tests="
                f"{solution.get('total_R_tests',0)} "
                f"aux="
                f"{solution.get('total_aux',0)} "
                f"K-tests="
                f"{solution.get('total_K_tests',0)} "
                f"time="
                f"{elapsed:.4f}s"
            )

    return results


# =============================================================================
# SUMMARY
# =============================================================================

def print_summary(
    all_results,
):
    print()
    print("=" * 100)
    print(
        "END-TO-END SOLVER SUMMARY"
    )
    print("=" * 100)

    print(
        "scale       cases   solved   "
        "valid   avg R-tests   "
        "avg aux   avg K-tests   avg time"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in all_results
            if r["scale"] == scale
        ]

        if not group:
            continue

        solved = sum(
            r["solution"].get(
                "success",
                False,
            )
            for r in group
        )

        valid = sum(
            r["recovered"]
            for r in group
        )

        avg_R = (
            sum(
                r["solution"].get(
                    "total_R_tests",
                    0,
                )
                for r in group
            )
            / len(group)
        )

        avg_aux = (
            sum(
                r["solution"].get(
                    "total_aux",
                    0,
                )
                for r in group
            )
            / len(group)
        )

        avg_K = (
            sum(
                r["solution"].get(
                    "total_K_tests",
                    0,
                )
                for r in group
            )
            / len(group)
        )

        avg_time = (
            sum(
                r["time"]
                for r in group
            )
            / len(group)
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):8d} "
            f"{solved:8d} "
            f"{valid:7d} "
            f"{avg_R:13.2f} "
            f"{avg_aux:10.2f} "
            f"{avg_K:13.2f} "
            f"{avg_time:10.4f}s"
        )


# =============================================================================
# SUCCESS DETAILS
# =============================================================================

def print_successes(
    all_results,
):
    print()
    print("=" * 100)
    print(
        "SUCCESSFUL RECOVERIES"
    )
    print("=" * 100)

    shown = 0

    for row in all_results:

        if not row["recovered"]:
            continue

        s = row[
            "solution"
        ]

        n = row[
            "n"
        ]

        print(
            f"scale={row['scale']:.0e} "
            f"n={n:,}"
        )

        print(
            f"    true=("
            f"{row['true_p']:,},"
            f"{row['true_q']:,})"
        )

        print(
            f"    recovered=("
            f"{s['p']:,},"
            f"{s['q']:,})"
        )

        print(
            f"    R={s['R']:,} "
            f"K={s['K']:,} "
            f"(k,l)=("
            f"{s['k']},"
            f"{s['l']})"
        )

        print(
            f"    R-tests="
            f"{s['total_R_tests']} "
            f"aux="
            f"{s['total_aux']} "
            f"K-tests="
            f"{s['total_K_tests']} "
            f"time="
            f"{row['time']:.4f}s"
        )

        print()

        shown += 1

        if shown >= MAX_EXAMPLES:
            break

    if shown == 0:
        print("none")


# =============================================================================
# INTERPRETATION
# =============================================================================

def print_interpretation():

    print()
    print("=" * 100)
    print(
        "INTERPRETATION"
    )
    print("=" * 100)

    print(
        r"""
This is the first experiment in this series that attempts the
complete proposed factoring mechanism.

The solver is given only:

    n

It does NOT receive:

    p
    q
    k
    l
    E.

For a selected modulus product:

    R = r1*r2

it computes:

    T = floor(n/R).

A deterministic auxiliary value is constructed:

    x = (-n) mod S

so that:

    n+x == 0 mod S.

The auxiliary integer is factored:

    n+x
        ->
    divisor pairs
        ->
    (kx,lx)
        ->
    Kx = kx*lx.

The small Kx is factored.

Its divisor structure generates candidate quotient products K.

For every candidate K:

    E = T-K.

The carry identity gives the necessary condition:

    E = c1+c2+c3

with:

    0 <= c1 < l
    0 <= c2 < k
    0 <= c3 <= 2.

Therefore:

    E <= k+l.

Candidates violating that condition are rejected immediately.

For every remaining divisor pair:

    K = k*l

the carry decomposition is used to produce residue cells:

    p = k*r1 + a
    q = l*r2 + b.

The residue intervals are derived from:

    c1 = floor(a*l/r1)
    c2 = floor(b*k/r2).

Only residue values inside those cells are tested.

The exact test is finally:

    p*q == n.

If a K candidate fails, it is added to an invalid-K set and is
not tested again.

If the current R produces no solution, the solver changes R and
repeats the complete recursive process.

The intended architecture is therefore:

    n
     |
     +--> choose R ~= n/1000
     |
     +--> deterministic n+x
     |        |
     |        +--> factor(n+x)
     |                 |
     |                 +--> Kx
     |                       |
     |                       +--> factor(Kx)
     |                               |
     |                               +--> candidate K
     |
     +----------------------------------------+
                                              |
                                              v
                                          E = T-K
                                              |
                                              v
                                          factor K
                                              |
                                              v
                                            (k,l)
                                              |
                                              v
                                        carry cells
                                              |
                                              v
                                        exact pq=n

The critical statistic is now:

    solved / tested.

More importantly, record:

    number of R values tested,
    number of auxiliary factorizations,
    number of distinct K candidates tested,
    number of invalid K values,
    residue-state work,
    total runtime.

This experiment therefore tests the actual recursive solver
hypothesis rather than merely measuring Kx collisions.

The factorization of n itself is NOT used by the solver.
The known p,q values are used only afterward to validate that a
reported solution is the original factorization.

SymPy is used for the auxiliary n+x factorizations and the small
K factorizations.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    rng = random.Random(
        SEED
    )

    all_results = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            rng,
        )

        all_results.extend(
            rows
        )

    print_summary(
        all_results
    )

    print_successes(
        all_results
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    print()
    print("=" * 100)
    print(
        "TOTAL TIMING"
    )
    print("=" * 100)

    print(
        f"total runtime = "
        f"{elapsed:.4f}s"
    )

    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 70")
    print("=" * 100)


if __name__ == "__main__":
    main()
