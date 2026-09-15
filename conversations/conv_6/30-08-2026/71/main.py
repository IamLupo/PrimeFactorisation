#!/usr/bin/env python3

import math
import random
import time
from collections import Counter

import sympy as sp


# =============================================================================
# START EXPERIMENT 71
# =============================================================================

print("=" * 100)
print("START EXPERIMENT 71")
print("RECURSIVE Kx DEBUG / STAGE-BY-STAGE BOTTLENECK EXPERIMENT")
print("=" * 100)


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 1_511_464_998

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 3

K_TARGET = 1000

R_OFFSETS = (
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
)

S_VALUES = (
    30,
    210,
    2310,
    30030,
)

PAIR_WINDOW = 3000

# Number of debug Kx values printed per R.
MAX_KX_DEBUG = 30

# Number of quotient candidates printed per R.
MAX_K_DEBUG = 30


# =============================================================================
# PRIME UTILITIES
# =============================================================================

def next_prime(x):
    return int(
        sp.nextprime(
            max(2, x - 1)
        )
    )


def choose_balanced_prime_pair(target_R):
    center = math.isqrt(
        max(target_R, 4)
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

                lo = min(
                    r1,
                    r2,
                )

                hi = max(
                    r1,
                    r2,
                )

                if hi / lo > 1.35:
                    continue

                R = r1 * r2

                score = (
                    abs(
                        R - target_R
                    ),
                    abs(
                        r1 - r2
                    ),
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
            f"Could not construct prime pair near R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# ANCHORS
# =============================================================================

def generate_anchor(
    scale,
    rng,
):
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
# FACTORIZATION + DIVISOR PAIRS
# =============================================================================

def factorint_timed(n):
    t0 = time.perf_counter()

    factors = dict(
        sp.factorint(n)
    )

    elapsed = (
        time.perf_counter()
        - t0
    )

    return factors, elapsed


def divisor_pairs_from_factors(
    factors
):
    t0 = time.perf_counter()

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)

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

    for prime, exponent in factors.items():
        N *= prime ** exponent

    root = math.isqrt(N)

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

    elapsed = (
        time.perf_counter()
        - t0
    )

    return pairs, elapsed


# =============================================================================
# RECONSTRUCTION
# =============================================================================

def quotient_residue_interval(
    quotient,
    radix,
    carry,
):
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


def test_quotient_candidate(
    n,
    r1,
    r2,
    k,
    l,
):
    """
    Exact reconstruction test.

    Returns:
        (p,q,tests)
    or:
        (None,None,tests)
    """

    if k <= 0 or l <= 0:
        return None, None, 0

    R = r1 * r2

    T = n // R

    K = k * l

    E = T - K

    if E < 0:
        return None, None, 0

    # Very cheap necessary condition.
    if E > k + l:
        return None, None, 0

    tests = 0

    for c3 in (
        0,
        1,
        2,
    ):

        remaining = E - c3

        if remaining < 0:
            continue

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

            c2 = remaining - c1

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

            aw = (
                a_hi
                - a_lo
                + 1
            )

            bw = (
                b_hi
                - b_lo
                + 1
            )

            if aw <= bw:

                for a in range(
                    a_lo,
                    a_hi + 1,
                ):

                    tests += 1

                    p = (
                        k * r1
                        + a
                    )

                    if p <= 1:
                        continue

                    if n % p != 0:
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

                    return p, q, tests

            else:

                for b in range(
                    b_lo,
                    b_hi + 1,
                ):

                    tests += 1

                    q = (
                        l * r2
                        + b
                    )

                    if q <= 1:
                        continue

                    if n % q != 0:
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

                    return p, q, tests

    return None, None, tests


# =============================================================================
# ONE AUXILIARY VALUE
# =============================================================================

def process_auxiliary(
    n,
    r1,
    r2,
    K_true,
    S,
):
    """
    Full timing breakdown for one n+x.
    """

    total_t0 = time.perf_counter()

    x = (-n) % S

    if x == 0:
        return {
            "S": S,
            "x": x,
            "nx": n,
            "factor_time": 0.0,
            "pair_time": 0.0,
            "Kx_factor_time": 0.0,
            "pairs": [],
            "Kx_records": [],
            "all_Kx": [],
            "total_time": 0.0,
        }

    nx = n + x

    factors, factor_time = (
        factorint_timed(
            nx
        )
    )

    pairs, pair_time = (
        divisor_pairs_from_factors(
            factors
        )
    )

    Kx_records = []
    seen_Kx = set()

    Kx_factor_time = 0.0

    for px, qx in pairs:

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

            t0 = time.perf_counter()

            Kx_factors = dict(
                sp.factorint(
                    Kx
                )
            )

            Kx_factor_time += (
                time.perf_counter()
                - t0
            )

            Kx_records.append(
                {
                    "Kx": Kx,
                    "kx": kx,
                    "lx": lx,
                    "factors": Kx_factors,
                    "match": (
                        Kx == K_true
                    ),
                }
            )

    total_time = (
        time.perf_counter()
        - total_t0
    )

    return {
        "S": S,
        "x": x,
        "nx": nx,
        "factors": factors,
        "factor_time": factor_time,
        "pair_time": pair_time,
        "pairs": pairs,
        "Kx_factor_time":
            Kx_factor_time,
        "Kx_records":
            Kx_records,
        "all_Kx":
            sorted(seen_Kx),
        "total_time":
            total_time,
    }


# =============================================================================
# ONE R
# =============================================================================

def run_R(
    n,
    p_true,
    q_true,
    r1,
    r2,
):
    total_t0 = time.perf_counter()

    R = r1 * r2

    T = n // R

    k_true = p_true // r1
    l_true = q_true // r2

    K_true = k_true * l_true

    E_true = T - K_true

    aux_results = []

    stage_factor = 0.0
    stage_pairs = 0.0
    stage_Kx_factor = 0.0

    all_Kx = set()

    for S in S_VALUES:

        aux = process_auxiliary(
            n,
            r1,
            r2,
            K_true,
            S,
        )

        aux_results.append(
            aux
        )

        stage_factor += (
            aux["factor_time"]
        )

        stage_pairs += (
            aux["pair_time"]
        )

        stage_Kx_factor += (
            aux["Kx_factor_time"]
        )

        all_Kx.update(
            aux["all_Kx"]
        )

    # Only Kx values themselves are candidate K values.
    candidate_Ks = sorted(
        all_Kx
    )

    # Keep testing statistics separate.
    candidate_results = []

    stage_reconstruct = 0.0

    total_residue_tests = 0

    solved = False
    solution = None

    tested_K = set()

    for K_candidate in candidate_Ks:

        if K_candidate <= 0:
            continue

        if K_candidate > T:
            continue

        if K_candidate in tested_K:
            continue

        tested_K.add(
            K_candidate
        )

        # Factor the small K separately for debug.
        t0 = time.perf_counter()

        K_factors = dict(
            sp.factorint(
                K_candidate
            )
        )

        factor_small_K_time = (
            time.perf_counter()
            - t0
        )

        K_pairs = (
            divisor_pairs_from_factors(
                K_factors
            )[0]
        )

        found_for_K = False

        for k, l in K_pairs:

            for kk, ll in (
                (k, l),
                (l, k),
            ):

                t0 = time.perf_counter()

                rp, rq, tests = (
                    test_quotient_candidate(
                        n,
                        r1,
                        r2,
                        kk,
                        ll,
                    )
                )

                reconstruct_time = (
                    time.perf_counter()
                    - t0
                )

                stage_reconstruct += (
                    reconstruct_time
                )

                total_residue_tests += (
                    tests
                )

                if rp is not None:

                    solution = {
                        "p": rp,
                        "q": rq,
                        "k": kk,
                        "l": ll,
                        "K": K_candidate,
                        "K_factors":
                            K_factors,
                    }

                    solved = True
                    found_for_K = True

                    break

            if found_for_K:
                break

        candidate_results.append(
            {
                "K": K_candidate,
                "matches_true":
                    K_candidate == K_true,
                "factors":
                    K_factors,
                "pairs":
                    K_pairs,
                "solved":
                    found_for_K,
                "small_factor_time":
                    factor_small_K_time,
            }
        )

        if solved:
            break

    total_time = (
        time.perf_counter()
        - total_t0
    )

    return {
        "R": R,
        "r1": r1,
        "r2": r2,
        "T": T,
        "E": E_true,
        "k": k_true,
        "l": l_true,
        "K": K_true,
        "aux_results":
            aux_results,
        "all_Kx":
            sorted(all_Kx),
        "candidate_results":
            candidate_results,
        "factor_time":
            stage_factor,
        "pair_time":
            stage_pairs,
        "Kx_factor_time":
            stage_Kx_factor,
        "reconstruct_time":
            stage_reconstruct,
        "residue_tests":
            total_residue_tests,
        "solution":
            solution,
        "total_time":
            total_time,
    }


# =============================================================================
# SCALE RUN
# =============================================================================

def run_scale(
    scale,
    rng,
):
    print()
    print("=" * 100)
    print(
        f"SCALE {scale:.0e}"
    )
    print("=" * 100)

    rows = []

    for index in range(
        ANCHORS_PER_SCALE
    ):

        p_true, q_true, n = (
            generate_anchor(
                scale,
                rng,
            )
        )

        print()
        print(
            f"anchor {index+1}/"
            f"{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        # Reference factorization is NOT used by the solver.
        # It is retained here solely for validation/debug.
        print(
            f"    TRUE p,q = "
            f"({p_true:,},"
            f"{q_true:,})"
        )

        base_R = (
            n // K_TARGET
        )

        regimes = []

        for offset in R_OFFSETS:

            target_R = int(
                round(
                    base_R * offset
                )
            )

            t0 = time.perf_counter()

            r1, r2, R = (
                choose_balanced_prime_pair(
                    target_R
                )
            )

            selection_time = (
                time.perf_counter()
                - t0
            )

            result = run_R(
                n,
                p_true,
                q_true,
                r1,
                r2,
            )

            result[
                "selection_time"
            ] = selection_time

            regimes.append(
                result
            )

        # Show every attempted R, but compactly.
        for r in regimes:

            R = r["R"]

            print()
            print(
                f"    R={R:,} "
                f"(r1,r2)=("
                f"{r['r1']},"
                f"{r['r2']})"
            )

            print(
                f"        TRUE: "
                f"k={r['k']} "
                f"l={r['l']} "
                f"K={r['K']:,} "
                f"T={r['T']:,} "
                f"E={r['E']:,}"
            )

            print(
                f"        times: "
                f"Rselect="
                f"{r['selection_time']:.4f}s "
                f"factor(n+x)="
                f"{r['factor_time']:.4f}s "
                f"divisors="
                f"{r['pair_time']:.4f}s "
                f"factor(Kx)="
                f"{r['Kx_factor_time']:.4f}s "
                f"reconstruct="
                f"{r['reconstruct_time']:.4f}s "
                f"TOTAL="
                f"{(
                    r['selection_time']
                    + r['total_time']
                ):.4f}s"
            )

            print(
                f"        Kx distinct="
                f"{len(r['all_Kx'])} "
                f"true-K found="
                f"{r['K'] in r['all_Kx']}"
            )

            # Explicit Kx debug.
            shown = 0

            for aux in r[
                "aux_results"
            ]:

                matches = [
                    rec
                    for rec in aux[
                        "Kx_records"
                    ]
                    if rec["match"]
                ]

                print(
                    f"        S={aux['S']} "
                    f"x={aux['x']:,} "
                    f"n+x={aux['nx']:,} "
                    f"divPairs="
                    f"{len(aux['pairs'])} "
                    f"Kx="
                    f"{len(aux['Kx_records'])} "
                    f"trueMatches="
                    f"{len(matches)}"
                )

                shown_here = 0

                for rec in aux[
                    "Kx_records"
                ]:

                    if shown >= MAX_KX_DEBUG:
                        break

                    print(
                        f"            "
                        f"Kx={rec['Kx']:,} "
                        f"(kx,lx)=("
                        f"{rec['kx']},"
                        f"{rec['lx']}) "
                        f"{'MATCH' if rec['match'] else ''}"
                    )

                    shown += 1
                    shown_here += 1

                if (
                    len(aux["Kx_records"])
                    > shown_here
                ):
                    print(
                        "            ..."
                    )

                if shown >= MAX_KX_DEBUG:
                    break

            print(
                "        K CANDIDATES TESTED:"
            )

            shown_K = 0

            for candidate in r[
                "candidate_results"
            ]:

                print(
                    f"            "
                    f"K={candidate['K']:,} "
                    f"{'TRUE-K' if candidate['matches_true'] else ''} "
                    f"{'SOLVED' if candidate['solved'] else 'rejected'}"
                )

                shown_K += 1

                if shown_K >= MAX_K_DEBUG:
                    print(
                        "            ..."
                    )
                    break

            if r["solution"]:

                s = r[
                    "solution"
                ]

                print(
                    f"        SOLUTION: "
                    f"({s['p']:,},"
                    f"{s['q']:,}) "
                    f"(k,l)=("
                    f"{s['k']},"
                    f"{s['l']}) "
                    f"K={s['K']:,}"
                )

            rows.append(
                {
                    "scale": scale,
                    "n": n,
                    "p_true": p_true,
                    "q_true": q_true,
                    **r,
                }
            )

    return rows


# =============================================================================
# BOTTLENECK SUMMARY
# =============================================================================

def print_bottleneck_summary(
    rows,
):
    print()
    print("=" * 100)
    print(
        "BOTTLENECK SUMMARY"
    )
    print("=" * 100)

    print(
        "scale       cases   "
        "factor(n+x)   divisors   "
        "factor(Kx)   reconstruct   "
        "other"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in rows
            if r["scale"] == scale
        ]

        if not group:
            continue

        factor_time = sum(
            r["factor_time"]
            for r in group
        )

        pair_time = sum(
            r["pair_time"]
            for r in group
        )

        Kx_time = sum(
            r["Kx_factor_time"]
            for r in group
        )

        reconstruct = sum(
            r["reconstruct_time"]
            for r in group
        )

        total = sum(
            r["selection_time"]
            + r["total_time"]
            for r in group
        )

        accounted = (
            factor_time
            + pair_time
            + Kx_time
            + reconstruct
        )

        other = max(
            0.0,
            total - accounted,
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):8d} "
            f"{factor_time:14.4f}s "
            f"{pair_time:11.4f}s "
            f"{Kx_time:12.4f}s "
            f"{reconstruct:14.4f}s "
            f"{other:10.4f}s"
        )


# =============================================================================
# FINAL SCALE SUMMARY
# =============================================================================

def print_final_summary(
    rows,
):
    print()
    print("=" * 100)
    print(
        "END-TO-END RESULT"
    )
    print("=" * 100)

    print(
        "scale       cases   solved   "
        "true-K-found   avg R-tests   "
        "avg Kx-count   avg time"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in rows
            if r["scale"] == scale
        ]

        solved = sum(
            r["solution"] is not None
            for r in group
        )

        true_K_found = sum(
            r["K"] in r["all_Kx"]
            for r in group
        )

        # Number of R values tried per anchor is fixed here.
        anchors = ANCHORS_PER_SCALE

        avg_R = (
            len(group) / anchors
        )

        avg_Kx = (
            sum(
                len(r["all_Kx"])
                for r in group
            )
            / len(group)
        )

        avg_time = (
            sum(
                r["selection_time"]
                + r["total_time"]
                for r in group
            )
            / len(group)
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):8d} "
            f"{solved:8d} "
            f"{true_K_found:13d} "
            f"{avg_R:13.2f} "
            f"{avg_Kx:14.2f} "
            f"{avg_time:10.4f}s"
        )


# =============================================================================
# INTERPRETATION
# =============================================================================

def print_interpretation():

    print()
    print("=" * 100)
    print(
        "MATHEMATICAL / COMPUTATIONAL INTERPRETATION"
    )
    print("=" * 100)

    print(
        r"""
Experiment 71 changes one important thing:

it does not merely ask whether the recursive solver succeeds.

It records the actual Kx information obtained before reconstruction.

For every R:

    K_true = floor(p/r1) * floor(q/r2)

is shown only as validation/debug information.

The solver itself obtains:

    n+x
      |
      v
    factor(n+x)
      |
      v
    divisor pairs
      |
      v
    (kx,lx)
      |
      v
    Kx = kx*lx.

Every distinct Kx is explicitly recorded.

The experiment then asks:

    Is K_true present among the Kx values?

This separates two failure modes.

A:

    K_true NOT FOUND

means the auxiliary trajectory did not expose the required
recursive object.

B:

    K_true FOUND
    but reconstruction fails

means the Kx bridge works, but the conversion:

    K -> (k,l) -> residue cells -> p,q

still needs improvement.

The timing instrumentation separates:

    factor(n+x)
    divisor enumeration
    factor(Kx)
    reconstruction.

The main expected bottleneck for large n is factor(n+x).

The Kx factorization is expected to be tiny because K is deliberately
kept near K_TARGET=1000.

If the timing confirms this, future optimization should target the
construction/factorization of n+x rather than the recursive K layer.

The experiment also records every Kx candidate so that a failed
10^12 or 10^16 case can be diagnosed directly.
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

    all_rows = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            rng,
        )

        all_rows.extend(
            rows
        )

    print_bottleneck_summary(
        all_rows
    )

    print_final_summary(
        all_rows
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
    print("FINISHED EXPERIMENT 71")
    print("=" * 100)


if __name__ == "__main__":
    main()
