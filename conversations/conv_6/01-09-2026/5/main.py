#!/usr/bin/env python3

import math
import random
import statistics
import time

# ============================================================
# START EXPERIMENT 78
# solve_ab TIMING AT 1e9 AND 1e16
# ============================================================

EXPERIMENT_ID = 78

SCALES = [
    10**9,
    10**16,
]

CASES_PER_SCALE = 5

# Keep the same style of modulus construction used
# throughout the previous experiments.
R_OFFSETS = (0.95, 1.00, 1.05)

SEED = 1511464998

# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

def fmt(n):
    return f"{n:,}"


def scale_label(n):
    if n == 10**9:
        return "1e+09"
    if n == 10**16:
        return "1e+16"
    return f"{n:.0e}"


# ------------------------------------------------------------
# Prime generation
# ------------------------------------------------------------

def is_prime(n):
    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Deterministic Miller-Rabin for the ranges used here.
    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    witnesses = (
        2, 325, 9375, 28178,
        450775, 9780504, 1795265022
    )

    for a in witnesses:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def next_prime(n):
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# ------------------------------------------------------------
# Generate balanced semiprime
# ------------------------------------------------------------

def make_semiprime(scale, rng):
    """
    Generate p*q ~= scale with p and q near sqrt(scale).
    """

    root = math.isqrt(scale)

    # Small random offsets around sqrt(scale).
    spread = max(100, root // 8)

    p = next_prime(
        root + rng.randint(-spread, spread)
    )

    q = next_prime(
        root + rng.randint(-spread, spread)
    )

    # Avoid equal primes.
    if p == q:
        q = next_prime(q + 2)

    n = p * q

    return n, p, q


# ------------------------------------------------------------
# R construction
# ------------------------------------------------------------

def make_R(n, offset, rng):
    """
    Choose R near sqrt(n).

    We deliberately search for a factorization
        R = r1*r2
    with r1 and r2 near sqrt(R).
    """

    root_n = math.isqrt(n)

    target_R = max(2, int(root_n * offset))

    root_R = math.isqrt(target_R)

    # Search around sqrt(R) for nearby factors.
    candidates = []

    for delta1 in range(-200, 201):
        r1 = root_R + delta1

        if r1 < 2:
            continue

        # Search a small neighborhood for r2.
        approx_r2 = max(2, target_R // r1)

        for delta2 in range(-200, 201):
            r2 = approx_r2 + delta2

            if r2 < 2:
                continue

            R = r1 * r2

            candidates.append(
                (abs(R - target_R), R, r1, r2)
            )

    if not candidates:
        raise RuntimeError("Could not construct R")

    candidates.sort()

    # Pick among the best few to avoid always getting
    # exactly the same construction.
    best = candidates[:20]
    _, R, r1, r2 = rng.choice(best)

    return R, r1, r2


# ------------------------------------------------------------
# Build K/k/l/a/b representation
# ------------------------------------------------------------

def representation_from_n(n, p, q, R, r1, r2):
    """
    For the experiment we know the hidden factorization and
    construct its exact quotient-cell representation.

    This is only the test setup.

        p = r1*k + a
        q = r2*l + b

    K = k*l
    T = floor(n/R)
    E = T-K
    """

    k = p // r1
    a = p - r1 * k

    l = q // r2
    b = q - r2 * l

    K = k * l

    T = n // R
    E = T - K

    # Sanity checks.
    assert p == r1 * k + a
    assert q == r2 * l + b
    assert n == (r1 * k + a) * (r2 * l + b)
    assert R == r1 * r2
    assert E >= 0

    return K, T, E, k, l, a, b


# ------------------------------------------------------------
# solve_ab
# ------------------------------------------------------------

def solve_ab(n, r1, r2, K, E, k, l):
    """
    Recover a,b from:

        p = r1*k + a
        q = r2*l + b
        p*q = n

    Search domain:

        0 <= a < r1
        0 <= b < r2

    For each candidate a:

        p = r1*k + a

    If p divides n, then:

        q = n/p

    and therefore:

        b = q - r2*l

    E is used as an exact consistency condition:

        E = floor(n/R) - K

    with:

        R = r1*r2
    """

    R = r1 * r2

    # E consistency check.
    expected_E = (n // R) - K

    if expected_E != E:
        return None

    base_p = r1 * k
    base_q = r2 * l

    tested = 0

    start = time.perf_counter()

    for a in range(r1):
        tested += 1

        p_candidate = base_p + a

        if p_candidate <= 1:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = n // p_candidate

        b = q_candidate - base_q

        if not (0 <= b < r2):
            continue

        # Final exact verification.
        if p_candidate * q_candidate != n:
            continue

        # Verify the original coordinates.
        if p_candidate != r1 * k + a:
            continue

        if q_candidate != r2 * l + b:
            continue

        elapsed = time.perf_counter() - start

        return {
            "a": a,
            "b": b,
            "p": p_candidate,
            "q": q_candidate,
            "tested_a": tested,
            "time": elapsed,
        }

    elapsed = time.perf_counter() - start

    return {
        "a": None,
        "b": None,
        "p": None,
        "q": None,
        "tested_a": tested,
        "time": elapsed,
    }


# ------------------------------------------------------------
# One case
# ------------------------------------------------------------

def run_case(n, p, q, offset, rng):
    R, r1, r2 = make_R(n, offset, rng)

    K, T, E, k, l, true_a, true_b = (
        representation_from_n(
            n,
            p,
            q,
            R,
            r1,
            r2,
        )
    )

    result = solve_ab(
        n,
        r1,
        r2,
        K,
        E,
        k,
        l,
    )

    solved = (
        result["p"] is not None
        and result["p"] * result["q"] == n
    )

    correct_ab = (
        solved
        and {
            result["a"], result["b"]
        } == {
            true_a, true_b
        }
    )

    return {
        "n": n,
        "p": p,
        "q": q,
        "R": R,
        "r1": r1,
        "r2": r2,
        "K": K,
        "T": T,
        "E": E,
        "k": k,
        "l": l,
        "true_a": true_a,
        "true_b": true_b,
        "result": result,
        "solved": solved,
        "correct_ab": correct_ab,
    }


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def main():
    print("=" * 92)
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print("solve_ab TIMING AT 1e9 AND 1e16")
    print("=" * 92)

    print()
    print("configuration")
    print(f"    scales              = {[scale_label(x) for x in SCALES]}")
    print(f"    cases / scale       = {CASES_PER_SCALE}")
    print(f"    R offsets            = {R_OFFSETS}")
    print(f"    seed                 = {SEED}")
    print()

    rng = random.Random(SEED)

    all_results = []

    for scale in SCALES:

        print("=" * 92)
        print(f"SCALE {scale_label(scale)}")
        print("=" * 92)

        scale_results = []

        for case_id in range(1, CASES_PER_SCALE + 1):

            n, p, q = make_semiprime(scale, rng)

            offset = rng.choice(R_OFFSETS)

            result = run_case(
                n,
                p,
                q,
                offset,
                rng,
            )

            scale_results.append(result)
            all_results.append(result)

            s = result["result"]

            print()
            print(
                f"case {case_id}/{CASES_PER_SCALE} "
                f"offset={offset:.2f}"
            )

            print(
                f"    n={fmt(n)}"
            )

            print(
                f"    true p,q=({fmt(p)},{fmt(q)})"
            )

            print(
                f"    R={fmt(result['R'])} "
                f"(r1,r2)=({fmt(result['r1'])},{fmt(result['r2'])})"
            )

            print(
                f"    K={fmt(result['K'])} "
                f"k,l=({fmt(result['k'])},{fmt(result['l'])}) "
                f"T={fmt(result['T'])} "
                f"E={fmt(result['E'])}"
            )

            print(
                f"    true a,b=({fmt(result['true_a'])},"
                f"{fmt(result['true_b'])})"
            )

            print(
                f"    recovered a,b=("
                f"{fmt(s['a']) if s['a'] is not None else '-'},"
                f"{fmt(s['b']) if s['b'] is not None else '-'}"
                f")"
            )

            print(
                f"    tested a values   = {fmt(s['tested_a'])} "
                f"/ {fmt(result['r1'])}"
            )

            print(
                f"    solve_ab time     = {s['time']:.8f}s"
            )

            print(
                f"    solved            = "
                f"{'YES' if result['solved'] else 'NO'}"
            )

            print(
                f"    correct (a,b)     = "
                f"{'YES' if result['correct_ab'] else 'NO'}"
            )

        # ----------------------------------------------------
        # Scale summary
        # ----------------------------------------------------

        times = [
            x["result"]["time"]
            for x in scale_results
        ]

        tested = [
            x["result"]["tested_a"]
            for x in scale_results
        ]

        solved = sum(
            x["solved"]
            for x in scale_results
        )

        correct = sum(
            x["correct_ab"]
            for x in scale_results
        )

        print()
        print("-" * 92)
        print(f"{scale_label(scale)} SUMMARY")
        print("-" * 92)

        print(
            f"cases                 = {len(scale_results)}"
        )

        print(
            f"solved                = "
            f"{solved}/{len(scale_results)}"
        )

        print(
            f"correct a,b           = "
            f"{correct}/{len(scale_results)}"
        )

        print(
            f"mean solve_ab time    = "
            f"{statistics.mean(times):.8f}s"
        )

        print(
            f"median solve_ab time  = "
            f"{statistics.median(times):.8f}s"
        )

        print(
            f"min solve_ab time     = "
            f"{min(times):.8f}s"
        )

        print(
            f"max solve_ab time     = "
            f"{max(times):.8f}s"
        )

        print(
            f"mean a values tested  = "
            f"{statistics.mean(tested):,.1f}"
        )

        print(
            f"max a values tested   = "
            f"{max(tested):,}"
        )

    # --------------------------------------------------------
    # Cross-scale comparison
    # --------------------------------------------------------

    print()
    print("=" * 92)
    print("CROSS-SCALE SOLVE_AB TIMING")
    print("=" * 92)

    print(
        "scale       cases   mean-time      median-time     "
        "min-time       max-time"
    )
    print("-" * 92)

    for scale in SCALES:

        rows = [
            x for x in all_results
            if (
                x["n"] >= scale // 2
                and x["n"] < scale * 2
            )
        ]

        times = [
            x["result"]["time"]
            for x in rows
        ]

        print(
            f"{scale_label(scale):<11}"
            f"{len(rows):<8}"
            f"{statistics.mean(times):<15.8f}"
            f"{statistics.median(times):<16.8f}"
            f"{min(times):<15.8f}"
            f"{max(times):<15.8f}"
        )

    # --------------------------------------------------------
    # Detailed timing data
    # --------------------------------------------------------

    print()
    print("=" * 92)
    print("DETAILED SOLVE_AB RESULTS")
    print("=" * 92)

    for idx, x in enumerate(all_results, 1):

        s = x["result"]

        print(
            f"case {idx:02d} "
            f"n={fmt(x['n'])} "
            f"R={fmt(x['R'])} "
            f"K={fmt(x['K'])} "
            f"E={fmt(x['E'])} "
            f"a-tested={fmt(s['tested_a'])} "
            f"time={s['time']:.8f}s "
            f"valid={x['correct_ab']}"
        )

    # --------------------------------------------------------
    # Important algebraic consistency
    # --------------------------------------------------------

    print()
    print("=" * 92)
    print("ALGEBRAIC CONSISTENCY")
    print("=" * 92)

    failures = 0

    for x in all_results:

        n = x["n"]
        R = x["R"]
        K = x["K"]
        E = x["E"]

        expected_E = n // R - K

        if expected_E != E:
            failures += 1

    print(
        f"E = floor(n/R)-K failures = {failures}"
    )

    print()
    print("=" * 92)
    print("INTERPRETATION TARGET")
    print("=" * 92)

    print(
        """
The experiment measures only the cost of:

    GIVEN:
        n
        R = r1*r2
        K
        k
        l
        E

    FIND:
        a,b

using:

    p = r1*k + a
    q = r2*l + b
    p*q = n

For every candidate a:

    p_candidate = r1*k + a

and if p_candidate divides n:

    q_candidate = n/p_candidate
    b = q_candidate - r2*l

The benchmark therefore isolates the remaining (a,b)
search after K, k and l are already known.

The experiment also verifies:

    E = floor(n/R) - K

and checks the recovered p*q against n.
"""
    )

    print("=" * 92)
    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")
    print("=" * 92)


if __name__ == "__main__":
    main()
