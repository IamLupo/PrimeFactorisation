#!/usr/bin/env python3

"""
============================================================================================
EXPERIMENT 82
C1/C2/C3 SOLVE_AB WITH R ~= n/100
K,E,k,l -> solve_ab() -> a,b -> p,q
============================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time


# ==========================================================================================
# PRIME UTILITIES
# ==========================================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n % p == 0:
            return n == p

    d = 41
    step = 2

    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step

    return True


def next_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


def prev_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n -= 1

    while n > 2 and not is_prime(n):
        n -= 2

    return n


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def generate_prime_near(
    center: int,
    rng: random.Random,
) -> int:

    if center < 100:
        center = 100

    spread = max(100, center // 100)

    for _ in range(10000):

        candidate = center + rng.randint(
            -spread,
            spread,
        )

        if candidate < 100:
            continue

        candidate |= 1

        if is_prime(candidate):
            return candidate

    return next_prime(center)


# ==========================================================================================
# CHOOSE R ~= n/100
# ==========================================================================================

def choose_r_pair(
    n: int,
    p: int,
    q: int,
    rng: random.Random,
):
    """
    Select prime r1,r2 such that:

        R = r1*r2 ~= n/100

    while enforcing:

        r1 < p
        r2 < q

    so that:

        k = p//r1 >= 1
        l = q//r2 >= 1
    """

    target_R = max(1, n // 100)
    root_R = math.isqrt(target_R)

    candidates = []

    # Search a relatively wide region around sqrt(n/100).
    search_radius = max(
        100,
        root_R // 4,
    )

    start = max(
        3,
        root_R - search_radius,
    )

    end = root_R + search_radius

    # Sample rather than testing every integer for very large values.
    step = max(
        1,
        (end - start) // 5000,
    )

    for center in range(
        start,
        end + 1,
        step,
    ):

        r1 = prev_prime(center)

        if r1 >= p:
            continue

        # Target r2 based on desired product.
        desired_r2 = max(
            3,
            target_R // r1,
        )

        # Check several nearby prime candidates.
        for delta in (
            0,
            -10,
            10,
            -50,
            50,
            -100,
            100,
            -250,
            250,
        ):

            candidate = desired_r2 + delta

            if candidate < 3:
                continue

            r2 = prev_prime(candidate)

            if r2 >= q:
                continue

            R = r1 * r2

            relative_error = abs(
                R - target_R
            ) / target_R

            candidates.append(
                (
                    relative_error,
                    abs(R - target_R),
                    r1,
                    r2,
                )
            )

    if not candidates:
        raise RuntimeError(
            "Could not construct valid r1,r2"
        )

    candidates.sort(
        key=lambda x: (x[0], x[1])
    )

    # Randomize slightly among very good candidates.
    best = candidates[:25]

    _, _, r1, r2 = rng.choice(best)

    if not (r1 < p and r2 < q):
        raise RuntimeError(
            "Internal R construction failure"
        )

    return r1, r2


# ==========================================================================================
# SOLVE_AB
# ==========================================================================================

def solve_ab(
    n: int,
    r1: int,
    r2: int,
    K: int,
    E: int,
    k: int,
    l: int,
):

    """
    Recover all (a,b) consistent with:

        n = (a + r1*k) * (b + r2*l)

        K = k*l

        E = c1 + c2 + c3

    where:

        c1 = floor(k*b / r2)
        c2 = floor(l*a / r1)

        beta  = (k*b) % r2
        alpha = (l*a) % r1

        c3 = floor(
            (r1*beta + r2*alpha + a*b)
            / (r1*r2)
        )
    """

    if k <= 0 or l <= 0:
        return []

    if k * l != K:
        return []

    if E < 0:
        return []

    R = r1 * r2

    solutions = []

    # ------------------------------------------------------------------
    # Carry states
    # ------------------------------------------------------------------

    for c1 in range(E + 1):

        for c2 in range(E - c1 + 1):

            c3 = E - c1 - c2

            # ----------------------------------------------------------
            # a-cell
            # ----------------------------------------------------------

            a_min = (
                c2 * r1 + l - 1
            ) // l

            a_max = (
                (c2 + 1) * r1 - 1
            ) // l

            a_min = max(
                a_min,
                0,
            )

            a_max = min(
                a_max,
                r1 - 1,
            )

            if a_min > a_max:
                continue

            # ----------------------------------------------------------
            # b-cell
            # ----------------------------------------------------------

            b_min = (
                c1 * r2 + k - 1
            ) // k

            b_max = (
                (c1 + 1) * r2 - 1
            ) // k

            b_min = max(
                b_min,
                0,
            )

            b_max = min(
                b_max,
                r2 - 1,
            )

            if b_min > b_max:
                continue

            width_a = (
                a_max - a_min + 1
            )

            width_b = (
                b_max - b_min + 1
            )

            # ----------------------------------------------------------
            # Scan smaller dimension
            # ----------------------------------------------------------

            if width_a <= width_b:

                for a in range(
                    a_min,
                    a_max + 1,
                ):

                    P = (
                        a
                        + r1 * k
                    )

                    if P <= 0:
                        continue

                    if n % P != 0:
                        continue

                    Q = n // P

                    b = (
                        Q
                        - r2 * l
                    )

                    if not (
                        b_min <= b <= b_max
                    ):
                        continue

                    c1_check = (
                        k * b
                    ) // r2

                    c2_check = (
                        l * a
                    ) // r1

                    beta = (
                        k * b
                    ) % r2

                    alpha = (
                        l * a
                    ) % r1

                    c3_check = (
                        r1 * beta
                        + r2 * alpha
                        + a * b
                    ) // R

                    if c1_check != c1:
                        continue

                    if c2_check != c2:
                        continue

                    if c3_check != c3:
                        continue

                    if (
                        c1_check
                        + c2_check
                        + c3_check
                        != E
                    ):
                        continue

                    if P * Q != n:
                        continue

                    solutions.append(
                        {
                            "a": a,
                            "b": b,
                            "c1": c1,
                            "c2": c2,
                            "c3": c3,
                            "E": E,
                        }
                    )

            else:

                for b in range(
                    b_min,
                    b_max + 1,
                ):

                    Q = (
                        b
                        + r2 * l
                    )

                    if Q <= 0:
                        continue

                    if n % Q != 0:
                        continue

                    P = n // Q

                    a = (
                        P
                        - r1 * k
                    )

                    if not (
                        a_min <= a <= a_max
                    ):
                        continue

                    c1_check = (
                        k * b
                    ) // r2

                    c2_check = (
                        l * a
                    ) // r1

                    beta = (
                        k * b
                    ) % r2

                    alpha = (
                        l * a
                    ) % r1

                    c3_check = (
                        r1 * beta
                        + r2 * alpha
                        + a * b
                    ) // R

                    if c1_check != c1:
                        continue

                    if c2_check != c2:
                        continue

                    if c3_check != c3:
                        continue

                    if (
                        c1_check
                        + c2_check
                        + c3_check
                        != E
                    ):
                        continue

                    if P * Q != n:
                        continue

                    solutions.append(
                        {
                            "a": a,
                            "b": b,
                            "c1": c1,
                            "c2": c2,
                            "c3": c3,
                            "E": E,
                        }
                    )

    # ------------------------------------------------------------------
    # Deduplicate
    # ------------------------------------------------------------------

    unique = {}

    for s in solutions:

        key = (
            s["a"],
            s["b"],
            s["c1"],
            s["c2"],
            s["c3"],
        )

        unique[key] = s

    return list(
        unique.values()
    )


# ==========================================================================================
# CALCULATE TRUE PARAMETERS
# ==========================================================================================

def calculate_parameters(
    p: int,
    q: int,
    r1: int,
    r2: int,
):

    n = p * q
    R = r1 * r2

    if p <= r1:
        raise ValueError(
            f"p <= r1: p={p} r1={r1}"
        )

    if q <= r2:
        raise ValueError(
            f"q <= r2: q={q} r2={r2}"
        )

    k, a = divmod(
        p,
        r1,
    )

    l, b = divmod(
        q,
        r2,
    )

    if k <= 0 or l <= 0:
        raise ValueError(
            f"k/l must be positive: "
            f"k={k} l={l}"
        )

    K = k * l

    quotient, remainder = divmod(
        n,
        R,
    )

    E = quotient - K

    if E < 0:
        raise ValueError(
            f"E < 0: "
            f"n={n} R={R} "
            f"quotient={quotient} K={K}"
        )

    # Exact carry decomposition.

    c1 = (
        k * b
    ) // r2

    c2 = (
        l * a
    ) // r1

    beta = (
        k * b
    ) % r2

    alpha = (
        l * a
    ) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    # Identity.

    if E != c1 + c2 + c3:
        raise AssertionError(
            "E != c1+c2+c3"
        )

    if (
        R * K
        + r1 * k * b
        + r2 * l * a
        + a * b
        != n
    ):
        raise AssertionError(
            "factor/carry identity failed"
        )

    return {
        "n": n,
        "R": R,
        "k": k,
        "l": l,
        "K": K,
        "E": E,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "remainder": remainder,
    }


# ==========================================================================================
# GENERATE ONE VALID CASE
# ==========================================================================================

def generate_case(
    n_target: int,
    rng: random.Random,
):

    sqrt_n = math.isqrt(
        n_target
    )

    # Keep p,q comfortably above sqrt(n/100),
    # so r1,r2 cannot accidentally consume them.

    for _ in range(1000):

        p = generate_prime_near(
            sqrt_n,
            rng,
        )

        q = generate_prime_near(
            sqrt_n,
            rng,
        )

        n = p * q

        r1, r2 = choose_r_pair(
            n,
            p,
            q,
            rng,
        )

        if r1 >= p or r2 >= q:
            continue

        try:

            params = calculate_parameters(
                p,
                q,
                r1,
                r2,
            )

        except ValueError:
            continue

        return (
            p,
            q,
            r1,
            r2,
            params,
        )

    raise RuntimeError(
        "Unable to generate valid case"
    )


# ==========================================================================================
# RUN SCALE
# ==========================================================================================

def run_scale(
    label: str,
    n_target: int,
    cases: int,
    rng: random.Random,
):

    print("=" * 92)
    print(f"SCALE {label}")
    print("=" * 92)

    records = []

    for case_index in range(
        1,
        cases + 1,
    ):

        (
            p,
            q,
            r1,
            r2,
            params,
        ) = generate_case(
            n_target,
            rng,
        )

        n = params["n"]
        R = params["R"]
        K = params["K"]
        E = params["E"]
        k = params["k"]
        l = params["l"]
        a_true = params["a"]
        b_true = params["b"]

        ratio = n / R

        # --------------------------------------------------------------
        # solve_ab benchmark
        # --------------------------------------------------------------

        t0 = time.perf_counter()

        solutions = solve_ab(
            n,
            r1,
            r2,
            K,
            E,
            k,
            l,
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        true_found = any(
            s["a"] == a_true
            and s["b"] == b_true
            for s in solutions
        )

        factor_correct = False

        recovered = []

        for s in solutions:

            p2 = (
                r1 * k
                + s["a"]
            )

            q2 = (
                r2 * l
                + s["b"]
            )

            recovered.append(
                (p2, q2)
            )

            if (
                p2 * q2 == n
                and (
                    (p2 == p and q2 == q)
                    or
                    (p2 == q and q2 == p)
                )
            ):
                factor_correct = True

        print()
        print(
            f"CASE {case_index}/{cases}"
        )

        print(
            f"    n={n:,}"
        )

        print(
            f"    true p,q=({p:,},{q:,})"
        )

        print(
            f"    r1,r2=({r1:,},{r2:,}) "
            f"R={R:,}"
        )

        print(
            f"    n/R={ratio:.12f}"
        )

        print()
        print("TRUE:")

        print(
            f"    K={K:,} E={E:,}"
        )

        print(
            f"    (k,l)=({k:,},{l:,})"
        )

        print(
            f"    (a,b)=({a_true:,},{b_true:,})"
        )

        print()
        print("CARRY:")

        print(
            f"    c1={params['c1']:,}"
            f" c2={params['c2']:,}"
            f" c3={params['c3']:,}"
        )

        print(
            f"    c1+c2+c3="
            f"{params['c1'] + params['c2'] + params['c3']:,}"
        )

        print()
        print("SOLVE_AB:")

        print(
            f"    solutions={len(solutions)}"
        )

        print(
            f"    true (a,b) found="
            f"{true_found}"
        )

        print(
            f"    correct factors="
            f"{factor_correct}"
        )

        print(
            f"    solve_ab time="
            f"{elapsed:.9f}s"
        )

        for idx, s in enumerate(
            solutions[:10],
            start=1,
        ):

            p2 = (
                r1 * k
                + s["a"]
            )

            q2 = (
                r2 * l
                + s["b"]
            )

            print(
                f"    solution {idx}: "
                f"(a,b)=("
                f"{s['a']:,},"
                f"{s['b']:,}) "
                f"(p,q)=("
                f"{p2:,},"
                f"{q2:,}) "
                f"(c1,c2,c3)=("
                f"{s['c1']:,},"
                f"{s['c2']:,},"
                f"{s['c3']:,})"
            )

        records.append(
            {
                "n": n,
                "R": R,
                "ratio": ratio,
                "p": p,
                "q": q,
                "r1": r1,
                "r2": r2,
                "K": K,
                "E": E,
                "k": k,
                "l": l,
                "a": a_true,
                "b": b_true,
                "solutions": len(solutions),
                "true_found": true_found,
                "factor_correct": factor_correct,
                "time": elapsed,
            }
        )

    # ======================================================================================
    # SUMMARY
    # ======================================================================================

    times = [
        x["time"]
        for x in records
    ]

    E_values = [
        x["E"]
        for x in records
    ]

    ratios = [
        x["ratio"]
        for x in records
    ]

    print()
    print("-" * 92)
    print(f"{label} SUMMARY")
    print("-" * 92)

    print(
        f"cases                  = {len(records)}"
    )

    print(
        f"true (a,b) recovered   = "
        f"{sum(x['true_found'] for x in records)}"
        f"/{len(records)}"
    )

    print(
        f"correct factors        = "
        f"{sum(x['factor_correct'] for x in records)}"
        f"/{len(records)}"
    )

    print(
        f"mean n/R               = "
        f"{statistics.mean(ratios):.9f}"
    )

    print(
        f"max |n/R-100|          = "
        f"{max(abs(x - 100) for x in ratios):.9f}"
    )

    print(
        f"mean E                 = "
        f"{statistics.mean(E_values):.3f}"
    )

    print(
        f"max E                  = "
        f"{max(E_values):,}"
    )

    print(
        f"mean solve_ab time     = "
        f"{statistics.mean(times):.9f}s"
    )

    print(
        f"median solve_ab time   = "
        f"{statistics.median(times):.9f}s"
    )

    print(
        f"min solve_ab time      = "
        f"{min(times):.9f}s"
    )

    print(
        f"max solve_ab time      = "
        f"{max(times):.9f}s"
    )

    return records


# ==========================================================================================
# MAIN
# ==========================================================================================

def main():

    SEED = 1511464998

    rng = random.Random(
        SEED
    )

    CASES_PER_SCALE = 10

    scales = (
        ("1e9", 10**9),
        ("1e12", 10**12),
        ("1e16", 10**16),
    )

    print("=" * 92)
    print("START EXPERIMENT 82")
    print("C1/C2/C3 SOLVE_AB WITH R ~= n/100")
    print("=" * 92)

    print()
    print("configuration")

    print(
        f"    scales = "
        f"{[x[0] for x in scales]}"
    )

    print(
        f"    cases / scale = "
        f"{CASES_PER_SCALE}"
    )

    print(
        f"    target n/R = 100.0"
    )

    print(
        f"    seed = {SEED}"
    )

    print(
        f"    solver = C1/C2/C3 carry-cell solver"
    )

    total_start = (
        time.perf_counter()
    )

    all_records = []

    for label, scale in scales:

        records = run_scale(
            label,
            scale,
            CASES_PER_SCALE,
            rng,
        )

        all_records.extend(
            records
        )

    total_time = (
        time.perf_counter()
        - total_start
    )

    # ======================================================================================
    # GLOBAL SUMMARY
    # ======================================================================================

    print()
    print("=" * 92)
    print("GLOBAL SUMMARY")
    print("=" * 92)

    total_cases = len(
        all_records
    )

    total_true = sum(
        x["true_found"]
        for x in all_records
    )

    total_correct = sum(
        x["factor_correct"]
        for x in all_records
    )

    all_times = [
        x["time"]
        for x in all_records
    ]

    all_E = [
        x["E"]
        for x in all_records
    ]

    all_ratios = [
        x["ratio"]
        for x in all_records
    ]

    print(
        f"total cases            = "
        f"{total_cases}"
    )

    print(
        f"true (a,b) recovered   = "
        f"{total_true}/{total_cases}"
    )

    print(
        f"correct factors        = "
        f"{total_correct}/{total_cases}"
    )

    print(
        f"maximum E              = "
        f"{max(all_E):,}"
    )

    print(
        f"mean E                 = "
        f"{statistics.mean(all_E):.3f}"
    )

    print(
        f"mean n/R               = "
        f"{statistics.mean(all_ratios):.9f}"
    )

    print(
        f"max |n/R-100|          = "
        f"{max(abs(x - 100) for x in all_ratios):.9f}"
    )

    print(
        f"mean solve_ab time     = "
        f"{statistics.mean(all_times):.9f}s"
    )

    print(
        f"median solve_ab time   = "
        f"{statistics.median(all_times):.9f}s"
    )

    print(
        f"max solve_ab time      = "
        f"{max(all_times):.9f}s"
    )

    # ======================================================================================
    # CONSISTENCY
    # ======================================================================================

    print()
    print("=" * 92)
    print("ALGEBRAIC CONSISTENCY")
    print("=" * 92)

    failures = 0

    for x in all_records:

        if (
            x["E"]
            != x["n"] // x["R"]
            - x["K"]
        ):
            failures += 1

    print(
        f"E = floor(n/R) - K failures = "
        f"{failures}"
    )

    print()
    print(
        f"total runtime = "
        f"{total_time:.6f}s"
    )

    print()
    print("=" * 92)
    print("FINISHED EXPERIMENT 82")
    print("=" * 92)


if __name__ == "__main__":
    main()