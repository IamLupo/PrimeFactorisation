#!/usr/bin/env python3

import math
import random
import time


# ============================================================================
# EXPERIMENT 81
# R ~= n / 100 -> BRUTE E -> K -> (k,l) -> solve_ab -> (p,q)
# ============================================================================

SEED = 1511464998

CASES_PER_SCALE = 20

# We want:
#
#     R = r1 * r2 ~= n / 100
#
# so that:
#
#     floor(n / R) ~= 100
#
# This makes K and E very small compared with the earlier experiments.

TARGET_RATIO = 100.0

# Brute-force E.
E_MAX = 1001

# Static K lookup.
K_MAX = 5000


# ============================================================================
# PRIME GENERATION
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step

    return True


def previous_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n -= 1

    while n >= 2:
        if is_prime(n):
            return n
        n -= 2

    raise RuntimeError("No previous prime found")


def next_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while True:
        if is_prime(n):
            return n
        n += 2


# ============================================================================
# STATIC K -> (k,l)
# ============================================================================

def build_static_k_table(k_max: int):
    table = {}

    for k in range(1, k_max + 1):
        for l in range(1, k_max // k + 1):
            K = k * l

            if K > k_max:
                break

            table.setdefault(K, []).append((k, l))

    pair_count = sum(len(v) for v in table.values())
    return table, pair_count


# ============================================================================
# solve_ab
#
# Given:
#
#     n
#     R = r1*r2
#     K
#     E
#     k
#     l
#
# solve:
#
#     (r1*k+a)(r2*l+b) = n
#
# Since:
#
#     E = floor(n/R) - K
#
# and
#
#     0 <= a < r1
#     0 <= b < r2
#
# we can narrow a using the quotient-cell constraint.
#
# For each candidate a:
#
#     p = r1*k + a
#
# and if p divides n:
#
#     q = n // p
#
# then:
#
#     b = q - r2*l
#
# ============================================================================

def solve_ab(n: int, r1: int, r2: int, K: int, E: int, k: int, l: int):
    if k <= 0 or l <= 0:
        return None, None, 0

    if k * l != K:
        return None, None, 0

    R = r1 * r2

    # Necessary quotient-cell condition:
    #
    # n = R*K + cross terms + a*b
    #
    # and:
    #
    # E = floor(n/R) - K
    #
    # Therefore:
    #
    # E*R <= cross < (E+1)*R
    #
    # We use the upper bound to obtain an a interval.
    #
    # cross = r1*k*b + r2*l*a + a*b
    #
    # Since b >= 0:
    #
    # cross >= r2*l*a
    #
    # so:
    #
    # a <= ((E+1)*R - 1) // (r2*l)
    #
    if r2 * l <= 0:
        return None, None, 0

    a_max_by_E = ((E + 1) * R - 1) // (r2 * l)
    a_max = min(r1 - 1, a_max_by_E)

    if a_max < 0:
        return None, None, 0

    tested = 0

    for a in range(a_max + 1):
        tested += 1

        p = r1 * k + a

        if p <= 0:
            continue

        if n % p != 0:
            continue

        q = n // p
        b = q - r2 * l

        if not (0 <= b < r2):
            continue

        if p * q != n:
            continue

        # Exact E verification.
        E_check = (n // R) - K
        if E_check != E:
            continue

        # Final representation verification.
        if p != r1 * k + a:
            continue

        if q != r2 * l + b:
            continue

        return p, q, tested

    return None, None, tested


# ============================================================================
# GENERATE CASE
#
# We first choose r1,r2 and therefore:
#
#     R = r1*r2
#
# Then choose p,q near sqrt(100R), so:
#
#     n = p*q ~= 100R
#
# This directly produces the desired scale.
# ============================================================================

def generate_case(target_n: int, rng: random.Random):
    # Desired:
    #
    #     R ~= target_n / 100
    #
    target_R = max(100, target_n // 100)

    root_R = math.isqrt(target_R)

    # Pick two nearby static-ish primes.
    r1 = previous_prime(root_R + rng.randint(-50, 50))
    r2 = next_prime(root_R + rng.randint(-50, 50))

    # Make sure R is sensible.
    R = r1 * r2

    # We want p,q around sqrt(100R).
    target_factor = math.isqrt(100 * R)

    p = previous_prime(
        target_factor + rng.randint(-target_factor // 1000, target_factor // 1000)
    )

    q = next_prime(
        target_factor + rng.randint(-target_factor // 1000, target_factor // 1000)
    )

    n = p * q

    return n, p, q, r1, r2


# ============================================================================
# REPRESENTATION
# ============================================================================

def calculate_parameters(n: int, p: int, q: int, r1: int, r2: int):
    R = r1 * r2

    k = p // r1
    l = q // r2

    a = p % r1
    b = q % r2

    K = k * l

    Q = n // R
    E = Q - K

    # Exact identity.
    assert p == r1 * k + a
    assert q == r2 * l + b
    assert n == (r1 * k + a) * (r2 * l + b)
    assert K == k * l
    assert E == (n // R) - K

    return R, K, E, k, l, a, b


# ============================================================================
# SEARCH
# ============================================================================

def search_factorization(n: int, r1: int, r2: int, static_table):
    R = r1 * r2
    Q = n // R

    E_tests = 0
    pair_tests = 0
    solve_calls = 0
    total_a_tests = 0

    # Since:
    #
    #     E = Q - K
    #
    # and we assume E is small:
    #
    #     K = Q - E
    #
    # Therefore test E directly.
    #
    for E in range(0, E_MAX + 1):
        E_tests += 1

        K = Q - E

        if K <= 0:
            continue

        pairs = static_table.get(K)
        if not pairs:
            continue

        for k, l in pairs:
            pair_tests += 1

            solve_calls += 1

            p, q, a_tests = solve_ab(
                n,
                r1,
                r2,
                K,
                E,
                k,
                l,
            )

            total_a_tests += a_tests

            if p is not None:
                # Final factor verification.
                if p * q != n:
                    continue

                return {
                    "solved": True,
                    "p": p,
                    "q": q,
                    "K": K,
                    "E": E,
                    "k": k,
                    "l": l,
                    "a_tests": total_a_tests,
                    "E_tests": E_tests,
                    "pair_tests": pair_tests,
                    "solve_calls": solve_calls,
                }

    return {
        "solved": False,
        "p": None,
        "q": None,
        "K": None,
        "E": None,
        "k": None,
        "l": None,
        "a_tests": total_a_tests,
        "E_tests": E_tests,
        "pair_tests": pair_tests,
        "solve_calls": solve_calls,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():
    rng = random.Random(SEED)

    print("=" * 100)
    print("START EXPERIMENT 81")
    print("R ~= n / 100 -> BRUTE E -> K -> (k,l) -> solve_ab -> (p,q)")
    print("=" * 100)
    print()
    print("configuration")
    print(f"    cases / scale       = {CASES_PER_SCALE}")
    print(f"    target n/R          = {TARGET_RATIO}")
    print(f"    E brute range       = 0..{E_MAX}")
    print(f"    static K table max  = {K_MAX}")
    print(f"    seed                = {SEED}")
    print()

    # ------------------------------------------------------------------------
    # Static table
    # ------------------------------------------------------------------------

    t0 = time.perf_counter()
    static_table, pair_count = build_static_k_table(K_MAX)
    build_time = time.perf_counter() - t0

    print("=" * 100)
    print("BUILDING STATIC K LOOKUP")
    print("=" * 100)
    print(f"    K entries          = {len(static_table):,}")
    print(f"    divisor pairs      = {pair_count:,}")
    print(f"    build time         = {build_time:.6f}s")
    print()

    scales = [
        ("1e9", 10**9),
        ("1e12", 10**12),
        ("1e16", 10**16),
    ]

    global_cases = []
    global_solved = 0
    global_correct = 0

    for scale_name, target_n in scales:
        print("=" * 100)
        print(f"SCALE {scale_name}")
        print("=" * 100)
        print()

        scale_solved = 0
        scale_correct = 0
        scale_E_small = 0

        sum_E = 0
        max_E = 0
        sum_ratio = 0.0
        max_ratio_error = 0.0
        sum_time = 0.0
        sum_E_tests = 0
        sum_pairs = 0
        sum_ab_calls = 0
        sum_a_tests = 0

        for case_idx in range(1, CASES_PER_SCALE + 1):
            n, p_true, q_true, r1, r2 = generate_case(target_n, rng)

            R, K_true, E_true, k_true, l_true, a_true, b_true = \
                calculate_parameters(
                    n,
                    p_true,
                    q_true,
                    r1,
                    r2,
                )

            ratio = n / R
            ratio_error = abs(ratio - TARGET_RATIO)

            t0 = time.perf_counter()

            result = search_factorization(
                n,
                r1,
                r2,
                static_table,
            )

            elapsed = time.perf_counter() - t0

            solved = result["solved"]

            correct = (
                solved
                and (
                    result["p"] == p_true
                    and result["q"] == q_true
                )
            )

            # Factor order may be reversed.
            correct_unordered = (
                solved
                and (
                    (result["p"] == p_true and result["q"] == q_true)
                    or
                    (result["p"] == q_true and result["q"] == p_true)
                )
            )

            if 0 <= E_true <= E_MAX:
                scale_E_small += 1

            if solved:
                scale_solved += 1

            if correct_unordered:
                scale_correct += 1

            global_cases.append(
                (
                    scale_name,
                    n,
                    p_true,
                    q_true,
                    r1,
                    r2,
                    R,
                    ratio,
                    K_true,
                    E_true,
                    k_true,
                    l_true,
                    a_true,
                    b_true,
                    result,
                    elapsed,
                )
            )

            sum_E += E_true
            max_E = max(max_E, E_true)
            sum_ratio += ratio
            max_ratio_error = max(max_ratio_error, ratio_error)
            sum_time += elapsed
            sum_E_tests += result["E_tests"]
            sum_pairs += result["pair_tests"]
            sum_ab_calls += result["solve_calls"]
            sum_a_tests += result["a_tests"]

            status = "YES" if solved else "NO"
            valid = (
                solved
                and result["p"] * result["q"] == n
            )

            print(f"CASE {case_idx}/{CASES_PER_SCALE}")
            print(f"    n={n:,}")
            print(f"    p,q=({p_true:,},{q_true:,})")
            print(
                f"    r1,r2=({r1:,},{r2:,}) "
                f"R={R:,}"
            )
            print(f"    n/R={ratio:.12f}")
            print()
            print("    TRUE:")
            print(
                f"        K={K_true:,} E={E_true:,}"
            )
            print(
                f"        (k,l)=({k_true:,},{l_true:,})"
            )
            print(
                f"        (a,b)=({a_true:,},{b_true:,})"
            )
            print()
            print("    SEARCH:")
            print(f"        solved={status}")
            print(
                f"        K={result['K'] if result['K'] is not None else '-'} "
                f"E={result['E'] if result['E'] is not None else '-'}"
            )
            print(
                f"        (k,l)=({result['k'] if result['k'] is not None else '-'},"
                f"{result['l'] if result['l'] is not None else '-'})"
            )

            if solved:
                print(
                    f"        recovered p,q=({result['p']:,},{result['q']:,})"
                )

                recovered_a = result["p"] - r1 * result["k"]
                recovered_b = result["q"] - r2 * result["l"]

                print(
                    f"        recovered (a,b)=({recovered_a:,},{recovered_b:,})"
                )

            print(
                f"        E tests={result['E_tests']:,}"
            )
            print(
                f"        divisor-pair tests={result['pair_tests']:,}"
            )
            print(
                f"        solve_ab calls={result['solve_calls']:,}"
            )
            print(
                f"        solve_ab a-tests={result['a_tests']:,}"
            )
            print(f"        valid factorization={valid}")
            print(
                f"        correct unordered factors={correct_unordered}"
            )
            print(f"        time={elapsed:.8f}s")
            print()

        print("-" * 100)
        print(f"{scale_name} SUMMARY")
        print("-" * 100)
        print(f"    cases                = {CASES_PER_SCALE}")
        print(
            f"    E <= {E_MAX:<4}             = "
            f"{scale_E_small}/{CASES_PER_SCALE}"
        )
        print(f"    solved               = {scale_solved}/{CASES_PER_SCALE}")
        print(f"    correct factors      = {scale_correct}/{CASES_PER_SCALE}")
        print(f"    maximum E            = {max_E:,}")
        print(f"    mean n/R             = {sum_ratio / CASES_PER_SCALE:.12f}")
        print(f"    max |n/R-100|        = {max_ratio_error:.12f}")
        print(
            f"    mean E tests        = "
            f"{sum_E_tests / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean pair tests     = "
            f"{sum_pairs / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean solve_ab calls = "
            f"{sum_ab_calls / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean a-tests        = "
            f"{sum_a_tests / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean runtime        = "
            f"{sum_time / CASES_PER_SCALE:.8f}s"
        )
        print()

    # ------------------------------------------------------------------------
    # Global summary
    # ------------------------------------------------------------------------

    global_solved = sum(
        1 for x in global_cases if x[-2]["solved"]
    )

    global_correct = sum(
        1
        for x in global_cases
        if (
            x[-2]["solved"]
            and (
                (x[-2]["p"] == x[2] and x[-2]["q"] == x[3])
                or
                (x[-2]["p"] == x[3] and x[-2]["q"] == x[2])
            )
        )
    )

    global_E_small = sum(
        1
        for x in global_cases
        if 0 <= x[9] <= E_MAX
    )

    max_E = max(x[9] for x in global_cases)

    print("=" * 100)
    print("GLOBAL SUMMARY")
    print("=" * 100)
    print(
        f"    total cases          = {len(global_cases)}"
    )
    print(
        f"    E <= {E_MAX:<4}             = "
        f"{global_E_small}/{len(global_cases)}"
    )
    print(
        f"    solved               = "
        f"{global_solved}/{len(global_cases)}"
    )
    print(
        f"    correct factors      = "
        f"{global_correct}/{len(global_cases)}"
    )
    print(
        f"    maximum E            = {max_E:,}"
    )

    # ------------------------------------------------------------------------
    # Show only cases where the solver did not return the true factor pair.
    # ------------------------------------------------------------------------

    failures = [
        x
        for x in global_cases
        if not (
            x[-2]["solved"]
            and (
                (x[-2]["p"] == x[2] and x[-2]["q"] == x[3])
                or
                (x[-2]["p"] == x[3] and x[-2]["q"] == x[2])
            )
        )
    ]

    print()
    print("=" * 100)
    print("NON-CORRECT CASES")
    print("=" * 100)

    if not failures:
        print("    none")
    else:
        for x in failures[:20]:
            (
                scale_name,
                n,
                p_true,
                q_true,
                r1,
                r2,
                R,
                ratio,
                K_true,
                E_true,
                k_true,
                l_true,
                a_true,
                b_true,
                result,
                elapsed,
            ) = x

            print(
                f"    scale={scale_name} "
                f"n={n:,} "
                f"R={R:,} "
                f"n/R={ratio:.9f} "
                f"true K={K_true} "
                f"true E={E_true}"
            )

            print(
                f"        true p,q=({p_true:,},{q_true:,})"
            )

            print(
                f"        search="
                f"{result['p']},{result['q']} "
                f"K={result['K']} "
                f"E={result['E']} "
                f"k,l=({result['k']},{result['l']})"
            )

    # ------------------------------------------------------------------------
    # Identity check
    # ------------------------------------------------------------------------

    identity_failures = 0

    for x in global_cases:
        (
            scale_name,
            n,
            p_true,
            q_true,
            r1,
            r2,
            R,
            ratio,
            K_true,
            E_true,
            k_true,
            l_true,
            a_true,
            b_true,
            result,
            elapsed,
        ) = x

        if E_true != (n // R) - K_true:
            identity_failures += 1

    print()
    print("=" * 100)
    print("ALGEBRAIC CONSISTENCY")
    print("=" * 100)
    print(
        f"    E = floor(n/R) - K failures = "
        f"{identity_failures}"
    )

    print()
    print("=" * 100)
    print("KEY EQUATIONS")
    print("=" * 100)
    print()
    print("    R = r1*r2")
    print()
    print("    p = r1*k + a")
    print("    q = r2*l + b")
    print()
    print("    K = k*l")
    print()
    print("    n = p*q")
    print()
    print("    E = floor(n/R) - K")
    print()
    print("    Therefore:")
    print()
    print("        K = floor(n/R) - E")
    print()
    print("    The search therefore becomes:")
    print()
    print("        Q = floor(n/R)")
    print("        for E = 0..E_MAX:")
    print("            K = Q-E")
    print("            lookup K -> (k,l)")
    print("            solve_ab(n,r1,r2,K,E,k,l)")
    print("            verify p*q == n")
    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 81")
    print("=" * 100)


if __name__ == "__main__":
    main()
