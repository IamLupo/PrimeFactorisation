#!/usr/bin/env python3

import math
import random
import time


# =============================================================================
# EXPERIMENT 79
#
# Test the idea:
#
#       Q = floor(n / R)
#       E = Q - K
#
# and brute-force K when E is assumed to be small.
#
# For:
#
#       n = p*q
#       p = r1*k + a
#       q = r2*l + b
#       R = r1*r2
#
#       K = k*l
#
# we have exactly:
#
#       n = R*K + r1*k*b + r2*l*a + a*b
#
# and therefore:
#
#       floor(n/R) = K + E
#
# where:
#
#       E = floor((r1*k*b + r2*l*a + a*b)/R)
#
# equivalently:
#
#       E = floor(n/R) - K
#
# =============================================================================


SEED = 1511464998
RNG = random.Random(SEED)

CASES_PER_SCALE = 10

SCALES = (
    10**9,
    10**12,
    10**16,
)

# Static prime r1/r2 pool.
STATIC_PRIMES = (
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)

# We investigate the hypothesis that K and E are small.
K_MIN = 0
K_MAX = 5000

E_MAX = 1001

# Static K -> (k,l) table.
K_TABLE_MAX = K_MAX


# =============================================================================
# PRIME HELPERS
# =============================================================================

def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3
    limit = math.isqrt(n)

    while d <= limit:
        if n % d == 0:
            return False
        d += 2

    return True


def next_prime(n):
    n = max(2, n)

    if n == 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# =============================================================================
# STATIC K TABLE
# =============================================================================

def build_k_lookup(max_k):
    """
    Build:

        K -> ((k,l), ...)

    for all positive divisor pairs satisfying:

        k*l = K

    K=0 is kept as an empty entry.
    """

    table = {}

    for K in range(max_k + 1):

        if K == 0:
            table[K] = ()
            continue

        pairs = []

        root = math.isqrt(K)

        for k in range(1, root + 1):

            if K % k != 0:
                continue

            l = K // k

            pairs.append((k, l))

            if k != l:
                pairs.append((l, k))

        table[K] = tuple(sorted(set(pairs)))

    return table


# =============================================================================
# GENERATE PRIME FACTORS
# =============================================================================

def generate_prime_pair(scale):
    """
    Generate two primes close to sqrt(scale), so n is approximately scale.
    """

    root = math.isqrt(scale)

    p = next_prime(
        root + RNG.randint(-5000, 5000)
    )

    q = next_prime(
        root + RNG.randint(-5000, 5000)
    )

    while q == p:
        q = next_prime(q + 2)

    return p, q


# =============================================================================
# STATIC r1/r2
# =============================================================================

def choose_r_values():
    r1 = RNG.choice(STATIC_PRIMES)
    r2 = RNG.choice(STATIC_PRIMES)

    return r1, r2


# =============================================================================
# EXACT REPRESENTATION
# =============================================================================

def calculate_parameters(p, q, r1, r2):
    """
    Exact quotient/remainder representation:

        p = r1*k + a
        q = r2*l + b

    with:

        0 <= a < r1
        0 <= b < r2
    """

    n = p * q
    R = r1 * r2

    k = p // r1
    a = p % r1

    l = q // r2
    b = q % r2

    K = k * l

    Q = n // R
    E = Q - K

    # Exact expanded identity.
    expanded = (
        R * K
        + r1 * k * b
        + r2 * l * a
        + a * b
    )

    assert expanded == n

    # Direct definition of E from the remainder contribution.
    remainder_contribution = (
        r1 * k * b
        + r2 * l * a
        + a * b
    )

    E_direct = remainder_contribution // R

    assert E == E_direct

    assert 0 <= a < r1
    assert 0 <= b < r2

    return {
        "n": n,
        "R": R,
        "Q": Q,
        "K": K,
        "E": E,
        "k": k,
        "l": l,
        "a": a,
        "b": b,
    }


# =============================================================================
# SOLVE a,b
# =============================================================================

def solve_ab(n, r1, r2, k, l):
    """
    Recover a,b after k,l are known.

    p = r1*k + a

    where:

        0 <= a < r1

    For each possible a, test whether p divides n.
    """

    base_p = r1 * k
    base_q = r2 * l

    tested = 0

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

        if p_candidate * q_candidate != n:
            continue

        return {
            "a": a,
            "b": b,
            "p": p_candidate,
            "q": q_candidate,
            "tested_a": tested,
        }

    return None


# =============================================================================
# SOLVE K/E
# =============================================================================

def solve_K_E(n, r1, r2, k_lookup):
    """
    Brute-force K using:

        Q = floor(n/R)
        E = Q-K

    We do NOT brute-force K and E independently.

    For each candidate K:

        E = Q-K

    and only candidates satisfying:

        0 <= E <= E_MAX

    are considered.

    Then use the static K lookup to obtain all possible (k,l).

    Finally solve a,b.
    """

    R = r1 * r2
    Q = n // R

    K_tests = 0
    pair_tests = 0
    solve_ab_calls = 0

    started = time.perf_counter()

    for K in range(K_MIN, K_MAX + 1):

        K_tests += 1

        E = Q - K

        # Since E is assumed to be small:
        if E < 0:
            continue

        if E > E_MAX:
            continue

        pairs = k_lookup.get(K, ())

        for k, l in pairs:

            pair_tests += 1

            if k * l != K:
                continue

            solve_ab_calls += 1

            result = solve_ab(
                n,
                r1,
                r2,
                k,
                l,
            )

            if result is None:
                continue

            elapsed = time.perf_counter() - started

            return {
                "solved": True,
                "K": K,
                "E": E,
                "k": k,
                "l": l,
                "a": result["a"],
                "b": result["b"],
                "p": result["p"],
                "q": result["q"],
                "tested_a": result["tested_a"],
                "K_tests": K_tests,
                "pair_tests": pair_tests,
                "solve_ab_calls": solve_ab_calls,
                "time": elapsed,
                "Q": Q,
            }

    elapsed = time.perf_counter() - started

    return {
        "solved": False,
        "K": None,
        "E": None,
        "k": None,
        "l": None,
        "a": None,
        "b": None,
        "p": None,
        "q": None,
        "tested_a": 0,
        "K_tests": K_tests,
        "pair_tests": pair_tests,
        "solve_ab_calls": solve_ab_calls,
        "time": elapsed,
        "Q": Q,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print("START EXPERIMENT 79")
    print("BRUTE-FORCE K USING E = FLOOR(n/R) - K")
    print("=" * 100)
    print()

    print("configuration")
    print(
        f"    scales             = "
        f"{[f'{x:.0e}' for x in SCALES]}"
    )
    print(
        f"    cases / scale      = {CASES_PER_SCALE}"
    )
    print(
        f"    K brute range      = "
        f"{K_MIN}..{K_MAX}"
    )
    print(
        f"    E allowed          = "
        f"0..{E_MAX}"
    )
    print(
        f"    static K table max = "
        f"{K_TABLE_MAX}"
    )
    print(
        f"    seed               = "
        f"{SEED}"
    )
    print()

    # =========================================================================
    # BUILD TABLE
    # =========================================================================

    print("=" * 100)
    print("BUILDING STATIC K LOOKUP")
    print("=" * 100)

    t0 = time.perf_counter()

    k_lookup = build_k_lookup(
        K_TABLE_MAX
    )

    lookup_time = time.perf_counter() - t0

    total_pairs = sum(
        len(v)
        for v in k_lookup.values()
    )

    print(
        f"    K entries          = "
        f"{len(k_lookup):,}"
    )
    print(
        f"    divisor pairs      = "
        f"{total_pairs:,}"
    )
    print(
        f"    build time         = "
        f"{lookup_time:.6f}s"
    )

    # =========================================================================
    # GLOBAL
    # =========================================================================

    global_cases = 0
    global_exposed = 0
    global_solved = 0
    global_correct = 0

    global_K_tests = 0
    global_pair_tests = 0
    global_ab_calls = 0

    global_E_values = []
    global_K_values = []

    global_runtime = 0.0

    # =========================================================================
    # SCALES
    # =========================================================================

    for scale in SCALES:

        print()
        print("=" * 100)
        print(f"SCALE {scale:.0e}")
        print("=" * 100)

        scale_cases = 0
        scale_exposed = 0
        scale_solved = 0
        scale_correct = 0

        scale_K_tests = 0
        scale_pair_tests = 0
        scale_ab_calls = 0

        scale_runtime = 0.0

        scale_E_values = []
        scale_K_values = []

        for case_no in range(
            1,
            CASES_PER_SCALE + 1
        ):

            p, q = generate_prime_pair(
                scale
            )

            r1, r2 = choose_r_values()

            actual = calculate_parameters(
                p,
                q,
                r1,
                r2,
            )

            n = actual["n"]
            R = actual["R"]
            Q = actual["Q"]

            K_true = actual["K"]
            E_true = actual["E"]

            k_true = actual["k"]
            l_true = actual["l"]

            a_true = actual["a"]
            b_true = actual["b"]

            scale_cases += 1
            global_cases += 1

            scale_E_values.append(E_true)
            scale_K_values.append(K_true)

            global_E_values.append(E_true)
            global_K_values.append(K_true)

            exposed = (
                K_MIN <= K_true <= K_MAX
                and
                0 <= E_true <= E_MAX
            )

            if exposed:
                scale_exposed += 1
                global_exposed += 1

            result = solve_K_E(
                n,
                r1,
                r2,
                k_lookup,
            )

            scale_runtime += result["time"]
            global_runtime += result["time"]

            scale_K_tests += result["K_tests"]
            scale_pair_tests += result["pair_tests"]
            scale_ab_calls += result["solve_ab_calls"]

            global_K_tests += result["K_tests"]
            global_pair_tests += result["pair_tests"]
            global_ab_calls += result["solve_ab_calls"]

            if result["solved"]:

                scale_solved += 1
                global_solved += 1

                correct = (
                    result["K"] == K_true
                    and result["E"] == E_true
                    and result["k"] == k_true
                    and result["l"] == l_true
                    and result["a"] == a_true
                    and result["b"] == b_true
                    and result["p"] == p
                    and result["q"] == q
                )

                if correct:
                    scale_correct += 1
                    global_correct += 1
                else:
                    correct = False

            else:
                correct = False

            # =================================================================
            # CASE OUTPUT
            # =================================================================

            print()
            print(
                f"CASE {case_no}/{CASES_PER_SCALE}"
            )

            print(
                f"    n={n:,}"
            )

            print(
                f"    p,q=({p:,},{q:,})"
            )

            print(
                f"    R={R:,} "
                f"(r1,r2)=({r1},{r2})"
            )

            print(
                f"    Q=floor(n/R)={Q:,}"
            )

            print(
                f"    TRUE:"
            )

            print(
                f"        K={K_true:,}"
                f" E={E_true:,}"
            )

            print(
                f"        (k,l)=({k_true},{l_true})"
            )

            print(
                f"        (a,b)=({a_true},{b_true})"
            )

            print(
                f"    K in brute range = "
                f"{K_MIN <= K_true <= K_MAX}"
            )

            print(
                f"    E in allowed range = "
                f"{0 <= E_true <= E_MAX}"
            )

            if result["solved"]:

                print()
                print("    RECOVERED:")

                print(
                    f"        K={result['K']:,}"
                    f" E={result['E']:,}"
                )

                print(
                    f"        (k,l)=("
                    f"{result['k']},"
                    f"{result['l']})"
                )

                print(
                    f"        (a,b)=("
                    f"{result['a']},"
                    f"{result['b']})"
                )

                print(
                    f"        p,q=("
                    f"{result['p']:,},"
                    f"{result['q']:,})"
                )

                print(
                    f"        tested a="
                    f"{result['tested_a']:,}"
                )

                print(
                    f"    solved = YES"
                )

                print(
                    f"    correct = "
                    f"{'YES' if correct else 'NO'}"
                )

            else:

                print()
                print(
                    "    solved = NO"
                )

            print(
                f"    K tests          = "
                f"{result['K_tests']:,}"
            )

            print(
                f"    divisor pairs    = "
                f"{result['pair_tests']:,}"
            )

            print(
                f"    solve_ab calls   = "
                f"{result['solve_ab_calls']:,}"
            )

            print(
                f"    time             = "
                f"{result['time']:.8f}s"
            )

        # =====================================================================
        # SCALE SUMMARY
        # =====================================================================

        print()
        print("-" * 100)
        print(
            f"{scale:.0e} SUMMARY"
        )
        print("-" * 100)

        print(
            f"    cases             = "
            f"{scale_cases}"
        )

        print(
            f"    K/E exposed       = "
            f"{scale_exposed}/{scale_cases}"
        )

        print(
            f"    solved            = "
            f"{scale_solved}/{scale_cases}"
        )

        print(
            f"    correct           = "
            f"{scale_correct}/{scale_cases}"
        )

        print(
            f"    max K             = "
            f"{max(scale_K_values):,}"
        )

        print(
            f"    max E             = "
            f"{max(scale_E_values):,}"
        )

        print(
            f"    mean K            = "
            f"{sum(scale_K_values) / len(scale_K_values):,.2f}"
        )

        print(
            f"    mean E            = "
            f"{sum(scale_E_values) / len(scale_E_values):,.2f}"
        )

        print(
            f"    mean K tests      = "
            f"{scale_K_tests / scale_cases:.2f}"
        )

        print(
            f"    mean pair tests   = "
            f"{scale_pair_tests / scale_cases:.2f}"
        )

        print(
            f"    mean solve_ab     = "
            f"{scale_ab_calls / scale_cases:.2f}"
        )

        print(
            f"    mean runtime      = "
            f"{scale_runtime / scale_cases:.8f}s"
        )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 100)
    print("GLOBAL SUMMARY")
    print("=" * 100)

    print(
        f"    total cases       = "
        f"{global_cases}"
    )

    print(
        f"    K/E exposed       = "
        f"{global_exposed}/{global_cases}"
    )

    print(
        f"    solved            = "
        f"{global_solved}/{global_cases}"
    )

    print(
        f"    correct           = "
        f"{global_correct}/{global_cases}"
    )

    print(
        f"    maximum K         = "
        f"{max(global_K_values):,}"
    )

    print(
        f"    maximum E         = "
        f"{max(global_E_values):,}"
    )

    print(
        f"    mean K            = "
        f"{sum(global_K_values) / len(global_K_values):,.2f}"
    )

    print(
        f"    mean E            = "
        f"{sum(global_E_values) / len(global_E_values):,.2f}"
    )

    print(
        f"    mean K tests      = "
        f"{global_K_tests / global_cases:.2f}"
    )

    print(
        f"    mean pair tests   = "
        f"{global_pair_tests / global_cases:.2f}"
    )

    print(
        f"    mean solve_ab     = "
        f"{global_ab_calls / global_cases:.2f}"
    )

    print(
        f"    mean runtime      = "
        f"{global_runtime / global_cases:.8f}s"
    )

    # =========================================================================
    # ALGEBRA
    # =========================================================================

    print()
    print("=" * 100)
    print("ALGEBRAIC IDENTITY")
    print("=" * 100)

    print()
    print("    p = r1*k + a")
    print("    q = r2*l + b")
    print()
    print("    n = p*q")
    print()
    print("      = (r1*k+a)(r2*l+b)")
    print()
    print("      = r1*r2*k*l")
    print("        + r1*k*b")
    print("        + r2*l*a")
    print("        + a*b")
    print()

    print("    Let:")
    print()
    print("        R = r1*r2")
    print("        K = k*l")
    print()

    print("    Then:")
    print()
    print(
        "        n = R*K + "
        "r1*k*b + r2*l*a + a*b"
    )
    print()

    print("    Therefore:")
    print()
    print(
        "        floor(n/R)"
        " = K + floor("
        "(r1*k*b + r2*l*a + a*b)/R"
        ")"
    )
    print()

    print("    Hence:")
    print()
    print(
        "        E = floor(n/R) - K"
    )
    print()

    print("    This identity is exact.")
    print()

    # =========================================================================
    # IMPORTANT OBSERVATION
    # =========================================================================

    print("=" * 100)
    print("IMPORTANT OBSERVATION")
    print("=" * 100)
    print()

    print(
        "The experiment does NOT prove that E <= 1001."
    )
    print()

    print(
        "It only tests whether E happens to remain in that"
    )
    print(
        "small range for the generated cases."
    )

    print()

    print(
        "The proposed search is therefore:"
    )

    print()

    print(
        "    Q = floor(n/(r1*r2))"
    )

    print(
        "    for candidate K:"
    )

    print(
        "        E = Q-K"
    )

    print(
        "        if 0 <= E <= E_MAX:"
    )

    print(
        "            lookup K -> (k,l)"
    )

    print(
        "            solve a,b"
    )

    print()

    print("=" * 100)
    print("FINISHED EXPERIMENT 79")
    print("=" * 100)


if __name__ == "__main__":
    main()