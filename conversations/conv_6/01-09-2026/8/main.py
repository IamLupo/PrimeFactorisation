import math
import random
import time


# ==========================================================================================
# START EXPERIMENT 80
# CORRECT-SCALE K/E SEARCH
# R = r1*r2 CHOSEN SO THAT n/R ~= 1000
# ==========================================================================================

SEED = 1511464998

SCALES = (
    10**9,
    10**12,
    10**16,
)

CASES_PER_SCALE = 20

TARGET_RATIO = 1000.0

# Proposed small-E search.
E_MAX = 1001

# Static K lookup.
K_TABLE_MAX = 5000

# Candidate primes for p,q.
# p and q are chosen so their product is close to the requested scale.
PRIME_MARGIN = 0.25

# For r1*r2 ~= n/1000:
# r1 and r2 are around sqrt(n/1000).
R_SEARCH_RADIUS = 5000


# ==========================================================================================
# PRIME UTILITIES
# ==========================================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3
    limit = math.isqrt(n)

    while d <= limit:
        if n % d == 0:
            return False
        d += 2

    return True


def previous_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n -= 1

    while n >= 2 and not is_prime(n):
        n -= 2

    return n


def next_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


def primes_near(center: int, radius: int = R_SEARCH_RADIUS):
    out = set()

    for d in range(radius + 1):
        lo = center - d
        hi = center + d

        if lo >= 2 and is_prime(lo):
            out.add(lo)

        if hi >= 2 and is_prime(hi):
            out.add(hi)

        if len(out) >= 32:
            break

    return sorted(out)


# ==========================================================================================
# STATIC K -> (k,l)
# ==========================================================================================

def build_k_lookup(k_max: int):
    table = {}

    total_pairs = 0

    for K in range(k_max + 1):
        pairs = []

        if K == 0:
            table[K] = [(0, 0)]
            continue

        limit = math.isqrt(K)

        for d in range(1, limit + 1):
            if K % d == 0:
                e = K // d

                pairs.append((d, e))

                if d != e:
                    pairs.append((e, d))

        pairs.sort()
        total_pairs += len(pairs)

        table[K] = pairs

    return table, total_pairs


# ==========================================================================================
# CONSTRUCT r1,r2
# ==========================================================================================

def choose_r1_r2(n: int):
    target_R = max(2, round(n / TARGET_RATIO))

    target_root = math.isqrt(target_R)

    candidates = primes_near(target_root)

    best = None

    for r1 in candidates:
        for r2 in candidates:
            R = r1 * r2

            ratio = n / R

            score = abs(ratio - TARGET_RATIO)

            item = (score, r1, r2, R, ratio)

            if best is None or item[0] < best[0]:
                best = item

    _, r1, r2, R, ratio = best

    return r1, r2, R, ratio


# ==========================================================================================
# FIND PRIME p,q
# ==========================================================================================

def choose_p_q(scale: int, rng: random.Random):
    root = math.isqrt(scale)

    # Pick nearby primes whose product is close to the requested scale.
    p_center = root + rng.randint(
        -int(root * PRIME_MARGIN),
        int(root * PRIME_MARGIN),
    )

    q_center = max(
        3,
        scale // max(3, p_center),
    )

    p = next_prime(max(3, p_center))
    q = next_prime(max(3, q_center))

    # Make sure they are distinct.
    if p == q:
        q = next_prime(q + 2)

    return p, q


# ==========================================================================================
# TRUE PARAMETERS
# ==========================================================================================

def derive_parameters(p: int, q: int, r1: int, r2: int):
    R = r1 * r2

    k = p // r1
    a = p % r1

    l = q // r2
    b = q % r2

    K = k * l

    Q = n_div_R = (p * q) // R
    E = Q - K

    return {
        "n": p * q,
        "R": R,
        "Q": Q,
        "K": K,
        "E": E,
        "k": k,
        "l": l,
        "a": a,
        "b": b,
    }


# ==========================================================================================
# solve_ab GIVEN K,k,l
# ==========================================================================================

def solve_ab(n: int, r1: int, r2: int, K: int, k: int, l: int):
    max_a = r1 - 1

    tested = 0

    base_p = r1 * k
    base_q = r2 * l

    for a in range(max_a + 1):
        tested += 1

        p_candidate = base_p + a

        if p_candidate <= 0:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = n // p_candidate

        b = q_candidate - base_q

        if 0 <= b < r2:
            if p_candidate * q_candidate == n:
                return a, b, tested

    return None, None, tested


# ==========================================================================================
# BRUTE-FORCE E
# ==========================================================================================

def solve_K_E(
    n: int,
    r1: int,
    r2: int,
    k_table,
    e_max: int = E_MAX,
):
    R = r1 * r2
    Q = n // R

    e_tests = 0
    pair_tests = 0
    solve_ab_calls = 0

    for E in range(e_max + 1):
        K = Q - E
        e_tests += 1

        if K < 0:
            continue

        pairs = k_table.get(K)

        if not pairs:
            continue

        for k, l in pairs:
            pair_tests += 1

            a, b, tested_a = solve_ab(
                n,
                r1,
                r2,
                K,
                k,
                l,
            )

            solve_ab_calls += 1

            if a is not None:
                return {
                    "solved": True,
                    "K": K,
                    "E": E,
                    "k": k,
                    "l": l,
                    "a": a,
                    "b": b,
                    "e_tests": e_tests,
                    "pair_tests": pair_tests,
                    "solve_ab_calls": solve_ab_calls,
                    "tested_a": tested_a,
                }

    return {
        "solved": False,
        "K": None,
        "E": None,
        "k": None,
        "l": None,
        "a": None,
        "b": None,
        "e_tests": e_tests,
        "pair_tests": pair_tests,
        "solve_ab_calls": solve_ab_calls,
        "tested_a": 0,
    }


# ==========================================================================================
# ONE CASE
# ==========================================================================================

def run_case(
    scale: int,
    rng: random.Random,
    k_table,
    case_index: int,
):
    p, q = choose_p_q(scale, rng)
    n = p * q

    r1, r2, R, ratio = choose_r1_r2(n)

    true_data = derive_parameters(p, q, r1, r2)

    # Exact sanity checks.
    assert true_data["K"] == true_data["k"] * true_data["l"]

    assert (
        true_data["n"]
        ==
        (
            r1 * true_data["k"] + true_data["a"]
        )
        *
        (
            r2 * true_data["l"] + true_data["b"]
        )
    )

    assert (
        true_data["E"]
        ==
        true_data["Q"] - true_data["K"]
    )

    t0 = time.perf_counter()

    result = solve_K_E(
        n,
        r1,
        r2,
        k_table,
        E_MAX,
    )

    elapsed = time.perf_counter() - t0

    correct = (
        result["solved"]
        and result["K"] == true_data["K"]
        and result["E"] == true_data["E"]
        and (
            result["k"] == true_data["k"]
            and result["l"] == true_data["l"]
        )
        and (
            result["a"] == true_data["a"]
            and result["b"] == true_data["b"]
        )
    )

    return {
        "scale": scale,
        "case": case_index,
        "p": p,
        "q": q,
        "n": n,
        "r1": r1,
        "r2": r2,
        "R": R,
        "ratio": ratio,
        "Q": true_data["Q"],
        "true_K": true_data["K"],
        "true_E": true_data["E"],
        "true_k": true_data["k"],
        "true_l": true_data["l"],
        "true_a": true_data["a"],
        "true_b": true_data["b"],
        "result": result,
        "correct": correct,
        "elapsed": elapsed,
    }


# ==========================================================================================
# MAIN
# ==========================================================================================

def main():
    rng = random.Random(SEED)

    print("=" * 100)
    print("START EXPERIMENT 80")
    print("CORRECT-SCALE K/E BRUTE-FORCE")
    print("=" * 100)

    print()
    print("configuration")
    print(f"    scales                = {[f'{x:.0e}' for x in SCALES]}")
    print(f"    cases / scale         = {CASES_PER_SCALE}")
    print(f"    target n/R            = {TARGET_RATIO}")
    print(f"    E brute range         = 0..{E_MAX}")
    print(f"    static K table max    = {K_TABLE_MAX}")
    print(f"    seed                  = {SEED}")

    print()
    print("=" * 100)
    print("BUILDING STATIC K LOOKUP")
    print("=" * 100)

    t_lookup = time.perf_counter()

    k_table, divisor_pairs = build_k_lookup(K_TABLE_MAX)

    lookup_time = time.perf_counter() - t_lookup

    print(f"    K entries             = {len(k_table):,}")
    print(f"    divisor pairs         = {divisor_pairs:,}")
    print(f"    build time             = {lookup_time:.6f}s")

    global_results = []

    for scale in SCALES:
        print()
        print("=" * 100)
        print(f"SCALE {scale:.0e}")
        print("=" * 100)

        scale_results = []

        exposed = 0
        solved = 0
        correct = 0

        max_E = 0
        max_K = 0
        max_ratio_error = 0.0

        total_e_tests = 0
        total_pair_tests = 0
        total_ab_calls = 0
        total_time = 0.0

        for i in range(1, CASES_PER_SCALE + 1):
            row = run_case(
                scale,
                rng,
                k_table,
                i,
            )

            scale_results.append(row)
            global_results.append(row)

            true_E = row["true_E"]
            true_K = row["true_K"]

            if 0 <= true_E <= E_MAX:
                exposed += 1

            if row["result"]["solved"]:
                solved += 1

            if row["correct"]:
                correct += 1

            max_E = max(max_E, true_E)
            max_K = max(max_K, true_K)

            max_ratio_error = max(
                max_ratio_error,
                abs(row["ratio"] - TARGET_RATIO),
            )

            total_e_tests += row["result"]["e_tests"]
            total_pair_tests += row["result"]["pair_tests"]
            total_ab_calls += row["result"]["solve_ab_calls"]
            total_time += row["elapsed"]

            print()
            print(f"CASE {i}/{CASES_PER_SCALE}")
            print(
                f"    n={row['n']:,}"
            )
            print(
                f"    p,q=({row['p']:,},{row['q']:,})"
            )
            print(
                f"    r1,r2=({row['r1']:,},{row['r2']:,})"
                f" R={row['R']:,}"
            )
            print(
                f"    n/R={row['ratio']:.9f}"
            )
            print(
                f"    Q=floor(n/R)={row['Q']:,}"
            )

            print()
            print("    TRUE:")
            print(
                f"        K={row['true_K']:,}"
                f" E={row['true_E']:,}"
            )
            print(
                f"        (k,l)=({row['true_k']:,},{row['true_l']:,})"
            )
            print(
                f"        (a,b)=({row['true_a']:,},{row['true_b']:,})"
            )

            print()
            print("    SEARCH:")
            print(
                f"        solved={row['result']['solved']}"
            )
            print(
                f"        K={row['result']['K']}"
                f" E={row['result']['E']}"
            )
            print(
                f"        (k,l)=({row['result']['k']},"
                f"{row['result']['l']})"
            )
            print(
                f"        (a,b)=({row['result']['a']},"
                f"{row['result']['b']})"
            )
            print(
                f"        E tests={row['result']['e_tests']:,}"
            )
            print(
                f"        divisor-pair tests={row['result']['pair_tests']:,}"
            )
            print(
                f"        solve_ab calls={row['result']['solve_ab_calls']:,}"
            )
            print(
                f"        solve_ab a-tests={row['result']['tested_a']:,}"
            )
            print(
                f"        correct={row['correct']}"
            )
            print(
                f"        time={row['elapsed']:.8f}s"
            )

        print()
        print("-" * 100)
        print(f"{scale:.0e} SUMMARY")
        print("-" * 100)

        print(f"    cases                 = {CASES_PER_SCALE}")
        print(f"    E <= {E_MAX}          = {exposed}/{CASES_PER_SCALE}")
        print(f"    solved                = {solved}/{CASES_PER_SCALE}")
        print(f"    correct               = {correct}/{CASES_PER_SCALE}")
        print(f"    maximum K             = {max_K:,}")
        print(f"    maximum E             = {max_E:,}")
        print(
            f"    mean n/R              = "
            f"{sum(x['ratio'] for x in scale_results) / len(scale_results):.6f}"
        )
        print(
            f"    max |n/R-1000|        = "
            f"{max_ratio_error:.6f}"
        )
        print(
            f"    mean E tests          = "
            f"{total_e_tests / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean pair tests       = "
            f"{total_pair_tests / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean solve_ab calls   = "
            f"{total_ab_calls / CASES_PER_SCALE:.2f}"
        )
        print(
            f"    mean runtime          = "
            f"{total_time / CASES_PER_SCALE:.8f}s"
        )

    # ======================================================================================
    # GLOBAL SUMMARY
    # ======================================================================================

    print()
    print("=" * 100)
    print("GLOBAL SUMMARY")
    print("=" * 100)

    total = len(global_results)
    exposed = sum(
        1
        for x in global_results
        if 0 <= x["true_E"] <= E_MAX
    )
    solved = sum(
        1
        for x in global_results
        if x["result"]["solved"]
    )
    correct = sum(
        1
        for x in global_results
        if x["correct"]
    )

    print(f"    total cases            = {total}")
    print(f"    E <= {E_MAX}            = {exposed}/{total}")
    print(f"    solved                 = {solved}/{total}")
    print(f"    correct                = {correct}/{total}")

    print(
        f"    maximum E             = "
        f"{max(x['true_E'] for x in global_results):,}"
    )

    print(
        f"    mean E                = "
        f"{sum(x['true_E'] for x in global_results) / total:.3f}"
    )

    print(
        f"    maximum |n/R-1000|    = "
        f"{max(abs(x['ratio'] - TARGET_RATIO) for x in global_results):.6f}"
    )

    print()
    print("=" * 100)
    print("ALGEBRAIC CHECKS")
    print("=" * 100)

    identity_failures = 0
    decomposition_failures = 0
    ratio_regime_failures = 0

    for x in global_results:
        n = x["n"]
        R = x["R"]

        # E = floor(n/R) - K
        if x["true_E"] != (n // R) - x["true_K"]:
            identity_failures += 1

        # p = r1*k+a and q = r2*l+b
        p_reconstructed = x["r1"] * x["true_k"] + x["true_a"]
        q_reconstructed = x["r2"] * x["true_l"] + x["true_b"]

        if p_reconstructed != x["p"] or q_reconstructed != x["q"]:
            decomposition_failures += 1

        # The important regime check.
        # We intentionally keep R close to n/1000.
        if abs(x["ratio"] - TARGET_RATIO) > 10.0:
            ratio_regime_failures += 1

    print(f"    E identity failures   = {identity_failures}")
    print(f"    decomposition failures= {decomposition_failures}")
    print(f"    ratio-regime failures = {ratio_regime_failures}")

    print()
    print("=" * 100)
    print("EXPERIMENT PURPOSE")
    print("=" * 100)

    print(
        """
The construction deliberately chooses:

    R = r1*r2 ~= n/1000

so that:

    n/R ~= 1000.

Then:

    Q = floor(n/R)

and because:

    E = Q - K

we test the proposed strategy:

    K = Q - E

with a small E search.

For every candidate K:

    static_lookup[K] -> possible (k,l)

and then:

    p = r1*k + a
    q = r2*l + b

with solve_ab() checking the remaining a,b space.

The experiment therefore isolates whether a small E range is sufficient
when r1 and r2 are chosen in the intended scale.
"""
    )

    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 80")
    print("=" * 100)


if __name__ == "__main__":
    main()