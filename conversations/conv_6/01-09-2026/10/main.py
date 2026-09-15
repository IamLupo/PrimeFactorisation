import math
import random
import time


# ==========================================================================================
# START EXPERIMENT 82
# C1/C2/C3 INFORMATION TEST FOR solve_ab
# ==========================================================================================

SEED = 1511464998

CASES_PER_SCALE = 10

SCALES = [
    10**9,
    10**12,
    10**16,
]

TARGET_N_OVER_R = 100.0

K_TABLE_MAX = 5000


# ==========================================================================================
# PRIMALITY
# ==========================================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def next_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


def previous_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n -= 1

    while n >= 2:
        if is_prime(n):
            return n
        n -= 2

    return 2


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def random_prime_near(target: int, rng: random.Random, spread: int) -> int:
    target = max(3, int(target))

    for _ in range(10000):
        candidate = target + rng.randint(-spread, spread)

        if candidate < 2:
            continue

        if candidate % 2 == 0:
            candidate += 1

        if is_prime(candidate):
            return candidate

    return next_prime(target)


def generate_semiprime(scale: int, rng: random.Random):
    root = math.isqrt(scale)

    spread = max(1000, root // 20)

    p = random_prime_near(root, rng, spread)
    q = random_prime_near(root, rng, spread)

    while q == p:
        q = random_prime_near(root + 100, rng, spread)

    return p, q


# ==========================================================================================
# CHOOSE r1,r2 SUCH THAT R ~= n/100
# ==========================================================================================

def choose_r1_r2(n: int, rng: random.Random):
    target_R = max(4, round(n / TARGET_N_OVER_R))

    target_root = math.isqrt(target_R)

    best = None
    best_error = None

    # Search around sqrt(R).
    for delta1 in range(-100, 101):
        t1 = target_root + delta1

        if t1 < 3:
            continue

        r1 = previous_prime(t1)

        if r1 < 3:
            continue

        target_r2 = max(3, target_R // r1)

        for delta2 in range(-30, 31):
            t2 = target_r2 + delta2

            if t2 < 3:
                continue

            r2 = previous_prime(t2)

            if r2 < 3:
                continue

            R = r1 * r2
            ratio = n / R
            error = abs(ratio - TARGET_N_OVER_R)

            if best_error is None or error < best_error:
                best_error = error
                best = (r1, r2)

    if best is None:
        r1 = previous_prime(target_root)
        r2 = previous_prime(max(3, target_R // r1))
        best = (r1, r2)

    return best


# ==========================================================================================
# PARAMETER CALCULATION
# ==========================================================================================

def calculate_parameters(p: int, q: int, r1: int, r2: int):
    n = p * q
    R = r1 * r2

    # Euclidean decompositions.
    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

    K = k * l

    Q, remainder = divmod(n, R)

    # Exact expanded cross-term.
    cross_term = (
        r1 * k * b
        + r2 * l * a
        + a * b
    )

    # This SHOULD equal the remainder because:
    #
    # n = R*K + cross_term.
    #
    # But we explicitly verify the complete identity instead of making
    # an assumption about how an earlier script constructed the values.
    expanded = (
        R * K
        + r1 * k * b
        + r2 * l * a
        + a * b
    )

    assert expanded == n
    assert 0 <= a < r1
    assert 0 <= b < r2
    assert Q == K + cross_term // R
    assert remainder == cross_term % R

    E = Q - K

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
        "p": p,
        "q": q,
        "r1": r1,
        "r2": r2,
        "remainder": remainder,
        "cross_term": cross_term,
    }


# ==========================================================================================
# STATIC K -> (k,l) TABLE
# ==========================================================================================

def build_k_table(max_K: int):
    table = {}

    divisor_pairs = 0

    for K in range(max_K + 1):
        pairs = []

        if K == 0:
            pairs.append((0, 0))
        else:
            d = 1

            while d * d <= K:
                if K % d == 0:
                    e = K // d

                    pairs.append((d, e))

                    if d != e:
                        pairs.append((e, d))

                    divisor_pairs += 1 if d == e else 2

                d += 1

        table[K] = tuple(sorted(set(pairs)))

    return table, divisor_pairs


# ==========================================================================================
# C VALUES
# ==========================================================================================

def c123(
    a: int,
    b: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
):
    R = r1 * r2

    # These are the natural normalized cross terms.
    #
    # c1 = k*b/r2
    # c2 = l*a/r1
    # c3 = a*b/(r1*r2)

    c1_num = k * b
    c1_den = r2

    c2_num = l * a
    c2_den = r1

    c3_num = a * b
    c3_den = R

    # Common denominator R.
    total_num = (
        r1 * k * b
        + r2 * l * a
        + a * b
    )

    total_den = R

    return (
        c1_num,
        c1_den,
        c2_num,
        c2_den,
        c3_num,
        c3_den,
        total_num,
        total_den,
    )


# ==========================================================================================
# REFERENCE solve_ab
# ==========================================================================================

def solve_ab_bruteforce(
    n: int,
    r1: int,
    r2: int,
    K: int,
    E: int,
    k: int,
    l: int,
):
    R = r1 * r2

    candidates = []
    tested = 0

    for a in range(r1):
        tested += 1

        p_candidate = r1 * k + a

        if p_candidate <= 1:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = n // p_candidate

        b_candidate = q_candidate - r2 * l

        if not (0 <= b_candidate < r2):
            continue

        if k * l != K:
            continue

        if n // R - K != E:
            continue

        if p_candidate * q_candidate != n:
            continue

        candidates.append(
            (
                a,
                b_candidate,
                p_candidate,
                q_candidate,
            )
        )

    return candidates, tested


# ==========================================================================================
# C1+C2+C3 TEST
# ==========================================================================================

def c123_check(
    n: int,
    r1: int,
    r2: int,
    K: int,
    E: int,
    k: int,
    l: int,
    a: int,
):
    R = r1 * r2

    p_candidate = r1 * k + a

    if p_candidate <= 1:
        return False, None

    if n % p_candidate != 0:
        return False, None

    q_candidate = n // p_candidate

    b_candidate = q_candidate - r2 * l

    if not (0 <= b_candidate < r2):
        return False, None

    (
        c1_num,
        c1_den,
        c2_num,
        c2_den,
        c3_num,
        c3_den,
        total_num,
        total_den,
    ) = c123(
        a,
        b_candidate,
        r1,
        r2,
        k,
        l,
    )

    # E condition.
    c_floor = total_num // total_den

    if c_floor != E:
        return False, None

    return True, (
        a,
        b_candidate,
        p_candidate,
        q_candidate,
    )


# ==========================================================================================
# CASE
# ==========================================================================================

def run_case(scale: int, case_no: int, rng: random.Random):
    p, q = generate_semiprime(scale, rng)

    n = p * q

    r1, r2 = choose_r1_r2(n, rng)

    params = calculate_parameters(
        p,
        q,
        r1,
        r2,
    )

    R = params["R"]
    K = params["K"]
    E = params["E"]
    k = params["k"]
    l = params["l"]
    a_true = params["a"]
    b_true = params["b"]

    print(f"CASE {case_no}/{CASES_PER_SCALE}")
    print(f"    n={n:,}")
    print(f"    p,q=({p:,},{q:,})")
    print(f"    r1,r2=({r1:,},{r2:,})")
    print(f"    R={R:,}")
    print(f"    n/R={n/R:.12f}")

    print()
    print("    TRUE:")
    print(f"        K={K:,}")
    print(f"        E={E:,}")
    print(f"        (k,l)=({k:,},{l:,})")
    print(f"        (a,b)=({a_true:,},{b_true:,})")

    # --------------------------------------------------------------------------
    # C values for true solution
    # --------------------------------------------------------------------------

    (
        c1_num,
        c1_den,
        c2_num,
        c2_den,
        c3_num,
        c3_den,
        total_num,
        total_den,
    ) = c123(
        a_true,
        b_true,
        r1,
        r2,
        k,
        l,
    )

    print()
    print("    C VALUES:")
    print(
        f"        c1 = {c1_num}/{c1_den}"
        f" = {c1_num/c1_den:.12f}"
    )
    print(
        f"        c2 = {c2_num}/{c2_den}"
        f" = {c2_num/c2_den:.12f}"
    )
    print(
        f"        c3 = {c3_num}/{c3_den}"
        f" = {c3_num/c3_den:.12f}"
    )
    print(
        f"        c1+c2+c3 = "
        f"{total_num}/{total_den}"
        f" = {total_num/total_den:.12f}"
    )
    print(
        f"        floor(c1+c2+c3) = "
        f"{total_num//total_den}"
    )

    # --------------------------------------------------------------------------
    # Brute reference
    # --------------------------------------------------------------------------

    t0 = time.perf_counter()

    brute_candidates, brute_tested = solve_ab_bruteforce(
        n,
        r1,
        r2,
        K,
        E,
        k,
        l,
    )

    brute_time = time.perf_counter() - t0

    # --------------------------------------------------------------------------
    # C123 pass
    # --------------------------------------------------------------------------

    t1 = time.perf_counter()

    c123_candidates = []
    c123_tested = 0

    for a in range(r1):
        c123_tested += 1

        ok, candidate = c123_check(
            n,
            r1,
            r2,
            K,
            E,
            k,
            l,
            a,
        )

        if ok:
            c123_candidates.append(candidate)

    c123_time = time.perf_counter() - t1

    true_brute = any(
        x[0] == a_true and x[1] == b_true
        for x in brute_candidates
    )

    true_c123 = any(
        x[0] == a_true and x[1] == b_true
        for x in c123_candidates
    )

    print()
    print("    SOLVER RESULTS:")
    print(
        f"        brute candidates       = "
        f"{len(brute_candidates)}"
    )
    print(
        f"        brute a tested        = "
        f"{brute_tested:,}"
    )
    print(
        f"        brute time            = "
        f"{brute_time:.8f}s"
    )

    print(
        f"        c123 candidates       = "
        f"{len(c123_candidates)}"
    )
    print(
        f"        c123 a tested        = "
        f"{c123_tested:,}"
    )
    print(
        f"        c123 time             = "
        f"{c123_time:.8f}s"
    )

    print(
        f"        true found brute      = "
        f"{true_brute}"
    )

    print(
        f"        true found c123       = "
        f"{true_c123}"
    )

    print(
        f"        factorization valid   = "
        f"{p*q == n}"
    )

    print()

    return {
        "scale": scale,
        "n": n,
        "R": R,
        "K": K,
        "E": E,
        "r1": r1,
        "r2": r2,
        "k": k,
        "l": l,
        "a": a_true,
        "b": b_true,
        "brute_candidates": len(brute_candidates),
        "c123_candidates": len(c123_candidates),
        "brute_time": brute_time,
        "c123_time": c123_time,
        "true_found": true_c123,
    }


# ==========================================================================================
# MAIN
# ==========================================================================================

def main():
    rng = random.Random(SEED)

    print("=" * 100)
    print("START EXPERIMENT 82")
    print("C1/C2/C3 INFORMATION TEST FOR solve_ab")
    print("=" * 100)

    print()
    print("configuration")
    print(
        f"    scales                = "
        f"{[f'{x:.0e}' for x in SCALES]}"
    )
    print(
        f"    cases / scale         = "
        f"{CASES_PER_SCALE}"
    )
    print(
        f"    target n/R            = "
        f"{TARGET_N_OVER_R}"
    )
    print(
        f"    seed                  = "
        f"{SEED}"
    )

    print()
    print("=" * 100)
    print("BUILDING STATIC K LOOKUP")
    print("=" * 100)

    t0 = time.perf_counter()

    k_table, divisor_pairs = build_k_table(
        K_TABLE_MAX
    )

    table_time = time.perf_counter() - t0

    print(
        f"K entries               = "
        f"{len(k_table):,}"
    )

    print(
        f"divisor pairs           = "
        f"{divisor_pairs:,}"
    )

    print(
        f"build time              = "
        f"{table_time:.6f}s"
    )

    all_results = []

    for scale in SCALES:
        print()
        print("=" * 100)
        print(f"SCALE {scale:.0e}")
        print("=" * 100)
        print()

        results = []

        for case_no in range(1, CASES_PER_SCALE + 1):
            result = run_case(
                scale,
                case_no,
                rng,
            )

            results.append(result)
            all_results.append(result)

        true_count = sum(
            1 for x in results
            if x["true_found"]
        )

        print("-" * 100)
        print(f"{scale:.0e} SUMMARY")
        print("-" * 100)

        print(
            f"    cases               = "
            f"{len(results)}"
        )

        print(
            f"    true c123 recovered = "
            f"{true_count}/{len(results)}"
        )

        print(
            f"    mean brute time     = "
            f"{sum(x['brute_time'] for x in results) / len(results):.8f}s"
        )

        print(
            f"    mean c123 time      = "
            f"{sum(x['c123_time'] for x in results) / len(results):.8f}s"
        )

        print()

    # ==========================================================================
    # GLOBAL
    # ==========================================================================

    true_count = sum(
        1 for x in all_results
        if x["true_found"]
    )

    print("=" * 100)
    print("GLOBAL SUMMARY")
    print("=" * 100)

    print(
        f"    total cases          = "
        f"{len(all_results)}"
    )

    print(
        f"    true c123 recovered  = "
        f"{true_count}/{len(all_results)}"
    )

    print()
    print("=" * 100)
    print("ALGEBRAIC IDENTITIES")
    print("=" * 100)
    print()

    print("    p = r1*k + a")
    print("    q = r2*l + b")
    print()
    print("    n = (r1*k+a)(r2*l+b)")
    print()
    print("    n = R*K + r1*k*b + r2*l*a + a*b")
    print()
    print("    c1 = k*b/r2")
    print("    c2 = l*a/r1")
    print("    c3 = a*b/R")
    print()
    print("    c1+c2+c3")
    print("        = (n-R*K)/R")
    print()
    print("    E = floor(c1+c2+c3)")
    print()
    print("The experiment therefore checks whether the C1/C2/C3")
    print("representation gives useful pruning beyond the known")
    print("values of K and E.")
    print()

    print("=" * 100)
    print("FINISHED EXPERIMENT 82")
    print("=" * 100)


if __name__ == "__main__":
    main()