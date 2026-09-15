import math
import time
import sympy


# ============================================================
# START EXPERIMENT 89
# ============================================================
#
# K/E RESIDUAL BOUND + MODULAR UNIQUENESS
#
# Main claim being tested:
#
#   If a candidate satisfies
#
#       floor(p*q / R) = floor(n / R)
#
#   then
#
#       |p*q - n| < R.
#
# Therefore, if we additionally know
#
#       p*q == n (mod M)
#
# and
#
#       M > R,
#
# then necessarily
#
#       p*q = n.
#
#
# This experiment tests that relationship over MANY random
# semiprimes and multiple small (r1,r2) pairs.
#
#
# We DO NOT use p*q == n during modular filtering.
#
# Exact equality is checked ONLY after the modular experiment
# has finished.
#
# ============================================================


# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

CASE_COUNT = 20

TARGETS = [
    10**9,
    10**10,
    10**11,
    10**12,
]

# Initial r values to test.
#
# We deliberately use several pairs rather than one.
R_PAIRS = [
    (5, 7),
    (7, 11),
    (11, 13),
    (13, 17),
    (17, 19),
]

# E search multiplier.
#
# This needs to be large enough to contain the true E.
E_LIMIT_MULTIPLIER = 3.0

# Small primes from which we construct M.
MODULI = [
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
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
]

# How many examples to print per case.
PRINT_CASES = 5

# How many surviving candidates to print.
PRINT_SURVIVORS = 20


# ------------------------------------------------------------
# Semiprime generation
# ------------------------------------------------------------

def generate_semiprime(target):

    root = math.isqrt(target)

    p = int(
        sympy.randprime(
            max(3, int(root * 0.85)),
            int(root * 1.05),
        )
    )

    q_est = target // p

    q = int(
        sympy.randprime(
            max(3, int(q_est * 0.90)),
            int(q_est * 1.10) + 100,
        )
    )

    return p, q, p * q


# ------------------------------------------------------------
# Carry calculation
# ------------------------------------------------------------

def compute_state(k, l, a, b, r1, r2):

    R = r1 * r2

    K = k * l

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    E = c1 + c2 + c3

    return K, E, c1, c2, c3


# ------------------------------------------------------------
# Initial candidate generation
# ------------------------------------------------------------

def generate_candidates(n, r1, r2):

    R = r1 * r2
    Q = n // R

    e_limit = (
        int(
            E_LIMIT_MULTIPLIER
            * math.isqrt(Q)
        )
        + 20
    )

    candidates = {}

    divisor_pairs = 0
    ab_tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        divisors = sympy.divisors(K)

        for k in divisors:

            l = K // k

            divisor_pairs += 1

            for a in range(r1):

                for b in range(r2):

                    ab_tests += 1

                    K2, E2, c1, c2, c3 = compute_state(
                        k,
                        l,
                        a,
                        b,
                        r1,
                        r2,
                    )

                    if K2 != K:
                        continue

                    if E2 != E:
                        continue

                    p = r1 * k + a
                    q = r2 * l + b

                    product = p * q

                    key = (p, q)

                    if key in candidates:
                        continue

                    candidates[key] = {
                        "p": p,
                        "q": q,
                        "product": product,
                        "D": product - n,

                        "k": k,
                        "l": l,
                        "a": a,
                        "b": b,

                        "K": K,
                        "E": E,

                        "c1": c1,
                        "c2": c2,
                        "c3": c3,
                    }

    elapsed = time.perf_counter() - start

    stats = {
        "Q": Q,
        "E_limit": e_limit,
        "divisor_pairs": divisor_pairs,
        "ab_tests": ab_tests,
        "time": elapsed,
    }

    return list(candidates.values()), stats


# ------------------------------------------------------------
# Verify the fundamental quotient bound
# ------------------------------------------------------------

def verify_bound(candidates, n, R):

    failures = []

    maximum = 0
    minimum = None

    for c in candidates:

        D = abs(c["D"])

        if D > maximum:
            maximum = D

        if minimum is None or D < minimum:
            minimum = D

        # Fundamental bound:
        #
        # |p*q-n| < R
        #
        if D >= R:
            failures.append(c)

    return {
        "ok": len(failures) == 0,
        "failures": failures,
        "min_abs_D": minimum,
        "max_abs_D": maximum,
    }


# ------------------------------------------------------------
# Build smallest modular prefix whose PRODUCT > R
# ------------------------------------------------------------

def build_modulus_prefix(R):

    M = 1
    used = []

    for prime in MODULI:

        M *= prime
        used.append(prime)

        if M > R:
            break

    return M, used


# ------------------------------------------------------------
# Modular filtering
# ------------------------------------------------------------

def modular_filter(candidates, moduli):

    survivors = []

    for c in candidates:

        p = c["p"]
        q = c["q"]

        ok = True

        for s in moduli:

            # DO NOT calculate p*q here.
            #
            # Use:
            #
            #     p*q mod s
            #
            # directly.
            #
            lhs = (
                (p % s)
                * (q % s)
            ) % s

            rhs = n_mod_global % s

            if lhs != rhs:

                ok = False
                break

        if ok:
            survivors.append(c)

    return survivors


# ------------------------------------------------------------
# Same filter but using D directly.
#
# This is only used to confirm equivalence.
# ------------------------------------------------------------

def modular_filter_residual(candidates, moduli):

    survivors = []

    for c in candidates:

        D = c["D"]

        ok = True

        for s in moduli:

            if D % s != 0:

                ok = False
                break

        if ok:
            survivors.append(c)

    return survivors


# ------------------------------------------------------------
# Determine minimum modulus prefix that forces D = 0
# ------------------------------------------------------------

def find_first_forcing_prefix(candidates, R):

    M = 1
    prefix = []

    for s in MODULI:

        M *= s
        prefix.append(s)

        survivors = [
            c
            for c in candidates
            if all(
                c["D"] % x == 0
                for x in prefix
            )
        ]

        if not survivors:
            return {
                "M": M,
                "moduli": prefix,
                "survivors": [],
            }

        all_zero = all(
            c["D"] == 0
            for c in survivors
        )

        if all_zero:

            return {
                "M": M,
                "moduli": prefix,
                "survivors": survivors,
            }

    return {
        "M": M,
        "moduli": prefix,
        "survivors": survivors,
    }


# ------------------------------------------------------------
# Exact check
# ------------------------------------------------------------

def exact_candidates(candidates):

    return [
        c
        for c in candidates
        if c["D"] == 0
    ]


# ------------------------------------------------------------
# True candidate check
# ------------------------------------------------------------

def true_alive(candidates, p, q):

    for c in candidates:

        if (
            (
                c["p"] == p
                and c["q"] == q
            )
            or
            (
                c["p"] == q
                and c["q"] == p
            )
        ):
            return True

    return False


# ------------------------------------------------------------
# Run one r pair
# ------------------------------------------------------------

def run_r_pair(n, p_true, q_true, r1, r2):

    global n_mod_global

    n_mod_global = n

    R = r1 * r2

    print()
    print(
        f"r1={r1} "
        f"r2={r2} "
        f"R={R}"
    )

    # --------------------------------------------------------
    # Generate cloud
    # --------------------------------------------------------

    start = time.perf_counter()

    candidates, stats = generate_candidates(
        n,
        r1,
        r2,
    )

    generation_time = (
        time.perf_counter()
        - start
    )

    print(
        f"    candidates = "
        f"{len(candidates):,}"
    )

    print(
        f"    E limit    = "
        f"{stats['E_limit']:,}"
    )

    print(
        f"    generation = "
        f"{generation_time:.6f} s"
    )

    # --------------------------------------------------------
    # Check |D| < R
    # --------------------------------------------------------

    bound = verify_bound(
        candidates,
        n,
        R,
    )

    print(
        f"    max |D|    = "
        f"{bound['max_abs_D']}"
    )

    print(
        f"    |D| < R    = "
        f"{bound['ok']}"
    )

    if not bound["ok"]:

        print()
        print(
            "    !!! BOUND FAILURE !!!"
        )

        return {
            "bound_ok": False,
            "forcing_ok": False,
            "initial": len(candidates),
        }

    # --------------------------------------------------------
    # Build modulus
    # --------------------------------------------------------

    M, used_moduli = build_modulus_prefix(R)

    print(
        f"    modulus M  = "
        f"{M}"
    )

    print(
        f"    moduli     = "
        f"{used_moduli}"
    )

    print(
        f"    M > R      = "
        f"{M > R}"
    )

    # --------------------------------------------------------
    # Modular filtering
    # --------------------------------------------------------

    start = time.perf_counter()

    modular_survivors = candidates

    for s in used_moduli:

        before = len(modular_survivors)

        modular_survivors = [
            c
            for c in modular_survivors
            if (
                (
                    (c["p"] % s)
                    *
                    (c["q"] % s)
                ) % s
            ) == (n % s)
        ]

        after = len(modular_survivors)

        print(
            f"        mod {s:<3} "
            f"{before:>6,} -> "
            f"{after:>6,}"
        )

    modular_time = (
        time.perf_counter()
        - start
    )

    # --------------------------------------------------------
    # Exact results
    # --------------------------------------------------------

    exact = exact_candidates(
        modular_survivors
    )

    exact_count = len(exact)

    true_survives = true_alive(
        modular_survivors,
        p_true,
        q_true,
    )

    forcing_ok = (
        len(modular_survivors) == exact_count
        and exact_count > 0
        and true_survives
    )

    print()
    print(
        f"    final survivors = "
        f"{len(modular_survivors):,}"
    )

    print(
        f"    exact survivors = "
        f"{exact_count:,}"
    )

    print(
        f"    true alive      = "
        f"{true_survives}"
    )

    print(
        f"    modular time    = "
        f"{modular_time:.9f} s"
    )

    # --------------------------------------------------------
    # Residuals
    # --------------------------------------------------------

    if modular_survivors:

        print()
        print(
            "    remaining residuals:"
        )

        for c in modular_survivors[
            :PRINT_SURVIVORS
        ]:

            print(
                f"        D={c['D']} "
                f"p={c['p']} "
                f"q={c['q']}"
            )

    # --------------------------------------------------------
    # Check direct theorem-style forcing
    # --------------------------------------------------------

    forcing = find_first_forcing_prefix(
        candidates,
        R,
    )

    print()
    print(
        "    first zero-forcing prefix:"
    )

    print(
        f"        M      = "
        f"{forcing['M']}"
    )

    print(
        f"        moduli = "
        f"{forcing['moduli']}"
    )

    print(
        f"        survivors = "
        f"{len(forcing['survivors']):,}"
    )

    all_zero = all(
        c["D"] == 0
        for c in forcing["survivors"]
    )

    print(
        f"        all D=0  = "
        f"{all_zero}"
    )

    return {
        "bound_ok": True,
        "forcing_ok": forcing_ok,
        "initial": len(candidates),
        "final": len(modular_survivors),
        "exact": exact_count,
        "true_alive": true_survives,
        "R": R,
        "M": M,
        "moduli": used_moduli,
        "max_abs_D": bound["max_abs_D"],
        "generation_time": generation_time,
        "modular_time": modular_time,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("K/E RESIDUAL BOUND + MODULAR UNIQUENESS")
    print("EXPERIMENT 89")
    print("=" * 78)

    total_cases = 0
    bound_failures = 0
    forcing_failures = 0

    for target in TARGETS:

        print()
        print("=" * 78)
        print(
            f"TARGET ~ {target:,}"
        )
        print("=" * 78)

        for case_index in range(CASE_COUNT):

            p_true, q_true, n = (
                generate_semiprime(target)
            )

            total_cases += 1

            print()
            print(
                f"CASE {case_index + 1}"
            )

            print(
                f"    n={n}"
            )

            print(
                f"    true p={p_true}"
            )

            print(
                f"    true q={q_true}"
            )

            print(
                f"    |p-q|="
                f"{abs(p_true-q_true)}"
            )

            # Only print a few pairs per case to prevent
            # gigantic output.
            pairs = R_PAIRS

            case_forcing_ok = True

            for r1, r2 in pairs:

                result = run_r_pair(
                    n,
                    p_true,
                    q_true,
                    r1,
                    r2,
                )

                if not result["bound_ok"]:

                    bound_failures += 1

                if not result["forcing_ok"]:

                    forcing_failures += 1
                    case_forcing_ok = False

            print()
            print(
                f"CASE RESULT: "
                f"{'PASS' if case_forcing_ok else 'FAIL'}"
            )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("EXPERIMENT 89 SUMMARY")
    print("=" * 78)

    print(
        f"total r-pair tests = "
        f"{total_cases * len(R_PAIRS):,}"
    )

    print(
        f"bound failures     = "
        f"{bound_failures:,}"
    )

    print(
        f"forcing failures   = "
        f"{forcing_failures:,}"
    )

    if (
        bound_failures == 0
        and forcing_failures == 0
    ):

        print()
        print(
            "RESULT: ALL TESTS PASSED"
        )

        print()
        print(
            "Every candidate satisfied:"
        )

        print(
            "    |p*q-n| < R"
        )

        print(
            "and the tested modular product M"
        )

        print(
            "was sufficient to leave only D=0."
        )

    else:

        print()
        print(
            "RESULT: SOME TESTS FAILED"
        )

    print()
    print("=" * 78)
    print("EXPERIMENT COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 89
# ============================================================
