import math
import time
import sympy


# ============================================================
# START EXPERIMENT 91
# ============================================================
#
# CORRECT K/E SEARCH OPTIMIZATION
#
# IMPORTANT CORRECTION:
#
# Experiment 90 incorrectly used
#
#     D = n - R*K
#
# as though
#
#     D = p*q - R*K.
#
# That is only true when p*q == n.
#
# The correct relationship is:
#
#     T = p*q - R*K
#
#       = r1*k*b + r2*l*a + a*b
#
# and
#
#     E = floor(T / R).
#
# Therefore:
#
#     R*E <= T < R*(E+1)
#
# For fixed a:
#
#     T = (r1*k+a)*b + r2*l*a
#
# so we can derive an EXACT interval for b:
#
#     ceil((R*E-r2*l*a)/(r1*k+a))
#
#        <= b <=
#
#     floor((R*(E+1)-1-r2*l*a)/(r1*k+a))
#
#
# This experiment compares:
#
#   METHOD A:
#       complete a*b scan
#
#   METHOD B:
#       interval-bounded b scan
#
#   METHOD C:
#       interval-bounded scan + direct E verification
#       only on surviving b values
#
# None of the methods require p*q == n.
#
# The candidate sets must agree.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGETS = [
    10**9,
    10**10,
    10**12,
]

R_PAIRS = [
    (5, 7),
    (7, 11),
    (11, 13),
    (13, 17),
]

E_LIMIT_MULTIPLIER = 3.0

PRINT_EXAMPLES = 10


# ------------------------------------------------------------
# Generate semiprime
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
# Correct mathematical ceil division
# ------------------------------------------------------------

def ceil_div(a, b):

    return -((-a) // b)


# ------------------------------------------------------------
# Compute E
# ------------------------------------------------------------

def compute_E(k, l, a, b, r1, r2):

    R = r1 * r2

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    return c1 + c2 + c3


# ------------------------------------------------------------
# Build candidate record
# ------------------------------------------------------------

def make_candidate(
    k,
    l,
    a,
    b,
    K,
    E,
    r1,
    r2,
):

    p = r1 * k + a
    q = r2 * l + b

    return {
        "p": p,
        "q": q,
        "k": k,
        "l": l,
        "a": a,
        "b": b,
        "K": K,
        "E": E,
    }


# ------------------------------------------------------------
# METHOD A
#
# Baseline complete a*b scan.
# ------------------------------------------------------------

def baseline_candidates(n, r1, r2):

    R = r1 * r2
    Q = n // R

    e_limit = (
        int(
            E_LIMIT_MULTIPLIER
            * math.isqrt(Q)
        )
        + 20
    )

    out = {}

    pair_tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        for k in sympy.divisors(K):

            l = K // k

            for a in range(r1):

                for b in range(r2):

                    pair_tests += 1

                    E2 = compute_E(
                        k,
                        l,
                        a,
                        b,
                        r1,
                        r2,
                    )

                    if E2 != E:
                        continue

                    c = make_candidate(
                        k,
                        l,
                        a,
                        b,
                        K,
                        E,
                        r1,
                        r2,
                    )

                    out[(c["p"], c["q"])] = c

    elapsed = (
        time.perf_counter()
        - start
    )

    return list(out.values()), {
        "time": elapsed,
        "pair_tests": pair_tests,
    }


# ------------------------------------------------------------
# METHOD B
#
# Exact interval-bounded search.
#
# We use:
#
#     R*E <= T < R*(E+1)
#
# without assuming pq=n.
# ------------------------------------------------------------

def bounded_candidates(n, r1, r2):

    R = r1 * r2
    Q = n // R

    e_limit = (
        int(
            E_LIMIT_MULTIPLIER
            * math.isqrt(Q)
        )
        + 20
    )

    out = {}

    b_tests = 0
    a_tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        lower_T = R * E
        upper_T = R * (E + 1) - 1

        for k in sympy.divisors(K):

            l = K // k

            for a in range(r1):

                a_tests += 1

                coefficient = (
                    r1 * k + a
                )

                constant = (
                    r2 * l * a
                )

                # ------------------------------------------------
                # Solve:
                #
                # lower_T <= coefficient*b + constant
                #
                # and
                #
                # coefficient*b + constant <= upper_T
                # ------------------------------------------------

                b_min = ceil_div(
                    lower_T - constant,
                    coefficient,
                )

                b_max = (
                    upper_T - constant
                ) // coefficient

                # Intersect with:
                #
                #     0 <= b < r2
                #
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

                for b in range(
                    b_min,
                    b_max + 1,
                ):

                    b_tests += 1

                    # This should already guarantee E,
                    # but verify it to make the result
                    # completely explicit.
                    E2 = compute_E(
                        k,
                        l,
                        a,
                        b,
                        r1,
                        r2,
                    )

                    if E2 != E:
                        continue

                    c = make_candidate(
                        k,
                        l,
                        a,
                        b,
                        K,
                        E,
                        r1,
                        r2,
                    )

                    out[(c["p"], c["q"])] = c

    elapsed = (
        time.perf_counter()
        - start
    )

    return list(out.values()), {
        "time": elapsed,
        "a_tests": a_tests,
        "b_tests": b_tests,
    }


# ------------------------------------------------------------
# METHOD C
#
# Same exact interval search, but exploit the fact that
# the interval calculation itself guarantees E.
#
# No compute_E() call is needed inside the b loop.
# ------------------------------------------------------------

def direct_interval_candidates(
    n,
    r1,
    r2,
):

    R = r1 * r2
    Q = n // R

    e_limit = (
        int(
            E_LIMIT_MULTIPLIER
            * math.isqrt(Q)
        )
        + 20
    )

    out = {}

    a_tests = 0
    b_tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        lower_T = R * E
        upper_T = R * (E + 1) - 1

        for k in sympy.divisors(K):

            l = K // k

            for a in range(r1):

                a_tests += 1

                coefficient = (
                    r1 * k + a
                )

                constant = (
                    r2 * l * a
                )

                b_min = ceil_div(
                    lower_T - constant,
                    coefficient,
                )

                b_max = (
                    upper_T - constant
                ) // coefficient

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

                for b in range(
                    b_min,
                    b_max + 1,
                ):

                    b_tests += 1

                    c = make_candidate(
                        k,
                        l,
                        a,
                        b,
                        K,
                        E,
                        r1,
                        r2,
                    )

                    out[(c["p"], c["q"])] = c

    elapsed = (
        time.perf_counter()
        - start
    )

    return list(out.values()), {
        "time": elapsed,
        "a_tests": a_tests,
        "b_tests": b_tests,
    }


# ------------------------------------------------------------
# Candidate keys
# ------------------------------------------------------------

def candidate_keys(candidates):

    return {
        (
            c["p"],
            c["q"],
        )
        for c in candidates
    }


# ------------------------------------------------------------
# Exact factors only for final reporting
# ------------------------------------------------------------

def exact_count(candidates, n):

    return sum(
        1
        for c in candidates
        if c["p"] * c["q"] == n
    )


# ------------------------------------------------------------
# Check true candidate
# ------------------------------------------------------------

def true_present(
    candidates,
    p,
    q,
):

    keys = candidate_keys(
        candidates
    )

    return (
        (p, q) in keys
        or
        (q, p) in keys
    )


# ------------------------------------------------------------
# Print examples
# ------------------------------------------------------------

def print_examples(
    candidates,
    n,
):

    for c in candidates[
        :PRINT_EXAMPLES
    ]:

        product = (
            c["p"] * c["q"]
        )

        print(
            "    "
            f"p={c['p']} "
            f"q={c['q']} "
            f"K={c['K']} "
            f"E={c['E']} "
            f"D=pq-n={product-n}"
        )


# ------------------------------------------------------------
# One comparison
# ------------------------------------------------------------

def run_comparison(
    n,
    p_true,
    q_true,
    r1,
    r2,
):

    print()
    print(
        f"r1={r1} "
        f"r2={r2} "
        f"R={r1*r2}"
    )

    # --------------------------------------------------------
    # Baseline
    # --------------------------------------------------------

    baseline, sb = baseline_candidates(
        n,
        r1,
        r2,
    )

    # --------------------------------------------------------
    # Correct bounded method
    # --------------------------------------------------------

    bounded, sbd = bounded_candidates(
        n,
        r1,
        r2,
    )

    # --------------------------------------------------------
    # Direct interval method
    # --------------------------------------------------------

    direct, sd = direct_interval_candidates(
        n,
        r1,
        r2,
    )

    kb = candidate_keys(
        baseline
    )

    kbd = candidate_keys(
        bounded
    )

    kd = candidate_keys(
        direct
    )

    # --------------------------------------------------------
    # Counts
    # --------------------------------------------------------

    print()
    print("CANDIDATE COUNTS")

    print(
        f"    baseline = "
        f"{len(baseline):,}"
    )

    print(
        f"    bounded  = "
        f"{len(bounded):,}"
    )

    print(
        f"    direct   = "
        f"{len(direct):,}"
    )

    # --------------------------------------------------------
    # Correctness
    # --------------------------------------------------------

    print()
    print("CORRECTNESS")

    print(
        f"    baseline == bounded = "
        f"{kb == kbd}"
    )

    print(
        f"    baseline == direct  = "
        f"{kb == kd}"
    )

    print(
        f"    true in baseline    = "
        f"{true_present(baseline,p_true,q_true)}"
    )

    print(
        f"    true in bounded     = "
        f"{true_present(bounded,p_true,q_true)}"
    )

    print(
        f"    true in direct      = "
        f"{true_present(direct,p_true,q_true)}"
    )

    # --------------------------------------------------------
    # Exact counts
    # --------------------------------------------------------

    print()
    print("EXACT FACTORS IN CANDIDATE CLOUD")

    print(
        f"    baseline = "
        f"{exact_count(baseline,n)}"
    )

    print(
        f"    bounded  = "
        f"{exact_count(bounded,n)}"
    )

    print(
        f"    direct   = "
        f"{exact_count(direct,n)}"
    )

    # --------------------------------------------------------
    # Timing
    # --------------------------------------------------------

    print()
    print("TIMING")

    print(
        f"    baseline = "
        f"{sb['time']:.6f} s"
    )

    print(
        f"    bounded  = "
        f"{sbd['time']:.6f} s"
    )

    print(
        f"    direct   = "
        f"{sd['time']:.6f} s"
    )

    print()
    print("SPEEDUP")

    print(
        f"    bounded vs baseline = "
        f"{sb['time'] / max(sbd['time'],1e-12):.3f}x"
    )

    print(
        f"    direct vs baseline  = "
        f"{sb['time'] / max(sd['time'],1e-12):.3f}x"
    )

    # --------------------------------------------------------
    # Search work
    # --------------------------------------------------------

    print()
    print("SEARCH WORK")

    print(
        f"    baseline pair tests = "
        f"{sb['pair_tests']:,}"
    )

    print(
        f"    bounded a tests     = "
        f"{sbd['a_tests']:,}"
    )

    print(
        f"    bounded b tests     = "
        f"{sbd['b_tests']:,}"
    )

    print(
        f"    direct a tests      = "
        f"{sd['a_tests']:,}"
    )

    print(
        f"    direct b tests      = "
        f"{sd['b_tests']:,}"
    )

    # --------------------------------------------------------
    # Direct examples
    # --------------------------------------------------------

    print()
    print("DIRECT EXAMPLES")

    print_examples(
        direct,
        n,
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("CORRECT RESIDUAL-INTERVAL K/E SEARCH")
    print("EXPERIMENT 91")
    print("=" * 78)

    for target in TARGETS:

        print()
        print("=" * 78)
        print(
            f"TARGET ~ {target:,}"
        )
        print("=" * 78)

        p_true, q_true, n = (
            generate_semiprime(target)
        )

        print()
        print(
            f"n      = {n}"
        )

        print(
            f"true p = {p_true}"
        )

        print(
            f"true q = {q_true}"
        )

        print(
            f"|p-q|  = "
            f"{abs(p_true-q_true)}"
        )

        for r1, r2 in R_PAIRS:

            run_comparison(
                n,
                p_true,
                q_true,
                r1,
                r2,
            )

    print()
    print("=" * 78)
    print("EXPERIMENT 91 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 91
# ============================================================
