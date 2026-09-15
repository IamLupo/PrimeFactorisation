import math
import time
import sympy


# ============================================================
# START EXPERIMENT 90
# ============================================================
#
# K/E CANDIDATE GENERATION:
#
#      BASELINE vs BOUNDED vs DIRECT-DIVISION
#
# We use:
#
#   n = (r1*k + a)(r2*l + b)
#
# and
#
#   K = k*l
#
#   E = c1+c2+c3
#
# with
#
#   floor(n / R) = K + E
#
# where
#
#   R = r1*r2
#
#
# For fixed K,k,l:
#
#   D = n - R*K
#
# and therefore
#
#   D = r1*k*b + r2*l*a + a*b
#
#            = (r1*k+a)b + r2*l*a
#
#
# Therefore:
#
#   b = (D-r2*l*a)/(r1*k+a)
#
#
# We compare:
#
#   METHOD A:
#       scan all a,b
#
#   METHOD B:
#       use residual bounds on a,b
#
#   METHOD C:
#       scan only a and derive b by division
#
#
# The candidate sets MUST agree.
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
# Carry E
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
# METHOD A
#
# Full a*b search.
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

    tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        for k in sympy.divisors(K):

            l = K // k

            for a in range(r1):

                for b in range(r2):

                    tests += 1

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

                    p = r1 * k + a
                    q = r2 * l + b

                    out[(p, q)] = {
                        "p": p,
                        "q": q,
                        "K": K,
                        "E": E,
                        "k": k,
                        "l": l,
                        "a": a,
                        "b": b,
                    }

    elapsed = time.perf_counter() - start

    return list(out.values()), {
        "time": elapsed,
        "tests": tests,
    }


# ------------------------------------------------------------
# METHOD B
#
# Residual-bounded search.
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

    tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        # Since the true candidate satisfies:
        #
        #     K = k*l
        #
        # try every divisor.
        for k in sympy.divisors(K):

            l = K // k

            # ------------------------------------------------
            # Residual:
            #
            # D = n - R*K
            #
            # D = r1*k*b + r2*l*a + a*b
            #
            # All terms are nonnegative.
            # ------------------------------------------------

            D = n - R * K

            if D < 0:
                continue

            # ------------------------------------------------
            # Since
            #
            #     D >= r2*l*a
            #
            # we have
            #
            #     a <= D/(r2*l)
            # ------------------------------------------------

            a_max = min(
                r1 - 1,
                D // (r2 * l)
                if l > 0
                else r1 - 1,
            )

            # ------------------------------------------------
            # Similarly:
            #
            #     D >= r1*k*b
            #
            # so
            #
            #     b <= D/(r1*k)
            # ------------------------------------------------

            b_max = min(
                r2 - 1,
                D // (r1 * k)
                if k > 0
                else r2 - 1,
            )

            for a in range(a_max + 1):

                for b in range(b_max + 1):

                    tests += 1

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

                    p = r1 * k + a
                    q = r2 * l + b

                    out[(p, q)] = {
                        "p": p,
                        "q": q,
                        "K": K,
                        "E": E,
                        "k": k,
                        "l": l,
                        "a": a,
                        "b": b,
                    }

    elapsed = time.perf_counter() - start

    return list(out.values()), {
        "time": elapsed,
        "tests": tests,
    }


# ------------------------------------------------------------
# METHOD C
#
# Direct division.
#
# For every a:
#
#     b = (D-r2*l*a)/(r1*k+a)
#
# ------------------------------------------------------------

def direct_candidates(n, r1, r2):

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

    tests = 0

    start = time.perf_counter()

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        for k in sympy.divisors(K):

            l = K // k

            D = n - R * K

            if D < 0:
                continue

            # ------------------------------------------------
            # We know:
            #
            #     r2*l*a <= D
            #
            # therefore:
            #
            #     a <= D/(r2*l)
            # ------------------------------------------------

            a_max = min(
                r1 - 1,
                D // (r2 * l)
            )

            for a in range(a_max + 1):

                tests += 1

                numerator = (
                    D
                    - r2 * l * a
                )

                denominator = (
                    r1 * k + a
                )

                if numerator < 0:
                    continue

                if denominator <= 0:
                    continue

                # ------------------------------------------------
                # Exact integrality condition.
                # ------------------------------------------------

                if numerator % denominator != 0:
                    continue

                b = numerator // denominator

                if not (0 <= b < r2):
                    continue

                # ------------------------------------------------
                # Verify E exactly.
                # ------------------------------------------------

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

                p = r1 * k + a
                q = r2 * l + b

                out[(p, q)] = {
                    "p": p,
                    "q": q,
                    "K": K,
                    "E": E,
                    "k": k,
                    "l": l,
                    "a": a,
                    "b": b,
                }

    elapsed = time.perf_counter() - start

    return list(out.values()), {
        "time": elapsed,
        "tests": tests,
    }


# ------------------------------------------------------------
# Normalize candidate sets
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
# Print examples
# ------------------------------------------------------------

def print_examples(candidates):

    if not candidates:

        print("    none")
        return

    for c in candidates[:PRINT_EXAMPLES]:

        print(
            "    "
            f"p={c['p']} "
            f"q={c['q']} "
            f"K={c['K']} "
            f"E={c['E']} "
            f"k={c['k']} "
            f"l={c['l']} "
            f"a={c['a']} "
            f"b={c['b']}"
        )


# ------------------------------------------------------------
# Run one comparison
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
    # Bounded
    # --------------------------------------------------------

    bounded, sbd = bounded_candidates(
        n,
        r1,
        r2,
    )

    # --------------------------------------------------------
    # Direct
    # --------------------------------------------------------

    direct, sd = direct_candidates(
        n,
        r1,
        r2,
    )

    kb = candidate_keys(baseline)
    kbd = candidate_keys(bounded)
    kd = candidate_keys(direct)

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

    true_key_1 = (
        p_true,
        q_true,
    )

    true_key_2 = (
        q_true,
        p_true,
    )

    print(
        f"    true in baseline    = "
        f"{true_key_1 in kb or true_key_2 in kb}"
    )

    print(
        f"    true in bounded     = "
        f"{true_key_1 in kbd or true_key_2 in kbd}"
    )

    print(
        f"    true in direct      = "
        f"{true_key_1 in kd or true_key_2 in kd}"
    )

    # --------------------------------------------------------
    # Speed
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

    # --------------------------------------------------------
    # Loop counts
    # --------------------------------------------------------

    print()
    print("SEARCH WORK")

    print(
        f"    baseline a,b tests = "
        f"{sb['tests']:,}"
    )

    print(
        f"    bounded  a,b tests = "
        f"{sbd['tests']:,}"
    )

    print(
        f"    direct   a tests   = "
        f"{sd['tests']:,}"
    )

    print()
    print("SPEEDUP")

    print(
        f"    bounded vs baseline = "
        f"{sb['time'] / max(sbd['time'], 1e-12):.3f}x"
    )

    print(
        f"    direct vs baseline  = "
        f"{sb['time'] / max(sd['time'], 1e-12):.3f}x"
    )

    # --------------------------------------------------------
    # Candidate examples
    # --------------------------------------------------------

    if direct:

        print()
        print("DIRECT EXAMPLES")

        print_examples(direct)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("K/E RESIDUAL-BOUNDED SEARCH")
    print("EXPERIMENT 90")
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
            f"n       = {n}"
        )

        print(
            f"true p  = {p_true}"
        )

        print(
            f"true q  = {q_true}"
        )

        print(
            f"|p-q|   = "
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
    print("EXPERIMENT 90 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 90
# ============================================================
