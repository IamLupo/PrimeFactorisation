import math
import time
import sympy


# ============================================================
# START EXPERIMENT 88
# ============================================================
#
# SMALL-r CANDIDATE CLOUD
#          +
# MODULAR RESIDUAL FILTERING
#
# We first generate candidates from:
#
#     Q = floor(n / (r1*r2))
#
#     K + E = Q
#
# using small r1,r2.
#
# We do NOT require:
#
#     p*q == n
#
# during candidate generation.
#
# Then each candidate has:
#
#     D = p*q - n
#
# A real factorization has:
#
#     D = 0
#
# therefore:
#
#     D == 0 mod s
#
# for EVERY prime s.
#
# We test increasing sets of small primes.
#
#
# IMPORTANT:
#
# This experiment is specifically testing whether the
# candidate cloud produced by the K/E construction can be
# collapsed using cheap independent modular constraints.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

INITIAL_R1 = 11
INITIAL_R2 = 13

# Initial E search range.
E_LIMIT_MULTIPLIER = 3.0

# Modular primes.
#
# We test cumulative prefixes:
#
#   first 1 prime
#   first 2 primes
#   first 4 primes
#   first 8 primes
#   ...
#
MODULI = [
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

# Print a few survivors.
PRINT_CANDIDATES = 20


# ------------------------------------------------------------
# Generate balanced semiprime
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
# Initial candidate cloud
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

        for k in sympy.divisors(K):

            l = K // k

            divisor_pairs += 1

            for a in range(r1):

                for b in range(r2):

                    ab_tests += 1

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

                    key = (p, q)

                    if key in candidates:
                        continue

                    candidates[key] = {
                        "p": p,
                        "q": q,
                        "product": p * q,
                        "residual": p * q - n,

                        # Save original state too.
                        "k0": k,
                        "l0": l,
                        "a0": a,
                        "b0": b,
                        "K0": K,
                        "E0": E,
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
# Modular filtering
# ------------------------------------------------------------

def filter_modulus(candidates, modulus):

    survivors = []

    for candidate in candidates:

        D = candidate["residual"]

        if D % modulus == 0:

            survivors.append(candidate)

    return survivors


# ------------------------------------------------------------
# Candidate statistics
# ------------------------------------------------------------

def residual_statistics(candidates, n):

    if not candidates:

        return None

    residuals = [
        abs(c["residual"])
        for c in candidates
    ]

    return {
        "min": min(residuals),
        "max": max(residuals),
        "sum": sum(residuals),
    }


# ------------------------------------------------------------
# Print examples
# ------------------------------------------------------------

def print_examples(candidates, n):

    if not candidates:

        print("    none")
        return

    for c in candidates[:PRINT_CANDIDATES]:

        D = c["residual"]

        print(
            "    "
            f"p={c['p']} "
            f"q={c['q']} "
            f"D=pq-n={D} "
            f"|D|={abs(D)}"
        )


# ------------------------------------------------------------
# Check true candidate
# ------------------------------------------------------------

def true_alive(candidates, p, q):

    for c in candidates:

        if (
            (c["p"] == p and c["q"] == q)
            or
            (c["p"] == q and c["q"] == p)
        ):

            return True

    return False


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("MODULAR RESIDUAL FILTERING")
    print("EXPERIMENT 88")
    print("=" * 78)

    # --------------------------------------------------------
    # Generate n
    # --------------------------------------------------------

    p_true, q_true, n = generate_semiprime(
        TARGET
    )

    print()
    print(f"n        = {n}")
    print(f"p        = {p_true}")
    print(f"q        = {q_true}")
    print(
        f"|p-q|    = "
        f"{abs(p_true-q_true)}"
    )

    # --------------------------------------------------------
    # Initial candidate generation
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("INITIAL SMALL-r CANDIDATE CLOUD")
    print("=" * 78)

    candidates, stats = generate_candidates(
        n,
        INITIAL_R1,
        INITIAL_R2,
    )

    print()
    print(
        f"r1                 = "
        f"{INITIAL_R1}"
    )

    print(
        f"r2                 = "
        f"{INITIAL_R2}"
    )

    print(
        f"Q                  = "
        f"{stats['Q']:,}"
    )

    print(
        f"E limit            = "
        f"{stats['E_limit']:,}"
    )

    print(
        f"divisor pairs      = "
        f"{stats['divisor_pairs']:,}"
    )

    print(
        f"a,b tests          = "
        f"{stats['ab_tests']:,}"
    )

    print(
        f"candidates         = "
        f"{len(candidates):,}"
    )

    print(
        f"true alive         = "
        f"{true_alive(candidates,p_true,q_true)}"
    )

    print(
        f"time               = "
        f"{stats['time']:.6f} s"
    )

    # --------------------------------------------------------
    # Residual distribution
    # --------------------------------------------------------

    s = residual_statistics(
        candidates,
        n,
    )

    if s:

        print()
        print("Residual statistics:")
        print(
            f"    min |D| = {s['min']}"
        )
        print(
            f"    max |D| = {s['max']}"
        )

    # --------------------------------------------------------
    # Print initial candidates
    # --------------------------------------------------------

    print()
    print("Initial examples:")

    print_examples(
        candidates,
        n,
    )

    # --------------------------------------------------------
    # Sequential modular filtering
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("MODULAR FILTERING")
    print("=" * 78)

    previous_count = len(candidates)

    for i, modulus in enumerate(MODULI, start=1):

        start = time.perf_counter()

        before = len(candidates)

        candidates = filter_modulus(
            candidates,
            modulus,
        )

        elapsed = time.perf_counter() - start

        after = len(candidates)

        alive = true_alive(
            candidates,
            p_true,
            q_true,
        )

        print()
        print(
            f"MODULUS {modulus}"
        )

        print(
            f"    before       = "
            f"{before:,}"
        )

        print(
            f"    after        = "
            f"{after:,}"
        )

        print(
            f"    reduction    = "
            f"{before / max(1,after):.3f}x"
        )

        print(
            f"    true alive   = "
            f"{alive}"
        )

        print(
            f"    time         = "
            f"{elapsed:.9f} s"
        )

        if candidates:

            residuals = [
                abs(c["residual"])
                for c in candidates
            ]

            print(
                f"    min |D|      = "
                f"{min(residuals)}"
            )

            print(
                f"    max |D|      = "
                f"{max(residuals)}"
            )

        if not candidates:

            print()
            print(
                "ALL CANDIDATES ELIMINATED."
            )

            break

        # ----------------------------------------------------
        # Print survivors when the cloud becomes small.
        # ----------------------------------------------------

        if (
            after <= PRINT_CANDIDATES
            or after < before / 5
        ):

            print()
            print("    survivors:")

            print_examples(
                candidates,
                n,
            )

        # ----------------------------------------------------
        # Detect accidental death of true candidate.
        # ----------------------------------------------------

        if not alive:

            print()
            print(
                "WARNING:"
            )

            print(
                "TRUE FACTORIZATION WAS "
                "ELIMINATED."
            )

            break

    # --------------------------------------------------------
    # Final exact check
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL RESULT")
    print("=" * 78)

    exact = []

    for c in candidates:

        if c["residual"] == 0:

            exact.append(c)

    print(
        f"remaining candidates = "
        f"{len(candidates):,}"
    )

    print(
        f"exact residual zero   = "
        f"{len(exact):,}"
    )

    if exact:

        for c in exact:

            print(
                f"p={c['p']} "
                f"q={c['q']}"
            )

    else:

        print()
        print(
            "No exact factorization "
            "remained."
        )

    print()
    print("=" * 78)
    print("EXPERIMENT COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 88
# ============================================================
