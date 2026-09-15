import math
import time
import sympy


# ============================================================
# START EXPERIMENT 87
# ============================================================
#
# REAL SMALL-r CANDIDATE CLOUD
#        ->
# MULTIPLE r1/r2 GRID FILTERS
#        ->
# CANDIDATE COLLAPSE TEST
#
#
# We first create the SAME KIND of candidate cloud as
# Experiment 85:
#
#     Q = floor(n / (r1*r2))
#
#     K + E = Q
#
# with small r1,r2.
#
# Importantly, we do NOT require p*q == n.
#
# Each surviving candidate implies:
#
#     p = r1*k + a
#     q = r2*l + b
#
# and therefore a candidate product:
#
#     m = p*q
#
# At later levels we test:
#
#     floor(m / (r1'*r2'))
#       ==
#     floor(n / (r1'*r2'))
#
# for MANY combinations of r1' and r2'.
#
# This is equivalent to the K+E condition at the new
# representation, but much cheaper to evaluate.
#
#
# Main hypothesis:
#
#     The TRUE candidate survives many independent
#     r1/r2 representations, while false candidates
#     gradually disappear.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

INITIAL_R1 = 11
INITIAL_R2 = 13

# Number of levels after the initial candidate generation.
LEVELS = 7

# Number of primes in each dimension of the grid.
#
# GRID_SIZE=5 -> at most 25 (r1,r2) combinations per level.
#
GRID_SIZE = 5

# Width of each grid around the nominal 2x scale.
#
# Example:
#
#   center = 100
#   spread = 0.20
#
# gives approximately 80..120.
#
GRID_SPREAD = 0.20

# Initial E search.
#
# For balanced factors this is normally comfortably above
# the actual E.
#
E_LIMIT_MULTIPLIER = 3.0

# Whether to also test the candidate against all previous
# grid levels cumulatively.
#
# TRUE means a candidate must have survived every previous
# level as well as the current level.
CUMULATIVE = True

# Print at most this many candidates.
PRINT_CANDIDATES = 15


# ------------------------------------------------------------
# Prime helpers
# ------------------------------------------------------------

def nearest_prime(x):
    x = max(2, int(x))
    return int(sympy.nextprime(x - 1))


def make_prime_grid(center, count, spread):
    """
    Return several primes around `center`.
    """

    if count <= 1:
        return [nearest_prime(center)]

    lo = max(3, int(center * (1.0 - spread)))
    hi = max(lo + 2, int(center * (1.0 + spread)))

    values = []

    for i in range(count):

        x = lo + (hi - lo) * i / (count - 1)

        values.append(
            nearest_prime(x)
        )

    return sorted(set(values))


# ------------------------------------------------------------
# Generate a balanced semiprime
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
# Candidate carry state
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
# Initial candidate generation
# ------------------------------------------------------------

def generate_initial_candidates(
    n,
    r1,
    r2,
):
    """
    Generate the actual candidate cloud.

    We use:

        Q = floor(n / R)
        K = Q - E

    and enumerate divisor pairs k*l=K.

    For each divisor pair we enumerate the small a,b space.

    This intentionally does NOT require p*q == n.
    """

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

    start = time.perf_counter()

    divisor_pairs = 0
    ab_tests = 0

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
# Grid test
# ------------------------------------------------------------

def grid_test(candidate, n, grid_r1, grid_r2):
    """
    Test one candidate against every (r1,r2) pair.

    The candidate survives a grid cell when:

        floor(p*q / R) == floor(n / R)

    for R = r1*r2.
    """

    m = candidate["product"]

    hits = 0
    total = 0

    for r1 in grid_r1:

        for r2 in grid_r2:

            total += 1

            R = r1 * r2

            if m // R == n // R:
                hits += 1

    return hits, total


# ------------------------------------------------------------
# Exact candidate test
# ------------------------------------------------------------

def is_exact(candidate, n):
    return candidate["product"] == n


# ------------------------------------------------------------
# Print candidates
# ------------------------------------------------------------

def print_candidate_examples(
    candidates,
    n,
    label,
):
    print()
    print(label)

    if not candidates:
        print("    none")
        return

    for c in candidates[:PRINT_CANDIDATES]:

        exact = is_exact(c, n)

        delta = c["product"] - n

        print(
            "    "
            f"p={c['p']} "
            f"q={c['q']} "
            f"|pq-n|={abs(delta)} "
            f"exact={exact}"
        )


# ------------------------------------------------------------
# Build all level grids
# ------------------------------------------------------------

def build_grids():

    grids = []

    for level in range(LEVELS + 1):

        scale = 2 ** level

        center_r1 = INITIAL_R1 * scale
        center_r2 = INITIAL_R2 * scale

        grid_r1 = make_prime_grid(
            center_r1,
            GRID_SIZE,
            GRID_SPREAD,
        )

        grid_r2 = make_prime_grid(
            center_r2,
            GRID_SIZE,
            GRID_SPREAD,
        )

        grids.append(
            (
                grid_r1,
                grid_r2,
            )
        )

    return grids


# ------------------------------------------------------------
# Run experiment
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("REAL SMALL-r CANDIDATE CLOUD")
    print("MULTI-R1 / MULTI-R2 GRID FILTER")
    print("EXPERIMENT 87")
    print("=" * 78)

    # --------------------------------------------------------
    # Generate test semiprime
    # --------------------------------------------------------

    p_true, q_true, n = generate_semiprime(
        TARGET
    )

    print()
    print(f"n        = {n}")
    print(f"p        = {p_true}")
    print(f"q        = {q_true}")
    print(f"|p-q|    = {abs(p_true - q_true)}")

    # --------------------------------------------------------
    # LEVEL 0: candidate cloud
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("LEVEL 0 — INITIAL SMALL-r SEARCH")
    print("=" * 78)

    candidates, stats = generate_initial_candidates(
        n,
        INITIAL_R1,
        INITIAL_R2,
    )

    true_alive = any(
        c["p"] == p_true
        and c["q"] == q_true
        for c in candidates
    )

    print()
    print(f"r1                 = {INITIAL_R1}")
    print(f"r2                 = {INITIAL_R2}")
    print(f"Q                  = {stats['Q']:,}")
    print(f"E limit            = {stats['E_limit']:,}")
    print(
        "divisor pairs      = "
        f"{stats['divisor_pairs']:,}"
    )
    print(
        "a,b tests          = "
        f"{stats['ab_tests']:,}"
    )
    print(
        "candidates         = "
        f"{len(candidates):,}"
    )
    print(
        "TRUE candidate     = "
        f"{true_alive}"
    )
    print(
        "search time        = "
        f"{stats['time']:.6f} s"
    )

    print_candidate_examples(
        candidates,
        n,
        "Initial candidate examples:",
    )

    # --------------------------------------------------------
    # Build grids
    # --------------------------------------------------------

    grids = build_grids()

    # --------------------------------------------------------
    # Cumulative candidate filtering
    # --------------------------------------------------------

    for level in range(1, LEVELS + 1):

        grid_r1, grid_r2 = grids[level]

        print()
        print("=" * 78)
        print(f"LEVEL {level}")
        print("=" * 78)

        print(
            f"r1 grid = {grid_r1}"
        )

        print(
            f"r2 grid = {grid_r2}"
        )

        grid_size = (
            len(grid_r1)
            * len(grid_r2)
        )

        print(
            f"grid cells = {grid_size}"
        )

        before = len(candidates)

        start = time.perf_counter()

        survivors = []

        for candidate in candidates:

            hits, total = grid_test(
                candidate,
                n,
                grid_r1,
                grid_r2,
            )

            # Strict intersection.
            #
            # Candidate must survive EVERY cell.
            if hits == total:

                survivor = dict(candidate)

                survivor[
                    f"hits_L{level}"
                ] = hits

                survivor[
                    f"grid_L{level}"
                ] = total

                survivors.append(
                    survivor
                )

        elapsed = time.perf_counter() - start

        candidates = survivors

        true_alive = any(
            c["p"] == p_true
            and c["q"] == q_true
            for c in candidates
        )

        exact_count = sum(
            1
            for c in candidates
            if is_exact(c, n)
        )

        print()
        print(
            f"before              = "
            f"{before:,}"
        )

        print(
            f"after               = "
            f"{len(candidates):,}"
        )

        print(
            f"reduction           = "
            f"{before / max(1, len(candidates)):.3f}x"
        )

        print(
            f"TRUE alive          = "
            f"{true_alive}"
        )

        print(
            f"exact p*q == n      = "
            f"{exact_count}"
        )

        print(
            f"time                = "
            f"{elapsed:.6f} s"
        )

        print_candidate_examples(
            candidates,
            n,
            "Survivor examples:",
        )

        if not candidates:

            print()
            print(
                "ALL CANDIDATES ELIMINATED."
            )

            break

        if not true_alive:

            print()
            print(
                "WARNING:"
            )

            print(
                "The TRUE candidate was "
                "eliminated."
            )

            break

    # --------------------------------------------------------
    # Final
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL RESULT")
    print("=" * 78)

    exact = [
        c
        for c in candidates
        if is_exact(c, n)
    ]

    print(
        f"remaining candidates = "
        f"{len(candidates):,}"
    )

    print(
        f"exact factorizations = "
        f"{len(exact):,}"
    )

    if exact:

        for c in exact:

            print(
                f"p={c['p']} "
                f"q={c['q']}"
            )

    elif candidates:

        print()
        print(
            "No exact factorization yet, "
            "but candidates remain."
        )

    print()
    print("=" * 78)
    print("EXPERIMENT COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 87
# ============================================================
