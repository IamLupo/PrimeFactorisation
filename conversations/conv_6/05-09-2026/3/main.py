import math
import random
import time
import sympy


# ============================================================
# START EXPERIMENT 86
# ============================================================
#
# MULTI-R1 / MULTI-R2 GRID PROPAGATION
#
# Hypothesis:
#
#   A candidate discovered at a small (r1,r2) pair can be
#   represented at many larger (r1,r2) pairs.
#
#   The different representations should provide multiple
#   constraints, potentially reducing false candidates.
#
#
# IMPORTANT:
#
#   Experiment 85 used only:
#
#       (r1,r2) -> (2r1,2r2) -> (4r1,4r2)
#
#   Experiment 86 instead uses a GRID:
#
#       r1 candidates x r2 candidates
#
#   at every level.
#
#
# Two modes:
#
#   synthetic:
#       Fast. Builds a candidate cloud containing the true
#       factorization + false factor pairs.
#
#   search:
#       Performs an actual small-r search.
#
# ============================================================


# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

MODE = "synthetic"

TARGET = 10**9

LEVELS = 6

# Number of prime choices on each side of each level.
#
# Example:
#
#   GRID_SIZE = 5
#
# gives 5 r1 values x 5 r2 values = 25 representations.
#
GRID_SIZE = 5

# Relative spread around the nominal doubled scale.
#
# 0.20 means +/-20%.
#
LEVEL_SPREAD = 0.20

# Number of synthetic false candidates.
SYNTHETIC_FALSE = 5000

# Initial primes for the search mode.
INITIAL_R1 = 11
INITIAL_R2 = 13

# Search mode E limit.
SEARCH_E_LIMIT = 2000

# Random seed
SEED = 86


# ------------------------------------------------------------
# Prime utilities
# ------------------------------------------------------------

def nearest_prime(x):
    """
    Return a nearby prime.
    """
    x = max(2, int(x))
    return int(sympy.nextprime(x - 1))


def prime_grid(center, count, spread):
    """
    Produce approximately `count` primes around

        center * [1-spread, 1+spread]

    """

    if count <= 1:
        return [nearest_prime(center)]

    lo = center * (1.0 - spread)
    hi = center * (1.0 + spread)

    if lo < 3:
        lo = 3

    values = []

    for i in range(count):

        t = i / (count - 1)

        x = lo + (hi - lo) * t

        p = nearest_prime(x)

        values.append(p)

    # Remove accidental duplicates.
    return sorted(set(values))


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

    q0 = target // p

    q = int(
        sympy.randprime(
            max(3, int(q0 * 0.85)),
            int(q0 * 1.15),
        )
    )

    return p, q, p * q


# ------------------------------------------------------------
# State calculation
# ------------------------------------------------------------

def state_from_factors(p, q, r1, r2):

    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

    K = k * l

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    R = r1 * r2

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    E = c1 + c2 + c3

    Q = None

    return {
        "p": p,
        "q": q,
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


# ------------------------------------------------------------
# Level consistency test
# ------------------------------------------------------------

def level_test(candidate, n, r1, r2):
    """
    Test

        K + E == floor(n / (r1*r2))

    for the candidate's p,q representation.
    """

    p = candidate["p"]
    q = candidate["q"]

    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

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

    Q = n // R

    return K + E == Q


# ------------------------------------------------------------
# Synthetic candidate generation
# ------------------------------------------------------------

def make_synthetic_cloud(
    p,
    q,
    n,
    count,
    seed,
):
    """
    Build:

        1 true candidate
        count false candidates

    False candidates are deliberately chosen with similar
    magnitude to the real factors.
    """

    rng = random.Random(seed)

    candidates = []

    # TRUE candidate.
    candidates.append({
        "p": p,
        "q": q,
        "true": True,
    })

    root = math.isqrt(n)

    while len(candidates) < count + 1:

        fp = rng.randint(
            max(2, int(root * 0.20)),
            int(root * 2.0),
        )

        fq = rng.randint(
            max(2, int(root * 0.20)),
            int(root * 2.0),
        )

        if fp == p and fq == q:
            continue

        candidates.append({
            "p": fp,
            "q": fq,
            "true": False,
        })

    return candidates


# ------------------------------------------------------------
# Search mode
# ------------------------------------------------------------

def search_small_r(n, r1, r2, e_limit):
    """
    Actual candidate search.

    This is intentionally bounded because the purpose of
    Experiment 86 is to study propagation, not optimize
    the original K/E search yet.
    """

    Q = n // (r1 * r2)

    candidates = {}

    for E in range(e_limit + 1):

        K = Q - E

        if K <= 0:
            continue

        divisors = sympy.divisors(K)

        for k in divisors:

            l = K // k

            for a in range(r1):

                # b is tested exhaustively at small r2.
                for b in range(r2):

                    K2 = k * l

                    c1 = (k * b) // r2
                    c2 = (l * a) // r1

                    beta = (k * b) % r2
                    alpha = (l * a) % r1

                    c3 = (
                        r1 * beta
                        + r2 * alpha
                        + a * b
                    ) // (r1 * r2)

                    E2 = c1 + c2 + c3

                    if E2 != E:
                        continue

                    p = r1 * k + a
                    q = r2 * l + b

                    key = (p, q)

                    candidates[key] = {
                        "p": p,
                        "q": q,
                        "true": False,
                    }

    return list(candidates.values())


# ------------------------------------------------------------
# Build level grids
# ------------------------------------------------------------

def build_level_grids():

    grids = []

    # Level 0.
    center_r1 = INITIAL_R1
    center_r2 = INITIAL_R2

    for level in range(LEVELS + 1):

        if level == 0:

            c1 = INITIAL_R1
            c2 = INITIAL_R2

        else:

            c1 = INITIAL_R1 * (2 ** level)
            c2 = INITIAL_R2 * (2 ** level)

        grid_r1 = prime_grid(
            c1,
            GRID_SIZE,
            LEVEL_SPREAD,
        )

        grid_r2 = prime_grid(
            c2,
            GRID_SIZE,
            LEVEL_SPREAD,
        )

        grids.append(
            (grid_r1, grid_r2)
        )

    return grids


# ------------------------------------------------------------
# Evaluate candidate against an entire grid
# ------------------------------------------------------------

def evaluate_grid(
    candidates,
    n,
    grid_r1,
    grid_r2,
):
    """
    A candidate survives a grid if there exists at least
    one (r1,r2) representation satisfying the level equation.
    """

    survivors = []

    for candidate in candidates:

        survive = False
        hits = 0

        for r1 in grid_r1:

            for r2 in grid_r2:

                if level_test(
                    candidate,
                    n,
                    r1,
                    r2,
                ):
                    survive = True
                    hits += 1

        if survive:

            item = dict(candidate)

            item["hits"] = hits

            survivors.append(item)

    return survivors


# ------------------------------------------------------------
# Strict grid intersection
# ------------------------------------------------------------

def evaluate_grid_strict(
    candidates,
    n,
    grid_r1,
    grid_r2,
):
    """
    Stronger version.

    Candidate must satisfy the level constraint for EVERY
    (r1,r2) pair in the grid.

    This is useful for measuring the strength of intersecting
    multiple representations.
    """

    survivors = []

    total_tests = (
        len(grid_r1) *
        len(grid_r2)
    )

    for candidate in candidates:

        hits = 0

        for r1 in grid_r1:

            for r2 in grid_r2:

                if level_test(
                    candidate,
                    n,
                    r1,
                    r2,
                ):
                    hits += 1

        if hits == total_tests:

            item = dict(candidate)

            item["hits"] = hits

            survivors.append(item)

    return survivors


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("MULTI-R1 / MULTI-R2 GRID PROPAGATION")
    print("EXPERIMENT 86")
    print("=" * 78)

    # --------------------------------------------------------
    # Generate test number
    # --------------------------------------------------------

    p, q, n = generate_semiprime(TARGET)

    print()
    print(f"MODE     = {MODE}")
    print(f"n        = {n}")
    print(f"p        = {p}")
    print(f"q        = {q}")
    print(f"|p-q|    = {abs(p-q)}")

    # --------------------------------------------------------
    # Candidate cloud
    # --------------------------------------------------------

    if MODE == "synthetic":

        candidates = make_synthetic_cloud(
            p,
            q,
            n,
            SYNTHETIC_FALSE,
            SEED,
        )

    elif MODE == "search":

        start = time.perf_counter()

        candidates = search_small_r(
            n,
            INITIAL_R1,
            INITIAL_R2,
            SEARCH_E_LIMIT,
        )

        elapsed = (
            time.perf_counter()
            - start
        )

        # Mark actual true candidate.
        for c in candidates:

            if (
                c["p"] == p
                and c["q"] == q
            ):
                c["true"] = True

        print()
        print(
            "Initial search time = "
            f"{elapsed:.6f} s"
        )

    else:

        raise ValueError(
            "MODE must be 'synthetic' or 'search'"
        )

    print()
    print(
        "Initial candidates = "
        f"{len(candidates):,}"
    )

    true_count = sum(
        1 for c in candidates
        if c.get("true", False)
    )

    print(
        "True candidates    = "
        f"{true_count:,}"
    )

    # --------------------------------------------------------
    # Build grids
    # --------------------------------------------------------

    grids = build_level_grids()

    # --------------------------------------------------------
    # Run levels
    # --------------------------------------------------------

    for level, (grid_r1, grid_r2) in enumerate(grids):

        print()
        print("-" * 78)
        print(f"LEVEL {level}")
        print("-" * 78)

        print(
            "r1 grid:",
            grid_r1,
        )

        print(
            "r2 grid:",
            grid_r2,
        )

        grid_count = (
            len(grid_r1)
            * len(grid_r2)
        )

        print(
            f"grid combinations = "
            f"{grid_count}"
        )

        # ----------------------------------------------------
        # OR test
        #
        # Candidate survives if at least ONE representation
        # survives.
        # ----------------------------------------------------

        start = time.perf_counter()

        survivors_or = evaluate_grid(
            candidates,
            n,
            grid_r1,
            grid_r2,
        )

        t_or = (
            time.perf_counter()
            - start
        )

        true_or = sum(
            1
            for c in survivors_or
            if c.get("true", False)
        )

        print()
        print("OR GRID")
        print(
            f"before              = "
            f"{len(candidates):,}"
        )
        print(
            f"after               = "
            f"{len(survivors_or):,}"
        )
        print(
            f"true alive          = "
            f"{true_or:,}"
        )
        print(
            f"time                = "
            f"{t_or:.6f} s"
        )

        if survivors_or:

            avg_hits = sum(
                c["hits"]
                for c in survivors_or
            ) / len(survivors_or)

            print(
                f"average grid hits   = "
                f"{avg_hits:.3f} / {grid_count}"
            )

        # ----------------------------------------------------
        # AND test
        #
        # Candidate must survive EVERY representation.
        # ----------------------------------------------------

        start = time.perf_counter()

        survivors_and = evaluate_grid_strict(
            candidates,
            n,
            grid_r1,
            grid_r2,
        )

        t_and = (
            time.perf_counter()
            - start
        )

        true_and = sum(
            1
            for c in survivors_and
            if c.get("true", False)
        )

        print()
        print("AND GRID")
        print(
            f"before              = "
            f"{len(candidates):,}"
        )
        print(
            f"after               = "
            f"{len(survivors_and):,}"
        )
        print(
            f"true alive          = "
            f"{true_and:,}"
        )
        print(
            f"time                = "
            f"{t_and:.6f} s"
        )

        # ----------------------------------------------------
        # Continue with the stronger AND intersection.
        # ----------------------------------------------------

        candidates = survivors_and

        if not candidates:

            print()
            print(
                "ALL CANDIDATES WERE ELIMINATED."
            )
            break

        # ----------------------------------------------------
        # Check actual factorization candidates.
        # ----------------------------------------------------

        exact = []

        for c in candidates:

            if c["p"] * c["q"] == n:

                exact.append(c)

        print()
        print(
            "Exact p*q == n survivors = "
            f"{len(exact):,}"
        )

        if exact:

            print()
            print("EXACT FACTORIZATION FOUND")

            for c in exact[:10]:

                print(
                    f"p={c['p']} "
                    f"q={c['q']}"
                )

    # --------------------------------------------------------
    # Final
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL RESULT")
    print("=" * 78)

    if candidates:

        exact = [
            c
            for c in candidates
            if c["p"] * c["q"] == n
        ]

        print(
            f"remaining candidates = "
            f"{len(candidates):,}"
        )

        print(
            f"exact factorizations  = "
            f"{len(exact):,}"
        )

        if exact:

            for c in exact[:10]:

                print(
                    f"p={c['p']}"
                )

                print(
                    f"q={c['q']}"
                )

    else:

        print("No candidates remain.")


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 86
# ============================================================
