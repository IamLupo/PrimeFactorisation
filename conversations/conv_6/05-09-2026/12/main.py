import math
import time
import itertools
import sympy


# ============================================================
# START EXPERIMENT 95
# ============================================================
#
# 2D GRID SCALING TEST
#
# Test whether increasing the number of r1/r2 values gives
# stronger K/E + carry collapse.
#
#
# Grid sizes tested:
#
#     2x2
#     3x3
#     4x4
#
#
# For every grid:
#
#     Q_ij = floor(n / (r1_i*r2_j))
#
#     Q_ij = k_i*l_j + E_ij
#
# with:
#
#     k_i = floor(p/r1_i)
#     l_j = floor(q/r2_j)
#
#
# Constraints:
#
#   1. K_ij = k_i*l_j
#
#   2. All rows must correspond to the SAME p:
#
#        r1_i*k_i <= p < r1_i*(k_i+1)
#
#   3. All columns must correspond to the SAME q.
#
#   4. The exact carry equation must reproduce E_ij.
#
#
# We DO NOT use p*q == n as a search filter.
#
# p*q == n is only used for final classification.
#
#
# IMPORTANT:
#
# The goal is not yet to solve huge n.
# The goal is to measure whether more 2D r1/r2 information
# systematically reduces the number of surviving structures.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

CASE_COUNT = 5

GRID_SIZES = [
    2,
    3,
    4,
]

# First r1/r2 primes.
#
# The grids are generated around these scales.
R1_CENTER = 7
R2_CENTER = 11

GRID_SPREAD = 0.60

# Search around sqrt(n)/r at the anchor cell.
#
# This must be wide enough to contain the true quotient.
ANCHOR_WINDOW = 1800

# We don't print every survivor.
PRINT_SURVIVORS = 10

# Safety limit for candidate structures.
#
# If exceeded, we stop that grid and report that the search
# became too large.
MAX_STRUCTURES = 2_000_000

# Random seed.
SEED = 9501


# ------------------------------------------------------------
# Prime helpers
# ------------------------------------------------------------

def nearest_prime(x):

    return int(
        sympy.nextprime(
            max(2, int(x)) - 1
        )
    )


def make_prime_grid(
    center,
    count,
    spread,
):

    if count == 1:
        return [
            nearest_prime(center)
        ]

    lo = max(
        3,
        int(center * (1.0 - spread))
    )

    hi = max(
        lo + 2,
        int(center * (1.0 + spread))
    )

    values = []

    for i in range(count):

        x = (
            lo
            + (hi - lo)
            * i
            / (count - 1)
        )

        values.append(
            nearest_prime(x)
        )

    return sorted(
        set(values)
    )


# ------------------------------------------------------------
# Semiprime generator
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
# Common interval
#
# Given:
#
#     x = r*k + a
#
# with:
#
#     0 <= a < r
#
# then:
#
#     r*k <= x <= r*(k+1)-1
# ------------------------------------------------------------

def common_interval(
    rs,
    ks,
):

    lo = 0
    hi = None

    for r, k in zip(rs, ks):

        this_lo = r * k
        this_hi = r * (k + 1) - 1

        lo = max(
            lo,
            this_lo,
        )

        if hi is None:
            hi = this_hi
        else:
            hi = min(
                hi,
                this_hi,
            )

    if hi is None:
        return None

    if lo > hi:
        return None

    return lo, hi


# ------------------------------------------------------------
# Exact E from p,q,r1,r2
# ------------------------------------------------------------

def exact_E(
    p,
    q,
    r1,
    r2,
):

    k, a = divmod(
        p,
        r1,
    )

    l, b = divmod(
        q,
        r2,
    )

    c1 = (
        k * b
    ) // r2

    c2 = (
        l * a
    ) // r1

    beta = (
        k * b
    ) % r2

    alpha = (
        l * a
    ) % r1

    R = r1 * r2

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    return c1 + c2 + c3


# ------------------------------------------------------------
# Validate complete candidate
# ------------------------------------------------------------

def validate_candidate(
    p,
    q,
    n,
    r1_grid,
    r2_grid,
    ks,
    ls,
    E_matrix,
):

    rows = len(r1_grid)
    cols = len(r2_grid)

    # --------------------------------------------------------
    # Quotient consistency.
    # --------------------------------------------------------

    for i, r1 in enumerate(
        r1_grid
    ):

        if p // r1 != ks[i]:
            return False

    for j, r2 in enumerate(
        r2_grid
    ):

        if q // r2 != ls[j]:
            return False

    # --------------------------------------------------------
    # Carry consistency.
    # --------------------------------------------------------

    for i, r1 in enumerate(
        r1_grid
    ):

        for j, r2 in enumerate(
            r2_grid
        ):

            E = exact_E(
                p,
                q,
                r1,
                r2,
            )

            if E != E_matrix[i][j]:
                return False

    return True


# ------------------------------------------------------------
# Build Q matrix.
# ------------------------------------------------------------

def build_Q_matrix(
    n,
    r1_grid,
    r2_grid,
):

    return [
        [
            n // (r1 * r2)
            for r2 in r2_grid
        ]
        for r1 in r1_grid
    ]


# ------------------------------------------------------------
# Search one grid
# ------------------------------------------------------------

def search_grid(
    n,
    p_true,
    q_true,
    r1_grid,
    r2_grid,
):

    rows = len(r1_grid)
    cols = len(r2_grid)

    Q = build_Q_matrix(
        n,
        r1_grid,
        r2_grid,
    )

    # --------------------------------------------------------
    # Anchor.
    # --------------------------------------------------------

    r1_0 = r1_grid[0]
    r2_0 = r2_grid[0]

    k_center = math.isqrt(n) // r1_0
    l_center = math.isqrt(n) // r2_0

    k_lo = max(
        1,
        k_center - ANCHOR_WINDOW,
    )

    k_hi = (
        k_center
        + ANCHOR_WINDOW
    )

    l_lo = max(
        1,
        l_center - ANCHOR_WINDOW,
    )

    l_hi = (
        l_center
        + ANCHOR_WINDOW
    )

    # --------------------------------------------------------
    # Anchor candidates.
    #
    # We don't independently enumerate E.
    #
    # E = Q00-k0*l0
    #
    # and require E to be compatible with the carry scale.
    # --------------------------------------------------------

    anchors = []

    Q00 = Q[0][0]

    for k0 in range(
        k_lo,
        k_hi + 1,
    ):

        l_max = min(
            l_hi,
            Q00 // k0,
        )

        for l0 in range(
            l_lo,
            l_max + 1,
        ):

            E0 = Q00 - k0 * l0

            if E0 < 0:
                continue

            # Conservative carry bound.
            #
            # c1 <= k-1
            # c2 <= l-1
            # c3 is small.
            if E0 > k0 + l0 + 10:
                continue

            anchors.append(
                (k0, l0)
            )

    # --------------------------------------------------------
    # Statistics.
    # --------------------------------------------------------

    structures_examined = 0
    pq_tests = 0
    carry_tests = 0

    survivors = []

    true_structure = None

    # --------------------------------------------------------
    # Search anchors.
    # --------------------------------------------------------

    for k0, l0 in anchors:

        # ----------------------------------------------------
        # p interval from anchor.
        # ----------------------------------------------------

        p_anchor = (
            r1_0 * k0,
            r1_0 * (k0 + 1) - 1,
        )

        # ----------------------------------------------------
        # q interval from anchor.
        # ----------------------------------------------------

        q_anchor = (
            r2_0 * l0,
            r2_0 * (l0 + 1) - 1,
        )

        # ----------------------------------------------------
        # Row options.
        # ----------------------------------------------------

        row_options = []

        failed = False

        for i, r1 in enumerate(
            r1_grid
        ):

            p_lo, p_hi = p_anchor

            k_min = p_lo // r1
            k_max = p_hi // r1

            values = []

            for k in range(
                max(1, k_min),
                k_max + 1,
            ):

                # Anchor column condition.
                E = Q[i][0] - k * l0

                if E < 0:
                    continue

                if E > k + l0 + 10:
                    continue

                values.append(k)

            if not values:

                failed = True
                break

            row_options.append(
                values
            )

        if failed:
            continue

        # ----------------------------------------------------
        # Column options.
        # ----------------------------------------------------

        col_options = []

        for j, r2 in enumerate(
            r2_grid
        ):

            q_lo, q_hi = q_anchor

            l_min = q_lo // r2
            l_max = q_hi // r2

            values = []

            for l in range(
                max(1, l_min),
                l_max + 1,
            ):

                E = Q[0][j] - k0 * l

                if E < 0:
                    continue

                if E > k0 + l + 10:
                    continue

                values.append(l)

            if not values:

                failed = True
                break

            col_options.append(
                values
            )

        if failed:
            continue

        # ----------------------------------------------------
        # Enumerate compatible row vectors.
        # ----------------------------------------------------

        for row_tail in itertools.product(
            *row_options[1:]
        ):

            ks = (
                k0,
                *row_tail,
            )

            p_interval = common_interval(
                r1_grid,
                ks,
            )

            if p_interval is None:
                continue

            # ------------------------------------------------
            # Enumerate compatible column vectors.
            # ------------------------------------------------

            for col_tail in itertools.product(
                *col_options[1:]
            ):

                ls = (
                    l0,
                    *col_tail,
                )

                q_interval = common_interval(
                    r2_grid,
                    ls,
                )

                if q_interval is None:
                    continue

                # ------------------------------------------------
                # Construct E matrix.
                # ------------------------------------------------

                E_matrix = []

                valid = True

                for i in range(rows):

                    E_row = []

                    for j in range(cols):

                        K = (
                            ks[i]
                            * ls[j]
                        )

                        E = (
                            Q[i][j]
                            - K
                        )

                        if E < 0:
                            valid = False
                            break

                        if E > (
                            ks[i]
                            + ls[j]
                            + 10
                        ):
                            valid = False
                            break

                        E_row.append(E)

                    if not valid:
                        break

                    E_matrix.append(
                        E_row
                    )

                if not valid:
                    continue

                structures_examined += 1

                if (
                    structures_examined
                    > MAX_STRUCTURES
                ):

                    return {
                        "status": "LIMIT",
                        "anchors": len(anchors),
                        "structures":
                            structures_examined,
                        "pq_tests":
                            pq_tests,
                        "carry_tests":
                            carry_tests,
                        "survivors":
                            survivors,
                        "true":
                            true_structure,
                    }

                # ------------------------------------------------
                # Enumerate actual p,q in the common intervals.
                # ------------------------------------------------

                p_lo, p_hi = p_interval
                q_lo, q_hi = q_interval

                for p in range(
                    p_lo,
                    p_hi + 1,
                ):

                    for q in range(
                        q_lo,
                        q_hi + 1,
                    ):

                        pq_tests += 1

                        # ----------------------------------------
                        # Exact carry structure.
                        # ----------------------------------------

                        ok = True

                        for i, r1 in enumerate(
                            r1_grid
                        ):

                            k = ks[i]

                            a = (
                                p
                                - r1 * k
                            )

                            if not (
                                0 <= a < r1
                            ):

                                ok = False
                                break

                        if not ok:
                            continue

                        for j, r2 in enumerate(
                            r2_grid
                        ):

                            l = ls[j]

                            b = (
                                q
                                - r2 * l
                            )

                            if not (
                                0 <= b < r2
                            ):

                                ok = False
                                break

                        if not ok:
                            continue

                        # ----------------------------------------
                        # All cells.
                        # ----------------------------------------

                        for i, r1 in enumerate(
                            r1_grid
                        ):

                            k = ks[i]

                            for j, r2 in enumerate(
                                r2_grid
                            ):

                                l = ls[j]

                                a = (
                                    p
                                    - r1 * k
                                )

                                b = (
                                    q
                                    - r2 * l
                                )

                                c1 = (
                                    k * b
                                ) // r2

                                c2 = (
                                    l * a
                                ) // r1

                                beta = (
                                    k * b
                                ) % r2

                                alpha = (
                                    l * a
                                ) % r1

                                R = (
                                    r1 * r2
                                )

                                c3 = (
                                    r1 * beta
                                    + r2 * alpha
                                    + a * b
                                ) // R

                                E_actual = (
                                    c1
                                    + c2
                                    + c3
                                )

                                carry_tests += 1

                                if (
                                    E_actual
                                    != E_matrix[i][j]
                                ):

                                    ok = False
                                    break

                            if not ok:
                                break

                        if not ok:
                            continue

                        # ----------------------------------------
                        # Complete structure survived.
                        # ----------------------------------------

                        exact = (
                            p * q == n
                        )

                        solution = {
                            "p": p,
                            "q": q,
                            "k": tuple(ks),
                            "l": tuple(ls),
                            "p_interval":
                                p_interval,
                            "q_interval":
                                q_interval,
                            "E":
                                E_matrix,
                            "exact":
                                exact,
                        }

                        survivors.append(
                            solution
                        )

                        # ------------------------------------------------
                        # Check true factorization.
                        # ------------------------------------------------

                        if (
                            (
                                p == p_true
                                and q == q_true
                            )
                            or
                            (
                                p == q_true
                                and q == p_true
                            )
                        ):

                            true_structure = (
                                solution
                            )

    return {
        "status": "OK",
        "anchors": len(anchors),
        "structures":
            structures_examined,
        "pq_tests":
            pq_tests,
        "carry_tests":
            carry_tests,
        "survivors":
            survivors,
        "true":
            true_structure,
    }


# ------------------------------------------------------------
# Print compact results
# ------------------------------------------------------------

def print_result(
    grid_size,
    result,
):

    print()
    print(
        "-" * 78
    )

    print(
        f"GRID {grid_size}x{grid_size}"
    )

    print(
        f"status               = "
        f"{result['status']}"
    )

    print(
        f"anchor candidates    = "
        f"{result['anchors']:,}"
    )

    print(
        f"structures examined  = "
        f"{result['structures']:,}"
    )

    print(
        f"p/q combinations     = "
        f"{result['pq_tests']:,}"
    )

    print(
        f"carry cell tests     = "
        f"{result['carry_tests']:,}"
    )

    print(
        f"surviving structures = "
        f"{len(result['survivors']):,}"
    )

    print(
        f"TRUE structure found = "
        f"{result['true'] is not None}"
    )

    exact = [
        s
        for s in result["survivors"]
        if s["exact"]
    ]

    print(
        f"exact p*q==n         = "
        f"{len(exact):,}"
    )

    if result["survivors"]:

        print()
        print(
            "survivor examples:"
        )

        for s in result[
            "survivors"
        ][:PRINT_SURVIVORS]:

            print(
                "    "
                f"p={s['p']} "
                f"q={s['q']} "
                f"exact={s['exact']}"
            )

    if result["status"] == "LIMIT":

        print()
        print(
            "WARNING:"
        )

        print(
            "This grid exceeded "
            "MAX_STRUCTURES."
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "2D GRID CARRY-COLLAPSE SCALING"
    )
    print(
        "EXPERIMENT 95"
    )
    print("=" * 78)

    rng = None

    # --------------------------------------------------------
    # Generate fixed test numbers so every grid size for a case
    # sees exactly the SAME n.
    # --------------------------------------------------------

    cases = []

    for index in range(
        CASE_COUNT
    ):

        p, q, n = (
            generate_semiprime(
                TARGET
            )
        )

        cases.append(
            (p, q, n)
        )

    # --------------------------------------------------------
    # Run.
    # --------------------------------------------------------

    summary = []

    for case_index, (
        p_true,
        q_true,
        n,
    ) in enumerate(
        cases,
        start=1,
    ):

        print()
        print("=" * 78)

        print(
            f"CASE {case_index}"
        )

        print(
            f"n      = {n}"
        )

        print(
            f"p      = {p_true}"
        )

        print(
            f"q      = {q_true}"
        )

        print(
            f"|p-q|  = "
            f"{abs(p_true-q_true)}"
        )

        print(
            "=" * 78
        )

        for grid_size in GRID_SIZES:

            # ------------------------------------------------
            # Construct grids.
            # ------------------------------------------------

            r1_grid = make_prime_grid(
                R1_CENTER,
                grid_size,
                GRID_SPREAD,
            )

            r2_grid = make_prime_grid(
                R2_CENTER,
                grid_size,
                GRID_SPREAD,
            )

            # Avoid accidental duplicate primes.
            if (
                len(r1_grid)
                < grid_size
                or
                len(r2_grid)
                < grid_size
            ):

                print()
                print(
                    f"GRID {grid_size}x{grid_size}"
                )

                print(
                    "Could not construct enough "
                    "distinct primes."
                )

                continue

            print()
            print(
                f"Testing grid "
                f"{grid_size}x{grid_size}"
            )

            print(
                f"r1 = {r1_grid}"
            )

            print(
                f"r2 = {r2_grid}"
            )

            start = time.perf_counter()

            result = search_grid(
                n,
                p_true,
                q_true,
                r1_grid,
                r2_grid,
            )

            elapsed = (
                time.perf_counter()
                - start
            )

            print_result(
                grid_size,
                result,
            )

            print(
                f"time                 = "
                f"{elapsed:.6f} s"
            )

            summary.append(
                (
                    case_index,
                    grid_size,
                    result,
                    elapsed,
                )
            )

    # --------------------------------------------------------
    # Final scaling table.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("SCALING SUMMARY")
    print("=" * 78)

    print()
    print(
        "case grid  anchors    structures    "
        "survivors   exact   true   time"
    )

    print(
        "-" * 78
    )

    for (
        case_index,
        grid_size,
        result,
        elapsed,
    ) in summary:

        exact_count = sum(
            1
            for s in result["survivors"]
            if s["exact"]
        )

        print(
            f"{case_index:>4} "
            f"{grid_size}x{grid_size} "
            f"{result['anchors']:>9,} "
            f"{result['structures']:>12,} "
            f"{len(result['survivors']):>10,} "
            f"{exact_count:>7,} "
            f"{str(result['true'] is not None):>6} "
            f"{elapsed:>8.3f}s"
        )

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 95 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 95
# ============================================================
