import math
import time
import itertools
import sympy


# ============================================================
# START EXPERIMENT 94
# ============================================================
#
# 2D K/E + SHARED-p + SHARED-q + EXACT CARRY CONSISTENCY
#
#
# For a grid:
#
#     r1_i
#     r2_j
#
# define:
#
#     Q_ij = floor(n / (r1_i*r2_j))
#
# and:
#
#     Q_ij = K_ij + E_ij
#
# with:
#
#     K_ij = k_i*l_j
#
#
# Experiment 93 tested:
#
#     K_ij = k_i*l_j
#
# and common p/q intervals.
#
# Experiment 94 additionally tests the ORIGINAL carry
# equations across every grid cell.
#
#
# For a candidate:
#
#     p = r1_i*k_i + a_i
#     q = r2_j*l_j + b_j
#
# where:
#
#     a_i = p - r1_i*k_i
#     b_j = q - r2_j*l_j
#
# We require:
#
#     c1 = floor(k_i*b_j/r2_j)
#     c2 = floor(l_j*a_i/r1_i)
#
#     beta  = (k_i*b_j) mod r2_j
#     alpha = (l_j*a_i) mod r1_i
#
#     c3 = floor(
#         (r1_i*beta + r2_j*alpha + a_i*b_j)
#         / (r1_i*r2_j)
#     )
#
# and:
#
#     E_ij = c1+c2+c3
#
#
# IMPORTANT:
#
# We do NOT use p*q == n as a search filter.
#
# The exact p*q == n check is ONLY performed at the end
# to classify the surviving mathematical structures.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

R1_GRID = [7, 11, 17]
R2_GRID = [11, 13, 17]

# Search window around sqrt(n)/r for anchor k,l.
#
# Increase this if the true structure falls outside.
ANCHOR_WINDOW = 1500

# Maximum number of structures PRINTED.
# Search itself continues.
PRINT_SOLUTIONS = 20

# Maximum number of p/q values explored per structure.
#
# This is only a safety guard for pathological intervals.
MAX_P_INTERVAL = 10000
MAX_Q_INTERVAL = 10000


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
# Interval intersection
#
# r*k <= x < r*(k+1)
#
# Therefore:
#
# r*k <= x <= r*(k+1)-1
# ------------------------------------------------------------

def common_interval(rs, ks):

    lo = 0
    hi = None

    for r, k in zip(rs, ks):

        this_lo = r * k
        this_hi = r * (k + 1) - 1

        lo = max(lo, this_lo)

        if hi is None:
            hi = this_hi
        else:
            hi = min(hi, this_hi)

    if hi is None:
        return None

    if lo > hi:
        return None

    return lo, hi


# ------------------------------------------------------------
# Compute exact carry state from p,q at one cell
# ------------------------------------------------------------

def exact_cell_state(
    n,
    p,
    q,
    r1,
    r2,
):

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

    Q = n // R

    return {
        "k": k,
        "l": l,
        "a": a,
        "b": b,
        "K": K,
        "E": E,
        "Q": Q,
        "c1": c1,
        "c2": c2,
        "c3": c3,
    }


# ------------------------------------------------------------
# Validate a proposed k/l structure against actual p,q
# ------------------------------------------------------------

def validate_pq_against_structure(
    n,
    p,
    q,
    r1_grid,
    r2_grid,
    expected_ks,
    expected_ls,
    expected_E,
):

    rows = len(r1_grid)
    cols = len(r2_grid)

    for i, r1 in enumerate(r1_grid):

        k, a = divmod(p, r1)

        if k != expected_ks[i]:
            return False

        if not (0 <= a < r1):
            return False

    for j, r2 in enumerate(r2_grid):

        l, b = divmod(q, r2)

        if l != expected_ls[j]:
            return False

        if not (0 <= b < r2):
            return False

    # --------------------------------------------------------
    # Exact carry validation for every cell.
    # --------------------------------------------------------

    for i, r1 in enumerate(r1_grid):

        k = expected_ks[i]
        a = p - r1 * k

        for j, r2 in enumerate(r2_grid):

            l = expected_ls[j]
            b = q - r2 * l

            state = exact_cell_state(
                n,
                p,
                q,
                r1,
                r2,
            )

            if state["K"] != k * l:
                return False

            if state["E"] != expected_E[i][j]:
                return False

    return True


# ------------------------------------------------------------
# Build Q matrix
# ------------------------------------------------------------

def build_Q(
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
# Build cell (k,l,E) candidate lists.
#
# We DON'T need to search arbitrary E independently forever.
#
# For a given Q:
#
#     K = Q-E
#
# and:
#
#     E = Q-k*l
#
# We search E up to a computed bound.
# ------------------------------------------------------------

def cell_candidates(
    n,
    r1,
    r2,
):

    Q = n // (r1 * r2)

    # Conservative experimental bound.
    #
    # k and l for balanced factors are on the order of sqrt(Q).
    #
    # This is deliberately generous.
    e_cap = (
        4 * math.isqrt(Q)
        + 100
    )

    result = []

    for E in range(
        min(e_cap, Q - 1) + 1
    ):

        K = Q - E

        if K <= 0:
            continue

        for k in sympy.divisors(K):

            l = K // k

            # Rough necessary carry bound.
            if E >= k + l + 5:
                continue

            result.append(
                {
                    "k": k,
                    "l": l,
                    "K": K,
                    "E": E,
                }
            )

    return result


# ------------------------------------------------------------
# Search global k/l structures
# ------------------------------------------------------------

def search_structures(
    n,
    r1_grid,
    r2_grid,
    p_true,
    q_true,
):

    rows = len(r1_grid)
    cols = len(r2_grid)

    Q = build_Q(
        n,
        r1_grid,
        r2_grid,
    )

    # --------------------------------------------------------
    # Anchor cell.
    # --------------------------------------------------------

    r1_0 = r1_grid[0]
    r2_0 = r2_grid[0]

    p_scale = math.isqrt(n) // r1_0
    q_scale = math.isqrt(n) // r2_0

    k_lo = max(
        1,
        p_scale - ANCHOR_WINDOW,
    )

    k_hi = (
        p_scale
        + ANCHOR_WINDOW
    )

    l_lo = max(
        1,
        q_scale - ANCHOR_WINDOW,
    )

    l_hi = (
        q_scale
        + ANCHOR_WINDOW
    )

    # --------------------------------------------------------
    # Anchor candidates.
    #
    # Instead of iterating every E and factoring Q-E,
    # simply search k,l in the anchor window.
    #
    # The anchor E is then determined.
    # --------------------------------------------------------

    anchors = []

    Q00 = Q[0][0]

    for k0 in range(
        k_lo,
        k_hi + 1,
    ):

        # Since E >= 0:
        #
        # k0*l0 <= Q00.
        #
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

            if E0 >= k0 + l0 + 5:
                continue

            anchors.append(
                (k0, l0)
            )

    print()
    print(
        f"anchor candidates = "
        f"{len(anchors):,}"
    )

    solutions = []

    structures_examined = 0
    pq_tests = 0
    carry_tests = 0

    # --------------------------------------------------------
    # Anchor expansion
    # --------------------------------------------------------

    for k0, l0 in anchors:

        p_anchor = (
            r1_0 * k0,
            r1_0 * (k0 + 1) - 1,
        )

        q_anchor = (
            r2_0 * l0,
            r2_0 * (l0 + 1) - 1,
        )

        # ----------------------------------------------------
        # Derive row quotient options.
        # ----------------------------------------------------

        row_options = []

        possible = True

        for i, r1 in enumerate(
            r1_grid
        ):

            values = []

            pmin, pmax = p_anchor

            kmin = pmin // r1
            kmax = pmax // r1

            for k in range(
                max(1, kmin),
                kmax + 1,
            ):

                # Necessary condition from cell (i,0):
                #
                # E = Q-k*l >= 0
                #
                l = l0

                E = Q[i][0] - k * l

                if E < 0:
                    continue

                if E >= k + l + 5:
                    continue

                values.append(k)

            if not values:

                possible = False
                break

            row_options.append(
                values
            )

        if not possible:
            continue

        # ----------------------------------------------------
        # Derive column quotient options.
        # ----------------------------------------------------

        col_options = []

        for j, r2 in enumerate(
            r2_grid
        ):

            values = []

            qmin, qmax = q_anchor

            lmin = qmin // r2
            lmax = qmax // r2

            for l in range(
                max(1, lmin),
                lmax + 1,
            ):

                k = k0

                E = Q[0][j] - k * l

                if E < 0:
                    continue

                if E >= k + l + 5:
                    continue

                values.append(l)

            if not values:

                possible = False
                break

            col_options.append(
                values
            )

        if not possible:
            continue

        # ----------------------------------------------------
        # Enumerate quotient vectors.
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

            p_lo, p_hi = p_interval

            if (
                p_hi - p_lo + 1
                > MAX_P_INTERVAL
            ):
                continue

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

                q_lo, q_hi = q_interval

                if (
                    q_hi - q_lo + 1
                    > MAX_Q_INTERVAL
                ):
                    continue

                # ------------------------------------------------
                # Calculate complete E matrix.
                #
                # This is the first strong grid-level filter.
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

                        E = Q[i][j] - K

                        if E < 0:
                            valid = False
                            break

                        if E >= (
                            ks[i]
                            + ls[j]
                            + 5
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

                # ------------------------------------------------
                # NOW test actual p/q values inside the common
                # intervals.
                #
                # This is where the complete residue/carry
                # structure is enforced.
                # ------------------------------------------------

                for p in range(
                    p_lo,
                    p_hi + 1,
                ):

                    for q in range(
                        q_lo,
                        q_hi + 1,
                    ):

                        pq_tests += 1

                        ok = True

                        # ------------------------------------------------
                        # Check every row/column quotient.
                        # ------------------------------------------------

                        for i, r1 in enumerate(
                            r1_grid
                        ):

                            if p // r1 != ks[i]:
                                ok = False
                                break

                        if not ok:
                            continue

                        for j, r2 in enumerate(
                            r2_grid
                        ):

                            if q // r2 != ls[j]:
                                ok = False
                                break

                        if not ok:
                            continue

                        # ------------------------------------------------
                        # Exact carry equations.
                        # ------------------------------------------------

                        for i, r1 in enumerate(
                            r1_grid
                        ):

                            k = ks[i]
                            a = p - r1 * k

                            for j, r2 in enumerate(
                                r2_grid
                            ):

                                l = ls[j]
                                b = q - r2 * l

                                # Exact E from the original carry
                                # definition.

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

                        # ------------------------------------------------
                        # We found a complete mathematical structure.
                        # ------------------------------------------------

                        exact = (
                            p * q == n
                        )

                        solution = {
                            "p": p,
                            "q": q,
                            "ks": tuple(ks),
                            "ls": tuple(ls),
                            "p_interval": p_interval,
                            "q_interval": q_interval,
                            "K": [
                                [
                                    ks[i] * ls[j]
                                    for j in range(cols)
                                ]
                                for i in range(rows)
                            ],
                            "E": E_matrix,
                            "exact": exact,
                        }

                        # Deduplicate.
                        key = (
                            p,
                            q,
                            tuple(ks),
                            tuple(ls),
                        )

                        already = any(
                            (
                                s["p"],
                                s["q"],
                                s["ks"],
                                s["ls"],
                            )
                            == key
                            for s in solutions
                        )

                        if not already:
                            solutions.append(
                                solution
                            )

                            # Keep searching.
                            # Only printing is limited.

    return (
        solutions,
        {
            "structures_examined":
                structures_examined,
            "pq_tests":
                pq_tests,
            "carry_tests":
                carry_tests,
        },
    )


# ------------------------------------------------------------
# Print matrix
# ------------------------------------------------------------

def print_matrix(
    matrix,
):

    for row in matrix:

        print(
            "      "
            + " ".join(
                f"{x:>10}"
                for x in row
            )
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("2D K/E + SHARED-p/q + EXACT CARRY CONSISTENCY")
    print("EXPERIMENT 94")
    print("=" * 78)

    # --------------------------------------------------------
    # Generate n.
    # --------------------------------------------------------

    p_true, q_true, n = (
        generate_semiprime(
            TARGET
        )
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

    # --------------------------------------------------------
    # Grids.
    # --------------------------------------------------------

    r1_grid = sorted(
        set(R1_GRID)
    )

    r2_grid = sorted(
        set(R2_GRID)
    )

    print()
    print(
        f"r1 grid = {r1_grid}"
    )

    print(
        f"r2 grid = {r2_grid}"
    )

    print(
        f"cells   = "
        f"{len(r1_grid)*len(r2_grid)}"
    )

    # --------------------------------------------------------
    # TRUE quotient structure.
    # --------------------------------------------------------

    true_ks = tuple(
        p_true // r
        for r in r1_grid
    )

    true_ls = tuple(
        q_true // r
        for r in r2_grid
    )

    print()
    print(
        f"TRUE k = {true_ks}"
    )

    print(
        f"TRUE l = {true_ls}"
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    start = time.perf_counter()

    solutions, stats = search_structures(
        n,
        r1_grid,
        r2_grid,
        p_true,
        q_true,
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    # --------------------------------------------------------
    # Results.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("RESULTS")
    print("=" * 78)

    print(
        f"structures examined = "
        f"{stats['structures_examined']:,}"
    )

    print(
        f"p/q combinations     = "
        f"{stats['pq_tests']:,}"
    )

    print(
        f"carry cell tests     = "
        f"{stats['carry_tests']:,}"
    )

    print(
        f"surviving structures = "
        f"{len(solutions):,}"
    )

    print(
        f"search time          = "
        f"{elapsed:.6f} s"
    )

    true_structure = None

    for s in solutions:

        if (
            s["p"] == p_true
            and s["q"] == q_true
        ) or (
            s["p"] == q_true
            and s["q"] == p_true
        ):

            true_structure = s
            break

    print(
        f"TRUE structure found = "
        f"{true_structure is not None}"
    )

    # --------------------------------------------------------
    # Exact count.
    # --------------------------------------------------------

    exact = [
        s
        for s in solutions
        if s["exact"]
    ]

    print(
        f"exact p*q==n        = "
        f"{len(exact):,}"
    )

    # --------------------------------------------------------
    # Print true structure.
    # --------------------------------------------------------

    if true_structure:

        print()
        print(
            "TRUE STRUCTURE"
        )

        print(
            f"    p = "
            f"{true_structure['p']}"
        )

        print(
            f"    q = "
            f"{true_structure['q']}"
        )

        print(
            f"    k = "
            f"{true_structure['ks']}"
        )

        print(
            f"    l = "
            f"{true_structure['ls']}"
        )

        print()
        print(
            "    K matrix:"
        )

        print_matrix(
            true_structure["K"]
        )

        print()
        print(
            "    E matrix:"
        )

        print_matrix(
            true_structure["E"]
        )

    # --------------------------------------------------------
    # Print surviving false structures.
    # --------------------------------------------------------

    false_count = 0

    for index, s in enumerate(
        solutions,
        start=1,
    ):

        if s is true_structure:
            continue

        if false_count >= PRINT_SOLUTIONS:
            break

        false_count += 1

        print()
        print(
            f"FALSE/OTHER STRUCTURE "
            f"{false_count}"
        )

        print(
            f"    p={s['p']} "
            f"q={s['q']} "
            f"exact={s['exact']}"
        )

        print(
            f"    k={s['ks']}"
        )

        print(
            f"    l={s['ls']}"
        )

        print(
            f"    p interval="
            f"{s['p_interval']}"
        )

        print(
            f"    q interval="
            f"{s['q_interval']}"
        )

    # --------------------------------------------------------
    # Final conclusion.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("EXPERIMENT 94 COMPLETE")
    print("=" * 78)

    if (
        true_structure
        and len(exact) <= 2
    ):

        print()
        print(
            "RESULT:"
        )

        print(
            "The full 2D carry consistency "
            "appears to collapse the tested"
        )

        print(
            "K/E structures to the true "
            "factorization (up to orientation)."
        )

    elif true_structure:

        print()
        print(
            "RESULT:"
        )

        print(
            "The TRUE structure survives the "
            "full carry constraints,"
        )

        print(
            "but additional false structures "
            "also survive."
        )

    else:

        print()
        print(
            "RESULT:"
        )

        print(
            "The TRUE structure was not found "
            "inside the current anchor/search window."
        )

        print(
            "Increase ANCHOR_WINDOW and rerun."
        )


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 94
# ============================================================
