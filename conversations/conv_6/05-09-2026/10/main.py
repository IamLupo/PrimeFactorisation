import math
import time
import itertools
import sympy


# ============================================================
# START EXPERIMENT 93
# ============================================================
#
# 2D GRID:
#
#       Q_ij = floor(n / (r1_i*r2_j))
#
#       Q_ij = k_i*l_j + E_ij
#
# with:
#
#       k_i = floor(p/r1_i)
#       l_j = floor(q/r2_j)
#
#
# We DO NOT search a,b.
#
# We use three structural constraints:
#
#   1. K_ij = k_i*l_j
#
#   2. Same p for every row:
#
#          r1_i*k_i <= p < r1_i*(k_i+1)
#
#      therefore all p-intervals must overlap.
#
#   3. Same q for every column:
#
#          r2_j*l_j <= q < r2_j*(l_j+1)
#
#      therefore all q-intervals must overlap.
#
#
# We then test the carry-size condition:
#
#       E_ij = Q_ij - k_i*l_j
#
# and require E_ij >= 0.
#
# We do NOT require p*q == n until FINAL REPORTING.
#
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

R1_GRID = [7, 11, 17]
R2_GRID = [11, 13, 17]

# How many possible k/l values around the natural scale
# are considered for the anchor cell.
#
# This is only a search-window parameter for the experiment.
ANCHOR_WINDOW = 1000

# Extra E allowance used while pruning.
#
# E is constrained by:
#
#     E = Q - k*l >= 0
#
# We additionally use the theoretical rough bound:
#
#     E < k + l + 3
#
# so a candidate cell must satisfy:
#
#     0 <= E < k+l+3
#
E_EXTRA = 10

PRINT_SOLUTIONS = 20


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
# Given values r_i and quotient k_i:
#
#     r_i*k_i <= p <= r_i*(k_i+1)-1
#
# Return the intersection.
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

        if this_lo > lo:
            lo = this_lo

        if hi is None or this_hi < hi:
            hi = this_hi

    if hi is None:
        return None

    if lo > hi:
        return None

    return lo, hi


# ------------------------------------------------------------
# E calculation
# ------------------------------------------------------------

def compute_E_from_kl(
    n,
    r1,
    r2,
    k,
    l,
):

    Q = n // (r1 * r2)

    K = k * l

    E = Q - K

    return K, E, Q


# ------------------------------------------------------------
# Carry bound
#
# Exact E is nonnegative.
#
# Since:
#
#     c1 <= k-1
#     c2 <= l-1
#
# and c3 is small, we use a conservative bound:
#
#     E < k+l+3
#
# ------------------------------------------------------------

def e_is_plausible(
    E,
    k,
    l,
):

    if E < 0:
        return False

    return E < (k + l + 3 + E_EXTRA)


# ------------------------------------------------------------
# Build Q matrix
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
# Find possible k values for one row from a chosen p interval.
#
# If we know possible p range [plo, phi], then:
#
#     k = floor(p/r)
#
# can only take a small interval of values.
# ------------------------------------------------------------

def k_values_from_p_interval(
    r,
    plo,
    phi,
):

    k_min = plo // r
    k_max = phi // r

    return range(
        k_min,
        k_max + 1,
    )


# ------------------------------------------------------------
# Find possible l values.
# ------------------------------------------------------------

def l_values_from_q_interval(
    r,
    qlo,
    qhi,
):

    l_min = qlo // r
    l_max = qhi // r

    return range(
        l_min,
        l_max + 1,
    )


# ------------------------------------------------------------
# Validate complete structure.
#
# This calculates E matrix and applies all constraints.
# ------------------------------------------------------------

def validate_structure(
    n,
    r1_grid,
    r2_grid,
    ks,
    ls,
):

    rows = len(r1_grid)
    cols = len(r2_grid)

    # --------------------------------------------------------
    # Shared p interval.
    # --------------------------------------------------------

    p_interval = common_interval(
        r1_grid,
        ks,
    )

    if p_interval is None:
        return None

    # --------------------------------------------------------
    # Shared q interval.
    # --------------------------------------------------------

    q_interval = common_interval(
        r2_grid,
        ls,
    )

    if q_interval is None:
        return None

    K_matrix = []
    E_matrix = {}

    # --------------------------------------------------------
    # Every cell.
    # --------------------------------------------------------

    for i in range(rows):

        K_row = []

        for j in range(cols):

            r1 = r1_grid[i]
            r2 = r2_grid[j]

            k = ks[i]
            l = ls[j]

            K = k * l

            Q = n // (r1 * r2)

            E = Q - K

            if not e_is_plausible(
                E,
                k,
                l,
            ):
                return None

            K_row.append(K)

            E_matrix[(i, j)] = E

        K_matrix.append(
            K_row
        )

    return {
        "ks": tuple(ks),
        "ls": tuple(ls),
        "p_interval": p_interval,
        "q_interval": q_interval,
        "K": K_matrix,
        "E": E_matrix,
    }


# ------------------------------------------------------------
# Search the structured grid.
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

    Q = build_Q_matrix(
        n,
        r1_grid,
        r2_grid,
    )

    # --------------------------------------------------------
    # Anchor cell.
    #
    # We search k0,l0 around sqrt(n)/(r1*r2 scale).
    #
    # IMPORTANT:
    # This does not assume k ~= l.
    # k0 and l0 are searched independently.
    # --------------------------------------------------------

    r1_0 = r1_grid[0]
    r2_0 = r2_grid[0]

    p_scale = math.isqrt(n) // r1_0
    q_scale = math.isqrt(n) // r2_0

    k_min = max(
        1,
        p_scale - ANCHOR_WINDOW,
    )

    k_max = (
        p_scale
        + ANCHOR_WINDOW
    )

    l_min = max(
        1,
        q_scale - ANCHOR_WINDOW,
    )

    l_max = (
        q_scale
        + ANCHOR_WINDOW
    )

    # --------------------------------------------------------
    # Build possible anchor (k,l) pairs.
    #
    # We immediately use Q00 and E bounds.
    # --------------------------------------------------------

    anchors = []

    Q00 = Q[0][0]

    for k0 in range(
        k_min,
        k_max + 1,
    ):

        # We can derive a rough l interval from:
        #
        #     E >= 0
        #
        # so:
        #
        #     k*l <= Q00
        #
        l_upper = Q00 // k0

        upper = min(
            l_max,
            l_upper,
        )

        for l0 in range(
            l_min,
            upper + 1,
        ):

            K0 = k0 * l0
            E0 = Q00 - K0

            if not e_is_plausible(
                E0,
                k0,
                l0,
            ):
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

    # --------------------------------------------------------
    # For each anchor, derive possible row k values and
    # column l values using interval overlap.
    # --------------------------------------------------------

    for anchor_index, (
        k0,
        l0,
    ) in enumerate(anchors):

        # ----------------------------------------------------
        # We don't know p exactly.
        #
        # Anchor gives:
        #
        # r1_0*k0 <= p < r1_0*(k0+1)
        #
        # and:
        #
        # r2_0*l0 <= q < r2_0*(l0+1)
        # ----------------------------------------------------

        p_anchor = (
            r1_0 * k0,
            r1_0 * (k0 + 1) - 1,
        )

        q_anchor = (
            r2_0 * l0,
            r2_0 * (l0 + 1) - 1,
        )

        # ----------------------------------------------------
        # Possible row k values.
        # ----------------------------------------------------

        row_options = [
            None
            for _ in range(rows)
        ]

        row_options[0] = [k0]

        row_failed = False

        for i in range(
            1,
            rows,
        ):

            vals = list(
                k_values_from_p_interval(
                    r1_grid[i],
                    p_anchor[0],
                    p_anchor[1],
                )
            )

            vals = [
                x
                for x in vals
                if x > 0
            ]

            if not vals:

                row_failed = True
                break

            row_options[i] = vals

        if row_failed:
            continue

        # ----------------------------------------------------
        # Possible column l values.
        # ----------------------------------------------------

        col_options = [
            None
            for _ in range(cols)
        ]

        col_options[0] = [l0]

        col_failed = False

        for j in range(
            1,
            cols,
        ):

            vals = list(
                l_values_from_q_interval(
                    r2_grid[j],
                    q_anchor[0],
                    q_anchor[1],
                )
            )

            vals = [
                x
                for x in vals
                if x > 0
            ]

            if not vals:

                col_failed = True
                break

            col_options[j] = vals

        if col_failed:
            continue

        # ----------------------------------------------------
        # Enumerate remaining possibilities.
        #
        # Since all r's are small and close, these option
        # lists should be tiny.
        # ----------------------------------------------------

        for ks_tail in itertools.product(
            *row_options[1:]
        ):

            ks = (
                k0,
                *ks_tail,
            )

            # Recompute common p interval.
            p_interval = common_interval(
                r1_grid,
                ks,
            )

            if p_interval is None:
                continue

            for ls_tail in itertools.product(
                *col_options[1:]
            ):

                ls = (
                    l0,
                    *ls_tail,
                )

                q_interval = common_interval(
                    r2_grid,
                    ls,
                )

                if q_interval is None:
                    continue

                result = validate_structure(
                    n,
                    r1_grid,
                    r2_grid,
                    ks,
                    ls,
                )

                if result is None:
                    continue

                # ------------------------------------------------
                # Record whether the TRUE quotient structure
                # is represented.
                # ------------------------------------------------

                true_ks = tuple(
                    p_true // r
                    for r in r1_grid
                )

                true_ls = tuple(
                    q_true // r
                    for r in r2_grid
                )

                result["true_structure"] = (
                    result["ks"] == true_ks
                    and result["ls"] == true_ls
                )

                solutions.append(
                    result
                )

                if len(solutions) >= PRINT_SOLUTIONS:
                    return solutions

    return solutions


# ------------------------------------------------------------
# Print matrix
# ------------------------------------------------------------

def print_K_matrix(matrix):

    for row in matrix:

        print(
            "    "
            + " ".join(
                f"{x:>10}"
                for x in row
            )
        )


# ------------------------------------------------------------
# Print E matrix
# ------------------------------------------------------------

def print_E_matrix(
    E,
    rows,
    cols,
):

    for i in range(rows):

        print(
            "    "
            + " ".join(
                f"{E[(i,j)]:>10}"
                for j in range(cols)
            )
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print("2D SHARED-p / SHARED-q K-E STRUCTURE")
    print("EXPERIMENT 93")
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
    # Grid.
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
    # TRUE structure.
    # --------------------------------------------------------

    true_ks = [
        p_true // r
        for r in r1_grid
    ]

    true_ls = [
        q_true // r
        for r in r2_grid
    ]

    print()
    print(
        f"TRUE k vector = {true_ks}"
    )

    print(
        f"TRUE l vector = {true_ls}"
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    start = time.perf_counter()

    solutions = search_structures(
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
        f"solutions found = "
        f"{len(solutions)}"
    )

    print(
        f"search time     = "
        f"{elapsed:.6f} s"
    )

    true_survives = any(
        s["true_structure"]
        for s in solutions
    )

    print(
        f"TRUE structure survived = "
        f"{true_survives}"
    )

    # --------------------------------------------------------
    # Print structures.
    # --------------------------------------------------------

    for index, s in enumerate(
        solutions,
        start=1,
    ):

        print()
        print(
            f"SOLUTION {index}"
        )

        print(
            f"    k = {s['ks']}"
        )

        print(
            f"    l = {s['ls']}"
        )

        print(
            f"    p interval = "
            f"{s['p_interval']}"
        )

        print(
            f"    q interval = "
            f"{s['q_interval']}"
        )

        print()
        print("    K matrix:")

        print_K_matrix(
            s["K"]
        )

        print()
        print("    E matrix:")

        print_E_matrix(
            s["E"],
            len(r1_grid),
            len(r2_grid),
        )

        print(
            f"    TRUE structure = "
            f"{s['true_structure']}"
        )

    # --------------------------------------------------------
    # Exact factor check ONLY as final measurement.
    #
    # We derive possible p and q from common intervals.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL p/q INTERVAL MEASUREMENT")
    print("=" * 78)

    for index, s in enumerate(
        solutions,
        start=1,
    ):

        plo, phi = s["p_interval"]
        qlo, qhi = s["q_interval"]

        print()
        print(
            f"SOLUTION {index}"
        )

        print(
            f"    possible p: "
            f"{plo} .. {phi}"
        )

        print(
            f"    possible q: "
            f"{qlo} .. {qhi}"
        )

        true_p_inside = (
            plo <= p_true <= phi
        )

        true_q_inside = (
            qlo <= q_true <= qhi
        )

        print(
            f"    true p inside = "
            f"{true_p_inside}"
        )

        print(
            f"    true q inside = "
            f"{true_q_inside}"
        )

        # Check whether the p/q intervals contain an exact
        # factorization pair. This is ONLY final reporting.

        exact_pairs = []

        for candidate_p in range(
            plo,
            min(
                phi,
                plo + 5000,
            ) + 1,
        ):

            if n % candidate_p != 0:
                continue

            candidate_q = n // candidate_p

            if (
                qlo
                <= candidate_q
                <= qhi
            ):

                exact_pairs.append(
                    (
                        candidate_p,
                        candidate_q,
                    )
                )

        if exact_pairs:

            print(
                "    exact pair(s) found "
                "inside intervals:"
            )

            for pair in exact_pairs[:10]:

                print(
                    f"        {pair}"
                )

    print()
    print("=" * 78)
    print("EXPERIMENT 93 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 93
# ============================================================
