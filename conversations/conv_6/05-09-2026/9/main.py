import math
import time
import sympy


# ============================================================
# START EXPERIMENT 92
# ============================================================
#
# 2D r1/r2 GRID -> K/E CONSISTENCY
#
# Main idea:
#
#     Q_ij = floor(n / (r1_i * r2_j))
#
#     Q_ij = K_ij + E_ij
#
# and
#
#     K_ij = k_i * l_j
#
# where:
#
#     k_i = floor(p / r1_i)
#     l_j = floor(q / r2_j)
#
#
# Therefore the whole K matrix has rank-1 multiplicative
# structure:
#
#     K_ij * K_mn = K_in * K_mj
#
#
# We try to discover that structure BEFORE touching a,b.
#
#
# For every grid cell:
#
#     K = Q - E
#
# with E in a configurable search range.
#
# For every candidate K we obtain its divisor pairs:
#
#     K = k*l
#
# A valid global solution must use the SAME k_i throughout
# row i and the SAME l_j throughout column j.
#
#
# We do NOT check:
#
#     p*q == n
#
# during the grid search.
#
# We only use the TRUE p,q at the end to measure whether
# the true K/E structure survived.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

# Number of r1 and r2 values.
#
# 3x3 = 9 grid cells.
#
R1_COUNT = 3
R2_COUNT = 3

# Initial scale.
R1_CENTER = 11
R2_CENTER = 13

# Spread around the center.
GRID_SPREAD = 0.35

# Search E from 0 through this value.
#
# IMPORTANT:
# This is the experimental search limit.
# It is NOT claimed to solve arbitrarily large n.
E_CAP = 500

# Maximum number of global solutions printed.
PRINT_SOLUTIONS = 10


# ------------------------------------------------------------
# Prime utilities
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
# TRUE K/E matrix
#
# Used only for measurement.
# ------------------------------------------------------------

def true_matrix(
    p,
    q,
    grid_r1,
    grid_r2,
):

    K = []
    E = []

    for r1 in grid_r1:

        K_row = []
        E_row = []

        k = p // r1

        for r2 in grid_r2:

            l = q // r2

            K_value = k * l

            Q = (
                N_GLOBAL
                // (r1 * r2)
            )

            E_value = Q - K_value

            K_row.append(
                K_value
            )

            E_row.append(
                E_value
            )

        K.append(K_row)
        E.append(E_row)

    return K, E


# ------------------------------------------------------------
# Build cell candidates
#
# For every cell:
#
#     K = Q - E
#
# and factor K.
#
# The result contains possible (k,l) pairs.
# ------------------------------------------------------------

def build_cell_candidates(
    n,
    r1,
    r2,
    e_cap,
):

    Q = n // (r1 * r2)

    cells = []

    for E in range(
        min(e_cap, Q - 1) + 1
    ):

        K = Q - E

        if K <= 0:
            break

        # Every divisor gives:
        #
        #     K = k*l
        #
        #
        # Store all possibilities.
        divisors = sympy.divisors(K)

        for k in divisors:

            l = K // k

            cells.append({
                "K": K,
                "E": E,
                "k": k,
                "l": l,
            })

    return Q, cells


# ------------------------------------------------------------
# Build complete 2D grid
# ------------------------------------------------------------

def build_grid_candidates(
    n,
    grid_r1,
    grid_r2,
    e_cap,
):

    grid = []

    total_cell_candidates = 0

    start = time.perf_counter()

    for r1 in grid_r1:

        row = []

        for r2 in grid_r2:

            Q, candidates = (
                build_cell_candidates(
                    n,
                    r1,
                    r2,
                    e_cap,
                )
            )

            total_cell_candidates += (
                len(candidates)
            )

            row.append({
                "r1": r1,
                "r2": r2,
                "Q": Q,
                "candidates": candidates,
            })

        grid.append(row)

    elapsed = (
        time.perf_counter()
        - start
    )

    return (
        grid,
        total_cell_candidates,
        elapsed,
    )


# ------------------------------------------------------------
# Solve global consistency
#
# Each row must use a single k_i.
# Each column must use a single l_j.
#
# We use the first cell of every row/column
# as an anchor and propagate.
# ------------------------------------------------------------

def solve_grid(
    grid,
    grid_r1,
    grid_r2,
):

    rows = len(grid_r1)
    cols = len(grid_r2)

    solutions = []

    # --------------------------------------------------------
    # For each possible k0,l0 from [0,0], propagate.
    # --------------------------------------------------------

    anchor_candidates = (
        grid[0][0]["candidates"]
    )

    for anchor in anchor_candidates:

        k0 = anchor["k"]
        l0 = anchor["l"]

        row_k = [
            None
            for _ in range(rows)
        ]

        col_l = [
            None
            for _ in range(cols)
        ]

        row_k[0] = k0
        col_l[0] = l0

        # ----------------------------------------------------
        # Derive candidate k values for all rows from
        # column 0.
        #
        # Derive candidate l values for all columns from
        # row 0.
        #
        # We keep sets initially.
        # ----------------------------------------------------

        row_options = [
            None
            for _ in range(rows)
        ]

        col_options = [
            None
            for _ in range(cols)
        ]

        # Row 0 is fixed.
        row_options[0] = {
            k0
        }

        # Column 0 is fixed.
        col_options[0] = {
            l0
        }

        # ----------------------------------------------------
        # Rows from column 0.
        # ----------------------------------------------------

        for i in range(1, rows):

            opts = set()

            for c in grid[i][0]["candidates"]:

                if c["l"] == l0:

                    opts.add(
                        c["k"]
                    )

            row_options[i] = opts

            if not opts:
                break

        else:

            # ------------------------------------------------
            # Columns from row 0.
            # ------------------------------------------------

            for j in range(1, cols):

                opts = set()

                for c in grid[0][j]["candidates"]:

                    if c["k"] == k0:

                        opts.add(
                            c["l"]
                        )

                col_options[j] = opts

                if not opts:
                    break

            else:

                # ------------------------------------------------
                # Cartesian combinations of the remaining
                # row/column possibilities.
                #
                # For each combination, every cell must contain
                # exactly the corresponding (k_i,l_j).
                # ------------------------------------------------

                row_lists = [
                    sorted(
                        row_options[i]
                    )
                    for i in range(rows)
                ]

                col_lists = [
                    sorted(
                        col_options[j]
                    )
                    for j in range(cols)
                ]

                import itertools

                for row_values in itertools.product(
                    *row_lists
                ):

                    for col_values in itertools.product(
                        *col_lists
                    ):

                        valid = True
                        E_matrix = []
                        K_matrix = []

                        for i in range(rows):

                            E_row = []
                            K_row = []

                            for j in range(cols):

                                expected_k = (
                                    row_values[i]
                                )

                                expected_l = (
                                    col_values[j]
                                )

                                matches = [
                                    c
                                    for c
                                    in grid[i][j][
                                        "candidates"
                                    ]
                                    if (
                                        c["k"]
                                        == expected_k
                                        and
                                        c["l"]
                                        == expected_l
                                    )
                                ]

                                if len(matches) != 1:

                                    valid = False
                                    break

                                match = matches[0]

                                K_row.append(
                                    match["K"]
                                )

                                E_row.append(
                                    match["E"]
                                )

                            if not valid:
                                break

                            K_matrix.append(
                                K_row
                            )

                            E_matrix.append(
                                E_row
                            )

                        if not valid:
                            continue

                        solutions.append({
                            "row_k": tuple(
                                row_values
                            ),
                            "col_l": tuple(
                                col_values
                            ),
                            "K": K_matrix,
                            "E": E_matrix,
                        })

                        if (
                            len(solutions)
                            >= PRINT_SOLUTIONS
                        ):

                            return solutions

    return solutions


# ------------------------------------------------------------
# Check true structure
# ------------------------------------------------------------

def true_solution_match(
    solution,
    true_K,
):

    return (
        solution["K"]
        == true_K
    )


# ------------------------------------------------------------
# Print matrix
# ------------------------------------------------------------

def print_matrix(
    name,
    matrix,
):

    print()
    print(name)

    for row in matrix:

        print(
            "    "
            + " ".join(
                f"{x:>10}"
                for x in row
            )
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    global N_GLOBAL

    print("=" * 78)
    print("2D r1/r2 GRID -> K/E CONSISTENCY")
    print("EXPERIMENT 92")
    print("=" * 78)

    # --------------------------------------------------------
    # Generate n
    # --------------------------------------------------------

    p_true, q_true, n = (
        generate_semiprime(
            TARGET
        )
    )

    N_GLOBAL = n

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
    # r grids
    # --------------------------------------------------------

    grid_r1 = make_prime_grid(
        R1_CENTER,
        R1_COUNT,
        GRID_SPREAD,
    )

    grid_r2 = make_prime_grid(
        R2_CENTER,
        R2_COUNT,
        GRID_SPREAD,
    )

    print()
    print(
        f"r1 grid = {grid_r1}"
    )

    print(
        f"r2 grid = {grid_r2}"
    )

    print(
        f"cells    = "
        f"{len(grid_r1) * len(grid_r2)}"
    )

    # --------------------------------------------------------
    # True matrix
    # --------------------------------------------------------

    true_K, true_E = true_matrix(
        p_true,
        q_true,
        grid_r1,
        grid_r2,
    )

    print_matrix(
        "TRUE K MATRIX",
        true_K,
    )

    print_matrix(
        "TRUE E MATRIX",
        true_E,
    )

    true_E_max = max(
        max(row)
        for row in true_E
    )

    print()
    print(
        f"true max E = {true_E_max}"
    )

    print(
        f"E_CAP      = {E_CAP}"
    )

    if true_E_max > E_CAP:

        print()
        print(
            "WARNING:"
        )

        print(
            "The true E is outside E_CAP."
        )

        print(
            "Increase E_CAP for this test."
        )

    # --------------------------------------------------------
    # Build candidate grid
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("BUILD CELL CANDIDATES")
    print("=" * 78)

    grid, total_cell_candidates, build_time = (
        build_grid_candidates(
            n,
            grid_r1,
            grid_r2,
            E_CAP,
        )
    )

    print()
    print(
        f"total cell candidates = "
        f"{total_cell_candidates:,}"
    )

    print(
        f"build time = "
        f"{build_time:.6f} s"
    )

    # --------------------------------------------------------
    # Cell statistics
    # --------------------------------------------------------

    print()
    print("CELL CANDIDATE COUNTS")

    for i, r1 in enumerate(grid_r1):

        counts = []

        for j, r2 in enumerate(grid_r2):

            counts.append(
                len(
                    grid[i][j][
                        "candidates"
                    ]
                )
            )

        print(
            f"    r1={r1:<5}: "
            + " ".join(
                f"{x:>7}"
                for x in counts
            )
        )

    # --------------------------------------------------------
    # Solve consistency
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("GLOBAL K = k_i * l_j CONSISTENCY")
    print("=" * 78)

    start = time.perf_counter()

    solutions = solve_grid(
        grid,
        grid_r1,
        grid_r2,
    )

    solve_time = (
        time.perf_counter()
        - start
    )

    print()
    print(
        f"global solutions found = "
        f"{len(solutions)}"
    )

    print(
        f"solve time              = "
        f"{solve_time:.6f} s"
    )

    # --------------------------------------------------------
    # Does true solution survive?
    # --------------------------------------------------------

    true_matches = [
        s
        for s in solutions
        if true_solution_match(
            s,
            true_K,
        )
    ]

    print(
        f"true K matrix survived  = "
        f"{len(true_matches) > 0}"
    )

    # --------------------------------------------------------
    # Print solutions
    # --------------------------------------------------------

    for index, solution in enumerate(
        solutions[:PRINT_SOLUTIONS],
        start=1,
    ):

        print()
        print(
            f"SOLUTION {index}"
        )

        print(
            "    row k = "
            f"{solution['row_k']}"
        )

        print(
            "    col l = "
            f"{solution['col_l']}"
        )

        print_matrix(
            "    K",
            solution["K"],
        )

        print_matrix(
            "    E",
            solution["E"],
        )

    # --------------------------------------------------------
    # Final
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("EXPERIMENT 92 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()


# ============================================================
# FINISHED EXPERIMENT 92
# ============================================================