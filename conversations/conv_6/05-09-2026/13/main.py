import math
import time
import sympy


# ============================================================
# START EXPERIMENT 96
# ============================================================
#
# 2x2 GRID ANCHOR COLLAPSE
#
# Goal:
#
#   Reduce the enormous (k0,l0) anchor search BEFORE
#   constructing full K/E structures.
#
#
# Grid:
#
#     r1 = [r10, r11]
#     r2 = [r20, r21]
#
#
# We have:
#
#     Qij = floor(n/(r1_i*r2_j))
#
# and
#
#     Qij = ki*lj + Eij
#
#
# where:
#
#     k0 = floor(p/r10)
#     k1 = floor(p/r11)
#
#     l0 = floor(q/r20)
#     l1 = floor(q/r21)
#
#
# The experiment progressively filters anchors:
#
#   Stage 1:
#       anchor cell only
#
#   Stage 2:
#       use row 0
#
#   Stage 3:
#       use column 0
#
#   Stage 4:
#       use the complete 2x2 K/E structure
#
#   Stage 5:
#       exact carry consistency
#
#
# We DO NOT use p*q == n for filtering.
#
# Exact multiplication is only reported at the end.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

# Strong 2x2 grid from Experiment 95.
R1_GRID = [3, 11]
R2_GRID = [5, 17]

# Search radius around the natural quotient estimates.
ANCHOR_WINDOW = 3000

# Conservative carry allowance.
E_EXTRA = 10

PRINT_SOLUTIONS = 20


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
# Exact E from p,q
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
# Carry plausibility
#
# Conservative:
#
#     E < k+l+3
# ------------------------------------------------------------

def plausible_E(E, k, l):

    if E < 0:
        return False

    return (
        E
        <= k + l + 3 + E_EXTRA
    )


# ------------------------------------------------------------
# Interval for quotient
# ------------------------------------------------------------

def quotient_interval(
    r,
    k,
):

    return (
        r * k,
        r * (k + 1) - 1,
    )


# ------------------------------------------------------------
# Intersect intervals
# ------------------------------------------------------------

def intersect_intervals(
    intervals,
):

    lo = max(
        x[0]
        for x in intervals
    )

    hi = min(
        x[1]
        for x in intervals
    )

    if lo > hi:
        return None

    return (
        lo,
        hi,
    )


# ------------------------------------------------------------
# Candidate E values for one K
#
# Q = K+E
# ------------------------------------------------------------

def cell_state(
    n,
    r1,
    r2,
    k,
    l,
):

    Q = n // (r1 * r2)

    K = k * l

    E = Q - K

    return Q, K, E


# ------------------------------------------------------------
# Full 2x2 structure validation
# ------------------------------------------------------------

def validate_2x2_structure(
    n,
    r1,
    r2,
    ks,
    ls,
):

    K = []
    E = []

    for i in range(2):

        Krow = []
        Erow = []

        for j in range(2):

            Q, kval, Eval = cell_state(
                n,
                r1[i],
                r2[j],
                ks[i],
                ls[j],
            )

            if Eval < 0:
                return None

            if not plausible_E(
                Eval,
                ks[i],
                ls[j],
            ):
                return None

            Krow.append(kval)
            Erow.append(Eval)

        K.append(Krow)
        E.append(Erow)

    # --------------------------------------------------------
    # Shared p interval.
    # --------------------------------------------------------

    p_interval = intersect_intervals(
        [
            quotient_interval(
                r1[0],
                ks[0],
            ),
            quotient_interval(
                r1[1],
                ks[1],
            ),
        ]
    )

    if p_interval is None:
        return None

    # --------------------------------------------------------
    # Shared q interval.
    # --------------------------------------------------------

    q_interval = intersect_intervals(
        [
            quotient_interval(
                r2[0],
                ls[0],
            ),
            quotient_interval(
                r2[1],
                ls[1],
            ),
        ]
    )

    if q_interval is None:
        return None

    return {
        "K": K,
        "E": E,
        "p_interval": p_interval,
        "q_interval": q_interval,
    }


# ------------------------------------------------------------
# Stage 1:
#
# Anchor cell only.
# ------------------------------------------------------------

def stage1_anchors(
    n,
    r10,
    r20,
    p_est,
    q_est,
):

    Q00 = n // (
        r10 * r20
    )

    k_lo = max(
        1,
        p_est - ANCHOR_WINDOW,
    )

    k_hi = (
        p_est
        + ANCHOR_WINDOW
    )

    l_lo = max(
        1,
        q_est - ANCHOR_WINDOW,
    )

    l_hi = (
        q_est
        + ANCHOR_WINDOW
    )

    anchors = []

    for k0 in range(
        k_lo,
        k_hi + 1,
    ):

        max_l = min(
            l_hi,
            Q00 // k0,
        )

        for l0 in range(
            l_lo,
            max_l + 1,
        ):

            E0 = (
                Q00
                - k0 * l0
            )

            if not plausible_E(
                E0,
                k0,
                l0,
            ):
                continue

            anchors.append(
                (k0, l0)
            )

    return anchors


# ------------------------------------------------------------
# Stage 2:
#
# Use row 0.
#
# Given k0,l0, derive possible l1 from:
#
#     E01 = Q01-k0*l1
#
# and carry plausibility.
# ------------------------------------------------------------

def stage2_row_filter(
    n,
    anchors,
    r20,
    r21,
):

    Q00_unused = None

    out = []

    Q01 = n // (
        R1_GLOBAL[0] * r21
    )

    for k0, l0 in anchors:

        # l1 near the quotient implied by q.
        # We derive a broad range from the anchor q interval.

        qlo = (
            r20 * l0
        )

        qhi = (
            r20 * (l0 + 1) - 1
        )

        l1_min = qlo // r21
        l1_max = qhi // r21

        for l1 in range(
            max(1, l1_min),
            l1_max + 1,
        ):

            E01 = (
                Q01
                - k0 * l1
            )

            if not plausible_E(
                E01,
                k0,
                l1,
            ):
                continue

            out.append(
                (
                    k0,
                    l0,
                    l1,
                )
            )

    return out


# ------------------------------------------------------------
# Stage 3:
#
# Add column 0 and derive k1.
# ------------------------------------------------------------

def stage3_column_filter(
    n,
    candidates,
    r10,
    r11,
    r20,
):

    Q10 = n // (
        r11 * r20
    )

    out = []

    for k0, l0, l1 in candidates:

        plo = (
            r10 * k0
        )

        phi = (
            r10 * (k0 + 1) - 1
        )

        k1_min = plo // r11
        k1_max = phi // r11

        for k1 in range(
            max(1, k1_min),
            k1_max + 1,
        ):

            E10 = (
                Q10
                - k1 * l0
            )

            if not plausible_E(
                E10,
                k1,
                l0,
            ):
                continue

            out.append(
                (
                    k0,
                    k1,
                    l0,
                    l1,
                )
            )

    return out


# ------------------------------------------------------------
# Stage 4:
#
# Complete 2x2 E matrix.
# ------------------------------------------------------------

def stage4_complete(
    n,
    candidates,
    r1_grid,
    r2_grid,
):

    out = []

    Q11 = n // (
        r1_grid[1]
        * r2_grid[1]
    )

    for k0, k1, l0, l1 in candidates:

        E11 = (
            Q11
            - k1 * l1
        )

        if not plausible_E(
            E11,
            k1,
            l1,
        ):
            continue

        ks = (
            k0,
            k1,
        )

        ls = (
            l0,
            l1,
        )

        structure = (
            validate_2x2_structure(
                n,
                r1_grid,
                r2_grid,
                ks,
                ls,
            )
        )

        if structure is None:
            continue

        out.append(
            (
                ks,
                ls,
                structure,
            )
        )

    return out


# ------------------------------------------------------------
# Stage 5:
#
# Exact carry consistency over the common p/q intervals.
# ------------------------------------------------------------

def stage5_carry(
    n,
    structures,
    r1_grid,
    r2_grid,
    p_true,
    q_true,
):

    survivors = []

    true_found = False

    for ks, ls, structure in structures:

        p_lo, p_hi = (
            structure["p_interval"]
        )

        q_lo, q_hi = (
            structure["q_interval"]
        )

        for p in range(
            p_lo,
            p_hi + 1,
        ):

            for q in range(
                q_lo,
                q_hi + 1,
            ):

                ok = True

                # --------------------------------------------
                # Verify quotient vectors.
                # --------------------------------------------

                for i, r1 in enumerate(
                    r1_grid
                ):

                    if (
                        p // r1
                        != ks[i]
                    ):
                        ok = False
                        break

                if not ok:
                    continue

                for j, r2 in enumerate(
                    r2_grid
                ):

                    if (
                        q // r2
                        != ls[j]
                    ):
                        ok = False
                        break

                if not ok:
                    continue

                # --------------------------------------------
                # Exact carry matrix.
                # --------------------------------------------

                for i, r1 in enumerate(
                    r1_grid
                ):

                    for j, r2 in enumerate(
                        r2_grid
                    ):

                        E_actual = exact_E(
                            p,
                            q,
                            r1,
                            r2,
                        )

                        if (
                            E_actual
                            != structure["E"][i][j]
                        ):
                            ok = False
                            break

                    if not ok:
                        break

                if not ok:
                    continue

                exact = (
                    p * q == n
                )

                item = {
                    "p": p,
                    "q": q,
                    "k": ks,
                    "l": ls,
                    "E": structure["E"],
                    "K": structure["K"],
                    "p_interval":
                        structure["p_interval"],
                    "q_interval":
                        structure["q_interval"],
                    "exact": exact,
                }

                survivors.append(
                    item
                )

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

                    true_found = True

    return survivors, true_found


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    global R1_GLOBAL

    print("=" * 78)
    print(
        "2x2 ANCHOR COLLAPSE"
    )
    print(
        "EXPERIMENT 96"
    )
    print("=" * 78)

    # --------------------------------------------------------
    # Generate semiprime.
    # --------------------------------------------------------

    p_true, q_true, n = (
        generate_semiprime(
            TARGET
        )
    )

    r1_grid = R1_GRID
    r2_grid = R2_GRID

    R1_GLOBAL = r1_grid

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

    print()
    print(
        f"r1 = {r1_grid}"
    )

    print(
        f"r2 = {r2_grid}"
    )

    # --------------------------------------------------------
    # Natural quotient estimates.
    # --------------------------------------------------------

    p_est = (
        math.isqrt(n)
        // r1_grid[0]
    )

    q_est = (
        math.isqrt(n)
        // r2_grid[0]
    )

    print()
    print(
        f"anchor p estimate = "
        f"{p_est}"
    )

    print(
        f"anchor q estimate = "
        f"{q_est}"
    )

    # --------------------------------------------------------
    # Stage 1.
    # --------------------------------------------------------

    start = time.perf_counter()

    a1 = stage1_anchors(
        n,
        r1_grid[0],
        r2_grid[0],
        p_est,
        q_est,
    )

    t1 = (
        time.perf_counter()
        - start
    )

    print()
    print(
        "STAGE 1 — ANCHOR"
    )

    print(
        f"anchors = "
        f"{len(a1):,}"
    )

    print(
        f"time    = "
        f"{t1:.6f} s"
    )

    # --------------------------------------------------------
    # Stage 2.
    # --------------------------------------------------------

    start = time.perf_counter()

    a2 = stage2_row_filter(
        n,
        a1,
        r2_grid[0],
        r2_grid[1],
    )

    t2 = (
        time.perf_counter()
        - start
    )

    print()
    print(
        "STAGE 2 — ROW FILTER"
    )

    print(
        f"candidates = "
        f"{len(a2):,}"
    )

    print(
        f"time       = "
        f"{t2:.6f} s"
    )

    # --------------------------------------------------------
    # Stage 3.
    # --------------------------------------------------------

    start = time.perf_counter()

    a3 = stage3_column_filter(
        n,
        a2,
        r1_grid[0],
        r1_grid[1],
        r2_grid[0],
    )

    t3 = (
        time.perf_counter()
        - start
    )

    print()
    print(
        "STAGE 3 — COLUMN FILTER"
    )

    print(
        f"candidates = "
        f"{len(a3):,}"
    )

    print(
        f"time       = "
        f"{t3:.6f} s"
    )

    # --------------------------------------------------------
    # Stage 4.
    # --------------------------------------------------------

    start = time.perf_counter()

    a4 = stage4_complete(
        n,
        a3,
        r1_grid,
        r2_grid,
    )

    t4 = (
        time.perf_counter()
        - start
    )

    print()
    print(
        "STAGE 4 — COMPLETE 2x2 K/E"
    )

    print(
        f"structures = "
        f"{len(a4):,}"
    )

    print(
        f"time       = "
        f"{t4:.6f} s"
    )

    # --------------------------------------------------------
    # Stage 5.
    # --------------------------------------------------------

    start = time.perf_counter()

    survivors, true_found = (
        stage5_carry(
            n,
            a4,
            r1_grid,
            r2_grid,
            p_true,
            q_true,
        )
    )

    t5 = (
        time.perf_counter()
        - start
    )

    print()
    print(
        "STAGE 5 — EXACT CARRY"
    )

    print(
        f"survivors  = "
        f"{len(survivors):,}"
    )

    print(
        f"true found = "
        f"{true_found}"
    )

    print(
        f"time       = "
        f"{t5:.6f} s"
    )

    # --------------------------------------------------------
    # Exact survivors.
    # --------------------------------------------------------

    exact = [
        s
        for s in survivors
        if s["exact"]
    ]

    print()
    print(
        "=" * 78
    )

    print(
        "FINAL"
    )

    print(
        f"exact p*q=n = "
        f"{len(exact):,}"
    )

    for s in survivors[
        :PRINT_SOLUTIONS
    ]:

        print()
        print(
            f"p={s['p']} "
            f"q={s['q']} "
            f"exact={s['exact']}"
        )

        print(
            f"    k={s['k']}"
        )

        print(
            f"    l={s['l']}"
        )

        print(
            f"    p interval="
            f"{s['p_interval']}"
        )

        print(
            f"    q interval="
            f"{s['q_interval']}"
        )

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 96 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":

    main()


# ============================================================
# FINISHED EXPERIMENT 96
# ============================================================
