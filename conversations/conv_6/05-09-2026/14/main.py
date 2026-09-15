import math
import time
import sympy


# ============================================================
# START EXPERIMENT 97
# ============================================================
#
# DIRECT 2x2 ANCHOR -> p,q -> FULL CARRY TEST
#
# Experiment 96 used:
#
#     anchor
#       ->
#     row options
#       ->
#     column options
#       ->
#     K/E structures
#       ->
#     p/q intervals
#       ->
#     carry test
#
# This experiment removes the middle stages.
#
# For an anchor:
#
#     k0 = floor(p/r10)
#     l0 = floor(q/r20)
#
# we immediately know:
#
#     r10*k0 <= p < r10*(k0+1)
#
#     r20*l0 <= q < r20*(l0+1)
#
# Therefore the anchor defines a SMALL rectangle of possible
# p,q values.
#
# We simply test every p,q in that rectangle against ALL
# four cells of the 2x2 grid.
#
#
# No p*q == n is used as a filter.
#
# We only use:
#
#     Qij = Kij + Eij
#
# and the exact carry equations.
#
#
# The purpose is to see whether we can replace the expensive
# 4.7 second anchor pipeline with one direct pass.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

R1 = [3, 11]
R2 = [5, 17]

ANCHOR_WINDOW = 3000

PRINT_SURVIVORS = 20


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
# Exact carry E
# ------------------------------------------------------------

def compute_E(
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
# Complete 2x2 test
#
# For a candidate p,q:
#
#     Kij = ki*lj
#     Eij = Qij-Kij
#
# and exact carry E must equal that value.
# ------------------------------------------------------------

def full_grid_test(
    n,
    p,
    q,
    r1_grid,
    r2_grid,
):

    E_matrix = []

    for r1 in r1_grid:

        row = []

        for r2 in r2_grid:

            k, a = divmod(
                p,
                r1,
            )

            l, b = divmod(
                q,
                r2,
            )

            K = k * l

            Q = (
                n
                // (r1 * r2)
            )

            E_expected = Q - K

            if E_expected < 0:
                return None

            E_actual = compute_E(
                p,
                q,
                r1,
                r2,
            )

            if E_actual != E_expected:
                return None

            row.append(
                E_actual
            )

        E_matrix.append(row)

    return E_matrix


# ------------------------------------------------------------
# Anchor E bound
#
# We know exactly:
#
#     c1 <= k-1
#     c2 <= l-1
#     c3 <= 2
#
# hence:
#
#     E <= k+l
#
# ------------------------------------------------------------

def anchor_is_possible(
    n,
    r1,
    r2,
    k,
    l,
):

    Q = n // (r1 * r2)

    E = Q - k * l

    if E < 0:
        return False

    return E <= k + l


# ------------------------------------------------------------
# Direct anchor search
# ------------------------------------------------------------

def direct_anchor_search(
    n,
    p_true,
    q_true,
    r1_grid,
    r2_grid,
):

    r10 = r1_grid[0]
    r20 = r2_grid[0]

    # Natural quotient estimates.
    p_center = math.isqrt(n) // r10
    q_center = math.isqrt(n) // r20

    k_min = max(
        1,
        p_center - ANCHOR_WINDOW,
    )

    k_max = (
        p_center
        + ANCHOR_WINDOW
    )

    l_min = max(
        1,
        q_center - ANCHOR_WINDOW,
    )

    l_max = (
        q_center
        + ANCHOR_WINDOW
    )

    anchors = 0
    pq_tests = 0
    carry_tests = 0

    survivors = []

    seen = set()

    start = time.perf_counter()

    # --------------------------------------------------------
    # Anchor enumeration.
    # --------------------------------------------------------

    for k0 in range(
        k_min,
        k_max + 1,
    ):

        # From:
        #
        #     k0*l0 <= Q00
        #
        # derive an upper l bound.
        #
        Q00 = n // (
            r10 * r20
        )

        l_upper_from_Q = (
            Q00 // k0
        )

        l_hi = min(
            l_max,
            l_upper_from_Q,
        )

        for l0 in range(
            l_min,
            l_hi + 1,
        ):

            if not anchor_is_possible(
                n,
                r10,
                r20,
                k0,
                l0,
            ):
                continue

            anchors += 1

            # ------------------------------------------------
            # Anchor p interval.
            # ------------------------------------------------

            p_lo = r10 * k0
            p_hi = (
                r10 * (k0 + 1)
                - 1
            )

            # ------------------------------------------------
            # Anchor q interval.
            # ------------------------------------------------

            q_lo = r20 * l0
            q_hi = (
                r20 * (l0 + 1)
                - 1
            )

            # ------------------------------------------------
            # Enumerate only the tiny anchor rectangle.
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

                    # ------------------------------------------------
                    # Verify that this p/q actually has the anchor
                    # quotient.
                    # ------------------------------------------------

                    if p // r10 != k0:
                        continue

                    if q // r20 != l0:
                        continue

                    # ------------------------------------------------
                    # Full 2x2 carry consistency.
                    # ------------------------------------------------

                    E_matrix = []

                    valid = True

                    for r1 in r1_grid:

                        E_row = []

                        for r2 in r2_grid:

                            k, a = divmod(
                                p,
                                r1,
                            )

                            l, b = divmod(
                                q,
                                r2,
                            )

                            K = k * l

                            Q = n // (
                                r1 * r2
                            )

                            E_expected = (
                                Q - K
                            )

                            if E_expected < 0:
                                valid = False
                                break

                            # Exact carry.
                            E_actual = compute_E(
                                p,
                                q,
                                r1,
                                r2,
                            )

                            carry_tests += 1

                            if (
                                E_actual
                                != E_expected
                            ):
                                valid = False
                                break

                            E_row.append(
                                E_actual
                            )

                        if not valid:
                            break

                        E_matrix.append(
                            E_row
                        )

                    if not valid:
                        continue

                    key = (
                        p,
                        q,
                    )

                    if key in seen:
                        continue

                    seen.add(key)

                    survivors.append({
                        "p": p,
                        "q": q,
                        "E": E_matrix,
                        "exact": (
                            p * q == n
                        ),
                    })

    elapsed = (
        time.perf_counter()
        - start
    )

    true_found = any(
        (
            s["p"] == p_true
            and s["q"] == q_true
        )
        or
        (
            s["p"] == q_true
            and s["q"] == p_true
        )
        for s in survivors
    )

    exact_count = sum(
        1
        for s in survivors
        if s["exact"]
    )

    return {
        "anchors": anchors,
        "pq_tests": pq_tests,
        "carry_tests": carry_tests,
        "survivors": survivors,
        "true_found": true_found,
        "exact_count": exact_count,
        "time": elapsed,
    }


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "DIRECT 2x2 ANCHOR -> p/q -> CARRY TEST"
    )
    print(
        "EXPERIMENT 97"
    )
    print("=" * 78)

    # --------------------------------------------------------
    # Test number.
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

    print()
    print(
        f"r1 = {R1}"
    )

    print(
        f"r2 = {R2}"
    )

    print(
        f"grid cells = "
        f"{len(R1)*len(R2)}"
    )

    # --------------------------------------------------------
    # Run.
    # --------------------------------------------------------

    result = direct_anchor_search(
        n,
        p_true,
        q_true,
        R1,
        R2,
    )

    print()
    print("=" * 78)
    print("RESULT")
    print("=" * 78)

    print(
        f"anchors              = "
        f"{result['anchors']:,}"
    )

    print(
        f"p/q tests            = "
        f"{result['pq_tests']:,}"
    )

    print(
        f"carry cell tests     = "
        f"{result['carry_tests']:,}"
    )

    print(
        f"survivors            = "
        f"{len(result['survivors']):,}"
    )

    print(
        f"TRUE found           = "
        f"{result['true_found']}"
    )

    print(
        f"exact p*q==n         = "
        f"{result['exact_count']:,}"
    )

    print(
        f"time                 = "
        f"{result['time']:.6f} s"
    )

    # --------------------------------------------------------
    # Survivor examples.
    # --------------------------------------------------------

    print()
    print(
        "SURVIVORS"
    )

    if not result["survivors"]:

        print(
            "    none"
        )

    else:

        for s in result[
            "survivors"
        ][:PRINT_SURVIVORS]:

            print(
                "    "
                f"p={s['p']} "
                f"q={s['q']} "
                f"exact={s['exact']}"
            )

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 97 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":

    main()


# ============================================================
# FINISHED EXPERIMENT 97
# ============================================================
