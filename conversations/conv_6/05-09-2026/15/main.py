import math
import time
import sympy


# ============================================================
# START EXPERIMENT 98
# ============================================================
#
# 1D ANCHOR REDUCTION USING THE E <= k+l BOUND
#
#
# At the anchor:
#
#     Q = floor(n / (r1*r2))
#
#     Q = k*l + E
#
# and:
#
#     E = c1+c2+c3
#
# with:
#
#     c1 <= k-1
#     c2 <= l-1
#     c3 <= 2
#
# therefore:
#
#     E <= k+l
#
# Thus:
#
#     0 <= Q-k*l <= k+l
#
# which gives:
#
#     l <= floor(Q/k)
#
# and:
#
#     l >= ceil((Q-k-C)/(k+1))
#
# for a small conservative C.
#
#
# Experiment 96/97 searched a large 2D (k,l) rectangle.
#
# Experiment 98 searches only k and DERIVES a tiny l interval.
#
#
# Then for each surviving (k,l):
#
#     p = r10*k + a
#     q = r20*l + b
#
# with small a,b.
#
# Every p,q is tested against the full 2x2 carry structure.
#
#
# NO p*q == n is used as a search filter.
#
# It is only reported at the end.
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

R1 = [3, 11]
R2 = [5, 17]

# Search range for k around sqrt(n)/r1.
#
# This remains one-dimensional.
ANCHOR_WINDOW = 3000

# Conservative constant in:
#
#     E <= k+l+C
#
E_C = 10

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
# Correct ceil division
# ------------------------------------------------------------

def ceil_div(a, b):

    return -((-a) // b)


# ------------------------------------------------------------
# Exact carry calculation
# ------------------------------------------------------------

def compute_state(
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

    E = (
        c1
        + c2
        + c3
    )

    K = k * l

    Q = None

    return {
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
# Full 2x2 carry test
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

            state = compute_state(
                p,
                q,
                r1,
                r2,
            )

            K = state["K"]
            E = state["E"]

            Q = n // (
                r1 * r2
            )

            # ------------------------------------------------
            # Core K/E identity.
            # ------------------------------------------------

            if K + E != Q:
                return None

            row.append(E)

        E_matrix.append(row)

    return E_matrix


# ------------------------------------------------------------
# Derive l interval for fixed k
#
# We know:
#
#     0 <= Q-k*l
#
# and:
#
#     Q-k*l <= k+l+C
#
# Therefore:
#
#     l <= Q/k
#
#     l >= (Q-k-C)/(k+1)
# ------------------------------------------------------------

def derive_l_interval(
    Q,
    k,
):

    # Upper:
    #
    #     k*l <= Q
    #
    l_max = (
        Q // k
    )

    # Lower:
    #
    #     Q-k*l <= k+l+C
    #
    #     Q-k-C <= l*(k+1)
    #
    #     l >= (Q-k-C)/(k+1)
    #

    l_min = ceil_div(
        Q - k - E_C,
        k + 1,
    )

    return (
        max(1, l_min),
        max(0, l_max),
    )


# ------------------------------------------------------------
# Direct 1D anchor search
# ------------------------------------------------------------

def search(
    n,
    p_true,
    q_true,
    r1_grid,
    r2_grid,
):

    r10 = r1_grid[0]
    r20 = r2_grid[0]

    Q00 = n // (
        r10 * r20
    )

    # --------------------------------------------------------
    # Estimate k from sqrt(n)/r1.
    # --------------------------------------------------------

    k_center = (
        math.isqrt(n)
        // r10
    )

    k_lo = max(
        1,
        k_center - ANCHOR_WINDOW,
    )

    k_hi = (
        k_center
        + ANCHOR_WINDOW
    )

    # --------------------------------------------------------
    # Statistics.
    # --------------------------------------------------------

    k_tests = 0
    l_tests = 0
    pq_tests = 0
    carry_tests = 0

    anchors = 0

    survivors = []

    seen = set()

    start = time.perf_counter()

    # --------------------------------------------------------
    # ONE-dimensional k search.
    # --------------------------------------------------------

    for k in range(
        k_lo,
        k_hi + 1,
    ):

        k_tests += 1

        l_min, l_max = derive_l_interval(
            Q00,
            k,
        )

        # ----------------------------------------------------
        # Intersect with a broad natural q estimate.
        #
        # This is NOT assuming k ~= l.
        # It is only limiting the anchor to the chosen
        # experimental search region.
        # ----------------------------------------------------

        l_center = (
            math.isqrt(n)
            // r20
        )

        l_search_min = max(
            1,
            l_center - ANCHOR_WINDOW,
        )

        l_search_max = (
            l_center
            + ANCHOR_WINDOW
        )

        l_min = max(
            l_min,
            l_search_min,
        )

        l_max = min(
            l_max,
            l_search_max,
        )

        if l_min > l_max:
            continue

        # ----------------------------------------------------
        # Usually this range should be very small.
        # ----------------------------------------------------

        for l in range(
            l_min,
            l_max + 1,
        ):

            l_tests += 1

            E00 = (
                Q00
                - k * l
            )

            if E00 < 0:
                continue

            if E00 > (
                k
                + l
                + E_C
            ):
                continue

            anchors += 1

            # ------------------------------------------------
            # Anchor p interval.
            # ------------------------------------------------

            p_lo = (
                r10 * k
            )

            p_hi = (
                r10 * (k + 1)
                - 1
            )

            # ------------------------------------------------
            # Anchor q interval.
            # ------------------------------------------------

            q_lo = (
                r20 * l
            )

            q_hi = (
                r20 * (l + 1)
                - 1
            )

            # ------------------------------------------------
            # Search the tiny rectangle.
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
                    # Full 2x2 structure.
                    # ------------------------------------------------

                    E_matrix = []

                    valid = True

                    for r1 in r1_grid:

                        E_row = []

                        for r2 in r2_grid:

                            state = compute_state(
                                p,
                                q,
                                r1,
                                r2,
                            )

                            Q = (
                                n
                                // (r1 * r2)
                            )

                            if (
                                state["K"]
                                + state["E"]
                                != Q
                            ):

                                valid = False
                                break

                            E_row.append(
                                state["E"]
                            )

                            carry_tests += 1

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

    exact = [
        s
        for s in survivors
        if s["exact"]
    ]

    return {
        "k_tests": k_tests,
        "l_tests": l_tests,
        "anchors": anchors,
        "pq_tests": pq_tests,
        "carry_tests": carry_tests,
        "survivors": survivors,
        "true_found": true_found,
        "exact": exact,
        "time": elapsed,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "1D ANCHOR REDUCTION USING E <= k+l"
    )
    print(
        "EXPERIMENT 98"
    )
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

    print()
    print(
        f"r1 = {R1}"
    )

    print(
        f"r2 = {R2}"
    )

    Q00 = n // (
        R1[0] * R2[0]
    )

    print()
    print(
        f"Q00 = {Q00}"
    )

    # --------------------------------------------------------
    # Show the TRUE anchor.
    # --------------------------------------------------------

    true_k = p_true // R1[0]
    true_l = q_true // R2[0]

    true_E = (
        Q00
        - true_k * true_l
    )

    true_l_min, true_l_max = (
        derive_l_interval(
            Q00,
            true_k,
        )
    )

    print()
    print(
        "TRUE ANCHOR"
    )

    print(
        f"    k = {true_k}"
    )

    print(
        f"    l = {true_l}"
    )

    print(
        f"    E = {true_E}"
    )

    print(
        f"    derived l interval = "
        f"{true_l_min} .. {true_l_max}"
    )

    print(
        f"    true l inside = "
        f"{true_l_min <= true_l <= true_l_max}"
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    result = search(
        n,
        p_true,
        q_true,
        R1,
        R2,
    )

    # --------------------------------------------------------
    # Results.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("RESULT")
    print("=" * 78)

    print(
        f"k values tested   = "
        f"{result['k_tests']:,}"
    )

    print(
        f"l values tested   = "
        f"{result['l_tests']:,}"
    )

    print(
        f"anchors           = "
        f"{result['anchors']:,}"
    )

    print(
        f"p/q tests         = "
        f"{result['pq_tests']:,}"
    )

    print(
        f"carry tests       = "
        f"{result['carry_tests']:,}"
    )

    print(
        f"survivors         = "
        f"{len(result['survivors']):,}"
    )

    print(
        f"TRUE found        = "
        f"{result['true_found']}"
    )

    print(
        f"exact p*q==n      = "
        f"{len(result['exact']):,}"
    )

    print(
        f"time              = "
        f"{result['time']:.6f} s"
    )

    # --------------------------------------------------------
    # Survivor details.
    # --------------------------------------------------------

    print()
    print(
        "SURVIVORS"
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

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 98 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":

    main()


# ============================================================
# FINISHED EXPERIMENT 98
# ============================================================
