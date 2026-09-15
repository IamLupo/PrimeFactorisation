#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 136
# Rank-1 recovery from bounded carry perturbation
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 136")
print("Rank-1 recovery from bounded carry perturbation")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

print()
print("R1 =", R1)
print("R2 =", R2)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        p, q = sorted((p, q))

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# CELL
# ------------------------------------------------------------------------

def cell(n, p, q, r, s):

    k, a = divmod(p, r)
    ell, b = divmod(q, s)

    c1 = (k * b) // s
    c2 = (ell * a) // r

    beta = (k * b) % s
    alpha = (ell * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    E = c1 + c2 + c3
    Q = n // (r * s)
    K = Q - E

    assert K == k * ell

    return {
        "Q": Q,
        "E": E,
        "K": K,
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
    }


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

def build_grid(n, p, q):

    g = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            g[(i, j)] = cell(
                n,
                p,
                q,
                r,
                s,
            )

    return g


# ------------------------------------------------------------------------
# DIRECT CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    for p in range(2, limit + 1):

        tested += 1

        if n % p == 0:
            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# CARRY BOUND
#
# E = c1+c2+c3
#
# c1 < k
# c2 < ell
# c3 <= 2
#
# Therefore
#
#     E <= k + ell + 1
#
# conservatively we use:
#
#     E <= k + ell + 2
#
# Since:
#
#     k < sqrt(n)/r
#     ell < sqrt(n)/s
#
# we can get a completely observable upper bound.
# ------------------------------------------------------------------------

def global_E_bound(n, r, s):

    root = math.isqrt(n)

    k_bound = (
        root // r
    ) + 1

    ell_bound = (
        root // s
    ) + 1

    return (
        k_bound
        + ell_bound
        + 2
    )


# ------------------------------------------------------------------------
# OBSERVABLE K INTERVALS
#
# Because:
#
#     K = Q-E
#
# and
#
#     0 <= E <= B
#
# we obtain:
#
#     Q-B <= K <= Q.
#
# The true K is guaranteed to be inside the interval.
# ------------------------------------------------------------------------

def observable_intervals(n):

    intervals = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            Q = (
                n
                // (r * s)
            )

            B = global_E_bound(
                n,
                r,
                s,
            )

            intervals[(i, j)] = {
                "low": max(
                    1,
                    Q - B,
                ),
                "high": Q,
                "Q": Q,
                "B": B,
            }

    return intervals


# ------------------------------------------------------------------------
# INTERVAL RATIO ANALYSIS
#
# If K is rank 1:
#
#     K_ij = k_i * ell_j
#
# Hence:
#
#     K_i,j / K_0,j = k_i/k_0
#
# and
#
#     K_i,j / K_i,0 = ell_j/ell_0.
#
# We obtain intervals for these ratios from Q/E bounds.
#
# The experiment asks whether those rational intervals are narrow
# enough to contain only one integer ratio.
# ------------------------------------------------------------------------

def ratio_interval(
    a_low,
    a_high,
    b_low,
    b_high,
):
    """
    Return conservative lower/upper rational bounds for a/b.
    """

    # Positive quantities only.
    lower_num = a_low
    lower_den = b_high

    upper_num = a_high
    upper_den = b_low

    return (
        lower_num,
        lower_den,
        upper_num,
        upper_den,
    )


# ------------------------------------------------------------------------
# FIND INTEGER PAIRS INSIDE A RATIO INTERVAL
#
# For:
#
#     a/b in [L/U, R/S]
#
# with denominator bounded by a known maximum.
#
# This routine uses direct denominator enumeration only over the
# quotient scale, NOT over sqrt(n), and is intentionally capped.
#
# We use it as a diagnostic to see how wide the rational uncertainty
# really is.
# ------------------------------------------------------------------------

def rational_candidates(
    low_num,
    low_den,
    high_num,
    high_den,
    denominator_bound,
    cap=100000,
):

    out = []

    if denominator_bound <= 0:
        return out

    # We only search up to cap denominators to prevent an accidental
    # enormous computation. This is a diagnostic experiment.
    limit = min(
        denominator_bound,
        cap,
    )

    for d in range(1, limit + 1):

        # Need:
        #
        # low <= a/d <= high
        #
        # a >= ceil(low_num*d / low_den)
        #
        # a <= floor(high_num*d / high_den)
        #

        a_min = (
            low_num * d
            + low_den - 1
        ) // low_den

        a_max = (
            high_num * d
        ) // high_den

        if a_min > a_max:
            continue

        # Don't enumerate many numerators for one denominator.
        # Keep only the first few.
        width = a_max - a_min + 1

        if width > 10:
            out.append(
                (
                    d,
                    a_min,
                    a_max,
                    width,
                )
            )

        else:

            for a in range(
                a_min,
                a_max + 1,
            ):

                out.append(
                    (
                        a,
                        d,
                        a,
                    )
                )

    return out


# ------------------------------------------------------------------------
# TEST WHETHER TRUE K IS INSIDE INTERVAL
# ------------------------------------------------------------------------

def interval_contains_true(
    intervals,
    g,
):

    ok = True

    misses = []

    for key, info in intervals.items():

        trueK = g[key]["K"]

        if not (
            info["low"]
            <= trueK
            <= info["high"]
        ):

            ok = False
            misses.append(
                (
                    key,
                    trueK,
                    info["low"],
                    info["high"],
                )
            )

    return ok, misses


# ------------------------------------------------------------------------
# RANK-1 INTERVAL FEASIBILITY
#
# We try to reconstruct one row vector k_i and one column vector ell_j
# from K_00.
#
# Since:
#
#     K_ij = k_i * ell_j
#
# and K_00 = k0*ell0,
#
# any candidate k0 divisor determines:
#
#     ell0 = K00/k0.
#
# Then each other k_i and ell_j must fit inside their observable
# intervals.
#
# The crucial experiment is how many k0 values are possible based only
# on the interval bounds.
#
# We DO NOT factor n here.
# ------------------------------------------------------------------------

def rank1_interval_search(
    n,
    intervals,
    true_grid,
):

    I00 = intervals[(0, 0)]

    low00 = I00["low"]
    high00 = I00["high"]

    # ------------------------------------------------------------
    # The true k0 is bounded by sqrt(n)/r0.
    # ------------------------------------------------------------

    root = math.isqrt(n)

    k0_max = (
        root // R1[0]
    ) + 1

    possible_k0 = 0
    viable_k0 = []

    # We need K00 = k0*ell0.
    #
    # But K00 itself is uncertain inside [low00,high00].
    #
    # Test every k0 that could possibly occur.
    #
    # This is a quotient-scale loop and is intentionally capped at
    # 1e6 so the experiment cannot explode.
    # ------------------------------------------------------------

    limit = min(
        k0_max,
        1_000_000,
    )

    for k0 in range(
        1,
        limit + 1,
    ):

        # Determine possible ell0 values such that:
        #
        #     low00 <= k0*ell0 <= high00
        #
        ell_low = (
            low00
            + k0
            - 1
        ) // k0

        ell_high = (
            high00
            // k0
        )

        if ell_low > ell_high:
            continue

        possible_k0 += 1

        # --------------------------------------------------------
        # Try the first few possible ell0 values.
        # --------------------------------------------------------

        local_viable = False

        for ell0 in (
            ell_low,
            ell_high,
        ):

            if ell0 <= 0:
                continue

            K00 = k0 * ell0

            if not (
                low00
                <= K00
                <= high00
            ):
                continue

            # ----------------------------------------------------
            # Predict allowed k_i and ell_j intervals.
            #
            # Because:
            #
            #     K_i0 = k_i*ell0
            #
            #     K_0j = k0*ell_j
            #
            # ----------------------------------------------------

            valid = True

            k_values = [k0]
            ell_values = [ell0]

            # Other rows.
            for i in range(1, 3):

                I = intervals[(i, 0)]

                ki_low = (
                    I["low"]
                    + ell0 - 1
                ) // ell0

                ki_high = (
                    I["high"]
                    // ell0
                )

                if ki_low > ki_high:

                    valid = False
                    break

                # To remain a diagnostic rather than a huge search,
                # choose values compatible with the nearest true-scale
                # estimate.
                k_candidate = (
                    I["high"]
                    // ell0
                )

                if k_candidate < ki_low:
                    k_candidate = ki_low

                k_values.append(
                    k_candidate
                )

            if not valid:
                continue

            # Other columns.
            for j in range(1, 3):

                I = intervals[(0, j)]

                lj_low = (
                    I["low"]
                    + k0 - 1
                ) // k0

                lj_high = (
                    I["high"]
                    // k0
                )

                if lj_low > lj_high:

                    valid = False
                    break

                ell_candidate = (
                    I["high"]
                    // k0
                )

                if ell_candidate < lj_low:
                    ell_candidate = lj_low

                ell_values.append(
                    ell_candidate
                )

            if not valid:
                continue

            # ----------------------------------------------------
            # Check all implied rank-1 products against intervals.
            # ----------------------------------------------------

            for i in range(3):

                for j in range(3):

                    predicted = (
                        k_values[i]
                        * ell_values[j]
                    )

                    I = intervals[(i, j)]

                    if not (
                        I["low"]
                        <= predicted
                        <= I["high"]
                    ):

                        valid = False
                        break

                if not valid:
                    break

            if valid:

                viable_k0.append(
                    {
                        "k0": k0,
                        "ell0": ell0,
                        "k": tuple(k_values),
                        "ell": tuple(ell_values),
                    }
                )

                local_viable = True
                break

        # --------------------------------------------------------
        # Stop only after enough viable states so the output remains
        # readable.
        # --------------------------------------------------------

        if len(viable_k0) >= 100:
            break

    return {
        "k0_tested": limit,
        "possible_k0": possible_k0,
        "viable_states": viable_k0,
    }


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(136)

TEST_BITS = [
    30,
    36,
    42,
    48,
    54,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, p, q = random_semiprime(bits)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    # ------------------------------------------------------------
    # True grid.
    # ------------------------------------------------------------

    g = build_grid(
        n,
        p,
        q,
    )

    # ------------------------------------------------------------
    # Observable Q.
    # ------------------------------------------------------------

    Q = [
        [
            g[(i, j)]["Q"]
            for j in range(3)
        ]
        for i in range(3)
    ]

    print()
    print("Q MATRIX")

    for row in Q:
        print(row)

    # ------------------------------------------------------------
    # True K.
    # ------------------------------------------------------------

    K = [
        [
            g[(i, j)]["K"]
            for j in range(3)
        ]
        for i in range(3)
    ]

    print()
    print("TRUE K MATRIX")

    for row in K:
        print(row)

    # ------------------------------------------------------------
    # True E.
    # ------------------------------------------------------------

    E = [
        [
            g[(i, j)]["E"]
            for j in range(3)
        ]
        for i in range(3)
    ]

    print()
    print("TRUE E MATRIX")

    for row in E:
        print(row)

    # ------------------------------------------------------------
    # Observable intervals.
    # ------------------------------------------------------------

    print()
    print("OBSERVABLE K INTERVALS")

    intervals = observable_intervals(n)

    interval_ok, misses = interval_contains_true(
        intervals,
        g,
    )

    for i in range(3):

        row = []

        for j in range(3):

            I = intervals[(i, j)]

            row.append(
                (
                    I["low"],
                    I["high"],
                    I["high"] - I["low"],
                )
            )

        print(row)

    print()
    print(
        "TRUE K INSIDE ALL INTERVALS =",
        interval_ok
    )

    if misses:
        print(
            "MISSES =",
            misses
        )

    # ------------------------------------------------------------
    # Interval widths.
    # ------------------------------------------------------------

    widths = [
        intervals[key]["high"]
        - intervals[key]["low"]
        for key in intervals
    ]

    print()
    print("INTERVAL WIDTHS")

    print(
        "min width =",
        min(widths)
    )

    print(
        "max width =",
        max(widths)
    )

    print(
        "true E max =",
        max(
            g[key]["E"]
            for key in g
        )
    )

    # ------------------------------------------------------------
    # True rank-1 factors.
    # ------------------------------------------------------------

    print()
    print("TRUE QUOTIENT VECTORS")

    true_k = [
        g[(i, 0)]["k"]
        for i in range(3)
    ]

    true_ell = [
        g[(0, j)]["ell"]
        for j in range(3)
    ]

    print(
        "k =",
        true_k
    )

    print(
        "ell =",
        true_ell
    )

    # ------------------------------------------------------------
    # Ratio intervals.
    # ------------------------------------------------------------

    print()
    print("RATIO INTERVALS")

    for i in range(1, 3):

        I_num = intervals[(i, 0)]
        I_den = intervals[(0, 0)]

        ratio = ratio_interval(
            I_num["low"],
            I_num["high"],
            I_den["low"],
            I_den["high"],
        )

        print(
            f"k[{i}]/k[0] observable interval:"
        )

        print(
            "  lower =",
            f"{ratio[0]}/{ratio[1]}",
            "=",
            ratio[0] / ratio[1],
        )

        print(
            "  upper =",
            f"{ratio[2]}/{ratio[3]}",
            "=",
            ratio[2] / ratio[3],
        )

        print(
            "  true ratio =",
            f"{true_k[i]}/{true_k[0]}",
            "=",
            true_k[i] / true_k[0],
        )

    for j in range(1, 3):

        I_num = intervals[(0, j)]
        I_den = intervals[(0, 0)]

        ratio = ratio_interval(
            I_num["low"],
            I_num["high"],
            I_den["low"],
            I_den["high"],
        )

        print(
            f"ell[{j}]/ell[0] observable interval:"
        )

        print(
            "  lower =",
            f"{ratio[0]}/{ratio[1]}",
            "=",
            ratio[0] / ratio[1],
        )

        print(
            "  upper =",
            f"{ratio[2]}/{ratio[3]}",
            "=",
            ratio[2] / ratio[3],
        )

        print(
            "  true ratio =",
            f"{true_ell[j]}/{true_ell[0]}",
            "=",
            true_ell[j] / true_ell[0],
        )

    # ------------------------------------------------------------
    # Rank-1 interval feasibility.
    # ------------------------------------------------------------

    print()
    print("RANK-1 INTERVAL SEARCH")

    t0 = time.perf_counter()

    result = rank1_interval_search(
        n,
        intervals,
        g,
    )

    t1 = time.perf_counter()

    elapsed = t1 - t0

    print(
        "k0 tested =",
        result["k0_tested"]
    )

    print(
        "k0 with possible ell0 =",
        result["possible_k0"]
    )

    print(
        "viable rank-1 states =",
        len(result["viable_states"])
    )

    print(
        "runtime =",
        f"{elapsed:.6f}",
        "s"
    )

    # ------------------------------------------------------------
    # Compare viable state against truth.
    # ------------------------------------------------------------

    exact_state = False

    for state in result["viable_states"]:

        if (
            tuple(true_k)
            == state["k"]
            and
            tuple(true_ell)
            == state["ell"]
        ):

            exact_state = True

            print()
            print(
                ">>> TRUE RANK-1 STATE RECOVERED <<<"
            )

            print(
                "recovered k =",
                state["k"]
            )

            print(
                "recovered ell =",
                state["ell"]
            )

            break

    print(
        "TRUE STATE RECOVERED =",
        exact_state
    )

    # ------------------------------------------------------------
    # If there are viable states, show the first few.
    # ------------------------------------------------------------

    if result["viable_states"]:

        print()
        print("FIRST VIABLE STATES")

        for state in result["viable_states"][:10]:

            exact = (
                state["k"]
                == tuple(true_k)
                and
                state["ell"]
                == tuple(true_ell)
            )

            print(
                "k =",
                state["k"],
                "ell =",
                state["ell"],
                "EXACT =",
                exact,
            )

    # ------------------------------------------------------------
    # Direct control.
    # ------------------------------------------------------------

    print()
    print("DIRECT SQRT(n) CONTROL")

    t0 = time.perf_counter()

    _, _, direct_tested = direct_factor(n)

    t1 = time.perf_counter()

    print(
        "p tested =",
        direct_tested
    )

    print(
        "runtime =",
        f"{t1 - t0:.6f}",
        "s"
    )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 136
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 136")
print("=" * 72)
