import math
import time
import sympy


# ============================================================
# START EXPERIMENT 99
# ============================================================
#
# 2D QUOTIENT-TRANSITION FINGERPRINT
#
# Goal:
#
#   Test whether a WIDELY SPACED 2D r1/r2 grid can reduce
#   the k,l search before touching p,q or a,b.
#
#
# For one anchor pair:
#
#     r10, r20
#
#     Q00 = floor(n/(r10*r20))
#
#     Q00 = k0*l0 + E00
#
#
# We use:
#
#     E00 <= k0+l0+C
#
# to derive a small l0 interval for each k0.
#
#
# Then a SECOND, MUCH LARGER r1 gives:
#
#     p in [r10*k0, r10*(k0+1)-1]
#
# so possible:
#
#     k1 = floor(p/r11)
#
# can be derived WITHOUT knowing p.
#
# Likewise a larger r2 gives possible l1.
#
#
# We then require:
#
#     E10 = Q10-k1*l0
#
#     E01 = Q01-k0*l1
#
#     E11 = Q11-k1*l1
#
# to satisfy the carry-scale bounds.
#
#
# NO a,b search.
# NO p,q search.
# NO p*q==n filter.
#
# At the end we only classify whether the TRUE quotient
# structure survived.
#
#
# Main question:
#
#     Can widely separated r values turn the enormous
#     k/l search into a much smaller transition search?
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

TARGET = 10**9

# Anchor values.
R1_SMALL = 3
R2_SMALL = 5

# Widely separated second level.
#
# Change these to test other scales.
R1_LARGE = 101
R2_LARGE = 103

# Search width for k0.
#
# This is deliberately only one-dimensional.
ANCHOR_K_WINDOW = 5000

# Conservative carry constant.
E_C = 10

# Print a few surviving structures.
PRINT_SURVIVORS = 20


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
# Correct ceil division
# ------------------------------------------------------------

def ceil_div(a, b):

    return -((-a) // b)


# ------------------------------------------------------------
# Carry upper bound
#
# Conservative:
#
#     E <= k+l+C
#
# ------------------------------------------------------------

def plausible_E(E, k, l):

    if E < 0:
        return False

    return E <= k + l + E_C


# ------------------------------------------------------------
# l interval from:
#
#     0 <= Q-k*l <= k+l+C
#
# ------------------------------------------------------------

def derive_l_interval(
    Q,
    k,
):

    # From:
    #
    #     k*l <= Q
    #
    l_max = Q // k

    # From:
    #
    #     Q-k*l <= k+l+C
    #
    #     Q-k-C <= l*(k+1)
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
# Quotient transition interval.
#
# If:
#
#     x in [r_small*k, r_small*(k+1)-1]
#
# then:
#
#     floor(x/r_large)
#
# lies in this interval.
# ------------------------------------------------------------

def quotient_transition_range(
    r_small,
    k_small,
    r_large,
):

    x_lo = (
        r_small
        * k_small
    )

    x_hi = (
        r_small
        * (k_small + 1)
        - 1
    )

    q_lo = x_lo // r_large
    q_hi = x_hi // r_large

    return (
        q_lo,
        q_hi,
    )


# ------------------------------------------------------------
# Build Q values
# ------------------------------------------------------------

def build_Q(
    n,
    r1,
    r2,
):

    return n // (
        r1 * r2
    )


# ------------------------------------------------------------
# True quotient structure
# ------------------------------------------------------------

def true_structure(
    p,
    q,
):

    return {
        "k0": p // R1_SMALL,
        "k1": p // R1_LARGE,
        "l0": q // R2_SMALL,
        "l1": q // R2_LARGE,
    }


# ------------------------------------------------------------
# Search transition structures
# ------------------------------------------------------------

def search(
    n,
    p_true,
    q_true,
):

    r10 = R1_SMALL
    r20 = R2_SMALL

    r11 = R1_LARGE
    r21 = R2_LARGE

    Q00 = build_Q(
        n,
        r10,
        r20,
    )

    Q10 = build_Q(
        n,
        r11,
        r20,
    )

    Q01 = build_Q(
        n,
        r10,
        r21,
    )

    Q11 = build_Q(
        n,
        r11,
        r21,
    )

    # --------------------------------------------------------
    # Search k0 only.
    # --------------------------------------------------------

    k_center = (
        math.isqrt(n)
        // r10
    )

    k_min = max(
        1,
        k_center - ANCHOR_K_WINDOW,
    )

    k_max = (
        k_center
        + ANCHOR_K_WINDOW
    )

    k_tests = 0
    l_tests = 0
    transition_tests = 0
    final_structures = 0

    survivors = []

    seen = set()

    start = time.perf_counter()

    # --------------------------------------------------------
    # MAIN 1D LOOP
    # --------------------------------------------------------

    for k0 in range(
        k_min,
        k_max + 1,
    ):

        k_tests += 1

        # ----------------------------------------------------
        # Derive l0 interval from Q00.
        # ----------------------------------------------------

        l0_min, l0_max = (
            derive_l_interval(
                Q00,
                k0,
            )
        )

        if l0_min > l0_max:
            continue

        # ----------------------------------------------------
        # For every possible l0.
        # ----------------------------------------------------

        for l0 in range(
            l0_min,
            l0_max + 1,
        ):

            l_tests += 1

            E00 = (
                Q00
                - k0 * l0
            )

            if not plausible_E(
                E00,
                k0,
                l0,
            ):
                continue

            # ------------------------------------------------
            # Transition k0 -> possible k1.
            # ------------------------------------------------

            k1_min, k1_max = (
                quotient_transition_range(
                    r10,
                    k0,
                    r11,
                )
            )

            # ------------------------------------------------
            # Transition l0 -> possible l1.
            # ------------------------------------------------

            l1_min, l1_max = (
                quotient_transition_range(
                    r20,
                    l0,
                    r21,
                )
            )

            # Discard zero quotients.
            k1_min = max(
                1,
                k1_min,
            )

            l1_min = max(
                1,
                l1_min,
            )

            if (
                k1_min > k1_max
                or
                l1_min > l1_max
            ):
                continue

            # ------------------------------------------------
            # Combine transitions.
            #
            # With widely separated r's these ranges should
            # often have width 0 or 1.
            # ------------------------------------------------

            for k1 in range(
                k1_min,
                k1_max + 1,
            ):

                for l1 in range(
                    l1_min,
                    l1_max + 1,
                ):

                    transition_tests += 1

                    # ----------------------------------------
                    # E10
                    # ----------------------------------------

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

                    # ----------------------------------------
                    # E01
                    # ----------------------------------------

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

                    # ----------------------------------------
                    # E11
                    # ----------------------------------------

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

                    final_structures += 1

                    key = (
                        k0,
                        l0,
                        k1,
                        l1,
                    )

                    if key in seen:
                        continue

                    seen.add(key)

                    survivors.append({
                        "k0": k0,
                        "l0": l0,
                        "k1": k1,
                        "l1": l1,

                        "E00": E00,
                        "E10": E10,
                        "E01": E01,
                        "E11": E11,

                        "K00": k0 * l0,
                        "K10": k1 * l0,
                        "K01": k0 * l1,
                        "K11": k1 * l1,
                    })

    elapsed = (
        time.perf_counter()
        - start
    )

    # --------------------------------------------------------
    # TRUE structure.
    # --------------------------------------------------------

    true_s = true_structure(
        p_true,
        q_true,
    )

    true_found = any(
        (
            s["k0"] == true_s["k0"]
            and
            s["l0"] == true_s["l0"]
            and
            s["k1"] == true_s["k1"]
            and
            s["l1"] == true_s["l1"]
        )
        for s in survivors
    )

    return {
        "Q00": Q00,
        "Q10": Q10,
        "Q01": Q01,
        "Q11": Q11,

        "k_tests": k_tests,
        "l_tests": l_tests,
        "transition_tests":
            transition_tests,
        "final_structures":
            final_structures,

        "survivors":
            survivors,

        "true":
            true_s,

        "true_found":
            true_found,

        "time":
            elapsed,
    }


# ------------------------------------------------------------
# Print survivor
# ------------------------------------------------------------

def print_survivors(
    survivors,
):

    if not survivors:

        print("    none")
        return

    for s in survivors[
        :PRINT_SURVIVORS
    ]:

        print(
            "    "
            f"k0={s['k0']} "
            f"l0={s['l0']} "
            f"k1={s['k1']} "
            f"l1={s['l1']} "
            f"E00={s['E00']} "
            f"E10={s['E10']} "
            f"E01={s['E01']} "
            f"E11={s['E11']}"
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "2D QUOTIENT-TRANSITION FINGERPRINT"
    )
    print(
        "EXPERIMENT 99"
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
        f"r1 = [{R1_SMALL}, {R1_LARGE}]"
    )

    print(
        f"r2 = [{R2_SMALL}, {R2_LARGE}]"
    )

    # --------------------------------------------------------
    # True quotient structure.
    # --------------------------------------------------------

    ts = true_structure(
        p_true,
        q_true,
    )

    print()
    print(
        "TRUE STRUCTURE"
    )

    for name, value in ts.items():

        print(
            f"    {name} = {value}"
        )

    # --------------------------------------------------------
    # Run.
    # --------------------------------------------------------

    result = search(
        n,
        p_true,
        q_true,
    )

    # --------------------------------------------------------
    # Results.
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("RESULT")
    print("=" * 78)

    print(
        f"Q00                = "
        f"{result['Q00']}"
    )

    print(
        f"Q10                = "
        f"{result['Q10']}"
    )

    print(
        f"Q01                = "
        f"{result['Q01']}"
    )

    print(
        f"Q11                = "
        f"{result['Q11']}"
    )

    print()
    print(
        f"k values tested    = "
        f"{result['k_tests']:,}"
    )

    print(
        f"l values tested    = "
        f"{result['l_tests']:,}"
    )

    print(
        f"transition tests   = "
        f"{result['transition_tests']:,}"
    )

    print(
        f"final structures   = "
        f"{result['final_structures']:,}"
    )

    print(
        f"unique survivors   = "
        f"{len(result['survivors']):,}"
    )

    print(
        f"TRUE found         = "
        f"{result['true_found']}"
    )

    print(
        f"time               = "
        f"{result['time']:.6f} s"
    )

    # --------------------------------------------------------
    # Survivor examples.
    # --------------------------------------------------------

    print()
    print(
        "SURVIVORS"
    )

    print_survivors(
        result["survivors"]
    )

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 99 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":

    main()


# ============================================================
# FINISHED EXPERIMENT 99
# ============================================================
