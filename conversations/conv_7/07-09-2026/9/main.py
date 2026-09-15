from collections import defaultdict
from sympy import primerange


# ============================================================
# SETTINGS
# ============================================================

TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
]

R1 = 11
R2 = 13


# ============================================================
# MACMAHON M2
# ============================================================

def macmahon_M2_terms(n):
    """
    Generate all states contributing to M2(n):

        m1*s1 + m2*s2 = n

    with

        0 < s1 < s2
        m1,m2 > 0

    Each state contributes weight m1*m2.
    """

    for s1 in range(1, n):

        for s2 in range(s1 + 1, n + 1):

            for m1 in range(1, n // s1 + 1):

                remainder = n - m1 * s1

                if remainder <= 0:
                    continue

                if remainder % s2 != 0:
                    continue

                m2 = remainder // s2

                if m2 <= 0:
                    continue

                yield s1, s2, m1, m2, m1 * m2


# ============================================================
# QUOTIENT / REMAINDER
# ============================================================

def qr(x, r):
    return x % r, x // r


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 6")
print()

print(
    "n\tp\tq\t"
    "a\tb\tk\tl\tK\tS\t"
    "M2\tstates\t"
    "factor_weight\t"
    "factor_weight_pct\t"
    "residue_match_weight\t"
    "residue_match_pct\t"
    "quotient_match_weight\t"
    "quotient_match_pct\t"
    "exact_factor_occurrences"
)

for p, q in TEST_CASES:

    n = p * q

    # Factor coordinates
    a, k = qr(p, R1)
    b, l = qr(q, R2)

    K = k * l
    S = p + q

    total_weight = 0
    state_count = 0

    factor_weight = 0
    residue_match_weight = 0
    quotient_match_weight = 0

    exact_factor_occurrences = 0

    # --------------------------------------------------------
    # Generate M2 states
    # --------------------------------------------------------

    for s1, s2, m1, m2, weight in macmahon_M2_terms(n):

        total_weight += weight
        state_count += 1

        # ----------------------------------------------------
        # 1. Does the M2 state literally contain p or q?
        # ----------------------------------------------------

        contains_factor = (
            p in (s1, s2, m1, m2)
            or q in (s1, s2, m1, m2)
        )

        if contains_factor:
            factor_weight += weight

        exact_count = (
            (s1 == p)
            or (s2 == p)
            or (m1 == p)
            or (m2 == p)
            or
            (s1 == q)
            or (s2 == q)
            or (m1 == q)
            or (m2 == q)
        )

        if exact_count:
            exact_factor_occurrences += 1

        # ----------------------------------------------------
        # 2. Residue signature match
        #
        # Check whether BOTH:
        #
        #   some partition quantity has residue a mod R1
        #   AND another has residue b mod R2
        #
        # This is deliberately symmetric.
        # ----------------------------------------------------

        values = (s1, s2, m1, m2)

        has_a = any(x % R1 == a for x in values)
        has_b = any(x % R2 == b for x in values)

        if has_a and has_b:
            residue_match_weight += weight

        # ----------------------------------------------------
        # 3. Quotient signature match
        #
        # Look for quotient coordinates k and l.
        # ----------------------------------------------------

        has_k = any(x // R1 == k for x in values)
        has_l = any(x // R2 == l for x in values)

        if has_k and has_l:
            quotient_match_weight += weight

    # --------------------------------------------------------
    # Percentages
    # --------------------------------------------------------

    factor_pct = (
        100.0 * factor_weight / total_weight
        if total_weight
        else 0.0
    )

    residue_pct = (
        100.0 * residue_match_weight / total_weight
        if total_weight
        else 0.0
    )

    quotient_pct = (
        100.0 * quotient_match_weight / total_weight
        if total_weight
        else 0.0
    )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{a}\t"
        f"{b}\t"
        f"{k}\t"
        f"{l}\t"
        f"{K}\t"
        f"{S}\t"
        f"{total_weight}\t"
        f"{state_count}\t"
        f"{factor_weight}\t"
        f"{factor_pct:.8f}\t"
        f"{residue_match_weight}\t"
        f"{residue_pct:.8f}\t"
        f"{quotient_match_weight}\t"
        f"{quotient_pct:.8f}\t"
        f"{exact_factor_occurrences}"
    )

print()
print("FINISHED EXPERIMENT 6")
