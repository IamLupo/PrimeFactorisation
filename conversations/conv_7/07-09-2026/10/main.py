from collections import defaultdict
from math import floor


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


def macmahon_M2_terms(n):
    """
    Generate the M2(n) states:

        n = m1*s1 + m2*s2
        s1 < s2
        m1,m2 > 0

    with weight m1*m2.
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


def coordinate(x, r):
    """
    Complete modular coordinate:

        x = residue + quotient*r

    Returns:
        (residue, quotient)
    """
    return x % r, x // r


def contains_coordinate(values, target, r):
    """
    True when any value in 'values' has exactly the target
    modular coordinate relative to modulus r.
    """
    return any(coordinate(x, r) == target for x in values)


print("START EXPERIMENT 7")
print()

print(
    "n\tp\tq\t"
    "a\tb\tk\tl\tK\tS\t"
    "M2\tstates\t"
    "joint_factor_weight\t"
    "joint_factor_pct\t"
    "p_coordinate_weight\t"
    "p_coordinate_pct\t"
    "q_coordinate_weight\t"
    "q_coordinate_pct\t"
    "paired_coordinate_weight\t"
    "paired_coordinate_pct\t"
    "exact_factor_weight\t"
    "exact_factor_pct\t"
    "factor_occurrences"
)


for p, q in TEST_CASES:

    n = p * q

    a = p % R1
    k = p // R1

    b = q % R2
    l = q // R2

    K = k * l
    S = p + q

    target_p = (a, k)
    target_q = (b, l)

    total_weight = 0
    state_count = 0

    p_coordinate_weight = 0
    q_coordinate_weight = 0

    paired_coordinate_weight = 0
    joint_factor_weight = 0

    exact_factor_weight = 0
    factor_occurrences = 0

    for s1, s2, m1, m2, weight in macmahon_M2_terms(n):

        values = (s1, s2, m1, m2)

        total_weight += weight
        state_count += 1

        #
        # Exact modular coordinate match for p.
        #
        has_p_coord = contains_coordinate(
            values,
            target_p,
            R1
        )

        #
        # Exact modular coordinate match for q.
        #
        has_q_coord = contains_coordinate(
            values,
            target_q,
            R2
        )

        if has_p_coord:
            p_coordinate_weight += weight

        if has_q_coord:
            q_coordinate_weight += weight

        #
        # Both factor coordinates occur somewhere
        # in the same M2 state.
        #
        if has_p_coord and has_q_coord:
            paired_coordinate_weight += weight

        #
        # Strongest version:
        #
        # A state must contain the actual value p AND
        # the actual value q.
        #
        if p in values and q in values:
            joint_factor_weight += weight

        #
        # Any state containing either exact factor.
        #
        if p in values or q in values:
            exact_factor_weight += weight

        #
        # Count occurrences rather than weight.
        #
        if p in values:
            factor_occurrences += 1

        if q in values:
            factor_occurrences += 1

    #
    # Percentages.
    #
    p_coord_pct = (
        100.0 * p_coordinate_weight / total_weight
        if total_weight else 0.0
    )

    q_coord_pct = (
        100.0 * q_coordinate_weight / total_weight
        if total_weight else 0.0
    )

    paired_coord_pct = (
        100.0 * paired_coordinate_weight / total_weight
        if total_weight else 0.0
    )

    joint_factor_pct = (
        100.0 * joint_factor_weight / total_weight
        if total_weight else 0.0
    )

    exact_factor_pct = (
        100.0 * exact_factor_weight / total_weight
        if total_weight else 0.0
    )

    print(
        f"{n}\t{p}\t{q}\t"
        f"{a}\t{b}\t{k}\t{l}\t{K}\t{S}\t"
        f"{total_weight}\t{state_count}\t"
        f"{paired_coordinate_weight}\t"
        f"{paired_coord_pct:.8f}\t"
        f"{p_coordinate_weight}\t"
        f"{p_coord_pct:.8f}\t"
        f"{q_coordinate_weight}\t"
        f"{q_coord_pct:.8f}\t"
        f"{paired_coordinate_weight}\t"
        f"{paired_coord_pct:.8f}\t"
        f"{joint_factor_weight}\t"
        f"{joint_factor_pct:.8f}\t"
        f"{exact_factor_weight}\t"
        f"{exact_factor_pct:.8f}\t"
        f"{factor_occurrences}"
    )


print()
print("FINISHED EXPERIMENT 7")

