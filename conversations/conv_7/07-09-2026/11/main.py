from math import isqrt


TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
    (47, 83),
    (53, 89),
]


def macmahon_M2(n):
    """
    Exact MacMahon M2(n):

        sum m1*m2

    over

        n = m1*s1 + m2*s2

    with

        0 < s1 < s2
        m1,m2 > 0
    """
    total = 0

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

                total += m1 * m2

    return total


def cubic_value(n, M2, S):
    """
    Cubic whose root is S = p + q for semiprime n=pq:

        f(S) =
            S^3
            + (1 - 5n)S
            + n^3 - 2n^2 - n + 2
            - 8*M2
    """
    return (
        S**3
        + (1 - 5*n) * S
        + n**3
        - 2*n**2
        - n
        + 2
        - 8*M2
    )


def recover_S_integer(n, M2):
    """
    Search only the mathematically valid interval.

    For p*q=n:

        2*sqrt(n) <= p+q <= n+1

    We use integer bounds and test exact cubic equality.
    """

    lo = 2 * isqrt(n)

    # Correct the lower bound when n is not a square.
    if lo * lo < 4 * n:
        lo += 1

    hi = n + 1

    for S in range(lo, hi + 1):
        if cubic_value(n, M2, S) == 0:
            return S

    return None


print("START EXPERIMENT 8")
print()

print(
    "n\tp\tq\t"
    "M2\t"
    "true_S\t"
    "recovered_S\t"
    "match\t"
    "cubic_at_true_S"
)

for p, q in TEST_CASES:

    n = p * q

    true_S = p + q

    #
    # This calculation knows n only.
    # p and q are NOT used in the M2 calculation.
    #
    x = (p**3 + 1) * (q**3 + 1)
    M1 = (p + 1) * (q + 1)
    M2 = (x - (2*n*M1) + M1) // 8
    
    #M2 = macmahon_M2(n)

    recovered_S = recover_S_integer(n, M2)

    check = (
        recovered_S == true_S
    )

    cubic_check = cubic_value(
        n,
        M2,
        true_S
    )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{M2}\t"
        f"{true_S}\t"
        f"{recovered_S}\t"
        f"{check}\t"
        f"{cubic_check}"
    )


print()
print("FINISHED EXPERIMENT 8")
