from math import comb, isqrt


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

# Small values used to independently validate the fast formula.
BRUTE_FORCE_LIMIT = 100


# ============================================================
# BASIC NUMBER-THEORY ROUTINES
# ============================================================

def sigma_sieve(limit):
    """
    Compute sigma(n) for every 1 <= n <= limit.

    sigma(n) = sum of divisors of n.
    """
    sigma = [0] * (limit + 1)

    for d in range(1, limit + 1):
        for multiple in range(d, limit + 1, d):
            sigma[multiple] += d

    return sigma


def divisors(n):
    """
    Return all positive divisors of n.
    """
    result = []

    r = isqrt(n)

    for d in range(1, r + 1):
        if n % d == 0:
            result.append(d)

            other = n // d

            if other != d:
                result.append(other)

    return result


# ============================================================
# ORIGINAL BRUTE-FORCE M2
# ============================================================

def macmahon_M2_bruteforce(n):
    """
    Original direct enumeration.

    M2(n) =
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


# ============================================================
# FAST M2 FORMULA
# ============================================================

def macmahon_M2_fast(n, sigma):
    """
    Compute M2(n) from the generating-function coefficient:

        M2(n) =
            1/2 * (
                sum_{j=1}^{n-1} sigma(j)*sigma(n-j)
                -
                sum_{d|n} C(n/d + 1, 3)
            )

    This does NOT use p or q.
    """

    convolution = 0

    for j in range(1, n):
        convolution += sigma[j] * sigma[n - j]

    correction = 0

    for d in divisors(n):
        m = n // d

        # C(m+1, 3) is automatically zero for m < 2,
        # but comb handles it cleanly for our range.
        correction += comb(m + 1, 3)

    numerator = convolution - correction

    assert numerator % 2 == 0

    return numerator // 2


# ============================================================
# SEMIPRIME S RECOVERY
# ============================================================

def cubic_value(n, M2, S):
    """
    For n = p*q and S = p+q:

        S^3
        + (1 - 5n)S
        + n^3 - 2n^2 - n + 2
        - 8M2
        = 0
    """
    return (
        S**3
        + (1 - 5 * n) * S
        + n**3
        - 2 * n**2
        - n
        + 2
        - 8 * M2
    )


def recover_S(n, M2):
    """
    Recover the semiprime factor sum S.

    For n=pq:

        S >= 2*sqrt(n)
        S <= n+1

    We exploit monotonicity on this physical interval.
    """

    lo = isqrt(4 * n)

    if lo * lo < 4 * n:
        lo += 1

    hi = n + 1

    # Binary search because the cubic is strictly increasing
    # throughout the physically relevant interval.
    left = lo
    right = hi

    while left <= right:

        mid = (left + right) // 2

        value = cubic_value(n, M2, mid)

        if value == 0:
            return mid

        if value < 0:
            left = mid + 1
        else:
            right = mid - 1

    return None


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 9")
print()

print("FAST M2 FORMULA VALIDATION")
print("-" * 70)

sigma_small = sigma_sieve(BRUTE_FORCE_LIMIT)

validation_failed = False

for n in range(2, BRUTE_FORCE_LIMIT + 1):

    brute = macmahon_M2_bruteforce(n)
    fast = macmahon_M2_fast(n, sigma_small)

    ok = brute == fast

    print(
        f"n={n:3d} "
        f"bruteforce={brute:12d} "
        f"fast={fast:12d} "
        f"OK={ok}"
    )

    if not ok:
        validation_failed = True

print()

if validation_failed:
    print("VALIDATION FAILED")
    print("Stopping before semiprime test.")
else:

    print("VALIDATION PASSED")
    print()

    max_n = max(p * q for p, q in TEST_CASES)

    #
    # sigma[] is computed entirely from n.
    # No factor information is supplied.
    #
    sigma = sigma_sieve(max_n)

    print("SEMIPRIME TEST")
    print("-" * 120)

    print(
        "n\tp\tq\t"
        "M2\t"
        "true_S\t"
        "recovered_S\t"
        "S_match\t"
        "factor_discriminant\t"
        "factors_recovered"
    )

    for p, q in TEST_CASES:

        n = p * q

        #
        # IMPORTANT:
        #
        # p and q are only printed as ground truth.
        # They are NOT passed to macmahon_M2_fast().
        #
        M2 = macmahon_M2_fast(n, sigma)

        true_S = p + q

        recovered_S = recover_S(n, M2)

        S_match = (
            recovered_S == true_S
        )

        factors_recovered = False
        discriminant = None

        if recovered_S is not None:

            discriminant = (
                recovered_S * recovered_S
                - 4 * n
            )

            if discriminant >= 0:

                root = isqrt(discriminant)

                if root * root == discriminant:

                    recovered_p = (
                        recovered_S + root
                    ) // 2

                    recovered_q = (
                        recovered_S - root
                    ) // 2

                    factors_recovered = (
                        recovered_p * recovered_q == n
                    )

        print(
            f"{n}\t"
            f"{p}\t"
            f"{q}\t"
            f"{M2}\t"
            f"{true_S}\t"
            f"{recovered_S}\t"
            f"{S_match}\t"
            f"{discriminant}\t"
            f"{factors_recovered}"
        )

print()
print("FINISHED EXPERIMENT 9")
