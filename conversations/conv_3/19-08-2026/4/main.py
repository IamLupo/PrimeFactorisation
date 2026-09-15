from fractions import Fraction
from math import isqrt


# ============================================================
# EXPERIMENT
#   EXACT MONOTONE INVERSION OF THE LAYER ORACLE
#
# Given:
#   N
#   R_TRUE = R(N / X_TRUE)
#
# recover X_TRUE without using p or q during the search.
#
# The reference R_TRUE is generated only once at the start
# from the known hidden test pair.
# ============================================================


# ============================================================
# EXACT LAYER POLYNOMIALS
# ============================================================

def h16(t):
    return (
        -9 * t**8
        -36 * t**7
        -84 * t**6
        -126 * t**5
        -126 * t**4
        -84 * t**3
        -36 * t**2
        -9 * t
        -1
    )


def h15(t):
    return (
        88 * t**9
        +396 * t**8
        +1164 * t**7
        +2226 * t**6
        +2898 * t**5
        +2604 * t**4
        +1596 * t**3
        +639 * t**2
        +151 * t
        +16
    )


def h14(t):
    return (
        -276 * t**10
        -1380 * t**9
        -5460 * t**8
        -13560 * t**7
        -23058 * t**6
        -27510 * t**5
        -23100 * t**4
        -13410 * t**3
        -5135 * t**2
        -1169 * t
        -120
    )


def R_from_t(t):
    a = h16(t)
    b = h15(t)
    c = h14(t)

    if a == 0 or c == 0:
        raise ZeroDivisionError("R(t) denominator vanished.")

    return (b * b) / (a * c)


def R_from_x(n, x):
    return R_from_t(Fraction(n, x))


# ============================================================
# KAPPA
# ============================================================

def kappa_from_x(n, x):
    s = x - 1

    return (
        -s * s
        + (n + 1) * s
        - n * n
        + n
    )


# ============================================================
# DISCRIMINANT
# ============================================================

def discriminant(n, x):
    s = x - 1
    return s * s - 4 * n


def is_square(n):
    if n < 0:
        return False

    r = isqrt(n)

    return r * r == n


# ============================================================
# TEST INSTANCES
# ============================================================

TEST_PAIRS = (
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
)


# ============================================================
# FEASIBLE X INTERVAL FOR ODD PRIME FACTORS
#
# For p,q odd primes:
#
#   p >= 3
#   q >= 3
#
# For fixed N=pq:
#
#   p+q >= 2 sqrt(N)
#
# and since p >= 3:
#
#   q <= N/3
#
#   p+q <= N/3 + 3
#
# hence:
#
#   X = p+q+1
#
# lies inside a finite interval.
#
# We deliberately make this only a SEARCH BOUND.
# No knowledge of p or q is used by the search itself.
# ============================================================

def search_bounds_for_odd_semiprime(n):
    root = isqrt(n)

    # Safe integer lower bound:
    # X >= ceil(2*sqrt(N)) + 1
    #
    # Since sqrt(N) may not be integral, use:
    lower_s = 2 * root

    if root * root < n:
        lower_s += 1

    lower_x = lower_s + 1

    # Safe upper bound for an odd semiprime:
    # p >= 3, q = N/p <= N/3
    upper_x = n // 3 + 4

    return lower_x, upper_x


# ============================================================
# EXACT BINARY SEARCH
#
# Since R(N/X) is strictly decreasing:
#
#   X < X_TRUE  => R(X) > R_TRUE
#   X > X_TRUE  => R(X) < R_TRUE
#
# ============================================================

def binary_search_x(n, r_true, lower_x, upper_x):
    lo = lower_x
    hi = upper_x

    queries = 0
    trace = []

    while lo <= hi:
        mid = (lo + hi) // 2

        r_mid = R_from_x(n, mid)

        queries += 1

        if r_mid == r_true:
            trace.append(
                ("EQUAL", lo, hi, mid)
            )
            return mid, queries, trace

        if r_mid > r_true:
            # Since R(X) is strictly decreasing:
            # R(mid) > R(TRUE) => mid < TRUE
            trace.append(
                ("TOO_SMALL", lo, hi, mid)
            )
            lo = mid + 1

        else:
            # R(mid) < R(TRUE) => mid > TRUE
            trace.append(
                ("TOO_LARGE", lo, hi, mid)
            )
            hi = mid - 1

    return None, queries, trace


# ============================================================
# LINEAR CHECK AROUND THE RECOVERED X
# ============================================================

def verify_recovered_x(n, x_true, r_true):
    checks = (
        -1000,
        -100,
        -10,
        -2,
        -1,
        0,
        1,
        2,
        10,
        100,
        1000,
    )

    failures = 0

    for delta in checks:
        x = x_true + delta

        if x <= 0:
            continue

        r = R_from_x(n, x)

        if delta < 0 and not (r > r_true):
            failures += 1

        elif delta > 0 and not (r < r_true):
            failures += 1

        elif delta == 0 and r != r_true:
            failures += 1

    return failures


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 110)
    print("EXPERIMENT: EXACT INVERSE OF THE MONOTONE X-LAYER ORACLE")
    print("=" * 110)

    print()
    print("OBJECTIVE")
    print("  Given N and the exact reference value R_TRUE, recover X_TRUE")
    print("  by monotone search only.")
    print()
    print("  During the search:")
    print("    - p is not used")
    print("    - q is not used")
    print("    - p+q is not used")
    print("    - the discriminant is not used")
    print("    - only N, R_TRUE, and monotonicity are used")
    print()

    total_instances = 0
    recovery_failures = 0
    total_queries = 0
    max_queries = 0
    verification_failures = 0

    for instance_id, (p, q) in enumerate(TEST_PAIRS, start=1):

        n = p * q
        x_true = p + q + 1

        # This is the hidden oracle value.
        # It is generated ONCE from the hidden X.
        r_true = R_from_x(n, x_true)

        lower_x, upper_x = search_bounds_for_odd_semiprime(n)

        recovered_x, queries, trace = binary_search_x(
            n,
            r_true,
            lower_x,
            upper_x
        )

        total_instances += 1
        total_queries += queries

        if queries > max_queries:
            max_queries = queries

        if recovered_x != x_true:
            recovery_failures += 1

        local_failures = verify_recovered_x(
            n,
            x_true,
            r_true
        )

        verification_failures += local_failures

        print()
        print("=" * 110)
        print(f"INSTANCE {instance_id}")
        print("=" * 110)

        print(f"  N                = {n}")
        print(f"  hidden X         = {x_true}")
        print(f"  search lower     = {lower_x}")
        print(f"  search upper     = {upper_x}")
        print(f"  search width     = {upper_x - lower_x + 1}")

        print()
        print(f"  recovered X      = {recovered_x}")
        print(f"  binary queries   = {queries}")
        print(f"  recovery correct = {recovered_x == x_true}")

        print()
        print(
            f"  local sign verification failures = "
            f"{local_failures}"
        )

        print()
        print("  SEARCH TRACE")

        # Print only the first few and final steps
        # so output remains manageable.

        if len(trace) <= 12:
            shown = trace
        else:
            shown = (
                trace[:6]
                + [("...", 0, 0, 0)]
                + trace[-6:]
            )

        for item in shown:
            status, lo, hi, mid = item

            if status == "...":
                print("    ...")
            else:
                print(
                    f"    {status:10s} "
                    f"lo={lo} "
                    f"hi={hi} "
                    f"mid={mid}"
                )

    # ========================================================
    # GLOBAL RESULT
    # ========================================================

    print()
    print("=" * 110)
    print("GLOBAL RESULT")
    print("=" * 110)

    print()
    print(f"  instances tested        = {total_instances}")
    print(f"  recovery failures       = {recovery_failures}")
    print(f"  local verification fail = {verification_failures}")

    if total_instances:
        average_queries = total_queries / total_instances
    else:
        average_queries = 0.0

    print(f"  total binary queries    = {total_queries}")
    print(f"  average binary queries  = {average_queries:.3f}")
    print(f"  maximum binary queries  = {max_queries}")

    print()
    print(
        "  EXACT INVERSE ORACLE PASS :",
        recovery_failures == 0
    )

    print(
        "  LOCAL SIGN PASS           :",
        verification_failures == 0
    )

    # ========================================================
    # INFORMATION-THEORETIC COMPARISON
    # ========================================================

    print()
    print("=" * 110)
    print("SEARCH COMPLEXITY COMPARISON")
    print("=" * 110)

    print()
    print("  This compares the number of monotone R-queries with")
    print("  the width of the X search interval.")
    print()

    for p, q in TEST_PAIRS:

        n = p * q
        x_true = p + q + 1

        lower_x, upper_x = search_bounds_for_odd_semiprime(n)

        width = upper_x - lower_x + 1

        theoretical = 0
        w = width

        while w > 1:
            w = (w + 1) // 2
            theoretical += 1

        print(
            f"  N={n:<16d} "
            f"width={width:<16d} "
            f"ceil(log2(width))={theoretical:<8d}"
        )

    # ========================================================
    # FINAL INTERPRETATION
    # ========================================================

    print()
    print("=" * 110)
    print("INTERPRETATION")
    print("=" * 110)

    print()
    print("  The experiment establishes the inverse side of the oracle:")
    print()
    print("      (N, R_TRUE) -> unique positive X")
    print()
    print("  by monotone binary search.")
    print()
    print("  The search never tests candidate primes.")
    print("  The discriminant is not used during the search.")
    print()
    print("  IMPORTANT:")
    print("    R_TRUE is still supplied from the hidden X.")
    print("    Therefore this does NOT yet provide a factoring algorithm.")
    print()
    print("  The remaining research question is:")
    print()
    print("      Can R_TRUE be generated from N")
    print("      without first knowing X or p,q?")
    print()

    print("=" * 110)
    print("END")
    print("=" * 110)


if __name__ == "__main__":
    main()
