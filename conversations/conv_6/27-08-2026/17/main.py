import sympy
import math
from collections import Counter


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

M = math.prod(MODULI)

TRIALS = 100


# ============================================================
# FIND FALSE COLLISIONS
# ============================================================

def find_collisions(p, q, n, primes):

    actual = tuple(sorted((p, q)))

    collisions = []

    for a in primes:

        if math.gcd(a, M) != 1:
            continue

        b = (n * pow(a, -1, M)) % M

        if not (
            FACTOR_LOW <= b <= FACTOR_HIGH
        ):
            continue

        if not sympy.isprime(b):
            continue

        pair = tuple(sorted((a, b)))

        if pair == actual:
            continue

        product = a * b

        difference = product - n

        if difference % M != 0:
            continue

        k = difference // M

        collisions.append(
            (
                pair[0],
                pair[1],
                k,
                product
            )
        )

    return sorted(set(collisions))


# ============================================================
# STATISTICS
# ============================================================

k_values = []

delta_p_values = []
delta_q_values = []

records = []


# ============================================================
# RUN
# ============================================================

for trial in range(1, TRIALS + 1):

    p = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    q = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    n = p * q

    actual_p, actual_q = sorted((p, q))

    sqrt_n = math.isqrt(n)

    primes = list(
        sympy.primerange(
            2,
            sqrt_n + 1
        )
    )

    collisions = find_collisions(
        p,
        q,
        n,
        primes
    )

    for a, b, k, product in collisions:

        dp = a - actual_p
        dq = b - actual_q

        # Because the pairs are sorted, this is not necessarily
        # the pairing of corresponding factors. So also examine
        # the cross pairing.

        cross_dp = a - actual_q
        cross_dq = b - actual_p

        k_values.append(k)

        records.append(
            {
                "trial": trial,
                "p": actual_p,
                "q": actual_q,
                "a": a,
                "b": b,
                "n": n,
                "k": k,
                "dp": dp,
                "dq": dq,
                "cross_dp": cross_dp,
                "cross_dq": cross_dq,
                "product": product,
            }
        )


# ============================================================
# K DISTRIBUTION
# ============================================================

print("=" * 90)
print("COLLISION K DISTRIBUTION")
print("=" * 90)

print(f"M = {M:,}")
print(f"Total false collisions = {len(k_values)}")
print()

for k, count in sorted(
    Counter(k_values).items()
):

    print(
        f"k = {k:+4} : {count:3}"
    )


# ============================================================
# COLLISION EQUATION
# ============================================================

print()
print("=" * 90)
print("COLLISION EQUATIONS")
print("=" * 90)

for r in records:

    target = r["n"] + r["k"] * M

    print()
    print(
        f"Trial {r['trial']}"
    )

    print(
        f"Actual: "
        f"{r['p']:,} * {r['q']:,}"
        f" = {r['n']:,}"
    )

    print(
        f"False:  "
        f"{r['a']:,} * {r['b']:,}"
        f" = {r['product']:,}"
    )

    print(
        f"k = {r['k']:+}"
    )

    print(
        f"n + k*M = {target:,}"
    )

    print(
        f"Equality check = "
        f"{target == r['product']}"
    )

    print(
        f"factor delta: "
        f"({r['a'] - r['p']:+,}, "
        f"{r['b'] - r['q']:+,})"
    )


# ============================================================
# LOOK FOR SYMMETRIC RELATIONSHIPS
# ============================================================

print()
print("=" * 90)
print("K VS FACTOR DELTAS")
print("=" * 90)

for r in records:

    print(
        f"k={r['k']:+4}  "
        f"dp={r['a']-r['p']:+7,}  "
        f"dq={r['b']-r['q']:+7,}  "
        f"|k|*M={abs(r['k'])*M:,}"
    )


# ============================================================
# FACTOR DIFFERENCE PRODUCT IDENTITY
# ============================================================

print()
print("=" * 90)
print("DIFFERENCE IDENTITY")
print("=" * 90)

for r in records:

    p = r["p"]
    q = r["q"]

    a = r["a"]
    b = r["b"]

    left = a * b - p * q

    # Exact expansion:
    #
    # (p + dp)(q + dq) - pq
    #
    # = p*dq + q*dp + dp*dq

    dp = a - p
    dq = b - q

    right = (
        p * dq
        + q * dp
        + dp * dq
    )

    print()
    print(
        f"k={r['k']:+4}"
    )

    print(
        f"ab-pq = {left:+,}"
    )

    print(
        f"p*dq + q*dp + dp*dq = {right:+,}"
    )

    print(
        f"Equals: {left == right}"
    )

    print(
        f"= k*M: "
        f"{r['k'] * M:+,}"
    )
