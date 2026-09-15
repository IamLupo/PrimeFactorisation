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

TWO_M = 2 * M

TRIALS = 100


# ============================================================
# FIND FALSE COLLISIONS
# ============================================================

def find_collisions(n, p, q, primes):

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

        if difference % TWO_M != 0:
            raise RuntimeError(
                "Expected difference to be divisible by 2M."
            )

        t = difference // TWO_M

        collisions.append(
            (
                pair[0],
                pair[1],
                t
            )
        )

    return sorted(set(collisions))


# ============================================================
# MAIN
# ============================================================

all_collisions = []

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
        n,
        p,
        q,
        primes
    )

    for a, b, t in collisions:

        dp = a - actual_p
        dq = b - actual_q

        # Exact expansion:
        #
        # (p+dp)(q+dq)-pq
        #
        # = p*dq + q*dp + dp*dq

        linear = (
            actual_p * dq
            + actual_q * dp
        )

        quadratic = dp * dq

        total = linear + quadratic

        assert total == 2 * M * t

        all_collisions.append(
            {
                "trial": trial,
                "p": actual_p,
                "q": actual_q,
                "a": a,
                "b": b,
                "dp": dp,
                "dq": dq,
                "t": t,
                "linear": linear,
                "quadratic": quadratic,
            }
        )


# ============================================================
# SUMMARY
# ============================================================

print("=" * 90)
print("NORMALIZED COLLISION EXPERIMENT")
print("=" * 90)

print(f"M    = {M:,}")
print(f"2M   = {TWO_M:,}")
print()

print(
    f"Total false collisions = "
    f"{len(all_collisions)}"
)

print()


# ============================================================
# t DISTRIBUTION
# ============================================================

print("=" * 90)
print("t DISTRIBUTION")
print("=" * 90)

counter = Counter(
    x["t"]
    for x in all_collisions
)

for t, count in sorted(counter.items()):

    print(
        f"t = {t:+4} : "
        f"{count:3}"
    )


# ============================================================
# COLLISION DETAILS
# ============================================================

print()
print("=" * 90)
print("COLLISION DETAILS")
print("=" * 90)

for x in all_collisions:

    print(
        f"trial={x['trial']:3} "
        f"t={x['t']:+4} "
        f"dp={x['dp']:+7,} "
        f"dq={x['dq']:+7,} "
        f"dp+dq={x['dp'] + x['dq']:+7,}"
    )

    print(
        f"    linear   = "
        f"{x['linear']:+,}"
    )

    print(
        f"    quadratic= "
        f"{x['quadratic']:+,}"
    )

    print(
        f"    total    = "
        f"{x['linear'] + x['quadratic']:+,}"
    )

    print(
        f"    2*M*t    = "
        f"{2 * M * x['t']:+,}"
    )


# ============================================================
# SEARCH FOR SPECIAL RELATIONSHIPS
# ============================================================

print()
print("=" * 90)
print("RELATIONSHIP TESTS")
print("=" * 90)

# Test whether t has a relationship with dp*dq.
same_sign = 0
opposite_sign = 0

for x in all_collisions:

    if x["dp"] * x["dq"] > 0:
        same_sign += 1

    elif x["dp"] * x["dq"] < 0:
        opposite_sign += 1


print(
    f"Same-sign dp,dq     : {same_sign}"
)

print(
    f"Opposite-sign dp,dq  : {opposite_sign}"
)


# ============================================================
# SMALLEST |t| COLLISIONS
# ============================================================

print()
print("=" * 90)
print("SMALLEST |t| COLLISIONS")
print("=" * 90)

smallest = sorted(
    all_collisions,
    key=lambda x: abs(x["t"])
)

for x in smallest[:30]:

    print(
        f"t={x['t']:+4} "
        f"actual=({x['p']:,},{x['q']:,}) "
        f"false=({x['a']:,},{x['b']:,}) "
        f"dp={x['dp']:+,} "
        f"dq={x['dq']:+,}"
    )
