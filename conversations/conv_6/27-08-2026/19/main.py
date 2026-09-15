import sympy
import math
from collections import Counter, defaultdict


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 500

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

M = math.prod(MODULI)
TWO_M = 2 * M


# ============================================================
# FIND COLLISIONS
# ============================================================

def find_collisions(n, p, q, primes):

    actual = tuple(sorted((p, q)))

    collisions = []

    for a in primes:

        if math.gcd(a, M) != 1:
            continue

        b = (n * pow(a, -1, M)) % M

        if not (FACTOR_LOW <= b <= FACTOR_HIGH):
            continue

        if not sympy.isprime(b):
            continue

        pair = tuple(sorted((a, b)))

        if pair == actual:
            continue

        product = a * b
        difference = product - n

        if difference % TWO_M != 0:
            continue

        t = difference // TWO_M

        collisions.append(
            (pair[0], pair[1], t)
        )

    return sorted(set(collisions))


# ============================================================
# STORAGE
# ============================================================

all_collisions = []

t_counter = Counter()

# Store collisions grouped by t
by_t = defaultdict(list)


# ============================================================
# MAIN EXPERIMENT
# ============================================================

print("=" * 90)
print("FIXED-t COLLISION EXPERIMENT")
print("=" * 90)

print(f"M  = {M:,}")
print(f"2M = {TWO_M:,}")
print()

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

        actual_p, actual_q = sorted((p, q))

        dp = a - actual_p
        dq = b - actual_q

        by_t[t].append(
            {
                "trial": trial,
                "p": actual_p,
                "q": actual_q,
                "a": a,
                "b": b,
                "n": n,
                "t": t,
                "dp": dp,
                "dq": dq,
            }
        )

        t_counter[t] += 1


# ============================================================
# t DISTRIBUTION
# ============================================================

print("=" * 90)
print("t DISTRIBUTION")
print("=" * 90)

for t, count in sorted(t_counter.items()):

    print(
        f"t={t:+4} : {count:4}"
    )


# ============================================================
# FOCUS ON SMALL |t|
# ============================================================

print()
print("=" * 90)
print("SMALL |t| COLLISIONS")
print("=" * 90)

small_t = sorted(
    [
        t
        for t in by_t
        if abs(t) <= 5
    ]
)

for t in small_t:

    records = by_t[t]

    print()
    print(
        f"t = {t:+d}   "
        f"count = {len(records)}"
    )

    print(
        f"Equation: "
        f"ab = n {'+' if t >= 0 else '-'} "
        f"{abs(2 * M * t):,}"
    )

    for r in records:

        print(
            f"  trial={r['trial']:3} "
            f"actual=({r['p']:,},{r['q']:,}) "
            f"false=({r['a']:,},{r['b']:,}) "
            f"dp={r['dp']:+,} "
            f"dq={r['dq']:+,}"
        )


# ============================================================
# ANALYZE t = +/-1
# ============================================================

print()
print("=" * 90)
print("t = +/-1 ANALYSIS")
print("=" * 90)

for t in [-1, 1]:

    records = by_t.get(t, [])

    print()
    print(
        f"t={t:+d}: "
        f"{len(records)} collisions"
    )

    if not records:
        continue

    # --------------------------------------------------------
    # Sum/difference of factor deltas
    # --------------------------------------------------------

    sums = [
        r["dp"] + r["dq"]
        for r in records
    ]

    products = [
        r["dp"] * r["dq"]
        for r in records
    ]

    print(
        f"dp+dq range: "
        f"{min(sums):+,} ... {max(sums):+,}"
    )

    print(
        f"dp*dq range: "
        f"{min(products):+,} ... {max(products):+,}"
    )

    # --------------------------------------------------------
    # Look for repeated delta patterns
    # --------------------------------------------------------

    sum_counter = Counter(sums)

    print()
    print("Most common dp+dq values:")

    for value, count in sum_counter.most_common(10):

        print(
            f"  {value:+,} : {count}"
        )


# ============================================================
# CHECK A POSSIBLE SYMMETRIC EQUATION
# ============================================================

print()
print("=" * 90)
print("NORMALIZED DELTA RELATION")
print("=" * 90)

for t in [-1, 1]:

    records = by_t.get(t, [])

    if not records:
        continue

    print()
    print(f"t = {t:+d}")

    for r in records:

        p = r["p"]
        q = r["q"]

        dp = r["dp"]
        dq = r["dq"]

        # From:
        #
        # p*dq + q*dp + dp*dq = 2*M*t
        #
        # Divide by p*q:
        #
        # dq/q + dp/p + dp*dq/(p*q)
        # = 2*M*t/(p*q)

        normalized = (
            dp / p
            + dq / q
            + (dp * dq) / (p * q)
        )

        rhs = (
            TWO_M * t
        ) / (p * q)

        print(
            f"actual=({p:,},{q:,}) "
            f"false=({r['a']:,},{r['b']:,}) "
            f"normalized={normalized:.12e} "
            f"rhs={rhs:.12e}"
        )


# ============================================================
# FACTORING n + 2Mt
# ============================================================

print()
print("=" * 90)
print("DIRECT FACTORIZATION OF n + 2Mt")
print("=" * 90)

tested_t = sorted(
    [
        t
        for t in by_t
        if abs(t) <= 5
    ]
)

seen = set()

for t in tested_t:

    for r in by_t[t]:

        key = (
            r["n"],
            t
        )

        if key in seen:
            continue

        seen.add(key)

        value = r["n"] + TWO_M * t

        factors = sympy.factorint(value)

        print()
        print(
            f"n = {r['n']:,}"
        )

        print(
            f"t = {t:+d}"
        )

        print(
            f"n + 2Mt = {value:,}"
        )

        print(
            f"factorization = {factors}"
        )


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 90)
print("SUMMARY")
print("=" * 90)

print(
    f"Total collisions: {len(all_collisions)}"
)

for t in [-3, -2, -1, 1, 2, 3]:

    print(
        f"t={t:+d}: "
        f"{len(by_t.get(t, []))} collisions"
    )
