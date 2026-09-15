import math
import secrets
import sympy
from collections import defaultdict, Counter


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 300

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

M = math.prod(MODULI)
STEP = 2 * M

# Only investigate small t.
T_LIMIT = 25


# ============================================================
# RANDOM PRIME
# ============================================================

def random_prime(low, high):
    """
    Force independent randomness instead of relying on the
    global RNG state used by sympy.randprime().
    """

    while True:

        x = secrets.randbelow(
            high - low + 1
        ) + low

        if x % 2 == 0:
            x += 1

        if x > high:
            continue

        if sympy.isprime(x):
            return x


# ============================================================
# FIND COLLISIONS
# ============================================================

def find_collisions(n, p, q):

    actual = tuple(sorted((p, q)))

    sqrt_n = math.isqrt(n)

    primes = sympy.primerange(
        2,
        sqrt_n + 1
    )

    collisions = []

    for a in primes:

        if math.gcd(a, M) != 1:
            continue

        b = (
            n * pow(a, -1, M)
        ) % M

        if not (
            FACTOR_LOW
            <= b
            <= FACTOR_HIGH
        ):
            continue

        if not sympy.isprime(b):
            continue

        pair = tuple(sorted((a, b)))

        if pair == actual:
            continue

        difference = a * b - n

        if difference % STEP != 0:
            continue

        t = difference // STEP

        if abs(t) > T_LIMIT:
            continue

        collisions.append(
            (pair[0], pair[1], t)
        )

    return sorted(set(collisions))


# ============================================================
# STORAGE
# ============================================================

seen_n = set()

events = []

by_t = defaultdict(list)

transformation_counter = Counter()


# ============================================================
# MAIN
# ============================================================

print("=" * 100)
print("UNIQUE-SEMIPRIME COLLISION EXPERIMENT")
print("=" * 100)

print(f"M       = {M:,}")
print(f"2M      = {STEP:,}")
print(
    f"Factor range = "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)
print(
    f"Unique trials requested = {TRIALS}"
)
print(
    f"|t| <= {T_LIMIT}"
)
print()


trial = 0

while trial < TRIALS:

    p = random_prime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    q = random_prime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    p, q = sorted((p, q))

    n = p * q

    # --------------------------------------------------------
    # Guarantee unique n
    # --------------------------------------------------------

    if n in seen_n:
        continue

    seen_n.add(n)

    trial += 1

    collisions = find_collisions(
        n,
        p,
        q
    )

    # --------------------------------------------------------
    # Record
    # --------------------------------------------------------

    print(
        f"trial={trial:3} "
        f"p={p:,} "
        f"q={q:,} "
        f"n={n:,} "
        f"collisions={len(collisions)}"
    )

    for a, b, t in collisions:

        dp = a - p
        dq = b - q

        ds = (a + b) - (p + q)
        dd = (b - a) - (q - p)

        # Transformation signature.
        signature = (
            t,
            dp,
            dq
        )

        transformation_counter[
            signature
        ] += 1

        record = {
            "n": n,
            "p": p,
            "q": q,
            "a": a,
            "b": b,
            "t": t,
            "dp": dp,
            "dq": dq,
            "ds": ds,
            "dd": dd,
        }

        events.append(record)

        by_t[t].append(record)

        print(
            f"    t={t:+3} "
            f"false=({a:,},{b:,}) "
            f"dp={dp:+,} "
            f"dq={dq:+,}"
        )


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    f"Unique semiprimes tested = {len(seen_n)}"
)

print(
    f"Total collision events    = {len(events)}"
)

print()


# ============================================================
# T DISTRIBUTION
# ============================================================

print("=" * 100)
print("t DISTRIBUTION")
print("=" * 100)

for t in sorted(by_t):

    print(
        f"t={t:+3} : "
        f"{len(by_t[t]):4}"
    )


# ============================================================
# UNIQUE COLLISION SIGNATURES
# ============================================================

print()
print("=" * 100)
print("REPEATED TRANSFORMATIONS")
print("=" * 100)

repeated = [
    (signature, count)
    for signature, count
    in transformation_counter.items()
    if count > 1
]

if not repeated:

    print(
        "No exact (t, dp, dq) transformation "
        "was repeated."
    )

else:

    for (
        (t, dp, dq),
        count
    ) in sorted(
        repeated,
        key=lambda x: -x[1]
    ):

        print(
            f"count={count:3} "
            f"t={t:+3} "
            f"dp={dp:+,} "
            f"dq={dq:+,}"
        )


# ============================================================
# PER-t DELTA STATISTICS
# ============================================================

print()
print("=" * 100)
print("PER-t DELTA RANGES")
print("=" * 100)

for t in sorted(by_t):

    records = by_t[t]

    dp_values = [
        r["dp"]
        for r in records
    ]

    dq_values = [
        r["dq"]
        for r in records
    ]

    ds_values = [
        r["ds"]
        for r in records
    ]

    dd_values = [
        r["dd"]
        for r in records
    ]

    print()
    print(
        f"t={t:+3} "
        f"count={len(records)}"
    )

    print(
        f"    dp range = "
        f"{min(dp_values):+,} ... "
        f"{max(dp_values):+,}"
    )

    print(
        f"    dq range = "
        f"{min(dq_values):+,} ... "
        f"{max(dq_values):+,}"
    )

    print(
        f"    dSum range = "
        f"{min(ds_values):+,} ... "
        f"{max(ds_values):+,}"
    )

    print(
        f"    dDiff range = "
        f"{min(dd_values):+,} ... "
        f"{max(dd_values):+,}"
    )


# ============================================================
# MOST INTERESTING: t = +/-1
# ============================================================

print()
print("=" * 100)
print("t = +/-1")
print("=" * 100)

for t in [-1, 1]:

    records = by_t.get(t, [])

    print()
    print(
        f"t={t:+d} "
        f"count={len(records)}"
    )

    for r in records:

        print(
            f"  n={r['n']:,} "
            f"actual=({r['p']:,},{r['q']:,}) "
            f"false=({r['a']:,},{r['b']:,})"
        )

        print(
            f"      "
            f"dp={r['dp']:+,} "
            f"dq={r['dq']:+,} "
            f"dSum={r['ds']:+,} "
            f"dDiff={r['dd']:+,}"
        )


# ============================================================
# CHECK SAME t FOR LINEAR RELATIONSHIPS
# ============================================================

print()
print("=" * 100)
print("NORMALIZED t RELATIONSHIPS")
print("=" * 100)

for t in sorted(by_t):

    records = by_t[t]

    print()
    print(
        f"t={t:+3}"
    )

    for r in records[:20]:

        p = r["p"]
        q = r["q"]

        dp = r["dp"]
        dq = r["dq"]

        # Equation:
        #
        # p*dq + q*dp + dp*dq = 2*M*t

        linear = (
            p * dq
            + q * dp
        )

        quadratic = dp * dq

        print(
            f"  dp={dp:+8,} "
            f"dq={dq:+8,} "
            f"linear={linear:+14,} "
            f"quadratic={quadratic:+14,}"
        )
