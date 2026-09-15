import sympy
import math
from collections import Counter


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

def find_collisions(n, actual_p, actual_q, primes):

    actual = tuple(sorted((actual_p, actual_q)))

    collisions = []

    for candidate_p in primes:

        if math.gcd(candidate_p, M) != 1:
            continue

        # p*q = n (mod M)
        candidate_q = (
            n * pow(candidate_p, -1, M)
        ) % M

        if not (
            FACTOR_LOW <= candidate_q <= FACTOR_HIGH
        ):
            continue

        if not sympy.isprime(candidate_q):
            continue

        pair = tuple(
            sorted((candidate_p, candidate_q))
        )

        if pair == actual:
            continue

        a, b = pair

        product = a * b

        difference = product - n

        if difference % TWO_M != 0:
            continue

        t = difference // TWO_M

        collisions.append(
            (a, b, t)
        )

    return sorted(set(collisions))


# ============================================================
# STORAGE
# ============================================================

all_collisions = []


# ============================================================
# MAIN
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
        n,
        actual_p,
        actual_q,
        primes
    )

    for a, b, t in collisions:

        dp = a - actual_p
        dq = b - actual_q

        # ----------------------------------------------------
        # Factor sums and differences
        # ----------------------------------------------------

        actual_sum = actual_p + actual_q
        false_sum = a + b

        actual_diff = actual_q - actual_p
        false_diff = b - a

        delta_sum = false_sum - actual_sum
        delta_diff = false_diff - actual_diff

        # ----------------------------------------------------
        # Product relation
        # ----------------------------------------------------

        product_difference = (
            a * b - actual_p * actual_q
        )

        # ----------------------------------------------------
        # Identity using sums/differences
        # ----------------------------------------------------

        left = (
            false_sum**2
            - actual_sum**2
            - false_diff**2
            + actual_diff**2
        )

        right = 4 * product_difference

        assert left == right
        assert right == 8 * M * t

        record = {
            "trial": trial,
            "p": actual_p,
            "q": actual_q,
            "a": a,
            "b": b,
            "t": t,
            "dp": dp,
            "dq": dq,
            "actual_sum": actual_sum,
            "false_sum": false_sum,
            "actual_diff": actual_diff,
            "false_diff": false_diff,
            "delta_sum": delta_sum,
            "delta_diff": delta_diff,
        }

        all_collisions.append(record)


# ============================================================
# SUMMARY
# ============================================================

print("=" * 90)
print("FACTOR-SHIFT GEOMETRY")
print("=" * 90)

print(f"M  = {M:,}")
print(f"2M = {TWO_M:,}")
print()

print(
    f"Total false collisions = "
    f"{len(all_collisions)}"
)


# ============================================================
# t DISTRIBUTION
# ============================================================

print()
print("=" * 90)
print("t DISTRIBUTION")
print("=" * 90)

t_counter = Counter(
    r["t"]
    for r in all_collisions
)

for t, count in sorted(t_counter.items()):

    print(
        f"t={t:+4}: {count:4}"
    )


# ============================================================
# SMALL t
# ============================================================

print()
print("=" * 90)
print("SMALL |t|")
print("=" * 90)

for r in sorted(
    all_collisions,
    key=lambda x: (abs(x["t"]), x["trial"])
):

    if abs(r["t"]) > 5:
        break

    print(
        f"t={r['t']:+3} "
        f"actual=({r['p']:,},{r['q']:,}) "
        f"false=({r['a']:,},{r['b']:,})"
    )

    print(
        f"    "
        f"dp={r['dp']:+,} "
        f"dq={r['dq']:+,}"
    )

    print(
        f"    "
        f"sum: "
        f"{r['actual_sum']:,}"
        f" -> "
        f"{r['false_sum']:,}"
        f" "
        f"delta={r['delta_sum']:+,}"
    )

    print(
        f"    "
        f"diff: "
        f"{r['actual_diff']:,}"
        f" -> "
        f"{r['false_diff']:,}"
        f" "
        f"delta={r['delta_diff']:+,}"
    )

    print()


# ============================================================
# TEST RELATIONSHIPS
# ============================================================

print("=" * 90)
print("RELATIONSHIP STATISTICS")
print("=" * 90)

sum_deltas = [
    r["delta_sum"]
    for r in all_collisions
]

diff_deltas = [
    r["delta_diff"]
    for r in all_collisions
]

print()
print(
    f"delta(sum) min = "
    f"{min(sum_deltas):+,}"
)

print(
    f"delta(sum) max = "
    f"{max(sum_deltas):+,}"
)

print(
    f"delta(diff) min = "
    f"{min(diff_deltas):+,}"
)

print(
    f"delta(diff) max = "
    f"{max(diff_deltas):+,}"
)


# ============================================================
# NORMALIZE BY t
# ============================================================

print()
print("=" * 90)
print("NORMALIZED SHIFT")
print("=" * 90)

for r in sorted(
    all_collisions,
    key=lambda x: abs(x["t"])
)[:50]:

    t = r["t"]

    print(
        f"t={t:+3} "
        f"dSum={r['delta_sum']:+7,} "
        f"dDiff={r['delta_diff']:+7,}"
    )


# ============================================================
# VERIFY THE CORE IDENTITY
# ============================================================

print()
print("=" * 90)
print("CORE IDENTITY")
print("=" * 90)

print(
    "(a+b)^2 - (p+q)^2"
    " - "
    "[(a-b)^2 - (p-q)^2]"
)

print(
    f"= 8*M*t"
)

print()

for r in all_collisions[:20]:

    lhs = (
        r["false_sum"]**2
        - r["actual_sum"]**2
        - r["false_diff"]**2
        + r["actual_diff"]**2
    )

    rhs = 8 * M * r["t"]

    print(
        f"t={r['t']:+3} "
        f"LHS={lhs:+,} "
        f"RHS={rhs:+,} "
        f"OK={lhs == rhs}"
    )
