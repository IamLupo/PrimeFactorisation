#!/usr/bin/env python3
"""
EXPERIMENT 506
======================================================================
2020 KERNEL: v/g CLASSIFICATION, CONIC REDUCTION, AND FACTOR RECOVERY
======================================================================

This experiment validates the historical construction supplied by the
research notes, then reduces the (x,y) equation to an explicit
difference-of-squares factorization of n.

No external factorization routine is used for the reconstruction tests.
The known p,q are used only to validate the identities.

For odd-prime products:

    n = p*q
    b = floor(n/2)
    v = b^2 mod n
    g = floor(b/2) + 1

Historical equation:

    y^2 - x^2 + 3x - 2 = v

Completed-square form:

    (2y)^2 - (2x-3)^2 = 4v - 1

For odd n:

    n == 3 (mod 4)  => v = (n+1)/4  => 4v-1 = n
    n == 1 (mod 4)  => v = (3n+1)/4 => 4v-1 = 3n

Thus the old (x,y) search is exactly a factor-pair search:
    A = 2y - 2x + 3
    B = 2y + 2x - 3

    A*B = n          if n == 3 (mod 4)
    A*B = 3n         if n == 1 (mod 4)

This experiment also audits the historical gcd claim and reports the
correct factor exposed by b-x-y+2.
"""

from __future__ import annotations

from math import gcd, isqrt
from itertools import combinations


# ----------------------------------------------------------------------
# Prime generation
# ----------------------------------------------------------------------

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []
    prime = bytearray(b"\x01") * (limit + 1)
    prime[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            start = p * p
            prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )
    return [i for i in range(2, limit + 1) if prime[i]]


def classify_residue(n: int) -> str:
    b = n // 2
    v = (b * b) % n
    g = b // 2 + 1
    return "v==g" if v == g else "v!=g"


def kernel_values(n: int) -> tuple[int, int, int]:
    b = n // 2
    v = (b * b) % n
    g = b // 2 + 1
    return b, v, g


def find_xy(n: int, v: int, max_x: int | None = None):
    """
    Solve y^2 - x^2 + 3x - 2 = v exactly by scanning x and checking
    whether the resulting y^2 is a square.

    This is deliberately a direct validation/search of the historical
    conic, not a use of p,q.
    """
    # Use an exact bound based on a known factor-pair solution:
    # A*B = n or 3n and A,B are positive in the normal branch.
    # For a general demonstration use a conservative finite bound.
    if max_x is None:
        max_x = max(100, 4 * isqrt(3 * n) + 20)

    for x in range(-max_x, max_x + 1):
        y2 = v + x * x - 3 * x + 2
        if y2 < 0:
            continue
        y = isqrt(y2)
        if y * y == y2:
            return x, y
    return None


def reconstructed_factors_from_xy(n: int, x: int, y: int, v: int):
    A = 2 * y - 2 * x + 3
    B = 2 * y + 2 * x - 3

    product = A * B
    if product == n:
        return [A, B], A, B, 1
    if product == 3 * n:
        # One of A,B is normally divisible by 3.
        if B % 3 == 0:
            return [A, B // 3], A, B, 3
        if A % 3 == 0:
            return [A // 3, B], A, B, 3
    return None, A, B, None


# ----------------------------------------------------------------------
# Closed-form historical x,y constructions
# ----------------------------------------------------------------------

def historical_xy(p: int, q: int):
    n = p * q
    b, v, g = kernel_values(n)

    if v == g:
        x_num = q - p + 6
        y_num = p + q
        if x_num % 4 != 0 or y_num % 4 != 0:
            return None
        return x_num // 4, y_num // 4, "v==g"

    # The supplied historical formulas split according to y ? 2x.
    # Rather than assuming which orientation p,q has, test both
    # algebraically valid oriented forms.
    candidates = []

    # Branch A
    x_num = 3 * p - q + 6
    y_num = 3 * p + q
    if x_num % 4 == 0 and y_num % 4 == 0:
        candidates.append((x_num // 4, y_num // 4, "v!=g:A"))

    # Branch B
    x_num = q - 3 * p + 6
    y_num = 3 * p + q
    if x_num % 4 == 0 and y_num % 4 == 0:
        candidates.append((x_num // 4, y_num // 4, "v!=g:B"))

    # Prefer the candidate satisfying the historical inequality.
    for x, y, tag in candidates:
        if y > 2 * x:
            return x, y, tag
    for x, y, tag in candidates:
        if y < 2 * x:
            return x, y, tag
    # Equality case y=2x occurs at a boundary and is worth retaining.
    if candidates:
        return candidates[0]
    return None


# ----------------------------------------------------------------------
# Main experiment
# ----------------------------------------------------------------------

print("=" * 78)
print("EXPERIMENT 506 START")
print("=" * 78)
print("2020 KERNEL: V/G CLASSIFICATION, CONIC REDUCTION, FACTOR EXPOSURE")
print("=" * 78)

# Test data used throughout the recent experiments.
TEST_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
    (3, 5),
    (5, 13),
    (13, 17),
    (3, 7),
    (7, 11),
]

print("\n[1] CLOSED FORM FOR v AND g")
print("-" * 78)

formula_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    b, v, g = kernel_values(n)

    if n % 4 == 1:
        v_expected = (3 * n + 1) // 4
        g_expected = (n + 3) // 4
    else:
        v_expected = (n + 1) // 4
        g_expected = (n + 1) // 4

    ok = (v == v_expected and g == g_expected)
    if not ok:
        formula_failures += 1

    print(
        f"  n={n:<15} n mod 4={n%4} "
        f"v={v:<15} g={g:<15} PASS={ok}"
    )

print(f"  CLOSED-FORM FAILURES = {formula_failures}")

print("\n[2] v==g CLASSIFICATION THEOREM")
print("-" * 78)

classification_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    _, v, g = kernel_values(n)
    predicted = (n % 4 == 3)
    observed = (v == g)
    ok = predicted == observed
    classification_failures += not ok
    print(
        f"  ({p},{q}) n mod 4={n%4} "
        f"v==g={observed} predicted={predicted} PASS={ok}"
    )

print(f"  CLASSIFICATION FAILURES = {classification_failures}")

print("\n[3] PRIME RESIDUE-CLASS INTERPRETATION")
print("-" * 78)

residue_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    opposite = (p % 4) != (q % 4)
    observed = classify_residue(n) == "v==g"
    ok = opposite == observed
    residue_failures += not ok
    print(
        f"  ({p},{q}) p mod4={p%4}, q mod4={q%4} "
        f"opposite={opposite} observed_v==g={observed} PASS={ok}"
    )

print(f"  RESIDUE-CLASS FAILURES = {residue_failures}")

print("\n[4] CONIC -> DIFFERENCE OF SQUARES")
print("-" * 78)

conic_failures = 0
factor_equation_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    b, v, g = kernel_values(n)

    # Historical x,y obtained directly from p,q only for validation.
    xy = historical_xy(p, q)
    if xy is None:
        conic_failures += 1
        print(f"  ({p},{q}) NO HISTORICAL x,y")
        continue

    x, y, tag = xy
    lhs = y * y - x * x + 3 * x - 2

    # Completed-square form.
    square_lhs = (2 * y) ** 2 - (2 * x - 3) ** 2
    expected_square_rhs = 4 * v - 1

    ok_conic = lhs == v and square_lhs == expected_square_rhs

    if n % 4 == 3:
        expected_product = n
    else:
        expected_product = 3 * n

    A = 2 * y - 2 * x + 3
    B = 2 * y + 2 * x - 3
    ok_factor_equation = A * B == expected_product

    conic_failures += not ok_conic
    factor_equation_failures += not ok_factor_equation

    print(
        f"  ({p},{q}) tag={tag:<8} x={x:<10} y={y:<10} "
        f"conic={ok_conic} factor-equation={ok_factor_equation}"
    )

print(f"  CONIC FAILURES = {conic_failures}")
print(f"  FACTOR-EQUATION FAILURES = {factor_equation_failures}")

print("\n[5] FACTOR RECOVERY FROM THE HISTORICAL x,y")
print("-" * 78)

recovery_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    b, v, g = kernel_values(n)
    xy = historical_xy(p, q)

    if xy is None:
        recovery_failures += 1
        continue

    x, y, tag = xy
    factors, A, B, scale = reconstructed_factors_from_xy(n, x, y, v)

    expected = sorted([p, q])
    ok = factors is not None and sorted(factors) == expected
    recovery_failures += not ok

    print(
        f"  ({p},{q}) "
        f"A={A} B={B} scale={scale} "
        f"recovered={factors} PASS={ok}"
    )

print(f"  RECOVERY FAILURES = {recovery_failures}")

print("\n[6] HISTORICAL GCD CLAIM AUDIT")
print("-" * 78)
print("Historical claim:")
print("  p = gcd(b - x - y + 2, n)")
print()
print("The algebra instead gives, in the v==g branch:")
print("  b - x - y + 2 = q*(p-1)/2")
print("so the gcd is q (for distinct odd primes), not p.")
print()

gcd_claim_failures = 0

for p, q in TEST_PAIRS:
    n = p * q
    _, v, g = kernel_values(n)

    if v != g:
        continue

    x, y, _ = historical_xy(p, q)
    value = n // 2 - x - y + 2
    got = gcd(value, n)

    claim_ok = got == p
    corrected_ok = got == q

    if not corrected_ok:
        gcd_claim_failures += 1

    print(
        f"  ({p},{q}) value={value} "
        f"gcd={got} historical_p_claim={claim_ok} "
        f"corrected_q_claim={corrected_ok}"
    )

print(f"  CORRECTED GCD FAILURES = {gcd_claim_failures}")

print("\n[7] SEARCH SIMPLE LINEAR gcd EXPRESSIONS")
print("-" * 78)
print("Search b + s_x*x + s_y*y + c for direct exposure of p or q.")
print()

candidate_forms = []
for sx in (-1, 0, 1):
    for sy in (-1, 0, 1):
        for c in range(-4, 5):
            candidate_forms.append((sx, sy, c))

hits = []

for sx, sy, c in candidate_forms:
    ok_p = True
    ok_q = True
    for p, q in TEST_PAIRS:
        n = p * q
        b, v, g = kernel_values(n)
        xy = historical_xy(p, q)
        if xy is None or v != g:
            continue
        x, y, _ = xy
        value = b + sx * x + sy * y + c
        gp = gcd(value, n)
        ok_p &= gp == p
        ok_q &= gp == q
    if ok_p:
        hits.append(("p", sx, sy, c))
    if ok_q:
        hits.append(("q", sx, sy, c))

for target, sx, sy, c in hits:
    print(
        f"  HIT {target}: gcd(b {sx:+d}x {sy:+d}y {c:+d}, n)"
    )

print(f"  TOTAL SIMPLE GCD HITS = {len(hits)}")

print("\n[8] N-ONLY CONIC REDUCTION")
print("-" * 78)
print(
    "The historical equation is already equivalent to a factor equation:"
)
print()
print("  y^2 - x^2 + 3x - 2 = v")
print("  => (2y)^2 - (2x-3)^2 = 4v-1")
print()

for residue in (1, 3):
    if residue == 3:
        print("  n == 3 (mod 4):")
        print("      v = (n+1)/4")
        print("      4v-1 = n")
        print("      (2y-2x+3)(2y+2x-3) = n")
    else:
        print("  n == 1 (mod 4):")
        print("      v = (3n+1)/4")
        print("      4v-1 = 3n")
        print("      (2y-2x+3)(2y+2x-3) = 3n")
    print()

print("\n[9] DIRECT x,y SEARCH WITHOUT p,q")
print("-" * 78)

search_failures = 0

for p, q in TEST_PAIRS[:8]:
    n = p * q
    b, v, g = kernel_values(n)

    xy = find_xy(n, v)

    if xy is None:
        search_failures += 1
        print(f"  ({p},{q}) no conic solution found in search range")
        continue

    x, y = xy
    factors, A, B, scale = reconstructed_factors_from_xy(n, x, y, v)
    ok = factors is not None and sorted(factors) == sorted([p, q])

    search_failures += not ok

    print(
        f"  ({p},{q}) x={x} y={y} A={A} B={B} "
        f"scale={scale} recovered={factors} PASS={ok}"
    )

print(f"  SOURCE-FREE CONIC SEARCH FAILURES = {search_failures}")

print("\n[10] EXHAUSTIVE SMALL PRIME-PAIR CLASSIFICATION")
print("-" * 78)

limit = 200
primes = [p for p in sieve(limit) if p != 2]

stats = {
    "pairs": 0,
    "v_eq_g": 0,
    "v_ne_g": 0,
    "conic_fail": 0,
    "factor_fail": 0,
}

for p, q in combinations(primes, 2):
    n = p * q
    stats["pairs"] += 1

    b, v, g = kernel_values(n)
    if v == g:
        stats["v_eq_g"] += 1
    else:
        stats["v_ne_g"] += 1

    xy = historical_xy(p, q)
    if xy is None:
        stats["conic_fail"] += 1
        continue

    x, y, _ = xy
    if y * y - x * x + 3 * x - 2 != v:
        stats["conic_fail"] += 1

    factors, A, B, scale = reconstructed_factors_from_xy(n, x, y, v)
    if factors is None or sorted(factors) != sorted([p, q]):
        stats["factor_fail"] += 1

print(f"  tested prime pairs = {stats['pairs']}")
print(f"  v==g pairs         = {stats['v_eq_g']}")
print(f"  v!=g pairs         = {stats['v_ne_g']}")
print(f"  conic failures     = {stats['conic_fail']}")
print(f"  factor failures    = {stats['factor_fail']}")

print("\n[11] RESEARCH INTERPRETATION")
print("-" * 78)
print("1. v==g is exactly an n mod 4 test for odd n.")
print("2. For odd prime products, it distinguishes opposite versus equal")
print("   prime residue classes modulo 4.")
print("3. The historical conic is not an independent mysterious relation:")
print("      y^2-x^2+3x-2=v")
print("   completes to")
print("      (2y-2x+3)(2y+2x-3)=n      if n==3 mod 4")
print("      (2y-2x+3)(2y+2x-3)=3n    if n==1 mod 4")
print("4. Therefore finding an integer (x,y) solution with the required")
print("   factor sign/positivity constraints is essentially equivalent")
print("   to finding a factor pair of n (or 3n).")
print("5. The older gcd formula contains a label reversal: in the")
print("   v==g branch, gcd(b-x-y+2,n)=q for the stated x,y definitions.")
print("6. This is still valuable for the current KAPPA research because")
print("   it gives a second historical coordinate system:")
print("      factor pair -> (x,y) -> difference-of-squares")
print("   that can be compared directly against the newer")
print("   S / Delta / shifted-product coordinates.")

print("\n[12] BRIDGE TO CURRENT RESEARCH")
print("-" * 78)
print("Current coordinates:")
print("  N = pq")
print("  S = p+q")
print("  Delta = S^2-4N = (p-q)^2")
print()
print("2020 coordinates:")
print("  A0 = 2y-2x+3")
print("  B0 = 2y+2x-3")
print()
print("Then:")
print("  A0*B0 = n or 3n")
print("  A0,B0 encode the factor pair.")
print()
print("This gives a concrete new research question:")
print("  Can an N-only object from the homogeneous-layer construction")
print("  reproduce the 2020 conic invariant v, or directly produce")
print("  the difference-of-squares factors A0,B0?")
print()
print("That question is structurally different from searching")
print("continued-fraction features for A=(p+q)/2.")

print("\n[13] FINAL STATUS")
print("-" * 78)
overall = all(
    stats[k] == 0
    for k in ("conic_fail", "factor_fail")
) and formula_failures == 0 and classification_failures == 0
print(f"  closed-form v/g formulas        = {formula_failures == 0}")
print(f"  v/g classification theorem      = {classification_failures == 0}")
print(f"  historical conic audit          = {stats['conic_fail'] == 0}")
print(f"  factor reconstruction           = {stats['factor_fail'] == 0}")
print(f"  corrected gcd identity          = {gcd_claim_failures == 0}")
print(f"  OVERALL VALIDATION               = {overall}")

print("\n" + "=" * 78)
print("EXPERIMENT 506 FINISHED")
print("=" * 78)
