#!/usr/bin/env python3

import math
from fractions import Fraction

import sympy as sp


print("EXPERIMENT 502 START")
print("=" * 78)
print("N-ONLY GAP-SQUARE / FERMAT / CONTINUED-FRACTION SIGNATURE SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLIC VARIABLES
# ============================================================================

p, q, N, S, Delta, a, k = sp.symbols(
    "p q N S Delta a k",
    integer=True
)


def fact(expr):
    return sp.factor(sp.expand(expr))


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


# ============================================================================
# TEST INSTANCES
# ============================================================================

INSTANCES = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]


# Additional same-N factor-pair controls.
SAME_N_CONTROLS = [
    (2, 6),
    (3, 4),
    (2, 9),
    (3, 6),
    (4, 5),
    (2, 10),
    (5, 6),
    (3, 10),
]


# ============================================================================
# 1. FUNDAMENTAL GAP IDENTITY
# ============================================================================

print("[1] FUNDAMENTAL GAP-SQUARE IDENTITY")
print("-" * 78)

S_expr = p + q
Delta_expr = (p - q) ** 2

identity = sp.expand(
    S_expr**2 - 4*p*q - Delta_expr
)

print("  (p+q)^2 - 4pq - (p-q)^2 =")
print("   ", fact(identity))

print("  PASS =", zero(identity))
print()


# ============================================================================
# 2. FERMAT REPRESENTATION
# ============================================================================

print("[2] FERMAT REPRESENTATION")
print("-" * 78)

a_expr = (p + q) / 2

fermat_identity = sp.expand(
    a_expr**2 - p*q - (p-q)**2 / 4
)

print("  a = (p+q)/2")
print()
print("  a^2 - N =")
print("   ", fact(a_expr**2 - p*q))
print()
print("  expected Delta/4 =")
print("   ", fact(Delta_expr / 4))
print()
print("  difference =")
print("   ", fact(fermat_identity))
print()
print("  PASS =", zero(fermat_identity))
print()


# ============================================================================
# 3. ODD PRIME SPECIALIZATION
# ============================================================================

print("[3] ODD-PRIME SPECIALIZATION")
print("-" * 78)

print(
    "For odd p,q:"
)
print()
print(
    "  A = (p+q)/2 is an integer"
)
print(
    "  B = (q-p)/2 is an integer"
)
print()
print(
    "  A^2 - B^2 = pq = N"
)
print()

A, B = sp.symbols("A B")

difference_of_squares = sp.expand(
    A**2 - B**2 - N
)

substitution_difference = sp.expand(
    difference_of_squares.subs({
        A: (p + q) / 2,
        B: (q - p) / 2,
    })
)

print(
    "  Symbolic certificate:",
    fact(substitution_difference),
)

print(
    "  PASS =",
    zero(substitution_difference),
)
print()


# ============================================================================
# 4. DIRECT N-ONLY FERMA T-SCAN
# ============================================================================

print("[4] DIRECT N-ONLY FERMAT SQUARE-DISTANCE SEARCH")
print("-" * 78)

print(
    "For each N, search a = ceil(sqrt(N)), ceil(sqrt(N))+1, ..."
)
print(
    "until a^2-N is a perfect square."
)
print()

print(
    "This does not use p or q during the search."
)
print()

fermat_results = []

for pp, qq in INSTANCES:
    nn = pp * qq

    start = math.isqrt(nn)

    if start * start < nn:
        start += 1

    found = None

    # Limit the direct scan so this experiment remains finite.
    # The limit is intentionally large enough for close factors.
    max_steps = 2_000_000

    for step in range(max_steps):
        aa = start + step
        bb2 = aa * aa - nn
        bb = math.isqrt(bb2)

        if bb * bb == bb2:
            found = (aa, bb, step)
            break

    true_A = (pp + qq) // 2
    true_B = abs(pp - qq) // 2

    if found is None:
        print(
            f"  ({pp},{qq}) "
            f"NO SQUARE FOUND within {max_steps} steps"
        )
        continue

    aa, bb, step = found

    passed = (
        aa == true_A
        and bb == true_B
    )

    fermat_results.append(
        (pp, qq, aa, bb, step, passed)
    )

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    N = {nn}"
    )
    print(
        f"    Fermat A = {aa}"
    )
    print(
        f"    true A   = {true_A}"
    )
    print(
        f"    Fermat B = {bb}"
    )
    print(
        f"    true B   = {true_B}"
    )
    print(
        f"    iterations from ceil(sqrt(N)) = {step}"
    )
    print(
        f"    PASS = {passed}"
    )

print()


# ============================================================================
# 5. GAP-SQUARE RECOVERY FROM FERMAT SEARCH
# ============================================================================

print("[5] GAP-SQUARE RECOVERY")
print("-" * 78)

gap_failures = 0

for pp, qq, aa, bb, step, _ in fermat_results:
    true_delta = (pp - qq) ** 2
    recovered_delta = 4 * bb * bb

    passed = recovered_delta == true_delta

    if not passed:
        gap_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"recovered Delta={recovered_delta} "
        f"expected={true_delta} "
        f"PASS={passed}"
    )

print()
print(
    "  GAP-SQUARE FAILURES =",
    gap_failures,
)
print()


# ============================================================================
# 6. N-ONLY QUANTITY: DISTANCE TO NEXT SQUARE
# ============================================================================

print("[6] DISTANCE-TO-NEXT-SQUARE SIGNATURE")
print("-" * 78)

for pp, qq in INSTANCES:
    nn = pp * qq

    r = math.isqrt(nn)

    if r * r < nn:
        r += 1

    distance = r * r - nn

    true_delta = (pp - qq) ** 2

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    ceil(sqrt(N)) = {r}"
    )
    print(
        f"    r^2-N         = {distance}"
    )
    print(
        f"    Delta         = {true_delta}"
    )

    if distance == 0:
        relation = "N is already a square"
    elif distance == true_delta // 4:
        relation = "EXACT Delta/4 MATCH"
    else:
        relation = "no direct Delta/4 match"

    print(
        f"    relation      = {relation}"
    )

print()


# ============================================================================
# 7. SAME-N FACTOR-PAIR CONTROL
# ============================================================================

print("[7] SAME-N FACTOR-PAIR CONTROL")
print("-" * 78)

grouped = {}

for pp, qq in SAME_N_CONTROLS:
    nn = pp * qq
    ss = pp + qq
    dd = abs(pp - qq)

    grouped.setdefault(nn, []).append(
        (pp, qq, ss, dd * dd)
    )

same_n_failures = 0

for nn, entries in grouped.items():
    if len(entries) < 2:
        continue

    print(f"  N = {nn}")

    values = set(
        entry[3]
        for entry in entries
    )

    for entry in entries:
        pp, qq, ss, dd2 = entry

        print(
            f"    ({pp},{qq}) "
            f"S={ss} "
            f"Delta={dd2}"
        )

    passed = len(values) > 1

    print(
        f"    distinct Delta values = {len(values)}"
    )
    print(
        f"    PASS = {passed}"
    )
    print()

    if not passed:
        same_n_failures += 1


# ============================================================================
# 8. CONTINUED FRACTION OF sqrt(N)
# ============================================================================

print("[8] CONTINUED-FRACTION ANALYSIS OF sqrt(N)")
print("-" * 78)

def continued_fraction_sqrt(n, max_terms=64):
    a0 = math.isqrt(n)

    if a0 * a0 == n:
        return [a0]

    m = 0
    d = 1
    a = a0

    period = []

    for _ in range(max_terms - 1):
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d

        period.append(a)

        if a == 2 * a0:
            break

    return [a0] + period


def convergents(cf):
    p0, p1 = 0, 1
    q0, q1 = 1, 0

    result = []

    for a in cf:
        p2 = a * p1 + p0
        q2 = a * q1 + q0

        result.append((p2, q2))

        p0, p1 = p1, p2
        q0, q1 = q1, q2

    return result


for pp, qq in INSTANCES:
    nn = pp * qq

    cf = continued_fraction_sqrt(nn)
    conv = convergents(cf)

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    sqrt(N) CF = {cf}"
    )

    # Search convergents for a Fermat-like relation.
    relation_hits = []

    for idx, (num, den) in enumerate(conv):
        # Difference between num^2 and N*den^2.
        residue = num * num - nn * den * den

        if residue != 0:
            abs_residue = abs(residue)
            rr = math.isqrt(abs_residue)

            if rr * rr == abs_residue:
                relation_hits.append(
                    (idx, num, den, residue, rr)
                )

    if relation_hits:
        for hit in relation_hits:
            print(
                "    square-residue convergent:",
                hit
            )
    else:
        print(
            "    no square-residue convergent found"
        )

print()


# ============================================================================
# 9. PELL-TYPE RESIDUE SEARCH
# ============================================================================

print("[9] PELL-TYPE RESIDUE SEARCH")
print("-" * 78)

print(
    "For convergents h/k of sqrt(N), inspect"
)
print(
    "    h^2 - N k^2."
)
print(
    "A small square residue would indicate a possible"
)
print(
    "Pell-like route to the factor gap."
)
print()

pell_hits = []

for pp, qq in INSTANCES:
    nn = pp * qq

    cf = continued_fraction_sqrt(nn, max_terms=128)
    conv = convergents(cf)

    local = []

    for idx, (h, kk) in enumerate(conv):
        residue = h*h - nn*kk*kk

        if residue == 0:
            continue

        absolute = abs(residue)

        # Search only small residues relative to N.
        if absolute <= 4 * math.isqrt(nn):
            local.append(
                (idx, h, kk, residue)
            )

    if local:
        print(
            f"  ({pp},{qq})"
        )
        for hit in local:
            print(
                f"    idx={hit[0]} "
                f"h={hit[1]} "
                f"k={hit[2]} "
                f"residue={hit[3]}"
            )

        pell_hits.extend(
            [(pp, qq, h) for h in local]
        )
    else:
        print(
            f"  ({pp},{qq}) no small Pell residues"
        )

print()


# ============================================================================
# 10. TARGETED FERMA GAP FORMULA
# ============================================================================

print("[10] TARGETED FERMA GAP FORMULA")
print("-" * 78)

print(
    "If A=(p+q)/2 is found from N alone, then:"
)
print()
print(
    "  Delta = 4*(A^2-N)"
)
print()
print(
    "  S = 2A"
)
print()
print(
    "  p,q = A +/- sqrt(A^2-N)"
)
print()

fermat_formula_pass = True

for pp, qq in INSTANCES:
    nn = pp * qq
    AA = (pp + qq) // 2

    recovered_delta = 4 * (AA * AA - nn)
    recovered_s = 2 * AA

    expected_delta = (pp - qq) ** 2
    expected_s = pp + qq

    local_pass = (
        recovered_delta == expected_delta
        and recovered_s == expected_s
    )

    if not local_pass:
        fermat_formula_pass = False

    print(
        f"  ({pp},{qq}) "
        f"A={AA} "
        f"S={recovered_s} "
        f"Delta={recovered_delta} "
        f"PASS={local_pass}"
    )

print()
print(
    "  FERMAT FORMULA PASS =",
    fermat_formula_pass,
)
print()


# ============================================================================
# 11. WHAT AN N-ONLY BRIDGE WOULD HAVE TO DO
# ============================================================================

print("[11] INFORMATION-THEORETIC TARGET")
print("-" * 78)

print(
    "An N-only bridge cannot use an arbitrary universal"
)
print(
    "function f(N)=S on the unrestricted factor space,"
)
print(
    "because different factor pairs can share the same N."
)
print()

print(
    "Therefore the useful target is not"
)
print(
    "    S = f(N)"
)
print()
print(
    "but rather an N-derived canonical object A(N)"
)
print(
    "satisfying"
)
print()
print(
    "    A(N) = (p+q)/2"
)
print()
print(
    "for the selected factorization class."
)
print()

print(
    "Then automatically:"
)
print(
    "    Delta = 4(A(N)^2-N)"
)
print(
    "    S     = 2A(N)"
)
print(
    "    p,q   = A(N) +/- sqrt(A(N)^2-N)"
)
print()


# ============================================================================
# 12. FINAL STATUS
# ============================================================================

print("[12] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  fundamental gap identity:",
    zero(identity),
)

print(
    "  Fermat representation:",
    zero(fermat_identity),
)

print(
    "  numerical Fermat searches:",
    len(fermat_results),
    "completed",
)

print(
    "  gap-square recovery failures:",
    gap_failures,
)

print(
    "  same-N control failures:",
    same_n_failures,
)

print(
    "  continued-fraction instances:",
    len(INSTANCES),
)

print(
    "  Pell-like hits:",
    len(pell_hits),
)

print(
    "  Fermat formula certificate:",
    fermat_formula_pass,
)

print()
print("MAIN RESULT")
print("-" * 78)

print(
    "The factorization problem can be rewritten as a square-distance"
)
print(
    "problem around sqrt(N):"
)
print()
print(
    "    A^2-N = Delta/4"
)
print()
print(
    "where A=(p+q)/2."
)
print()

print(
    "Thus a successful N-only construction of A would immediately"
)
print(
    "produce both the symmetric sum and the factor gap:"
)
print()
print(
    "    A(N)"
)
print(
    "      -> Delta = 4(A(N)^2-N)"
)
print(
    "      -> S = 2A(N)"
)
print(
    "      -> z^2-Sz+N"
)
print(
    "      -> p,q"
)
print()

print(
    "This gives a concrete upstream target that is independent of"
)
print(
    "the shifted KAPPA kernel."
)
print()

print("NEXT TARGET")
print("-" * 78)

print(
    "If the continued-fraction/Pell/Fermat channels show a structured"
)
print(
    "relationship with the gap square, Experiment 503 should search"
)
print(
    "for an exact map from the homogeneous-layer observable to"
)
print(
    "the Fermat quantity A=(p+q)/2."
)
print()

print(
    "If they do not, the next step should inspect the existing"
)
print(
    "N-only homogeneous-layer values directly for square-distance"
)
print(
    "or discriminant signatures, rather than generating more"
)
print(
    "shifted-kernel identities."
)

print()
print("=" * 78)
print("EXPERIMENT 502 FINISHED")
print("=" * 78)
