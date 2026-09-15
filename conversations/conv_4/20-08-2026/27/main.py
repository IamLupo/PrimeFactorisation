#!/usr/bin/env python3

import math
import sympy as sp


print("EXPERIMENT 503 START")
print("=" * 78)
print("CONTINUED-FRACTION SQUARE-RESIDUE / FACTOR-EXPOSURE SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q, N = sp.symbols("p q N", integer=True)
h, k, r = sp.symbols("h k r", integer=True)


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


# ============================================================================
# PARAMETERS
# ============================================================================

MAX_CF_TERMS = 128
MAX_SEMICONVERGENT_MULTIPLIER = 32


# ============================================================================
# SAFE SYMBOLIC HELPERS
# ============================================================================

def exact_zero(expr):
    """
    Exact symbolic zero test.
    """
    return sp.simplify(expr) == 0


def exact_equal(lhs, rhs):
    """
    Exact symbolic equality test.
    """
    return sp.simplify(lhs - rhs) == 0


def is_square(n):
    """
    Return sqrt(n) if n is a non-negative perfect square,
    otherwise return None.
    """
    if n < 0:
        return None

    root = math.isqrt(n)

    if root * root == n:
        return root

    return None


# ============================================================================
# CONTINUED FRACTION OF sqrt(N)
# ============================================================================

def continued_fraction_sqrt(n, max_terms=128):
    """
    Exact continued fraction expansion of sqrt(n).

    For nonsquare n:
        sqrt(n) = [a0; a1, a2, ...]

    Returns at most max_terms entries.
    """
    a0 = math.isqrt(n)

    if a0 * a0 == n:
        return [a0]

    m = 0
    d = 1
    a = a0

    result = [a0]

    for _ in range(max_terms - 1):
        m = d * a - m

        numerator = n - m * m
        if numerator <= 0:
            break

        d = numerator // d
        a = (a0 + m) // d

        result.append(a)

        # End of one period.
        if a == 2 * a0:
            break

    return result


# ============================================================================
# CONVERGENTS
# ============================================================================

def convergents(cf):
    """
    Return exact convergents (h, k).
    """
    h_m2, h_m1 = 0, 1
    k_m2, k_m1 = 1, 0

    result = []

    for a in cf:
        h_curr = a * h_m1 + h_m2
        k_curr = a * k_m1 + k_m2

        result.append((h_curr, k_curr))

        h_m2, h_m1 = h_m1, h_curr
        k_m2, k_m1 = k_m1, k_curr

    return result


# ============================================================================
# SEMICONVERGENTS
# ============================================================================

def semiconvergents(cf, max_multiplier=32):
    """
    Generate bounded intermediate convergents.

    For each continued-fraction coefficient a_n, inspect

        (t*h_{n-1} + h_{n-2}) /
        (t*k_{n-1} + k_{n-2})

    for 1 <= t <= min(a_n, max_multiplier).
    """
    result = []

    h_m2, h_m1 = 0, 1
    k_m2, k_m1 = 1, 0

    for idx, a in enumerate(cf):
        limit = min(a, max_multiplier)

        for t in range(1, limit + 1):
            hh = t * h_m1 + h_m2
            kk = t * k_m1 + k_m2

            result.append((idx, t, hh, kk))

        h_curr = a * h_m1 + h_m2
        k_curr = a * k_m1 + k_m2

        h_m2, h_m1 = h_m1, h_curr
        k_m2, k_m1 = k_m1, k_curr

    return result


# ============================================================================
# GCD CHANNEL
# ============================================================================

def factor_gcd_channels(n, hh, rr):
    """
    Given

        hh^2 - n*kk^2 = rr^2,

    the factorization identity is

        (hh-rr)(hh+rr) = n*kk^2.

    Return gcd(h-r,n), gcd(h+r,n).
    """
    g_minus = math.gcd(abs(hh - rr), n)
    g_plus = math.gcd(abs(hh + rr), n)

    return g_minus, g_plus


def is_proper_factor(g, n):
    return 1 < g < n


# ============================================================================
# 1. SYMBOLIC SQUARE-RESIDUE FACTORIZATION
# ============================================================================

print("[1] SYMBOLIC SQUARE-RESIDUE FACTORIZATION")
print("-" * 78)

identity_1 = sp.expand(
    (h - r) * (h + r) - (h**2 - r**2)
)

print("  (h-r)(h+r) - (h^2-r^2) =")
print("   ", sp.factor(identity_1))
print("  PASS =", exact_zero(identity_1))
print()


# ============================================================================
# 2. CONDITIONAL FACTORIZATION WITH N*k^2
# ============================================================================

print("[2] CONDITIONAL FACTORIZATION")
print("-" * 78)

conditional_identity = sp.expand(
    (h - r) * (h + r)
    - (N * k**2)
)

print("  Under h^2-N*k^2=r^2:")
print("    (h-r)(h+r)=N*k^2")
print("  Residual before imposing condition =")
print("   ", conditional_identity)
print()


# ============================================================================
# 3. FERMAT REPRESENTATION
# ============================================================================

print("[3] FERMAT REPRESENTATION")
print("-" * 78)

A = (p + q) / 2

fermat_identity = sp.expand(
    A**2
    - p*q
    - ((p - q) / 2)**2
)

print("  A = (p+q)/2")
print("  A^2 - pq - (p-q)^2/4 =")
print("   ", sp.factor(fermat_identity))
print("  PASS =", exact_zero(fermat_identity))
print()


# ============================================================================
# 4. ODD-PRIME / INTEGER FERΜAT COORDINATES
# ============================================================================

print("[4] FERMAT COORDINATES")
print("-" * 78)

A_sym, B_sym = sp.symbols("A B", integer=True)

fermat_coordinate_identity = sp.expand(
    A_sym**2 - B_sym**2 - p*q
).subs(
    {
        A_sym: (p + q) / 2,
        B_sym: (q - p) / 2,
    }
)

fermat_coordinate_identity = sp.expand(
    fermat_coordinate_identity
)

print("  A=(p+q)/2")
print("  B=(q-p)/2")
print("  A^2-B^2-pq =")
print("   ", sp.factor(fermat_coordinate_identity))
print(
    "  PASS =",
    exact_zero(fermat_coordinate_identity),
)
print()


# ============================================================================
# 5. EXACT CF SQUARE-RESIDUE SEARCH
# ============================================================================

print("[5] EXACT CF SQUARE-RESIDUE SEARCH")
print("-" * 78)

all_square_hits = []

for pp, qq in INSTANCES:
    n = pp * qq

    cf = continued_fraction_sqrt(
        n,
        MAX_CF_TERMS,
    )

    conv = convergents(cf)

    square_hits = []

    print(f"  ({pp},{qq}) N={n}")

    for idx, (hh, kk) in enumerate(conv):
        residue = hh * hh - n * kk * kk
        rr = is_square(residue)

        if rr is None:
            continue

        square_hits.append(
            {
                "idx": idx,
                "h": hh,
                "k": kk,
                "residue": residue,
                "r": rr,
            }
        )

        print(
            f"    square residue: "
            f"idx={idx} h={hh} k={kk} "
            f"r={rr}"
        )

    print(
        f"    exact square-residue hits = "
        f"{len(square_hits)}"
    )

    all_square_hits.append(
        (pp, qq, square_hits)
    )

print()


# ============================================================================
# 6. GCD FACTOR-EXPOSURE SEARCH
# ============================================================================

print("[6] GCD FACTOR-EXPOSURE SEARCH")
print("-" * 78)

convergent_factor_hits = []

for pp, qq, square_hits in all_square_hits:
    n = pp * qq

    print(f"  ({pp},{qq})")

    local_hits = 0

    for hit in square_hits:
        idx = hit["idx"]
        hh = hit["h"]
        kk = hit["k"]
        rr = hit["r"]

        g_minus, g_plus = factor_gcd_channels(
            n,
            hh,
            rr,
        )

        proper_minus = is_proper_factor(g_minus, n)
        proper_plus = is_proper_factor(g_plus, n)

        print(
            f"    idx={idx}: "
            f"gcd(h-r,N)={g_minus}, "
            f"gcd(h+r,N)={g_plus}"
        )

        if proper_minus or proper_plus:
            local_hits += 1

            convergent_factor_hits.append(
                (
                    pp,
                    qq,
                    idx,
                    hh,
                    kk,
                    rr,
                    g_minus,
                    g_plus,
                )
            )

    print(
        f"    nontrivial factor hits = {local_hits}"
    )

print()


# ============================================================================
# 7. SEMICONVERGENT SEARCH
# ============================================================================

print("[7] SEMICONVERGENT SQUARE-RESIDUE SEARCH")
print("-" * 78)

semiconv_square_hits = []
semiconv_factor_hits = []

for pp, qq in INSTANCES:
    n = pp * qq

    cf = continued_fraction_sqrt(
        n,
        MAX_CF_TERMS,
    )

    semi = semiconvergents(
        cf,
        MAX_SEMICONVERGENT_MULTIPLIER,
    )

    local_square_hits = 0
    local_factor_hits = 0

    print(f"  ({pp},{qq})")

    for idx, t, hh, kk in semi:
        residue = hh * hh - n * kk * kk
        rr = is_square(residue)

        if rr is None:
            continue

        local_square_hits += 1

        semiconv_square_hits.append(
            (
                pp,
                qq,
                idx,
                t,
                hh,
                kk,
                residue,
                rr,
            )
        )

        g_minus, g_plus = factor_gcd_channels(
            n,
            hh,
            rr,
        )

        proper_minus = is_proper_factor(
            g_minus,
            n,
        )

        proper_plus = is_proper_factor(
            g_plus,
            n,
        )

        if proper_minus or proper_plus:
            local_factor_hits += 1

            semiconv_factor_hits.append(
                (
                    pp,
                    qq,
                    idx,
                    t,
                    hh,
                    kk,
                    rr,
                    g_minus,
                    g_plus,
                )
            )

    print(
        f"    exact square residues = "
        f"{local_square_hits}"
    )
    print(
        f"    nontrivial factor hits = "
        f"{local_factor_hits}"
    )

print()


# ============================================================================
# 8. FERΜAT SOLUTION VISIBILITY
# ============================================================================

print("[8] FERMAT SOLUTION VISIBILITY")
print("-" * 78)

fermat_visibility_failures = 0

for pp, qq, square_hits in all_square_hits:
    n = pp * qq

    true_A = (pp + qq) // 2
    true_B = abs(pp - qq) // 2

    found = False

    for hit in square_hits:
        if (
            hit["h"] == true_A
            and hit["k"] == 1
            and hit["r"] == true_B
        ):
            found = True
            break

    if not found:
        fermat_visibility_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"A={true_A} "
        f"B={true_B} "
        f"CF convergent A/1 = {found}"
    )

print(
    "  FERMAT VISIBILITY FAILURES =",
    fermat_visibility_failures,
)
print()


# ============================================================================
# 9. DIRECT FERMAT BASELINE
# ============================================================================

print("[9] DIRECT FERMAT BASELINE")
print("-" * 78)

fermat_failures = 0

for pp, qq in INSTANCES:
    n = pp * qq

    a = math.isqrt(n)

    if a * a < n:
        a += 1

    iterations = 0

    while True:
        diff = a * a - n
        b = math.isqrt(diff)

        if b * b == diff:
            break

        a += 1
        iterations += 1

    recovered_p = a - b
    recovered_q = a + b

    passed = (
        {recovered_p, recovered_q}
        ==
        {pp, qq}
    )

    if not passed:
        fermat_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"A={a} "
        f"B={b} "
        f"iterations={iterations} "
        f"roots=({recovered_p},{recovered_q}) "
        f"PASS={passed}"
    )

print(
    "  FERMAT FAILURES =",
    fermat_failures,
)
print()


# ============================================================================
# 10. GAP-SQUARE PROFILE
# ============================================================================

print("[10] GAP-SQUARE PROFILE")
print("-" * 78)

for pp, qq in INSTANCES:
    n = pp * qq

    A_true = (pp + qq) // 2
    B_true = abs(pp - qq) // 2
    delta_true = (pp - qq) ** 2

    floor_sqrt = math.isqrt(n)

    ceil_sqrt = (
        floor_sqrt
        if floor_sqrt * floor_sqrt == n
        else floor_sqrt + 1
    )

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    N              = {n}"
    )
    print(
        f"    floor(sqrt(N)) = {floor_sqrt}"
    )
    print(
        f"    ceil(sqrt(N))  = {ceil_sqrt}"
    )
    print(
        f"    A              = {A_true}"
    )
    print(
        f"    B              = {B_true}"
    )
    print(
        f"    Delta          = {delta_true}"
    )
    print(
        f"    A-ceil(sqrtN)  = {A_true - ceil_sqrt}"
    )

print()


# ============================================================================
# 11. SAME-N CONTROL
# ============================================================================

print("[11] SAME-N CONTROL")
print("-" * 78)

same_n_cases = [
    (12, [(2, 6), (3, 4)]),
    (18, [(2, 9), (3, 6)]),
    (20, [(4, 5), (2, 10)]),
    (30, [(5, 6), (3, 10)]),
]

same_n_failures = 0

for n, pairs in same_n_cases:
    deltas = []

    print(f"  N={n}")

    for aa, bb in pairs:
        S_val = aa + bb
        delta_val = (aa - bb) ** 2

        deltas.append(delta_val)

        print(
            f"    ({aa},{bb}) "
            f"S={S_val} "
            f"Delta={delta_val}"
        )

    passed = len(set(deltas)) == len(deltas)

    if not passed:
        same_n_failures += 1

    print(
        f"    distinct Delta values = "
        f"{len(set(deltas))} "
        f"PASS={passed}"
    )

print(
    "  SAME-N CONTROL FAILURES =",
    same_n_failures,
)
print()


# ============================================================================
# 12. EXACT GCD IDENTITY AUDIT
# ============================================================================

print("[12] EXACT GCD IDENTITY AUDIT")
print("-" * 78)

gcd_identity_failures = 0

for pp, qq, square_hits in all_square_hits:
    n = pp * qq

    for hit in square_hits:
        hh = hit["h"]
        kk = hit["k"]
        rr = hit["r"]

        lhs = (hh - rr) * (hh + rr)
        rhs = n * kk * kk

        if lhs != rhs:
            gcd_identity_failures += 1

print(
    "  (h-r)(h+r)=N*k^2 failures =",
    gcd_identity_failures,
)
print()


# ============================================================================
# 13. FACTOR-HIT DETAIL
# ============================================================================

print("[13] NONTRIVIAL FACTOR-HIT DETAIL")
print("-" * 78)

if not convergent_factor_hits:
    print(
        "  No nontrivial factor exposures from exact-square "
        "convergent residues."
    )
else:
    for (
        pp,
        qq,
        idx,
        hh,
        kk,
        rr,
        gm,
        gp,
    ) in convergent_factor_hits:
        print(
            f"  ({pp},{qq}) "
            f"idx={idx} "
            f"h={hh} "
            f"k={kk} "
            f"r={rr} "
            f"gcd(h-r,N)={gm} "
            f"gcd(h+r,N)={gp}"
        )

print()

if not semiconv_factor_hits:
    print(
        "  No nontrivial factor exposures from exact-square "
        "semiconvergent residues."
    )
else:
    for (
        pp,
        qq,
        idx,
        t,
        hh,
        kk,
        rr,
        gm,
        gp,
    ) in semiconv_factor_hits:
        print(
            f"  ({pp},{qq}) "
            f"idx={idx} t={t} "
            f"h={hh} k={kk} r={rr} "
            f"gcd(h-r,N)={gm} "
            f"gcd(h+r,N)={gp}"
        )

print()


# ============================================================================
# 14. FINAL SYMBOLIC CERTIFICATES
# ============================================================================

print("[14] FINAL SYMBOLIC CERTIFICATES")
print("-" * 78)

gap_identity = sp.expand(
    (p + q) ** 2
    - 4 * p * q
    - (p - q) ** 2
)

fermat_identity_final = sp.expand(
    ((p + q) / 2) ** 2
    - p * q
    - ((p - q) / 2) ** 2
)

coordinate_identity = sp.expand(
    ((p + q) / 2) ** 2
    - ((q - p) / 2) ** 2
    - p * q
)

print(
    "  gap-square identity PASS =",
    exact_zero(gap_identity),
)

print(
    "  Fermat identity PASS =",
    exact_zero(fermat_identity_final),
)

print(
    "  A^2-B^2=N identity PASS =",
    exact_zero(coordinate_identity),
)

print(
    "  gcd factorization PASS =",
    gcd_identity_failures == 0,
)

print()


# ============================================================================
# 15. FINAL STATUS
# ============================================================================

total_cf_square_hits = sum(
    len(hits)
    for _, _, hits in all_square_hits
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  fundamental gap-square identity =",
    exact_zero(gap_identity),
)

print(
    "  Fermat representation =",
    exact_zero(fermat_identity_final),
)

print(
    "  A^2-B^2=N representation =",
    exact_zero(coordinate_identity),
)

print(
    "  exact CF square-residue hits =",
    total_cf_square_hits,
)

print(
    "  convergent factor-exposure hits =",
    len(convergent_factor_hits),
)

print(
    "  semiconvergent square-residue hits =",
    len(semiconv_square_hits),
)

print(
    "  semiconvergent factor-exposure hits =",
    len(semiconv_factor_hits),
)

print(
    "  Fermat visibility failures =",
    fermat_visibility_failures,
)

print(
    "  same-N control failures =",
    same_n_failures,
)

print(
    "  Fermat baseline failures =",
    fermat_failures,
)

print()


# ============================================================================
# 16. INTERPRETATION
# ============================================================================

print("MAIN RESULT")
print("-" * 78)

print(
    "The exact upstream Fermat identity is:"
)
print()
print(
    "    A=(p+q)/2"
)
print(
    "    A^2-N=((p-q)/2)^2"
)
print()

print(
    "An exact continued-fraction square residue"
)
print()
print(
    "    h^2-N*k^2=r^2"
)
print()
print(
    "implies:"
)
print()
print(
    "    (h-r)(h+r)=N*k^2."
)
print()

print(
    "Therefore the strongest useful test is not whether"
)
print(
    "a Pell residue merely becomes small."
)
print(
    "It is whether the residue is an exact square and"
)
print(
    "whether gcd(h-r,N) or gcd(h+r,N) yields a proper"
)
print(
    "nontrivial divisor of N."
)
print()

if convergent_factor_hits or semiconv_factor_hits:
    print(
        "A genuine factor-exposure event was observed."
    )
    print(
        "This deserves a dedicated follow-up experiment"
    )
    print(
        "to determine whether its index and residue can"
    )
    print(
        "be predicted directly from N."
    )
else:
    print(
        "No nontrivial factor exposure was found in the"
    )
    print(
        "tested exact-square convergent/semiconvergent family."
    )
    print(
        "The tested Pell-style channel therefore does not,"
    )
    print(
        "by itself, provide a factor leak."
    )

print()
print(
    "The central unresolved upstream target remains:"
)
print()
print(
    "    N"
)
print(
    "     -> canonical A=(p+q)/2"
)
print(
    "     -> Delta=4(A^2-N)"
)
print(
    "     -> S=2A"
)
print(
    "     -> z^2-Sz+N"
)
print(
    "     -> p,q"
)
print()

print(
    "Do not interpret the existence of many small Pell residues"
)
print(
    "as evidence of a new factorization mechanism unless they"
)
print(
    "survive the exact-square and proper-gcd tests."
)
print()

print("NEXT TARGET")
print("-" * 78)

print(
    "If this experiment produces no factor hits, the productive"
)
print(
    "next step is to stop enumerating Pell residues and instead"
)
print(
    "search for an exact N-only observable whose value is the"
)
print(
    "Fermat midpoint A or the gap-square Delta."
)

print()
print("=" * 78)
print("EXPERIMENT 503 FINISHED")
print("=" * 78)
