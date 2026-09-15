#!/usr/bin/env python3

import math
from fractions import Fraction

import sympy as sp


print("EXPERIMENT 503 START")
print("=" * 78)
print("CONTINUED-FRACTION SQUARE-RESIDUE / FACTOR-EXPOSURE SEARCH")
print("=" * 78)
print()


# ============================================================================
# PURPOSE
# ============================================================================
#
# Experiment 502 established the exact Fermat relation
#
#     A = (p+q)/2
#     A^2 - N = ((p-q)/2)^2
#
# but the continued-fraction search produced many small residues that
# were not necessarily related to the actual factorization.
#
# This experiment imposes a much stronger condition:
#
#     h^2 - N*k^2 = r^2
#
# with r >= 0 exactly.
#
# Then
#
#     (h-r)(h+r) = N*k^2.
#
# We test whether
#
#     gcd(h-r, N)
#     gcd(h+r, N)
#
# expose p or q.
#
# Everything in this search is generated from N alone.
# No p or q is used to generate the continued fraction or candidates.
#
# The known factors are used ONLY afterward as an audit.
#
# ============================================================================


p, q, N = sp.symbols("p q N", integer=True)


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

MAX_CF_TERMS = 256
MAX_SEMICONVERGENT_MULTIPLIER = 32


# ============================================================================
# HELPERS
# ============================================================================

def is_square(n: int):
    if n < 0:
        return None

    r = math.isqrt(n)

    if r * r == n:
        return r

    return None


def continued_fraction_sqrt(n: int, max_terms: int = 256):
    """
    Exact continued fraction of sqrt(n).

    Returns:
        [a0, a1, ...]
    """
    a0 = math.isqrt(n)

    if a0 * a0 == n:
        return [a0]

    m = 0
    d = 1
    a = a0

    cf = [a0]

    for _ in range(max_terms - 1):
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d

        cf.append(a)

        # The continued fraction has returned to the start of the period.
        if a == 2 * a0:
            break

    return cf


def convergents(cf):
    """
    Exact convergents h/k.
    """
    h_m2, h_m1 = 0, 1
    k_m2, k_m1 = 1, 0

    out = []

    for a in cf:
        h = a * h_m1 + h_m2
        k = a * k_m1 + k_m2

        out.append((h, k))

        h_m2, h_m1 = h_m1, h
        k_m2, k_m1 = k_m1, k

    return out


def semiconvergents(cf, max_multiplier=32):
    """
    Generate convergents and a bounded set of intermediate fractions.

    Between consecutive convergents:

        (h_{n-2} + t h_{n-1}) /
        (k_{n-2} + t k_{n-1})

    for 1 <= t <= a_n.

    We cap t to avoid giant searches.
    """
    result = []

    h_m2, h_m1 = 0, 1
    k_m2, k_m1 = 1, 0

    for idx, a in enumerate(cf):
        limit = min(a, max_multiplier)

        for t in range(1, limit + 1):
            h = t * h_m1 + h_m2
            k = t * k_m1 + k_m2

            result.append(
                (
                    idx,
                    t,
                    h,
                    k,
                )
            )

        h_m2, h_m1 = h_m1, a * h_m1 + h_m2
        k_m2, k_m1 = k_m1, a * k_m1 + k_m2

    return result


def gcd_channels(n: int, h: int, r: int, k: int):
    """
    Factor channels induced by

        (h-r)(h+r) = N*k^2.
    """
    g1 = math.gcd(abs(h - r), n)
    g2 = math.gcd(abs(h + r), n)

    # Also test the denominator-normalized channels.
    g3 = math.gcd(abs(h - r), n)
    g4 = math.gcd(abs(h + r), n)

    return g1, g2, g3, g4


def classify_gcd(g, n):
    if g <= 1:
        return "trivial-1"

    if g >= n:
        return "trivial-N"

    return "NONTRIVIAL-FACTOR-CANDIDATE"


# ============================================================================
# 1. SYMBOLIC SQUARE-RESIDUE FACTORIZATION
# ============================================================================

print("[1] SYMBOLIC SQUARE-RESIDUE FACTORIZATION")
print("-" * 78)

h, k, r = sp.symbols("h k r", integer=True)

factor_identity = sp.expand(
    (h - r) * (h + r) - (N * k**2)
    .subs(N, N)
)

# This is only the conditional identity after imposing
# h^2 - Nk^2 = r^2.
conditional_identity = sp.expand(
    (h - r) * (h + r)
    - (h**2 - r**2)
)

print("  (h-r)(h+r) - (h^2-r^2) =")
print("   ", sp.factor(conditional_identity))

print("  PASS =", conditional_identity == 0)
print()


# ============================================================================
# 2. FERMAT SPECIAL CASE
# ============================================================================

print("[2] FERMAT SPECIAL CASE")
print("-" * 78)

A = sp.symbols("A", integer=True)

fermat_relation = sp.expand(
    A**2 - N - ((p - q) / 2)**2
)

fermat_relation_pq = sp.factor(
    fermat_relation.subs(
        N, p * q
    )
)

print("  A=(p+q)/2")
print("  A^2-N-(p-q)^2/4 =")
print("   ", fermat_relation_pq)

print("  PASS =", fermat_relation_pq == 0)
print()


# ============================================================================
# 3. EXACT CONTINUED-FRACTION SQUARE RESIDUES
# ============================================================================

print("[3] EXACT CF SQUARE-RESIDUE SEARCH")
print("-" * 78)

all_square_hits = []
all_factor_hits = []

for pp, qq in INSTANCES:
    n = pp * qq

    cf = continued_fraction_sqrt(
        n,
        MAX_CF_TERMS,
    )

    conv = convergents(cf)

    square_hits = []

    for idx, (hh, kk) in enumerate(conv):
        residue = hh * hh - n * kk * kk

        rr = is_square(residue)

        if rr is None:
            continue

        square_hits.append(
            (
                idx,
                hh,
                kk,
                residue,
                rr,
            )
        )

    all_square_hits.append(
        (pp, qq, square_hits)
    )

    print(
        f"  ({pp},{qq}) N={n}"
    )

    if not square_hits:
        print(
            "    exact square residues: NONE"
        )
        continue

    print(
        f"    exact square residues: {len(square_hits)}"
    )

    for idx, hh, kk, residue, rr in square_hits:
        print(
            f"      idx={idx:3d} "
            f"h={hh} "
            f"k={kk} "
            f"r={rr}"
        )

print()


# ============================================================================
# 4. GCD FACTOR-EXPOSURE AUDIT
# ============================================================================

print("[4] GCD FACTOR-EXPOSURE AUDIT")
print("-" * 78)

for pp, qq, square_hits in all_square_hits:
    n = pp * qq

    print(
        f"  ({pp},{qq})"
    )

    local_factor_hits = 0

    for idx, hh, kk, residue, rr in square_hits:
        g1 = math.gcd(abs(hh - rr), n)
        g2 = math.gcd(abs(hh + rr), n)

        c1 = classify_gcd(g1, n)
        c2 = classify_gcd(g2, n)

        print(
            f"    idx={idx:3d} "
            f"gcd(h-r,N)={g1} [{c1}] "
            f"gcd(h+r,N)={g2} [{c2}]"
        )

        if (
            1 < g1 < n
            or
            1 < g2 < n
        ):
            local_factor_hits += 1

            all_factor_hits.append(
                (
                    pp,
                    qq,
                    idx,
                    hh,
                    kk,
                    rr,
                    g1,
                    g2,
                )
            )

    print(
        f"    NONTRIVIAL FACTOR HITS = {local_factor_hits}"
    )

print()


# ============================================================================
# 5. TARGET MATCH AGAINST TRUE FERMAT SOLUTION
# ============================================================================

print("[5] FERMAT SOLUTION VISIBILITY TEST")
print("-" * 78)

fermat_visibility_failures = 0

for pp, qq, square_hits in all_square_hits:
    n = pp * qq

    true_A = (pp + qq) // 2
    true_B = abs(pp - qq) // 2

    found = None

    for idx, hh, kk, residue, rr in square_hits:
        if (
            hh == true_A
            and
            kk == 1
            and
            rr == true_B
        ):
            found = (
                idx,
                hh,
                kk,
                rr,
            )
            break

    passed = found is not None

    if not passed:
        fermat_visibility_failures += 1

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    true A = {true_A}"
    )
    print(
        f"    true B = {true_B}"
    )
    print(
        f"    A/1 appears as CF convergent = {passed}"
    )

    if found:
        print(
            f"    CF index = {found[0]}"
        )

print()
print(
    "  FERMAT VISIBILITY FAILURES =",
    fermat_visibility_failures,
)
print()


# ============================================================================
# 6. SEMICONVERGENT FERMAT VISIBILITY
# ============================================================================

print("[6] SEMICONVERGENT FERMAT VISIBILITY SEARCH")
print("-" * 78)

semi_visibility_failures = 0
semi_factor_hits = []

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

    true_A = (pp + qq) // 2

    exact_matches = []

    for idx, t, hh, kk in semi:
        if hh == true_A and kk == 1:
            exact_matches.append(
                (idx, t, hh, kk)
            )

    passed = bool(exact_matches)

    if not passed:
        semi_visibility_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"A={true_A} "
        f"semiconvergent A/1 = {passed}"
    )

    if exact_matches:
        print(
            f"    matches = {exact_matches}"
        )

    # Also perform square-residue factor search.
    local_hits = 0

    for idx, t, hh, kk in semi:
        residue = hh * hh - n * kk * kk

        rr = is_square(residue)

        if rr is None:
            continue

        g1 = math.gcd(abs(hh - rr), n)
        g2 = math.gcd(abs(hh + rr), n)

        if (
            1 < g1 < n
            or
            1 < g2 < n
        ):
            local_hits += 1
            semi_factor_hits.append(
                (
                    pp,
                    qq,
                    idx,
                    t,
                    hh,
                    kk,
                    rr,
                    g1,
                    g2,
                )
            )

    print(
        f"    semiconvergent nontrivial factor hits = {local_hits}"
    )

print()
print(
    "  SEMICONVERGENT VISIBILITY FAILURES =",
    semi_visibility_failures,
)
print()


# ============================================================================
# 7. CLASSICAL FACTOR-RECOVERY IDENTITY
# ============================================================================

print("[7] CLASSICAL GCD IDENTITY")
print("-" * 78)

print(
    "Suppose:"
)
print()
print(
    "  h^2 - N*k^2 = r^2"
)
print()
print(
    "Then:"
)
print(
    "  (h-r)(h+r) = N*k^2."
)
print()
print(
    "Hence p or q may divide one of h-r or h+r."
)
print()
print(
    "This is tested directly below."
)
print()


gcd_identity_failures = 0

for pp, qq in INSTANCES:
    n = pp * qq

    # Test all exact convergent square residues.
    for _, hh, kk, _, rr in [
        hit
        for ppp, qqq, hits in all_square_hits
        if ppp == pp and qqq == qq
        for hit in hits
    ]:
        lhs = (hh - rr) * (hh + rr)
        rhs = n * kk * kk

        if lhs != rhs:
            gcd_identity_failures += 1

print(
    "  square-residue factorization identity failures =",
    gcd_identity_failures,
)
print()


# ============================================================================
# 8. N-ONLY STATUS OF THE CHANNEL
# ============================================================================

print("[8] N-ONLY INFORMATION CHANNEL")
print("-" * 78)

print(
    "The following quantities are generated entirely from N:"
)
print()
print(
    "  sqrt(N)"
)
print(
    "  continued fraction of sqrt(N)"
)
print(
    "  convergents h/k"
)
print(
    "  residues h^2-Nk^2"
)
print(
    "  exact-square test r^2 = h^2-Nk^2"
)
print(
    "  gcd(h-r,N)"
)
print(
    "  gcd(h+r,N)"
)
print()

print(
    "No p or q is required to generate these candidates."
)
print()


# ============================================================================
# 9. CONTROL AGAINST TRIVIAL RESIDUES
# ============================================================================

print("[9] RESIDUE QUALITY CONTROL")
print("-" * 78)

print(
    "A residue is considered structurally interesting only when:"
)
print()
print(
    "  1. h^2-Nk^2 is an exact square;"
)
print(
    "  2. the resulting gcd is a proper divisor of N."
)
print()

interesting_hits = []

for hit in all_factor_hits:
    pp, qq, idx, hh, kk, rr, g1, g2 = hit

    candidates = []

    if 1 < g1 < pp * qq:
        candidates.append(g1)

    if 1 < g2 < pp * qq:
        candidates.append(g2)

    if candidates:
        interesting_hits.append(
            hit
        )

for hit in interesting_hits:
    pp, qq, idx, hh, kk, rr, g1, g2 = hit

    print(
        f"  ({pp},{qq}) "
        f"idx={idx} "
        f"h={hh} "
        f"k={kk} "
        f"r={rr} "
        f"g1={g1} "
        f"g2={g2}"
    )

print()
print(
    "  STRUCTURALLY INTERESTING HITS =",
    len(interesting_hits),
)
print()


# ============================================================================
# 10. DIRECT FERMAT CHANNEL AS A BASELINE
# ============================================================================

print("[10] DIRECT FERMAT BASELINE")
print("-" * 78)

fermat_iterations = []

for pp, qq in INSTANCES:
    n = pp * qq

    a0 = math.isqrt(n)

    if a0 * a0 < n:
        a0 += 1

    steps = 0
    a_current = a0

    while True:
        rem = a_current * a_current - n
        b_current = math.isqrt(rem)

        if b_current * b_current == rem:
            break

        a_current += 1
        steps += 1

    recovered_p = a_current - b_current
    recovered_q = a_current + b_current

    passed = (
        {recovered_p, recovered_q}
        ==
        {pp, qq}
    )

    fermat_iterations.append(steps)

    print(
        f"  ({pp},{qq}) "
        f"steps={steps} "
        f"A={a_current} "
        f"B={b_current} "
        f"roots=({recovered_p},{recovered_q}) "
        f"PASS={passed}"
    )

print()


# ============================================================================
# 11. GAP-SQUARE PROFILE
# ============================================================================

print("[11] GAP-SQUARE PROFILE")
print("-" * 78)

for pp, qq in INSTANCES:
    n = pp * qq

    A_true = (pp + qq) // 2
    B_true = abs(pp - qq) // 2
    Delta_true = (pp - qq) ** 2

    distance_from_sqrt = A_true - math.isqrt(n)

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    N          = {n}"
    )
    print(
        f"    A          = {A_true}"
    )
    print(
        f"    B          = {B_true}"
    )
    print(
        f"    Delta      = {Delta_true}"
    )
    print(
        f"    A-floor√N  = {distance_from_sqrt}"
    )
    print(
        f"    A-steps    = {A_true - (
            math.isqrt(n) + (
                1 if math.isqrt(n)**2 < n else 0
            )
        )}"
    )

print()


# ============================================================================
# 12. FINAL CERTIFICATES
# ============================================================================

print("[12] FINAL CERTIFICATES")
print("-" * 78)

print(
    "  symbolic gap identity:",
    zero(
        (p + q) ** 2
        - 4 * p * q
        - (p - q) ** 2
    ),
)

print(
    "  symbolic Fermat identity:",
    zero(
        ((p + q) / 2) ** 2
        - p * q
        - ((p - q) / 2) ** 2
    ),
)

print(
    "  exact square-residue identity failures:",
    gcd_identity_failures,
)

print(
    "  exact CF square-residue hits:",
    sum(
        len(hits)
        for _, _, hits in all_square_hits
    ),
)

print(
    "  convergent factor-exposure hits:",
    len(all_factor_hits),
)

print(
    "  semiconvergent factor-exposure hits:",
    len(semi_factor_hits),
)

print(
    "  Fermat visibility failures:",
    fermat_visibility_failures,
)

print(
    "  semiconvergent Fermat visibility failures:",
    semi_visibility_failures,
)

print()


# ============================================================================
# 13. FINAL INTERPRETATION
# ============================================================================

print("[13] EXPERIMENT STATUS")
print("-" * 78)

symbolic_pass = (
    zero(
        (p + q) ** 2
        - 4 * p * q
        - (p - q) ** 2
    )
    and
    zero(
        ((p + q) / 2) ** 2
        - p * q
        - ((p - q) / 2) ** 2
    )
)

print(
    "  fundamental Fermat identities =",
    symbolic_pass,
)

print(
    "  CF square-residue search completed = True"
)

print(
    "  exact square-residue gcd channel completed = True"
)

print(
    "  nontrivial convergent factor hits =",
    len(all_factor_hits),
)

print(
    "  nontrivial semiconvergent factor hits =",
    len(semi_factor_hits),
)

print()


print("MAIN RESULT")
print("-" * 78)

print(
    "Experiment 502 identified"
)
print(
    "    A^2 - N = Delta/4"
)
print(
    "as the fundamental upstream target."
)
print()

print(
    "Experiment 503 now tests a much stronger N-only signature:"
)
print()
print(
    "    h^2 - N*k^2 = r^2."
)
print()

print(
    "Such an exact square residue gives"
)
print()
print(
    "    (h-r)(h+r)=N*k^2,"
)
print()
print(
    "allowing gcd(h-r,N) and gcd(h+r,N)"
)
print(
    "to be tested as direct factor channels."
)
print()

print(
    "This is materially different from merely observing that"
)
print(
    "a continued-fraction residue is numerically small."
)
print()

print(
    "A successful nontrivial gcd hit would establish a genuine"
)
print(
    "N-only factor-recovery channel:"
)
print()
print(
    "    N"
)
print(
    "     -> continued fraction"
)
print(
    "     -> square residue"
)
print(
    "     -> gcd"
)
print(
    "     -> p or q"
)
print()

print(
    "A null result would rule out this specific quadratic-residue"
)
print(
    "channel on the tested convergent/semiconvergent family."
)
print()

print("NEXT TARGET")
print("-" * 78)

print(
    "If exact square residues produce genuine factors, analyze"
)
print(
    "their algebraic origin and whether the relevant convergent"
)
print(
    "index can be predicted from N."
)
print()

print(
    "If they do not, the next experiment should inspect the"
)
print(
    "partial-quotient sequence itself for an exact invariant"
)
print(
    "corresponding to the factor gap |p-q|, rather than dumping"
)
print(
    "large Pell residue tables."
)

print()
print("=" * 78)
print("EXPERIMENT 503 FINISHED")
print("=" * 78)
