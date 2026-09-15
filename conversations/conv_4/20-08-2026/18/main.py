#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 495 START")
print("=" * 78)
print("SHIFTED-KERNEL GENERATING FUNCTION AND S-ELIMINATION SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q, x, z, t = sp.symbols("p q x z t")
N, S = sp.symbols("N S")


def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def cancel(expr):
    return sp.cancel(sp.together(expr))


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(expr):
    """
    Reduce a symmetric polynomial in p,q to N=pq, S=p+q.
    """
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ValueError(
            f"Expression was not symmetric. Remainder={remainder}"
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    reduced = reduced.subs({
        s1: S,
        s2: N,
    })

    return sp.expand(reduced)


# ============================================================================
# SHIFTED KERNEL
# ============================================================================

def F_shifted(k, ell):
    return sp.expand(
        p**k * (q + x)**ell
        + q**k * (p + x)**ell
        - p**ell * (q + x)**k
        - q**ell * (p + x)**k
    )


def H_shifted(k, ell):
    F = F_shifted(k, ell)
    return sp.expand(sp.cancel(F / x))


# ============================================================================
# [1] BUILD H_(1,ell)
# ============================================================================

print("[1] BUILD k=1 SHIFTED KERNEL")
print("-" * 78)

H = {}

for ell in range(2, 11):
    raw = H_shifted(1, ell)
    H[ell] = symmetric_reduce(raw)

    print(
        f"  ell={ell}: degree={sp.Poly(H[ell], x).degree()}"
    )

print()


# ============================================================================
# [2] FOUR-EXPONENTIAL DECOMPOSITION IN ell
# ============================================================================

print("[2] FOUR-EXPONENTIAL STRUCTURE")
print("-" * 78)

print(
    "For fixed x, the shifted kernel has bases:"
)

print(
    "  p+x, q+x, p, q"
)

print(
    "and therefore admits a quartic characteristic polynomial."
)

print()

lam = [
    p + x,
    q + x,
    p,
    q,
]

weights = [
    p / x,
    q / x,
    -(q + x) / x,
    -(p + x) / x,
]

# Check the decomposition numerically/symbolically.
decomp_failures = 0

for ell in range(2, 11):

    reconstructed = sum(
        weights[i] * lam[i]**ell
        for i in range(4)
    )

    difference = sp.expand(
        reconstructed
        - F_shifted(1, ell) / x
    )

    if not zero(difference):
        decomp_failures += 1

    print(
        f"  ell={ell}: PASS={zero(difference)}"
    )

print(
    "  DECOMPOSITION FAILURES =",
    decomp_failures
)

print()


# ============================================================================
# [3] CHARACTERISTIC POLYNOMIAL
# ============================================================================

print("[3] SHIFTED CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi_x = sp.expand(
    (z - p - x)
    * (z - q - x)
    * (z - p)
    * (z - q)
)

chi_ns = symmetric_reduce(
    chi_x
)

print(
    "  chi_x(z) =",
    chi_ns
)

poly_chi = sp.Poly(
    chi_ns,
    z
)

coeffs = poly_chi.all_coeffs()

print()
print("  coefficients:")

for i, coeff in enumerate(coeffs):
    print(
        f"    z^{4-i}: {sp.factor(coeff)}"
    )

print()


# ============================================================================
# [4] EXPLICIT SYMMETRIC CHARACTERISTIC COEFFICIENTS
# ============================================================================

print("[4] CHARACTERISTIC COEFFICIENT STRUCTURE")
print("-" * 78)

a3 = sp.expand(
    2 * S + 2 * x
)

a2 = sp.expand(
    2 * N
    + S**2
    + 3 * S * x
    + x**2
)

a1 = sp.expand(
    (2 * N + S**2) * x
    + 2 * N * S
)

a0 = sp.expand(
    N * (N + S * x + x**2)
)

expected_chi = sp.expand(
    z**4
    - a3 * z**3
    + a2 * z**2
    - a1 * z
    + a0
)

chi_difference = sp.expand(
    chi_ns - expected_chi
)

print(
    "  expected form =",
    expected_chi
)

print(
    "  difference =",
    chi_difference
)

print(
    "  PASS =",
    zero(chi_difference)
)

print()


# ============================================================================
# [5] SPECIAL VALUES x = 0,1,-1
# ============================================================================

print("[5] SPECIAL SHIFT CHARACTERISTICS")
print("-" * 78)

for xv in [0, 1, -1, 2, -2]:

    specialized = sp.expand(
        chi_ns.subs(x, xv)
    )

    print(
        f"  x={xv}:"
    )

    print(
        f"    chi(z) = {fact(specialized)}"
    )

print()


# ============================================================================
# [6] GENERATING FUNCTION
# ============================================================================

print("[6] ORDINARY GENERATING FUNCTION IN ell")
print("-" * 78)

print(
    "Define"
)

print(
    "  G_x(t) = sum_{ell>=1} F_ell(x) t^ell."
)

print(
    "The four-exponential decomposition gives"
)

print(
    "  G_x(t) = sum_i c_i * lambda_i*t/(1-lambda_i*t)."
)

print()

G = 0

for i in range(4):

    G += (
        weights[i]
        * lam[i]
        * t
        / (1 - lam[i] * t)
    )

G = sp.factor(
    sp.together(G)
)

print(
    "  G_x(t) =",
    G
)

print()


# ============================================================================
# [7] GENERATING FUNCTION NUMERATOR / DENOMINATOR
# ============================================================================

print("[7] GENERATING FUNCTION DENOMINATOR")
print("-" * 78)

G_together = sp.together(G)

G_num, G_den = sp.fraction(
    G_together
)

G_num = sp.expand(G_num)
G_den = sp.expand(G_den)

G_num_ns = symmetric_reduce(
    G_num
)

G_den_ns = symmetric_reduce(
    G_den
)

print(
    "  numerator ="
)

print(
    f"    {sp.factor(G_num_ns)}"
)

print(
    "  denominator ="
)

print(
    f"    {sp.factor(G_den_ns)}"
)

print()


# ============================================================================
# [8] DENOMINATOR / CHARACTERISTIC CONNECTION
# ============================================================================

print("[8] DENOMINATOR / CHARACTERISTIC CONNECTION")
print("-" * 78)

expected_den = sp.expand(
    1
    - a3 * t
    + a2 * t**2
    - a1 * t**3
    + a0 * t**4
)

den_difference = sp.expand(
    G_den_ns - expected_den
)

print(
    "  denominator difference =",
    den_difference
)

print(
    "  PASS =",
    zero(den_difference)
)

print()


# ============================================================================
# [9] GENERATING-FUNCTION NUMERATOR AUDIT
# ============================================================================

print("[9] NUMERATOR AUDIT")
print("-" * 78)

print(
    "The numerator is inspected for simple factors"
)

num_factor = fact(G_num_ns)

print(
    "  numerator factorization ="
)

print(
    f"    {num_factor}"
)

print()


# ============================================================================
# [10] x -> -x COMPARISON
# ============================================================================

print("[10] x -> -x GENERATING-FUNCTION COMPARISON")
print("-" * 78)

G_plus_num, G_plus_den = sp.fraction(
    sp.together(
        G.subs(x, 1)
    )
)

G_minus_num, G_minus_den = sp.fraction(
    sp.together(
        G.subs(x, -1)
    )
)

G_plus_num_ns = symmetric_reduce(
    sp.expand(G_plus_num)
)

G_minus_num_ns = symmetric_reduce(
    sp.expand(G_minus_num)
)

G_plus_den_ns = symmetric_reduce(
    sp.expand(G_plus_den)
)

G_minus_den_ns = symmetric_reduce(
    sp.expand(G_minus_den)
)

print(
    "  x=+1 numerator =",
    fact(G_plus_num_ns)
)

print(
    "  x=-1 numerator =",
    fact(G_minus_num_ns)
)

print(
    "  x=+1 denominator =",
    fact(G_plus_den_ns)
)

print(
    "  x=-1 denominator =",
    fact(G_minus_den_ns)
)

print()


# ============================================================================
# [11] DENOMINATOR SUM / DIFFERENCE
# ============================================================================

print("[11] x=+1 / x=-1 DENOMINATOR COMBINATIONS")
print("-" * 78)

den_sum = sp.expand(
    G_plus_den_ns + G_minus_den_ns
)

den_diff = sp.expand(
    G_plus_den_ns - G_minus_den_ns
)

print(
    "  D(+1)+D(-1) =",
    fact(den_sum)
)

print(
    "  D(+1)-D(-1) =",
    fact(den_diff)
)

print()


# ============================================================================
# [12] NUMERATOR SUM / DIFFERENCE
# ============================================================================

print("[12] x=+1 / x=-1 NUMERATOR COMBINATIONS")
print("-" * 78)

num_sum = sp.expand(
    G_plus_num_ns + G_minus_num_ns
)

num_diff = sp.expand(
    G_plus_num_ns - G_minus_num_ns
)

print(
    "  Num(+1)+Num(-1) =",
    fact(num_sum)
)

print(
    "  Num(+1)-Num(-1) =",
    fact(num_diff)
)

print()


# ============================================================================
# [13] S-ELIMINATION FROM CHARACTERISTIC COEFFICIENTS
# ============================================================================

print("[13] S-ELIMINATION SEARCH")
print("-" * 78)

print(
    "Treat the characteristic coefficients"
)

print(
    "  a3, a2, a1, a0"
)

print(
    "as observable symbolic quantities."
)

print(
    "Test simple combinations for exact cancellation of S."
)

print()

targets = {
    "a3": a3,
    "a2": a2,
    "a1": a1,
    "a0": a0,
}

for name, expr in targets.items():

    print(
        f"  {name} = {sp.factor(expr)}"
    )

print()


# Simple discriminant-like quantities.

disc01 = sp.expand(
    a3**2 - 4*a2
)

disc02 = sp.expand(
    a3 * a1 - 4 * a0
)

mixed1 = sp.expand(
    a2 - sp.Rational(3, 4) * a3**2
)

print(
    "  a3^2 - 4*a2 =",
    fact(disc01)
)

print(
    "  a3*a1 - 4*a0 =",
    fact(disc02)
)

print(
    "  a2 - 3/4*a3^2 =",
    fact(mixed1)
)

print()


# ============================================================================
# [14] NORMALIZE BY x
# ============================================================================

print("[14] x-NORMALIZED COEFFICIENTS")
print("-" * 78)

for name, expr in targets.items():

    # Even/odd decomposition with respect to x.
    plus = sp.expand(
        expr.subs(x, x) + expr.subs(x, -x)
    )

    minus = sp.expand(
        expr.subs(x, x) - expr.subs(x, -x)
    )

    print(
        f"  {name}:"
    )

    print(
        f"    even = {fact(plus / 2)}"
    )

    print(
        f"    odd/x = {fact(minus / (2*x))}"
    )

print()


# ============================================================================
# [15] SPECTRAL NULLSPACE RECHECK
# ============================================================================

print("[15] INTERPRETATION OF 494 SPECTRAL NULLSPACE")
print("-" * 78)

print(
    "Experiment 494 found a 2-dimensional spectral nullspace."
)

print(
    "We now test elementary operator identities that can explain"
)

print(
    "such a nullspace without introducing a new invariant."
)

print()

# Euler operator E = x d/dx.
Euler = lambda expr: sp.expand(
    x * sp.diff(expr, x)
)

# Identity:
# x H' - E(H) = 0
# This is deliberately tautological and serves as a control.

tautological_control = sp.expand(
    Euler(H[4]) - x * sp.diff(H[4], x)
)

print(
    "  Euler identity control =",
    tautological_control
)

print(
    "  PASS =",
    zero(tautological_control)
)

# Test whether H has a universal homogeneity relation in x,N,S.
homogeneity_candidates = [
    x * sp.diff(H[4], x),
    N * sp.diff(H[4], N),
    S * sp.diff(H[4], S),
]

for idx, expr in enumerate(
    homogeneity_candidates,
    start=1
):

    print(
        f"  homogeneous candidate {idx} degree =",
        sp.Poly(
            expr,
            x,
            N,
            S
        ).total_degree()
    )

print()


# ============================================================================
# [16] DIRECT N-ONLY CANDIDATE SEARCH IN GENERATING-FUNCTION COEFFICIENTS
# ============================================================================

print("[16] N-ONLY GENERATING-FUNCTION COMBINATION SEARCH")
print("-" * 78)

# Search a small fixed family of combinations of characteristic coefficients:
#
# a2
# a3^2
# a3*x
# a1/x
# a0/x^2
#
# and their simple integer linear combinations.

feature_names = [
    "a2",
    "a3^2",
    "a3*x",
    "a1/x",
    "a0/x^2",
]

features = [
    a2,
    sp.expand(a3**2),
    sp.expand(a3*x),
    sp.expand(a1 / x),
    sp.expand(a0 / x**2),
]

# Remove denominators safely.
cleared_features = []

for f in features:

    numerator, denominator = sp.fraction(
        sp.cancel(f)
    )

    cleared_features.append(
        (
            sp.expand(numerator),
            sp.expand(denominator)
        )
    )

print(
    "  fixed feature family size =",
    len(features)
)

print()


# Test pairwise differences/sums for S independence.
s_independent_hits = []

for i in range(len(features)):

    for j in range(
        i + 1,
        len(features)
    ):

        expr = sp.cancel(
            features[i] - features[j]
        )

        num, den = sp.fraction(
            sp.together(expr)
        )

        num = sp.expand(num)

        dS = sp.diff(
            num,
            S
        )

        # A numerator independent of S is a useful candidate.
        if zero(dS):

            s_independent_hits.append(
                (
                    feature_names[i],
                    feature_names[j],
                    fact(num),
                    fact(den)
                )
            )

print(
    "  pairwise S-independent hits =",
    len(s_independent_hits)
)

for hit in s_independent_hits:

    left, right, numerator, denominator = hit

    print(
        f"  {left} - {right}:"
    )

    print(
        f"    numerator   = {numerator}"
    )

    print(
        f"    denominator = {denominator}"
    )

print()


# ============================================================================
# [17] NUMERICAL SANITY
# ============================================================================

print("[17] NUMERICAL SANITY")
print("-" * 78)

instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]

numeric_failures = 0

for pp, qq in instances:

    nn = pp * qq
    ss = pp + qq

    local_pass = True

    for ell in range(2, 8):

        raw = F_shifted(
            1,
            ell
        )

        h_numeric = sp.expand(
            sp.cancel(
                raw.subs(
                    {
                        p: pp,
                        q: qq
                    }
                ) / x
            )
        )

        top = sp.Poly(
            h_numeric,
            x
        ).coeff_monomial(
            x**(ell - 1)
        )

        if top != ss:
            local_pass = False

    if not local_pass:
        numeric_failures += 1

    print(
        f"  ({pp},{qq}) PASS={local_pass}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [18] FINAL STATUS
# ============================================================================

print("[18] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  shifted kernel decomposition =",
    decomp_failures == 0
)

print(
    "  characteristic polynomial =",
    zero(chi_difference)
)

print(
    "  generating denominator =",
    zero(den_difference)
)

print(
    "  spectral nullspace from 494 is explained only by explicit controls"
)

print(
    "  N-only coefficient-combination hits =",
    len(s_independent_hits)
)

print(
    "  numerical sanity =",
    numeric_failures == 0
)

print()

print(
    "INTERPRETATION"
)

print(
    "----------------------------------------------------------------------------"
)

print(
    "The shifted kernel is a four-exponential sequence in ell with bases"
)

print(
    "  p, q, p+x, q+x."
)

print()

print(
    "Its generating function therefore has an exact quartic denominator."
)

print()

print(
    "The important question is whether combinations of the generating"
)

print(
    "function coefficients remove S while retaining nontrivial information"
)

print(
    "about N or the shifted product."
)

print()

print(
    "A positive result would give a new transfer mechanism:"
)

print(
    "  kernel generating data -> N-only invariant -> M_x -> S -> p,q"
)

print()

print(
    "A negative result means the search should move upstream toward"
)

print(
    "the original N-only homogeneous-layer construction."
)

print()


print("=" * 78)
print("EXPERIMENT 495 FINISHED")
print("=" * 78)
