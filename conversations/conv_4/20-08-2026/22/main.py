#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 499 START")
print("=" * 78)
print("UPSTREAM SYMMETRIC-BRIDGE TARGET: H2(x), SHIFTED PRODUCTS, AND DISCRIMINANT")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S, Delta = sp.symbols("N S Delta")
x, a, b = sp.symbols("x a b")

z = sp.symbols("z")


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def symmetrize(expr):
    """
    Rewrite a symmetric p,q expression in N=pq and S=p+q.
    """
    expr = sp.expand(expr)

    sym_expr, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ValueError(
            f"Expression is not symmetric in p,q: {remainder}"
        )

    u1 = mapping[0][0]
    u2 = mapping[1][0]

    return sp.expand(
        sym_expr.subs(
            {
                u1: S,
                u2: N,
            }
        )
    )


# ============================================================================
# BASIC DEFINITIONS
# ============================================================================

def M(xv):
    return sp.expand(
        N + S*xv + xv**2
    )


def H2(xv):
    return sp.expand(
        6*N - S**2 + S*xv
    )


Delta_from_NS = sp.expand(
    S**2 - 4*N
)

print("[1] BASIC SYMMETRIC VARIABLES")
print("-" * 78)

print("  N     = p*q")
print("  S     = p+q")
print("  Delta = S^2 - 4*N")
print("  Delta = (p-q)^2")
print()

print(
    "  Verified Delta difference =",
    fact(
        Delta_from_NS
        - (
            (p + q)**2
            - 4*p*q
        )
    )
)

print()


# ============================================================================
# H2 REWRITE IN TERMS OF DELTA
# ============================================================================

print("[2] H2(x) IN DISCRIMINANT COORDINATES")
print("-" * 78)

h2 = H2(x)

h2_delta = sp.expand(
    2*N - Delta + S*x
).subs(
    Delta,
    Delta_from_NS
)

print(
    "  H2(x)          =",
    fact(h2)
)

print(
    "  2N-Delta+S*x   =",
    fact(h2_delta)
)

print(
    "  difference      =",
    fact(
        h2 - h2_delta
    )
)

print(
    "  PASS =",
    zero(
        h2 - h2_delta
    )
)

print()


# ============================================================================
# H2 IN TERMS OF SHIFTED PRODUCT M_x
# ============================================================================

print("[3] H2(x) <-> SHIFTED PRODUCT M_x")
print("-" * 78)

mx = M(x)

candidate_h2 = sp.expand(
    mx + N - Delta - x**2
).subs(
    Delta,
    Delta_from_NS
)

print(
    "  M_x =",
    fact(mx)
)

print(
    "  Candidate:"
)

print(
    "    M_x + N - Delta - x^2"
)

print(
    "  =", 
    fact(candidate_h2)
)

print(
    "  H2(x) - candidate =",
    fact(
        h2 - candidate_h2
    )
)

print(
    "  PASS =",
    zero(
        h2 - candidate_h2
    )
)

print()


# ============================================================================
# SHIFTED PRODUCT DIFFERENCE
# ============================================================================

print("[4] SHIFTED PRODUCT PAIR")
print("-" * 78)

Ma = M(a)
Mb = M(b)

print(
    "  M_a =", fact(Ma)
)

print(
    "  M_b =", fact(Mb)
)

print(
    "  M_a - M_b =",
    fact(
        Ma - Mb
    )
)

print(
    "  Expected   =",
    fact(
        S*(a-b) + a**2 - b**2
    )
)

print(
    "  PASS =",
    zero(
        Ma - Mb
        - (
            S*(a-b)
            + a**2
            - b**2
        )
    )
)

print()


# ============================================================================
# RECOVER S FROM TWO SHIFTED PRODUCTS
# ============================================================================

print("[5] S FROM TWO SHIFTED PRODUCTS")
print("-" * 78)

S_from_M = sp.cancel(
    (
        Ma
        - Mb
        - a**2
        + b**2
    )
    / (a-b)
)

print(
    "  recovered S =",
    fact(S_from_M)
)

print(
    "  difference   =",
    fact(
        S_from_M - S
    )
)

print(
    "  PASS =",
    zero(
        S_from_M - S
    )
)

print()


# ============================================================================
# DISCRIMINANT FROM SHIFTED PRODUCTS
# ============================================================================

print("[6] DISCRIMINANT FROM SHIFTED PRODUCTS")
print("-" * 78)

# M_a = N + aS + a^2
# M_b = N + bS + b^2
#
# First recover S, then Delta=S^2-4N.

Delta_from_M = sp.cancel(
    S_from_M**2 - 4*N
)

print(
    "  Delta(M_a,M_b) ="
)

print(
    "   ",
    fact(Delta_from_M)
)

print(
    "  Delta - Delta(M_a,M_b) =",
    fact(
        Delta_from_NS - Delta_from_M
    )
)

print(
    "  PASS =",
    zero(
        Delta_from_NS - Delta_from_M
    )
)

print()


# ============================================================================
# SPECIAL +1 / -1 SHIFTED PRODUCTS
# ============================================================================

print("[7] SPECIAL SHIFTED PRODUCTS M_+1 AND M_-1")
print("-" * 78)

Mp = sp.expand(
    M(1)
)

Mm = sp.expand(
    M(-1)
)

print(
    "  M_+1 =", fact(Mp)
)

print(
    "  M_-1 =", fact(Mm)
)

print(
    "  M_+1 + M_-1 =",
    fact(
        Mp + Mm
    )
)

print(
    "  M_+1 - M_-1 =",
    fact(
        Mp - Mm
    )
)

print(
    "  M_+1*M_-1 =",
    fact(
        Mp*Mm
    )
)

print()


# ============================================================================
# SPECIAL DISCRIMINANT IDENTITY
# ============================================================================

print("[8] PRODUCT IDENTITY FOR M_+1 M_-1")
print("-" * 78)

special_delta = sp.expand(
    (N - 1)**2
    - Mp*Mm
)

print(
    "  (N-1)^2 - M_+1*M_-1 =",
    fact(special_delta)
)

print(
    "  Delta =",
    fact(Delta_from_NS)
)

print(
    "  difference =",
    fact(
        special_delta
        - Delta_from_NS
    )
)

print(
    "  PASS =",
    zero(
        special_delta
        - Delta_from_NS
    )
)

print()


# ============================================================================
# H2 FROM M+ AND M-
# ============================================================================

print("[9] H2(x) FROM SHIFTED-PRODUCT DATA")
print("-" * 78)

# For x=1:
#
# H2(1)=M1+N-Delta-1
#
# For x=-1:
#
# H2(-1)=M_-1+N-Delta-1

H2_plus_from_M = sp.expand(
    Mp
    + N
    - Delta
    - 1
).subs(
    Delta,
    Delta_from_NS
)

H2_minus_from_M = sp.expand(
    Mm
    + N
    - Delta
    - 1
).subs(
    Delta,
    Delta_from_NS
)

print(
    "  H2(+1) reconstructed =",
    fact(H2_plus_from_M)
)

print(
    "  target H2(+1)       =",
    fact(H2(1))
)

print(
    "  PASS =",
    zero(
        H2_plus_from_M
        - H2(1)
    )
)

print()

print(
    "  H2(-1) reconstructed =",
    fact(H2_minus_from_M)
)

print(
    "  target H2(-1)       =",
    fact(H2(-1))
)

print(
    "  PASS =",
    zero(
        H2_minus_from_M
        - H2(-1)
    )
)

print()


# ============================================================================
# H2 PRODUCT AND DISCRIMINANT
# ============================================================================

print("[10] H2(+1)H2(-1) STRUCTURE")
print("-" * 78)

h2_product = sp.expand(
    H2(1)*H2(-1)
)

print(
    "  H2(+1)H2(-1) =",
    fact(h2_product)
)

print()

print(
    "  Alternative form:"
)

alternative_product = sp.expand(
    (
        6*N - S**2
    )**2
    - S**2
)

print(
    "  (6N-S^2)^2-S^2 =",
    fact(alternative_product)
)

print(
    "  Difference =",
    fact(
        h2_product
        - alternative_product
    )
)

print(
    "  PASS =",
    zero(
        h2_product
        - alternative_product
    )
)

print()


# ============================================================================
# EXPRESS H2 PRODUCT THROUGH DELTA
# ============================================================================

print("[11] H2 PRODUCT IN (N,S,DELTA)")
print("-" * 78)

h2_product_delta = sp.expand(
    (
        2*N - Delta
    )**2
    - S**2
)

print(
    "  Candidate = (2N-Delta)^2-S^2"
)

print(
    "  =", 
    fact(
        h2_product_delta.subs(
            Delta,
            Delta_from_NS
        )
    )
)

print(
    "  Difference =",
    fact(
        h2_product
        - h2_product_delta.subs(
            Delta,
            Delta_from_NS
        )
    )
)

print(
    "  PASS =",
    zero(
        h2_product
        - h2_product_delta.subs(
            Delta,
            Delta_from_NS
        )
    )
)

print()


# ============================================================================
# SEARCH FOR LOW-DEGREE EXPRESSIONS IN N AND SHIFTED PRODUCTS
# ============================================================================

print("[12] LOW-DEGREE SHIFTED-PRODUCT BRIDGE SEARCH")
print("-" * 78)

print(
    "Search whether H2(+1) or H2(-1) can be expressed"
)

print(
    "using N, M_+1, M_-1 with small polynomial degree."
)

basis_vars = [
    1,
    N,
    Mp,
    Mm,
    N**2,
    N*Mp,
    N*Mm,
    Mp**2,
    Mp*Mm,
    Mm**2,
]

targets = {
    "H2(+1)": sp.expand(H2(1)),
    "H2(-1)": sp.expand(H2(-1)),
}

for target_name, target in targets.items():

    coeffs = sp.symbols(
        "c0:" + str(len(basis_vars))
    )

    candidate = sum(
        coeff*term
        for coeff, term
        in zip(
            coeffs,
            basis_vars
        )
    )

    poly = sp.Poly(
        sp.expand(
            candidate - target
        ),
        N,
        S
    )

    equations = [
        coeff
        for _, coeff
        in poly.terms()
    ]

    solution_set = sp.linsolve(
        equations,
        coeffs
    )

    hits = []

    for solution in solution_set:

        if any(
            value != 0
            for value in solution
        ):
            hits.append(solution)

    print(
        f"  {target_name}: "
        f"hits = {len(hits)}"
    )

    for sol in hits[:3]:

        print(
            "    ",
            sol
        )

print()


# ============================================================================
# SEARCH FOR H2 AS FUNCTION OF M_x AND DELTA
# ============================================================================

print("[13] MINIMAL H2 REPRESENTATION")
print("-" * 78)

candidate_minimal = sp.expand(
    M(x)
    + N
    - Delta
    - x**2
).subs(
    Delta,
    Delta_from_NS
)

print(
    "  Candidate:"
)

print(
    "    H2(x) = M_x + N - Delta - x^2"
)

print(
    "  Candidate expanded =",
    fact(candidate_minimal)
)

print(
    "  Target              =",
    fact(H2(x))
)

print(
    "  Difference           =",
    fact(
        candidate_minimal
        - H2(x)
    )
)

print(
    "  PASS =",
    zero(
        candidate_minimal
        - H2(x)
    )
)

print()


# ============================================================================
# NUMERICAL AUDIT
# ============================================================================

print("[14] NUMERICAL AUDIT")
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
    dd = (pp - qq)**2

    for xx in [1, -1, 2, -2]:

        true_m = (
            pp + xx
        ) * (
            qq + xx
        )

        true_h = (
            6*nn
            - ss**2
            + ss*xx
        )

        reconstructed_h = (
            true_m
            + nn
            - dd
            - xx**2
        )

        passed = (
            reconstructed_h
            == true_h
        )

        if not passed:
            numeric_failures += 1

        print(
            f"  ({pp},{qq}), x={xx}: "
            f"H2={true_h}, "
            f"Mx+N-Delta-x^2={reconstructed_h}, "
            f"PASS={passed}"
        )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# INFORMATION-TARGET ANALYSIS
# ============================================================================

print("[15] INFORMATION-TARGET ANALYSIS")
print("-" * 78)

print(
    "The identity"
)

print(
    "  H2(x) = M_x + N - Delta - x^2"
)

print(
    "shows that the missing scalar H2(x)"
)

print(
    "is equivalent to the pair"
)

print(
    "  (M_x, Delta)"
)

print(
    "once N and x are known."
)

print()

print(
    "Conversely:"
)

print(
    "  M_x = H2(x) - N + Delta + x^2"
)

print(
    "so H2(x) and M_x carry the same"
)

print(
    "symmetric information once Delta is available."
)

print()

print(
    "For x=1 and x=-1:"
)

print(
    "  Delta = (N-1)^2 - M_1*M_-1"
)

print(
    "and therefore the old shifted-product route"
)

print(
    "can be expressed entirely through"
)

print(
    "the pair (M_1,M_-1) together with N."
)

print()


# ============================================================================
# CENTRAL UPSTREAM QUESTION
# ============================================================================

print("[16] CENTRAL UPSTREAM QUESTION")
print("-" * 78)

print(
    "The downstream algebra is now completely characterized:"
)

print(
    "  M_x + N - Delta - x^2"
)

print(
    "       -> H2(x)"
)

print(
    "       -> S"
)

print(
    "       -> p,q."
)

print()

print(
    "Therefore the upstream research should search for"
)

print(
    "either:"
)

print(
    "  (A) M_x and Delta,"
)

print(
    "or:"
)

print(
    "  (B) one scalar equal to H2(x),"
)

print(
    "or:"
)

print(
    "  (C) M_1*M_-1 and another independent"
)

print(
    "      shifted-product quantity."
)

print()

print(
    "These are now equivalent bridge targets."
)

print()


# ============================================================================
# STATUS
# ============================================================================

two_shift_pass = zero(
    S_from_M - S
)

delta_shift_pass = zero(
    special_delta
    - Delta_from_NS
)

h2_mx_pass = zero(
    candidate_minimal
    - H2(x)
)

h2_plus_pass = zero(
    H2_plus_from_M
    - H2(1)
)

h2_minus_pass = zero(
    H2_minus_from_M
    - H2(-1)
)

product_pass = zero(
    h2_product
    - alternative_product
)

overall = (
    two_shift_pass
    and delta_shift_pass
    and h2_mx_pass
    and h2_plus_pass
    and h2_minus_pass
    and product_pass
    and numeric_failures == 0
)

print("[17] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  two-shift S recovery =",
    two_shift_pass
)

print(
    "  Delta from shifted products =",
    delta_shift_pass
)

print(
    "  H2 = Mx+N-Delta-x^2 =",
    h2_mx_pass
)

print(
    "  H2(+1) reconstruction =",
    h2_plus_pass
)

print(
    "  H2(-1) reconstruction =",
    h2_minus_pass
)

print(
    "  H2(+1)H2(-1) identity =",
    product_pass
)

print(
    "  numerical audit =",
    numeric_failures == 0
)

print(
    "  OVERALL EXACT AUDIT =",
    overall
)

print()

print("MAIN RESULT")
print("-" * 78)

print(
    "The experiment establishes the exact identity"
)

print(
    "  H2(x) = M_x + N - (p-q)^2 - x^2."
)

print(
    "Since"
)

print(
    "  M_x = (p+x)(q+x),"
)

print(
    "and"
)

print(
    "  Delta = (p-q)^2,"
)

print(
    "the shifted-kernel target is algebraically equivalent"
)

print(
    "to recovering the shifted symmetric product together"
)

print(
    "with the factor gap-square."
)

print()

print(
    "For the special shifts ±1:"
)

print(
    "  Delta = (N-1)^2 - M_1 M_-1."
)

print()

print(
    "Thus the old 2023 shifted-product idea and the modern"
)

print(
    "KAPPA kernel are connected through the discriminant"
)

print(
    "Delta = (p-q)^2."
)

print()

print(
    "NEXT TARGET"
)

print(
    "-" * 78
)

print(
    "Search the existing N-only homogeneous-layer construction"
)

print(
    "for an exact object corresponding to either:"
)

print(
    "  Delta = (p-q)^2,"
)

print(
    "  M_1*M_-1,"
)

print(
    "  M_x,"
)

print(
    "or directly H2(x)."
)

print()

print(
    "The most promising upstream candidate is now the"
)

print(
    "factor-gap/discriminant channel:"
)

print(
    "  N-only construction"
)

print(
    "       -> (p-q)^2"
)

print(
    "       -> M_x or H2(x)"
)

print(
    "       -> S"
)

print(
    "       -> p,q."
)

print()

print("=" * 78)
print("EXPERIMENT 499 FINISHED")
print("=" * 78)
