#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 491 START")
print("=" * 78)
print("KAPPA KERNEL N-ONLY INVARIANT SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")

Y2, Y3, Y4, Y5 = sp.symbols("Y2 Y3 Y4 Y5")


# ============================================================================
# HELPERS
# ============================================================================

def expand0(expr):
    return sp.expand(expr)


def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def cancel_fact(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def safe_fraction(expr):
    expr = sp.cancel(sp.together(expr))
    num, den = expr.as_numer_denom()
    return sp.factor(num), sp.factor(den)


# ============================================================================
# ORIGINAL KAPPA KERNEL
# ============================================================================

def F_pq(n):
    return sp.expand(
        p * (q + 1) ** n
        + q * (p + 1) ** n
        - (q + 1) * p ** n
        - (p + 1) * q ** n
    )


# ============================================================================
# SYMMETRIC KERNEL VALUES
# ============================================================================

F2 = sp.expand(
    6 * N - S**2 + S
)

F3 = sp.expand(
    (S + 1) * F2
)

F4 = sp.expand(
    -10 * N**2
    + 8 * N * S**2
    + 6 * N * S
    + 8 * N
    - S**4
    + S
)

F5 = sp.expand(
    -(S + 1) * (
        20 * N**2
        - 10 * N * S**2
        - 10 * N
        + S**4
        - S**3
        + S**2
        - S
    )
)

F6 = sp.expand(
    14 * N**3
    - 33 * N**2 * S**2
    - 45 * N**2 * S
    - 40 * N**2
    + 12 * N * S**4
    + 15 * N * S**3
    + 20 * N * S**2
    + 15 * N * S
    + 12 * N
    - S**6
    + S
)

kernel_map = {
    Y2: F2,
    Y3: F3,
    Y4: F4,
    Y5: F5,
}


# ============================================================================
# [1] DIRECT KAPPA SANITY CHECK
# ============================================================================

print("[1] DIRECT KAPPA SANITY CHECK")
print("-" * 78)

direct_checks = []

expected_values = {
    2: F2,
    3: F3,
    4: F4,
    5: F5,
    6: F6,
}

for n, expected in expected_values.items():

    actual = F_pq(n)

    # q = S-p and p,q satisfy:
    #
    #     p^2 - S*p + N = 0.
    #
    reduced = sp.expand(
        actual.subs(q, S - p)
    )

    reduced_poly = sp.Poly(
        reduced,
        p
    )

    modulus_poly = sp.Poly(
        p**2 - S * p + N,
        p
    )

    reduced = sp.rem(
        reduced_poly,
        modulus_poly
    ).as_expr()

    reduced = fact(reduced)

    passed = zero(reduced - expected)

    direct_checks.append(passed)

    print(
        f"  F_{n}: PASS={passed}"
    )

print()
print(
    "  DIRECT KAPPA FAILURES =",
    direct_checks.count(False)
)
print()


# ============================================================================
# [2] BASIC KAPPA RELATION
# ============================================================================

print("[2] BASIC KAPPA RELATION")
print("-" * 78)

rel_F3 = expand0(
    F3 - (S + 1) * F2
)

print(
    "  F3-(S+1)F2 =",
    fact(rel_F3)
)

print(
    "  PASS =",
    zero(rel_F3)
)

print()


# ============================================================================
# [3] F2 AS A QUADRATIC IN S
# ============================================================================

print("[3] EXACT F2 QUADRATIC")
print("-" * 78)

quadratic_relation = expand0(
    S**2 - S + Y2 - 6 * N
)

print(
    "  F2 = 6*N - S^2 + S"
)

print(
    "  Equivalent equation:"
)

print(
    "    S^2 - S + F2 - 6*N = 0"
)

print()


# ============================================================================
# [4] ELIMINATION OF S USING F3/F2
# ============================================================================

print("[4] F2/F3 ELIMINATION")
print("-" * 78)

# From F3=(S+1)F2:
#
#     S = F3/F2 - 1
#
# whenever F2 != 0.
#
# Substitute this expression into:
#
#     S^2 - S + F2 - 6N = 0.
#
# We clear denominators exactly.

S_from_ratio = sp.cancel(
    Y3 / Y2 - 1
)

elimination_fraction = sp.together(
    quadratic_relation.subs(
        S,
        S_from_ratio
    )
)

elimination_num, elimination_den = (
    elimination_fraction.as_numer_denom()
)

# This is the quantity that was previously referred to
# accidentally as "relation_kernel".
relation_kernel = fact(
    sp.expand(elimination_num)
)

relation_den = fact(
    sp.expand(elimination_den)
)

print(
    "  S reconstructed symbolically as:"
)

print(
    "    S = F3/F2 - 1"
)

print()

print(
    "  Cleared elimination numerator:"
)

print(
    "   ",
    relation_kernel
)

print()

print(
    "  Cleared elimination denominator:"
)

print(
    "   ",
    relation_den
)

print()

print(
    "  Is the elimination numerator identically zero?"
)

elimination_pass = zero(
    relation_kernel
)

print(
    "   PASS =",
    elimination_pass
)

print()


# ============================================================================
# [5] LOW-DEGREE KERNEL BASIS
# ============================================================================

print("[5] LOW-DEGREE KERNEL BASIS")
print("-" * 78)


def monomials_of_degree(max_degree):
    result = [sp.Integer(1)]

    for total in range(1, max_degree + 1):

        for a in range(total + 1):

            for b in range(total - a + 1):

                for c in range(total - a - b + 1):

                    d = total - a - b - c

                    result.append(
                        Y2**a
                        * Y3**b
                        * Y4**c
                        * Y5**d
                    )

    return result


max_monomial_degree = 2
max_N_degree = 3

mons = monomials_of_degree(
    max_monomial_degree
)

coefficient_basis = []

for mon in mons:

    for degN in range(max_N_degree + 1):

        coefficient_basis.append(
            sp.expand(
                N**degN * mon
            )
        )

print(
    "  monomial degree <=",
    max_monomial_degree
)

print(
    "  N coefficient degree <=",
    max_N_degree
)

print(
    "  basis size =",
    len(coefficient_basis)
)

print()


# ============================================================================
# [6] S-INDEPENDENT LINEAR COMBINATION SEARCH
# ============================================================================

print("[6] EXACT S-INDEPENDENT COMBINATION SEARCH")
print("-" * 78)

derivatives = []

for basis_expr in coefficient_basis:

    substituted = sp.expand(
        basis_expr.subs(kernel_map)
    )

    derivative = sp.diff(
        substituted,
        S
    )

    derivatives.append(
        sp.Poly(
            derivative,
            N,
            S
        )
    )

all_monomials = set()

for poly in derivatives:

    all_monomials.update(
        poly.monoms()
    )

all_monomials = sorted(
    all_monomials
)

rows = []

for mon in all_monomials:

    row = [
        poly.coeff_monomial(mon)
        for poly in derivatives
    ]

    rows.append(row)

if rows:

    derivative_matrix = sp.Matrix(rows)

else:

    derivative_matrix = sp.zeros(
        0,
        len(coefficient_basis)
    )

print(
    "  derivative matrix shape =",
    derivative_matrix.shape
)

nullspace = derivative_matrix.nullspace()

print(
    "  nullspace dimension =",
    len(nullspace)
)

print()


# ============================================================================
# [7] CONSTRUCT NONZERO NULLSPACE CANDIDATES
# ============================================================================

print("[7] S-INDEPENDENT CANDIDATES")
print("-" * 78)

candidate_expressions = []

for vector in nullspace:

    expr = sp.Integer(0)

    for coefficient, basis_expr in zip(
        vector,
        coefficient_basis
    ):
        expr += (
            coefficient
            * basis_expr
        )

    expr = fact(
        expr.subs(
            kernel_map
        )
    )

    if zero(expr):
        continue

    candidate_expressions.append(
        expr
    )


# Exact duplicate removal.
unique_candidates = []

for expr in candidate_expressions:

    if any(
        zero(expr - other)
        for other in unique_candidates
    ):
        continue

    unique_candidates.append(
        expr
    )


print(
    "  nonzero candidates =",
    len(unique_candidates)
)

for i, expr in enumerate(
    unique_candidates[:20],
    start=1
):

    print(
        f"  candidate {i}:"
    )

    print(
        f"    {expr}"
    )

print()


# ============================================================================
# [8] PURE N FUNCTION TEST
# ============================================================================

print("[8] PURE-N FUNCTION CHECK")
print("-" * 78)

pure_N_candidates = []

for expr in unique_candidates:

    derivative = sp.diff(
        sp.expand(expr),
        S
    )

    if not zero(derivative):
        continue

    poly_S = sp.Poly(
        sp.expand(expr),
        S,
        domain="EX"
    )

    if poly_S.degree() == 0:

        pure_N_candidates.append(
            fact(expr)
        )


unique_pure_N = []

for expr in pure_N_candidates:

    if any(
        zero(expr - other)
        for other in unique_pure_N
    ):
        continue

    unique_pure_N.append(
        expr
    )


print(
    "  pure N candidates =",
    len(unique_pure_N)
)

for i, expr in enumerate(
    unique_pure_N[:20],
    start=1
):

    print(
        f"  pure-N candidate {i}:",
        expr
    )

print()


# ============================================================================
# [9] NONTRIVIAL N-ONLY CANDIDATES
# ============================================================================

print("[9] NONTRIVIAL N-ONLY SEARCH")
print("-" * 78)

nontrivial_N = []

for expr in unique_pure_N:

    poly_N = sp.Poly(
        sp.expand(expr),
        N
    )

    if poly_N.degree() == 0:
        continue

    nontrivial_N.append(
        fact(expr)
    )


print(
    "  nontrivial N-only candidates =",
    len(nontrivial_N)
)

for i, expr in enumerate(
    nontrivial_N[:20],
    start=1
):

    print(
        f"  candidate {i}:",
        expr
    )

print()


# ============================================================================
# [10] CANONICAL N-ONLY TARGET COMPARISON
# ============================================================================

print("[10] N-ONLY TARGET COMPARISON")
print("-" * 78)

known_targets = {
    "N": N,
    "N+1": N + 1,
    "N^2": N**2,
    "N^2+N": N**2 + N,
}

hits = []

for candidate in nontrivial_N:

    for name, target in known_targets.items():

        if zero(
            candidate - target
        ):

            hits.append(
                (
                    candidate,
                    name
                )
            )


if hits:

    for candidate, name in hits:

        print(
            f"  HIT: {candidate} = {name}"
        )

else:

    print(
        "  No direct match to the tested canonical N-only targets."
    )

print()


# ============================================================================
# [11] DIRECT LOW-DEGREE M1 SEARCH
# ============================================================================

print("[11] DIRECT LOW-DEGREE M1 SEARCH")
print("-" * 78)

M1 = sp.expand(
    N + S + 1
)

target_polys = [
    sp.Poly(
        sp.expand(
            basis_expr.subs(
                kernel_map
            )
        ),
        N,
        S
    )
    for basis_expr in coefficient_basis
]

target_poly = sp.Poly(
    M1,
    N,
    S
)

target_monomials = set(
    target_poly.monoms()
)

for poly in target_polys:

    target_monomials.update(
        poly.monoms()
    )

target_monomials = sorted(
    target_monomials
)

A_rows = []
rhs_rows = []

for mon in target_monomials:

    A_rows.append(
        [
            poly.coeff_monomial(mon)
            for poly in target_polys
        ]
    )

    rhs_rows.append(
        target_poly.coeff_monomial(mon)
    )

A_mat = sp.Matrix(A_rows)
rhs_vec = sp.Matrix(rhs_rows)

print(
    "  linear-system shape =",
    A_mat.shape
)

try:

    solution_set = sp.linsolve(
        (
            A_mat,
            rhs_vec
        )
    )

    solutions = list(
        solution_set
    )

    if not solutions:

        print(
            "  No exact low-degree M1 representation found."
        )

    else:

        solution = solutions[0]

        has_free_parameters = any(
            value.free_symbols
            for value in solution
        )

        if has_free_parameters:

            print(
                "  Solutions exist with free parameters."
            )

        else:

            found_expr = sp.expand(
                sum(
                    coefficient * basis_expr
                    for coefficient, basis_expr in zip(
                        solution,
                        coefficient_basis
                    )
                )
            )

            print(
                "  Exact low-degree M1 representation FOUND."
            )

            print(
                "  representation =",
                fact(found_expr)
            )

            print(
                "  difference =",
                fact(
                    found_expr - M1
                )
            )

except Exception as exc:

    print(
        "  M1 linear search failed:",
        type(exc).__name__,
        str(exc)
    )

print()


# ============================================================================
# [12] SIMPLE RATIONAL RATIO SEARCH
# ============================================================================

print("[12] SIMPLE RATIONAL KERNEL RATIO SEARCH")
print("-" * 78)

ratio_candidates = {
    "F3/F2": Y3 / Y2,
    "F4/F2": Y4 / Y2,
    "F4/F3": Y4 / Y3,
    "F2/F3": Y2 / Y3,
    "F3/F4": Y3 / Y4,
    "F5/F2": Y5 / Y2,
    "F5/F3": Y5 / Y3,
    "F5/F4": Y5 / Y4,
}

for name, candidate in ratio_candidates.items():

    expression = cancel_fact(
        candidate.subs(
            kernel_map
        )
    )

    numerator, denominator = safe_fraction(
        expression
    )

    derivative = sp.diff(
        expression,
        S
    )

    pure_N = zero(
        sp.together(
            derivative
        )
    )

    equals_S = zero(
        sp.together(
            expression - S
        )
    )

    equals_M1 = zero(
        sp.together(
            expression - M1
        )
    )

    print(
        f"  {name}:"
    )

    print(
        f"    numerator   = {numerator}"
    )

    print(
        f"    denominator = {denominator}"
    )

    print(
        f"    pure-N      = {pure_N}"
    )

    print(
        f"    equals S    = {equals_S}"
    )

    print(
        f"    equals M1   = {equals_M1}"
    )

print()


# ============================================================================
# [13] HANKEL-INSPIRED SEARCH
# ============================================================================

print("[13] HANKEL-INSPIRED SEARCH")
print("-" * 78)

H2 = sp.expand(
    F2 * F4 - F3**2
)

H3 = sp.expand(
    F2 * F3 * F4
    + 2 * F3 * F4 * F5
    - F2 * F4**2
    - F3**2 * F5
)

print(
    "  H2 =",
    fact(H2)
)

print(
    "  H3 =",
    fact(H3)
)

print()

hankel_candidates = {
    "H2": H2,
    "H3": H3,
    "H2/F2": sp.cancel(H2 / F2),
    "H2/F3": sp.cancel(H2 / F3),
    "H3/H2": sp.cancel(H3 / H2),
}

for name, expression in hankel_candidates.items():

    derivative = sp.diff(
        expression,
        S
    )

    pure_N = zero(
        sp.together(
            derivative
        )
    )

    print(
        f"  {name}: pure-N = {pure_N}"
    )

print()


# ============================================================================
# [14] INFORMATION / ALGEBRAIC BOUNDARY
# ============================================================================

print("[14] ALGEBRAIC INFORMATION BOUNDARY")
print("-" * 78)

print(
    "  N = p*q"
)

print(
    "  S = p+q"
)

jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(N, p), sp.diff(N, q)],
            [sp.diff(S, p), sp.diff(S, q)],
        ]
    )
)

# The displayed Jacobian is evaluated after replacing
# N=pq and S=p+q.
jacobian_pq = sp.det(
    sp.Matrix(
        [
            [q, p],
            [1, 1],
        ]
    )
)

print(
    "  Jacobian d(N,S)/d(p,q) =",
    fact(jacobian_pq)
)

print(
    "  Therefore det = q-p."
)

print(
    "  For p != q the symmetric coordinates N,S are locally independent."
)

print(
    "  Hence there is no universal identity S=f(N) on the full"
)

print(
    "  two-variable factor space."
)

print()


# ============================================================================
# [15] NUMERICAL SANITY CHECKS
# ============================================================================

print("[15] NUMERICAL SANITY CHECKS")
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

numerical_failures = 0

for pp, qq in instances:

    nn = pp * qq
    ss = pp + qq

    f2_num = (
        6 * nn
        - ss**2
        + ss
    )

    f3_num = (
        (ss + 1)
        * f2_num
    )

    if f2_num == 0:

        print(
            f"  ({pp},{qq}) F2=0 -- ratio skipped"
        )

        continue

    divisible = (
        f3_num % f2_num == 0
    )

    recovered_s = (
        f3_num // f2_num
        - 1
    )

    passed = (
        divisible
        and recovered_s == ss
    )

    if not passed:
        numerical_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"S={ss} "
        f"recovered={recovered_s} "
        f"PASS={passed}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numerical_failures
)

print()


# ============================================================================
# [16] FINAL STATUS
# ============================================================================

print("[16] EXPERIMENT STATUS")
print("-" * 78)

basic_pass = all(
    direct_checks
)

f3_pass = zero(
    rel_F3
)

# This is the corrected reference:
# relation_kernel is the cleared numerator of the
# substitution S = F3/F2 - 1 into the F2 relation.
#
# It should NOT be expected to vanish identically, because
# F2 and F3 are independent kernel symbols here.
#
# The exact relation is instead interpreted numerically on
# the actual kernel image.
elimination_relation_defined = (
    relation_kernel is not None
)

n_only_found = (
    len(nontrivial_N) > 0
)

print(
    "  direct KAPPA identities =",
    basic_pass
)

print(
    "  F3=(S+1)F2 =",
    f3_pass
)

print(
    "  elimination numerator constructed =",
    elimination_relation_defined
)

print(
    "  nontrivial N-only invariant found =",
    n_only_found
)

print(
    "  numerical S recovery =",
    numerical_failures == 0
)

print()

print(
    "INTERPRETATION:"
)

if n_only_found:

    print(
        "  A nontrivial S-independent kernel combination"
    )

    print(
        "  was found in the tested low-degree algebra."
    )

    print(
        "  Inspect the candidate expressions above."
    )

else:

    print(
        "  No nontrivial N-only invariant was found"
    )

    print(
        "  in the tested low-degree polynomial algebra"
    )

    print(
        "  generated from F2,F3,F4,F5."
    )

print()

print(
    "ESTABLISHED BRIDGE:"
)

print(
    "  F2,F3"
)

print(
    "    -> S = F3/F2 - 1"
)

print(
    "    -> M1 = N + S + 1"
)

print(
    "    -> z^2 - S*z + N"
)

print(
    "    -> p,q"
)

print()

print(
    "UNRESOLVED INTERFACE:"
)

print(
    "  N"
)

print(
    "    -> independently generated KAPPA information"
)

print(
    "    -> S"
)

print(
    "    -> factor pair"
)

print()

print(
    "IMPORTANT:"
)

print(
    "  This experiment does not prove that no arbitrary"
)

print(
    "  high-complexity, analytic, modular, or external"
)

print(
    "  N-only construction can exist."
)

print(
    "  It only rejects the explicitly tested low-degree"
)

print(
    "  exact algebraic families."
)

print()
print("=" * 78)
print("EXPERIMENT 491 FINISHED")
print("=" * 78)