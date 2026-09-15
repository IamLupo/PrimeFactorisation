#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 494 START")
print("=" * 78)
print("SHIFTED-KERNEL OPERATOR EQUATION / SPECTRAL-PARAMETER SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S, x, z = sp.symbols("N S x z")


def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def together(expr):
    return sp.cancel(sp.together(expr))


# ============================================================================
# SYMMETRIC SUBSTITUTION
# ============================================================================

S_pq = p + q
N_pq = p * q


# ============================================================================
# KAPPA / SHIFTED KERNEL
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

    # F is divisible by x for all k < ell.
    quotient = sp.cancel(F / x)

    return sp.expand(quotient)


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(expr):
    """
    Replace symmetric expressions in p,q by N=pq and S=p+q.

    We use SymPy's symmetric reduction rather than ad-hoc substitution.
    """
    expr = sp.expand(expr)

    result = sp.symmetrize(
        expr,
        [p, q],
        formal=True
    )

    reduced, remainder, mapping = result

    if remainder != 0:
        raise ValueError(
            f"Expression is not symmetric: remainder={remainder}"
        )

    # mapping looks like [(s1, p+q), (s2, p*q)]
    s1, s2 = mapping[0][0], mapping[1][0]

    reduced = reduced.subs(
        {
            s1: S,
            s2: N
        }
    )

    return sp.expand(reduced)


# ============================================================================
# [1] BUILD SHIFTED KERNEL FAMILY
# ============================================================================

print("[1] BUILD k=1 SHIFTED KERNEL FAMILY")
print("-" * 78)

H = {}

for ell in range(2, 11):

    H[ell] = H_shifted(1, ell)

    H_ns = symmetric_reduce(H[ell])

    H[ell] = H_ns

    print(
        f"  ell={ell}:"
    )

    print(
        f"    H(x) = {H[ell]}"
    )

print()


# ============================================================================
# [2] BASIC STRUCTURAL AUDIT
# ============================================================================

print("[2] BASIC STRUCTURAL AUDIT")
print("-" * 78)

failures = 0

for ell in range(2, 11):

    expr = H[ell]

    # Degree must be ell-1.
    degree = sp.Poly(
        expr,
        x
    ).degree()

    expected = ell - 1

    degree_pass = degree == expected

    # Highest coefficient must be S.
    top = sp.Poly(
        expr,
        x
    ).coeff_monomial(
        x**(ell - 1)
    )

    top_pass = zero(top - S)

    # Restore F=xH and compare to p,q definition.
    original = F_shifted(1, ell)

    reconstructed = sp.expand(
        x * expr
    )

    reconstruction = symmetric_reduce(
        sp.expand(
            reconstructed - original
        )
    )

    reconstruction_pass = zero(
        reconstruction
    )

    passed = (
        degree_pass
        and top_pass
        and reconstruction_pass
    )

    print(
        f"  ell={ell}: "
        f"degree PASS={degree_pass}, "
        f"top PASS={top_pass}, "
        f"reconstruction PASS={reconstruction_pass}"
    )

    if not passed:
        failures += 1

print(
    "  STRUCTURAL FAILURES =",
    failures
)

print()


# ============================================================================
# [3] OPERATOR FAMILY
# ============================================================================

print("[3] SECOND-ORDER OPERATOR SEARCH")
print("-" * 78)

print(
    "Search for operators of the form"
)

print(
    "  A2(x,N) D^2 + A1(x,N) D + A0(x,N)"
)

print(
    "whose coefficients contain N but NOT S,"
)

print(
    "and for which"
)

print(
    "  L[H_(1,ell)] = 0"
)

print(
    "or"
)

print(
    "  L[H_(1,ell)] = lambda(S,N) H_(1,ell)"
)

print(
    "over several independent ell."
)

print()


D = lambda expr, r=1: sp.diff(
    expr,
    x,
    r
)


# ============================================================================
# [4] MONOMIAL N-ONLY OPERATOR BASIS
# ============================================================================

print("[4] N-ONLY OPERATOR BASIS")
print("-" * 78)

# A small but meaningful basis:
#
# A_j(x,N) consists of x^a N^b with
# a <= 2
# b <= 2
#
# We search operators whose coefficient polynomials are built
# from this common basis.

operator_basis = []

for a in range(3):
    for b in range(3):
        operator_basis.append(
            x**a * N**b
        )

print(
    "  coefficient basis size =",
    len(operator_basis)
)

print()


# ============================================================================
# [5] ANNIHILATING OPERATOR SEARCH
# ============================================================================

print("[5] COMMON ANNIHILATING OPERATOR SEARCH")
print("-" * 78)

# Candidate operator:
#
# L = A2 D2 + A1 D + A0
#
# Each Aj is an unknown linear combination of operator_basis.
#
# Search requires L(H)=0 for ell=2..10.

unknowns = []

A2_coeffs = sp.symbols(
    "a2_0:" + str(len(operator_basis))
)

A1_coeffs = sp.symbols(
    "a1_0:" + str(len(operator_basis))
)

A0_coeffs = sp.symbols(
    "a0_0:" + str(len(operator_basis))
)

unknowns.extend(A2_coeffs)
unknowns.extend(A1_coeffs)
unknowns.extend(A0_coeffs)


A2 = sum(
    c * b
    for c, b in zip(
        A2_coeffs,
        operator_basis
    )
)

A1 = sum(
    c * b
    for c, b in zip(
        A1_coeffs,
        operator_basis
    )
)

A0 = sum(
    c * b
    for c, b in zip(
        A0_coeffs,
        operator_basis
    )
)


equations = []


for ell in range(2, 11):

    h = H[ell]

    expression = sp.expand(
        A2 * D(h, 2)
        + A1 * D(h, 1)
        + A0 * h
    )

    poly = sp.Poly(
        expression,
        x,
        S,
        N
    )

    for coeff in poly.coeffs():

        equations.append(
            coeff
        )


print(
    "  unknown operator coefficients =",
    len(unknowns)
)

print(
    "  exact scalar equations =",
    len(equations)
)

# Convert equations to a linear matrix system.
M = sp.linear_eq_to_matrix(
    equations,
    unknowns
)[0]

print(
    "  operator matrix shape =",
    M.shape
)

nullspace = M.nullspace()

print(
    "  nullspace dimension =",
    len(nullspace)
)

print()


# ============================================================================
# [6] INTERPRET OPERATOR CANDIDATES
# ============================================================================

print("[6] OPERATOR CANDIDATE AUDIT")
print("-" * 78)

operator_hits = []

for idx, vec in enumerate(
    nullspace[:20],
    start=1
):

    A2_candidate = sp.expand(
        sum(
            vec[i] * operator_basis[i]
            for i in range(
                len(operator_basis)
            )
        )
    )

    offset = len(operator_basis)

    A1_candidate = sp.expand(
        sum(
            vec[offset + i]
            * operator_basis[i]
            for i in range(
                len(operator_basis)
            )
        )
    )

    offset *= 2

    A0_candidate = sp.expand(
        sum(
            vec[offset + i]
            * operator_basis[i]
            for i in range(
                len(operator_basis)
            )
        )
    )

    # Normalize by first nonzero coefficient.
    coeffs = list(vec)

    nonzero_coeffs = [
        c
        for c in coeffs
        if c != 0
    ]

    if nonzero_coeffs:

        scale = nonzero_coeffs[0]

        A2_candidate = sp.factor(
            A2_candidate / scale
        )

        A1_candidate = sp.factor(
            A1_candidate / scale
        )

        A0_candidate = sp.factor(
            A0_candidate / scale
        )

    print(
        f"  candidate {idx}:"
    )

    print(
        f"    A2 = {A2_candidate}"
    )

    print(
        f"    A1 = {A1_candidate}"
    )

    print(
        f"    A0 = {A0_candidate}"
    )

    # Test whether the operator is genuinely nonzero.
    nonzero_operator = not (
        zero(A2_candidate)
        and zero(A1_candidate)
        and zero(A0_candidate)
    )

    if nonzero_operator:

        operator_hits.append(
            (
                A2_candidate,
                A1_candidate,
                A0_candidate
            )
        )

print()


# ============================================================================
# [7] FIRST-ORDER OPERATOR SEARCH
# ============================================================================

print("[7] FIRST-ORDER N-ONLY OPERATOR SEARCH")
print("-" * 78)

unknown1 = sp.symbols(
    "b0:" + str(len(operator_basis))
)

unknown0 = sp.symbols(
    "c0:" + str(len(operator_basis))
)

A1_first = sum(
    c * b
    for c, b in zip(
        unknown1,
        operator_basis
    )
)

A0_first = sum(
    c * b
    for c, b in zip(
        unknown0,
        operator_basis
    )
)

eq_first = []

for ell in range(2, 11):

    expr = sp.expand(
        A1_first * D(H[ell], 1)
        + A0_first * H[ell]
    )

    poly = sp.Poly(
        expr,
        x,
        S,
        N
    )

    eq_first.extend(
        poly.coeffs()
    )

M_first = sp.linear_eq_to_matrix(
    eq_first,
    list(unknown1) + list(unknown0)
)[0]

ns_first = M_first.nullspace()

print(
    "  matrix shape =",
    M_first.shape
)

print(
    "  nullspace dimension =",
    len(ns_first)
)

print()


# ============================================================================
# [8] EIGENVALUE-TYPE OPERATOR SEARCH
# ============================================================================

print("[8] SPECTRAL-PARAMETER SEARCH")
print("-" * 78)

print(
    "Now test whether an N-only first-order operator"
)

print(
    "can have H_(1,ell) as an eigenfunction with"
)

print(
    "eigenvalue depending only on S and N."
)

print(
    "Candidate relation:"
)

print(
    "  A1(N,x) H'(x) + A0(N,x) H(x)"
)

print(
    "       = Lambda(N,S) H(x)"
)

print()


# Use simple scalar spectral ansatz:
#
# Lambda = l0(N) + l1(N) S
#
# where l0,l1 have degree <= 2 in N.

lambda_basis = [
    N**0,
    N,
    N**2
]

lambda0_coeffs = sp.symbols(
    "u0:3"
)

lambda1_coeffs = sp.symbols(
    "v0:3"
)

Lambda0 = sum(
    lambda0_coeffs[i] * lambda_basis[i]
    for i in range(3)
)

Lambda1 = sum(
    lambda1_coeffs[i] * lambda_basis[i]
    for i in range(3)
)

Lambda = (
    Lambda0
    + Lambda1*S
)


# Coefficient basis for A1 and A0:
small_basis = [
    1,
    x,
    x**2,
    N,
    N*x,
    N*x**2
]

alpha = sp.symbols(
    "alpha0:" + str(len(small_basis))
)

beta = sp.symbols(
    "beta0:" + str(len(small_basis))
)

A1_small = sum(
    alpha[i] * small_basis[i]
    for i in range(len(small_basis))
)

A0_small = sum(
    beta[i] * small_basis[i]
    for i in range(len(small_basis))
)

spec_unknowns = list(alpha) + list(beta)
spec_unknowns.extend(lambda0_coeffs)
spec_unknowns.extend(lambda1_coeffs)

spec_eq = []

for ell in range(2, 9):

    h = H[ell]

    residual = sp.expand(
        A1_small * D(h)
        + A0_small * h
        - Lambda * h
    )

    poly = sp.Poly(
        residual,
        x,
        S,
        N
    )

    spec_eq.extend(
        poly.coeffs()
    )

M_spec = sp.linear_eq_to_matrix(
    spec_eq,
    spec_unknowns
)[0]

ns_spec = M_spec.nullspace()

print(
    "  spectral matrix shape =",
    M_spec.shape
)

print(
    "  spectral nullspace dimension =",
    len(ns_spec)
)

print()


# ============================================================================
# [9] KNOWN TOP-COEFFICIENT CONNECTION
# ============================================================================

print("[9] TOP-COEFFICIENT CHECK")
print("-" * 78)

for ell in range(2, 11):

    top = sp.Poly(
        H[ell],
        x
    ).coeff_monomial(
        x**(ell - 1)
    )

    print(
        f"  ell={ell}: top={top}, PASS={zero(top-S)}"
    )

print()


# ============================================================================
# [10] DIRECT CONNECTION TO SHIFTED PRODUCT
# ============================================================================

print("[10] SHIFTED-PRODUCT COMPARISON")
print("-" * 78)

for ell in range(2, 7):

    h2 = H[ell]

    # Compare top coefficient with S.
    top = sp.Poly(
        h2,
        x
    ).coeff_monomial(
        x**(ell - 1)
    )

    # Compare against M_x = N + S*x + x^2.
    Mx = N + S*x + x**2

    # Search whether H can be written as a low-degree polynomial
    # in Mx with coefficients independent of S.
    a0, a1, a2 = sp.symbols(
        f"a0_{ell} a1_{ell} a2_{ell}"
    )

    candidate = (
        a2 * Mx**2
        + a1 * Mx
        + a0
    )

    diff = sp.expand(
        h2 - candidate
    )

    equations_mx = sp.Poly(
        diff,
        x,
        S,
        N
    ).coeffs()

    if equations_mx:

        matrix_mx = sp.linear_eq_to_matrix(
            equations_mx,
            [a0, a1, a2]
        )[0]

        ns_mx = matrix_mx.nullspace()

    else:

        ns_mx = []

    print(
        f"  ell={ell}: "
        f"top=S PASS={zero(top-S)}, "
        f"quadratic-Mx representations={len(ns_mx)}"
    )

print()


# ============================================================================
# [11] NUMERICAL SANITY CHECKS
# ============================================================================

print("[11] NUMERICAL SANITY CHECKS")
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

    local_fail = False

    for ell in range(2, 7):

        expr = H[ell].subs(
            {
                N: nn,
                S: ss
            }
        )

        top = sp.Poly(
            expr,
            x
        ).coeff_monomial(
            x**(ell - 1)
        )

        if top != ss:

            local_fail = True

    if local_fail:

        numeric_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"top-coefficient PASS={not local_fail}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [12] FINAL STATUS
# ============================================================================

print("[12] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  shifted-kernel construction =",
    failures == 0
)

print(
    "  second-order annihilating operators found =",
    len(operator_hits)
)

print(
    "  first-order nullspace dimension =",
    len(ns_first)
)

print(
    "  spectral nullspace dimension =",
    len(ns_spec)
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
    "This experiment deliberately attacks the source interface at the"
)

print(
    "operator level rather than by adding another finite collection of"
)

print(
    "algebraic feature combinations."
)

print()

if len(operator_hits) == 0:

    print(
        "No common low-complexity second-order N-only annihilating"
    )

    print(
        "operator was found for H_(1,ell), ell=2..10."
    )

else:

    print(
        "At least one common N-only second-order operator exists."
    )

print()

print(
    "The decisive unresolved interface remains:"
)

print(
    "  N"
)

print(
    "   -> independently generated shifted/original kernel"
)

print(
    "   -> spectral/operator extraction"
)

print(
    "   -> S"
)

print(
    "   -> p,q"
)

print()

print(
    "A surviving operator would be a qualitatively new bridge."
)

print(
    "A null result would further support the conclusion that the"
)

print(
    "missing information must enter before the shifted kernel is formed."
)

print()


print("=" * 78)
print("EXPERIMENT 494 FINISHED")
print("=" * 78)
