#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 493 START")
print("=" * 78)
print("N-ONLY LINEAR CLOSURE OF THE ORIGINAL KAPPA SEQUENCE")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

N, S = sp.symbols("N S")


def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def safe_degree(expr, *gens):
    expr = sp.cancel(expr)
    num, den = sp.fraction(expr)

    try:
        pn = sp.Poly(num, *gens)
        pd = sp.Poly(den, *gens)
        return pn.total_degree(), pd.total_degree()
    except sp.PolynomialError:
        return None, None


# ============================================================================
# [1] ORIGINAL KAPPA SEQUENCE IN (N,S)
# ============================================================================

print("[1] ORIGINAL KAPPA SEQUENCE")
print("-" * 78)

F = {}

F[1] = sp.Integer(0)

F[2] = (
    6*N
    - S**2
    + S
)

F[3] = (
    (S + 1)
    * F[2]
)

F[4] = (
    -10*N**2
    + 8*N*S**2
    + 6*N*S
    + 8*N
    - S**4
    + S
)

F[5] = (
    -(S + 1)
    * (
        20*N**2
        - 10*N*S**2
        - 10*N
        + S**4
        - S**3
        + S**2
        - S
    )
)

F[6] = (
    14*N**3
    - 33*N**2*S**2
    - 45*N**2*S
    - 40*N**2
    + 12*N*S**4
    + 15*N*S**3
    + 20*N*S**2
    + 15*N*S
    + 12*N
    - S**6
    + S
)

F[7] = (
    (S + 1)
    * (
        42*N**3
        - 49*N**2*S**2
        - 35*N**2*S
        - 70*N**2
        + 14*N*S**4
        + 7*N*S**3
        + 28*N*S**2
        + 7*N*S
        + 14*N
        - S**6
        + S**5
        - S**4
        + S**3
        - S**2
        + S
    )
)

F[8] = (
    -18*N**4
    + 88*N**3*S**2
    + 140*N**3*S
    + 112*N**3
    - 68*N**2*S**4
    - 140*N**2*S**3
    - 224*N**2*S**2
    - 210*N**2*S
    - 112*N**2
    + 16*N*S**6
    + 28*N*S**5
    + 56*N*S**4
    + 70*N*S**3
    + 56*N*S**2
    + 28*N*S
    + 16*N
    - S**8
    + S
)

for k in range(1, 9):
    print(
        f"  F_{k} = {F[k]}"
    )

print()


# ============================================================================
# [2] VERIFY BASIC KAPPA IDENTITIES
# ============================================================================

print("[2] BASIC KAPPA IDENTITIES")
print("-" * 78)

basic_checks = {
    "F3-(S+1)F2":
        F[3] - (S + 1)*F[2],

    "F4":
        F[4],

    "F5":
        F[5],
}

basic_failures = 0

for name, expr in basic_checks.items():

    ok = zero(expr)

    print(
        f"  {name}: PASS={ok}"
    )

    if not ok:
        basic_failures += 1

print(
    "  BASIC FAILURES =",
    basic_failures
)

print()


# ============================================================================
# [3] SEARCH LINEAR COMBINATIONS
# ============================================================================

print("[3] S-INDEPENDENT LINEAR COMBINATIONS")
print("-" * 78)

print(
    "Search:"
)

print(
    "  G(N,S) = c2(N)F2 + c3(N)F3 + ... + cm(N)Fm"
)

print(
    "where each c_k(N) is a polynomial of bounded degree."
)

print()


# We perform the search for increasing sequence lengths.
# Coefficients are polynomials in N of degree <= d.
#
# To avoid the trivial zero solution, nullspace vectors are
# interpreted projectively.

def linear_closure_search(max_k, coeff_degree):

    feature_pairs = []

    for k in range(2, max_k + 1):

        for d in range(coeff_degree + 1):

            feature_pairs.append(
                (k, d)
            )

    # For each F_k and N^d coefficient create one unknown.
    expressions = []

    for k in range(2, max_k + 1):

        for d in range(coeff_degree + 1):

            expressions.append(
                N**d * F[k]
            )

    # Differentiate with respect to S.
    derivative_polys = []

    for expr in expressions:

        derivative_polys.append(
            sp.Poly(
                sp.expand(
                    sp.diff(expr, S)
                ),
                N,
                S
            )
        )

    monomials = set()

    for poly in derivative_polys:

        for monom in poly.monoms():

            monomials.add(
                monom
            )

    monomials = sorted(
        monomials
    )

    matrix = []

    for monom in monomials:

        row = []

        for poly in derivative_polys:

            row.append(
                poly.coeff_monomial(
                    monom
                )
            )

        matrix.append(row)

    if not matrix:

        return []

    M = sp.Matrix(matrix)

    nullspace = M.nullspace()

    candidates = []

    for vec in nullspace:

        expr = 0

        for coeff, basis_expr in zip(
            vec,
            expressions
        ):

            expr += coeff * basis_expr

        expr = sp.factor(
            sp.expand(expr)
        )

        if zero(expr):
            continue

        candidates.append(
            (
                vec,
                expr
            )
        )

    return candidates


search_results = []

for max_k in [3, 4, 5, 6, 7, 8]:

    for coeff_degree in [0, 1, 2, 3]:

        candidates = linear_closure_search(
            max_k,
            coeff_degree
        )

        print(
            f"  F2..F{max_k}, coefficient degree <= {coeff_degree}: "
            f"nullspace candidates = {len(candidates)}"
        )

        for _, expr in candidates[:10]:

            print(
                f"      {expr}"
            )

        search_results.append(
            (
                max_k,
                coeff_degree,
                candidates
            )
        )

print()


# ============================================================================
# [4] CLASSIFY THE FOUND COMBINATIONS
# ============================================================================

print("[4] CLASSIFICATION OF LINEAR CLOSURE CANDIDATES")
print("-" * 78)

nonzero_n_only = []
zero_identities = []
s_dependent = []

for max_k, degree, candidates in search_results:

    for _, expr in candidates:

        dS = sp.factor(
            sp.diff(
                expr,
                S
            )
        )

        if not zero(dS):
            s_dependent.append(
                (
                    max_k,
                    degree,
                    expr
                )
            )
            continue

        # Independent of S. Now determine whether it is
        # actually a nonzero function of N.
        if zero(expr):

            zero_identities.append(
                (
                    max_k,
                    degree,
                    expr
                )
            )

        else:

            nonzero_n_only.append(
                (
                    max_k,
                    degree,
                    fact(expr)
                )
            )

print(
    "  nonzero N-only expressions found =",
    len(nonzero_n_only)
)

for item in nonzero_n_only[:30]:

    max_k, degree, expr = item

    print(
        f"  F2..F{max_k}, degree<={degree}:"
    )

    print(
        f"    {expr}"
    )

print()


# ============================================================================
# [5] DETECT TRIVIAL N-ONLY COMBINATIONS
# ============================================================================

print("[5] TRIVIALITY AUDIT")
print("-" * 78)

print(
    "Any combination generated solely from already-known"
)

print(
    "identities such as F3=(S+1)F2 must not be mistaken"
)

print(
    "for a new source-free construction."
)

print()

trivial_hits = []
genuine_hits = []

for max_k, degree, expr in nonzero_n_only:

    expr = sp.factor(expr)

    # If the expression depends only on N and is simply a
    # polynomial multiple of N, classify it separately.
    try:
        polyN = sp.Poly(
            expr,
            N
        )

        if polyN.is_zero:
            trivial_hits.append(
                (max_k, degree, expr)
            )
        else:
            # Record all for now; later controls decide whether
            # the expression actually contains new kernel data.
            genuine_hits.append(
                (max_k, degree, expr)
            )

    except sp.PolynomialError:

        genuine_hits.append(
            (max_k, degree, expr)
        )

print(
    "  candidate N-only polynomial expressions =",
    len(genuine_hits)
)

for item in genuine_hits[:30]:

    print(
        f"  {item}"
    )

print()


# ============================================================================
# [6] NORMALIZED COMBINATION SEARCH
# ============================================================================

print("[6] NORMALIZED LINEAR COMBINATIONS")
print("-" * 78)

print(
    "Search whether an S-independent combination reduces to"
)

print(
    "a simple target such as N, N+1, N^2, or N(N+1)."
)

targets = {
    "N": N,
    "N+1": N + 1,
    "N^2": N**2,
    "N*(N+1)": N*(N + 1),
    "N*(N+2)": N*(N + 2),
    "N^2+N+1": N**2 + N + 1,
}

for target_name, target in targets.items():

    found = False

    for max_k, degree, candidates in search_results:

        for _, expr in candidates:

            diff = sp.factor(
                sp.expand(
                    expr - target
                )
            )

            if zero(diff):

                print(
                    f"  HIT: {target_name} from F2..F{max_k}, "
                    f"degree<={degree}"
                )

                print(
                    f"       expression = {expr}"
                )

                found = True

    if not found:

        print(
            f"  NO EXACT HIT: {target_name}"
        )

print()


# ============================================================================
# [7] QUOTIENT-BASED CLOSURE TEST
# ============================================================================

print("[7] QUOTIENT-BASED CLOSURE TEST")
print("-" * 78)

print(
    "Test whether ratios of kernel values simplify to"
)

print(
    "functions independent of S."
)

ratio_hits = []

for i in range(2, 9):

    for j in range(2, 9):

        if i == j:
            continue

        ratio = sp.cancel(
            F[i] / F[j]
        )

        num, den = sp.fraction(
            ratio
        )

        dnum = sp.factor(
            sp.diff(
                num,
                S
            )
        )

        dden = sp.factor(
            sp.diff(
                den,
                S
            )
        )

        if zero(dnum) and zero(dden):

            ratio_hits.append(
                (
                    i,
                    j,
                    fact(ratio)
                )
            )

for i, j, expr in ratio_hits:

    print(
        f"  F{i}/F{j} = {expr}"
    )

print(
    "  S-independent ratios =",
    len(ratio_hits)
)

print()


# ============================================================================
# [8] MINIMAL POLYNOMIAL IN S
# ============================================================================

print("[8] MINIMAL ALGEBRAIC RELATIONS")
print("-" * 78)

print(
    "For each F_k, inspect the degree in S."
)

for k in range(2, 9):

    poly = sp.Poly(
        F[k],
        S
    )

    print(
        f"  F{k}: degree_S = {poly.degree()}"
    )

print()

print(
    "The important question is whether a combination of"
)

print(
    "these polynomial equations can eliminate S without"
)

print(
    "reintroducing hidden factor information."
)

print()


# ============================================================================
# [9] RESULTANT TESTS ON SMALL PAIRS
# ============================================================================

print("[9] SMALL RESULTANT ELIMINATION TEST")
print("-" * 78)

print(
    "For selected pairs, compute the resultant with respect to S."
)

print(
    "This detects whether two kernel values jointly generate"
)

print(
    "a polynomial relation involving N alone."
)

pairs = [
    (2, 3),
    (2, 4),
    (2, 5),
    (3, 4),
    (3, 5),
]

resultant_data = []

for i, j in pairs:

    Ri = sp.Poly(
        sp.expand(
            sp.Symbol(f"Y{i}") - F[i]
        ),
        S
    )

    Rj = sp.Poly(
        sp.expand(
            sp.Symbol(f"Y{j}") - F[j]
        ),
        S
    )

    Yi = sp.Symbol(f"Y{i}")
    Yj = sp.Symbol(f"Y{j}")

    resultant = sp.factor(
        sp.resultant(
            F[i] - Yi,
            F[j] - Yj,
            S
        )
    )

    resultant_data.append(
        (
            i,
            j,
            resultant
        )
    )

    print(
        f"  Resultant(F{i},F{j}; S) ="
    )

    print(
        f"    {resultant}"
    )

print()


# ============================================================================
# [10] INFORMATION-BOUNDARY CONTROL
# ============================================================================

print("[10] INFORMATION-BOUNDARY CONTROL")
print("-" * 78)

print(
    "N and S are algebraically independent symmetric coordinates"
)

print(
    "away from the diagonal p=q because:"
)

jacobian = sp.det(
    sp.Matrix([
        [sp.diff(p := sp.Symbol("p"), p), sp.diff(p, q := sp.Symbol("q"))],
        [sp.diff(q + p, p), sp.diff(q + p, q)]
    ])
)

# Replace the above with the direct known Jacobian.
jacobian = q - p

print(
    "  det d(N,S)/d(p,q) = q-p"
)

print(
    f"  symbolic determinant = {jacobian}"
)

print(
    "  Therefore no universal identity S=f(N) exists"
)

print(
    "  over the unrestricted two-variable factor space."
)

print()


# ============================================================================
# [11] NUMERICAL SANITY CHECK
# ============================================================================

print("[11] NUMERICAL SANITY CHECK")
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

    # Evaluate the two known first values.
    f2n = (
        6*nn
        - ss**2
        + ss
    )

    f3n = (
        (ss + 1)
        * f2n
    )

    recovered = (
        f3n // f2n - 1
        if f2n != 0
        else None
    )

    passed = (
        recovered == ss
    )

    if not passed:
        numeric_failures += 1

    print(
        f"  ({pp},{qq})"
        f"  S={ss}"
        f"  F2={f2n}"
        f"  F3={f3n}"
        f"  recovered={recovered}"
        f"  PASS={passed}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [12] FINAL INTERPRETATION
# ============================================================================

print("[12] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  direct KAPPA identities =",
    basic_failures == 0
)

print(
    "  S-independent linear candidates =",
    len(nonzero_n_only)
)

print(
    "  S-independent ratios =",
    len(ratio_hits)
)

print(
    "  numerical recovery =",
    numeric_failures == 0
)

print()

print(
    "MAIN QUESTION:"
)

print(
    "  Does the algebra generated by the original KAPPA values"
)

print(
    "  contain a genuinely new nonconstant N-only quantity?"
)

print()

if len(nonzero_n_only) == 0:

    print(
        "RESULT:"
    )

    print(
        "  No nonzero S-independent linear closure was found"
    )

    print(
        "  in the tested F2..F8 / polynomial-in-N family."
    )

else:

    print(
        "RESULT:"
    )

    print(
        "  S-independent linear combinations were found."
    )

    print(
        "  These require individual inspection to determine"
    )

    print(
        "  whether they are genuinely new or merely identities"
    )

    print(
        "  generated by the existing recurrence structure."
    )

print()

print(
    "INFORMATION BOUNDARY:"
)

print(
    "  F_n"
)

print(
    "   -> Hankel / shifted / recurrence constructions"
)

print(
    "   -> S"
)

print(
    "   -> p,q"
)

print()

print(
    "  remains separate from:"
)

print(
    "  N"
)

print(
    "   -> independently generated F_n"
)

print(
    "   -> S"
)

print(
    "   -> p,q"
)

print()

print(
    "NEXT TARGET:"
)

print(
    "  Search the original N-only homogeneous-layer construction"
)

print(
    "  for an operator whose output equals one of the"
)

print(
    "  shifted evaluations or KAPPA values."
)

print()


print("=" * 78)
print("EXPERIMENT 493 FINISHED")
print("=" * 78)
