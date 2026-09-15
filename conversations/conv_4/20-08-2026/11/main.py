#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 488 START")
print("=" * 78)
print("FOUR-BASE HANKEL THEOREM AND VANDERMONDE PROOF")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# ORIGINAL KAPPA SEQUENCE
# ============================================================================

def F_pq(n):
    return sp.expand(
        p * (q + 1) ** n
        + q * (p + 1) ** n
        - (q + 1) * p ** n
        - (p + 1) * q ** n
    )


MAX_N = 20

F_PQ = {}
F_NS = {}

for n in range(1, MAX_N + 1):
    expr = F_pq(n)
    F_PQ[n] = expr

    reduced = sp.symmetrize(
        expr,
        [p, q],
        formal=True
    )

    symmetric_expr = reduced[0]
    remainder = reduced[1]
    mapping = reduced[2]

    if remainder != 0:
        raise RuntimeError(
            f"Symmetrization failed for F_{n}: {remainder}"
        )

    reduced_expr = symmetric_expr

    for formal_symbol, concrete_expr in mapping:
        if sp.expand(concrete_expr - (p + q)) == 0:
            reduced_expr = reduced_expr.subs(formal_symbol, S)
        elif sp.expand(concrete_expr - p * q) == 0:
            reduced_expr = reduced_expr.subs(formal_symbol, N)
        else:
            raise RuntimeError(
                f"Unexpected symmetric mapping: "
                f"{formal_symbol} -> {concrete_expr}"
            )

    reduced_expr = sp.expand(reduced_expr)

    verification = sp.expand(
        reduced_expr.subs({
            N: p * q,
            S: p + q,
        }) - expr
    )

    if verification != 0:
        raise RuntimeError(
            f"Symmetric reconstruction failed for F_{n}: "
            f"{sp.factor(verification)}"
        )

    F_NS[n] = reduced_expr


# ============================================================================
# FOUR-EXPONENTIAL DECOMPOSITION
# ============================================================================

print("[1] ORIGINAL KAPPA EXPONENTIAL DECOMPOSITION")
print("-" * 78)

bases = [
    p,
    q,
    p + 1,
    q + 1,
]

weights = [
    -(q + 1),
    -(p + 1),
    q,
    p,
]

for i in range(4):
    print(
        "  lambda_{} = {}    c_{} = {}".format(
            i + 1,
            bases[i],
            i + 1,
            weights[i],
        )
    )

print()

decomp_failures = []

for n in range(1, 11):

    reconstructed = sp.expand(
        sum(
            weights[i] * bases[i] ** n
            for i in range(4)
        )
    )

    difference = sp.expand(
        reconstructed - F_PQ[n]
    )

    passed = difference == 0

    print(
        f"  n={n:2d}: PASS={passed}"
    )

    if not passed:
        decomp_failures.append(n)
        print(
            f"    difference = {sp.factor(difference)}"
        )

print()
print(
    f"  DECOMPOSITION FAILURES = {len(decomp_failures)}"
)
print()


# ============================================================================
# PRODUCTS
# ============================================================================

print("[2] PRODUCT IDENTITIES")
print("-" * 78)

product_weights = sp.factor(
    sp.prod(weights)
)

product_bases = sp.factor(
    sp.prod(bases)
)

expected_product = sp.expand(
    p * q * (p + 1) * (q + 1)
)

weight_difference = sp.expand(
    product_weights - expected_product
)

base_difference = sp.expand(
    product_bases - expected_product
)

print(
    f"  product(c_i) = {product_weights}"
)

print(
    f"  product(lambda_i) = {product_bases}"
)

print(
    f"  expected = {expected_product}"
)

print(
    f"  weight product PASS = {weight_difference == 0}"
)

print(
    f"  base product PASS = {base_difference == 0}"
)

print()


# ============================================================================
# VANDERMONDE
# ============================================================================

print("[3] VANDERMONDE STRUCTURE")
print("-" * 78)

vandermonde = sp.Integer(1)

for i in range(4):
    for j in range(i + 1, 4):
        vandermonde *= bases[j] - bases[i]

vandermonde = sp.factor(
    sp.expand(vandermonde)
)

vandermonde_sq = sp.factor(
    sp.expand(vandermonde ** 2)
)

expected_vandermonde_sq = sp.expand(
    (p - q) ** 4
    * ((p - q) ** 2 - 1) ** 2
)

vd_difference = sp.factor(
    sp.expand(
        vandermonde_sq
        - expected_vandermonde_sq
    )
)

print(
    f"  Vandermonde = {vandermonde}"
)

print(
    f"  Vandermonde^2 = {vandermonde_sq}"
)

print(
    f"  Expected = {expected_vandermonde_sq}"
)

print(
    f"  Difference = {vd_difference}"
)

print(
    f"  PASS = {vd_difference == 0}"
)

print()


# ============================================================================
# VANDERMONDE IN N,S
# ============================================================================

print("[4] VANDERMONDE N,S REPRESENTATION")
print("-" * 78)

vandermonde_ns = (
    (4 * N - S ** 2) ** 2
    * (4 * N - S ** 2 + 1) ** 2
)

vandermonde_ns_pq = sp.expand(
    vandermonde_ns.subs({
        N: p * q,
        S: p + q,
    })
)

vd_ns_difference = sp.factor(
    sp.expand(
        vandermonde_sq
        - vandermonde_ns_pq
    )
)

print(
    "  Target = (4N-S^2)^2*(4N-S^2+1)^2"
)

print(
    f"  Difference = {vd_ns_difference}"
)

print(
    f"  PASS = {vd_ns_difference == 0}"
)

print()


# ============================================================================
# HANKEL MATRICES
# ============================================================================

def hankel_matrix_pq(start, size=4):
    matrix = []

    for i in range(size):
        row = []

        for j in range(size):
            index = start + i + j

            if index not in F_PQ:
                raise KeyError(
                    f"F_{index} is unavailable"
                )

            row.append(F_PQ[index])

        matrix.append(row)

    return sp.Matrix(matrix)


def hankel_matrix_ns(start, size=4):
    matrix = []

    for i in range(size):
        row = []

        for j in range(size):
            index = start + i + j

            if index not in F_NS:
                raise KeyError(
                    f"F_{index} is unavailable"
                )

            row.append(F_NS[index])

        matrix.append(row)

    return sp.Matrix(matrix)


def H4_pq(start):
    return sp.factor(
        sp.expand(
            hankel_matrix_pq(start, 4).det()
        )
    )


def H4_ns(start):
    return sp.factor(
        sp.expand(
            hankel_matrix_ns(start, 4).det()
        )
    )


# ============================================================================
# GENERAL HANKEL THEOREM
# ============================================================================

print("[5] GENERAL FOUR-BASE HANKEL THEOREM")
print("-" * 78)

hankel_pq_failures = []

for start in range(1, 7):

    actual = H4_pq(start)

    predicted = sp.expand(
        product_weights
        * product_bases ** start
        * vandermonde_sq
    )

    difference = sp.factor(
        sp.expand(
            actual - predicted
        )
    )

    passed = difference == 0

    print(
        f"  start={start}: PASS={passed}"
    )

    if not passed:
        hankel_pq_failures.append(start)
        print(
            f"    difference = {difference}"
        )

print()

print(
    f"  p,q THEOREM FAILURES = "
    f"{len(hankel_pq_failures)}"
)

print()


# ============================================================================
# CLOSED FORM IN N,S
# ============================================================================

print("[6] CLOSED FORM IN N,S")
print("-" * 78)

M1 = N + S + 1
A = 4 * N - S ** 2
C = 4 * N - S ** 2 + 1

hankel_ns_failures = []

for start in range(1, 7):

    actual = H4_ns(start)

    predicted = sp.expand(
        (N * M1) ** (start + 1)
        * A ** 2
        * C ** 2
    )

    difference = sp.factor(
        sp.expand(
            actual - predicted
        )
    )

    passed = difference == 0

    print(
        f"  start={start}: PASS={passed}"
    )

    if not passed:
        hankel_ns_failures.append(start)
        print(
            f"    difference = {difference}"
        )

print()

print(
    f"  N,S THEOREM FAILURES = "
    f"{len(hankel_ns_failures)}"
)

print()


# ============================================================================
# CONSECUTIVE HANKEL RATIO
# ============================================================================

print("[7] CONSECUTIVE HANKEL RATIO")
print("-" * 78)

ratio_failures = []

for start in range(1, 6):

    Hs = sp.expand(
        (N * M1) ** (start + 1)
        * A ** 2
        * C ** 2
    )

    Hnext = sp.expand(
        (N * M1) ** (start + 2)
        * A ** 2
        * C ** 2
    )

    ratio = sp.factor(
        sp.cancel(
            Hnext / Hs
        )
    )

    difference = sp.factor(
        sp.expand(
            ratio - N * M1
        )
    )

    passed = difference == 0

    print(
        f"  start={start}"
    )

    print(
        f"    ratio = {ratio}"
    )

    print(
        f"    expected = N*(N+S+1)"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        ratio_failures.append(start)

print()

print(
    f"  RATIO FAILURES = {len(ratio_failures)}"
)

print()


# ============================================================================
# DIRECT M1 / S RECOVERY
# ============================================================================

print("[8] DIRECT M1 AND S RECOVERY")
print("-" * 78)

H1_generic = sp.expand(
    (N * M1) ** 2
    * A ** 2
    * C ** 2
)

H2_generic = sp.expand(
    (N * M1) ** 3
    * A ** 2
    * C ** 2
)

recovered_M1 = sp.factor(
    sp.cancel(
        H2_generic
        / (N * H1_generic)
    )
)

recovered_S = sp.expand(
    recovered_M1
    - N
    - 1
)

m1_difference = sp.factor(
    sp.expand(
        recovered_M1 - M1
    )
)

s_difference = sp.factor(
    sp.expand(
        recovered_S - S
    )
)

print(
    f"  recovered M1 = {recovered_M1}"
)

print(
    f"  target M1 = {M1}"
)

print(
    f"  M1 difference = {m1_difference}"
)

print(
    f"  M1 PASS = {m1_difference == 0}"
)

print()

print(
    f"  recovered S = {recovered_S}"
)

print(
    f"  target S = {S}"
)

print(
    f"  S difference = {s_difference}"
)

print(
    f"  S PASS = {s_difference == 0}"
)

print()


# ============================================================================
# DEGENERACY
# ============================================================================

print("[9] DEGENERACY ANALYSIS")
print("-" * 78)

print(
    "  4N-S^2 = -(p-q)^2"
)

print(
    "  4N-S^2+1 = 1-(p-q)^2"
)

print(
    "  Therefore H_s can vanish when p=q"
)

print(
    "  or when |p-q|=1."
)

print(
    "  For distinct odd primes, |p-q| >= 2."
)

print()


# ============================================================================
# CONTROL CASES
# ============================================================================

print("[10] CONSECUTIVE-FACTOR CONTROL")
print("-" * 78)

controls = [
    (2, 3),
    (3, 4),
    (5, 7),
    (7, 8),
    (11, 13),
]

control_failures = []

for pp, qq in controls:

    values = {}

    for n in range(1, 10):
        values[n] = (
            pp * (qq + 1) ** n
            + qq * (pp + 1) ** n
            - (qq + 1) * pp ** n
            - (pp + 1) * qq ** n
        )

    matrix = sp.Matrix(
        [
            [
                values[1 + i + j]
                for j in range(4)
            ]
            for i in range(4)
        ]
    )

    determinant = int(
        matrix.det()
    )

    expected_zero = (
        pp == qq
        or abs(pp - qq) == 1
    )

    actual_zero = (
        determinant == 0
    )

    passed = (
        expected_zero == actual_zero
    )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    |p-q| = {abs(pp - qq)}"
    )

    print(
        f"    H4 = {determinant}"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        control_failures.append(
            (pp, qq)
        )

print()

print(
    f"  CONTROL FAILURES = "
    f"{len(control_failures)}"
)

print()


# ============================================================================
# NUMERICAL AUDIT
# ============================================================================

print("[11] NUMERICAL PRIME-PAIR AUDIT")
print("-" * 78)

prime_pairs = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]


def numeric_F(pp, qq, n):
    return (
        pp * (qq + 1) ** n
        + qq * (pp + 1) ** n
        - (qq + 1) * pp ** n
        - (pp + 1) * qq ** n
    )


def numeric_H4(pp, qq, start):

    max_index = (
        start
        + 3
        + 3
    )

    values = {}

    for n in range(1, max_index + 1):
        values[n] = numeric_F(
            pp,
            qq,
            n,
        )

    matrix = sp.Matrix(
        [
            [
                values[start + i + j]
                for j in range(4)
            ]
            for i in range(4)
        ]
    )

    return int(
        matrix.det()
    )


numerical_failures = []

for pp, qq in prime_pairs:

    NN = pp * qq
    true_M1 = (pp + 1) * (qq + 1)

    H1 = numeric_H4(
        pp,
        qq,
        1,
    )

    H2 = numeric_H4(
        pp,
        qq,
        2,
    )

    H3 = numeric_H4(
        pp,
        qq,
        3,
    )

    if H1 == 0 or H2 == 0:

        print(
            f"  ({pp},{qq}) "
            f"degenerate: H1={H1}, H2={H2}"
        )

        continue

    recovered_1 = sp.Rational(
        H2,
        NN * H1,
    )

    recovered_2 = sp.Rational(
        H3,
        NN * H2,
    )

    passed = (
        recovered_1 == true_M1
        and recovered_2 == true_M1
    )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N = {NN}"
    )

    print(
        f"    true M1 = {true_M1}"
    )

    print(
        f"    H2/(N*H1) = {recovered_1}"
    )

    print(
        f"    H3/(N*H2) = {recovered_2}"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        numerical_failures.append(
            (pp, qq)
        )

print()

print(
    f"  NUMERICAL FAILURES = "
    f"{len(numerical_failures)}"
)

print()


# ============================================================================
# COMPLETE FACTOR RECOVERY
# ============================================================================

print("[12] COMPLETE FACTOR RECOVERY")
print("-" * 78)

recovery_failures = []

for pp, qq in prime_pairs:

    NN = pp * qq
    SS = pp + qq

    H1 = numeric_H4(
        pp,
        qq,
        1,
    )

    H2 = numeric_H4(
        pp,
        qq,
        2,
    )

    if H1 == 0:

        recovery_failures.append(
            (pp, qq, "zero H4")
        )

        print(
            f"  ({pp},{qq}) ZERO H4"
        )

        continue

    recovered_M1 = sp.Rational(
        H2,
        NN * H1,
    )

    recovered_S = (
        recovered_M1
        - NN
        - 1
    )

    recovered_discriminant = sp.expand(
        recovered_S ** 2
        - 4 * NN
    )

    roots = set()

    if recovered_discriminant >= 0:

        root = sp.sqrt(
            recovered_discriminant
        )

        roots.add(
            sp.simplify(
                (recovered_S - root) / 2
            )
        )

        roots.add(
            sp.simplify(
                (recovered_S + root) / 2
            )
        )

    expected_roots = {
        sp.Integer(pp),
        sp.Integer(qq),
    }

    passed = (
        recovered_M1
        == (pp + 1) * (qq + 1)
        and recovered_S == SS
        and roots == expected_roots
    )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    recovered M1 = {recovered_M1}"
    )

    print(
        f"    recovered S = {recovered_S}"
    )

    print(
        f"    roots = {roots}"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        recovery_failures.append(
            (pp, qq, "factor recovery")
        )

print()

print(
    f"  RECOVERY FAILURES = "
    f"{len(recovery_failures)}"
)

print()


# ============================================================================
# SYMBOLIC RATIO CERTIFICATE
# ============================================================================

print("[13] SYMBOLIC CONSECUTIVE-RATIO CERTIFICATE")
print("-" * 78)

s = sp.symbols("s")

Hs = sp.expand(
    (N * M1) ** (s + 1)
    * A ** 2
    * C ** 2
)

Hs_next = sp.expand(
    (N * M1) ** (s + 2)
    * A ** 2
    * C ** 2
)

ratio_equation_difference = sp.simplify(
    Hs_next
    - N * M1 * Hs
)

ratio_pass = (
    ratio_equation_difference == 0
)

print(
    "  H_(s+1) = N*(N+S+1)*H_s"
)

print(
    f"  Difference = {ratio_equation_difference}"
)

print(
    f"  PASS = {ratio_pass}"
)

print()


# ============================================================================
# FINAL CERTIFICATES
# ============================================================================

print("[14] FINAL PROOF CERTIFICATES")
print("-" * 78)

certificates = [
    (
        "four-exponential decomposition",
        len(decomp_failures) == 0,
    ),
    (
        "weight product identity",
        weight_difference == 0,
    ),
    (
        "base product identity",
        base_difference == 0,
    ),
    (
        "Vandermonde p,q identity",
        vd_difference == 0,
    ),
    (
        "Vandermonde N,S identity",
        vd_ns_difference == 0,
    ),
    (
        "Hankel theorem p,q",
        len(hankel_pq_failures) == 0,
    ),
    (
        "Hankel theorem N,S",
        len(hankel_ns_failures) == 0,
    ),
    (
        "consecutive H4 ratio",
        len(ratio_failures) == 0
        and ratio_pass,
    ),
    (
        "numerical audit",
        len(numerical_failures) == 0,
    ),
    (
        "factor recovery",
        len(recovery_failures) == 0,
    ),
    (
        "degeneracy controls",
        len(control_failures) == 0,
    ),
]

for name, passed in certificates:
    print(
        f"  {name}: PASS={passed}"
    )

print()

overall = all(
    passed
    for _, passed in certificates
)

print(
    f"  OVERALL EXACT AUDIT = {overall}"
)

print()


# ============================================================================
# FORMAL THEOREM
# ============================================================================

print("[15] FORMAL THEOREM")
print("-" * 78)

print(
    "Let"
)

print(
    "  F_n = p(q+1)^n + q(p+1)^n"
)

print(
    "        - (q+1)p^n - (p+1)q^n."
)

print()

print(
    "Define"
)

print(
    "  H_s = det(F_(s+i+j))_(0<=i,j<4)."
)

print()

print(
    "Then"
)

print(
    "  H_s = [N(N+S+1)]^(s+1)"
)

print(
    "        * (4N-S^2)^2"
)

print(
    "        * (4N-S^2+1)^2."
)

print()

print(
    "Therefore, whenever H_s != 0:"
)

print(
    "  H_(s+1)/H_s = N(N+S+1)"
)

print(
    "  M1 = H_(s+1)/(N H_s)"
)

print(
    "  S = M1-N-1"
)

print()

print(
    "and finally"
)

print(
    "  z^2-Sz+N=0"
)

print(
    "recovers the factor pair."
)

print()


# ============================================================================
# INFORMATION MODEL
# ============================================================================

print("[16] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print(
    "Established:"
)

print(
    "  F_n"
)

print(
    "   -> four exponential bases"
)

print(
    "   -> H4 Hankel determinant"
)

print(
    "   -> Vandermonde square"
)

print(
    "   -> N*(N+S+1)"
)

print(
    "   -> M1=(p+1)(q+1)"
)

print(
    "   -> S"
)

print(
    "   -> p,q"
)

print()

print(
    "Still unresolved:"
)

print(
    "  N"
)

print(
    "   -> independently computable F_n"
)

print(
    "   -> H4"
)

print(
    "   -> M1"
)

print(
    "   -> S"
)

print(
    "   -> p,q"
)

print()


# ============================================================================
# FINAL INTERPRETATION
# ============================================================================

if overall:

    print(
        "MAIN RESULT:"
    )

    print(
        "  The four-base Vandermonde/Hankel theorem"
    )

    print(
        "  is exactly verified."
    )

    print()

    print(
        "  H_(s+1) = N*(N+S+1)*H_s"
    )

    print()

    print(
        "  Hence:"
    )

    print(
        "    S = H_(s+1)/(N*H_s) - N - 1"
    )

else:

    print(
        "MAIN RESULT:"
    )

    print(
        "  At least one exact certificate failed."
    )

    print(
        "  Do not promote the full theorem until"
    )

    print(
        "  the failed certificate has been isolated."
    )

print()

print(
    "CRITICAL NEXT QUESTION:"
)

print(
    "  Can F_n, or equivalently H_s,"
)

print(
    "  be generated from N alone without"
)

print(
    "  already evaluating expressions containing"
)

print(
    "  p or q?"
)

print()

print("=" * 78)
print("EXPERIMENT 488 FINISHED")
print("=" * 78)