#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp


print("EXPERIMENT 487 START")
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

def Fpq(ell: int):
    """
    Original k=1 KAPPA sequence:
        F_l = p(q+1)^l + q(p+1)^l
              - (q+1)p^l - (p+1)q^l
    """
    return sp.expand(
        p * (q + 1) ** ell
        + q * (p + 1) ** ell
        - (q + 1) * p ** ell
        - (p + 1) * q ** ell
    )


# ============================================================================
# EXACT SYMMETRIC REDUCTION p,q -> N,S
# ============================================================================

def to_NS(expr):
    """
    Rewrite a symmetric polynomial in p,q using
        N = pq
        S = p+q
    """
    expr = sp.expand(expr)

    result, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Unexpected symmetrization remainder: {remainder}"
        )

    for formal_var, concrete in mapping:

        if sp.expand(concrete - (p + q)) == 0:
            result = result.subs(formal_var, S)

        elif sp.expand(concrete - p * q) == 0:
            result = result.subs(formal_var, N)

        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: "
                f"{formal_var} -> {concrete}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs(
            {
                N: p * q,
                S: p + q,
            }
        )
        - expr
    )

    if check != 0:
        raise AssertionError(
            f"Symmetric reconstruction failed:\n{check}"
        )

    return result


# ============================================================================
# BUILD SYMBOLIC SEQUENCE
# ============================================================================

MAX_ELL = 20

F = {}

for ell in range(1, MAX_ELL + 1):
    F[ell] = to_NS(Fpq(ell))


# ============================================================================
# HANKEL MATRIX
# ============================================================================

def hankel(start: int, size: int = 4):
    """
    H[start] = (F[start+i+j])_(0<=i,j<size)
    """
    return sp.Matrix(
        [
            [
                F[start + i + j]
                for j in range(size)
            ]
            for i in range(size)
        ]
    )


def h4_symbolic(start: int):
    return sp.factor(
        hankel(start, 4).det()
    )


# ============================================================================
# [1] FOUR-EXPONENTIAL DECOMPOSITION
# ============================================================================

print("[1] ORIGINAL KAPPA EXPONENTIAL DECOMPOSITION")
print("-" * 78)

print(
    "F_l = p*(q+1)^l + q*(p+1)^l"
)
print(
    "      - (q+1)*p^l - (p+1)*q^l"
)
print()

bases = [
    p,
    q,
    p + 1,
    q + 1,
]

weights = [
    -(q + 1),
    -(p + 1),
    p,
    q,
]

print("Bases:")
for i, base in enumerate(bases, start=1):
    print(
        f"  lambda_{i} = {base}"
    )

print()
print("Weights:")
for i, weight in enumerate(weights, start=1):
    print(
        f"  c_{i} = {weight}"
    )

print()


# ============================================================================
# [2] VERIFY FOUR-EXPONENTIAL DECOMPOSITION
# ============================================================================

print("[2] FOUR-EXPONENTIAL DECOMPOSITION AUDIT")
print("-" * 78)

decomposition_failures = []

for ell in range(1, 11):

    reconstructed = sp.expand(
        sum(
            weights[i] * bases[i] ** ell
            for i in range(4)
        )
    )

    actual = Fpq(ell)

    passed = (
        sp.expand(reconstructed - actual) == 0
    )

    print(
        f"  ell={ell:2d}: PASS={passed}"
    )

    if not passed:
        decomposition_failures.append(ell)

print()
print(
    f"  DECOMPOSITION FAILURES = "
    f"{len(decomposition_failures)}"
)
print()


# ============================================================================
# [3] PRODUCT OF WEIGHTS
# ============================================================================

print("[3] PRODUCT OF EXPONENTIAL WEIGHTS")
print("-" * 78)

weight_product = sp.factor(
    sp.prod(weights)
)

expected_weight_product = sp.expand(
    N * (N + S + 1)
)

weight_product_difference = sp.factor(
    sp.expand(
        weight_product
        - expected_weight_product.subs(
            {
                N: p * q,
                S: p + q,
            }
        )
    )
)

print(
    f"  product(c_i) = {weight_product}"
)
print(
    f"  expected     = "
    f"{expected_weight_product.subs({N: p*q, S: p+q})}"
)
print(
    f"  difference   = {weight_product_difference}"
)
print(
    f"  PASS = {weight_product_difference == 0}"
)
print()


# ============================================================================
# [4] PRODUCT OF BASES
# ============================================================================

print("[4] PRODUCT OF EXPONENTIAL BASES")
print("-" * 78)

base_product = sp.factor(
    sp.prod(bases)
)

base_product_difference = sp.factor(
    sp.expand(
        base_product
        - expected_weight_product.subs(
            {
                N: p * q,
                S: p + q,
            }
        )
    )
)

print(
    f"  product(lambda_i) = {base_product}"
)
print(
    "  expected          = p*q*(p+1)*(q+1)"
)
print(
    f"  difference        = {base_product_difference}"
)
print(
    f"  PASS = {base_product_difference == 0}"
)
print()


# ============================================================================
# [5] VANDERMONDE PRODUCT
# ============================================================================

print("[5] VANDERMONDE STRUCTURE")
print("-" * 78)

pairwise_differences = []

for i in range(4):
    for j in range(i + 1, 4):

        diff = sp.expand(
            bases[j] - bases[i]
        )

        pairwise_differences.append(
            diff
        )

        print(
            f"  lambda_{j+1} - lambda_{i+1}"
            f" = {diff}"
        )

vandermonde = sp.factor(
    sp.prod(pairwise_differences)
)

print()
print(
    f"  Vandermonde = {vandermonde}"
)

vandermonde_square = sp.factor(
    sp.expand(
        vandermonde ** 2
    )
)

print(
    f"  Vandermonde^2 = {vandermonde_square}"
)
print()


# ============================================================================
# [6] VANDERMONDE SQUARE IN TERMS OF p,q
# ============================================================================

print("[6] VANDERMONDE SQUARE p,q IDENTITY")
print("-" * 78)

expected_vd_pq = sp.expand(
    (p - q) ** 4
    * (
        (p - q) ** 2 - 1
    ) ** 2
)

vd_pq_difference = sp.factor(
    sp.expand(
        vandermonde_square
        - expected_vd_pq
    )
)

print(
    "Expected:"
)
print(
    "  (p-q)^4 * ((p-q)^2-1)^2"
)
print(
    f"Difference = {vd_pq_difference}"
)
print(
    f"PASS = {vd_pq_difference == 0}"
)
print()


# ============================================================================
# [7] VANDERMONDE SQUARE IN N,S
# ============================================================================

print("[7] VANDERMONDE SQUARE N,S IDENTITY")
print("-" * 78)

Delta = 4 * N - S ** 2
Delta_plus = 4 * N - S ** 2 + 1

expected_vd_NS = sp.expand(
    Delta ** 2
    * Delta_plus ** 2
)

vd_NS_difference = sp.factor(
    sp.expand(
        vandermonde_square
        - expected_vd_NS.subs(
            {
                N: p * q,
                S: p + q,
            }
        )
    )
)

print(
    "Expected:"
)
print(
    "  (4N-S^2)^2 * (4N-S^2+1)^2"
)
print(
    f"Difference = {vd_NS_difference}"
)
print(
    f"PASS = {vd_NS_difference == 0}"
)
print()


# ============================================================================
# [8] GENERAL HANKEL THEOREM
# ============================================================================

print("[8] GENERAL 4x4 HANKEL THEOREM")
print("-" * 78)

theorem_failures = []

for start in range(1, 7):

    # General four-exponential identity:
    #
    # F_n = sum_i c_i lambda_i^n
    #
    # H_s = det(F_{s+i+j})_{0<=i,j<4}
    #
    #      = product(c_i lambda_i^s)
    #        * Vandermonde(lambda)^2
    #
    # Here:
    #
    # product(c_i)       = N*M1
    # product(lambda_i)  = N*M1
    #
    # therefore:
    #
    # H_s = [N*M1]^(s+1) * Vandermonde^2.

    M1 = N + S + 1

    predicted = sp.expand(
        (
            N * M1
        ) ** (start + 1)
        * Delta ** 2
        * Delta_plus ** 2
    )

    actual = h4_symbolic(start)

    predicted_pq = predicted.subs(
        {
            N: p * q,
            S: p + q,
        }
    )

    difference = sp.factor(
        sp.expand(
            actual
            - predicted_pq
        )
    )

    passed = (
        difference == 0
    )

    print(
        f"  start={start}: PASS={passed}"
    )

    if not passed:
        print(
            f"    difference = {difference}"
        )
        theorem_failures.append(
            start
        )

print()

print(
    f"  THEOREM FAILURES = "
    f"{len(theorem_failures)}"
)

print()


# ============================================================================
# [9] DIRECT RATIO THEOREM
# ============================================================================

print("[9] CONSECUTIVE HANKEL RATIO")
print("-" * 78)

M1 = N + S + 1

generic_hs = (
    N * M1
) ** 2 * Delta ** 2 * Delta_plus ** 2

generic_hs1 = (
    N * M1
) ** 3 * Delta ** 2 * Delta_plus ** 2

generic_ratio = sp.factor(
    sp.cancel(
        generic_hs1
        /
        generic_hs
    )
)

expected_ratio = sp.expand(
    N * M1
)

ratio_difference = sp.factor(
    sp.expand(
        generic_ratio
        - expected_ratio
    )
)

print(
    f"  H_(s+1)/H_s = {generic_ratio}"
)

print(
    f"  expected     = {expected_ratio}"
)

print(
    f"  difference   = {ratio_difference}"
)

print(
    f"  PASS = {ratio_difference == 0}"
)

print()


# ============================================================================
# [10] DIRECT M1 AND S RECOVERY
# ============================================================================

print("[10] SYMBOLIC M1 AND S RECOVERY")
print("-" * 78)

recovered_M1 = sp.cancel(
    generic_ratio / N
)

recovered_S = sp.expand(
    recovered_M1 - N - 1
)

M1_difference = sp.factor(
    sp.expand(
        recovered_M1 - M1
    )
)

S_difference = sp.factor(
    sp.expand(
        recovered_S - S
    )
)

print(
    f"  recovered M1 = {recovered_M1}"
)
print(
    f"  target M1    = {M1}"
)
print(
    f"  M1 difference = {M1_difference}"
)
print(
    f"  M1 PASS = {M1_difference == 0}"
)

print()

print(
    f"  recovered S = {recovered_S}"
)
print(
    f"  target S    = {S}"
)
print(
    f"  S difference = {S_difference}"
)
print(
    f"  S PASS = {S_difference == 0}"
)

print()


# ============================================================================
# [11] DEGENERACY ANALYSIS
# ============================================================================

print("[11] DEGENERACY ANALYSIS")
print("-" * 78)

print(
    "The general theorem contains the factor:"
)
print(
    "  (4N-S^2)^2 * (4N-S^2+1)^2"
)
print()

print(
    "For an integer factor pair:"
)
print(
    "  S^2-4N = (p-q)^2"
)
print()

print(
    "Therefore:"
)
print(
    "  4N-S^2 = -(p-q)^2"
)
print(
    "  4N-S^2+1 = 1-(p-q)^2"
)
print()

print(
    "Thus H_s can vanish only when:"
)
print(
    "  p=q"
)
print(
    "or"
)
print(
    "  |p-q|=1."
)
print()

print(
    "For distinct odd primes:"
)
print(
    "  |p-q| >= 2"
)
print(
    "so the Hankel bridge is nonzero."
)

print()


# ============================================================================
# [12] CONSECUTIVE-FACTOR CONTROLS
# ============================================================================

print("[12] CONSECUTIVE-FACTOR CONTROL")
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

    vals = {
        ell: (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )
        for ell in range(1, 9)
    }

    matrix = sp.Matrix(
        [
            [
                vals[1 + i + j]
                for j in range(4)
            ]
            for i in range(4)
        ]
    )

    h = int(
        matrix.det()
    )

    expected_zero = (
        abs(pp - qq) == 1
    )

    actual_zero = (
        h == 0
    )

    passed = (
        expected_zero == actual_zero
    )

    print(
        f"  ({pp},{qq})"
        f"  |p-q|={abs(pp-qq)}"
        f"  H4={h}"
        f"  expected_zero={expected_zero}"
        f"  PASS={passed}"
    )

    if not passed:
        control_failures.append(
            (
                pp,
                qq,
            )
        )

print()

print(
    f"  CONTROL FAILURES = "
    f"{len(control_failures)}"
)

print()


# ============================================================================
# [13] NUMERICAL PRIME-PAIR AUDIT
# ============================================================================

print("[13] NUMERICAL PRIME-PAIR HANKEL RATIO AUDIT")
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

numerical_failures = []

for pp, qq in prime_pairs:

    NN = pp * qq
    SS = pp + qq
    true_M1 = (
        pp + 1
    ) * (
        qq + 1
    )

    # IMPORTANT:
    # We need values through index 2+3+3 = 8.
    #
    # Build a complete range once, rather than the old
    # erroneous range(1, 8), which stopped at 7.

    vals = {
        ell: (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )
        for ell in range(1, 10)
    }

    def numeric_h4(start: int):

        required_max = (
            start
            + 3
            + 3
        )

        if required_max not in vals:
            raise KeyError(
                f"Missing sequence value F_{required_max}"
            )

        matrix = sp.Matrix(
            [
                [
                    vals[start + i + j]
                    for j in range(4)
                ]
                for i in range(4)
            ]
        )

        return int(
            matrix.det()
        )

    h1 = numeric_h4(1)
    h2 = numeric_h4(2)
    h3 = numeric_h4(3)

    if h1 == 0 or h2 == 0:
        print(
            f"  ({pp},{qq}) DEGENERATE"
        )
        numerical_failures.append(
            (
                pp,
                qq,
                "zero H4",
            )
        )
        continue

    recovered_1 = sp.Rational(
        h2,
        NN * h1,
    )

    recovered_2 = sp.Rational(
        h3,
        NN * h2,
    )

    ok = (
        recovered_1 == true_M1
        and recovered_2 == true_M1
    )

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    true M1      = {true_M1}"
    )
    print(
        f"    H4[2]/(N*H4[1]) = {recovered_1}"
    )
    print(
        f"    H4[3]/(N*H4[2]) = {recovered_2}"
    )
    print(
        f"    PASS = {ok}"
    )

    if not ok:
        numerical_failures.append(
            (
                pp,
                qq,
                "M1 ratio mismatch",
            )
        )

print()

print(
    f"  NUMERICAL FAILURES = "
    f"{len(numerical_failures)}"
)

print()


# ============================================================================
# [14] COMPLETE NUMERICAL RECOVERY CHAIN
# ============================================================================

print("[14] COMPLETE NUMERICAL RECOVERY CHAIN")
print("-" * 78)

chain_failures = []

for pp, qq in prime_pairs:

    NN = pp * qq
    SS = pp + qq

    true_M1 = (
        pp + 1
    ) * (
        qq + 1
    )

    vals = {
        ell: (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )
        for ell in range(1, 10)
    }

    h1 = int(
        sp.Matrix(
            [
                [
                    vals[1 + i + j]
                    for j in range(4)
                ]
                for i in range(4)
            ]
        ).det()
    )

    h2 = int(
        sp.Matrix(
            [
                [
                    vals[2 + i + j]
                    for j in range(4)
                ]
                for i in range(4)
            ]
        ).det()
    )

    if h1 == 0:
        print(
            f"  ({pp},{qq}) zero H4 -- skipped"
        )
        chain_failures.append(
            (
                pp,
                qq,
                "zero H4",
            )
        )
        continue

    recovered_M1 = sp.Rational(
        h2,
        NN * h1,
    )

    recovered_S = (
        recovered_M1
        - NN
        - 1
    )

    discriminant = sp.expand(
        recovered_S ** 2
        - 4 * NN
    )

    root = sp.sqrt(
        discriminant
    )

    roots = {
        sp.simplify(
            (recovered_S - root) / 2
        ),
        sp.simplify(
            (recovered_S + root) / 2
        ),
    }

    expected_roots = {
        sp.Integer(pp),
        sp.Integer(qq),
    }

    roots_ok = (
        roots == expected_roots
    )

    ok = (
        recovered_M1 == true_M1
        and recovered_S == SS
        and roots_ok
    )

    print(
        f"  ({pp},{qq})"
    )
    print(
        f"    recovered M1 = {recovered_M1}"
    )
    print(
        f"    recovered S  = {recovered_S}"
    )
    print(
        f"    recovered roots = {roots}"
    )
    print(
        f"    PASS = {ok}"
    )

    if not ok:
        chain_failures.append(
            (
                pp,
                qq,
                "recovery mismatch",
            )
        )

print()

print(
    f"  COMPLETE CHAIN FAILURES = "
    f"{len(chain_failures)}"
)

print()


# ============================================================================
# [15] DIRECT SYMBOLIC VANDERMONDE PROOF CERTIFICATE
# ============================================================================

print("[15] SYMBOLIC VANDERMONDE PROOF CERTIFICATE")
print("-" * 78)

certificate_checks = []

# Weight product:
certificate_checks.append(
    weight_product_difference == 0
)

# Vandermonde:
certificate_checks.append(
    vd_pq_difference == 0
)

# N,S Vandermonde:
certificate_checks.append(
    vd_NS_difference == 0
)

# General H4:
certificate_checks.append(
    len(theorem_failures) == 0
)

# Ratio:
certificate_checks.append(
    ratio_difference == 0
)

# Numerical:
certificate_checks.append(
    len(numerical_failures) == 0
)

# Recovery:
certificate_checks.append(
    len(chain_failures) == 0
)

# Controls:
certificate_checks.append(
    len(control_failures) == 0
)

certificate_names = [
    "weight product identity",
    "Vandermonde p,q identity",
    "Vandermonde N,S identity",
    "general H4 theorem",
    "consecutive ratio theorem",
    "numerical H4 audit",
    "complete recovery chain",
    "degeneracy controls",
]

for name, passed in zip(
    certificate_names,
    certificate_checks,
):
    print(
        f"  {name}: PASS={passed}"
    )

print()


# ============================================================================
# [16] FORMAL THEOREM
# ============================================================================

print("[16] FORMAL THEOREM")
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
    "Write"
)

print(
    "  H_s = det(F_{s+i+j})_(0<=i,j<4)."
)

print()

print(
    "Then F_n has four exponential bases"
)

print(
    "  lambda = {p,q,p+1,q+1},"
)

print(
    "with weights"
)

print(
    "  c = {-(q+1),-(p+1),p,q}."
)

print()

print(
    "Hence"
)

print(
    "  H_s"
)

print(
    "  = [product(c_i)]"
)

print(
    "    [product(lambda_i)]^s"
)

print(
    "    [Vandermonde(lambda)]^2."
)

print()

print(
    "For this kernel:"
)

print(
    "  product(c_i)"
)

print(
    "  = product(lambda_i)"
)

print(
    "  = pq(p+1)(q+1)"
)

print(
    "  = N(N+S+1)."
)

print()

print(
    "Furthermore:"
)

print(
    "  Vandermonde(lambda)^2"
)

print(
    "  = (p-q)^4 * ((p-q)^2-1)^2"
)

print(
    "  = (4N-S^2)^2"
)

print(
    "    * (4N-S^2+1)^2."
)

print()

print(
    "Therefore:"
)

print(
    "  H_s"
)

print(
    "  = [N(N+S+1)]^(s+1)"
)

print(
    "    * (4N-S^2)^2"
)

print(
    "    * (4N-S^2+1)^2."
)

print()

print(
    "Whenever H_s != 0:"
)

print(
    "  M1 = H_(s+1)/(N*H_s)"
)

print(
    "  S  = M1 - N - 1."
)

print()

print(
    "Finally:"
)

print(
    "  z^2 - S*z + N = 0"
)

print(
    "recovers p and q."
)

print()


# ============================================================================
# [17] INFORMATION-MODEL BOUNDARY
# ============================================================================

print("[17] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print(
    "The exact algebraic chain is:"
)

print(
    "  F-sequence"
)

print(
    "      -> four exponential bases"
)

print(
    "      -> 4x4 Hankel determinant"
)

print(
    "      -> Vandermonde square"
)

print(
    "      -> N*(N+S+1)"
)

print(
    "      -> M1=(p+1)(q+1)"
)

print(
    "      -> S"
)

print(
    "      -> p,q"
)

print()

print(
    "The unresolved source interface remains:"
)

print(
    "  N"
)

print(
    "   -> independently computable F_n"
)

print(
    "   -> H_s"
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
# [18] FINAL STATUS
# ============================================================================

overall = all(
    certificate_checks
)

print("[18] EXPERIMENT STATUS")
print("-" * 78)

for name, passed in zip(
    certificate_names,
    certificate_checks,
):
    print(
        f"  {name}: {passed}"
    )

print()

print(
    f"  OVERALL EXACT AUDIT = {overall}"
)

print()

if overall:
    print(
        "MAIN RESULT:"
    )

    print(
        "  The consecutive H4 bridge is now derived"
    )

    print(
        "  directly from the four-exponential/Vandermonde"
    )

    print(
        "  structure."
    )

    print()

    print(
        "  H_(s+1) = N*(N+S+1)*H_s"
    )

    print()

    print(
        "  Therefore:"
    )

    print(
        "    M1 = H_(s+1)/(N*H_s)"
    )

    print(
        "    S  = H_(s+1)/(N*H_s) - N - 1"
    )

    print()

    print(
        "  and finally:"
    )

    print(
        "    z^2 - S*z + N = 0"
    )

    print(
        "  recovers the factor pair."
    )

else:
    print(
        "MAIN RESULT:"
    )

    print(
        "  At least one certificate failed."
    )

    print(
        "  Inspect the failed section before"
    )

    print(
        "  promoting the theorem to the research record."
    )

print()

print(
    "CRITICAL NEXT QUESTION:"
)

print(
    "  Can the F_n sequence itself be generated"
)

print(
    "  from N alone without evaluating expressions"
)

print(
    "  that already encode p and q?"
)

print()

print("=" * 78)
print("EXPERIMENT 487 FINISHED")
print("=" * 78)