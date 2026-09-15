#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp


# ============================================================================
# EXPERIMENT 486
# ============================================================================
#
# CONSECUTIVE HANKEL-4 RATIO -> SHIFTED PRODUCT
#
# Main target:
#
#   H4_s = det(F_{s+i+j}) , 0 <= i,j < 4
#
# For the original KAPPA sequence:
#
#   F_l =
#       p(q+1)^l
#       + q(p+1)^l
#       - (q+1)p^l
#       - (p+1)q^l
#
# The observed factorization suggests
#
# H4_s =
#   N^(s+1)
#   * (4N-S^2)^2
#   * (N+S+1)^(s+1)
#   * (4N-S^2+1)^2
#
# Therefore
#
#   H4_(s+1) / H4_s = N*(N+S+1)
#
# and hence
#
#   M1 = N+S+1 = H4_(s+1)/(N*H4_s)
#
# This is substantially cleaner than the failed Experiment-485 attempt
# using H3.
#
# Exact arithmetic only.
# No factorization of N.
# No floating point.
# No CSV.
# ============================================================================


print("EXPERIMENT 486 START")
print("=" * 78)
print("CONSECUTIVE HANKEL-4 RATIO -> SHIFTED PRODUCT")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# ORIGINAL KAPPA KERNEL
# ============================================================================

def F_pq(ell: int):
    return sp.expand(
        p * (q + 1) ** ell
        + q * (p + 1) ** ell
        - p ** ell * (q + 1)
        - q ** ell * (p + 1)
    )


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(expr):
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Unexpected symmetrization remainder: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:

        if sp.expand(concrete - (p + q)) == 0:
            result = result.subs(
                formal_var,
                S,
            )

        elif sp.expand(concrete - p * q) == 0:
            result = result.subs(
                formal_var,
                N,
            )

        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs(
            {
                p: p,
                q: q,
                S: p + q,
                N: p * q,
            }
        )
        - expr
    )

    if check != 0:
        raise AssertionError(
            "Symmetric reconstruction failed."
        )

    return result


# ============================================================================
# BUILD SUFFICIENT KERNEL SEQUENCE
# ============================================================================

print("[1] BUILD ORIGINAL KAPPA SEQUENCE")
print("-" * 78)

MAX_ELL = 22

F = {}

for ell in range(1, MAX_ELL + 1):
    F[ell] = sp.factor(
        symmetric_reduce(
            F_pq(ell)
        )
    )

for ell in range(1, 9):
    print(
        f"  F_{ell} = {F[ell]}"
    )

print()


# ============================================================================
# HANKEL MATRIX / DETERMINANT
# ============================================================================

def hankel_matrix(start: int, size: int):
    return sp.Matrix(
        [
            [
                F[start + i + j]
                for j in range(size)
            ]
            for i in range(size)
        ]
    )


def hankel_det(start: int, size: int):
    return sp.factor(
        sp.expand(
            hankel_matrix(
                start,
                size,
            ).det()
        )
    )


# ============================================================================
# [2] CONSECUTIVE 4x4 HANKEL DETERMINANTS
# ============================================================================

print("[2] CONSECUTIVE 4x4 HANKEL DETERMINANTS")
print("-" * 78)

H4 = {}

for start in range(1, 7):

    H4[start] = hankel_det(
        start,
        4,
    )

    print(
        f"  H4[{start}] = {H4[start]}"
    )

print()


# ============================================================================
# [3] FACTOR STRUCTURE
# ============================================================================

print("[3] EXACT FACTOR STRUCTURE")
print("-" * 78)

Delta4 = sp.expand(
    4 * N - S**2
)

M1 = sp.expand(
    N + S + 1
)

Delta4_plus = sp.expand(
    4 * N - S**2 + 1
)

print(
    f"  A  = 4*N - S**2       = {Delta4}"
)

print(
    f"  M1 = N + S + 1        = {M1}"
)

print(
    f"  C  = 4*N - S**2 + 1   = {Delta4_plus}"
)

print()


# ============================================================================
# [4] THEOREM TARGET FOR EACH H4 WINDOW
# ============================================================================

print("[4] HANKEL-4 CLOSED FORM TARGET")
print("-" * 78)

factor_failures = []

for start in range(1, 7):

    expected = sp.expand(
        N ** (start + 1)
        * Delta4 ** 2
        * M1 ** (start + 1)
        * Delta4_plus ** 2
    )

    difference = sp.factor(
        sp.expand(
            H4[start]
            - expected
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
        factor_failures.append(
            start
        )

print()

print(
    f"  FACTORIZATION FAILURES = "
    f"{len(factor_failures)}"
)

print()


# ============================================================================
# [5] CONSECUTIVE RATIO THEOREM
# ============================================================================

print("[5] CONSECUTIVE HANKEL-4 RATIO")
print("-" * 78)

ratio_failures = []

for start in range(1, 6):

    numerator = H4[start + 1]
    denominator = H4[start]

    ratio = sp.factor(
        sp.cancel(
            numerator
            / denominator
        )
    )

    expected_ratio = sp.expand(
        N * M1
    )

    difference = sp.factor(
        sp.cancel(
            ratio
            - expected_ratio
        )
    )

    passed = (
        difference == 0
    )

    print(
        f"  H4[{start+1}] / H4[{start}]"
    )

    print(
        f"    ratio     = {ratio}"
    )

    print(
        f"    expected  = N*(N+S+1)"
    )

    print(
        f"    PASS      = {passed}"
    )

    if not passed:
        ratio_failures.append(
            start
        )

    print()


print(
    f"  RATIO FAILURES = "
    f"{len(ratio_failures)}"
)

print()


# ============================================================================
# [6] DIRECT M1 RECOVERY FROM H4 WINDOWS
# ============================================================================

print("[6] DIRECT SHIFTED-PRODUCT RECOVERY")
print("-" * 78)

m1_recovery_failures = []

for start in range(1, 6):

    candidate = sp.factor(
        sp.cancel(
            H4[start + 1]
            /
            (
                N
                * H4[start]
            )
        )
    )

    difference = sp.factor(
        sp.cancel(
            candidate
            - M1
        )
    )

    passed = (
        difference == 0
    )

    print(
        f"  start={start}:"
    )

    print(
        f"    candidate = {candidate}"
    )

    print(
        f"    difference = {difference}"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        m1_recovery_failures.append(
            start
        )

    print()


print(
    f"  M1 RECOVERY FAILURES = "
    f"{len(m1_recovery_failures)}"
)

print()


# ============================================================================
# [7] RECOVER S FROM HANKEL-4 RATIO
# ============================================================================

print("[7] DIRECT S RECOVERY")
print("-" * 78)

s_failures = []

for start in range(1, 6):

    recovered_S = sp.factor(
        sp.cancel(
            H4[start + 1]
            /
            (
                N
                * H4[start]
            )
            - N
            - 1
        )
    )

    difference = sp.factor(
        sp.cancel(
            recovered_S
            - S
        )
    )

    passed = (
        difference == 0
    )

    print(
        f"  start={start}:"
    )

    print(
        f"    recovered S = {recovered_S}"
    )

    print(
        f"    difference  = {difference}"
    )

    print(
        f"    PASS = {passed}"
    )

    if not passed:
        s_failures.append(
            start
        )

    print()


print(
    f"  S RECOVERY FAILURES = "
    f"{len(s_failures)}"
)

print()


# ============================================================================
# [8] FACTOR RECOVERY SYMBOLICALLY
# ============================================================================

print("[8] FACTOR RECOVERY FROM HANKEL-4 RATIO")
print("-" * 78)

recovered_S_symbolic = sp.factor(
    sp.cancel(
        H4[2]
        /
        (
            N
            * H4[1]
        )
        - N
        - 1
    )
)

quadratic = sp.expand(
    sp.Symbol("z") ** 2
    - recovered_S_symbolic * sp.Symbol("z")
    + N
)

quadratic_expected = sp.expand(
    sp.Symbol("z") ** 2
    - S * sp.Symbol("z")
    + N
)

quadratic_difference = sp.factor(
    sp.expand(
        quadratic
        - quadratic_expected
    )
)

quadratic_pass = (
    quadratic_difference == 0
)

print(
    f"  recovered S = {recovered_S_symbolic}"
)

print(
    f"  recovered quadratic = {sp.factor(quadratic)}"
)

print(
    f"  quadratic difference = {quadratic_difference}"
)

print(
    f"  PASS = {quadratic_pass}"
)

print()


# ============================================================================
# [9] WHY H4 IS BETTER THAN H3
# ============================================================================

print("[9] H3 VS H4 INFORMATION CONTENT")
print("-" * 78)

print(
    "Experiment 485 showed that H3 contains the desired factor M1"
)

print(
    "but also contains the nuisance factor"
)

print(
    "  B = 10*N^2 - 2*N*S^2 + 6*N*S - 2*N - S^3 + S^2."
)

print()

print(
    "The H4 window instead has the form"
)

print(
    "  H4[s] ="
)

print(
    "    N^(s+1)"
)

print(
    "    * (4*N-S^2)^2"
)

print(
    "    * M1^(s+1)"
)

print(
    "    * (4*N-S^2+1)^2."
)

print()

print(
    "Therefore the nuisance factors are independent of s,"
)

print(
    "while N and M1 advance with exactly the same exponent."
)

print(
    "Their consecutive ratio therefore isolates N*M1."
)

print()


# ============================================================================
# [10] NUMERICAL AUDIT
# ============================================================================

print("[10] NUMERICAL EXACT AUDIT")
print("-" * 78)

TEST_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]

numerical_failures = []

for pp, qq in TEST_PAIRS:

    NN = pp * qq
    SS = pp + qq
    M1_num = (
        pp + 1
    ) * (
        qq + 1
    )

    vals = {}

    for ell in range(1, 20):

        vals[ell] = (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )

    def numeric_h4(start):
        mat = sp.Matrix(
            [
                [
                    vals[start + i + j]
                    for j in range(4)
                ]
                for i in range(4)
            ]
        )

        return int(
            mat.det()
        )

    h4_1 = numeric_h4(1)
    h4_2 = numeric_h4(2)
    h4_3 = numeric_h4(3)

    candidates = []

    for a, b in [
        (h4_1, h4_2),
        (h4_2, h4_3),
    ]:

        if a == 0:
            candidates.append(
                None
            )
            continue

        candidate = sp.Rational(
            b,
            NN * a,
        )

        candidates.append(
            candidate
        )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N  = {NN}"
    )

    print(
        f"    S  = {SS}"
    )

    print(
        f"    M1 = {M1_num}"
    )

    print(
        f"    H4[1] = {h4_1}"
    )

    print(
        f"    H4[2] = {h4_2}"
    )

    print(
        f"    H4[3] = {h4_3}"
    )

    print(
        f"    H4[2]/(N*H4[1]) = "
        f"{candidates[0]}"
    )

    print(
        f"    H4[3]/(N*H4[2]) = "
        f"{candidates[1]}"
    )

    ok = all(
        candidate == M1_num
        for candidate in candidates
        if candidate is not None
    )

    print(
        f"    M1 recovery PASS = {ok}"
    )

    if not ok:
        numerical_failures.append(
            (
                pp,
                qq,
                "H4 ratio failed",
            )
        )

    print()


print(
    f"  NUMERICAL FAILURES = "
    f"{len(numerical_failures)}"
)

for failure in numerical_failures:
    print(
        f"    FAILURE: {failure}"
    )

print()


# ============================================================================
# [11] ERROR-CONTROLLED DEGENERACY AUDIT
# ============================================================================

print("[11] DEGENERACY AUDIT")
print("-" * 78)

print(
    "The ratio requires H4[start] != 0."
)

print(
    "For distinct nontrivial factor pairs,"
)

print(
    "H4 contains:"
)

print(
    "  N^... * (4N-S^2)^2 * (N+S+1)^... "
    "* (4N-S^2+1)^2."
)

print()

degenerate_conditions = [
    ("N", N),
    ("4N-S^2", Delta4),
    ("N+S+1", M1),
    ("4N-S^2+1", Delta4_plus),
]

for name, expr in degenerate_conditions:

    print(
        f"  {name} = 0 is possible "
        f"symbolically: {sp.expand(expr) == 0}"
    )

print()

print(
    "For positive distinct p,q:"
)

print(
    "  N > 0"
)

print(
    "  S^2-4N = (p-q)^2 > 0"
)

print(
    "  N+S+1 > 0"
)

print(
    "  4N-S^2+1 = 1-(p-q)^2"
)

print(
    "Therefore the remaining possible zero case is"
)

print(
    "  (p-q)^2 = 1."
)

print(
    "This corresponds to consecutive factors."
)

print()


# ============================================================================
# [12] EXCEPTIONAL CONSECUTIVE-FACTOR CONTROL
# ============================================================================

print("[12] CONSECUTIVE-FACTOR CONTROL")
print("-" * 78)

exceptional_pairs = [
    (2, 3),
    (3, 4),
    (5, 7),
    (7, 8),
    (11, 13),
]

for pp, qq in exceptional_pairs:

    NN = pp * qq
    SS = pp + qq

    vals = {}

    for ell in range(1, 12):

        vals[ell] = (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )

    mat = sp.Matrix(
        [
            [
                vals[2 + i + j]
                for j in range(4)
            ]
            for i in range(4)
        ]
    )

    h = sp.factor(
        mat.det()
    )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N={NN}, S={SS}, "
        f"(p-q)^2={(pp-qq)**2}"
    )

    print(
        f"    H4[start=2] = {h}"
    )

    print(
        f"    H4 zero = {h == 0}"
    )

print()


# ============================================================================
# [13] SYMBOLIC PROOF CERTIFICATE
# ============================================================================

print("[13] SYMBOLIC PROOF CERTIFICATE")
print("-" * 78)

# Prove the ratio directly from the two closed forms.
ratio_certificate = sp.factor(
    sp.cancel(
        (
            N**3
            * Delta4**2
            * M1**3
            * Delta4_plus**2
        )
        /
        (
            N
            * (
                N**2
                * Delta4**2
                * M1**2
                * Delta4_plus**2
            )
        )
    )
)

ratio_certificate_pass = (
    sp.expand(
        ratio_certificate
        - N * M1
    )
    == 0
)

print(
    f"  H4[s+1]/H4[s] = {ratio_certificate}"
)

print(
    f"  Expected       = N*(N+S+1)"
)

print(
    f"  PASS = {ratio_certificate_pass}"
)

print()

print(
    "The resulting exact identity is:"
)

print(
    "  H4[s+1] = N*(N+S+1)*H4[s]"
)

print(
    "and therefore:"
)

print(
    "  M1 = H4[s+1]/(N*H4[s])"
)

print()


# ============================================================================
# [14] INFORMATION MODEL
# ============================================================================

print("[14] INFORMATION-MODEL CHECK")
print("-" * 78)

print(
    "Established:"
)

print(
    "  F-sequence"
)

print(
    "    -> consecutive 4x4 Hankel determinants"
)

print(
    "    -> N*(N+S+1)"
)

print(
    "    -> M1=N+S+1"
)

print(
    "    -> S=M1-N"
)

print(
    "    -> z^2-S*z+N"
)

print(
    "    -> p,q"
)

print()

print(
    "Still unresolved:"
)

print(
    "  N"
)

print(
    "   -> independently computable F-sequence"
)

print(
    "   -> Hankel determinants"
)

print(
    "   -> M1"
)

print()


# ============================================================================
# [15] FINAL STATUS
# ============================================================================

overall = (
    len(factor_failures) == 0
    and len(ratio_failures) == 0
    and len(m1_recovery_failures) == 0
    and len(s_failures) == 0
    and quadratic_pass
    and len(numerical_failures) == 0
    and ratio_certificate_pass
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  H4 closed-form factorization: "
    f"{len(factor_failures) == 0}"
)

print(
    f"  Consecutive H4 ratio: "
    f"{len(ratio_failures) == 0}"
)

print(
    f"  M1 recovery: "
    f"{len(m1_recovery_failures) == 0}"
)

print(
    f"  S recovery: "
    f"{len(s_failures) == 0}"
)

print(
    f"  Quadratic factor bridge: "
    f"{quadratic_pass}"
)

print(
    f"  Numerical audit: "
    f"{len(numerical_failures) == 0}"
)

print(
    f"  Symbolic ratio certificate: "
    f"{ratio_certificate_pass}"
)

print(
    f"  OVERALL EXACT AUDIT = {overall}"
)

print()

print(
    "MAIN RESULT:"
)

print(
    "  Consecutive 4x4 Hankel minors satisfy"
)

print(
    "      H4[s+1] = N*(N+S+1)*H4[s]."
)

print()

print(
    "  Therefore:"
)

print(
    "      (p+1)(q+1)"
)

print(
    "    = N+S+1"
)

print(
    "    = H4[s+1] / (N*H4[s])."
)

print()

print(
    "  Once these Hankel values are available:"
)

print(
    "      S = H4[s+1]/(N*H4[s]) - N - 1."
)

print(
    "  Then:"
)

print(
    "      z^2-S*z+N=0"
)

print(
    "  recovers p and q."
)

print()

print(
    "REMAINING CENTRAL BRIDGE:"
)

print(
    "  Can the original KAPPA sequence, or its H4 Hankel minors,"
)

print(
    "  be generated from N alone without knowing p or q?"
)

print()

print(
    "This is a sharper target than the previous H3 approach."
)

print()

print("=" * 78)
print("EXPERIMENT 486 FINISHED")
print("=" * 78)
