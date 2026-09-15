#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp


# ============================================================================
# EXPERIMENT 485
# ============================================================================
#
# HANKEL-DETERMINANT BRIDGE TO THE SHIFTED PRODUCT
#
# This version avoids passing rational expressions to sp.Poly().
#
# Main targets:
#
#   M1    = (p+1)(q+1) = N+S+1
#   Delta = (p-q)^2     = S^2-4N
#
# The k=1 KAPPA sequence is:
#
#   F_l =
#       p(q+1)^l
#       + q(p+1)^l
#       - (q+1)p^l
#       - (p+1)q^l
#
# It is a four-exponential sequence with bases
#
#   p, q, p+1, q+1.
#
# Therefore the Hankel rank is at most four and every 5x5 Hankel
# determinant must vanish.
#
# Experiment 484 found:
#
#   H2 = -F2 * B
#
#   H3 = N*Delta*(N+S+1)*(Delta+1)*B
#
# with a common nuisance factor B.
#
# The present experiment audits that structure carefully and searches
# for cancellation patterns without assuming rational expressions are
# polynomials.
#
# Exact arithmetic only.
# No factorization of N.
# No floating point.
# No CSV.
# ============================================================================


print("EXPERIMENT 485 START")
print("=" * 78)
print("HANKEL-DETERMINANT BRIDGE TO THE SHIFTED PRODUCT")
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
        p * (1 + q) ** ell
        + q * (1 + p) ** ell
        - p ** ell * (1 + q)
        - q ** ell * (1 + p)
    )


def symmetric_reduce(expr):
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization remainder: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if sp.expand(concrete - (p + q)) == 0:
            result = result.subs(formal_var, S)

        elif sp.expand(concrete - p * q) == 0:
            result = result.subs(formal_var, N)

        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs(
            {
                S: p + q,
                N: p * q,
            }
        ) - expr
    )

    if check != 0:
        raise AssertionError(
            "Symmetric reconstruction failed."
        )

    return result


# ============================================================================
# BUILD SEQUENCE
# ============================================================================

print("[1] BUILD ORIGINAL KAPPA SEQUENCE")
print("-" * 78)

F = {
    ell: sp.factor(
        symmetric_reduce(
            F_pq(ell)
        )
    )
    for ell in range(1, 15)
}

for ell in range(1, 9):
    print(
        f"  F_{ell} = {F[ell]}"
    )

print()


# ============================================================================
# HANKEL HELPER
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
    M = hankel_matrix(
        start,
        size,
    )

    return sp.factor(
        sp.expand(
            M.det()
        )
    )


# ============================================================================
# [2] BASIC HANKEL DETERMINANTS
# ============================================================================

print("[2] BASIC HANKEL DETERMINANTS")
print("-" * 78)

H2 = hankel_det(2, 2)
H3 = hankel_det(2, 3)
H4 = hankel_det(2, 4)
H5 = hankel_det(1, 5)

print(
    f"  H2 = {H2}"
)

print(
    f"  H3 = {H3}"
)

print(
    f"  H4 = {H4}"
)

print(
    f"  H5 = {H5}"
)

print()


# ============================================================================
# [3] FOUR-BASE RANK CERTIFICATE
# ============================================================================

print("[3] FOUR-BASE HANKEL RANK CERTIFICATE")
print("-" * 78)

print(
    "Bases:"
)
print(
    "  p, q, p+1, q+1"
)

print(
    "Expected Hankel rank <= 4."
)

print(
    f"  H5 = {H5}"
)

H5_pass = sp.expand(H5) == 0

print(
    f"  H5 == 0 PASS = {H5_pass}"
)

print()


# ============================================================================
# [4] SHIFTED HANKEL WINDOWS
# ============================================================================

print("[4] SHIFTED HANKEL WINDOWS")
print("-" * 78)

H2_shifted = {}
H3_shifted = {}
H4_shifted = {}
H5_shifted = {}

for start in range(1, 8):

    H2_shifted[start] = hankel_det(
        start,
        2,
    )

    H3_shifted[start] = hankel_det(
        start,
        3,
    )

    H4_shifted[start] = hankel_det(
        start,
        4,
    )

    # Only compute 5x5 where enough F_l values exist.
    if start + 8 <= max(F):
        H5_shifted[start] = hankel_det(
            start,
            5,
        )

    print(
        f"  start={start}:"
    )

    print(
        f"    H2 = {H2_shifted[start]}"
    )

    print(
        f"    H3 = {H3_shifted[start]}"
    )

    print(
        f"    H4 = {H4_shifted[start]}"
    )

    if start in H5_shifted:
        print(
            f"    H5 = {H5_shifted[start]}"
        )

print()


# ============================================================================
# [5] KNOWN STRUCTURAL TARGETS
# ============================================================================

print("[5] KNOWN STRUCTURAL TARGETS")
print("-" * 78)

M1 = sp.expand(
    N + S + 1
)

Delta = sp.expand(
    S**2 - 4*N
)

Delta_plus_1 = sp.expand(
    Delta + 1
)

B484 = sp.expand(
    10*N**2
    - 2*N*S**2
    + 6*N*S
    - 2*N
    - S**3
    + S**2
)

F2 = F[2]

print(
    f"  M1       = {M1}"
)

print(
    f"  Delta    = {Delta}"
)

print(
    f"  Delta+1  = {Delta_plus_1}"
)

print(
    f"  F2       = {F2}"
)

print(
    f"  B484     = {B484}"
)

print()


# ============================================================================
# [6] H2 EXACT FACTORIZATION
# ============================================================================

print("[6] H2 EXACT FACTORIZATION")
print("-" * 78)

H2_expected = sp.expand(
    -F2 * B484
)

H2_difference = sp.factor(
    sp.expand(
        H2 - H2_expected
    )
)

H2_pass = (
    H2_difference == 0
)

print(
    f"  H2 + F2*B = {H2_difference}"
)

print(
    f"  H2 = -F2*B PASS = {H2_pass}"
)

print()


# ============================================================================
# [7] H3 EXACT FACTORIZATION
# ============================================================================

print("[7] H3 EXACT FACTORIZATION")
print("-" * 78)

H3_expected = sp.expand(
    N
    * Delta
    * M1
    * Delta_plus_1
    * B484
)

H3_difference = sp.factor(
    sp.expand(
        H3 - H3_expected
    )
)

H3_pass = (
    H3_difference == 0
)

print(
    f"  H3 difference = {H3_difference}"
)

print(
    f"  H3 factorization PASS = {H3_pass}"
)

print()


# ============================================================================
# [8] H3/H2 NUISANCE CANCELLATION
# ============================================================================

print("[8] H3/H2 NUISANCE CANCELLATION")
print("-" * 78)

ratio_H3_H2 = sp.factor(
    sp.cancel(
        H3 / H2
    )
)

expected_ratio = sp.factor(
    sp.cancel(
        -N
        * Delta
        * M1
        * Delta_plus_1
        / F2
    )
)

ratio_difference = sp.factor(
    sp.cancel(
        ratio_H3_H2
        - expected_ratio
    )
)

ratio_pass = (
    ratio_difference == 0
)

print(
    f"  H3/H2 = {ratio_H3_H2}"
)

print(
    f"  Expected = {expected_ratio}"
)

print(
    f"  Difference = {ratio_difference}"
)

print(
    f"  PASS = {ratio_pass}"
)

print()


# ============================================================================
# [9] M1 RECONSTRUCTION
# ============================================================================

print("[9] DIRECT M1 RECONSTRUCTION")
print("-" * 78)

candidate_M1 = sp.factor(
    sp.cancel(
        -H3
        * F2
        / (
            H2
            * N
            * Delta
            * Delta_plus_1
        )
    )
)

candidate_M1_difference = sp.factor(
    sp.cancel(
        candidate_M1
        - M1
    )
)

candidate_M1_pass = (
    candidate_M1_difference == 0
)

print(
    "Candidate:"
)

print(
    "  -H3*F2 / [H2*N*Delta*(Delta+1)]"
)

print(
    f"  = {candidate_M1}"
)

print()

print(
    f"  Difference from M1 = "
    f"{candidate_M1_difference}"
)

print(
    f"  PASS = {candidate_M1_pass}"
)

print()


# ============================================================================
# [10] H4 FACTORIZATION EXPLORATION
# ============================================================================

print("[10] H4 FACTORIZATION EXPLORATION")
print("-" * 78)

print(
    f"  H4 = {H4}"
)

print()

# Test obvious factors.
h4_targets = {
    "N": N,
    "Delta": Delta,
    "Delta+1": Delta_plus_1,
    "M1": M1,
    "B484": B484,
    "F2": F2,
}

for name, target in h4_targets.items():

    poly_H4 = sp.Poly(
        H4,
        S,
        domain=sp.QQ.frac_field(N),
    )

    poly_target = sp.Poly(
        target,
        S,
        domain=sp.QQ.frac_field(N),
    )

    _, rem = sp.div(
        poly_H4,
        poly_target,
    )

    print(
        f"  H4 divisible by {name}: "
        f"{rem.as_expr() == 0}"
    )

print()


# ============================================================================
# [11] SAFE RATIONAL COMPLEXITY AUDIT
# ============================================================================

print("[11] SAFE RATIONAL COMPLEXITY AUDIT")
print("-" * 78)

print(
    "Rational expressions are split into numerator and denominator"
)
print(
    "before polynomial degree inspection."
)
print()


def rational_degree_info(expr):
    expr = sp.cancel(expr)

    numerator, denominator = sp.fraction(
        expr
    )

    numerator = sp.expand(
        numerator
    )

    denominator = sp.expand(
        denominator
    )

    numerator_poly = sp.Poly(
        numerator,
        N,
        S,
    )

    denominator_poly = sp.Poly(
        denominator,
        N,
        S,
    )

    return {
        "numerator": numerator,
        "denominator": denominator,
        "num_degree": numerator_poly.total_degree(),
        "den_degree": denominator_poly.total_degree(),
    }


# Selected shifted-window ratios.
for start in range(1, 7):

    if H2_shifted[start] == 0:
        continue

    for other_start in range(start + 1, 7):

        if other_start not in H3_shifted:
            continue

        expr = sp.cancel(
            H3_shifted[other_start]
            / H2_shifted[start]
            / M1
        )

        info = rational_degree_info(
            expr
        )

        print(
            f"  H3(start={other_start}) / "
            f"H2(start={start}) / M1:"
        )

        print(
            f"    numerator degree   = "
            f"{info['num_degree']}"
        )

        print(
            f"    denominator degree = "
            f"{info['den_degree']}"
        )

print()


# ============================================================================
# [12] SEARCH FOR CLEAN SHIFTED-WINDOW CANCELLATIONS
# ============================================================================

print("[12] CLEAN SHIFTED-WINDOW CANCELLATION SEARCH")
print("-" * 78)

clean_hits = []

for a in range(1, 7):

    for b in range(a + 1, 7):

        if H2_shifted[a] == 0:
            continue

        expr = sp.cancel(
            H3_shifted[b]
            / H2_shifted[a]
        )

        # Check whether M1 divides the numerator or denominator
        # after rational normalization.
        num, den = sp.fraction(expr)

        num_poly = sp.Poly(
            sp.expand(num),
            S,
            domain=sp.QQ.frac_field(N),
        )

        M1_poly = sp.Poly(
            M1,
            S,
            domain=sp.QQ.frac_field(N),
        )

        _, rem_num = sp.div(
            num_poly,
            M1_poly,
        )

        if rem_num.as_expr() == 0:

            clean_hits.append(
                (
                    a,
                    b,
                    "numerator",
                )
            )

        den_poly = sp.Poly(
            sp.expand(den),
            S,
            domain=sp.QQ.frac_field(N),
        )

        _, rem_den = sp.div(
            den_poly,
            M1_poly,
        )

        if rem_den.as_expr() == 0:

            clean_hits.append(
                (
                    a,
                    b,
                    "denominator",
                )
            )

if clean_hits:

    for hit in clean_hits:
        print(
            f"  H3(start={hit[1]}) / "
            f"H2(start={hit[0]}) has M1 in "
            f"{hit[2]}"
        )

else:

    print(
        "  No additional shifted-window M1 cancellation found."
    )

print()


# ============================================================================
# [13] FOUR-BASE ELEMENTARY SYMMETRIC DATA
# ============================================================================

print("[13] FOUR-BASE ELEMENTARY SYMMETRIC DATA")
print("-" * 78)

bases = [
    p,
    q,
    p + 1,
    q + 1,
]

e1 = sp.expand(
    sum(bases)
)

e2 = sp.expand(
    sum(
        bases[i] * bases[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
)

e3 = sp.expand(
    sum(
        bases[i] * bases[j] * bases[k]
        for i in range(4)
        for j in range(i + 1, 4)
        for k in range(j + 1, 4)
    )
)

e4 = sp.expand(
    sp.prod(
        bases
    )
)

print(
    f"  e1 = {symmetric_reduce(e1)}"
)

print(
    f"  e2 = {symmetric_reduce(e2)}"
)

print(
    f"  e3 = {symmetric_reduce(e3)}"
)

print(
    f"  e4 = {symmetric_reduce(e4)}"
)

print()


# ============================================================================
# [14] CHARACTERISTIC POLYNOMIAL
# ============================================================================

print("[14] CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

z = sp.symbols("z")

chi_factorized = sp.expand(
    (z - p)
    * (z - q)
    * (z - p - 1)
    * (z - q - 1)
)

chi_NS = symmetric_reduce(
    chi_factorized
)

print(
    f"  chi(z) = {sp.factor(chi_NS)}"
)

expected_chi = sp.expand(
    z**4
    - 2*(S + 1)*z**3
    + (2*N + S**2 + 3*S + 1)*z**2
    - (2*N + S)*(S + 1)*z
    + N*(N + S + 1)
)

chi_pass = (
    sp.expand(
        chi_NS - expected_chi
    )
    == 0
)

print(
    f"  characteristic identity PASS = "
    f"{chi_pass}"
)

print()


# ============================================================================
# [15] SECOND STRATEGY: SHIFTED PRODUCT AS A ROOT PRODUCT
# ============================================================================

print("[15] SHIFTED PRODUCT AS ROOT PRODUCT")
print("-" * 78)

print(
    "The two shifted roots p+1 and q+1 have:"
)

print(
    "  sum     = S+2"
)

print(
    "  product = N+S+1 = M1"
)

print()

root_sum = sp.expand(
    (p + 1) + (q + 1)
)

root_product = sp.expand(
    (p + 1) * (q + 1)
)

root_sum_NS = symmetric_reduce(
    root_sum
)

root_product_NS = symmetric_reduce(
    root_product
)

print(
    f"  root sum     = {root_sum_NS}"
)

print(
    f"  root product = {root_product_NS}"
)

root_product_pass = (
    sp.expand(
        root_product_NS
        - M1
    )
    == 0
)

print(
    f"  M1 identity PASS = "
    f"{root_product_pass}"
)

print()


# ============================================================================
# [16] NUMERICAL EXACT AUDIT
# ============================================================================

print("[16] NUMERICAL EXACT AUDIT")
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
        (pp + 1)
        * (qq + 1)
    )

    vals = {}

    for ell in range(1, 10):
        vals[ell] = (
            pp * (qq + 1) ** ell
            + qq * (pp + 1) ** ell
            - pp ** ell * (qq + 1)
            - qq ** ell * (pp + 1)
        )

    h2_num = (
        vals[2] * vals[4]
        - vals[3] ** 2
    )

    h3_num = int(
        sp.Matrix(
            [
                [vals[2], vals[3], vals[4]],
                [vals[3], vals[4], vals[5]],
                [vals[4], vals[5], vals[6]],
            ]
        ).det()
    )

    delta_num = (
        SS ** 2
        - 4 * NN
    )

    F2_num = vals[2]

    if delta_num == 0:
        numerical_failures.append(
            (
                pp,
                qq,
                "Delta=0",
            )
        )
        continue

    denominator = (
        h2_num
        * NN
        * delta_num
        * (delta_num + 1)
    )

    if denominator == 0:
        numerical_failures.append(
            (
                pp,
                qq,
                "zero reconstruction denominator",
            )
        )
        continue

    recovered_M1 = sp.Rational(
        -h3_num * F2_num,
        denominator,
    )

    ok = (
        recovered_M1 == M1_num
    )

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N = {NN}"
    )

    print(
        f"    S = {SS}"
    )

    print(
        f"    M1 = {M1_num}"
    )

    print(
        f"    Delta = {delta_num}"
    )

    print(
        f"    recovered M1 = {recovered_M1}"
    )

    print(
        f"    PASS = {ok}"
    )

    if not ok:
        numerical_failures.append(
            (
                pp,
                qq,
                "M1 recovery failed",
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
# [17] INFORMATION MODEL
# ============================================================================

print("[17] INFORMATION-MODEL CHECK")
print("-" * 78)

print(
    "Established:"
)

print(
    "  F2,F3 -> S+1"
)

print(
    "  F2,F3,N -> M1"
)

print(
    "  M1-N -> S"
)

print(
    "  S,N -> p,q"
)

print(
    "  H3 contains M1 exactly."
)

print()

print(
    "Still unresolved:"
)

print(
    "  N -> F2,F3"
)

print(
    "  N -> Hankel minors"
)

print(
    "  N -> M1"
)

print()


# ============================================================================
# [18] PROOF CERTIFICATES
# ============================================================================

print("[18] PROOF CERTIFICATES")
print("-" * 78)

certificate_1 = (
    H2_pass
)

certificate_2 = (
    H3_pass
)

certificate_3 = (
    candidate_M1_pass
)

certificate_4 = (
    H5_pass
)

certificate_5 = (
    chi_pass
)

certificate_6 = (
    ratio_pass
)

certificate_7 = (
    root_product_pass
)

print(
    f"  H2 factorization = {certificate_1}"
)

print(
    f"  H3 factorization = {certificate_2}"
)

print(
    f"  M1 reconstruction = {certificate_3}"
)

print(
    f"  Hankel rank certificate = {certificate_4}"
)

print(
    f"  characteristic polynomial = {certificate_5}"
)

print(
    f"  H3/H2 cancellation = {certificate_6}"
)

print(
    f"  shifted root product = {certificate_7}"
)

print()


# ============================================================================
# [19] FINAL STATUS
# ============================================================================

overall = (
    certificate_1
    and certificate_2
    and certificate_3
    and certificate_4
    and certificate_5
    and certificate_6
    and certificate_7
    and len(numerical_failures) == 0
)

print("[19] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  H2 exact factorization: "
    f"{certificate_1}"
)

print(
    f"  H3 exact factorization: "
    f"{certificate_2}"
)

print(
    f"  M1 exact reconstruction: "
    f"{certificate_3}"
)

print(
    f"  Hankel rank <= 4: "
    f"{certificate_4}"
)

print(
    f"  Four-base characteristic polynomial: "
    f"{certificate_5}"
)

print(
    f"  H3/H2 nuisance cancellation: "
    f"{certificate_6}"
)

print(
    f"  (p+1)(q+1)=N+S+1: "
    f"{certificate_7}"
)

print(
    f"  Numerical audit: "
    f"{len(numerical_failures) == 0}"
)

print(
    f"  OVERALL EXACT AUDIT = {overall}"
)

print()

print(
    "MAIN RESULT:"
)

print(
    "  The Hankel minors of the original KAPPA sequence"
)

print(
    "  contain the shifted product M1=(p+1)(q+1)"
)

print(
    "  as an exact algebraic factor."
)

print()

print(
    "The structural chain is:"
)

print(
    "  F_l sequence"
)

print(
    "     -> Hankel minors"
)

print(
    "     -> M1=N+S+1"
)

print(
    "     -> S=M1-N"
)

print(
    "     -> z^2-S*z+N"
)

print(
    "     -> p,q"
)

print()

print(
    "REMAINING BRIDGE:"
)

print(
    "  Determine whether the Hankel minors, or M1 itself,"
)

print(
    "  admit a construction from N alone."
)

print()

print("=" * 78)
print("EXPERIMENT 485 FINISHED")
print("=" * 78)