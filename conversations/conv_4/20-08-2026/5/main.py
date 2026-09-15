#!/usr/bin/env python3

"""
EXPERIMENT 482
==============================================================================
CORRECTED RECURRENCE ROOT ISOLATION FOR S = p+q

Experiment 481 established that the original k=1 KAPPA kernel sequence

    F_l =
        p(q+1)^l
        + q(p+1)^l
        - (q+1)p^l
        - (p+1)q^l

is a four-exponential sequence.

Two bookkeeping errors were found in Experiment 481:

    1. The z^2 characteristic coefficient was printed as
           2N + S^2 + 2S + 1
       instead of
           2N + S^2 + 3S + 1.

    2. The unknown-T numerical recurrence used the same missing-S error.

This experiment corrects those errors.

Primary target
--------------

Let

    N = pq
    S = p+q.

The exact characteristic polynomial is

    chi(z)
      = (z-p)(z-q)(z-p-1)(z-q-1)

      = z^4
        - 2(S+1) z^3
        + (2N+S^2+3S+1) z^2
        - (2NS+2N+S^2+S) z
        + N(N+S+1).

Replace S by an independent variable T.

For each recurrence window define the residual R_l(T).

At the true value T=S,

    R_l(S)=0.

Therefore every R_l(T) must contain the factor (T-S).

The experiment asks:

    1. Does the corrected recurrence hold symbolically?
    2. Does every residual contain T-S?
    3. Is gcd(R_1,R_2,...) exactly T-S?
    4. Does this isolate S algebraically?
    5. Does numerical exact recovery work for independent semiprimes?
    6. What is the remaining algebraic factor in a single recurrence?
    7. Does combining two recurrence windows remove the extraneous root?

This is still NOT an N-only factorization algorithm.

The information model being tested is:

    N + consecutive exact F_l values
        -> S
        -> p,q.

No factorization of N is used.
No floating point.
No CSV.
No interpolation.
"""

from __future__ import annotations

import sympy as sp


# ============================================================================
# START
# ============================================================================

print("EXPERIMENT 482 START")
print("=" * 78)
print("CORRECTED RECURRENCE ROOT ISOLATION FOR S = p+q")
print("=" * 78)
print()


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
N, S, T, z = sp.symbols("N S T z")


# ============================================================================
# Test configuration
# ============================================================================

ELL_MAX = 18

TEST_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
]


# ============================================================================
# Helpers
# ============================================================================

def exact_zero(expr) -> bool:
    return sp.cancel(sp.expand(expr)) == 0


def symmetric_reduce(expr):
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization failed: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if exact_zero(concrete - (p + q)):
            result = result.subs(formal_var, S)
        elif exact_zero(concrete - p*q):
            result = result.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs({
            S: p + q,
            N: p*q,
        }) - expr
    )

    if not exact_zero(check):
        raise AssertionError(
            "Symmetric reconstruction failed."
        )

    return result


def F_raw(ell: int):
    return sp.expand(
        p * (1 + q)**ell
        + q * (1 + p)**ell
        - p**ell * (1 + q)
        - q**ell * (1 + p)
    )


def F_NS(ell: int):
    return symmetric_reduce(F_raw(ell))


def characteristic_polynomial():
    raw = sp.expand(
        (z - p)
        * (z - q)
        * (z - p - 1)
        * (z - q - 1)
    )

    return sp.expand(
        symmetric_reduce(raw)
    )


def corrected_coefficients():
    """
    chi(z) =
        z^4
        - a3 z^3
        + a2 z^2
        - a1 z
        + a0
    """

    a3 = sp.expand(
        2*(S + 1)
    )

    a2 = sp.expand(
        2*N + S**2 + 3*S + 1
    )

    a1 = sp.expand(
        (2*N + S) * (S + 1)
    )

    a0 = sp.expand(
        N*(N + S + 1)
    )

    return a3, a2, a1, a0


def unknown_T_coefficients():
    """
    Exactly the same recurrence coefficients with S -> T.
    """

    a3 = sp.expand(
        2*(T + 1)
    )

    a2 = sp.expand(
        2*N + T**2 + 3*T + 1
    )

    a1 = sp.expand(
        (2*N + T) * (T + 1)
    )

    a0 = sp.expand(
        N*(N + T + 1)
    )

    return a3, a2, a1, a0


def recurrence_residual(
    ell: int,
    sequence: dict[int, sp.Expr],
    coefficients,
):
    a3, a2, a1, a0 = coefficients

    return sp.expand(
        sequence[ell + 4]
        - a3 * sequence[ell + 3]
        + a2 * sequence[ell + 2]
        - a1 * sequence[ell + 1]
        + a0 * sequence[ell]
    )


# ============================================================================
# [1] Correct characteristic polynomial
# ============================================================================

print("[1] CORRECT CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi_factorized = sp.expand(
    (z - p)
    * (z - q)
    * (z - p - 1)
    * (z - q - 1)
)

chi_NS = characteristic_polynomial()

chi_expected = sp.expand(
    z**4
    - 2*(S + 1)*z**3
    + (2*N + S**2 + 3*S + 1)*z**2
    - (2*N*S + 2*N + S**2 + S)*z
    + N*(N + S + 1)
)

print("Factorized:")
print(
    "  (z-p)(z-q)(z-p-1)(z-q-1)"
)
print()

print("Exact N,S form:")
print(
    f"  {sp.factor(chi_NS)}"
)
print()

print("Expected:")
print(
    f"  {sp.factor(chi_expected)}"
)
print()

chi_pass = exact_zero(
    chi_NS - chi_expected
)

print(
    f"CHARACTERISTIC IDENTITY PASS = {chi_pass}"
)
print()


# ============================================================================
# [2] Correct recurrence coefficients
# ============================================================================

print("[2] CORRECT RECURRENCE COEFFICIENTS")
print("-" * 78)

a3, a2, a1, a0 = corrected_coefficients()

print(
    f"  a3 = {sp.factor(a3)}"
)
print(
    f"  a2 = {sp.factor(a2)}"
)
print(
    f"  a1 = {sp.factor(a1)}"
)
print(
    f"  a0 = {sp.factor(a0)}"
)
print()

coefficient_check = (
    exact_zero(
        a3 - 2*(S + 1)
    )
    and exact_zero(
        a2 - (2*N + S**2 + 3*S + 1)
    )
    and exact_zero(
        a1 - (2*N + S)*(S + 1)
    )
    and exact_zero(
        a0 - N*(N + S + 1)
    )
)

print(
    f"RECURRENCE COEFFICIENT PASS = {coefficient_check}"
)
print()


# ============================================================================
# [3] Build original sequence
# ============================================================================

print("[3] BUILD ORIGINAL KAPPA SEQUENCE")
print("-" * 78)

sequence = {
    ell: F_NS(ell)
    for ell in range(1, ELL_MAX + 1)
}

for ell in range(1, 8):
    print(
        f"  F_{ell} = {sp.factor(sequence[ell])}"
    )

print()


# ============================================================================
# [4] Symbolic recurrence
# ============================================================================

print("[4] SYMBOLIC RECURRENCE AUDIT")
print("-" * 78)

recurrence_failures = []

for ell in range(1, ELL_MAX - 3):
    residual = recurrence_residual(
        ell,
        sequence,
        (a3, a2, a1, a0),
    )

    ok = exact_zero(residual)

    print(
        f"  ell={ell:2d}: PASS={ok}"
    )

    if not ok:
        recurrence_failures.append(ell)
        print(
            f"    residual = {sp.factor(residual)}"
        )

print(
    f"  RECURRENCE FAILURES = {len(recurrence_failures)}"
)
print()


# ============================================================================
# [5] Unknown-T recurrence
# ============================================================================

print("[5] UNKNOWN-T RECURRENCE")
print("-" * 78)

T_coeffs = unknown_T_coefficients()

print(
    "Unknown-T coefficients:"
)
print(
    f"  a3(T) = {sp.factor(T_coeffs[0])}"
)
print(
    f"  a2(T) = {sp.factor(T_coeffs[1])}"
)
print(
    f"  a1(T) = {sp.factor(T_coeffs[2])}"
)
print(
    f"  a0(T) = {sp.factor(T_coeffs[3])}"
)
print()


residuals = {}

for ell in range(1, 7):

    residual = sp.factor(
        recurrence_residual(
            ell,
            sequence,
            T_coeffs,
        )
    )

    residuals[ell] = residual

    print(
        f"  ell={ell}:"
    )
    print(
        f"    R_{ell}(T) = {residual}"
    )

print()


# ============================================================================
# [6] Exact T-S factor
# ============================================================================

print("[6] EXACT T-S FACTORIZATION")
print("-" * 78)

factor_failures = []

quotients = {}

for ell, residual in residuals.items():

    divisor = sp.Poly(
        T - S,
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    polynomial = sp.Poly(
        sp.expand(residual),
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    quotient, remainder = sp.div(
        polynomial,
        divisor,
    )

    quotient_expr = sp.factor(
        quotient.as_expr()
    )

    remainder_expr = sp.factor(
        remainder.as_expr()
    )

    ok = exact_zero(remainder_expr)

    print(
        f"  ell={ell}: "
        f"(T-S) divides R_{ell} -> {ok}"
    )

    if ok:
        quotients[ell] = quotient_expr
        print(
            f"    quotient = {quotient_expr}"
        )
    else:
        factor_failures.append(ell)
        print(
            f"    remainder = {remainder_expr}"
        )

print(
    f"  T-S FACTORIZATION FAILURES = "
    f"{len(factor_failures)}"
)
print()


# ============================================================================
# [7] Common gcd of residuals
# ============================================================================

print("[7] COMMON GCD OF RECURRENCE RESIDUALS")
print("-" * 78)

residual_polys = [
    sp.Poly(
        sp.expand(residuals[ell]),
        T,
        domain=sp.QQ.frac_field(N, S),
    )
    for ell in sorted(residuals)
]

common_gcd = residual_polys[0]

for polynomial in residual_polys[1:]:
    common_gcd = sp.gcd(
        common_gcd,
        polynomial,
    )

common_gcd_expr = sp.factor(
    common_gcd.as_expr()
)

print(
    f"  gcd(R_1,...,R_6) = {common_gcd_expr}"
)
print()

gcd_minus = sp.cancel(
    common_gcd_expr / (T - S)
)

print(
    f"  gcd/(T-S) = {sp.factor(gcd_minus)}"
)

common_gcd_is_T_minus_S = (
    common_gcd.degree() == 1
    and exact_zero(
        sp.rem(
            common_gcd,
            sp.Poly(
                T - S,
                T,
                domain=sp.QQ.frac_field(N, S),
            ),
        ).as_expr()
    )
)

print(
    f"  GCD IS EXACTLY T-S UP TO A UNIT = "
    f"{common_gcd_is_T_minus_S}"
)
print()


# ============================================================================
# [8] Single-equation extraneous root analysis
# ============================================================================

print("[8] SINGLE-RECURRENCE EXTRANEOUS ROOTS")
print("-" * 78)

for ell in [1, 2, 3]:
    polynomial = sp.Poly(
        residuals[ell],
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    print(
        f"  R_{ell}(T) factorization:"
    )

    print(
        f"    {sp.factor(polynomial.as_expr())}"
    )

    print()

    print(
        f"  quotient by T-S:"
    )

    print(
        f"    {sp.factor(quotients[ell])}"
    )

    print()


# ============================================================================
# [9] Cross-window gcds
# ============================================================================

print("[9] CROSS-WINDOW ROOT INTERSECTION")
print("-" * 78)

for ell in range(1, 6):

    left = sp.Poly(
        residuals[ell],
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    right = sp.Poly(
        residuals[ell + 1],
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    pair_gcd = sp.gcd(
        left,
        right,
    )

    print(
        f"  gcd(R_{ell},R_{ell+1}) = "
        f"{sp.factor(pair_gcd.as_expr())}"
    )

print()


# ============================================================================
# [10] Exact numerical reconstruction
# ============================================================================

print("[10] NUMERICAL EXACT S RECOVERY")
print("-" * 78)

numerical_failures = []

for pp, qq in TEST_PAIRS:

    NN = pp * qq
    SS = pp + qq

    print(
        f"--- p={pp}, q={qq} ---"
    )
    print(
        f"  N = {NN}"
    )
    print(
        f"  true S = {SS}"
    )

    Fvals = {}

    for ell in range(1, 10):
        Fvals[ell] = (
            pp * (qq + 1)**ell
            + qq * (pp + 1)**ell
            - pp**ell * (qq + 1)
            - qq**ell * (pp + 1)
        )

    T_num = sp.symbols("T_num")

    def numeric_residual(ell):
        a3n = 2*(T_num + 1)

        a2n = (
            2*NN
            + T_num**2
            + 3*T_num
            + 1
        )

        a1n = (
            2*NN*T_num
            + 2*NN
            + T_num**2
            + T_num
        )

        a0n = NN*(NN + T_num + 1)

        return sp.expand(
            Fvals[ell + 4]
            - a3n * Fvals[ell + 3]
            + a2n * Fvals[ell + 2]
            - a1n * Fvals[ell + 1]
            + a0n * Fvals[ell]
        )

    R1 = sp.factor(
        numeric_residual(1)
    )

    R2 = sp.factor(
        numeric_residual(2)
    )

    print(
        f"  R_1(T) = {R1}"
    )

    print(
        f"  R_1(S) = "
        f"{sp.simplify(R1.subs(T_num, SS))}"
    )

    true_pass = exact_zero(
        R1.subs(T_num, SS)
    )

    print(
        f"  true S is a root = {true_pass}"
    )

    gcd_num = sp.gcd(
        sp.Poly(R1, T_num),
        sp.Poly(R2, T_num),
    )

    gcd_num_expr = sp.factor(
        gcd_num.as_expr()
    )

    print(
        f"  gcd(R_1,R_2) = {gcd_num_expr}"
    )

    gcd_factor = sp.cancel(
        gcd_num_expr / (T_num - SS)
    )

    print(
        f"  gcd/(T-S) = {sp.factor(gcd_factor)}"
    )

    gcd_pass = (
        gcd_num.degree() == 1
        and exact_zero(
            gcd_num.eval(SS)
        )
    )

    print(
        f"  joint gcd isolates true S = {gcd_pass}"
    )

    if not true_pass:
        numerical_failures.append(
            (pp, qq, "true S failed")
        )

    if not gcd_pass:
        numerical_failures.append(
            (pp, qq, "joint gcd failed")
        )

    print()


print(
    f"  NUMERICAL FAILURES = "
    f"{len(numerical_failures)}"
)

if numerical_failures:
    for failure in numerical_failures:
        print(
            f"    {failure}"
        )

print()


# ============================================================================
# [11] Direct algebraic recovery from two windows
# ============================================================================

print("[11] TWO-WINDOW ALGEBRAIC RECOVERY")
print("-" * 78)

print(
    "The central question is whether two recurrence windows"
)
print(
    "remove every root other than T=S."
)
print()

two_window_gcd = sp.gcd(
    sp.Poly(
        residuals[1],
        T,
        domain=sp.QQ.frac_field(N, S),
    ),
    sp.Poly(
        residuals[2],
        T,
        domain=sp.QQ.frac_field(N, S),
    ),
)

print(
    f"  gcd(R_1,R_2) = "
    f"{sp.factor(two_window_gcd.as_expr())}"
)

print()


# ============================================================================
# [12] Important distinction: this still requires kernel values
# ============================================================================

print("[12] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print(
    "If N and F_1,...,F_6 are available:"
)
print()
print(
    "  recurrence residuals -> gcd in T -> S"
)
print(
    "                         -> z^2-S*z+N"
)
print(
    "                         -> p,q"
)
print()

print(
    "But this experiment does NOT establish:"
)
print()
print(
    "  N -> F_1,...,F_6"
)
print()

print(
    "Thus the remaining bridge is now:"
)
print()
print(
    "  N"
)
print(
    "   -> independently computable original KAPPA sequence"
)
print(
    "   -> recurrence equations"
)
print(
    "   -> S"
)
print(
    "   -> factor pair"
)
print()


# ============================================================================
# [13] Proof certificate
# ============================================================================

print("[13] SYMBOLIC PROOF CERTIFICATE")
print("-" * 78)

print(
    "The characteristic polynomial follows directly from"
)
print()
print(
    "  chi(z)"
)
print(
    "    = (z^2-Sz+N)"
)
print(
    "      * (z^2-(S+2)z+(N+S+1))."
)
print()

print(
    "The second factor is the characteristic polynomial of"
)
print(
    "the shifted roots p+1 and q+1 because"
)
print()
print(
    "  (p+1)+(q+1)=S+2"
)
print(
    "  (p+1)(q+1)=N+S+1."
)
print()

print(
    "Therefore the order-four recurrence is exact."
)
print()

certificate_pass = (
    chi_pass
    and coefficient_check
    and len(recurrence_failures) == 0
    and len(factor_failures) == 0
    and common_gcd_is_T_minus_S
)

print(
    f"  SYMBOLIC CERTIFICATE PASS = {certificate_pass}"
)
print()


# ============================================================================
# [14] Final status
# ============================================================================

print("[14] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  Correct characteristic polynomial: {chi_pass}"
)
print(
    f"  Exact order-four recurrence: "
    f"{len(recurrence_failures) == 0}"
)
print(
    f"  Every residual divisible by T-S: "
    f"{len(factor_failures) == 0}"
)
print(
    f"  Common gcd equals T-S: "
    f"{common_gcd_is_T_minus_S}"
)
print(
    f"  Numerical recovery: "
    f"{len(numerical_failures) == 0}"
)
print(
    f"  Overall symbolic certificate: "
    f"{certificate_pass}"
)

print()

print(
    "MAIN TARGET:"
)
print(
    "  Determine whether the original KAPPA sequence itself"
)
print(
    "  provides an algebraically identifiable S when N is known."
)
print()

print(
    "NEXT RESEARCH TARGET:"
)
print(
    "  Connect the independently generated sequence values F_l"
)
print(
    "  to the existing N,S homogeneous-layer representation."
)
print()

print("=" * 78)
print("EXPERIMENT 482 FINISHED")
print("=" * 78)
