#!/usr/bin/env python3

"""
EXPERIMENT 481
==============================================================================
ORDER-FOUR RECURRENCE OF THE ORIGINAL KAPPA KERNEL

Purpose
-------

Experiment 480 established the exact shifted-kernel coefficient theorem

    [x^(ell-1)] H_(k,ell)(p,q;x) = p^k + q^k,

and therefore, for k=1,

    [x^(ell-1)] H_(1,ell) = S = p+q.

The unresolved question is whether the original x=1 KAPPA kernel itself,
without introducing an explicit x-variable, contains enough information to
recover S.

For k=1 define

    F_l
      =
      p(1+q)^l + q(1+p)^l
      - p^l(1+q) - q^l(1+p).

Each term is a linear combination of the four exponential sequences

    p^l,
    q^l,
    (p+1)^l,
    (q+1)^l.

Therefore F_l must satisfy an exact order-four recurrence whose characteristic
polynomial is

    (z-p)(z-q)(z-p-1)(z-q-1).

Because

    N = p*q
    S = p+q,

the characteristic polynomial is entirely expressible in N and S.

This experiment:

    1. Derives the exact characteristic polynomial.
    2. Derives the exact order-four recurrence for F_l.
    3. Verifies it symbolically for many l.
    4. Checks whether the sequence has smaller recurrence order.
    5. Treats S as an unknown T and constructs recurrence residuals.
    6. Tests whether those residuals have the exact factor T-S.
    7. Computes gcds of multiple residual equations to determine whether
       S is algebraically isolated by the original kernel sequence.
    8. Tests exact numerical semiprime instances.
    9. Records the resulting bridge:

          N + consecutive original KAPPA values
              -> recurrence equations
              -> S
              -> p,q.

The experiment does NOT:
    - factor N,
    - use floating point,
    - interpolate from numerical data,
    - enumerate candidate factors,
    - write CSV files,
    - use statistical fitting.

Everything is exact symbolic arithmetic.
"""

from __future__ import annotations

import sympy as sp


# ============================================================================
# START
# ============================================================================

print("EXPERIMENT 481 START")
print("=" * 78)
print("ORDER-FOUR RECURRENCE OF THE ORIGINAL KAPPA KERNEL")
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
    """
    Exact p,q -> S,N reduction.
    """
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization failed: remainder={remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if exact_zero(concrete - (p + q)):
            result = result.subs(formal_var, S)
        elif exact_zero(concrete - p*q):
            result = result.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
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


def F_original(k: int, ell: int):
    """
    Original KAPPA kernel at k=1.

        F_l =
            p(1+q)^l + q(1+p)^l
            - p^l(1+q) - q^l(1+p)
    """
    if k != 1:
        raise ValueError("Experiment 481 is specialized to k=1.")

    return sp.expand(
        p * (1 + q)**ell
        + q * (1 + p)**ell
        - p**ell * (1 + q)
        - q**ell * (1 + p)
    )


def F_original_NS(ell: int):
    return symmetric_reduce(F_original(1, ell))


def characteristic_polynomial_pq():
    """
    Characteristic polynomial

        (z-p)(z-q)(z-(p+1))(z-(q+1))

    reduced to N,S.
    """
    raw = sp.expand(
        (z - p)
        * (z - q)
        * (z - p - 1)
        * (z - q - 1)
    )

    return sp.expand(symmetric_reduce(raw))


def recurrence_coefficients():
    """
    If

        chi(z) = z^4 - a3 z^3 + a2 z^2 - a1 z + a0,

    then

        F_(l+4)
          = a3 F_(l+3)
          - a2 F_(l+2)
          + a1 F_(l+1)
          - a0 F_l.
    """
    chi = characteristic_polynomial_pq()
    poly = sp.Poly(chi, z)

    a4 = poly.coeff_monomial(z**4)
    a3 = poly.coeff_monomial(z**3)
    a2 = poly.coeff_monomial(z**2)
    a1 = poly.coeff_monomial(z)
    a0 = poly.coeff_monomial(1)

    if a4 != 1:
        raise AssertionError("Characteristic polynomial not monic.")

    return (
        sp.expand(-a3),
        sp.expand(a2),
        sp.expand(-a1),
        sp.expand(a0),
    )


def recurrence_residual(
    ell: int,
    sequence: dict[int, sp.Expr],
    coeffs,
):
    """
    Compute

      F_(ell+4)
      - a3 F_(ell+3)
      + a2 F_(ell+2)
      - a1 F_(ell+1)
      + a0 F_ell
    """
    a3, a2, a1, a0 = coeffs

    return sp.expand(
        sequence[ell + 4]
        - a3 * sequence[ell + 3]
        + a2 * sequence[ell + 2]
        - a1 * sequence[ell + 1]
        + a0 * sequence[ell]
    )


# ============================================================================
# [1] Characteristic polynomial
# ============================================================================

print("[1] CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi_raw = sp.expand(
    (z - p)
    * (z - q)
    * (z - p - 1)
    * (z - q - 1)
)

chi_NS = characteristic_polynomial_pq()

print("In factorized p,q form:")
print(
    "  chi(z) = "
    "(z-p)(z-q)(z-p-1)(z-q-1)"
)
print()

print("In N,S coordinates:")
print(
    f"  chi(z) = {sp.factor(chi_NS)}"
)
print()

expected_chi = sp.expand(
    z**4
    - (2*S + 2)*z**3
    + (2*N + S**2 + 2*S + 1)*z**2
    - (2*N*S + S**2 + S + 2*N)*z
    + N*(N + S + 1)
)

chi_check = exact_zero(
    chi_NS - expected_chi
)

print(
    "Expected closed form:"
)
print(
    f"  {sp.factor(expected_chi)}"
)
print(
    f"CHARACTERISTIC POLYNOMIAL PASS = {chi_check}"
)
print()


# ============================================================================
# [2] Exact recurrence
# ============================================================================

print("[2] EXACT ORDER-FOUR RECURRENCE")
print("-" * 78)

a3, a2, a1, a0 = recurrence_coefficients()

print("The recurrence is")
print()
print(
    "  F_(l+4)"
)
print(
    "    = a3*F_(l+3)"
)
print(
    "      - a2*F_(l+2)"
)
print(
    "      + a1*F_(l+1)"
)
print(
    "      - a0*F_l"
)
print()

print(f"  a3 = {sp.factor(a3)}")
print(f"  a2 = {sp.factor(a2)}")
print(f"  a1 = {sp.factor(a1)}")
print(f"  a0 = {sp.factor(a0)}")
print()


# ============================================================================
# [3] Build original kernel sequence in N,S coordinates
# ============================================================================

print("[3] ORIGINAL KAPPA SEQUENCE")
print("-" * 78)

sequence = {}

for ell in range(1, ELL_MAX + 1):
    sequence[ell] = F_original_NS(ell)

for ell in range(1, min(8, ELL_MAX + 1)):
    print(
        f"  F_{ell} = {sp.factor(sequence[ell])}"
    )

print()


# ============================================================================
# [4] Symbolic recurrence verification
# ============================================================================

print("[4] SYMBOLIC RECURRENCE VERIFICATION")
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
# [5] Direct exponential proof
# ============================================================================

print("[5] DIRECT EXPONENTIAL-DECOMPOSITION PROOF")
print("-" * 78)

print(
    "Rewrite F_l as a linear combination of four exponential sequences:"
)
print()
print(
    "  F_l = p*(q+1)^l"
)
print(
    "        + q*(p+1)^l"
)
print(
    "        - (q+1)*p^l"
)
print(
    "        - (p+1)*q^l"
)
print()

print(
    "The four bases are:"
)
print(
    "  p, q, p+1, q+1."
)
print()

print(
    "Therefore the characteristic polynomial must be"
)
print(
    "  (z-p)(z-q)(z-p-1)(z-q-1)."
)
print()

# Verify exact equality of the rewritten expression for representative ell.
decomposition_failures = []

for ell in range(1, 10):
    direct = F_original(1, ell)

    decomposed = sp.expand(
        p * (q + 1)**ell
        + q * (p + 1)**ell
        - (q + 1) * p**ell
        - (p + 1) * q**ell
    )

    ok = exact_zero(direct - decomposed)

    print(
        f"  ell={ell:2d}: decomposition PASS={ok}"
    )

    if not ok:
        decomposition_failures.append(ell)

print(
    f"  DECOMPOSITION FAILURES = "
    f"{len(decomposition_failures)}"
)
print()


# ============================================================================
# [6] Treat S as unknown T
# ============================================================================

print("[6] UNKNOWN-S RECURRENCE RESIDUALS")
print("-" * 78)

print(
    "Now replace the true S in the recurrence coefficients by an"
)
print(
    "independent symbol T."
)
print()
print(
    "If the observed sequence is generated by the actual p,q pair,"
)
print(
    "the true value T=S must make every recurrence residual vanish."
)
print()

unknown_coefficients = (
    2*T + 2,
    2*N + T**2 + 2*T + 1,
    2*N*T + T**2 + T + 2*N,
    N*(N + T + 1),
)

residuals_T = {}

for ell in range(1, 7):
    residual = sp.expand(
        sequence[ell + 4]
        - unknown_coefficients[0] * sequence[ell + 3]
        + unknown_coefficients[1] * sequence[ell + 2]
        - unknown_coefficients[2] * sequence[ell + 1]
        + unknown_coefficients[3] * sequence[ell]
    )

    residual = sp.factor(residual)

    residuals_T[ell] = residual

    print(
        f"  ell={ell}:"
    )
    print(
        f"    residual = {residual}"
    )

print()


# ============================================================================
# [7] Verify T-S divisibility
# ============================================================================

print("[7] EXACT T-S FACTOR AUDIT")
print("-" * 78)

factor_failures = []

for ell, residual in residuals_T.items():
    quotient, remainder = sp.div(
        sp.Poly(
            sp.expand(residual),
            T,
            domain=sp.QQ.frac_field(N, S),
        ),
        sp.Poly(
            T - S,
            T,
            domain=sp.QQ.frac_field(N, S),
        ),
    )

    quotient = sp.factor(quotient.as_expr())
    remainder = sp.factor(remainder.as_expr())

    ok = exact_zero(remainder)

    print(
        f"  ell={ell}: "
        f"(T-S) divides residual -> {ok}"
    )

    if ok:
        print(
            f"    quotient = {quotient}"
        )
    else:
        print(
            f"    remainder = {remainder}"
        )
        factor_failures.append(ell)

print(
    f"  T-S DIVISIBILITY FAILURES = "
    f"{len(factor_failures)}"
)
print()


# ============================================================================
# [8] Common gcd of multiple residual equations
# ============================================================================

print("[8] COMMON RESIDUAL GCD")
print("-" * 78)

residual_polys = []

for ell, residual in residuals_T.items():
    residual_poly = sp.Poly(
        sp.expand(residual),
        T,
        domain=sp.QQ.frac_field(N, S),
    )

    residual_polys.append(residual_poly)

if residual_polys:
    gcd_poly = residual_polys[0]

    for poly in residual_polys[1:]:
        gcd_poly = sp.gcd(gcd_poly, poly)

    gcd_expr = sp.factor(gcd_poly.as_expr())

    print(
        "Common gcd of all tested recurrence residuals:"
    )
    print(
        f"  gcd = {gcd_expr}"
    )

    gcd_expected = T - S

    # Normalize by checking whether gcd/(T-S) is nonzero constant
    # in the field Q(N,S).
    quotient = sp.cancel(
        gcd_expr / gcd_expected
    )

    print(
        f"  gcd/(T-S) = {sp.factor(quotient)}"
    )

    gcd_exact = (
        sp.Poly(
            gcd_expr,
            T,
            domain=sp.QQ.frac_field(N, S),
        ).degree() == 1
        and exact_zero(
            sp.rem(
                sp.Poly(
                    gcd_expr,
                    T,
                    domain=sp.QQ.frac_field(N, S),
                ),
                sp.Poly(
                    T - S,
                    T,
                    domain=sp.QQ.frac_field(N, S),
                ),
            ).as_expr()
        )
    )

    print(
        f"  COMMON GCD IS EXACTLY LINEAR = {gcd_exact}"
    )

print()


# ============================================================================
# [9] Candidate S recovery from a single recurrence equation
# ============================================================================

print("[9] ALGEBRAIC RECOVERY OF S FROM ONE RECURRENCE")
print("-" * 78)

print(
    "Each recurrence equation is a polynomial equation in T."
)
print(
    "The experiment factors the first residual explicitly."
)
print()

first_residual = residuals_T[1]

print(
    f"  R_1(T) = {sp.factor(first_residual)}"
)
print()

quotient_R1, remainder_R1 = sp.div(
    sp.Poly(
        first_residual,
        T,
        domain=sp.QQ.frac_field(N, S),
    ),
    sp.Poly(
        T - S,
        T,
        domain=sp.QQ.frac_field(N, S),
    ),
)

print(
    "  R_1(T)/(T-S) ="
)
print(
    f"    {sp.factor(quotient_R1.as_expr())}"
)
print(
    f"  remainder = {sp.factor(remainder_R1.as_expr())}"
)
print()


# ============================================================================
# [10] Is the alternate root useful?
# ============================================================================

print("[10] ALTERNATE-ROOT ANALYSIS")
print("-" * 78)

print(
    "The first recurrence residual may contain additional roots in T."
)
print(
    "We factor R_1(T) and inspect its remaining factor."
)
print()

R1_factor = sp.factor(first_residual)

print(
    f"  R_1(T) = {R1_factor}"
)
print()

factor_list_R1 = sp.factor_list(first_residual)

print(
    "  Factor list:"
)

for factor_expr, multiplicity in factor_list_R1[1]:
    print(
        f"    multiplicity={multiplicity}: "
        f"{sp.factor(factor_expr)}"
    )

print()


# ============================================================================
# [11] Cross-equation uniqueness
# ============================================================================

print("[11] CROSS-EQUATION ROOT INTERSECTION")
print("-" * 78)

print(
    "If a single recurrence equation contains multiple algebraic T roots,"
)
print(
    "common gcds across different ell values determine which roots survive."
)
print()

for i in range(len(residual_polys) - 1):
    pair_gcd = sp.gcd(
        residual_polys[i],
        residual_polys[i + 1],
    )

    print(
        f"  residuals ell={i+1} and ell={i+2}:"
    )
    print(
        f"    gcd = {sp.factor(pair_gcd.as_expr())}"
    )

print()


# ============================================================================
# [12] Numerical exact recovery
# ============================================================================

print("[12] NUMERICAL EXACT RECOVERY OF S")
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

    # Construct exact integer kernel values.
    Fvals = {}

    for ell in range(1, 10):
        Fvals[ell] = (
            pp * (1 + qq)**ell
            + qq * (1 + pp)**ell
            - pp**ell * (1 + qq)
            - qq**ell * (1 + pp)
        )

    # Numerical unknown T polynomial for ell=1.
    T_num = sp.symbols("T_num")

    R_num = sp.expand(
        Fvals[5]
        - (2*T_num + 2) * Fvals[4]
        + (2*NN + T_num**2 + 2*T_num + 1) * Fvals[3]
        - (2*NN*T_num + T_num**2 + T_num + 2*NN) * Fvals[2]
        + NN*(NN + T_num + 1) * Fvals[1]
    )

    factored_num = sp.factor(R_num)

    print(
        f"  recurrence residual polynomial:"
    )
    print(
        f"    {factored_num}"
    )

    roots = sp.solve(
        sp.Eq(R_num, 0),
        T_num,
    )

    print(
        f"  candidate roots from ell=1 recurrence:"
    )

    for root in roots:
        print(
            f"    T = {root}"
        )

    # Determine whether true S is a root.
    true_ok = exact_zero(
        R_num.subs(T_num, SS)
    )

    if not true_ok:
        numerical_failures.append(
            (pp, qq, "true S is not a recurrence root")
        )

    print(
        f"  true S satisfies recurrence = {true_ok}"
    )

    # Cross-check with second recurrence.
    R2_num = sp.expand(
        Fvals[6]
        - (2*T_num + 2) * Fvals[5]
        + (2*NN + T_num**2 + 2*T_num + 1) * Fvals[4]
        - (2*NN*T_num + T_num**2 + T_num + 2*NN) * Fvals[3]
        + NN*(NN + T_num + 1) * Fvals[2]
    )

    joint_gcd = sp.gcd(
        sp.Poly(R_num, T_num),
        sp.Poly(R2_num, T_num),
    )

    print(
        f"  gcd(first,second recurrence) = "
        f"{sp.factor(joint_gcd.as_expr())}"
    )

    if joint_gcd.degree() != 1:
        numerical_failures.append(
            (pp, qq, "joint recurrence gcd not linear")
        )
    else:
        # Check whether joint gcd is proportional to T-trueS.
        normalized = sp.factor(
            joint_gcd.as_expr() / (T_num - SS)
        )

        print(
            f"  gcd/(T-S) = {normalized}"
        )

    print()


print(
    f"NUMERICAL RECOVERY FAILURES = "
    f"{len(numerical_failures)}"
)
if numerical_failures:
    for failure in numerical_failures:
        print(
            f"  FAILURE: {failure}"
        )

print()


# ============================================================================
# [13] Minimal-recurrence audit
# ============================================================================

print("[13] SMALLER RECURRENCE AUDIT")
print("-" * 78)

print(
    "Test whether the k=1 original kernel sequence satisfies a universal"
)
print(
    "recurrence of order 1, 2, or 3 over Q(N,S)."
)
print()

smaller_failures = []

# Build a symbolic sequence long enough.
SEQ = {
    ell: sequence[ell]
    for ell in range(1, ELL_MAX + 1)
}

# Order 1: F_(l+1) = c0 F_l with c0 in Q(N,S)
#
# A universal order-1 recurrence would require
# F_(l+1)/F_l to be independent of ell.
def test_order_one():
    ratios = []

    for ell in range(1, 6):
        ratio = sp.cancel(
            SEQ[ell + 1] / SEQ[ell]
        )
        ratios.append(ratio)

    for i in range(1, len(ratios)):
        if not exact_zero(ratios[i] - ratios[0]):
            return False

    return True


order1 = test_order_one()

print(
    f"  Order 1 recurrence exists = {order1}"
)

# Order 2 / order 3 are tested through exact linear-system compatibility.
#
# We require constants c_j in Q(N,S), independent of ell.
# Solve using the first equations and verify all remaining ones.

def test_constant_coefficient_order(order: int):
    coeff_symbols = sp.symbols(
        f"c0:{order}"
    )

    equations = []

    # F_(l+order) =
    #   c0 F_(l+order-1) + ... + c_(order-1) F_l
    for ell in range(
        1,
        min(ELL_MAX - order + 1, 8)
    ):
        lhs = SEQ[ell + order]

        rhs = sum(
            coeff_symbols[j]
            * SEQ[ell + order - 1 - j]
            for j in range(order)
        )

        equations.append(
            sp.Eq(lhs, rhs)
        )

    solved = sp.solve(
        equations,
        coeff_symbols,
        dict=True,
    )

    if not solved:
        return False, None

    candidate = solved[0]

    for ell in range(
        1,
        ELL_MAX - order + 1
    ):
        lhs = sp.expand(SEQ[ell + order])

        rhs = sp.expand(
            sum(
                candidate[coeff_symbols[j]]
                * SEQ[ell + order - 1 - j]
                for j in range(order)
            )
        )

        if not exact_zero(lhs - rhs):
            return False, candidate

    return True, candidate


for order in [2, 3]:
    exists, candidate = test_constant_coefficient_order(order)

    print(
        f"  Order {order} recurrence exists = {exists}"
    )

    if exists:
        print(
            f"    coefficients = {candidate}"
        )

print()


# ============================================================================
# [14] Bridge theorem
# ============================================================================

print("[14] BRIDGE THEOREM CANDIDATE")
print("-" * 78)

print(
    "The exact original-kernel sequence"
)
print(
    "  F_l = p(1+q)^l + q(1+p)^l"
)
print(
    "        - p^l(1+q) - q^l(1+p)"
)
print(
    "is a linear combination of four exponentials."
)
print()
print(
    "Its characteristic polynomial is"
)
print(
    "  (z-p)(z-q)(z-p-1)(z-q-1)."
)
print()
print(
    "This polynomial is expressible entirely through N=pq and S=p+q."
)
print()
print(
    "Therefore consecutive original-kernel values satisfy an exact"
)
print(
    "order-four recurrence whose coefficients are functions of N and S."
)
print()
print(
    "If N is known and enough consecutive F_l values are independently"
)
print(
    "available, the unknown S can be recovered algebraically from the"
)
print(
    "recurrence residual equations."
)
print()


# ============================================================================
# [15] Information-model distinction
# ============================================================================

print("[15] INFORMATION-MODEL DISTINCTION")
print("-" * 78)

print(
    "Established by this experiment:"
)
print(
    "  F_l sequence -> exact order-four recurrence -> S."
)
print()
print(
    "Still unresolved:"
)
print(
    "  N alone -> independently computable F_l sequence."
)
print()
print(
    "This is now the precise remaining bridge."
)
print()


# ============================================================================
# [16] Final status
# ============================================================================

all_exact = (
    chi_check
    and len(recurrence_failures) == 0
    and len(decomposition_failures) == 0
    and len(factor_failures) == 0
    and len(numerical_failures) == 0
)

print("[16] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  Characteristic polynomial identity: {chi_check}"
)
print(
    f"  Symbolic recurrence: "
    f"{len(recurrence_failures) == 0}"
)
print(
    f"  Four-exponential decomposition: "
    f"{len(decomposition_failures) == 0}"
)
print(
    f"  Exact T-S residual factor: "
    f"{len(factor_failures) == 0}"
)
print(
    f"  Numerical recovery checks: "
    f"{len(numerical_failures) == 0}"
)
print(
    f"  OVERALL EXACT AUDIT: {all_exact}"
)

print()
print(
    "MAIN RESULT:"
)
print(
    "  The original x=1 KAPPA kernel is itself a finite exponential"
)
print(
    "  sequence in ell with a characteristic polynomial depending only"
)
print(
    "  on N and S."
)
print()
print(
    "  This gives a new exact bridge:"
)
print(
    "      original KAPPA kernel sequence"
)
print(
    "            -> order-four characteristic recurrence"
)
print(
    "            -> S=p+q"
)
print(
    "            -> p,q."
)
print()
print(
    "  The remaining question is whether the required consecutive kernel"
)
print(
    "  values can be generated from N without prior access to p,q."
)
print()


# ============================================================================
# FINISH
# ============================================================================

print("=" * 78)
print("EXPERIMENT 481 FINISHED")
print("=" * 78)
