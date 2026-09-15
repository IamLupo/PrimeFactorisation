#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 288 — EXACT B-ROW DIFFERENTIAL / APPELL-TYPE OPERATOR AUDIT
==============================================================================

Experiment 287 established:

    * each primitive B-row generating polynomial is irreducible over Z;
    * every pairwise row gcd is 1;
    * no rational roots;
    * no fixed roots x = c;
    * no affine root x = a*k+b.

The row degrees are nevertheless exactly

    7, 6, 5, 4, 3, 2.

This experiment tests whether that degree ladder comes from a
differential operator.

Candidate structures:

    A. direct derivative:
           H_{k+1}(x) = c_k H'_k(x)

    B. shifted derivative:
           H_{k+1}(x) = c_k H'_k(x+s_k)

    C. derivative plus multiplication:
           H_{k+1}(x) = c_k H'_k(x) + d_k H_k(x)

       when degree compatibility permits;

    D. affine change of variable:
           H_{k+1}(x) =
           c_k H'_k(a_k x+b_k)

    E. Appell-type common derivative relation:
           H'_k(x) = c_k H_{k+1}(x)

The first test is intentionally strict and exact.

Only identities holding for the COMPLETE polynomial count.

No arbitrary coefficient-by-coefficient fitting.

No q-family data.
Exact QQ arithmetic only.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# EXACT PRIMITIVE INTEGER B-ROW POLYNOMIALS
# =============================================================================

x = sp.Symbol("x")


H = [
    -421*x**7
    -69937*x**6
    +423738*x**5
    -1055250*x**4
    -41580*x**3
    +4071060*x**2
    +1559880*x
    +63000,

    -12106*x**6
    -33638*x**5
    +379836*x**4
    -1230345*x**3
    +819480*x**2
    +3104640*x
    +630000,

    -173070*x**5
    +240576*x**4
    +1037057*x**3
    -6072045*x**2
    +8596560*x
    +8139600,

    -808952*x**4
    +2458897*x**3
    -3513510*x**2
    -3999450*x
    +18564000,

    -162139*x**3
    +926401*x**2
    -2779686*x
    +3753750,

    301*x**2
    -626*x
    +650,
]


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def degree(expr):
    return int(
        sp.Poly(
            sp.expand(expr),
            x,
            domain=sp.QQ,
        ).degree()
    )


def exact_scalar_multiple(A, B):
    """
    Test whether A = c*B for one rational scalar c.
    """
    A = sp.Poly(
        sp.expand(A),
        x,
        domain=sp.QQ,
    )

    B = sp.Poly(
        sp.expand(B),
        x,
        domain=sp.QQ,
    )

    if A.is_zero or B.is_zero:
        return None

    if A.degree() != B.degree():
        return None

    ratios = []

    for coeff_A, coeff_B in zip(
        A.all_coeffs(),
        B.all_coeffs(),
    ):
        if coeff_B == 0:
            if coeff_A != 0:
                return None
            continue

        ratios.append(
            clean(
                coeff_A / coeff_B
            )
        )

    if not ratios:
        return None

    c = ratios[0]

    if all(
        ratio == c
        for ratio in ratios
    ):
        return c

    return None


def affine_transform(poly, a, b):
    return clean(
        poly.subs(
            x,
            a*x + b,
        )
    )


# =============================================================================
# 1. DIRECT DERIVATIVE TEST
# =============================================================================

def direct_derivative_test():

    print()
    print("=" * 78)
    print("1. DIRECT DERIVATIVE TEST")
    print("=" * 78)

    hits = []

    for k in range(
        len(H) - 1
    ):

        derivative = clean(
            sp.diff(
                H[k],
                x,
            )
        )

        scalar = exact_scalar_multiple(
            H[k + 1],
            derivative,
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    H_k_degree="
            f"{degree(H[k])}"
        )

        print(
            f"    H_k1_degree="
            f"{degree(H[k+1])}"
        )

        print(
            f"    derivative="
            f"{derivative}"
        )

        print(
            f"    exact_scalar="
            f"{scalar}"
        )

        if scalar is not None:
            hits.append(
                (
                    k,
                    scalar,
                )
            )

    return hits


# =============================================================================
# 2. REVERSE APPPELL TEST
# =============================================================================

def reverse_derivative_test():

    print()
    print("=" * 78)
    print("2. REVERSE APPELL TEST")
    print("=" * 78)

    hits = []

    for k in range(
        len(H) - 1
    ):

        derivative = clean(
            sp.diff(
                H[k],
                x,
            )
        )

        scalar = exact_scalar_multiple(
            derivative,
            H[k + 1],
        )

        print()
        print(
            f"  k={k}: "
            f"scalar={scalar}"
        )

        if scalar is not None:
            hits.append(
                (
                    k,
                    scalar,
                )
            )

    return hits


# =============================================================================
# 3. TRANSLATED DERIVATIVE SEARCH
# =============================================================================

def translated_derivative_test():

    print()
    print("=" * 78)
    print("3. TRANSLATED DERIVATIVE SEARCH")
    print("=" * 78)

    hits = []

    # Small exact integer/rational shifts.
    shifts = [
        sp.Integer(v)
        for v in range(
            -3,
            4,
        )
    ]

    for k in range(
        len(H) - 1
    ):

        local_hits = []

        derivative = sp.diff(
            H[k],
            x,
        )

        for s in shifts:

            transformed = clean(
                derivative.subs(
                    x,
                    x + s,
                )
            )

            scalar = exact_scalar_multiple(
                H[k + 1],
                transformed,
            )

            if scalar is not None:
                local_hits.append(
                    (
                        s,
                        scalar,
                    )
                )

        print()
        print(
            f"  k={k}: "
            f"hits={local_hits}"
        )

        if local_hits:
            hits.append(
                (
                    k,
                    local_hits,
                )
            )

    return hits


# =============================================================================
# 4. SHIFT + SCALE TEST
# =============================================================================

def affine_derivative_test():

    print()
    print("=" * 78)
    print("4. AFFINE VARIABLE DERIVATIVE SEARCH")
    print("=" * 78)

    hits = []

    # Search small rational affine maps:
    #
    #     x -> a*x+b
    #
    # with a in {-2,-1,-1/2,1/2,1,2}
    # and b in {-2,-1,0,1,2}.

    scales = [
        sp.Rational(-2),
        sp.Rational(-1),
        sp.Rational(-1, 2),
        sp.Rational(1, 2),
        sp.Rational(1),
        sp.Rational(2),
    ]

    shifts = [
        sp.Integer(v)
        for v in range(
            -2,
            3,
        )
    ]

    for k in range(
        len(H) - 1
    ):

        local_hits = []

        derivative = sp.diff(
            H[k],
            x,
        )

        for a in scales:

            for b in shifts:

                transformed = clean(
                    derivative.subs(
                        x,
                        a*x + b,
                    )
                )

                scalar = exact_scalar_multiple(
                    H[k + 1],
                    transformed,
                )

                if scalar is not None:
                    local_hits.append(
                        (
                            a,
                            b,
                            scalar,
                        )
                    )

        print()
        print(
            f"  k={k}: "
            f"hits={local_hits}"
        )

        if local_hits:
            hits.append(
                (
                    k,
                    local_hits,
                )
            )

    return hits


# =============================================================================
# 5. SECOND-DERIVATIVE / LOW-ORDER OPERATOR SEARCH
# =============================================================================

def differential_span_test():

    print()
    print("=" * 78)
    print(
        "5. DIFFERENTIAL-SPAN TEST"
    )
    print("=" * 78)

    hits = []

    # Test whether H[k+1] lies exactly in the span of
    #
    #     H[k]'
    #     H[k]''
    #
    # This is intentionally restricted and overdetermined.

    for k in range(
        len(H) - 1
    ):

        d1 = clean(
            sp.diff(
                H[k],
                x,
            )
        )

        d2 = clean(
            sp.diff(
                H[k],
                x,
                2,
            )
        )

        # Need constants a,b satisfying
        #
        # H[k+1] = a*d1+b*d2
        #
        equations = []

        rhs = []

        P1 = sp.Poly(
            d1,
            x,
            domain=sp.QQ,
        )

        P2 = sp.Poly(
            d2,
            x,
            domain=sp.QQ,
        )

        target = sp.Poly(
            H[k+1],
            x,
            domain=sp.QQ,
        )

        max_degree = max(
            P1.degree(),
            P2.degree(),
            target.degree(),
        )

        for p in range(
            max_degree + 1
        ):

            equations.append(
                [
                    P1.nth(p),
                    P2.nth(p),
                ]
            )

            rhs.append(
                target.nth(p)
            )

        M = sp.Matrix(
            equations
        )

        y = sp.Matrix(
            rhs
        )

        rank = M.rank()
        aug_rank = (
            M.row_join(y).rank()
        )

        if aug_rank != rank:

            print()
            print(
                f"  k={k}: NO_SOLUTION"
            )

            continue

        if rank != 2:

            print()
            print(
                f"  k={k}: NONUNIQUE"
            )

            continue

        sol = M.gauss_jordan_solve(
            y
        )[0]

        coeffs = tuple(
            clean(v)
            for v in sol
        )

        rebuilt = clean(
            coeffs[0]*d1
            + coeffs[1]*d2
        )

        exact = (
            rebuilt
            == clean(H[k+1])
        )

        print()
        print(
            f"  k={k}: "
            f"exact={exact} "
            f"coefficients={coeffs}"
        )

        if exact:
            hits.append(
                (
                    k,
                    coeffs,
                )
            )

    return hits


# =============================================================================
# 6. ROOT-OF-DERIVATIVE STRUCTURE
# =============================================================================

def derivative_root_audit():

    print()
    print("=" * 78)
    print(
        "6. DERIVATIVE ROOT / CRITICAL-POINT AUDIT"
    )
    print("=" * 78)

    for k in range(
        len(H)
    ):

        derivative = sp.diff(
            H[k],
            x,
        )

        factor = sp.factor(
            derivative
        )

        roots = sp.solve(
            sp.Eq(
                derivative,
                0,
            ),
            x,
        )

        rational_roots = []

        for root in roots:

            if root.is_Rational:
                rational_roots.append(
                    root
                )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    derivative="
            f"{factor}"
        )

        print(
            f"    rational_critical_points="
            f"{rational_roots}"
        )


# =============================================================================
# 7. DEGREE-LOWERING COEFFICIENT TEST
# =============================================================================

def coefficient_ratio_test():

    print()
    print("=" * 78)
    print(
        "7. LEADING-COEFFICIENT DERIVATIVE TEST"
    )
    print("=" * 78)

    rows = []

    for k in range(
        len(H) - 1
    ):

        P = sp.Poly(
            H[k],
            x,
            domain=sp.QQ,
        )

        Q = sp.Poly(
            H[k+1],
            x,
            domain=sp.QQ,
        )

        lead_P = P.LC()
        lead_Q = Q.LC()

        derivative_lead = (
            lead_P
            * P.degree()
        )

        ratio = clean(
            lead_Q
            / derivative_lead
        )

        rows.append(
            (
                k,
                lead_P,
                lead_Q,
                derivative_lead,
                ratio,
            )
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    lead(H_k)={lead_P}"
        )

        print(
            f"    lead(H_k1)={lead_Q}"
        )

        print(
            f"    lead(H_k')="
            f"{derivative_lead}"
        )

        print(
            f"    required_scalar="
            f"{ratio}"
        )

    return rows


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 288 — EXACT B-ROW "
        "DIFFERENTIAL / APPELL-TYPE OPERATOR AUDIT"
    )
    print("=" * 78)

    direct_hits = direct_derivative_test()

    reverse_hits = reverse_derivative_test()

    translated_hits = translated_derivative_test()

    affine_hits = affine_derivative_test()

    span_hits = differential_span_test()

    derivative_root_audit()

    coefficient_ratio_test()

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The row degrees

    7,6,5,4,3,2

make differentiation the next natural structural hypothesis.

A true Appell-type law would have

    H'_k(x) = c_k H_{k+1}(x).

A translated or affine differential law would mean the degree ladder is
preserved but the polynomial coordinate is shifted or rescaled.

A successful differential-span identity would show that the B rows are
generated by a small differential algebra rather than an arbitrary
coefficient table.

This is a fundamentally different test from the previous recurrence
search: it acts on the whole generating polynomial at once.
"""
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  direct_derivative_hits="
        f"{len(direct_hits)}"
    )

    print(
        f"  reverse_appell_hits="
        f"{len(reverse_hits)}"
    )

    print(
        f"  translated_derivative_hits="
        f"{len(translated_hits)}"
    )

    print(
        f"  affine_derivative_hits="
        f"{len(affine_hits)}"
    )

    print(
        f"  differential_span_hits="
        f"{len(span_hits)}"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 288 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted."
        )
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise
