#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 149
EXACT DETECTOR RECURRENCE IN ELL
N-ONLY RECURRENCE COEFFICIENT SEARCH

PURPOSE
-------
Experiments 147R and 148 ruled out low-degree S-independent
linear/quadratic detector invariants.

We now test a different mechanism:

    Can the detector sequence itself satisfy a recurrence in ell

        Q_(1,ell+r)
          + A_(r-1)(N) Q_(1,ell+r-1)
          + ...
          + A_0(N) Q_(1,ell)
        = 0

with coefficients depending ONLY on N?

If YES, then detector values might be generated recursively from N
without first knowing S.

This is the first experiment aimed directly at:

        N -> detector sequence -> hidden trace S.

We use k=1 initially and odd ell only:

    ell = 3,5,7,9,...

Therefore the natural sequence index is

    m = (ell-3)/2.

The recurrence is searched in this m-index.

SEARCH:
    recurrence order r = 1,...,4
    degree in N <= 3

For each candidate order/degree:

    Q_(m+r)
      + A_(r-1)(N) Q_(m+r-1)
      + ...
      + A_0(N) Q_m
      = 0

is imposed as an EXACT polynomial identity in S.

The leading coefficient is fixed to 1, so there is no scaling
ambiguity.

Then every discovered recurrence is tested on future ell values
not used during fitting.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# Exact finite quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q_pq(k, ell):
    F = kernel_F(k, ell)

    poly = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    divisor = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(
        poly,
        divisor
    )

    rem = sp.expand(
        R.as_expr()
    )

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(rem)}"
        )

    return sp.expand(
        Q.as_expr()
    )


# ============================================================================
# Symmetric reduction
# ============================================================================

def symmetric_to_NS(expr):
    """
    q = S-p

    reduce modulo

        p^2 - S p + N.

    The result must be independent of p.
    """

    substituted = sp.expand(
        expr.subs(
            q,
            S - p
        )
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        substituted,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.expand(
        sp.rem(
            poly,
            modulus
        ).as_expr()
    )

    p_part = sp.expand(
        rem.coeff(p, 1)
    )

    if p_part != 0:
        raise ArithmeticError(
            "symmetric reduction retained p-dependence"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


def detector(k, ell):
    return symmetric_to_NS(
        quotient_Q_pq(k, ell)
    )


# ============================================================================
# Cache detector sequence
# ============================================================================

def build_cache(max_ell):
    cache = {}

    print(
        "building detector cache..."
    )

    for ell in range(
        3,
        max_ell + 1,
        2
    ):
        print(
            f"  ell={ell}"
        )

        cache[ell] = sp.expand(
            detector(
                1,
                ell
            )
        )

    return cache


# ============================================================================
# Coefficient vector in S
# ============================================================================

def coeff_in_S(expr, power):
    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if power > poly.degree():
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            S**power
        )
    )


# ============================================================================
# Build exact recurrence system
# ============================================================================

def build_recurrence_matrix(
    cache,
    order,
    degree_N,
    training_ell
):
    """
    Search:

        Q_(m+r) + sum_{i=0}^{r-1} A_i(N) Q_(m+i) = 0

    where each

        A_i(N) = sum_{d=0}^{degree_N} a_(i,d) N^d.

    The coefficient of the highest recurrence term is fixed to 1.
    """

    block = degree_N + 1
    n_unknowns = order * block

    rows = []

    for ell0 in training_ell:

        target_ell = ell0 + 2*order

        if target_ell not in cache:
            continue

        # Sequence terms:
        #
        # ell0, ell0+2, ..., ell0+2r

        terms = [
            cache[
                ell0 + 2*i
            ]
            for i in range(
                order + 1
            )
        ]

        max_S_degree = max(
            sp.Poly(
                term,
                S,
                domain=sp.QQ.frac_field(N)
            ).degree()
            for term in terms
        )

        for s_power in range(
            max_S_degree + 1
        ):

            # Build coefficient of S^s as a polynomial in N.
            #
            # Highest term Q_(m+r) contributes with coefficient 1.
            # Unknowns multiply Q_(m+i), i=0,...,r-1.

            known = coeff_in_S(
                terms[order],
                s_power
            )

            basis_polys = []

            for i in range(order):

                q_coeff = coeff_in_S(
                    terms[i],
                    s_power
                )

                for dN in range(
                    degree_N + 1
                ):

                    basis_polys.append(
                        sp.expand(
                            N**dN * q_coeff
                        )
                    )

            # The equation is:
            #
            # known + sum c_j*basis_j = 0
            #
            # Produce one row for every N-power.

            max_N_degree = 0

            for expr in basis_polys:

                poly_N = sp.Poly(
                    sp.expand(expr),
                    N,
                    domain=sp.QQ,
                )

                if not poly_N.is_zero:
                    max_N_degree = max(
                        max_N_degree,
                        poly_N.degree()
                    )

            known_poly = sp.Poly(
                sp.expand(known),
                N,
                domain=sp.QQ,
            )

            if not known_poly.is_zero:
                max_N_degree = max(
                    max_N_degree,
                    known_poly.degree()
                )

            for n_power in range(
                max_N_degree + 1
            ):

                row = []

                for expr in basis_polys:

                    poly_N = sp.Poly(
                        sp.expand(expr),
                        N,
                        domain=sp.QQ,
                    )

                    row.append(
                        sp.Rational(
                            poly_N.coeff_monomial(
                                N**n_power
                            )
                        )
                    )

                known_value = sp.Rational(
                    known_poly.coeff_monomial(
                        N**n_power
                    )
                )

                # Equation:
                #
                # row*c + known_value = 0.
                #
                # Therefore matrix*c = -known_value.

                if (
                    any(
                        value != 0
                        for value in row
                    )
                    or known_value != 0
                ):

                    rows.append(
                        (
                            row,
                            -known_value
                        )
                    )

    if not rows:
        return (
            sp.zeros(0, n_unknowns),
            sp.zeros(0, 1)
        )

    A = sp.Matrix(
        [
            row for row, _ in rows
        ]
    )

    b = sp.Matrix(
        [
            rhs for _, rhs in rows
        ]
    )

    return A, b


# ============================================================================
# Solve exact recurrence system
# ============================================================================

def solve_recurrence(
    cache,
    order,
    degree_N,
    training_ell
):
    A, b = build_recurrence_matrix(
        cache,
        order,
        degree_N,
        training_ell
    )

    print(
        f"      matrix = {A.rows} x {A.cols}"
    )

    if A.cols == 0:
        return []

    solution_set = sp.linsolve(
        (
            A,
            b
        )
    )

    if solution_set == sp.EmptySet:
        return []

    solutions = list(
        solution_set
    )

    if not solutions:
        return []

    solution = solutions[0]

    # Reject parametric solutions for this experiment.
    free_symbols = set()

    for value in solution:
        free_symbols |= value.free_symbols

    free_symbols -= {
        N,
        S
    }

    if free_symbols:
        return [
            {
                "parametric": True,
                "solution": solution
            }
        ]

    return [
        {
            "parametric": False,
            "solution": solution
        }
    ]


# ============================================================================
# Convert recurrence coefficients to polynomials
# ============================================================================

def solution_to_coefficients(
    solution,
    order,
    degree_N
):
    coeffs = list(
        solution["solution"]
    )

    block = degree_N + 1

    polys = []

    for i in range(order):

        A_i = sp.expand(
            sum(
                coeffs[
                    i*block + dN
                ] * N**dN
                for dN in range(
                    degree_N + 1
                )
            )
        )

        polys.append(
            sp.factor(A_i)
        )

    return polys


# ============================================================================
# Verify recurrence exactly
# ============================================================================

def verify_recurrence(
    cache,
    order,
    coeffs,
    test_ell
):
    """
    Check:

        Q_(ell+2r)
        + sum_i A_i(N) Q_(ell+2i)
        = 0

    exactly as a polynomial in N,S.
    """

    failures = []

    for ell0 in test_ell:

        target_ell = ell0 + 2*order

        if target_ell not in cache:
            continue

        expr = cache[target_ell]

        for i, A_i in enumerate(
            coeffs
        ):
            expr += (
                A_i
                * cache[
                    ell0 + 2*i
                ]
            )

        expr = sp.expand(
            expr
        )

        if expr != 0:

            failures.append(
                (
                    ell0,
                    sp.factor(expr)
                )
            )

    return failures


# ============================================================================
# Pretty print recurrence
# ============================================================================

def print_recurrence(
    order,
    coeffs
):
    print()

    print(
        "      recurrence:"
    )

    text = "      Q_(m+%d)" % order

    for i, A_i in enumerate(
        coeffs
    ):

        text += (
            f" + ({sp.factor(A_i)})"
            f"*Q_(m+{i})"
        )

    text += " = 0"

    print(
        text
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 149")
    print("EXACT DETECTOR RECURRENCE IN ELL")
    print("N-ONLY RECURRENCE COEFFICIENT SEARCH")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 15
    FORWARD_MAX_ELL = 21

    cache = build_cache(
        FORWARD_MAX_ELL
    )

    print()

    print("=" * 78)
    print("1. DETECTOR SEQUENCE")
    print("=" * 78)

    for ell in range(
        3,
        12,
        2
    ):
        print(
            f"Q_(1,{ell}) = "
            f"{sp.factor(cache[ell])}"
        )

    print()

    # Training windows.
    #
    # Each starting ell consumes 2*order additional ell-units.
    #
    training_starts = [
        3,
        5,
        7,
    ]

    forward_starts = [
        9,
        11,
        13,
    ]

    MAX_ORDER = 4
    MAX_DEGREE_N = 3

    winners = []

    print("=" * 78)
    print("2. EXACT RECURRENCE SEARCH")
    print("=" * 78)

    for order in range(
        1,
        MAX_ORDER + 1
    ):

        for degree_N in range(
            MAX_DEGREE_N + 1
        ):

            print()
            print(
                f"order={order}, "
                f"degree_N<={degree_N}"
            )

            solutions = solve_recurrence(
                cache,
                order,
                degree_N,
                training_starts
            )

            print(
                f"      solutions = "
                f"{len(solutions)}"
            )

            for solution in solutions:

                if solution["parametric"]:

                    print(
                        "      parametric solution: "
                        "SKIPPED"
                    )

                    continue

                coeffs = solution_to_coefficients(
                    solution,
                    order,
                    degree_N
                )

                print_recurrence(
                    order,
                    coeffs
                )

                train_failures = verify_recurrence(
                    cache,
                    order,
                    coeffs,
                    training_starts
                )

                forward_failures = verify_recurrence(
                    cache,
                    order,
                    coeffs,
                    forward_starts
                )

                print(
                    "      training failures =",
                    len(train_failures)
                )

                print(
                    "      forward failures =",
                    len(forward_failures)
                )

                if (
                    not train_failures
                    and not forward_failures
                ):

                    winners.append(
                        {
                            "order": order,
                            "degree_N": degree_N,
                            "coeffs": coeffs,
                        }
                    )

    # ------------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FINAL DIAGNOSTIC")
    print("=" * 78)

    if winners:

        best = min(
            winners,
            key=lambda item: (
                item["order"],
                item["degree_N"]
            )
        )

        print(
            "STATUS = PASS"
        )

        print()

        print(
            "A genuine N-only recurrence for the detector"
        )

        print(
            "sequence was found."
        )

        print()

        print(
            "best order =",
            best["order"]
        )

        print(
            "N-degree =",
            best["degree_N"]
        )

        print_recurrence(
            best["order"],
            best["coeffs"]
        )

        print()

        print(
            "This is the first direct candidate for generating"
        )

        print(
            "detector values from N without explicitly knowing S."
        )

        print()

        print(
            "NEXT TARGET:"
        )

        print(
            "derive the recurrence symbolically from the finite"
        )

        print(
            "binomial kernel and determine its initial values."
        )

    else:

        print(
            "STATUS = NO-RECURRENCE"
        )

        print()

        print(
            "No N-only detector recurrence was found through"
        )

        print(
            f"order <= {MAX_ORDER}, "
            f"degree_N <= {MAX_DEGREE_N}."
        )

        print()

        print(
            "The recurrence/generating-function route is therefore"
        )

        print(
            "not supported at this tested complexity."
        )

        print()

        print(
            "NEXT TARGET:"
        )

        print(
            "test a recurrence whose coefficients depend on ell,"
        )

        print(
            "or move to a generating-function / holonomic analysis."
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

