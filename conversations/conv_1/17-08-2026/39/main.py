#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 107
BIVARIATE BINOMIAL / CYCLOTOMIC COEFFICIENT LAW
m=(ell-1)/2, r=(k-1)/2
NEWTON-MOMENT COEFFICIENT TENSOR
LEAVE-ONE-ELL-OUT SYMBOLIC VALIDATION
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================

Goal
----
Experiment 106 found:

  1. exact Newton-moment coefficient decompositions;
  2. strong N-power structure;
  3. cyclotomic factors in the k-generating polynomials;
  4. very regular terminal coefficients.

The wrong next move would be another unconstrained polynomial fit.

This experiment instead asks whether the coefficient tensor can be
represented in a LOW-COMPLEXITY BINOMIAL BASIS in BOTH

    r = (k-1)/2
    m = (ell-1)/2

with coefficients polynomial in N.

We explicitly use leave-one-ell-out validation.

A candidate law must:

  * be exact over QQ[N];
  * use a bounded binomial/falling-factorial basis;
  * survive prediction of an unseen ell;
  * reproduce every supplied coefficient exactly.

This is a structural test, not an N-only factoring algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Dict, List, Tuple, Optional

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

N = sp.Symbol("N")
m_sym = sp.Symbol("m")
r_sym = sp.Symbol("r")
x = sp.Symbol("x")


# ============================================================================
# DATA
# ============================================================================
#
# Coefficient data from Experiment 106.
#
# Each detector is (k, ell), and the dictionary contains
# Newton coefficients c_j(N), j=1,...,ell-1.
#
# We include enough rows for leave-one-ell-out tests:
# ell = 5,7,9,11,13
#
# The data are entered symbolically, not numerically.
# ============================================================================

C = {
    (1, 3): {
        1: sp.Integer(1),
        2: sp.Integer(-1),
    },

    (1, 5): {
        1: 3*N + 1,
        2: 6*N - 1,
        3: 1,
        4: -1,
    },

    (3, 5): {
        1: 6*N**2,
        2: -N*(3*N - 4),
        3: 1 - 3*N,
        4: -1,
    },

    (1, 7): {
        1: -4*N**2 + 10*N + 1,
        2: -8*N**2 + 24*N - 1,
        3: 12*N + 1,
        4: 8*N - 1,
        5: 1,
        6: -1,
    },

    (3, 7): {
        1: 5*N**2*(N + 3),
        2: 2*N*(5*N**2 - 5*N + 3),
        3: (2*N - 1)*(3*N - 1),
        4: -(N - 1)*(3*N - 1),
        5: 1 - 3*N,
        6: -1,
    },

    (5, 7): {
        1: 15*N**4,
        2: -5*N**3*(N - 4),
        3: -5*N**2*(2*N - 3),
        4: -2*N*(5*N - 3),
        5: 1 - 5*N,
        6: -1,
    },

    (1, 9): {
        1: 5*N**3 - 20*N**2 + 21*N + 1,
        2: 10*N**3 - 45*N**2 + 62*N - 1,
        3: -15*N**2 + 65*N + 1,
        4: -10*N**2 + 60*N - 1,
        5: 25*N + 1,
        6: 10*N - 1,
        7: 1,
        8: -1,
    },

    (3, 9): {
        1: -N**2*(6*N**2 - 35*N - 28),
        2: -N*(12*N**3 - 76*N**2 + 21*N - 8),
        3: 18*N**3 + 15*N**2 - 7*N + 1,
        4: 12*N**3 - 10*N**2 + 6*N - 1,
        5: (2*N - 1)*(3*N - 1),
        6: -(N - 1)*(3*N - 1),
        7: 1 - 3*N,
        8: -1,
    },

    (5, 9): {
        1: 7*N**4*(N + 10),
        2: 7*N**3*(2*N**2 - 5*N + 8),
        3: N**2*(15*N**2 - 35*N + 28),
        4: -N*(5*N**3 - 20*N**2 + 21*N - 8),
        5: -(2*N - 1)*(5*N**2 - 5*N + 1),
        6: -10*N**2 + 6*N - 1,
        7: 1 - 5*N,
        8: -1,
    },

    (7, 9): {
        1: 28*N**6,
        2: -7*N**5*(N - 8),
        3: -7*N**4*(3*N - 10),
        4: -7*N**3*(5*N - 8),
        5: -7*N**2*(5*N - 4),
        6: -N*(21*N - 8),
        7: 1 - 7*N,
        8: -1,
    },

    (1, 11): {
        1: -6*N**4 + 35*N**3 - 56*N**2 + 36*N + 1,
        2: -12*N**4 + 76*N**3 - 147*N**2 + 128*N - 1,
        3: 18*N**3 - 111*N**2 + 203*N + 1,
        4: 12*N**3 - 94*N**2 + 258*N - 1,
        5: -30*N**2 + 205*N + 1,
        6: -12*N**2 + 124*N - 1,
        7: 42*N + 1,
        8: 12*N - 1,
        9: 1,
        10: -1,
    },

    (3, 11): {
        1: N**2*(7*N**3 - 56*N**2 + 126*N + 45),
        2: N*(14*N**4 - 119*N**3 + 308*N**2 - 36*N + 10),
        3: -21*N**4 + 175*N**3 + 28*N**2 - 9*N + 1,
        4: -14*N**4 + 140*N**3 - 21*N**2 + 8*N - 1,
        5: 35*N**3 + 15*N**2 - 7*N + 1,
        6: 14*N**3 - 10*N**2 + 6*N - 1,
        7: (2*N - 1)*(3*N - 1),
        8: -(N - 1)*(3*N - 1),
        9: 1 - 3*N,
        10: -1,
    },

    (5, 11): {
        1: -2*N**4*(4*N**2 - 42*N - 105),
        2: -2*N**3*(8*N**3 - 88*N**2 + 63*N - 60),
        3: N**2*(24*N**3 + 70*N**2 - 84*N + 45),
        4: N*(16*N**4 - 35*N**3 + 56*N**2 - 36*N + 10),
        5: (N - 1)*(3*N - 1)*(5*N**2 - 5*N + 1),
        6: -(N**2 - 3*N + 1)*(5*N**2 - 5*N + 1),
        7: -(2*N - 1)*(5*N**2 - 5*N + 1),
        8: -10*N**2 + 6*N - 1,
        9: 1 - 5*N,
        10: -1,
    },

    (7, 11): {
        1: 3*N**6*(3*N + 70),
        2: 6*N**5*(3*N**2 - 14*N + 42),
        3: 14*N**4*(2*N**2 - 9*N + 15),
        4: -N**3*(7*N**3 - 56*N**2 + 126*N - 120),
        5: -N**2*(21*N**3 - 70*N**2 + 84*N - 45),
        6: -N*(35*N**3 - 56*N**2 + 36*N - 10),
        7: -35*N**3 + 28*N**2 - 9*N + 1,
        8: -21*N**2 + 8*N - 1,
        9: 1 - 7*N,
        10: -1,
    },

    (9, 11): {
        1: 45*N**8,
        2: -3*N**7*(3*N - 40),
        3: -6*N**6*(6*N - 35),
        4: -84*N**5*(N - 3),
        5: -42*N**4*(3*N - 5),
        6: -6*N**3*(21*N - 20),
        7: -3*N**2*(28*N - 15),
        8: -2*N*(18*N - 5),
        9: 1 - 9*N,
        10: -1,
    },

    (1, 13): {
        1: 7*N**5 - 56*N**4 + 126*N**3 - 120*N**2 + 55*N + 1,
        2: 14*N**5 - 119*N**4 + 308*N**3 - 366*N**2 + 230*N - 1,
        3: -21*N**4 + 175*N**3 - 434*N**2 + 486*N + 1,
        4: -14*N**4 + 140*N**3 - 483*N**2 + 800*N - 1,
        5: 35*N**3 - 315*N**2 + 917*N + 1,
        6: 14*N**3 - 175*N**2 + 798*N - 1,
        7: -49*N**2 + 490*N + 1,
        8: -14*N**2 + 224*N - 1,
        9: 63*N + 1,
        10: 14*N - 1,
        11: 1,
        12: -1,
    },

    (3, 13): {
        1: -2*N**2*(4*N**4 - 42*N**3 + 126*N**2 - 165*N - 33),
        2: -N*(16*N**5 - 176*N**4 + 588*N**3 - 912*N**2 + 55*N - 12),
        3: 24*N**5 - 260*N**4 + 840*N**3 + 45*N**2 - 11*N + 1,
        4: 16*N**5 - 200*N**4 + 848*N**3 - 36*N**2 + 10*N - 1,
        5: -40*N**4 + 460*N**3 + 28*N**2 - 9*N + 1,
        6: -16*N**4 + 240*N**3 - 21*N**2 + 8*N - 1,
        7: 56*N**3 + 15*N**2 - 7*N + 1,
        8: 16*N**3 - 10*N**2 + 6*N - 1,
        9: (2*N - 1)*(3*N - 1),
        10: -(N - 1)*(3*N - 1),
        11: 1 - 3*N,
        12: -1,
    },

    (5, 13): {
        1: 3*N**4*(3*N**3 - 40*N**2 + 154*N + 165),
        2: N**3*(18*N**4 - 249*N**3 + 1044*N**2 - 330*N + 220),
        3: -3*N**2*(9*N**4 - 123*N**3 - 70*N**2 + 55*N - 22),
        4: -N*(18*N**5 - 276*N**4 + 126*N**3 - 120*N**2 + 55*N - 12),
        5: 45*N**5 + 70*N**4 - 84*N**3 + 45*N**2 - 11*N + 1,
        6: 18*N**5 - 35*N**4 + 56*N**3 - 36*N**2 + 10*N - 1,
        7: (N - 1)*(3*N - 1)*(5*N**2 - 5*N + 1),
        8: -(N**2 - 3*N + 1)*(5*N**2 - 5*N + 1),
        9: -(2*N - 1)*(5*N**2 - 5*N + 1),
        10: -10*N**2 + 6*N - 1,
        11: 1 - 5*N,
        12: -1,
    },

    (7, 13): {
        1: -N**6*(10*N**2 - 165*N - 924),
        2: -2*N**5*(10*N**3 - 170*N**2 + 231*N - 396),
        3: 3*N**4*(10*N**3 + 70*N**2 - 154*N + 165),
        4: 2*N**3*(10*N**4 - 42*N**3 + 126*N**2 - 165*N + 110),
        5: N**2*(28*N**4 - 126*N**3 + 210*N**2 - 165*N + 66),
        6: -N*(N - 1)*(7*N**4 - 49*N**3 + 77*N**2 - 43*N + 12),
        7: -(N - 1)*(3*N - 1)*(7*N**3 - 14*N**2 + 7*N - 1),
        8: -35*N**4 + 56*N**3 - 36*N**2 + 10*N - 1,
        9: -35*N**3 + 28*N**2 - 9*N + 1,
        10: -21*N**2 + 8*N - 1,
        11: 1 - 7*N,
        12: -1,
    },

    (9, 13): {
        1: 11*N**8*(N + 45),
        2: 11*N**7*(2*N**2 - 15*N + 72),
        3: 3*N**6*(15*N**2 - 110*N + 308),
        4: -3*N**5*(3*N**3 - 40*N**2 + 154*N - 264),
        5: -3*N**4*(12*N**3 - 70*N**2 + 154*N - 165),
        6: -2*N**3*(42*N**3 - 126*N**2 + 165*N - 110),
        7: -3*N**2*(42*N**3 - 70*N**2 + 55*N - 22),
        8: -N*(126*N**3 - 120*N**2 + 55*N - 12),
        9: -84*N**3 + 45*N**2 - 11*N + 1,
        10: -36*N**2 + 10*N - 1,
        11: 1 - 9*N,
        12: -1,
    },

    (11, 13): {
        1: 66*N**10,
        2: -11*N**9*(N - 20),
        3: -55*N**8*(N - 9),
        4: -33*N**7*(5*N - 24),
        5: -66*N**6*(5*N - 14),
        6: -66*N**5*(7*N - 12),
        7: -33*N**4*(14*N - 15),
        8: -110*N**3*(3*N - 2),
        9: -33*N**2*(5*N - 2),
        10: -N*(55*N - 12),
        11: 1 - 11*N,
        12: -1,
    },
}


# ============================================================================
# BASIC HELPERS
# ============================================================================

def falling(xv: sp.Expr, a: int) -> sp.Expr:
    """Falling factorial x(x-1)...(x-a+1)."""
    if a == 0:
        return sp.Integer(1)
    return sp.prod(xv - t for t in range(a))


def binom_poly(xv: sp.Expr, a: int) -> sp.Expr:
    """Polynomial binomial basis C(x,a)."""
    return sp.expand(falling(xv, a) / sp.factorial(a))


def m_of_ell(ell: int) -> int:
    assert ell % 2 == 1 and ell >= 3
    return (ell - 1) // 2


def r_of_k(k: int) -> int:
    assert k % 2 == 1 and k >= 1
    return (k - 1) // 2


def coeff_degree_in_N(expr: sp.Expr) -> int:
    expr = sp.Poly(sp.expand(expr), N)
    return expr.degree()


def leading_coeff_in_N(expr: sp.Expr) -> sp.Expr:
    poly = sp.Poly(sp.expand(expr), N)
    return poly.LC()


# ============================================================================
# NORMALIZATION
# ============================================================================
#
# For each coefficient c_{k,l,j}, use the highest observed N-degree
# on that (l,j) column, independent of k.
#
# Then subtract the terminal +/-1 part for j close to l-1.
#
# We DO NOT assume the normalization is correct. It is tested.
# ============================================================================

def column_degree(ell: int, j: int) -> int:
    rows = [(k, C[(k, ell)][j]) for (k, e) in C if e == ell and j in C[(k, e)]]
    return max(coeff_degree_in_N(v) for _, v in rows)


def normalize_coeff(ell: int, j: int, expr: sp.Expr) -> sp.Expr:
    d = column_degree(ell, j)
    if d <= 0:
        return sp.expand(expr)
    return sp.expand(expr / N**d)


# ============================================================================
# BINOMIAL-BASIS RECOGNITION
# ============================================================================
#
# For fixed ell and j, test whether the k dependence can be represented by
#
#   sum_{a=0}^A B_a(N) * C(r,a)
#
# with small A.
#
# This is not ordinary interpolation: we bound A and independently bound
# the degree of B_a(N).
# ============================================================================

def fit_binomial_in_r(rows: List[Tuple[int, sp.Expr]],
                      max_order: int,
                      max_n_degree: int):
    """
    Fit:
        f(r,N) = sum_a sum_b alpha[a,b] C(r,a) N^b

    exactly over QQ.

    Unknown count = (max_order+1)*(max_n_degree+1).

    We construct equations coefficient-wise in N.
    """
    alpha = {}

    unknowns = []
    for a in range(max_order + 1):
        for b in range(max_n_degree + 1):
            z = sp.Symbol(f"A_{a}_{b}")
            alpha[(a, b)] = z
            unknowns.append(z)

    equations = []

    for k, expr in rows:
        r = r_of_k(k)
        candidate = 0
        for a in range(max_order + 1):
            for b in range(max_n_degree + 1):
                candidate += alpha[(a, b)] * binom_poly(r, a) * N**b

        diff = sp.Poly(sp.expand(candidate - expr), N)

        for coeff in diff.all_coeffs():
            equations.append(sp.Eq(coeff, 0))

    if not equations:
        return None

    sol = sp.linsolve([eq.lhs for eq in equations], unknowns)

    if sol == sp.EmptySet:
        return None

    sols = list(sol)
    if len(sols) != 1:
        return None

    vals = sols[0]

    # Reject free parameters.
    if any(v.free_symbols.intersection(set(unknowns)) for v in vals):
        return None

    formula = 0
    for (a, b), v in zip(unknowns, vals):
        formula += v * binom_poly(r_sym, a) * N**b

    return sp.factor(formula)


# ============================================================================
# LEAVE-ONE-ELL-OUT
# ============================================================================

def fit_and_predict_ell(held_out_ell: int,
                        orders=(1, 2, 3, 4, 5),
                        n_degrees=(0, 1, 2, 3, 4, 5, 6)):
    training_ells = sorted({
        ell for (_, ell) in C
        if ell != held_out_ell
    })

    test_rows = [(k, expr)
                 for (k, ell) in C
                 if ell == held_out_ell
                 for expr in [C[(k, ell)][1]]]

    best = None

    for j in range(1, held_out_ell):
        # Find smallest exact binomial model among the training ells.
        train_rows = []
        for k, ell in C:
            if ell == held_out_ell:
                continue
            if j in C[(k, ell)]:
                train_rows.append((k, C[(k, ell)][j]))

        found = None

        for order in orders:
            for degN in n_degrees:
                # A single coefficient function is NOT enough across varying ell.
                # This section intentionally reports per-ell fits only.
                #
                # We therefore skip cross-ell guessing here and perform
                # consistency checks below.
                found = fit_binomial_in_r(
                    [(k, v) for k, v in train_rows],
                    max_order=order,
                    max_n_degree=degN,
                )
                if found is not None:
                    break
            if found is not None:
                break

        # Per-ell rows do not share all k across ell, so this test is only
        # diagnostic. The real cross-ell test is handled later.
        best = best or {}
        best[j] = found

    return best


# ============================================================================
# CROSS-ELL MODEL
# ============================================================================
#
# The real test:
#
# f(r,m,N) =
#   sum_{a=0}^A sum_{b=0}^B sum_{c=0}^D
#      gamma[a,b,c] C(r,a) C(m,b) N^c.
#
# This is deliberately bounded.
#
# We remove the most obvious terminal coefficients first:
# j = 2m and j = 2m-1 are checked separately.
# ============================================================================

def cross_ell_rows(j: int,
                   held_out_ell: Optional[int] = None):
    rows = []
    for (k, ell), coeffs in C.items():
        if held_out_ell is not None and ell == held_out_ell:
            continue
        if j in coeffs:
            rows.append((k, ell, coeffs[j]))
    return rows


def fit_bivariate(j: int,
                  max_r_order: int,
                  max_m_order: int,
                  max_n_degree: int,
                  held_out_ell: Optional[int] = None):
    rows = cross_ell_rows(j, held_out_ell=held_out_ell)

    unknowns = []
    basis = {}

    for a in range(max_r_order + 1):
        for b in range(max_m_order + 1):
            for c in range(max_n_degree + 1):
                z = sp.Symbol(f"G_{a}_{b}_{c}")
                unknowns.append(z)
                basis[(a, b, c)] = z

    equations = []

    for k, ell, expr in rows:
        rv = r_of_k(k)
        mv = m_of_ell(ell)

        pred = 0
        for a in range(max_r_order + 1):
            for b in range(max_m_order + 1):
                for c in range(max_n_degree + 1):
                    pred += (
                        basis[(a, b, c)]
                        * binom_poly(rv, a)
                        * binom_poly(mv, b)
                        * N**c
                    )

        poly = sp.Poly(sp.expand(pred - expr), N)

        for coeff in poly.all_coeffs():
            equations.append(coeff)

    if not equations:
        return None

    try:
        solution = sp.linsolve(equations, unknowns)
    except Exception:
        return None

    if solution == sp.EmptySet:
        return None

    tuples = list(solution)
    if len(tuples) != 1:
        return None

    vals = tuples[0]

    # Free parameters mean underdetermination, not evidence.
    unknown_set = set(unknowns)
    if any(
        bool(v.free_symbols.intersection(unknown_set))
        for v in vals
    ):
        return None

    formula = 0
    for (a, b, c), v in zip(basis, vals):
        formula += (
            v
            * binom_poly(r_sym, a)
            * binom_poly(m_sym, b)
            * N**c
        )

    return sp.factor(formula)


def evaluate_formula(formula: sp.Expr, k: int, ell: int) -> sp.Expr:
    return sp.expand(
        formula.subs({
            r_sym: r_of_k(k),
            m_sym: m_of_ell(ell),
        })
    )


# ============================================================================
# TERMINAL COEFFICIENT TEST
# ============================================================================

def terminal_law_test():
    print("\n6. TERMINAL COEFFICIENT LAW")
    print("-" * 78)

    failures = 0
    total = 0

    for (k, ell), coeffs in sorted(C.items(), key=lambda t: (t[0][1], t[0][0])):
        m = m_of_ell(ell)

        expected_top = sp.Integer(-1)
        expected_prev = 1 - ell*N

        for j, expected in [
            (2*m, expected_top),
            (2*m - 1, expected_prev),
        ]:
            total += 1
            got = sp.expand(coeffs[j])

            if sp.expand(got - expected) != 0:
                failures += 1
                print(
                    f"FAIL detector=({k},{ell}) j={j}: "
                    f"got={got}, expected={expected}"
                )

    print(f"terminal-law failures = {failures}/{total}")
    print("status =", "PASS" if failures == 0 else "FAIL")


# ============================================================================
# CYCLOTOMIC TERMINAL GENERATOR
# ============================================================================

def cyclotomic_terminal_test():
    print("\n7. CYCLOTOMIC TERMINAL GENERATOR")
    print("-" * 78)

    for ell in sorted({ell for _, ell in C}):
        m = m_of_ell(ell)
        poly = -x * sum(x**(2*t) for t in range(m))

        # Factorization over Z.
        fac = sp.factor(poly)

        print(f"ell={ell}")
        print(f"  terminal generator = {fac}")

        expected_terminal = -x * sum(x**(2*t) for t in range(m))

        if sp.expand(poly - expected_terminal) != 0:
            raise ArithmeticError("Terminal generator identity failed.")

    print("all terminal cyclotomic identities = PASS")


# ============================================================================
# LEAVE-ONE-ELL-OUT TEST
# ============================================================================

def leave_one_ell_out_test():
    print("\n8. LEAVE-ONE-ELL-OUT BIVARIATE TEST")
    print("-" * 78)

    ells = sorted({ell for _, ell in C})

    # Candidate model sizes.
    candidates = [
        (1, 1, 1),
        (2, 1, 1),
        (2, 2, 1),
        (3, 2, 2),
        (3, 3, 2),
        (4, 3, 3),
        (5, 4, 4),
    ]

    for held_out in ells:
        print(f"\nheld-out ell={held_out}")

        total = 0
        best_found = None

        # We fit each j independently, but force a common model class.
        for j in range(1, held_out):
            possible = None

            for a, b, c in candidates:
                formula = fit_bivariate(
                    j,
                    max_r_order=a,
                    max_m_order=b,
                    max_n_degree=c,
                    held_out_ell=held_out,
                )

                if formula is not None:
                    possible = (a, b, c, formula)
                    break

            if possible is None:
                continue

            a, b, c, formula = possible

            # Validate on every row of held-out ell that contains j.
            for k, coeffs in C.items():
                if coeffs[1] is not None:
                    pass

            for (k, ell), coeffs in C.items():
                if ell != held_out or j not in coeffs:
                    continue

                total += 1

                pred = evaluate_formula(formula, k, ell)
                actual = sp.expand(coeffs[j])

                if sp.expand(pred - actual) != 0:
                    print(
                        f"  FAIL j={j}, detector=({k},{ell})"
                    )
                    break
            else:
                if best_found is None:
                    best_found = (a, b, c)

        if best_found is None:
            print("  no bounded bivariate law found")
        else:
            print(
                "  smallest model class with at least one successful column:",
                best_found
            )


# ============================================================================
# STRONGER TEST:
# fit using ell <= 11, predict all ell=13.
# Then rotate.
# ============================================================================

def forward_holdout_test():
    print("\n9. FORWARD HOLDOUT: TRAIN ell<=11, TEST ell=13")
    print("-" * 78)

    train_ells = {5, 7, 9, 11}
    test_ell = 13

    total = 0
    failures = 0
    discovered = 0

    for j in range(1, test_ell):
        formula = None

        for a, b, c in [
            (1, 1, 1),
            (2, 1, 1),
            (2, 2, 1),
            (3, 2, 2),
            (3, 3, 2),
            (4, 3, 3),
            (5, 4, 4),
        ]:
            formula = fit_bivariate(
                j,
                a,
                b,
                c,
                held_out_ell=test_ell,
            )
            if formula is not None:
                discovered += 1
                break

        if formula is None:
            continue

        for (k, ell), coeffs in C.items():
            if ell != test_ell or j not in coeffs:
                continue

            total += 1
            pred = evaluate_formula(formula, k, ell)
            actual = sp.expand(coeffs[j])

            if sp.expand(pred - actual) != 0:
                failures += 1
                print(
                    f"  FAIL j={j}, k={k}, "
                    f"pred={pred}, actual={actual}"
                )

    print(f"discovered columns = {discovered}/{test_ell-1}")
    print(f"prediction failures = {failures}/{total}")

    if failures == 0 and total > 0:
        print("STATUS = POSITIVE CANDIDATE")
    else:
        print("STATUS = NO UNIVERSAL LAW YET")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("KAPPA EXPERIMENT 107")
    print("BIVARIATE BINOMIAL / CYCLOTOMIC COEFFICIENT LAW")
    print("m=(ell-1)/2, r=(k-1)/2")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    print("\n1. DATASET")
    ells = sorted({ell for _, ell in C})
    ks_by_ell = {
        ell: sorted(k for (k, e) in C if e == ell)
        for ell in ells
    }

    for ell in ells:
        print(
            f"ell={ell}: m={m_of_ell(ell)} "
            f"k={ks_by_ell[ell]}"
        )

    print("\n2. EXACT INPUT CONSISTENCY")
    failures = 0

    for (k, ell), coeffs in C.items():
        expected_len = ell - 1
        if sorted(coeffs.keys()) != list(range(1, expected_len + 1)):
            failures += 1
            print(
                f"bad coefficient support for ({k},{ell}): "
                f"{sorted(coeffs.keys())}"
            )

    print(f"support failures = {failures}")
    if failures:
        raise ArithmeticError("Malformed coefficient tensor.")

    print("\n3. TERMINAL STRUCTURE")
    terminal_law_test()

    print("\n4. DEGREE PROFILE")
    print("-" * 78)

    for ell in ells:
        print(f"ell={ell}")
        for j in range(1, ell):
            vals = [
                C[(k, ell)][j]
                for k in ks_by_ell[ell]
                if j in C[(k, ell)]
            ]
            d = max(coeff_degree_in_N(v) for v in vals)
            print(f"  j={j:2d}: max deg_N={d}")

    print("\n5. LEADING-TERM CONSISTENCY")
    print("-" * 78)

    for ell in ells:
        m = m_of_ell(ell)
        print(f"ell={ell}")

        for j in range(1, 2*m + 1):
            vals = [
                C[(k, ell)][j]
                for k in ks_by_ell[ell]
            ]
            d = max(coeff_degree_in_N(v) for v in vals)
            lead = sorted({
                sp.Poly(v, N).LC()
                for v in vals
                if sp.Poly(v, N).degree() == d
            }, key=str)

            print(
                f"  j={j:2d}: degree={d}, "
                f"leading_coeffs={lead}"
            )

    cyclotomic_terminal_test()

    leave_one_ell_out_test()
    forward_holdout_test()

    print("\n10. FINAL DIAGNOSTIC")
    print("=" * 78)
    print(
        """
Experiment 107 asks a narrower question than Experiment 106:

    Is the coefficient tensor generated by a bounded bivariate
    binomial/cyclotomic law in

        r=(k-1)/2
        m=(ell-1)/2

    rather than by unrelated polynomial interpolation?

The experiment is successful only if:

  * the model is exact over QQ[N];
  * the same bounded basis works across multiple ell;
  * an entire ell value is predicted exactly;
  * the terminal cyclotomic structure is reproduced.

A positive result would give a genuine closed-form coefficient law.

A negative result would mean the observed coefficient regularity is
real but needs a different basis.

Either result is useful because it tells us whether the next step
should be symbolic derivation or abandonment of this detector-tensor
route.

This remains a structural/oracle experiment.
It does not assume that detector values are computable from N alone.
"""
    )

    print("=" * 78)
    print("EXPERIMENT 107 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

