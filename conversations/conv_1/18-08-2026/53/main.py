#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
R=5 EXPERIMENT 53
=================

Exact cross-k / cross-D boundary factor experiment.

The supplied discrepancies are treated as exact input data.

No old-law subtraction is performed.
No r=6 analysis.
No full pq-kernel expansion.
No replacement universal law.

Main goals
----------
1. Represent E_j(k,D) exactly.
2. Recover fixed-D cross-k polynomials P_{j,D}(K).
3. Compute exact gcds between nonzero D-slices.
4. Extract primitive integer forms of those gcds robustly.
5. Track factor persistence and primitive parts across D.
6. Reconstruct and verify all supplied discrepancies.
7. Keep j=5 terminal behaviour separate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

import sympy as sp


# ============================================================================
# SYMBOLS / EXACT DOMAINS
# ============================================================================

D = sp.Symbol("D")
K = sp.Symbol("K")

KD = sp.QQ.frac_field(K)
POLY_D_DOMAIN = sp.QQ


# ============================================================================
# EXPERIMENT GRID
# ============================================================================

K_VALUES: List[int] = [3, 5, 7, 9, 11, 13]
J_VALUES: List[int] = [0, 1, 2, 3, 4, 5]
D_VALUES: List[int] = [6, 8, 10, 12, 14, 16]

EXPECTED_POINTS = len(K_VALUES) * len(J_VALUES) * len(D_VALUES)


# ============================================================================
# EXACT INPUT POLYNOMIALS E_j(k,D)
#
# These are the exact D-polynomials established in the previous experiment.
# Each entry is a SymPy expression in D.
# ============================================================================

E_POLYS: Dict[Tuple[int, int], sp.Expr] = {
    # ----------------------------------------------------------------------
    # j = 0
    # ----------------------------------------------------------------------
    (3, 0):
        -(7363*D**5
          - 395215*D**4
          + 8629700*D**3
          - 91891700*D**2
          + 472005792*D
          - 932561280) / 960,

    (5, 0):
        (76393*D**5
         - 4305640*D**4
         + 90582740*D**3
         - 916615760*D**2
         + 4483939392*D
         - 8489617920) / 1920,

    (7, 0):
        -(15541*D**5
          - 632044*D**4
          + 11202092*D**3
          - 98824400*D**2
          + 426557568*D
          - 720147456) / 384,

    (9, 0):
        -(615341*D**5
          - 28009020*D**4
          + 533637820*D**3
          - 5012129040*D**2
          + 22981059264*D
          - 41124902400) / 3840,

    (11, 0):
        -13*(9677*D**5
             - 485334*D**4
             + 9685172*D**3
             - 94225800*D**2
             + 445818272*D
             - 820598016) / 128,

    (13, 0):
        -119*(25623*D**5
              - 1350100*D**4
              + 27793940*D**3
              - 276874160*D**2
              + 1335166912*D
              - 2495823360) / 960,

    # ----------------------------------------------------------------------
    # j = 1
    # ----------------------------------------------------------------------
    (3, 1):
        -(D - 6) * (
            149*D**4
            - 2683*D**3
            - 2402*D**2
            + 398140*D
            - 2180640
        ) / 96,

    (5, 1):
        -7*(D - 6) * (
            181*D**4
            - 8406*D**3
            + 150268*D**2
            - 1168344*D
            + 3311680
        ) / 48,

    (7, 1):
        (D - 6) * (
            1679*D**4
            - 96396*D**3
            + 1753904*D**2
            - 13562106*D
            + 38170160
        ) / 80,

    (9, 1):
        -11*(D - 6) * (
            11371*D**4
            - 508124*D**3
            + 8504336*D**2
            - 62038984*D
            + 166242240
        ) / 480,

    (11, 1):
        -91*(D - 6) * (
            2863*D**4
            - 136182*D**3
            + 2382328*D**2
            - 17966172*D
            + 49397920
        ) / 480,

    (13, 1):
        -17*(D - 6) * (
            531*D**4
            - 79544*D**3
            + 2118616*D**2
            - 19978264*D
            + 63103040
        ) / 120,

    # ----------------------------------------------------------------------
    # j = 2
    # ----------------------------------------------------------------------
    (3, 2):
        (D - 6) * (
            22707*D**4
            - 1110703*D**3
            + 19642662*D**2
            - 149850488*D
            + 416996160
        ) / 960,

    (5, 2):
        -(D - 6) * (
            2312*D**4
            - 103197*D**3
            + 1707416*D**2
            - 12332148*D
            + 32788640
        ) / 48,

    (7, 2):
        -(D - 6) * (
            15669*D**4
            - 732214*D**3
            + 12707648*D**2
            - 95638016*D
            + 263119360
        ) / 384,

    (9, 2):
        -(D - 6) * (
            21077*D**4
            - 1071324*D**3
            + 19742740*D**2
            - 155018688*D
            + 439799040
        ) / 384,

    (11, 2):
        7*(D - 6) * (
            4955*D**4
            - 182498*D**3
            + 2446376*D**2
            - 14262136*D
            + 30560960
        ) / 192,

    (13, 2):
        17*(D - 6) * (
            94421*D**4
            - 3956744*D**3
            + 61568116*D**2
            - 421903984*D
            + 1073639040
        ) / 960,

    # ----------------------------------------------------------------------
    # j = 3
    # ----------------------------------------------------------------------
    (3, 3):
        -7*(D - 8)*(D - 6) * (
            131*D**3
            - 5546*D**2
            + 77468*D
            - 350330
        ) / 80,

    (5, 3):
        (D - 8)*(D - 6) * (
            4051*D**3
            - 157178*D**2
            + 1973512*D
            - 8082000
        ) / 128,

    (7, 3):
        -11*(D - 8)*(D - 6) * (
            2463*D**3
            - 93638*D**2
            + 1173424*D
            - 4819640
        ) / 960,

    (9, 3):
        13*(D - 8)*(D - 6) * (
            595*D**3
            - 20198*D**2
            + 226612*D
            - 842940
        ) / 96,

    (11, 3):
        (D - 8)*(D - 6) * (
            3358*D**3
            - 117757*D**2
            + 1366672*D
            - 5251345
        ) / 8,

    (13, 3):
        17*(D - 8)*(D - 6) * (
            163617*D**3
            - 5798122*D**2
            + 67971936*D
            - 263566160
        ) / 1920,

    # ----------------------------------------------------------------------
    # j = 4
    # ----------------------------------------------------------------------
    (3, 4):
        -(D - 16)*(D - 8)*(D - 6) * (
            1892*D**2
            - 40889*D
            + 219490
        ) / 64,

    (5, 4):
        (D - 8)*(D - 6) * (
            8315*D**3
            - 272506*D**2
            + 2955352*D
            - 10636160
        ) / 768,

    (7, 4):
        (D - 8)*(D - 6) * (
            124993*D**3
            - 4453618*D**2
            + 52402784*D
            - 203769440
        ) / 3840,

    (9, 4):
        (D - 8)*(D - 6) * (
            5649*D**3
            - 201758*D**2
            + 2380888*D
            - 9284240
        ) / 64,

    (11, 4):
        (D - 8)*(D - 6) * (
            487679*D**3
            - 17467064*D**2
            + 206706812*D
            - 808115120
        ) / 1920,

    (13, 4):
        (D - 8)*(D - 6) * (
            2350941*D**3
            - 84327646*D**2
            + 999345048*D
            - 3911794880
        ) / 3840,

    # ----------------------------------------------------------------------
    # j = 5
    # ----------------------------------------------------------------------
    (3, 5):
        (D - 12)*(D - 10)*(D - 8)*(D - 6) * (
            22277*D - 371378
        ) / 960,

    (5, 5):
        -1295*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 96,

    (7, 5):
        7007*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 480,

    (9, 5):
        1377*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 40,

    (11, 5):
        46189*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 640,

    (13, 5):
        -7*(D - 12)*(D - 10)*(D - 8)*(D - 6) * (
            37989*D - 626891
        ) / 480,
}


# ============================================================================
# DATA ACCESS
# ============================================================================

def discrepancy(k: int, j: int, d: int) -> sp.Rational:
    """Return the exact supplied discrepancy E_j(k,D)."""
    if (k, j) not in E_POLYS:
        raise KeyError(f"Missing polynomial for k={k}, j={j}")

    value = sp.factor(
        sp.expand(E_POLYS[(k, j)].subs(D, d))
    )

    if not value.is_Rational:
        value = sp.cancel(value)

    return sp.Rational(value)


def make_D_poly(k: int, j: int) -> sp.Poly:
    """Return E_j(k,D) as an exact polynomial over QQ."""
    return sp.Poly(
        sp.expand(E_POLYS[(k, j)]),
        D,
        domain=sp.QQ,
    )


# ============================================================================
# BASIC UTILITIES
# ============================================================================

def clean_expr(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.expand(expr)))


def poly_zero(poly: sp.Poly) -> bool:
    return poly.is_zero


def polynomial_from_k_values(
    values: Sequence[sp.Expr],
    k_values: Sequence[int],
) -> sp.Poly:
    """
    Exact interpolation in K through the supplied k-grid.
    """
    if len(values) != len(k_values):
        raise ValueError("values and k_values must have the same length")

    expr = sp.interpolate(
        [
            (sp.Integer(k), sp.Rational(v))
            for k, v in zip(k_values, values)
        ],
        K,
    )

    return sp.Poly(
        sp.expand(expr),
        K,
        domain=sp.QQ,
    )


def fixed_D_cross_k_poly(j: int, d: int) -> sp.Poly:
    """
    Exact fixed-D cross-k interpolation polynomial.
    Six k-values are available.
    """
    values = [
        discrepancy(k, j, d)
        for k in K_VALUES
    ]
    return polynomial_from_k_values(values, K_VALUES)


def normalized_poly(poly: sp.Poly) -> sp.Poly:
    """
    Monic normalization only for gcd comparison.

    Zero remains zero.
    """
    if poly.is_zero:
        return poly

    lc = poly.LC()
    return sp.Poly(
        sp.expand(poly.as_expr() / lc),
        K,
        domain=sp.QQ,
    )


def gcd_poly(p: sp.Poly, q: sp.Poly) -> sp.Poly:
    """
    Exact gcd over QQ[K], normalized to positive leading coefficient.
    """
    if p.is_zero and q.is_zero:
        return sp.Poly(0, K, domain=sp.QQ)

    if p.is_zero:
        return normalized_poly(q)

    if q.is_zero:
        return normalized_poly(p)

    g = sp.gcd(p, q)

    if g.is_zero:
        return g

    return normalized_poly(g)


def polynomial_degree(poly: sp.Poly) -> int:
    if poly.is_zero:
        return -sp.oo
    return int(poly.degree())


# ============================================================================
# ROBUST PRIMITIVE INTEGER FORM
# ============================================================================

def primitive_integer_form(
    poly: sp.Poly,
) -> Tuple[sp.Poly, Fraction]:
    """
    Convert a rational polynomial into

        poly = scale * primitive_integer_polynomial

    where the primitive polynomial has integer coefficients whose gcd is 1.

    Robust for:
      * zero polynomial
      * nonzero constant polynomial
      * one nonzero coefficient
      * arbitrary rational coefficients

    The primitive polynomial is normalized to have positive leading
    coefficient.  The returned rational scale preserves the original sign.
    """

    if poly.is_zero:
        return (
            sp.Poly(0, K, domain=sp.ZZ),
            Fraction(0, 1),
        )

    coeffs = [
        sp.Rational(c)
        for c in poly.all_coeffs()
    ]

    # ------------------------------------------------------------
    # Clear rational denominators.
    # ------------------------------------------------------------
    den_lcm = 1
    for c in coeffs:
        den = int(sp.denom(c))
        den_lcm = math.lcm(den_lcm, den)

    int_coeffs = [
        int(c * den_lcm)
        for c in coeffs
    ]

    nonzero_abs = [
        abs(v)
        for v in int_coeffs
        if v != 0
    ]

    if not nonzero_abs:
        return (
            sp.Poly(0, K, domain=sp.ZZ),
            Fraction(0, 1),
        )

    # ------------------------------------------------------------
    # Integer content.
    #
    # Never call sp.igcd(*list) directly because a one-element list
    # is rejected by the SymPy version used in the environment.
    # ------------------------------------------------------------
    content = nonzero_abs[0]

    for value in nonzero_abs[1:]:
        content = math.gcd(content, value)

    if content <= 0:
        raise ArithmeticError(
            f"Invalid integer content: {content}"
        )

    # ------------------------------------------------------------
    # Divide content.
    # ------------------------------------------------------------
    primitive_coeffs = [
        value // content
        for value in int_coeffs
    ]

    # ------------------------------------------------------------
    # Normalize primitive leading coefficient to positive.
    # ------------------------------------------------------------
    sign = 1

    if primitive_coeffs[0] < 0:
        primitive_coeffs = [
            -value
            for value in primitive_coeffs
        ]
        sign = -1

    primitive = sp.Poly.from_list(
        primitive_coeffs,
        gens=K,
        domain=sp.ZZ,
    )

    scale = Fraction(
        sign * content,
        den_lcm,
    )

    return primitive, scale


def primitive_integer_expression(poly: sp.Poly) -> Tuple[sp.Expr, Fraction]:
    """
    Convenience wrapper returning the primitive polynomial as an expression.
    """
    prim, scale = primitive_integer_form(poly)
    return prim.as_expr(), scale


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class Transition:
    j: int
    d1: int
    d2: int
    p1: sp.Poly
    p2: sp.Poly
    gcd: sp.Poly
    primitive: sp.Poly
    scale: Fraction

    @property
    def gcd_degree(self) -> int:
        return polynomial_degree(self.gcd)

    @property
    def primitive_degree(self) -> int:
        return polynomial_degree(self.primitive)


# ============================================================================
# SECTION 0
# ============================================================================

def section_0_data_validation() -> None:
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    print(f"k values = {K_VALUES}")
    print(f"j values = {J_VALUES}")
    print(f"D values = {D_VALUES}")
    print(f"points   = {EXPECTED_POINTS}")

    expected_keys = {
        (k, j)
        for k in K_VALUES
        for j in J_VALUES
    }

    actual_keys = set(E_POLYS.keys())

    missing = sorted(expected_keys - actual_keys)
    extra = sorted(actual_keys - expected_keys)

    if missing:
        raise RuntimeError(f"missing E-polynomials: {missing}")

    if extra:
        raise RuntimeError(f"extra E-polynomials: {extra}")

    for j in J_VALUES:
        for d in D_VALUES:
            for k in K_VALUES:
                discrepancy(k, j, d)

    print("grid status = OK")
    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_fixed_D_factorization() -> None:
    print("=" * 78)
    print("1. EXACT FIXED-D CROSS-k FACTORIZATION")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        for d in D_VALUES:
            p = fixed_D_cross_k_poly(j, d)

            if p.is_zero:
                print(f"  D={d}: ZERO POLYNOMIAL")
                continue

            print(f"  D={d}: degree={p.degree()}")
            print(f"       P(K)={sp.factor(p.as_expr())}")

            try:
                _, roots = sp.factor_list(p.as_expr())
                root_list = []

                for factor, multiplicity in roots:
                    if sp.degree(factor, K) == 1:
                        root_list.append(
                            sp.solve(sp.Eq(factor, 0), K)[0]
                        )

                print(f"       roots={root_list}")

            except Exception:
                print("       roots=<factorization unavailable>")

        print()


# ============================================================================
# SECTION 2
# ============================================================================

def build_neighbor_transitions() -> List[Transition]:
    transitions: List[Transition] = []

    for j in J_VALUES:
        for d1, d2 in zip(D_VALUES[:-1], D_VALUES[1:]):
            p1 = fixed_D_cross_k_poly(j, d1)
            p2 = fixed_D_cross_k_poly(j, d2)

            # Skip if either side is zero.
            if p1.is_zero or p2.is_zero:
                continue

            g = gcd_poly(p1, p2)

            if g.is_zero:
                primitive = sp.Poly(0, K, domain=sp.ZZ)
                scale = Fraction(0, 1)
            else:
                primitive, scale = primitive_integer_form(g)

            transitions.append(
                Transition(
                    j=j,
                    d1=d1,
                    d2=d2,
                    p1=p1,
                    p2=p2,
                    gcd=g,
                    primitive=primitive,
                    scale=scale,
                )
            )

    return transitions


def section_2_neighbor_gcds(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("2. GCD OF NEIGHBORING NONZERO D-SLICES")
    print("=" * 78)

    by_j: Dict[int, List[Transition]] = {
        j: []
        for j in J_VALUES
    }

    for t in transitions:
        by_j[t.j].append(t)

    for j in J_VALUES:
        print(f"j={j}")

        for t in by_j[j]:
            print(
                f"  D=({t.d1},{t.d2}): "
                f"degree={t.gcd_degree} "
                f"gcd={sp.factor(t.gcd.as_expr())}"
            )

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_primitive_parts(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("3. PRIMITIVE INTEGER PARTS OF NEIGHBOR GCDs")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        ts = [
            t for t in transitions
            if t.j == j
        ]

        if not ts:
            print("  no nonzero neighboring transitions")
            print()
            continue

        for t in ts:
            print(
                f"  D=({t.d1},{t.d2})"
            )
            print(
                f"    gcd={sp.factor(t.gcd.as_expr())}"
            )
            print(
                f"    primitive={sp.factor(t.primitive.as_expr())}"
            )
            print(
                f"    scale={t.scale}"
            )

        print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_primitive_gcd_across_transitions(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("4. GCD OF PRIMITIVE PARTS ACROSS D-TRANSITIONS")
    print("=" * 78)

    for j in J_VALUES:
        ts = [
            t for t in transitions
            if t.j == j
            and not t.primitive.is_zero
        ]

        print(f"j={j}")

        if not ts:
            print("  no primitive transition parts")
            print()
            continue

        current = ts[0].primitive

        for t in ts[1:]:
            current = gcd_poly(
                current,
                t.primitive,
            )

        if current.is_zero:
            print("  common primitive gcd=0")
            print()
            continue

        primitive_common, common_scale = primitive_integer_form(
            current
        )

        print(
            f"  common primitive gcd="
            f"{sp.factor(primitive_common.as_expr())}"
        )
        print(
            f"  scale={common_scale}"
        )

        for t in ts:
            g = gcd_poly(
                t.primitive,
                primitive_common,
            )

            print(
                f"    transition ({t.d1},{t.d2}): "
                f"quotient="
                f"{sp.factor((t.primitive.exquo(g)).as_expr())}"
            )

        print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_transition_factor_patterns(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("5. TRANSITION FACTOR PATTERNS")
    print("=" * 78)

    for j in J_VALUES:
        print(f"j={j}")

        ts = [
            t for t in transitions
            if t.j == j
            and not t.gcd.is_zero
        ]

        if not ts:
            print("  none")
            print()
            continue

        for t in ts:
            g = t.gcd

            factor_list = sp.factor_list(g.as_expr())[1]

            factors = []

            for factor, mult in factor_list:
                factors.append(
                    f"({sp.factor(factor)})^{mult}"
                )

            print(
                f"  ({t.d1},{t.d2}) "
                f"degree={t.gcd_degree}"
            )

            print(
                "    factors="
                + (
                    " * ".join(factors)
                    if factors
                    else "1"
                )
            )

        print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_transition_survival_matrix(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("6. TRANSITION SURVIVAL MATRIX")
    print("=" * 78)

    candidate_factors: Dict[int, List[sp.Expr]] = {
        0: [
            K + 1,
            K + 2,
            K + 3,
            K + 4,
            2*K + 5,
        ],
        1: [
            K + 2,
            K + 3,
            K + 4,
            2*K + 7,
        ],
        2: [
            K + 3,
            K + 4,
            2*K + 7,
        ],
        3: [
            K + 4,
            2*K + 9,
        ],
        4: [
            2*K + 9,
        ],
        5: [
            K - 5,
            K - 7,
            K - 9,
            K - 11,
        ],
    }

    for j in J_VALUES:
        print(f"j={j}")

        ts = {
            (t.d1, t.d2): t
            for t in transitions
            if t.j == j
        }

        if not candidate_factors[j]:
            print("  no candidate factors")
            print()
            continue

        labels = [
            f"{d1}->{d2}"
            for d1, d2 in zip(
                D_VALUES[:-1],
                D_VALUES[1:],
            )
        ]

        print("  D-transition:", " ".join(f"{x:>10}" for x in labels))

        for factor in candidate_factors[j]:
            row = []

            for d1, d2 in zip(
                D_VALUES[:-1],
                D_VALUES[1:],
            ):
                t = ts.get((d1, d2))

                if t is None or t.gcd.is_zero:
                    row.append("ZERO")
                    continue

                remainder = sp.rem(
                    t.gcd,
                    sp.Poly(factor, K, domain=sp.QQ),
                )

                survives = remainder.is_zero

                row.append(
                    "YES" if survives else "NO"
                )

            print(
                f"  factor={sp.factor(factor)}"
            )
            print(
                "    "
                + " ".join(
                    f"{x:>10}"
                    for x in row
                )
            )

        print()


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_terminal_j5() -> None:
    print("=" * 78)
    print("7. TERMINAL j=5 TRANSITION AUDIT")
    print("=" * 78)

    ladder = (
        (K - 5)
        * (K - 7)
        * (K - 9)
        * (K - 11)
    )

    print(
        "terminal ladder="
        f"{sp.factor(ladder)}"
    )

    for d in D_VALUES:
        p = fixed_D_cross_k_poly(5, d)

        if p.is_zero:
            print(f"D={d}: ZERO")
            continue

        remainder = sp.rem(
            p,
            sp.Poly(ladder, K, domain=sp.QQ),
        )

        survives = remainder.is_zero

        print(
            f"D={d}: "
            f"ladder_survives={survives}"
        )

        if survives:
            quotient = p.exquo(
                sp.Poly(ladder, K, domain=sp.QQ)
            )
            print(
                f"  quotient={sp.factor(quotient.as_expr())}"
            )

    print()


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_onset_check() -> None:
    print("=" * 78)
    print("8. EXACT ONSET RECONSTRUCTION")
    print("=" * 78)

    # Established onset structure for j=0,...,4.
    onset_D = {
        0: 6,
        1: 8,
        2: 8,
        3: 10,
        4: 10,
    }

    onset_constants = {
        0: sp.Rational(-1, 20),
        1: sp.Rational(-2, 3),
        2: sp.Rational(-4, 3),
        3: sp.Rational(-5, 1),
        4: sp.Rational(-10, 1),
    }

    for j in range(5):
        d0 = onset_D[j]
        c = onset_constants[j]

        ladder = sp.prod(
            K + m
            for m in range(j + 1, 5)
        )

        candidate = sp.factor(
            c
            * ladder
            * (2*K + d0 - 1)
        )

        actual = fixed_D_cross_k_poly(
            j,
            d0,
        ).as_expr()

        actual = sp.factor(actual)

        print(f"j={j}")
        print(f"  D0={d0}")
        print(f"  candidate={candidate}")
        print(f"  actual={actual}")
        print(
            "  exact match="
            f"{sp.simplify(candidate - actual) == 0}"
        )
        print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_reconstruction_check() -> None:
    print("=" * 78)
    print("9. EXACT INPUT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for j in J_VALUES:
        for k in K_VALUES:
            poly = make_D_poly(k, j)

            for d in D_VALUES:
                tested += 1

                predicted = sp.factor(
                    poly.as_expr().subs(D, d)
                )
                supplied = discrepancy(k, j, d)

                if sp.simplify(
                    predicted - supplied
                ) != 0:
                    failures.append(
                        (k, j, d, predicted, supplied)
                    )

    print(f"tested={tested}")
    print(
        f"reconstruction failures={len(failures)}"
    )

    if failures:
        for failure in failures[:10]:
            print(
                "  failure="
                f"{failure}"
            )
        raise RuntimeError(
            "Exact input reconstruction failed."
        )

    print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_support_summary() -> None:
    print("=" * 78)
    print("10. SUPPORT / ONSET STRUCTURAL SUMMARY")
    print("=" * 78)

    support = {
        0: sp.Integer(1),
        1: D - 6,
        2: D - 6,
        3: (D - 8)*(D - 6),
        4: (D - 8)*(D - 6),
        5: (D - 12)*(D - 10)*(D - 8)*(D - 6),
    }

    onset_D = {
        0: 6,
        1: 8,
        2: 8,
        3: 10,
        4: 10,
        5: 14,
    }

    for j in J_VALUES:
        s = sp.factor(
            support[j]
        )

        print(
            f"j={j}: "
            f"support={s} "
            f"support_degree={sp.degree(s, D) if s != 1 else 0} "
            f"D0={onset_D[j]}"
        )

    print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_diagnostics(
    transitions: Sequence[Transition],
) -> None:
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The experiment uses the exact supplied "
        "R=5 discrepancy polynomials."
    )
    print()

    print(
        "The primary cross-D object is:"
    )
    print(
        "    gcd(P_{j,D1}(K), P_{j,D2}(K))"
    )
    print(
        "for neighboring nonzero D-slices."
    )
    print()

    print(
        "Primitive integer normalization is applied "
        "after exact polynomial gcd computation."
    )
    print()

    print(
        "The primitive normalization is robust for:"
    )
    print(
        "  * zero polynomials;"
    )
    print(
        "  * nonzero constants;"
    )
    print(
        "  * single nonzero coefficients;"
    )
    print(
        "  * general rational coefficient lists."
    )
    print()

    print(
        f"neighboring nonzero transitions={len(transitions)}"
    )

    degree_drops = []

    for t in transitions:
        if t.gcd_degree < 5:
            degree_drops.append(
                (
                    t.j,
                    t.d1,
                    t.d2,
                    t.gcd_degree,
                    sp.factor(t.gcd.as_expr()),
                )
            )

    print(
        f"gcd degree-drop transitions={len(degree_drops)}"
    )

    for item in degree_drops:
        print(
            f"  j={item[0]}, "
            f"D=({item[1]},{item[2]}), "
            f"degree={item[3]}, "
            f"gcd={item[4]}"
        )

    print()

    print("Important previously observed onset structure:")
    print(
        "  j=0: -(K+1)(K+2)(K+3)(K+4)(2K+5)/20"
    )
    print(
        "  j=1: -2(K+2)(K+3)(K+4)(2K+7)/3"
    )
    print(
        "  j=2: -4(K+3)(K+4)(2K+7)/3"
    )
    print(
        "  j=3: -5(K+4)(2K+9)"
    )
    print(
        "  j=4: -10(2K+9)"
    )
    print()

    print(
        "j=5 remains terminal and is not folded into "
        "the j=0,...,4 onset template."
    )
    print()

    print("No r=6.")
    print("No full pq-kernel expansion.")
    print("No replacement universal r,j law.")
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    section_0_data_validation()

    section_1_fixed_D_factorization()

    transitions = build_neighbor_transitions()

    section_2_neighbor_gcds(
        transitions
    )

    section_3_primitive_parts(
        transitions
    )

    section_4_primitive_gcd_across_transitions(
        transitions
    )

    section_5_transition_factor_patterns(
        transitions
    )

    section_6_transition_survival_matrix(
        transitions
    )

    section_7_terminal_j5()

    section_8_onset_check()

    section_9_reconstruction_check()

    section_10_support_summary()

    section_11_diagnostics(
        transitions
    )


if __name__ == "__main__":
    main()