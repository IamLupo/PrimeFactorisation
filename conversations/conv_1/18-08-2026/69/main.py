#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 71 — FULL SYMBOLIC LAYER DEPENDENCY AUDIT
==============================================================================

Purpose
-------
Determine whether the 17 homogeneous layers of the exact G_{9,16}(N,X)
contain independent information, or whether all layers are generated from
the same dimensionless quantity

    t = N/X

with explicit homogeneous scaling

    L_d(N,X) = X^d h_d(t).

Everything is exact SymPy arithmetic.

NO floating point.
NO fitting.
NO interpolation.
NO regression.
NO statistics.
==============================================================================
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Iterable, Sequence

import sympy as sp


# ============================================================================
# GLOBAL SYMBOLS
# ============================================================================

N, X, t = sp.symbols("N X t")
p, q = sp.symbols("p q")

K = 9
ELL = 16


# ============================================================================
# EXACT G_{9,16}(N,X)
# ============================================================================

def exact_G_9_16(Nv: sp.Expr, Xv: sp.Expr) -> sp.Expr:
    Nl = Nv
    Xl = Xv

    return sp.expand(
        -50*Nl**12
        + 288*Nl**11*Xl**2
        - 300*Nl**11*Xl
        + 1300*Nl**11
        - 276*Nl**10*Xl**4
        + 1584*Nl**10*Xl**3
        - 7436*Nl**10*Xl**2
        + 7150*Nl**10*Xl
        - 10010*Nl**10
        + 88*Nl**9*Xl**6
        - 1380*Nl**9*Xl**5
        + 10340*Nl**9*Xl**4
        - 34430*Nl**9*Xl**3
        + 62634*Nl**9*Xl**2
        - 50050*Nl**9*Xl
        + 35750*Nl**9
        - 9*Nl**8*Xl**8
        + 396*Nl**8*Xl**7
        - 5460*Nl**8*Xl**6
        + 34650*Nl**8*Xl**5
        - 117810*Nl**8*Xl**4
        + 228228*Nl**8*Xl**3
        - 252252*Nl**8*Xl**2
        + 160875*Nl**8*Xl
        - 71500*Nl**8
        - 36*Nl**7*Xl**9
        + 1164*Nl**7*Xl**8
        - 13560*Nl**7*Xl**7
        + 79464*Nl**7*Xl**6
        - 267960*Nl**7*Xl**5
        + 552552*Nl**7*Xl**4
        - 708708*Nl**7*Xl**3
        + 563420*Nl**7*Xl**2
        - 286000*Nl**7*Xl
        + 88400*Nl**7
        - 84*Nl**6*Xl**10
        + 2226*Nl**6*Xl**9
        - 23058*Nl**6*Xl**8
        + 127512*Nl**6*Xl**7
        - 426888*Nl**6*Xl**6
        + 918918*Nl**6*Xl**5
        - 1303302*Nl**6*Xl**4
        + 1221220*Nl**6*Xl**3
        - 755664*Nl**6*Xl**2
        + 309400*Nl**6*Xl
        - 71400*Nl**6
        - 126*Nl**5*Xl**11
        + 2898*Nl**5*Xl**10
        - 27510*Nl**5*Xl**9
        + 145530*Nl**5*Xl**8
        - 484110*Nl**5*Xl**7
        + 1075074*Nl**5*Xl**6
        - 1639638*Nl**5*Xl**5
        + 1732458*Nl**5*Xl**4
        - 1265992*Nl**5*Xl**3
        + 636888*Nl**5*Xl**2
        - 214200*Nl**5*Xl
        + 38760*Nl**5
        - 126*Nl**4*Xl**12
        + 2604*Nl**4*Xl**11
        - 23100*Nl**4*Xl**10
        + 117975*Nl**4*Xl**9
        - 390225*Nl**4*Xl**8
        + 887172*Nl**4*Xl**7
        - 1429428*Nl**4*Xl**6
        + 1653470*Nl**4*Xl**5
        - 1375360*Nl**4*Xl**4
        + 818720*Nl**4*Xl**3
        - 344352*Nl**4*Xl**2
        + 96900*Nl**4*Xl
        - 14250*Nl**4
        - 84*Nl**3*Xl**13
        + 1596*Nl**3*Xl**12
        - 13410*Nl**3*Xl**11
        + 66550*Nl**3*Xl**10
        - 219010*Nl**3*Xl**9
        + 507078*Nl**3*Xl**8
        - 852852*Nl**3*Xl**7
        + 1058540*Nl**3*Xl**6
        - 974400*Nl**3*Xl**5
        + 663680*Nl**3*Xl**4
        - 331704*Nl**3*Xl**3
        + 119016*Nl**3*Xl**2
        - 28500*Nl**3*Xl
        + 3500*Nl**3
        - 36*Nl**2*Xl**14
        + 639*Nl**2*Xl**13
        - 5135*Nl**2*Xl**12
        + 24882*Nl**2*Xl**11
        - 81510*Nl**2*Xl**10
        + 191477*Nl**2*Xl**9
        - 333333*Nl**2*Xl**8
        + 437700*Nl**2*Xl**7
        - 436832*Nl**2*Xl**6
        + 331500*Nl**2*Xl**5
        - 190296*Nl**2*Xl**4
        + 81624*Nl**2*Xl**3
        - 25380*Nl**2*Xl**2
        + 5250*Nl**2*Xl
        - 550*Nl**2
        - 9*Nl*Xl**15
        + 151*Nl*Xl**14
        - 1169*Nl*Xl**13
        + 5551*Nl*Xl**12
        - 18109*Nl*Xl**11
        + 43043*Nl*Xl**10
        - 77077*Nl*Xl**9
        + 105979*Nl*Xl**8
        - 112964*Nl*Xl**7
        + 93604*Nl*Xl**6
        - 60144*Nl*Xl**5
        + 29736*Nl*Xl**4
        - 11130*Nl*Xl**3
        + 3038*Nl*Xl**2
        - 550*Nl*Xl
        + 50*Nl
        - Xl**16
        + 16*Xl**15
        - 120*Xl**14
        + 560*Xl**13
        - 1820*Xl**12
        + 4368*Xl**11
        - 8008*Xl**10
        + 11441*Xl**9
        - 12879*Xl**8
        + 11476*Xl**7
        - 8092*Xl**6
        + 4494*Xl**5
        - 1946*Xl**4
        + 644*Xl**3
        - 156*Xl**2
        + 25*Xl
        - 2
    )


# ============================================================================
# POLYNOMIAL CONSTRUCTION
# ============================================================================

def poly_NX(expr: sp.Expr) -> sp.Poly:
    """
    Exact bivariate polynomial in N and X.
    """
    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain="QQ",
    )


def poly_t(expr: sp.Expr) -> sp.Poly:
    """
    Exact univariate polynomial in t.
    """
    return sp.Poly(
        sp.expand(expr),
        t,
        domain="QQ",
    )


# ============================================================================
# CORRECT MULTIVARIATE DEGREE FUNCTIONS
# ============================================================================

def degree_N(expr: sp.Expr) -> int | sp.Integer:
    """
    Degree in N, while allowing X to remain a genuine second variable.
    """
    poly = poly_NX(expr)

    if poly.is_zero:
        return -sp.oo

    return int(
        max(
            monom[0]
            for monom, coeff in poly.terms()
            if coeff != 0
        )
    )


def degree_X(expr: sp.Expr) -> int | sp.Integer:
    """
    Degree in X, while allowing N to remain a genuine second variable.
    """
    poly = poly_NX(expr)

    if poly.is_zero:
        return -sp.oo

    return int(
        max(
            monom[1]
            for monom, coeff in poly.terms()
            if coeff != 0
        )
    )


def total_degree(expr: sp.Expr) -> int | sp.Integer:
    poly = poly_NX(expr)

    if poly.is_zero:
        return -sp.oo

    return int(
        max(
            monom[0] + monom[1]
            for monom, coeff in poly.terms()
            if coeff != 0
        )
    )


def t_degree(expr: sp.Expr) -> int | sp.Integer:
    poly = poly_t(expr)

    if poly.is_zero:
        return -sp.oo

    return int(poly.degree())


def term_count_NX(expr: sp.Expr) -> int:
    poly = poly_NX(expr)
    return len(poly.terms())


def term_count_t(expr: sp.Expr) -> int:
    poly = poly_t(expr)
    return len(poly.terms())


# ============================================================================
# HOMOGENEOUS LAYERS
# ============================================================================

def homogeneous_layers(expr: sp.Expr) -> dict[int, sp.Expr]:
    poly = poly_NX(expr)

    result = {
        d: sp.Integer(0)
        for d in range(ELL + 1)
    }

    for (n_exp, x_exp), coeff in poly.terms():
        d = n_exp + x_exp
        result[d] += (
            coeff
            * N**n_exp
            * X**x_exp
        )

    return {
        d: sp.expand(value)
        for d, value in result.items()
    }


# ============================================================================
# EXACT X-VALUATION
# ============================================================================

def x_valuation(expr: sp.Expr) -> int | sp.Integer:
    """
    Minimal exponent of X among monomials of expr.

    This is the correct polynomial notion of v_X for these layers.
    """
    poly = poly_NX(expr)

    if poly.is_zero:
        return sp.oo

    return int(
        min(
            monom[1]
            for monom, coeff in poly.terms()
            if coeff != 0
        )
    )


# ============================================================================
# EXACT N-VALUATION
# ============================================================================

def n_valuation(expr: sp.Expr) -> int | sp.Integer:
    poly = poly_NX(expr)

    if poly.is_zero:
        return sp.oo

    return int(
        min(
            monom[0]
            for monom, coeff in poly.terms()
            if coeff != 0
        )
    )


# ============================================================================
# RATIONAL CONTENT
# ============================================================================

def rational_content(
    expr: sp.Expr,
    gens: Sequence[sp.Symbol],
) -> sp.Rational:
    """
    Exact content for arbitrary polynomial variables.

    Examples:
        rational_content(expr, (N, X))
        rational_content(h, (t,))
    """
    expanded = sp.expand(expr)

    if expanded == 0:
        return sp.Integer(0)

    poly = sp.Poly(
        expanded,
        *gens,
        domain="QQ",
    )

    coeffs = poly.coeffs()

    if not coeffs:
        return sp.Integer(0)

    numerator_gcd = 0
    denominator_lcm = 1

    for c in coeffs:
        num = int(sp.numer(c))
        den = int(sp.denom(c))

        numerator_gcd = math.gcd(
            numerator_gcd,
            abs(num),
        )

        denominator_lcm = math.lcm(
            denominator_lcm,
            den,
        )

    return sp.Rational(
        numerator_gcd,
        denominator_lcm,
    )


# ============================================================================
# PRIMITIVE COEFFICIENTS
# ============================================================================

def primitive_integer_coefficients(
    expr: sp.Expr,
    gens: Sequence[sp.Symbol],
) -> list[int]:

    expanded = sp.expand(expr)

    if expanded == 0:
        return [0]

    poly = sp.Poly(
        expanded,
        *gens,
        domain="QQ",
    )

    coeffs = poly.coeffs()

    den_lcm = 1

    for c in coeffs:
        den_lcm = math.lcm(
            den_lcm,
            int(sp.denom(c)),
        )

    ints: list[int] = []

    for c in coeffs:
        scaled = (
            sp.Rational(c)
            * sp.Integer(den_lcm)
        )
        ints.append(int(scaled))

    content = 0

    for value in ints:
        content = math.gcd(
            content,
            abs(value),
        )

    if content == 0:
        return [0]

    ints = [
        value // content
        for value in ints
    ]

    for value in ints:
        if value != 0:
            if value < 0:
                ints = [-v for v in ints]
            break

    return ints


# ============================================================================
# NORMALIZED LAYER
# ============================================================================

@dataclass
class NormalizedLayer:
    d: int
    h: sp.Expr
    success: bool
    reconstruction: bool


def normalize_layer(
    layer: sp.Expr,
    d: int,
) -> NormalizedLayer:

    poly = poly_NX(layer)

    if poly.is_zero:
        return NormalizedLayer(
            d=d,
            h=sp.Integer(0),
            success=False,
            reconstruction=True,
        )

    h = sp.Integer(0)

    for (n_exp, x_exp), coeff in poly.terms():

        if n_exp + x_exp != d:
            raise ValueError(
                f"Layer d={d} is not homogeneous: "
                f"monomial {(n_exp, x_exp)}"
            )

        # Since
        #
        #   X^d * (N/X)^n = N^n X^(d-n)
        #
        # the coefficient is copied directly to t^n.
        h += coeff * t**n_exp

    h = sp.expand(h)

    reconstructed = sp.expand(
        X**d * h.subs(
            t,
            N / X,
        )
    )

    return NormalizedLayer(
        d=d,
        h=h,
        success=True,
        reconstruction=(
            sp.expand(reconstructed - layer) == 0
        ),
    )


# ============================================================================
# GCD
# ============================================================================

def gcd_polynomials(
    expressions: Iterable[sp.Expr],
    gens: Sequence[sp.Symbol],
) -> sp.Expr:

    items = [
        sp.expand(x)
        for x in expressions
        if sp.expand(x) != 0
    ]

    if not items:
        return sp.Integer(0)

    g = sp.Poly(
        items[0],
        *gens,
        domain="QQ",
    )

    for expr in items[1:]:
        g = sp.gcd(
            g,
            sp.Poly(
                expr,
                *gens,
                domain="QQ",
            ),
        )

    return sp.factor(g.as_expr())


# ============================================================================
# POLYNOMIAL RANK
# ============================================================================

def coefficient_row_t(expr: sp.Expr) -> list[sp.Rational]:
    poly = poly_t(expr)

    if poly.is_zero:
        return [sp.Integer(0)]

    return [
        sp.Rational(c)
        for c in poly.all_coeffs()
    ]


def coefficient_matrix_t(
    expressions: Sequence[sp.Expr],
) -> sp.Matrix:

    if not expressions:
        return sp.Matrix([])

    max_deg = max(
        t_degree(expr)
        for expr in expressions
    )

    rows = []

    for expr in expressions:

        coeffs = coefficient_row_t(expr)

        padding = (
            max_deg + 1 - len(coeffs)
        )

        rows.append(
            [sp.Integer(0)] * padding
            + coeffs
        )

    return sp.Matrix(rows)


def polynomial_rank_t(
    expressions: Sequence[sp.Expr],
) -> int:

    if not expressions:
        return 0

    M = coefficient_matrix_t(
        expressions
    )

    return int(M.rank())


# ============================================================================
# RANDOM PRIME
# ============================================================================

def random_prime(
    rng: random.Random,
    low: int,
    high: int,
) -> int:

    while True:
        candidate = rng.randrange(
            low | 1,
            high + 1,
            2,
        )

        if sp.isprime(candidate):
            return candidate


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("EXPERIMENT 71 — FULL SYMBOLIC LAYER DEPENDENCY AUDIT")
    print("=" * 78)
    print()

    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)
    print(
        f"  k = {K}\n"
        f"  ell = {ELL}\n"
        "  kernel = self-contained exact G_{9,16}(N,X)\n"
        "  floating point = forbidden"
    )
    print()

    # ----------------------------------------------------------------------
    # Build kernel and layers.
    # ----------------------------------------------------------------------

    G = exact_G_9_16(N, X)
    layers = homogeneous_layers(G)

    # ----------------------------------------------------------------------
    # 1
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("1. COMPLETE HOMOGENEOUS LAYER PROFILE")
    print("=" * 78)

    for d in range(ELL, -1, -1):

        Ld = layers[d]

        print(
            f"  degree={d:2d}: "
            f"present={Ld != 0} "
            f"deg_N={degree_N(Ld)} "
            f"deg_X={degree_X(Ld)} "
            f"terms={term_count_NX(Ld)}"
        )

    nonzero_count = sum(
        1
        for d in range(ELL + 1)
        if layers[d] != 0
    )

    print(
        f"  nonzero layers = "
        f"{nonzero_count} / {ELL + 1}"
    )

    print()

    # ----------------------------------------------------------------------
    # 2
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT X-DIVISIBILITY PROFILE")
    print("=" * 78)
    print()
    print(
        "  d | v_X(L_d) | deg_N | deg_X | normalized_terms"
    )
    print(
        "  -------------------------------------------------"
    )

    normalized: dict[int, NormalizedLayer] = {}

    for d in range(ELL, -1, -1):

        Ld = layers[d]
        info = normalize_layer(
            Ld,
            d,
        )
        normalized[d] = info

        print(
            f" {d:2d} |"
            f" {str(x_valuation(Ld)):>8} |"
            f" {str(degree_N(Ld)):>6} |"
            f" {str(degree_X(Ld)):>6} |"
            f" {term_count_t(info.h):>16}"
        )

    print()

    # ----------------------------------------------------------------------
    # 3
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("3. EXACT NORMALIZED REPRESENTATION")
    print("=" * 78)
    print()
    print("Testing")
    print()
    print("    L_d(N,X) = X^d h_d(N/X)")
    print()

    normalization_ok = True

    for d in range(ELL, -1, -1):

        info = normalized[d]

        if not info.success:
            print(
                f"  degree={d}: zero layer"
            )
            continue

        normalization_ok &= info.reconstruction

        print(
            f"  degree={d}: "
            f"alpha={d} "
            f"degree(h)={t_degree(info.h)} "
            f"terms(h)={term_count_t(info.h)} "
            f"reconstruction={info.reconstruction}"
        )

    print()
    print(
        f"  ALL NORMALIZATIONS EXACT = "
        f"{normalization_ok}"
    )
    print()

    # ----------------------------------------------------------------------
    # 4
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("4. ALL NORMALIZED LAYER POLYNOMIALS h_d(t)")
    print("=" * 78)

    for d in range(ELL, -1, -1):

        h = normalized[d].h

        if h == 0:
            continue

        print()
        print(
            f"  h_{d}(t) ="
        )
        print(
            f"    {sp.factor(h)}"
        )

    print()

    # ----------------------------------------------------------------------
    # 5
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("5. GCD STRUCTURE OF NORMALIZED LAYER POLYNOMIALS")
    print("=" * 78)

    h_list = [
        normalized[d].h
        for d in range(ELL, -1, -1)
        if normalized[d].h != 0
    ]

    gcd_all = gcd_polynomials(
        h_list,
        (t,),
    )

    print(
        f"  gcd(h_16,...,h_0) = "
        f"{sp.factor(gcd_all)}"
    )

    print(
        f"  gcd degree = "
        f"{t_degree(gcd_all)}"
    )

    print()
    print("  ADJACENT GCDS")

    for d in range(ELL, 0, -1):

        h1 = normalized[d].h
        h2 = normalized[d - 1].h

        g = gcd_polynomials(
            [h1, h2],
            (t,),
        )

        print(
            f"    gcd(h_{d},h_{d-1}) = "
            f"{sp.factor(g)}"
        )

    print()

    # ----------------------------------------------------------------------
    # 6
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("6. SCALE STRUCTURE")
    print("=" * 78)
    print()

    print(
        "  d | X-power | alpha=d | content(L_d) | content(h_d)"
    )
    print(
        "  ------------------------------------------------------"
    )

    for d in range(ELL, -1, -1):

        Ld = layers[d]
        info = normalized[d]

        content_L = rational_content(
            Ld,
            (N, X),
        )

        content_h = rational_content(
            info.h,
            (t,),
        )

        print(
            f" {d:2d} |"
            f" {str(x_valuation(Ld)):>7} |"
            f" {str(d):>8} |"
            f" {str(content_L):>14} |"
            f" {str(content_h):>14}"
        )

    print()

    # ----------------------------------------------------------------------
    # 7
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("7. PRIMITIVE INTEGER SIGNATURES")
    print("=" * 78)

    for d in range(ELL, -1, -1):

        h = normalized[d].h

        print(
            f"  h_{d}: "
            f"{primitive_integer_coefficients(h, (t,))}"
        )

    print()

    # ----------------------------------------------------------------------
    # 8
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("8. CUMULATIVE NORMALIZED-LAYER RANK")
    print("=" * 78)

    cumulative: list[sp.Expr] = []

    for d in range(ELL, -1, -1):

        cumulative.append(
            normalized[d].h
        )

        rank = polynomial_rank_t(
            cumulative
        )

        print(
            f"  layers {ELL}..{d:2d}: "
            f"count={len(cumulative):2d} "
            f"rank={rank:2d}"
        )

    print()

    # ----------------------------------------------------------------------
    # 9
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("9. TOP-THREE GENERATION TEST")
    print("=" * 78)

    top_three = [
        normalized[16].h,
        normalized[15].h,
        normalized[14].h,
    ]

    top_rank = polynomial_rank_t(
        top_three
    )

    print(
        f"  top-three rank = {top_rank}"
    )
    print()

    for d in range(13, -1, -1):

        hd = normalized[d].h

        before = polynomial_rank_t(
            top_three
        )

        after = polynomial_rank_t(
            top_three + [hd]
        )

        print(
            f"  h_{d:2d}: "
            f"rank_before={before} "
            f"rank_after={after} "
            f"new_independent_direction="
            f"{after > before}"
        )

    print()

    # ----------------------------------------------------------------------
    # 10
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("10. ADJACENT NORMALIZED-LAYER RATIO TEST")
    print("=" * 78)

    for d in range(ELL, 0, -1):

        upper = normalized[d].h
        lower = normalized[d - 1].h

        ratio = sp.factor(
            sp.cancel(
                lower / upper
            )
        )

        numerator, denominator = sp.fraction(
            ratio
        )

        print()
        print(
            f"  h_{d-1}/h_{d} = {ratio}"
        )
        print(
            f"    numerator_degree="
            f"{t_degree(numerator)}"
        )
        print(
            f"    denominator_degree="
            f"{t_degree(denominator)}"
        )

    print()

    # ----------------------------------------------------------------------
    # 11
    # Fresh semiprime check
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("11. FRESH SEMIPRIME VERIFICATION")
    print("=" * 78)

    rng = random.Random(20260818)

    pv = random_prime(
        rng,
        50000,
        300000,
    )

    qv = random_prime(
        rng,
        50000,
        300000,
    )

    while qv == pv:
        qv = random_prime(
            rng,
            50000,
            300000,
        )

    Nv = sp.Integer(pv * qv)
    Sv = sp.Integer(pv + qv)
    Xv = Sv + 1

    print(
        f"  p = {pv}\n"
        f"  q = {qv}\n"
        f"  N = {Nv}\n"
        f"  S = {Sv}\n"
        f"  X = {Xv}"
    )
    print()

    fresh_G = exact_G_9_16(
        Nv,
        Xv,
    )

    substituted_G = exact_G_9_16(
        sp.Integer(pv * qv),
        sp.Integer(pv + qv + 1),
    )

    print(
        "  G(N,X) == G(pq,p+q+1): "
        f"{sp.expand(fresh_G - substituted_G) == 0}"
    )

    # ----------------------------------------------------------------------
    # 12
    # Fresh layer checks
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. FRESH LAYER CHECKS")
    print("=" * 78)

    fresh_layer_values: dict[int, sp.Expr] = {}

    for d in range(ELL, -1, -1):

        value = sp.expand(
            layers[d].subs(
                {
                    N: Nv,
                    X: Xv,
                }
            )
        )

        fresh_layer_values[d] = value

        print(
            f"  L_{d:2d} = {value}"
        )

    layer_sum = sp.expand(
        sum(
            fresh_layer_values.values()
        )
    )

    print()
    print(
        "  sum(L_d) == G(N,X): "
        f"{sp.expand(layer_sum - fresh_G) == 0}"
    )

    print()

    # ----------------------------------------------------------------------
    # 13
    # Fresh normalized h_d reconstruction
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("13. FRESH NORMALIZED-LAYER RECONSTRUCTION")
    print("=" * 78)

    fresh_normalization_ok = True

    for d in range(ELL, -1, -1):

        h = normalized[d].h

        predicted = sp.expand(
            Xv**d
            * h.subs(
                t,
                Nv / Xv,
            )
        )

        observed = fresh_layer_values[d]

        ok = (
            sp.expand(predicted - observed)
            == 0
        )

        fresh_normalization_ok &= ok

        print(
            f"  degree={d:2d}: "
            f"predicted==observed={ok}"
        )

    print()
    print(
        f"  ALL FRESH NORMALIZATIONS EXACT = "
        f"{fresh_normalization_ok}"
    )
    print()

    # ----------------------------------------------------------------------
    # 14
    # Symbolic p,q substitution
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("14. SYMBOLIC p,q -> N,X SUBSTITUTION")
    print("=" * 78)

    symbolic_check = sp.expand(
        exact_G_9_16(
            p * q,
            p + q + 1,
        )
        -
        G.subs(
            {
                N: p * q,
                X: p + q + 1,
            }
        )
    )

    print(
        "  exact symbolic substitution = "
        f"{symbolic_check == 0}"
    )
    print()

    # ----------------------------------------------------------------------
    # 15
    # Key interpretation
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("15. FINAL STRUCTURAL SUMMARY")
    print("=" * 78)
    print()

    print(
        "  EXACT RESULT:"
    )
    print(
        "    Every homogeneous layer satisfies"
    )
    print(
        "        L_d(N,X) = X^d h_d(N/X)"
    )
    print()

    print(
        "  Therefore the entire kernel can be written as"
    )
    print(
        "        G(N,X) = sum_d X^d h_d(N/X)"
    )
    print()

    print(
        "  This does NOT mean that every lower layer is redundant."
    )
    print(
        "  The h_d(t) themselves may be algebraically independent."
    )
    print()

    print(
        f"  gcd(all h_d) = {sp.factor(gcd_all)}"
    )
    print()

    print(
        "  The critical unresolved question is:"
    )
    print()
    print(
        "    Are h_16,h_15,h_14 sufficient to determine t=N/X,"
    )
    print(
        "    or do lower h_d contain genuinely new invariants?"
    )
    print()

    print(
        "  If lower layers increase the exact QQ-rank, they contain"
    )
    print(
        "  additional information about the same hidden parameter t."
    )
    print()

    print(
        "  That would make them candidates for a stronger inverse"
    )
    print(
        "  reconstruction using more than three observed layers."
    )
    print()

    print(
        "  No factorization algorithm is inferred."
    )
    print(
        "  No fitting, interpolation, regression, or statistics."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 71 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()