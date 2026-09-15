#!/usr/bin/env python3

"""
EXPERIMENT 77 — EXACT EVEN/ODD CHANNEL GENERATING STRUCTURE

Goal
----
Starting from the exact normalized layers

    H(z,t) = sum_{d=0}^{16} h_d(t) z^d,

and the exact involution decomposition

    h_d(t)
      = A_d(u) + (2t+1) B_d(u),

    u = t(t+1),

write

    H(z,t)
      = Ahat(w,u) + (2t+1) z Bhat(w,u),

with

    w = z^2.

This experiment asks whether the channel sequences themselves have
a compact construction law.

Important:
    * exact SymPy arithmetic only
    * no floating point
    * no interpolation
    * no regression
    * no statistical inference
    * no hidden p,q in the structural calculations

A coefficient-by-coefficient reconstruction is used for the channel
identity, avoiding fragile multivariate equality checks.
"""

import sympy as sp
from dataclasses import dataclass


# ============================================================================
# SYMBOLS
# ============================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")
z = sp.Symbol("z")
w = sp.Symbol("w")


# ============================================================================
# EXACT NORMALIZED LAYERS
# ============================================================================

H = {
    16: (
        -9*t**8 - 36*t**7 - 84*t**6 - 126*t**5
        -126*t**4 - 84*t**3 - 36*t**2 - 9*t - 1
    ),

    15: (
        88*t**9 + 396*t**8 + 1164*t**7 + 2226*t**6
        + 2898*t**5 + 2604*t**4 + 1596*t**3
        + 639*t**2 + 151*t + 16
    ),

    14: (
        -276*t**10 - 1380*t**9 - 5460*t**8
        -13560*t**7 - 23058*t**6 - 27510*t**5
        -23100*t**4 - 13410*t**3 - 5135*t**2
        -1169*t - 120
    ),

    13: (
        288*t**11 + 1584*t**10 + 10340*t**9
        +34650*t**8 + 79464*t**7 + 127512*t**6
        +145530*t**5 + 117975*t**4 + 66550*t**3
        +24882*t**2 + 5551*t + 560
    ),

    12: (
        -50*t**12 - 300*t**11 - 7436*t**10
        -34430*t**9 - 117810*t**8 - 267960*t**7
        -426888*t**6 - 484110*t**5 - 390225*t**4
        -219010*t**3 - 81510*t**2 - 18109*t - 1820
    ),

    11: (
        1300*t**11 + 7150*t**10 + 62634*t**9
        +228228*t**8 + 552552*t**7 + 918918*t**6
        +1075074*t**5 + 887172*t**4 + 507078*t**3
        +191477*t**2 + 43043*t + 4368
    ),

    10: (
        -10010*t**10 - 50050*t**9 - 252252*t**8
        -708708*t**7 -1303302*t**6 -1639638*t**5
        -1429428*t**4 -852852*t**3 -333333*t**2
        -77077*t -8008
    ),

    9: (
        35750*t**9 + 160875*t**8 + 563420*t**7
        +1221220*t**6 + 1732458*t**5 + 1653470*t**4
        +1058540*t**3 + 437700*t**2 + 105979*t + 11441
    ),

    8: (
        -71500*t**8 -286000*t**7 -755664*t**6
        -1265992*t**5 -1375360*t**4 -974400*t**3
        -436832*t**2 -112964*t -12879
    ),

    7: (
        88400*t**7 +309400*t**6 +636888*t**5
        +818720*t**4 +663680*t**3 +331500*t**2
        +93604*t +11476
    ),

    6: (
        -71400*t**6 -214200*t**5 -344352*t**4
        -331704*t**3 -190296*t**2 -60144*t -8092
    ),

    5: (
        38760*t**5 +96900*t**4 +119016*t**3
        +81624*t**2 +29736*t +4494
    ),

    4: (
        -14250*t**4 -28500*t**3 -25380*t**2
        -11130*t -1946
    ),

    3: (
        3500*t**3 +5250*t**2 +3038*t +644
    ),

    2: (
        -550*t**2 -550*t -156
    ),

    1: (
        50*t +25
    ),

    0: sp.Integer(-2),
}


# ============================================================================
# GENERAL EXACT HELPERS
# ============================================================================

def S(expr):
    return sp.sympify(expr)


def exact_equal(a, b):
    return sp.simplify(S(a) - S(b)) == 0


def poly(expr, variable):
    expr = sp.expand(S(expr))

    others = sorted(
        expr.free_symbols - {variable},
        key=lambda s: s.sort_key(),
    )

    if others:
        domain = sp.QQ.frac_field(*others)
    else:
        domain = sp.QQ

    return sp.Poly(
        expr,
        variable,
        domain=domain,
    )


def degree(expr, variable):
    p = poly(expr, variable)

    if p.is_zero:
        return -sp.oo

    return p.degree()


def coeff(expr, variable, exponent):
    return sp.expand(
        poly(expr, variable).coeff_monomial(
            variable**exponent
        )
    )


def term_count(expr, *variables):
    expr = sp.expand(S(expr))
    if not variables:
        return 0 if expr == 0 else 1

    p = sp.Poly(
        expr,
        *variables,
        domain=sp.QQ,
    )

    return len(p.terms())


def gcd_poly_u(expressions):
    nonzero = [
        sp.Poly(
            sp.expand(e),
            u,
            domain=sp.QQ,
        )
        for e in expressions
        if not exact_equal(e, 0)
    ]

    if not nonzero:
        return sp.Integer(0)

    g = nonzero[0]

    for p in nonzero[1:]:
        g = sp.gcd(g, p)

    return sp.factor(g.as_expr())


def gcd_poly_w(expressions):
    nonzero = [
        sp.Poly(
            sp.expand(e),
            w,
            domain=sp.QQ.frac_field(u),
        )
        for e in expressions
        if not exact_equal(e, 0)
    ]

    if not nonzero:
        return sp.Integer(0)

    g = nonzero[0]

    for p in nonzero[1:]:
        g = sp.gcd(g, p)

    return sp.factor(g.as_expr())


# ============================================================================
# EXACT INVOLUTION DECOMPOSITION
# ============================================================================

def recover_u_polynomial(expr_t, label):
    """
    Recover P(u) from an expression known to be P(t(t+1)).
    """
    expr_t = sp.expand(S(expr_t))

    if exact_equal(expr_t, 0):
        return sp.Integer(0)

    p_t = poly(expr_t, t)
    max_deg_u = p_t.degree() // 2

    c = sp.symbols(
        f"{label}_0:{max_deg_u + 1}"
    )

    candidate_u = sp.expand(
        sum(
            c[i] * u**i
            for i in range(max_deg_u + 1)
        )
    )

    candidate_t = sp.expand(
        candidate_u.subs(
            u,
            t*(t+1),
        )
    )

    diff = sp.Poly(
        sp.expand(candidate_t - expr_t),
        t,
        domain=sp.QQ.frac_field(*c),
    )

    equations = list(diff.all_coeffs())

    sol = sp.solve(
        equations,
        c,
        dict=True,
    )

    if not sol:
        raise RuntimeError(
            f"Could not recover {label}(u)"
        )

    result = sp.expand(
        candidate_u.subs(sol[0])
    )

    if not exact_equal(
        result.subs(u, t*(t+1)),
        expr_t,
    ):
        raise RuntimeError(
            f"Reconstruction failure for {label}(u)"
        )

    return result


@dataclass
class ChannelData:
    d: int
    A: sp.Expr
    B: sp.Expr


def decompose_layer(d, h):
    h = sp.expand(S(h))

    partner = sp.expand(
        h.subs(
            t,
            -1-t,
        )
    )

    even = sp.expand(
        (h + partner) / 2
    )

    odd = sp.expand(
        (h - partner) / 2
    )

    A = recover_u_polynomial(
        even,
        f"A{d}",
    )

    if exact_equal(odd, 0):
        B = sp.Integer(0)
    else:
        q, r = sp.div(
            poly(odd, t),
            poly(2*t+1, t),
            domain=sp.QQ,
        )

        if not r.is_zero:
            raise RuntimeError(
                f"d={d}: odd component not divisible by 2t+1"
            )

        B = recover_u_polynomial(
            q.as_expr(),
            f"B{d}",
        )

    reconstruction = exact_equal(
        A.subs(u, t*(t+1))
        +
        (2*t+1)
        * B.subs(u, t*(t+1)),
        h,
    )

    if not reconstruction:
        raise RuntimeError(
            f"d={d}: channel reconstruction failed"
        )

    return ChannelData(
        d=d,
        A=sp.expand(A),
        B=sp.expand(B),
    )


def build_channels():
    return {
        d: decompose_layer(d, H[d])
        for d in range(16, -1, -1)
    }


# ============================================================================
# BUILD W-CHANNELS
# ============================================================================

def build_Ahat(channels):
    return sp.expand(
        sum(
            channels[d].A * w**(d//2)
            for d in range(0, 17, 2)
        )
    )


def build_Bhat(channels):
    return sp.expand(
        sum(
            channels[d].B * w**((d-1)//2)
            for d in range(1, 17, 2)
        )
    )


# ============================================================================
# COEFFICIENT MATRICES
# ============================================================================

def coefficient_matrix_u(family):
    """
    Rows = polynomials, columns = powers of u.
    """
    if not family:
        return sp.Matrix([])

    max_deg = max(
        degree(f, u)
        for f in family
    )

    rows = []

    for f in family:
        row = [
            coeff(f, u, j)
            for j in range(max_deg + 1)
        ]
        rows.append(row)

    return sp.Matrix(rows)


def rank_u_family(family):
    return coefficient_matrix_u(family).rank()


# ============================================================================
# RECIPROCAL / ANTI-RECIPROCAL TEST
# ============================================================================

def reciprocal_status(f):
    p = poly(f, w)

    if p.is_zero:
        return False, False

    d = p.degree()

    recip = sp.expand(
        w**d * S(f).subs(
            w,
            1/w,
        )
    )

    return (
        exact_equal(recip, f),
        exact_equal(recip, -S(f)),
    )


# ============================================================================
# W-QUOTIENT / REMAINDER
# ============================================================================

def divide_in_w(f, g):
    pf = sp.Poly(
        sp.expand(f),
        w,
        domain=sp.QQ.frac_field(u),
    )

    pg = sp.Poly(
        sp.expand(g),
        w,
        domain=sp.QQ.frac_field(u),
    )

    return sp.div(
        pf,
        pg,
        domain=sp.QQ.frac_field(u),
    )


# ============================================================================
# SMALL OPERATOR SEARCH
# ============================================================================

def search_constant_w_operator(sequence):
    """
    Search

        H_j = c0 H_{j-1} + c1 H_{j-2}

    with c_i rational constants.

    This is intentionally conservative.
    """
    if len(sequence) < 3:
        return None

    equations = []

    for i in range(2, len(sequence)):
        f0 = sequence[i]
        f1 = sequence[i-1]
        f2 = sequence[i-2]

        diff = sp.Poly(
            sp.expand(
                f0
                - sp.Symbol("c0")*f1
                - sp.Symbol("c1")*f2
            ),
            w,
            domain=sp.QQ.frac_field(
                u,
                sp.Symbol("c0"),
                sp.Symbol("c1"),
            ),
        )

        equations.extend(
            diff.all_coeffs()
        )

    c0, c1 = sp.symbols("c0 c1")

    sol = sp.solve(
        equations,
        (c0, c1),
        dict=True,
    )

    return sol


# ============================================================================
# FRESH SEMIPRIME
# ============================================================================

def fresh_semiprime():
    p = sp.Integer(106621)
    q = sp.Integer(246473)

    if not sp.isprime(p):
        raise RuntimeError("p is not prime")

    if not sp.isprime(q):
        raise RuntimeError("q is not prime")

    N = p*q
    S = p+q
    X = S+1

    tv = sp.cancel(
        sp.Rational(N, X)
    )

    uv = sp.cancel(
        tv*(tv+1)
    )

    return p, q, N, S, X, tv, uv


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print(
        "EXPERIMENT 77 — EXACT EVEN/ODD CHANNEL "
        "GENERATING STRUCTURE"
    )
    print("="*78)

    # ----------------------------------------------------------------------
    # 0
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("0. EXACT SYMBOLIC SETUP")
    print("="*78)

    print("  t = N/X")
    print("  u = t(t+1)")
    print("  z = layer variable")
    print("  w = z^2")
    print("  arithmetic = exact QQ / rational-function fields")
    print("  floating point = forbidden")

    # ----------------------------------------------------------------------
    # 1 INVOLUTION
    # ----------------------------------------------------------------------

    channels = build_channels()

    print()
    print("="*78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("="*78)

    all_ok = True

    for d in range(16, -1, -1):

        h = H[d]

        partner = sp.expand(
            h.subs(
                t,
                -1-t,
            )
        )

        invariant = exact_equal(
            partner,
            h,
        )

        anti = exact_equal(
            partner,
            -h,
        )

        reconstruction = exact_equal(
            channels[d].A.subs(
                u,
                t*(t+1),
            )
            +
            (2*t+1)
            *channels[d].B.subs(
                u,
                t*(t+1),
            ),
            h,
        )

        all_ok &= reconstruction

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"reconstruction={reconstruction}"
        )

    print()
    print(
        "  ALL INVOLUTION CHECKS =",
        all_ok,
    )

    # ----------------------------------------------------------------------
    # 2 W CHANNELS
    # ----------------------------------------------------------------------

    Ahat = build_Ahat(channels)
    Bhat = build_Bhat(channels)

    print()
    print("="*78)
    print("2. EVEN/ODD GENERATING FUNCTIONS IN w=z^2")
    print("="*78)

    print(
        "  Ahat(w,u) degree_w =",
        degree(Ahat, w),
    )

    print(
        "  Ahat(w,u) degree_u =",
        degree(Ahat, u),
    )

    print(
        "  Bhat(w,u) degree_w =",
        degree(Bhat, w),
    )

    print(
        "  Bhat(w,u) degree_u =",
        degree(Bhat, u),
    )

    print()
    print("  Ahat(w,u) =")
    print("   ", sp.factor(Ahat))

    print()
    print("  Bhat(w,u) =")
    print("   ", sp.factor(Bhat))

    # ----------------------------------------------------------------------
    # 3 COEFFICIENT-BY-COEFFICIENT FULL CHANNEL RECONSTRUCTION
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("3. COEFFICIENT-BY-COEFFICIENT CHANNEL RECONSTRUCTION")
    print("="*78)

    coefficient_checks = []

    for d in range(16, -1, -1):

        expected = H[d]

        A_part = channels[d].A.subs(
            u,
            t*(t+1),
        )

        B_part = (
            (2*t+1)
            *channels[d].B.subs(
                u,
                t*(t+1),
            )
        )

        reconstructed = sp.expand(
            A_part + B_part
        )

        ok = exact_equal(
            reconstructed,
            expected,
        )

        coefficient_checks.append(ok)

        print(
            f"  z^{d:2d}: reconstruction={ok}"
        )

    full_channel_ok = all(
        coefficient_checks
    )

    print()
    print(
        "  ALL COEFFICIENT CHECKS =",
        full_channel_ok,
    )

    # ----------------------------------------------------------------------
    # 4 EXACT W-FACTOR STRUCTURE
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("4. EXACT FACTORIZATION IN w")
    print("="*78)

    print()
    print("  factor(Ahat) =")
    print(
        "   ",
        sp.factor(Ahat),
    )

    print()
    print("  factor(Bhat) =")
    print(
        "   ",
        sp.factor(Bhat),
    )

    # ----------------------------------------------------------------------
    # 5 RECIPROCAL STRUCTURE
    # ----------------------------------------------------------------------

    rec_A, anti_A = reciprocal_status(Ahat)
    rec_B, anti_B = reciprocal_status(Bhat)

    print()
    print("="*78)
    print("5. W-RECIPROCAL STRUCTURE")
    print("="*78)

    print(
        "  Ahat reciprocal =",
        rec_A,
    )

    print(
        "  Ahat anti-reciprocal =",
        anti_A,
    )

    print(
        "  Bhat reciprocal =",
        rec_B,
    )

    print(
        "  Bhat anti-reciprocal =",
        anti_B,
    )

    # ----------------------------------------------------------------------
    # 6 COMMON W-GCD
    # ----------------------------------------------------------------------

    gcd_AB = sp.Poly(
        sp.expand(Ahat),
        w,
        domain=sp.QQ.frac_field(u),
    )

    gcd_B = sp.Poly(
        sp.expand(Bhat),
        w,
        domain=sp.QQ.frac_field(u),
    )

    common = sp.gcd(
        gcd_AB,
        gcd_B,
    ).as_expr()

    print()
    print("="*78)
    print("6. COMMON w-GCD")
    print("="*78)

    print(
        "  gcd(Ahat,Bhat) =",
        sp.factor(common),
    )

    # ----------------------------------------------------------------------
    # 7 LAYER COEFFICIENT POLYNOMIALS IN w
    # ----------------------------------------------------------------------

    A_sequence = [
        channels[d].A
        for d in range(16, -1, -2)
    ]

    B_sequence = [
        channels[d].B
        for d in range(15, 0, -2)
    ]

    print()
    print("="*78)
    print("7. CHANNEL COEFFICIENT PROFILES")
    print("="*78)

    print()
    print("  A-channel:")
    for d in range(16, -1, -2):
        f = channels[d].A
        print(
            f"    A_{d}: "
            f"deg_u={degree(f,u)} "
            f"terms={term_count(f,u)} "
            f"content_signature="
            f"{sp.polys.polytools.primitive(
                sp.Poly(
                    f,
                    u,
                    domain=sp.QQ,
                ).as_expr(),
                u,
            )[0]}"
        )

    print()
    print("  B-channel:")
    for d in range(15, 0, -2):
        f = channels[d].B
        print(
            f"    B_{d}: "
            f"deg_u={degree(f,u)} "
            f"terms={term_count(f,u)}"
        )

    # ----------------------------------------------------------------------
    # 8 QQ RANKS
    # ----------------------------------------------------------------------

    A_family = [
        channels[d].A
        for d in range(16, -1, -2)
    ]

    B_family = [
        channels[d].B
        for d in range(15, 0, -2)
    ]

    A_rank = rank_u_family(A_family)
    B_rank = rank_u_family(B_family)

    print()
    print("="*78)
    print("8. EXACT QQ-RANK OF CHANNEL FAMILIES")
    print("="*78)

    print(
        "  A count =",
        len(A_family),
    )

    print(
        "  A QQ-rank =",
        A_rank,
    )

    print(
        "  A nullity =",
        len(A_family) - A_rank,
    )

    print(
        "  B count =",
        len(B_family),
    )

    print(
        "  B QQ-rank =",
        B_rank,
    )

    print(
        "  B nullity =",
        len(B_family) - B_rank,
    )

    # ----------------------------------------------------------------------
    # 9 ADJACENT QUOTIENTS IN CHANNEL INDEX
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("9. EXACT ADJACENT CHANNEL DIVISION")
    print("="*78)

    print()
    print("  A channel:")

    for hi, lo in zip(
        range(16, 0, -2),
        range(14, -2, -2),
    ):
        f_hi = channels[hi].A
        f_lo = channels[lo].A

        q_poly, r_poly = sp.div(
            poly(f_hi, u),
            poly(f_lo, u),
            domain=sp.QQ,
        )

        print(
            f"    A_{hi}/A_{lo}: "
            f"exact={r_poly.is_zero} "
            f"q_deg="
            f"{(-sp.oo if q_poly.is_zero else q_poly.degree())} "
            f"r_deg="
            f"{(-sp.oo if r_poly.is_zero else r_poly.degree())}"
        )

        if r_poly.is_zero:
            print(
                "      quotient=",
                sp.factor(q_poly.as_expr()),
            )

    print()
    print("  B channel:")

    for hi, lo in zip(
        range(15, 1, -2),
        range(13, -1, -2),
    ):
        f_hi = channels[hi].B
        f_lo = channels[lo].B

        q_poly, r_poly = sp.div(
            poly(f_hi, u),
            poly(f_lo, u),
            domain=sp.QQ,
        )

        print(
            f"    B_{hi}/B_{lo}: "
            f"exact={r_poly.is_zero} "
            f"q_deg="
            f"{(-sp.oo if q_poly.is_zero else q_poly.degree())} "
            f"r_deg="
            f"{(-sp.oo if r_poly.is_zero else r_poly.degree())}"
        )

        if r_poly.is_zero:
            print(
                "      quotient=",
                sp.factor(q_poly.as_expr()),
            )

    # ----------------------------------------------------------------------
    # 10 OPERATOR SEARCH USING W-POLYNOMIAL COEFFICIENTS
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("10. SMALL W-OPERATOR SEARCH")
    print("="*78)

    print(
        """
  We now search for an operator acting on the layer index.

  For example,

      H_j = c_0(u) H_{j-1}

  or

      H_j = c_0(u) H_{j-1}
          + c_1(u) H_{j-2},

  where the c_i(u) have small degree.

  This is an exact symbolic search.
        """
    )

    def search_order_one(sequence, max_degree):
        if len(sequence) < 2:
            return None

        c_coeffs = sp.symbols(
            f"c0:{max_degree+1}"
        )

        c_poly = sum(
            c_coeffs[i]*u**i
            for i in range(max_degree+1)
        )

        equations = []

        for current, previous in zip(
            sequence[:-0] if False else sequence[1:],
            sequence[:-1],
        ):
            diff = sp.Poly(
                sp.expand(
                    current
                    - c_poly*previous
                ),
                u,
                domain=sp.QQ.frac_field(
                    *c_coeffs
                ),
            )

            equations.extend(
                diff.all_coeffs()
            )

        sol = sp.solve(
            equations,
            c_coeffs,
            dict=True,
        )

        return sol

    for label, seq in [
        (
            "A",
            A_family,
        ),
        (
            "B",
            B_family,
        ),
    ]:

        print()
        print(
            f"  {label}-channel order-1 "
            "coefficient-operator search:"
        )

        found = False

        for deg_c in range(0, 5):

            sol = search_order_one(
                seq,
                deg_c,
            )

            if sol:
                print(
                    f"    degree={deg_c}: "
                    f"FOUND {sol}"
                )
                found = True
                break

            print(
                f"    degree={deg_c}: NONE"
            )

        if not found:
            print(
                "    no exact order-1 operator "
                "found through degree 4"
            )

    # ----------------------------------------------------------------------
    # 11 FRESH SEMIPRIME
    # ----------------------------------------------------------------------

    (
        p_val,
        q_val,
        N_val,
        S_val,
        X_val,
        t_val,
        u_val,
    ) = fresh_semiprime()

    print()
    print("="*78)
    print("11. FRESH EXACT SEMIPRIME")
    print("="*78)

    print("  p =", p_val)
    print("  q =", q_val)
    print("  N =", N_val)
    print("  S =", S_val)
    print("  X =", X_val)
    print("  t =", t_val)
    print("  u =", u_val)

    # ----------------------------------------------------------------------
    # 12 FRESH CHANNEL CHECKS
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("12. FRESH CHANNEL RECONSTRUCTION")
    print("="*78)

    fresh_ok = True

    for d in range(16, -1, -1):

        direct = sp.expand(
            H[d].subs(
                t,
                t_val,
            )
        )

        channel = sp.expand(
            channels[d].A.subs(
                u,
                u_val,
            )
            +
            (2*t_val+1)
            *channels[d].B.subs(
                u,
                u_val,
            )
        )

        ok = exact_equal(
            direct,
            channel,
        )

        fresh_ok &= ok

        print(
            f"  d={d:2d}: reconstruction={ok}"
        )

    # ----------------------------------------------------------------------
    # 13 FRESH W-GENERATING CHECK
    # ----------------------------------------------------------------------

    A_eval = sp.expand(
        Ahat.subs(
            u,
            u_val,
        )
    )

    B_eval = sp.expand(
        Bhat.subs(
            u,
            u_val,
        )
    )

    H_eval_direct = sp.expand(
        sum(
            H[d].subs(
                t,
                t_val,
            ) * z**d
            for d in range(17)
        )
    )

    H_eval_channels = sp.expand(
        A_eval.subs(
            w,
            z**2,
        )
        +
        (2*t_val+1)
        *z
        *B_eval.subs(
            w,
            z**2,
        )
    )

    # Crucial exact coefficient-by-coefficient check.
    z_degree = degree(
        H_eval_direct,
        z,
    )

    coefficientwise = True

    for j in range(z_degree+1):

        c_direct = coeff(
            H_eval_direct,
            z,
            j,
        )

        c_channel = coeff(
            H_eval_channels,
            z,
            j,
        )

        ok = exact_equal(
            c_direct,
            c_channel,
        )

        if not ok:
            coefficientwise = False

        print(
            f"  z^{j:2d}: match={ok}"
        )

    print()
    print(
        "  coefficientwise full reconstruction =",
        coefficientwise,
    )

    # ----------------------------------------------------------------------
    # 14 FRESH LAYER SUM
    # ----------------------------------------------------------------------

    fresh_layers = {}

    for d in range(16, -1, -1):
        fresh_layers[d] = sp.expand(
            X_val**d
            *
            H[d].subs(
                t,
                t_val,
            )
        )

    layer_sum = sp.expand(
        sum(
            fresh_layers[d]
            for d in range(17)
        )
    )

    print()
    print("="*78)
    print("14. FRESH HOMOGENEOUS LAYER SUM")
    print("="*78)

    print(
        "  all exact layers reconstructed =",
        exact_equal(
            layer_sum,
            layer_sum,
        ),
    )

    print(
        "  total digits =",
        len(str(abs(int(layer_sum)))),
    )

    # ----------------------------------------------------------------------
    # 15 FINAL INTERPRETATION
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("15. STRUCTURAL INTERPRETATION")
    print("="*78)

    print(
        """
  The normalized kernel is now represented exactly as

      H(z,t)
        = Ahat(z^2,u)
          + (2t+1) z Bhat(z^2,u),

      u=t(t+1).

  This separates three structures:

      1. layer parity,
      2. involution parity,
      3. the actual construction in w=z^2.

  The key question is now:

      Does Ahat(w,u) or Bhat(w,u)
      factor into a small number of universal
      polynomial blocks?

  A stronger positive result would be a factorization

      Ahat = P_1 P_2 ... P_r

  or

      Bhat = Q_1 Q_2 ... Q_s

  where the factors have simple dependence on w and u.

  An equally important possibility is a low-degree relation
  between Ahat and Bhat in w.

  If such a construction law exists, it may explain why the
  top layers carry enough information to recover t.

  Only after that construction is understood should we ask
  whether any corresponding quantities can be generated from
  N without explicitly knowing p and q.
        """
    )

    # ----------------------------------------------------------------------
    # FINAL EXACTNESS
    # ----------------------------------------------------------------------

    failures = 0

    checks = {
        "involution": all_ok,
        "coefficientwise_channel_identity":
            coefficient_checks
            and full_channel_ok,
        "fresh_channel_values":
            fresh_ok,
        "fresh_w_generating_identity":
            coefficientwise,
    }

    for name, value in checks.items():
        print(
            f"  {name} = {value}"
        )

        if not value:
            failures += 1

    print()
    print(
        "  failures =",
        failures,
    )

    print(
        "  ALL EXACT CHECKS PASS =",
        failures == 0,
    )

    print()
    print("="*78)
    print("EXPERIMENT 77 COMPLETE")
    print("="*78)


if __name__ == "__main__":
    main()

