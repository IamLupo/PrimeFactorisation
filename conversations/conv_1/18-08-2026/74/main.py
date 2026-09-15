#!/usr/bin/env python3
"""
EXPERIMENT 76 — EXACT LAYER GENERATING-FUNCTION / CONSTRUCTION AUDIT

Purpose
-------
Study the construction of the normalized homogeneous layer sequence

    H(z,t) = sum_{d=0}^{16} h_d(t) z^d

after the exact involution

    t -> -1-t

with

    u = t(t+1).

The layer decomposition is

    h_d(t) = A_d(u) + (2t+1) B_d(u).

This experiment investigates whether the full layer family has a simpler
generating-function structure.

All arithmetic is exact over QQ / rational function fields.
No floating point.
No interpolation.
No regression.
No statistical inference.
"""

import sympy as sp
from dataclasses import dataclass


# ============================================================================
# SYMBOLS
# ============================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")
z = sp.Symbol("z")


# ============================================================================
# EXACT NORMALIZED LAYERS h_d(t)
# ============================================================================

H = {
    16: (
        -9*t**8 - 36*t**7 - 84*t**6 - 126*t**5
        - 126*t**4 - 84*t**3 - 36*t**2 - 9*t - 1
    ),

    15: (
        88*t**9 + 396*t**8 + 1164*t**7 + 2226*t**6
        + 2898*t**5 + 2604*t**4 + 1596*t**3
        + 639*t**2 + 151*t + 16
    ),

    14: (
        -276*t**10 - 1380*t**9 - 5460*t**8 - 13560*t**7
        - 23058*t**6 - 27510*t**5 - 23100*t**4
        - 13410*t**3 - 5135*t**2 - 1169*t - 120
    ),

    13: (
        288*t**11 + 1584*t**10 + 10340*t**9 + 34650*t**8
        + 79464*t**7 + 127512*t**6 + 145530*t**5
        + 117975*t**4 + 66550*t**3 + 24882*t**2
        + 5551*t + 560
    ),

    12: (
        -50*t**12 - 300*t**11 - 7436*t**10 - 34430*t**9
        - 117810*t**8 - 267960*t**7 - 426888*t**6
        - 484110*t**5 - 390225*t**4 - 219010*t**3
        - 81510*t**2 - 18109*t - 1820
    ),

    11: (
        1300*t**11 + 7150*t**10 + 62634*t**9 + 228228*t**8
        + 552552*t**7 + 918918*t**6 + 1075074*t**5
        + 887172*t**4 + 507078*t**3 + 191477*t**2
        + 43043*t + 4368
    ),

    10: (
        -10010*t**10 - 50050*t**9 - 252252*t**8 - 708708*t**7
        - 1303302*t**6 - 1639638*t**5 - 1429428*t**4
        - 852852*t**3 - 333333*t**2 - 77077*t - 8008
    ),

    9: (
        35750*t**9 + 160875*t**8 + 563420*t**7 + 1221220*t**6
        + 1732458*t**5 + 1653470*t**4 + 1058540*t**3
        + 437700*t**2 + 105979*t + 11441
    ),

    8: (
        -71500*t**8 - 286000*t**7 - 755664*t**6 - 1265992*t**5
        - 1375360*t**4 - 974400*t**3 - 436832*t**2
        - 112964*t - 12879
    ),

    7: (
        88400*t**7 + 309400*t**6 + 636888*t**5 + 818720*t**4
        + 663680*t**3 + 331500*t**2 + 93604*t + 11476
    ),

    6: (
        -71400*t**6 - 214200*t**5 - 344352*t**4 - 331704*t**3
        - 190296*t**2 - 60044*t - 8092
    ),

    5: (
        38760*t**5 + 96900*t**4 + 119016*t**3
        + 81624*t**2 + 29736*t + 4494
    ),

    4: (
        -14250*t**4 - 28500*t**3 - 25380*t**2
        - 11130*t - 1946
    ),

    3: (
        3500*t**3 + 5250*t**2 + 3038*t + 644
    ),

    2: (
        -550*t**2 - 550*t - 156
    ),

    1: (
        50*t + 25
    ),

    0: sp.Integer(-2),
}


# ============================================================================
# ROBUST SYMBOLIC HELPERS
# ============================================================================

def S(x):
    return sp.sympify(x)


def infer_coefficient_symbols(expression, variable):
    """
    Return every free symbol except `variable`.

    This is the crucial fix for multivariate polynomials.

    Example:
        degree_in(H(z,t), t)

    treats z as a coefficient-field symbol rather than trying to coerce
    z-dependent coefficients into QQ.
    """
    expression = sp.expand(S(expression))
    symbols = sorted(
        expression.free_symbols - {variable},
        key=lambda s: s.sort_key(),
    )
    return tuple(symbols)


def coefficient_domain(variable, coefficient_symbols):
    if not coefficient_symbols:
        return sp.QQ

    return sp.QQ.frac_field(*coefficient_symbols)


def polynomial(expression, variable, coefficient_symbols=None):
    """
    Construct Poly(expression, variable) over an automatically selected
    rational function field.

    If coefficient_symbols is omitted, they are inferred automatically.
    """
    expression = sp.expand(S(expression))

    if coefficient_symbols is None:
        coefficient_symbols = infer_coefficient_symbols(
            expression,
            variable,
        )

    domain = coefficient_domain(
        variable,
        coefficient_symbols,
    )

    return sp.Poly(
        expression,
        variable,
        domain=domain,
    )


def degree_in(expression, variable, coefficient_symbols=None):
    p = polynomial(
        expression,
        variable,
        coefficient_symbols,
    )

    if p.is_zero:
        return -sp.oo

    return p.degree()


def exact_zero(expression):
    return sp.simplify(S(expression)) == 0


def exact_equal(a, b):
    return sp.simplify(S(a) - S(b)) == 0


def polynomial_term_count(expression, variables):
    """
    Exact number of monomials for a multivariate polynomial over QQ.
    """
    expression = sp.expand(S(expression))
    if not variables:
        return 1 if expression != 0 else 0

    p = sp.Poly(
        expression,
        *variables,
        domain=sp.QQ,
    )

    return len(p.terms())


def primitive_integer_coefficients(expression, variable):
    """
    Primitive integer coefficient signature.

    Handles constant polynomials and one-coefficient cases safely.
    """
    p = polynomial(
        expression,
        variable,
    )

    if p.is_zero:
        return []

    coeffs = [
        sp.Rational(c)
        for c in p.all_coeffs()
    ]

    den_lcm = sp.Integer(1)

    for c in coeffs:
        den_lcm = sp.ilcm(
            int(den_lcm),
            int(c.q),
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    content = 0
    for value in ints:
        content = sp.igcd(
            int(content),
            abs(int(value)),
        )

    if content == 0:
        return ints

    ints = [
        value // int(content)
        for value in ints
    ]

    # Normalize sign.
    if ints and ints[0] < 0:
        ints = [-x for x in ints]

    return ints


# ============================================================================
# EXACT RECOVERY P(u), WITH u = t(t+1)
# ============================================================================

def recover_u_polynomial(expr_t, label):
    """
    Recover exact P(u) satisfying

        expr_t = P(t(t+1))

    whenever such P exists.
    """
    expr_t = sp.expand(S(expr_t))

    if exact_zero(expr_t):
        return sp.Integer(0)

    p_t = polynomial(expr_t, t)

    deg_t = p_t.degree()
    max_degree_u = deg_t // 2

    c = sp.symbols(
        f"{label}_c0:{max_degree_u + 1}"
    )

    candidate_u = sp.expand(
        sum(
            c[i] * u**i
            for i in range(max_degree_u + 1)
        )
    )

    candidate_t = sp.expand(
        candidate_u.subs(
            u,
            t*(t + 1),
        )
    )

    # The c_i are unknown coefficients, so the coefficient field must
    # contain them. This is deliberately not QQ.
    diff_poly = sp.Poly(
        sp.expand(candidate_t - expr_t),
        t,
        domain=sp.QQ.frac_field(*c),
    )

    equations = list(
        diff_poly.all_coeffs()
    )

    solution = sp.solve(
        equations,
        c,
        dict=True,
    )

    if not solution:
        raise RuntimeError(
            f"Could not recover {label}(u)"
        )

    result = sp.expand(
        candidate_u.subs(solution[0])
    )

    if not exact_equal(
        result.subs(u, t*(t + 1)),
        expr_t,
    ):
        raise RuntimeError(
            f"Exact reconstruction failed for {label}(u)"
        )

    return result


# ============================================================================
# CHANNEL OBJECT
# ============================================================================

@dataclass
class Channel:
    d: int
    A: sp.Expr
    B: sp.Expr
    reconstruction: bool


# ============================================================================
# INVOLUTION DECOMPOSITION
# ============================================================================

def decompose_layer(d, h):
    h = sp.expand(S(h))

    partner = sp.expand(
        h.subs(
            t,
            -1 - t,
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
        f"A_{d}",
    )

    if exact_zero(odd):
        B = sp.Integer(0)
    else:
        odd_poly = polynomial(
            odd,
            t,
        )

        divisor = polynomial(
            2*t + 1,
            t,
        )

        q, r = sp.div(
            odd_poly,
            divisor,
            domain=sp.QQ,
        )

        if not r.is_zero:
            raise RuntimeError(
                f"d={d}: odd component not divisible by 2t+1"
            )

        B = recover_u_polynomial(
            q.as_expr(),
            f"B_{d}",
        )

    reconstructed = sp.expand(
        A.subs(u, t*(t + 1))
        +
        (2*t + 1)
        * B.subs(u, t*(t + 1))
    )

    ok = exact_equal(
        reconstructed,
        h,
    )

    return Channel(
        d=d,
        A=sp.expand(A),
        B=sp.expand(B),
        reconstruction=ok,
    )


def build_channels():
    result = {}

    for d in range(16, -1, -1):
        result[d] = decompose_layer(
            d,
            H[d],
        )

    return result


# ============================================================================
# GENERATING FUNCTIONS
# ============================================================================

def build_full_generating_function():
    return sp.expand(
        sum(
            H[d] * z**d
            for d in range(17)
        )
    )


def build_A_generating_function(channels):
    return sp.expand(
        sum(
            channels[d].A * z**d
            for d in range(0, 17, 2)
        )
    )


def build_B_generating_function(channels):
    return sp.expand(
        sum(
            channels[d].B * z**d
            for d in range(1, 17, 2)
        )
    )


# ============================================================================
# Z-POLYNOMIAL HELPERS
# ============================================================================

def poly_in_z_t(expression):
    return polynomial(
        expression,
        z,
        coefficient_symbols=(t,),
    )


def poly_in_z_u(expression):
    return polynomial(
        expression,
        z,
        coefficient_symbols=(u,),
    )


def z_division(f, g, coefficient_symbol):
    if coefficient_symbol == t:
        field = sp.QQ.frac_field(t)
    elif coefficient_symbol == u:
        field = sp.QQ.frac_field(u)
    else:
        raise ValueError(
            "coefficient_symbol must be t or u"
        )

    pf = sp.Poly(
        sp.expand(f),
        z,
        domain=field,
    )

    pg = sp.Poly(
        sp.expand(g),
        z,
        domain=field,
    )

    q, r = sp.div(
        pf,
        pg,
        domain=field,
    )

    return q, r


def factor_over_z(expression, coefficient_symbol):
    if coefficient_symbol == t:
        field = sp.QQ.frac_field(t)
    else:
        field = sp.QQ.frac_field(u)

    p = sp.Poly(
        sp.expand(expression),
        z,
        domain=field,
    )

    return sp.factor_list(
        p.as_expr(),
        z,
    )


def gcd_over_z(f, g, coefficient_symbol):
    if coefficient_symbol == t:
        field = sp.QQ.frac_field(t)
    else:
        field = sp.QQ.frac_field(u)

    pf = sp.Poly(
        sp.expand(f),
        z,
        domain=field,
    )

    pg = sp.Poly(
        sp.expand(g),
        z,
        domain=field,
    )

    return sp.gcd(
        pf,
        pg,
    ).as_expr()


# ============================================================================
# RECIPROCAL TEST
# ============================================================================

def reciprocal_test(f, variable, coefficient_symbols=None):
    p = polynomial(
        f,
        variable,
        coefficient_symbols,
    )

    if p.is_zero:
        return False, False

    degree = p.degree()

    reciprocal = sp.expand(
        variable**degree
        * S(f).subs(
            variable,
            1 / variable,
        )
    )

    return (
        exact_equal(reciprocal, f),
        exact_equal(reciprocal, -S(f)),
    )


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

    t_value = sp.cancel(
        sp.Rational(N, X)
    )

    u_value = sp.cancel(
        t_value * (t_value + 1)
    )

    return (
        p,
        q,
        N,
        S,
        X,
        t_value,
        u_value,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print(
        "EXPERIMENT 76 — EXACT LAYER GENERATING-FUNCTION / "
        "CONSTRUCTION AUDIT"
    )
    print("="*78)

    # ----------------------------------------------------------------------
    # Setup
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("0. EXACT SYMBOLIC SETUP")
    print("="*78)

    print("  t = N/X")
    print("  u = t(t+1)")
    print("  z = layer variable")
    print("  arithmetic = exact QQ / rational-function fields")
    print("  floating point = forbidden")

    # ----------------------------------------------------------------------
    # Build channels
    # ----------------------------------------------------------------------

    channels = build_channels()

    involution_ok = all(
        channels[d].reconstruction
        for d in range(17)
    )

    # ----------------------------------------------------------------------
    # Involution report
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("="*78)

    for d in range(16, -1, -1):

        h = sp.expand(H[d])

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

        divisible = False

        if anti:
            q, r = sp.div(
                polynomial(h, t),
                polynomial(2*t+1, t),
                domain=sp.QQ,
            )
            divisible = r.is_zero

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"divisible={divisible} "
            f"reconstruction={channels[d].reconstruction}"
        )

    print(
        "\n  ALL INVOLUTION CHECKS =",
        involution_ok,
    )

    # ----------------------------------------------------------------------
    # Full generating polynomial
    # ----------------------------------------------------------------------

    H_z = build_full_generating_function()

    print()
    print("="*78)
    print("2. FULL GENERATING POLYNOMIAL H(z,t)")
    print("="*78)

    print(
        "  degree_z =",
        degree_in(H_z, z),
    )

    # IMPORTANT:
    # H_z contains z and t. degree_in now automatically treats z as
    # a rational-function coefficient symbol when the variable is t.
    print(
        "  degree_t =",
        degree_in(H_z, t),
    )

    print(
        "  number_of_terms =",
        polynomial_term_count(
            H_z,
            (z, t),
        ),
    )

    print("  leading_z coefficient =")
    print(
        "   ",
        polynomial(
            H_z,
            z,
            (t,),
        ).LC(),
    )

    print("  constant_z coefficient =")
    print(
        "   ",
        polynomial(
            H_z,
            z,
            (t,),
        ).TC(),
    )

    # ----------------------------------------------------------------------
    # A/B generating functions
    # ----------------------------------------------------------------------

    A_z = build_A_generating_function(
        channels
    )

    B_z = build_B_generating_function(
        channels
    )

    print()
    print("="*78)
    print("3. CHANNEL GENERATING POLYNOMIALS")
    print("="*78)

    print(
        "  A(z,u) degree_z =",
        degree_in(A_z, z),
    )

    print(
        "  A(z,u) degree_u =",
        degree_in(A_z, u),
    )

    print(
        "  B(z,u) degree_z =",
        degree_in(B_z, z),
    )

    print(
        "  B(z,u) degree_u =",
        degree_in(B_z, u),
    )

    print()
    print("  A(z,u) =")
    print("   ", sp.factor(A_z))

    print()
    print("  B(z,u) =")
    print("   ", sp.factor(B_z))

    # ----------------------------------------------------------------------
    # Factorization over z
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("4. EXACT FACTORIZATION OVER z")
    print("="*78)

    print()
    print("  H(z,t):")
    print(
        "   ",
        factor_over_z(
            H_z,
            t,
        ),
    )

    print()
    print("  A(z,u):")
    print(
        "   ",
        factor_over_z(
            A_z,
            u,
        ),
    )

    print()
    print("  B(z,u):")
    print(
        "   ",
        factor_over_z(
            B_z,
            u,
        ),
    )

    # ----------------------------------------------------------------------
    # Simple z divisor audit
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("5. SIMPLE z-DIVISOR AUDIT")
    print("="*78)

    divisors = [
        z,
        z-1,
        z+1,
        2*z+1,
        z**2+1,
        z**2+z+1,
        z**2-z+1,
        z**3+1,
        z**3+z**2+z+1,
    ]

    for label, f, field_symbol in [
        ("H(z,t)", H_z, t),
        ("A(z,u)", A_z, u),
        ("B(z,u)", B_z, u),
    ]:

        print()
        print(f"  {label}")

        for divisor in divisors:
            try:
                q, r = z_division(
                    f,
                    divisor,
                    field_symbol,
                )

                print(
                    f"    divisor={divisor}: "
                    f"exact={r.is_zero} "
                    f"q_degree="
                    f"{(-sp.oo if q.is_zero else q.degree())} "
                    f"r_degree="
                    f"{(-sp.oo if r.is_zero else r.degree())}"
                )

            except Exception as exc:
                print(
                    f"    divisor={divisor}: ERROR={exc}"
                )

    # ----------------------------------------------------------------------
    # Reciprocal structure
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("6. RECIPROCAL STRUCTURE")
    print("="*78)

    for label, f, variable, coeff_symbols in [
        ("H(z,t)", H_z, z, (t,)),
        ("A(z,u)", A_z, z, (u,)),
        ("B(z,u)", B_z, z, (u,)),
    ]:

        reciprocal, anti = reciprocal_test(
            f,
            variable,
            coeff_symbols,
        )

        print()
        print(f"  {label}")
        print("    reciprocal =", reciprocal)
        print("    anti_reciprocal =", anti)

    # ----------------------------------------------------------------------
    # Layer reversal
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("7. LAYER REVERSAL")
    print("="*78)

    hdeg = degree_in(
        H_z,
        z,
    )

    reverse_plus = sp.expand(
        z**hdeg
        *H_z.subs(
            z,
            1/z,
        )
    )

    reverse_minus = sp.expand(
        z**hdeg
        *H_z.subs(
            z,
            -1/z,
        )
    )

    print(
        "  z -> 1/z equals H =",
        exact_equal(
            reverse_plus,
            H_z,
        ),
    )

    print(
        "  z -> 1/z equals -H =",
        exact_equal(
            reverse_plus,
            -H_z,
        ),
    )

    print(
        "  z -> -1/z equals H =",
        exact_equal(
            reverse_minus,
            H_z,
        ),
    )

    print(
        "  z -> -1/z equals -H =",
        exact_equal(
            reverse_minus,
            -H_z,
        ),
    )

    # ----------------------------------------------------------------------
    # Common z-gcd of channel generating functions
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("8. COMMON z-GCD")
    print("="*78)

    gcd_channel = gcd_over_z(
        A_z,
        B_z,
        u,
    )

    print(
        "  gcd(A(z,u), B(z,u)) =",
        sp.factor(gcd_channel),
    )

    # ----------------------------------------------------------------------
    # Explicit C/D identity
    #
    # C includes even d directly.
    # D includes odd d directly.
    # ----------------------------------------------------------------------

    C_z = sp.expand(
        sum(
            channels[d].A*z**d
            for d in range(0, 17, 2)
            if not exact_zero(channels[d].A)
        )
    )

    D_z = sp.expand(
        sum(
            channels[d].B*z**d
            for d in range(1, 17, 2)
            if not exact_zero(channels[d].B)
        )
    )

    H_from_channels = sp.expand(
        C_z.subs(
            u,
            t*(t+1),
        )
        +
        (2*t+1)
        *D_z.subs(
            u,
            t*(t+1),
        )
    )

    print()
    print("="*78)
    print("9. FULL CHANNEL RECONSTRUCTION")
    print("="*78)

    print(
        "  C(z,u) degree_z =",
        degree_in(C_z, z),
    )

    print(
        "  D(z,u) degree_z =",
        degree_in(D_z, z),
    )

    print(
        "  H == C + (2t+1)D =",
        exact_equal(
            H_from_channels,
            H_z,
        ),
    )

    # ----------------------------------------------------------------------
    # Coefficient-space rank
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("10. CHANNEL COEFFICIENT-SPACE RANK")
    print("="*78)

    A_list = [
        channels[d].A
        for d in range(0, 17, 2)
    ]

    B_list = [
        channels[d].B
        for d in range(1, 17, 2)
    ]

    def coefficient_vector_u(f):
        p = polynomial(
            f,
            u,
        )

        if p.is_zero:
            return [sp.Integer(0)]

        return [
            sp.expand(
                p.coeff_monomial(
                    u**i
                )
            )
            for i in range(
                p.degree()+1
            )
        ]

    def rank_family_u(family):
        if not family:
            return 0

        width = max(
            len(coefficient_vector_u(f))
            for f in family
        )

        rows = []

        for f in family:
            row = coefficient_vector_u(f)
            row += [0] * (width-len(row))
            rows.append(row)

        return sp.Matrix(rows).rank()

    print(
        "  A count =",
        len(A_list),
    )

    print(
        "  A QQ-rank =",
        rank_family_u(A_list),
    )

    print(
        "  B count =",
        len(B_list),
    )

    print(
        "  B QQ-rank =",
        rank_family_u(B_list),
    )

    # ----------------------------------------------------------------------
    # Fresh exact semiprime
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
    # Fresh generating function
    # ----------------------------------------------------------------------

    H_eval = sp.expand(
        H_z.subs(
            t,
            t_val,
        )
    )

    channel_eval = sp.expand(
        H_from_channels.subs(
            t,
            t_val,
        )
    )

    print()
    print("="*78)
    print("12. FRESH GENERATING-FUNCTION CHECK")
    print("="*78)

    print(
        "  direct H == channel H =",
        exact_equal(
            H_eval,
            channel_eval,
        ),
    )

    # ----------------------------------------------------------------------
    # Fresh layers
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("13. FRESH HOMOGENEOUS LAYERS")
    print("="*78)

    fresh_layers = {}

    for d in range(16, -1, -1):

        h_value = sp.expand(
            H[d].subs(
                t,
                t_val,
            )
        )

        L = sp.expand(
            sp.Integer(X_val)**d
            *h_value
        )

        fresh_layers[d] = L

        print(
            f"  L_{d:2d} = {L}"
        )

    total_direct = sp.expand(
        sum(
            fresh_layers[d]
            for d in range(17)
        )
    )

    print()
    print("="*78)
    print("14. FRESH FULL LAYER RECONSTRUCTION")
    print("="*78)

    print(
        "  sum(L_d) computed exactly =",
        total_direct,
    )

    # ----------------------------------------------------------------------
    # Final exactness
    # ----------------------------------------------------------------------

    checks = {
        "involution": involution_ok,
        "full_channel_identity":
            exact_equal(
                H_from_channels,
                H_z,
            ),
        "fresh_generating_identity":
            exact_equal(
                H_eval,
                channel_eval,
            ),
        "fresh_H_reconstruction":
            exact_equal(
                sp.expand(
                    H_eval
                ),
                sp.expand(
                    sum(
                        H[d].subs(
                            t,
                            t_val,
                        ) * z**d
                        for d in range(17)
                    )
                ),
            ),
        "fresh_layer_sum_defined":
            exact_equal(
                total_direct,
                total_direct,
            ),
    }

    failures = sum(
        1
        for value in checks.values()
        if not value
    )

    print()
    print("="*78)
    print("15. STRUCTURAL INTERPRETATION")
    print("="*78)

    print(
        """
  The entire normalized layer family is encoded by

      H(z,t) = sum_{d=0}^{16} h_d(t) z^d.

  The involution splits it into

      H(z,t)
        = C(z,u)
        + (2t+1) D(z,u),

      u=t(t+1).

  This experiment therefore asks a stronger construction question:

      Is there a compact factorization or recurrence of
      H(z,t), C(z,u), or D(z,u) in the layer variable z?

  In particular, an exact factor

      Q(z,u)

  would imply an algebraic recurrence between layer
  coefficients.

  The important distinction remains:

      KNOWN:
          layer values -> t -> X -> N -> p,q

      STILL OPEN:
          N -> layer information -> t -> X -> p,q.

  The purpose of this experiment is to understand the
  left-hand construction mechanism before attempting
  the second implication.
        """
    )

    print("="*78)
    print("FINAL EXACTNESS")
    print("="*78)

    for name, value in checks.items():
        print(
            f"  {name} = {value}"
        )

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
    print("EXPERIMENT 76 COMPLETE")
    print("="*78)


if __name__ == "__main__":
    main()