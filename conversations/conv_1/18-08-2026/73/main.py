#!/usr/bin/env python3
"""
EXPERIMENT 75 — EXACT MULTIPLICATIVE NESTING / LAYER-BUILDING AUDIT

Robust version.

Important implementation rule:
    every symbolic layer is explicitly converted with sp.sympify(...)
    before .subs(), .expand(), Poly(), etc.

This prevents failures such as:

    AttributeError: 'int' object has no attribute 'subs'

All arithmetic is exact over QQ.
No floating point is used.
"""

import sympy as sp
from dataclasses import dataclass
from typing import Dict, List, Tuple


# ======================================================================
# SYMBOLS
# ======================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")
N = sp.Symbol("N")
X = sp.Symbol("X")


def S(expr):
    """Force any Python int/bool/etc. into a SymPy expression."""
    return sp.sympify(expr)


# ======================================================================
# EXACT NORMALIZED LAYER POLYNOMIALS h_d(t)
# ======================================================================

H_T: Dict[int, sp.Expr] = {
    16: S(
        -9*t**8 - 36*t**7 - 84*t**6 - 126*t**5
        - 126*t**4 - 84*t**3 - 36*t**2 - 9*t - 1
    ),

    15: S(
        88*t**9 + 396*t**8 + 1164*t**7 + 2226*t**6
        + 2898*t**5 + 2604*t**4 + 1596*t**3
        + 639*t**2 + 151*t + 16
    ),

    14: S(
        -276*t**10 - 1380*t**9 - 5460*t**8 - 13560*t**7
        - 23058*t**6 - 27510*t**5 - 23100*t**4
        - 13410*t**3 - 5135*t**2 - 1169*t - 120
    ),

    13: S(
        288*t**11 + 1584*t**10 + 10340*t**9 + 34650*t**8
        + 79464*t**7 + 127512*t**6 + 145530*t**5
        + 117975*t**4 + 66550*t**3 + 24882*t**2
        + 5551*t + 560
    ),

    12: S(
        -50*t**12 - 300*t**11 - 7436*t**10 - 34430*t**9
        - 117810*t**8 - 267960*t**7 - 426888*t**6
        - 484110*t**5 - 390225*t**4 - 219010*t**3
        - 81510*t**2 - 18109*t - 1820
    ),

    11: S(
        1300*t**11 + 7150*t**10 + 62634*t**9 + 228228*t**8
        + 552552*t**7 + 918918*t**6 + 1075074*t**5
        + 887172*t**4 + 507078*t**3 + 191477*t**2
        + 43043*t + 4368
    ),

    10: S(
        -10010*t**10 - 50050*t**9 - 252252*t**8 - 708708*t**7
        - 1303302*t**6 - 1639638*t**5 - 1429428*t**4
        - 852852*t**3 - 333333*t**2 - 77077*t - 8008
    ),

    9: S(
        35750*t**9 + 160875*t**8 + 563420*t**7 + 1221220*t**6
        + 1732458*t**5 + 1653470*t**4 + 1058540*t**3
        + 437700*t**2 + 105979*t + 11441
    ),

    8: S(
        -71500*t**8 - 286000*t**7 - 755664*t**6 - 1265992*t**5
        - 1375360*t**4 - 974400*t**3 - 436832*t**2
        - 112964*t - 12879
    ),

    7: S(
        88400*t**7 + 309400*t**6 + 636888*t**5 + 818720*t**4
        + 663680*t**3 + 331500*t**2 + 93604*t + 11476
    ),

    6: S(
        -71400*t**6 - 214200*t**5 - 344352*t**4
        - 331704*t**3 - 190296*t**2 - 60144*t - 8092
    ),

    5: S(
        38760*t**5 + 96900*t**4 + 119016*t**3
        + 81624*t**2 + 29736*t + 4494
    ),

    4: S(
        -14250*t**4 - 28500*t**3 - 25380*t**2
        - 11130*t - 1946
    ),

    3: S(
        3500*t**3 + 5250*t**2 + 3038*t + 644
    ),

    2: S(
        -550*t**2 - 550*t - 156
    ),

    1: S(
        50*t + 25
    ),

    0: S(-2),
}


# ======================================================================
# POLYNOMIAL HELPERS
# ======================================================================

def poly_t(expr) -> sp.Poly:
    return sp.Poly(sp.expand(S(expr)), t, domain=sp.QQ)


def poly_u(expr) -> sp.Poly:
    return sp.Poly(sp.expand(S(expr)), u, domain=sp.QQ)


def degree_u(expr) -> int:
    expr = S(expr)
    p = poly_u(expr)
    return p.degree() if not p.is_zero else -sp.oo


def primitive_integer_coefficients(expr) -> List[int]:
    """
    Primitive integer coefficient vector in descending degree.
    """
    p = poly_u(expr)

    if p.is_zero:
        return [0]

    coeffs = [sp.Rational(c) for c in p.all_coeffs()]

    den = 1
    for c in coeffs:
        den = sp.ilcm(den, int(c.q))

    ints = [int(c * den) for c in coeffs]

    g = 0
    for z in ints:
        g = sp.igcd(g, abs(z))

    if g:
        ints = [z // g for z in ints]

    if ints and ints[0] < 0:
        ints = [-z for z in ints]

    return ints


def exact_equal(a, b) -> bool:
    return sp.simplify(S(a) - S(b)) == 0


# ======================================================================
# EXACT INVOLUTION DECOMPOSITION
# ======================================================================

@dataclass
class ChannelInfo:
    d: int
    A: sp.Expr
    B: sp.Expr
    invariant: bool
    anti_invariant: bool
    divisible_by_linear: bool
    reconstructed: bool


def recover_polynomial_in_u(expr_t, variable_name: str) -> sp.Expr:
    """
    Given a polynomial known to be in u=t(t+1), recover its exact
    polynomial in u.

    No interpolation is used.
    """
    expr_t = S(expr_t)

    if exact_equal(expr_t, 0):
        return S(0)

    p = poly_t(expr_t)
    degree = int(p.degree())

    max_degree_u = degree // 2

    coeffs = sp.symbols(
        f"{variable_name}_0:{max_degree_u + 1}"
    )

    candidate = sum(
        coeffs[i] * u**i
        for i in range(max_degree_u + 1)
    )

    candidate_t = sp.expand(
        candidate.subs(u, t * (t + 1))
    )

    difference = sp.expand(candidate_t - expr_t)

    equations = sp.Poly(
        difference,
        t,
        domain=sp.QQ.frac_field(*coeffs),
    ).all_coeffs()

    equations = [sp.Eq(eq, 0) for eq in equations]

    solution = sp.solve(
        equations,
        coeffs,
        dict=True,
    )

    if not solution:
        raise RuntimeError(
            f"Could not recover polynomial in u for {variable_name}: "
            f"{expr_t}"
        )

    result = sp.expand(
        candidate.subs(solution[0])
    )

    # Exact verification
    verify = sp.expand(
        result.subs(u, t * (t + 1)) - expr_t
    )

    if verify != 0:
        raise RuntimeError(
            f"u-reconstruction failed for {variable_name}"
        )

    return result


def decompose_involution(expr, d: int) -> ChannelInfo:
    """
    Exact decomposition

        h_d(t) = A_d(u) + (2t+1) B_d(u),

        u=t(t+1).

    """
    h = S(expr)

    partner = sp.expand(
        h.subs(t, -1 - t)
    )

    even_part = sp.expand(
        sp.cancel((h + partner) / 2)
    )

    odd_part = sp.expand(
        sp.cancel((h - partner) / 2)
    )

    A = recover_polynomial_in_u(
        even_part,
        f"A{d}"
    )

    if exact_equal(odd_part, 0):
        B = S(0)
    else:
        quotient = sp.cancel(
            odd_part / (2*t + 1)
        )

        # Verify exact divisibility first.
        rem = sp.rem(
            poly_t(odd_part),
            poly_t(2*t + 1),
            domain=sp.QQ,
        )

        if not rem.is_zero:
            raise RuntimeError(
                f"d={d}: odd part not divisible by 2t+1"
            )

        B = recover_polynomial_in_u(
            quotient,
            f"B{d}"
        )

    reconstructed = exact_equal(
        h,
        A.subs(u, t * (t + 1))
        + (2*t + 1) * B.subs(u, t * (t + 1))
    )

    invariant = exact_equal(
        partner,
        h
    )

    anti_invariant = exact_equal(
        partner,
        -h
    )

    divisible = False

    p_h = poly_t(h)
    p_linear = poly_t(2*t + 1)

    if not p_h.is_zero:
        divisible = sp.rem(
            p_h,
            p_linear,
            domain=sp.QQ,
        ).is_zero

    return ChannelInfo(
        d=d,
        A=sp.expand(A),
        B=sp.expand(B),
        invariant=invariant,
        anti_invariant=anti_invariant,
        divisible_by_linear=divisible,
        reconstructed=reconstructed,
    )


# ======================================================================
# EXACT DIVISION
# ======================================================================

@dataclass
class DivisionResult:
    hi: int
    lo: int
    divides: bool
    quotient: sp.Expr
    remainder: sp.Expr
    quotient_degree: int
    remainder_degree: int


def exact_division(
    high_expr,
    low_expr,
    hi: int,
    lo: int,
) -> DivisionResult:

    high_expr = S(high_expr)
    low_expr = S(low_expr)

    if exact_equal(low_expr, 0):
        raise ZeroDivisionError(
            f"Cannot divide by zero channel polynomial A_{lo}/B_{lo}"
        )

    high_poly = poly_u(high_expr)
    low_poly = poly_u(low_expr)

    q, r = sp.div(
        high_poly,
        low_poly,
        domain=sp.QQ,
    )

    q_expr = sp.expand(q.as_expr())
    r_expr = sp.expand(r.as_expr())

    return DivisionResult(
        hi=hi,
        lo=lo,
        divides=r.is_zero,
        quotient=q_expr,
        remainder=r_expr,
        quotient_degree=(
            q.degree()
            if not q.is_zero
            else -sp.oo
        ),
        remainder_degree=(
            r.degree()
            if not r.is_zero
            else -sp.oo
        ),
    )


# ======================================================================
# CHANNEL DATA
# ======================================================================

def build_channels(
    infos: Dict[int, ChannelInfo],
):
    A = {}
    B = {}

    for d in range(16, -1, -1):
        info = infos[d]

        if not exact_equal(info.A, 0):
            A[d] = S(info.A)

        if not exact_equal(info.B, 0):
            B[d] = S(info.B)

    return A, B


# ======================================================================
# DIVISIBILITY GRAPH
# ======================================================================

def divisibility_graph(
    channel: Dict[int, sp.Expr]
):
    keys = sorted(channel.keys(), reverse=True)
    graph = {}

    for hi in keys:
        for lo in keys:
            if hi <= lo:
                continue

            graph[(hi, lo)] = exact_division(
                channel[hi],
                channel[lo],
                hi,
                lo,
            )

    return graph


# ======================================================================
# QUOTIENT RECOGNITION
# ======================================================================

def scalar_multiple_of(
    qexpr,
    candidate,
) -> bool:
    qexpr = S(qexpr)
    candidate = S(candidate)

    if exact_equal(candidate, 0):
        return exact_equal(qexpr, 0)

    pq = poly_u(qexpr)
    pc = poly_u(candidate)

    if pq.degree() != pc.degree():
        return False

    ratio = sp.cancel(
        sp.Rational(pq.LC())
        / sp.Rational(pc.LC())
    )

    return exact_equal(
        qexpr,
        ratio * candidate
    )


def recognize_existing(
    quotient,
    channel: Dict[int, sp.Expr]
) -> List[int]:

    result = []

    for d, candidate in sorted(channel.items()):
        if scalar_multiple_of(
            quotient,
            candidate
        ):
            result.append(d)

    return result


def recognize_products(
    quotient,
    channel: Dict[int, sp.Expr]
):
    result = []

    keys = sorted(channel.keys())

    for i, di in enumerate(keys):
        for dj in keys[i:]:
            product = sp.expand(
                channel[di] * channel[dj]
            )

            if scalar_multiple_of(
                quotient,
                product
            ):
                pq = poly_u(quotient)
                pp = poly_u(product)

                ratio = sp.cancel(
                    sp.Rational(pq.LC())
                    / sp.Rational(pp.LC())
                )

                result.append(
                    (di, dj, ratio)
                )

    return result


# ======================================================================
# GCD OF EXPRESSIONS
# ======================================================================

def gcd_channel(expressions: List[sp.Expr]) -> sp.Expr:
    if not expressions:
        return S(0)

    current = poly_u(expressions[0])

    for expr in expressions[1:]:
        current = sp.gcd(
            current,
            poly_u(expr),
        )

    return sp.expand(
        current.as_expr()
    )


# ======================================================================
# LAYER EVALUATION
# ======================================================================

def layer_value(
    d: int,
    n: int,
    x: int,
) -> sp.Expr:

    x = sp.Integer(x)
    n = sp.Integer(n)

    ratio = sp.cancel(
        n / x
    )

    h = S(H_T[d])

    return sp.expand(
        x**d * h.subs(
            t,
            ratio
        )
    )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 75 — EXACT MULTIPLICATIVE NESTING / "
        "LAYER-BUILDING AUDIT"
    )
    print("=" * 78)

    # --------------------------------------------------------------
    # 0
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)

    print("  t = N/X")
    print("  u = t(t+1)")
    print("  involution = t -> -1-t")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # --------------------------------------------------------------
    # 1. involution
    # --------------------------------------------------------------

    infos = {}

    all_involution_ok = True

    print()
    print("=" * 78)
    print("1. EXACT INVOLUTION DECOMPOSITION")
    print("=" * 78)

    for d in range(16, -1, -1):
        info = decompose_involution(
            H_T[d],
            d,
        )

        infos[d] = info

        ok = (
            info.reconstructed
            and (
                info.invariant
                or info.anti_invariant
            )
        )

        all_involution_ok &= ok

        print(
            f"  d={d:2d}: "
            f"invariant={info.invariant} "
            f"anti_invariant={info.anti_invariant} "
            f"divisible={info.divisible_by_linear} "
            f"reconstruction={info.reconstructed}"
        )

    print()
    print(
        f"  ALL INVOLUTION CHECKS = "
        f"{all_involution_ok}"
    )

    # --------------------------------------------------------------
    # channels
    # --------------------------------------------------------------

    A, B = build_channels(infos)

    A_even = {
        d: A[d]
        for d in sorted(A)
        if d % 2 == 0
    }

    B_odd = {
        d: B[d]
        for d in sorted(B)
        if d % 2 == 1
    }

    print()
    print("-" * 78)
    print("A CHANNEL")
    print("-" * 78)

    for d in sorted(A_even, reverse=True):
        expr = S(A_even[d])

        print(
            f"  A_{d}: "
            f"degree={degree_u(expr)} "
            f"primitive="
            f"{primitive_integer_coefficients(expr)}"
        )

    print()
    print("-" * 78)
    print("B CHANNEL")
    print("-" * 78)

    for d in sorted(B_odd, reverse=True):
        expr = S(B_odd[d])

        print(
            f"  B_{d}: "
            f"degree={degree_u(expr)} "
            f"primitive="
            f"{primitive_integer_coefficients(expr)}"
        )

    # --------------------------------------------------------------
    # 2. complete divisibility
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. COMPLETE SAME-CHANNEL DIVISIBILITY AUDIT")
    print("=" * 78)

    A_graph = divisibility_graph(A_even)
    B_graph = divisibility_graph(B_odd)

    exact_A_edges = []
    exact_B_edges = []

    print()
    print("  A CHANNEL EXACT DIVISIONS")

    for (hi, lo), result in sorted(
        A_graph.items(),
        reverse=True,
    ):
        if result.divides:
            exact_A_edges.append(
                (hi, lo, result.quotient)
            )

            print(
                f"    A_{hi}/A_{lo} = "
                f"{sp.factor(result.quotient)}"
            )

    if not exact_A_edges:
        print("    none")

    print()
    print("  B CHANNEL EXACT DIVISIONS")

    for (hi, lo), result in sorted(
        B_graph.items(),
        reverse=True,
    ):
        if result.divides:
            exact_B_edges.append(
                (hi, lo, result.quotient)
            )

            print(
                f"    B_{hi}/B_{lo} = "
                f"{sp.factor(result.quotient)}"
            )

    if not exact_B_edges:
        print("    none")

    # --------------------------------------------------------------
    # 3. adjacent same parity
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ADJACENT SAME-PARITY NESTING")
    print("=" * 78)

    A_keys = sorted(A_even.keys())
    B_keys = sorted(B_odd.keys())

    print()
    print("  A CHANNEL")

    for lo, hi in zip(
        A_keys,
        A_keys[1:],
    ):
        result = A_graph[(hi, lo)]

        print(
            f"    A_{hi}/A_{lo}: "
            f"exact={result.divides} "
            f"q_degree={result.quotient_degree} "
            f"r_degree={result.remainder_degree}"
        )

        if result.divides:
            print(
                "      quotient="
                f"{sp.factor(result.quotient)}"
            )

    print()
    print("  B CHANNEL")

    for lo, hi in zip(
        B_keys,
        B_keys[1:],
    ):
        result = B_graph[(hi, lo)]

        print(
            f"    B_{hi}/B_{lo}: "
            f"exact={result.divides} "
            f"q_degree={result.quotient_degree} "
            f"r_degree={result.remainder_degree}"
        )

        if result.divides:
            print(
                "      quotient="
                f"{sp.factor(result.quotient)}"
            )

    # --------------------------------------------------------------
    # 4. quotient recognition
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. QUOTIENT RECOGNITION")
    print("=" * 78)

    print()
    print("  A CHANNEL")

    for hi, lo, quotient in exact_A_edges:

        direct = recognize_existing(
            quotient,
            A_even,
        )

        products = recognize_products(
            quotient,
            A_even,
        )

        print(
            f"    A_{hi}/A_{lo}:"
        )
        print(
            f"      quotient="
            f"{sp.factor(quotient)}"
        )
        print(
            f"      scalar*existing={direct}"
        )
        print(
            f"      scalar*product-existing={products}"
        )

    print()
    print("  B CHANNEL")

    for hi, lo, quotient in exact_B_edges:

        direct = recognize_existing(
            quotient,
            B_odd,
        )

        products = recognize_products(
            quotient,
            B_odd,
        )

        print(
            f"    B_{hi}/B_{lo}:"
        )
        print(
            f"      quotient="
            f"{sp.factor(quotient)}"
        )
        print(
            f"      scalar*existing={direct}"
        )
        print(
            f"      scalar*product-existing={products}"
        )

    # --------------------------------------------------------------
    # 5. quotient/remainder profile
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LOW-DEGREE QUOTIENT / REMAINDER PROFILE")
    print("=" * 78)

    print()
    print("  A CHANNEL")

    for (hi, lo), result in sorted(
        A_graph.items(),
        reverse=True,
    ):
        print(
            f"    A_{hi}/A_{lo}: "
            f"q_deg={result.quotient_degree} "
            f"r_deg={result.remainder_degree} "
            f"exact={result.divides}"
        )

    print()
    print("  B CHANNEL")

    for (hi, lo), result in sorted(
        B_graph.items(),
        reverse=True,
    ):
        print(
            f"    B_{hi}/B_{lo}: "
            f"q_deg={result.quotient_degree} "
            f"r_deg={result.remainder_degree} "
            f"exact={result.divides}"
        )

    # --------------------------------------------------------------
    # 6. factorization
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FACTORIZATION OF CHANNEL POLYNOMIALS")
    print("=" * 78)

    print()
    print("  A CHANNEL")

    for d in sorted(
        A_even,
        reverse=True,
    ):
        print(
            f"    A_{d} = "
            f"{sp.factor(A_even[d])}"
        )

    print()
    print("  B CHANNEL")

    for d in sorted(
        B_odd,
        reverse=True,
    ):
        print(
            f"    B_{d} = "
            f"{sp.factor(B_odd[d])}"
        )

    # --------------------------------------------------------------
    # 7. gcd of exact quotients
    # --------------------------------------------------------------

    A_quotients = [
        S(qexpr)
        for _, _, qexpr in exact_A_edges
        if not exact_equal(qexpr, 0)
    ]

    B_quotients = [
        S(qexpr)
        for _, _, qexpr in exact_B_edges
        if not exact_equal(qexpr, 0)
    ]

    print()
    print("=" * 78)
    print("7. GCD OF EXACT QUOTIENTS")
    print("=" * 78)

    if A_quotients:
        print(
            "  gcd(A exact quotients) = "
            f"{sp.factor(gcd_channel(A_quotients))}"
        )
    else:
        print(
            "  gcd(A exact quotients) = none"
        )

    if B_quotients:
        print(
            "  gcd(B exact quotients) = "
            f"{sp.factor(gcd_channel(B_quotients))}"
        )
    else:
        print(
            "  gcd(B exact quotients) = none"
        )

    # --------------------------------------------------------------
    # 8. quotient seed summary
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. QUOTIENT SEED SEARCH")
    print("=" * 78)

    unique_A = []

    for _, _, qexpr in exact_A_edges:
        if not any(
            exact_equal(qexpr, existing)
            for existing in unique_A
        ):
            unique_A.append(
                S(qexpr)
            )

    unique_B = []

    for _, _, qexpr in exact_B_edges:
        if not any(
            exact_equal(qexpr, existing)
            for existing in unique_B
        ):
            unique_B.append(
                S(qexpr)
            )

    print()
    print(
        f"  unique exact A quotients = "
        f"{len(unique_A)}"
    )

    for qexpr in unique_A:
        print(
            f"    {sp.factor(qexpr)}"
        )

    print()
    print(
        f"  unique exact B quotients = "
        f"{len(unique_B)}"
    )

    for qexpr in unique_B:
        print(
            f"    {sp.factor(qexpr)}"
        )

    # --------------------------------------------------------------
    # 9. Fresh semiprime
    # --------------------------------------------------------------

    p = 106621
    q = 246473

    if not sp.isprime(p):
        raise RuntimeError(f"p={p} is not prime")

    if not sp.isprime(q):
        raise RuntimeError(f"q={q} is not prime")

    if p == q:
        raise RuntimeError("Fresh prime pair must be distinct")

    n_val = p * q
    s_val = p + q
    x_val = s_val + 1

    t_val = sp.cancel(
        sp.Rational(n_val, x_val)
    )

    u_val = sp.cancel(
        t_val * (t_val + 1)
    )

    print()
    print("=" * 78)
    print("9. FRESH EXACT SEMIPRIME VALIDATION")
    print("=" * 78)

    print(f"  p={p}")
    print(f"  q={q}")
    print(f"  N={n_val}")
    print(f"  S={s_val}")
    print(f"  X={x_val}")
    print(f"  t={t_val}")
    print(f"  u={u_val}")

    fresh_channel_ok = True

    for d in range(16, -1, -1):

        h = S(H_T[d])

        observed = sp.expand(
            h.subs(t, t_val)
        )

        info = infos[d]

        predicted = sp.expand(
            info.A.subs(u, u_val)
            +
            (2*t_val + 1)
            *
            info.B.subs(u, u_val)
        )

        ok = exact_equal(
            observed,
            predicted,
        )

        fresh_channel_ok &= ok

        print(
            f"  d={d:2d}: "
            f"channel reconstruction={ok}"
        )

    print()
    print(
        f"  ALL FRESH CHANNEL CHECKS = "
        f"{fresh_channel_ok}"
    )

    # --------------------------------------------------------------
    # 10. Fresh actual layers
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FRESH LAYER RECONSTRUCTION FROM A/B CHANNELS")
    print("=" * 78)

    fresh_layer_ok = True

    for d in range(16, -1, -1):

        direct = layer_value(
            d,
            n_val,
            x_val,
        )

        info = infos[d]

        channel_value = sp.expand(
            info.A.subs(u, u_val)
            +
            (2*t_val + 1)
            *
            info.B.subs(u, u_val)
        )

        reconstructed = sp.expand(
            sp.Integer(x_val)**d
            * channel_value
        )

        ok = exact_equal(
            direct,
            reconstructed,
        )

        fresh_layer_ok &= ok

        print(
            f"  L_{d:2d}: "
            f"direct=channel={ok}"
        )

    print()
    print(
        f"  ALL LAYER RECONSTRUCTIONS = "
        f"{fresh_layer_ok}"
    )

    # --------------------------------------------------------------
    # 11. Construction interpretation
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. CONSTRUCTION INTERPRETATION")
    print("=" * 78)

    print(
        """
  The layer family is split exactly into two channels:

      even d:
          h_d(t) = A_d(u)

      odd d:
          h_d(t) = (2t+1) B_d(u),

      u=t(t+1).

  The present experiment therefore asks whether consecutive layers
  are multiplicatively nested inside their own channel.

  An exact relation

      A_hi = Q(u) A_lo

  or

      B_hi = Q(u) B_lo

  would identify an exact construction operator.

  Quotient recognition then tests whether Q(u) is itself already
  present among the channel polynomials, or is a product of existing
  channel elements.

  This produces a directed algebraic construction graph.

  A small graph would be strong evidence that the layer family is
  generated from a small seed set.

  Conversely, if almost no exact divisibility occurs, the channel
  decomposition is useful for understanding symmetry but does not
  by itself reveal a multiplicative construction law.

  This experiment intentionally does NOT infer a factorization
  algorithm.  It only studies how the already-known layer family
  is constructed.
        """
    )

    # --------------------------------------------------------------
    # 12. Exactness
    # --------------------------------------------------------------

    failures = 0

    if not all_involution_ok:
        failures += 1

    if not fresh_channel_ok:
        failures += 1

    if not fresh_layer_ok:
        failures += 1

    print()
    print("=" * 78)
    print("FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  involution_decomposition = "
        f"{all_involution_ok}"
    )

    print(
        f"  fresh_channel_validation = "
        f"{fresh_channel_ok}"
    )

    print(
        f"  fresh_layer_validation = "
        f"{fresh_layer_ok}"
    )

    print(
        f"  exact_A_divisions = "
        f"{len(exact_A_edges)}"
    )

    print(
        f"  exact_B_divisions = "
        f"{len(exact_B_edges)}"
    )

    print(
        f"  failures = {failures}"
    )

    print(
        f"  ALL BASIC CHECKS PASS = "
        f"{failures == 0}"
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 75 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()