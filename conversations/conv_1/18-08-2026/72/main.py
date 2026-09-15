import sympy as sp
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# =============================================================================
# EXPERIMENT 74
#
# EXACT LAYER-CONSTRUCTION / CHANNEL-OPERATOR SEARCH
#
# We know:
#
#     L_d(N,X) = X^d h_d(t),     t=N/X
#
# and
#
#     h_d(t) = A_d(u) + (2t+1) B_d(u),
#     u = t(t+1).
#
# This experiment asks:
#
#     HOW ARE THE A_d AND B_d BUILT AS d CHANGES?
#
# We test:
#
#   1. exact involution reconstruction;
#   2. polynomial-coefficient recurrences in u;
#   3. neighboring-layer operator forms;
#   4. simple-factor normalized recurrences;
#   5. constant-rational relations;
#   6. fresh exact validation.
#
# No floating point.
# No regression.
# No numerical fitting.
# Everything exact over QQ / ZZ.
# =============================================================================


# =============================================================================
# SYMBOLS
# =============================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")


# =============================================================================
# NORMALIZED LAYER POLYNOMIALS h_d(t)
# =============================================================================

H_RAW = {
    16: -(3*t**2 + 3*t + 1) * (
        3*t**6 + 9*t**5 + 18*t**4 + 21*t**3 +
        15*t**2 + 6*t + 1
    ),

    15: (2*t + 1) * (
        44*t**8 + 176*t**7 + 494*t**6 + 866*t**5 +
        1016*t**4 + 794*t**3 + 401*t**2 + 119*t + 16
    ),

    14: -(
        276*t**10 + 1380*t**9 + 5460*t**8 + 13560*t**7 +
        23058*t**6 + 27510*t**5 + 23100*t**4 +
        13410*t**3 + 5135*t**2 + 1169*t + 120
    ),

    13: (2*t + 1) * (
        144*t**10 + 720*t**9 + 4810*t**8 + 14920*t**7 +
        32272*t**6 + 47620*t**5 + 48955*t**4 +
        34510*t**3 + 16020*t**2 + 4431*t + 560
    ),

    12: -(
        50*t**12 + 300*t**11 + 7436*t**10 + 34430*t**9 +
        117810*t**8 + 267960*t**7 + 426888*t**6 +
        484110*t**5 + 390225*t**4 + 219010*t**3 +
        81510*t**2 + 18109*t + 1820
    ),

    11: 13*(2*t + 1) * (
        50*t**10 + 250*t**9 + 2284*t**8 + 7636*t**7 +
        17434*t**6 + 26626*t**5 + 28036*t**4 +
        20104*t**3 + 9451*t**2 + 2639*t + 336
    ),

    10: -1001*(
        10*t**10 + 50*t**9 + 252*t**8 + 708*t**7 +
        1302*t**6 + 1638*t**5 + 1428*t**4 +
        852*t**3 + 333*t**2 + 77*t + 8
    ),

    9: (2*t + 1) * (
        17875*t**8 + 71500*t**7 + 245960*t**6 +
        487630*t**5 + 622414*t**4 + 515528*t**3 +
        271506*t**2 + 83097*t + 11441
    ),

    8: -(
        71500*t**8 + 286000*t**7 + 755664*t**6 +
        1265992*t**5 + 1375360*t**4 + 974400*t**3 +
        436832*t**2 + 112964*t + 12879
    ),

    7: 4*(2*t + 1) * (
        11050*t**6 + 33150*t**5 + 63036*t**4 +
        70822*t**3 + 47549*t**2 + 17663*t + 2869
    ),

    6: -4*(
        17850*t**6 + 53550*t**5 + 86088*t**4 +
        82926*t**3 + 47574*t**2 + 15036*t + 2023
    ),

    5: 6*(2*t + 1) * (
        3230*t**4 + 6460*t**3 + 6688*t**2 +
        3458*t + 749
    ),

    4: -2*(
        7125*t**4 + 14250*t**3 + 12690*t**2 +
        5565*t + 973
    ),

    3: 14*(2*t + 1) * (
        125*t**2 + 125*t + 46
    ),

    2: -2*(275*t**2 + 275*t + 78),

    1: 25*(2*t + 1),

    0: -2,
}

H: Dict[int, sp.Expr] = {
    d: sp.expand(sp.sympify(expr))
    for d, expr in H_RAW.items()
}


# =============================================================================
# EXACT POLYNOMIAL HELPERS
# =============================================================================

def simplify_zero(expr) -> bool:
    return sp.expand(sp.cancel(sp.sympify(expr))) == 0


def poly_t(expr: sp.Expr) -> sp.Poly:
    return sp.Poly(
        sp.expand(sp.sympify(expr)),
        t,
        domain=sp.QQ,
    )


def poly_u(expr: sp.Expr) -> sp.Poly:
    return sp.Poly(
        sp.expand(sp.sympify(expr)),
        u,
        domain=sp.QQ,
    )


def degree_u(expr):
    p = poly_u(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def primitive_signature(
    expr: sp.Expr,
    variable: sp.Symbol,
) -> Tuple[int, List[int]]:

    p = sp.Poly(
        sp.expand(sp.sympify(expr)),
        variable,
        domain=sp.QQ,
    )

    if p.is_zero:
        return 0, [0]

    coeffs = p.all_coeffs()

    den_lcm = 1
    for c in coeffs:
        den_lcm = int(
            sp.ilcm(
                den_lcm,
                int(sp.denom(c))
            )
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    g = 0
    for z in ints:
        g = sp.igcd(g, abs(z))

    if g == 0:
        return 0, [0]

    ints = [z // g for z in ints]

    if ints[0] < 0:
        ints = [-z for z in ints]

    return int(g), ints


# =============================================================================
# CANONICAL INVOLUTION DECOMPOSITION
# =============================================================================

def canonical_involution_decomposition(
    h_expr: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr]:

    h = sp.expand(
        sp.sympify(h_expr)
    )

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

    # -------------------------------------------------------------------------
    # A(u)
    # -------------------------------------------------------------------------

    y = sp.Symbol("y")

    even_y = sp.expand(
        even.subs(
            t,
            (y - 1)/2
        )
    )

    py = sp.Poly(
        even_y,
        y,
        domain=sp.QQ,
    )

    A = sp.Integer(0)

    for (power,), coeff in py.terms():

        if power % 2:
            raise RuntimeError(
                "Even channel unexpectedly contains odd powers of y."
            )

        A += coeff * (4*u + 1)**(power // 2)

    A = sp.expand(A)

    # -------------------------------------------------------------------------
    # B(u)
    # -------------------------------------------------------------------------

    divisor = sp.Poly(
        2*t + 1,
        t,
        domain=sp.QQ,
    )

    p_odd = sp.Poly(
        odd,
        t,
        domain=sp.QQ,
    )

    q, r = sp.div(
        p_odd,
        divisor,
        domain=sp.QQ,
    )

    if not r.is_zero:
        raise RuntimeError(
            "Odd channel is not divisible by 2t+1."
        )

    B_t = sp.expand(
        q.as_expr()
    )

    B_y = sp.expand(
        B_t.subs(
            t,
            (y - 1)/2
        )
    )

    pB = sp.Poly(
        B_y,
        y,
        domain=sp.QQ,
    )

    B = sp.Integer(0)

    for (power,), coeff in pB.terms():

        if power % 2:
            raise RuntimeError(
                "B channel unexpectedly contains odd powers of y."
            )

        B += coeff * (4*u + 1)**(power // 2)

    B = sp.expand(B)

    return A, B


A: Dict[int, sp.Expr] = {}
B: Dict[int, sp.Expr] = {}

for d in range(17):
    A[d], B[d] = canonical_involution_decomposition(H[d])


# =============================================================================
# ROBUST RECONSTRUCTION CHECK
# =============================================================================

def reconstruct_h(
    d: int,
) -> sp.Expr:

    return sp.expand(
        A[d].subs(
            u,
            t*(t + 1)
        )
        +
        (2*t + 1) *
        B[d].subs(
            u,
            t*(t + 1)
        )
    )


def reconstruction_check(d: int) -> bool:

    lhs = sp.expand(
        H[d]
    )

    rhs = reconstruct_h(d)

    return simplify_zero(
        lhs - rhs
    )


# =============================================================================
# POLYNOMIAL SEQUENCE UTILITIES
# =============================================================================

def sequence_gcd(
    seq: List[sp.Expr],
) -> sp.Expr:

    nonzero = [
        poly_u(x)
        for x in seq
        if not poly_u(x).is_zero
    ]

    if not nonzero:
        return sp.Integer(0)

    g = nonzero[0]

    for p in nonzero[1:]:
        g = sp.gcd(g, p)

    return sp.factor(
        g.as_expr()
    )


def constant_rational_relations(
    seq: List[sp.Expr],
) -> Tuple[int, int, List[List[sp.Expr]]]:

    polys = [
        poly_u(x)
        for x in seq
    ]

    max_degree = max(
        [
            -1 if p.is_zero else p.degree()
            for p in polys
        ]
    )

    if max_degree < 0:
        return 0, len(seq), []

    rows = []

    for power in range(
        max_degree + 1
    ):

        rows.append(
            [
                p.coeff_monomial(u**power)
                for p in polys
            ]
        )

    M = sp.Matrix(rows)

    rank = M.rank()

    nullspace = M.nullspace()

    return (
        rank,
        len(nullspace),
        [
            list(v)
            for v in nullspace
        ],
    )


# =============================================================================
# GENERAL EXACT RECURRENCE SEARCH
# =============================================================================

@dataclass
class RecurrenceCandidate:
    order: int
    degree: int
    coefficients: List[sp.Expr]
    rank: int
    unknowns: int


def recurrence_system(
    sequence: List[sp.Expr],
    order: int,
    coefficient_degree: int,
):

    unknowns = []

    polys = []

    for j in range(
        order + 1
    ):

        c = sp.Integer(0)

        for k in range(
            coefficient_degree + 1
        ):

            z = sp.Symbol(
                f"c_{j}_{k}"
            )

            unknowns.append(z)

            c += z*u**k

        polys.append(
            sp.expand(c)
        )

    equations = []

    for r in range(
        order,
        len(sequence)
    ):

        expr = sp.Integer(0)

        for j in range(
            order + 1
        ):

            expr += (
                polys[j] *
                sequence[r - j]
            )

        expr = sp.expand(expr)

        pu = sp.Poly(
            expr,
            u,
            domain=sp.EX,
        )

        for (_,), coeff in pu.terms():

            coeff = sp.expand(coeff)

            for variable in unknowns:

                equations.append(
                    coeff.coeff(variable)
                )

            # remove the linear part and verify nothing nonlinear survives
            linear_part = sp.Integer(0)

            for variable in unknowns:

                linear_part += (
                    variable *
                    sp.expand(
                        coeff.coeff(variable)
                    )
                )

            nonlinear = sp.expand(
                coeff -
                linear_part
            )

            if nonlinear != 0:
                raise RuntimeError(
                    "Recurrence system became nonlinear."
                )

    M = sp.Matrix(
        sp.linear_eq_to_matrix(
            equations,
            unknowns,
        )[0]
    )

    return (
        unknowns,
        M,
    )


def recurrence_search(
    sequence: List[sp.Expr],
    orders=(1, 2, 3, 4),
    degrees=(0, 1, 2, 3, 4),
) -> List[RecurrenceCandidate]:

    found = []

    for order in orders:

        for degree in degrees:

            unknowns, M = recurrence_system(
                sequence,
                order,
                degree,
            )

            rank = M.rank()
            nullity = len(unknowns) - rank

            if nullity <= 0:
                continue

            vectors = M.nullspace()

            for vector in vectors:

                candidate = []

                idx = 0

                for j in range(
                    order + 1
                ):

                    c = sp.Integer(0)

                    for k in range(
                        degree + 1
                    ):

                        c += (
                            vector[idx] *
                            u**k
                        )

                        idx += 1

                    candidate.append(
                        sp.expand(c)
                    )

                # Verify exact recurrence.
                valid = True

                for r in range(
                    order,
                    len(sequence)
                ):

                    expr = sp.Integer(0)

                    for j in range(
                        order + 1
                    ):

                        expr += (
                            candidate[j] *
                            sequence[r-j]
                        )

                    if not simplify_zero(expr):
                        valid = False
                        break

                if valid:
                    found.append(
                        RecurrenceCandidate(
                            order,
                            degree,
                            candidate,
                            rank,
                            len(unknowns),
                        )
                    )

    return found


# =============================================================================
# OPERATOR SEARCH
# =============================================================================

def operator_search(
    name: str,
    sequence: List[sp.Expr],
):

    print()
    print("=" * 78)
    print(
        f"{name} SIMPLE NEIGHBORING-LAYER OPERATOR SEARCH"
    )
    print("=" * 78)

    factors = [
        sp.Integer(1),
        u,
        u + 1,
        4*u + 1,
        2*u + 1,
        3*u + 1,
        3*u + 3,
        5*u + 1,
    ]

    for f in factors:

        transformed = [
            sp.expand(
                f * x
            )
            for x in sequence
        ]

        g = sequence_gcd(
            transformed
        )

        print(
            f"  multiplier={sp.factor(f)} "
            f"common_gcd={sp.factor(g)}"
        )

    print()
    print("  Neighbor quotient/remainder:")

    for i in range(
        len(sequence) - 1
    ):

        hi = poly_u(
            sequence[i]
        )

        lo = poly_u(
            sequence[i + 1]
        )

        if lo.is_zero:
            print(
                f"    step {i}: denominator zero"
            )
            continue

        q, r = sp.div(
            hi,
            lo,
            domain=sp.QQ,
        )

        print(
            f"    step {i}: "
            f"q_deg={(-sp.oo if q.is_zero else q.degree())}, "
            f"r_deg={(-sp.oo if r.is_zero else r.degree())}"
        )


# =============================================================================
# NORMALIZED-DIVISION SEARCH
# =============================================================================

def normalized_division_search(
    name: str,
    sequence: List[sp.Expr],
):

    print()
    print("=" * 78)
    print(
        f"{name} EXACT NORMALIZED-DIVISION SEARCH"
    )
    print("=" * 78)

    for shift in range(1, 5):

        if len(sequence) <= shift:
            continue

        print(
            f"  shift={shift}"
        )

        for i in range(
            len(sequence) - shift
        ):

            a = poly_u(
                sequence[i]
            )

            b = poly_u(
                sequence[i + shift]
            )

            if b.is_zero:
                continue

            q, r = sp.div(
                a,
                b,
                domain=sp.QQ,
            )

            if r.is_zero:

                print(
                    f"    EXACT DIVISION "
                    f"step={i}->{i+shift}: "
                    f"quotient={sp.factor(q.as_expr())}"
                )


# =============================================================================
# CHANNEL PROFILE
# =============================================================================

def profile(
    name: str,
    degrees: List[int],
    channel: Dict[int, sp.Expr],
):

    print()
    print("-" * 78)
    print(
        f"{name} CHANNEL"
    )
    print("-" * 78)

    for d in degrees:

        content, primitive = (
            primitive_signature(
                channel[d],
                u
            )
        )

        print(
            f"  {name}_{d}: "
            f"degree={degree_u(channel[d])} "
            f"terms={len(poly_u(channel[d]).terms())} "
            f"content={content}"
        )

        print(
            f"    {sp.factor(channel[d])}"
        )

        print(
            f"    primitive={primitive}"
        )


# =============================================================================
# FRESH SEMIPRIME CHECK
# =============================================================================

def fresh_check() -> bool:

    p = sp.Integer(106621)
    q = sp.Integer(246473)

    N = p*q
    X = p + q + 1

    tv = sp.cancel(
        N / X
    )

    uv = sp.cancel(
        tv*(tv + 1)
    )

    ok_all = True

    print()
    print("=" * 78)
    print("FRESH EXACT SEMIPRIME CHECK")
    print("=" * 78)

    print(
        f"  p={p}"
    )
    print(
        f"  q={q}"
    )
    print(
        f"  N={N}"
    )
    print(
        f"  X={X}"
    )
    print(
        f"  t={tv}"
    )
    print(
        f"  u={uv}"
    )

    for d in range(
        16,
        -1,
        -1
    ):

        observed = sp.expand(
            H[d].subs(
                t,
                tv
            )
        )

        reconstructed = sp.expand(
            A[d].subs(
                u,
                uv
            )
            +
            (2*tv + 1) *
            B[d].subs(
                u,
                uv
            )
        )

        ok = simplify_zero(
            observed - reconstructed
        )

        ok_all = ok_all and ok

        print(
            f"  d={d:2d}: "
            f"reconstruction={ok}"
        )

    return ok_all


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 74 — EXACT LAYER-CONSTRUCTION / "
        "CHANNEL-OPERATOR SEARCH"
    )
    print("=" * 78)

    print()
    print("0. EXACT SYMBOLIC SETUP")
    print("  t = N/X")
    print("  u = t(t+1)")
    print("  involution = t -> -1-t")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # -------------------------------------------------------------------------
    # 1. INVOLUTION VALIDATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("=" * 78)

    involution_ok = True

    for d in range(
        16,
        -1,
        -1
    ):

        h = H[d]

        partner = sp.expand(
            h.subs(
                t,
                -1 - t
            )
        )

        invariant = simplify_zero(
            partner - h
        )

        anti = simplify_zero(
            partner + h
        )

        reconstructed = reconstruct_h(
            d
        )

        ok = simplify_zero(
            reconstructed - h
        )

        if not ok:
            involution_ok = False

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"reconstruction={ok}"
        )

    print(
        f"\n  ALL INVOLUTION CHECKS = "
        f"{involution_ok}"
    )

    # -------------------------------------------------------------------------
    # 2. CHANNEL PROFILES
    # -------------------------------------------------------------------------

    A_degrees = [
        16, 14, 12, 10,
        8, 6, 4, 2, 0
    ]

    B_degrees = [
        15, 13, 11, 9,
        7, 5, 3, 1
    ]

    profile(
        "A",
        A_degrees,
        A,
    )

    profile(
        "B",
        B_degrees,
        B,
    )

    # -------------------------------------------------------------------------
    # 3. CHANNEL GCDS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CHANNEL GCD STRUCTURE")
    print("=" * 78)

    gA = sequence_gcd(
        [
            A[d]
            for d in A_degrees
        ]
    )

    gB = sequence_gcd(
        [
            B[d]
            for d in B_degrees
        ]
    )

    print(
        f"  gcd(A_16,A_14,...,A_0) = "
        f"{sp.factor(gA)}"
    )

    print(
        f"  gcd(B_15,B_13,...,B_1) = "
        f"{sp.factor(gB)}"
    )

    # -------------------------------------------------------------------------
    # 4. CONSTANT-RATIONAL RELATIONS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CONSTANT-RATIONAL CHANNEL RELATIONS")
    print("=" * 78)

    A_rank, A_nullity, A_relations = (
        constant_rational_relations(
            [
                A[d]
                for d in A_degrees
            ]
        )
    )

    B_rank, B_nullity, B_relations = (
        constant_rational_relations(
            [
                B[d]
                for d in B_degrees
            ]
        )
    )

    print(
        f"  A rank={A_rank} "
        f"nullity={A_nullity}"
    )

    for idx, rel in enumerate(
        A_relations,
        1
    ):

        expr = sp.Integer(0)

        for coeff, d in zip(
            rel,
            A_degrees
        ):

            expr += (
                coeff *
                A[d]
            )

        print(
            f"  A relation #{idx}: "
            f"zero={simplify_zero(expr)}"
        )
        print(
            f"    {sp.factor(expr)}"
        )

    print(
        f"\n  B rank={B_rank} "
        f"nullity={B_nullity}"
    )

    for idx, rel in enumerate(
        B_relations,
        1
    ):

        expr = sp.Integer(0)

        for coeff, d in zip(
            rel,
            B_degrees
        ):

            expr += (
                coeff *
                B[d]
            )

        print(
            f"  B relation #{idx}: "
            f"zero={simplify_zero(expr)}"
        )

    # -------------------------------------------------------------------------
    # 5. RECURRENCE SEARCH
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. POLYNOMIAL-COEFFICIENT RECURRENCE SEARCH")
    print("=" * 78)

    for name, degrees, channel in (
        ("A", A_degrees, A),
        ("B", B_degrees, B),
    ):

        sequence = [
            channel[d]
            for d in degrees
        ]

        print()
        print(
            f"  {name}-channel:"
        )

        found = recurrence_search(
            sequence,
            orders=(1, 2, 3, 4),
            degrees=(0, 1, 2, 3, 4),
        )

        if not found:

            print(
                "    No exact recurrence in scanned range."
            )

        else:

            for item in found:

                print(
                    f"    FOUND order={item.order} "
                    f"degree={item.degree} "
                    f"rank={item.rank} "
                    f"unknowns={item.unknowns}"
                )

                for j, coeff in enumerate(
                    item.coefficients
                ):

                    print(
                        f"      c_{j}(u) = "
                        f"{sp.factor(coeff)}"
                    )

    # -------------------------------------------------------------------------
    # 6. SIMPLE OPERATOR SEARCH
    # -------------------------------------------------------------------------

    operator_search(
        "A",
        [
            A[d]
            for d in A_degrees
        ]
    )

    operator_search(
        "B",
        [
            B[d]
            for d in B_degrees
        ]
    )

    # -------------------------------------------------------------------------
    # 7. NORMALIZED DIVISION SEARCH
    # -------------------------------------------------------------------------

    normalized_division_search(
        "A",
        [
            A[d]
            for d in A_degrees
        ]
    )

    normalized_division_search(
        "B",
        [
            B[d]
            for d in B_degrees
        ]
    )

    # -------------------------------------------------------------------------
    # 8. CROSS-CHANNEL CORRESPONDENCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. CROSS-CHANNEL NEIGHBOR TEST")
    print("=" * 78)

    pairs = [
        (16, 15),
        (14, 13),
        (12, 11),
        (10, 9),
        (8, 7),
        (6, 5),
        (4, 3),
        (2, 1),
    ]

    for dA, dB in pairs:

        a = A[dA]
        b = B[dB]

        print(
            f"\n  pair A_{dA}, B_{dB}"
        )

        ga = poly_u(a)
        gb = poly_u(b)

        g = sp.gcd(
            ga,
            gb
        )

        print(
            f"    gcd={sp.factor(g.as_expr())}"
        )

        # Test simple linear relation:
        # A_d = (alpha*u + beta) B_d
        if not gb.is_zero:

            q, r = sp.div(
                ga,
                gb,
                domain=sp.QQ
            )

            if r.is_zero:
                print(
                    f"    exact quotient="
                    f"{sp.factor(q.as_expr())}"
                )
            else:
                print(
                    f"    quotient_degree="
                    f"{(-sp.oo if q.is_zero else q.degree())}"
                )
                print(
                    f"    remainder_degree="
                    f"{(-sp.oo if r.is_zero else r.degree())}"
                )

    # -------------------------------------------------------------------------
    # 9. FRESH VALIDATION
    # -------------------------------------------------------------------------

    fresh_ok = fresh_check()

    # -------------------------------------------------------------------------
    # 10. FINAL INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The exact structural identity remains

    L_d(N,X) = X^d h_d(N/X)

with

    h_d(t)
      = A_d(u)
      + (2t+1)B_d(u),

    u=t(t+1).

The repaired reconstruction test now compares symbolic expressions
through exact cancellation rather than raw Python/SymPy expression
identity.

The main new question is construction rather than inversion:

    Are the sequences

        A_16,A_14,...,A_0

    and

        B_15,B_13,...,B_1

    generated by a small exact operator in d?

There are several possible outcomes.

1. A short polynomial-coefficient recurrence exists.

   Then the whole layer family may be recursively generated from a
   small number of seed polynomials.

2. No short recurrence exists, but exact quotient/remainder patterns
   stabilize.

   Then the construction may involve a variable-coefficient operator.

3. Constant-rational relations exist without a local recurrence.

   Then the family may lie in a low-dimensional algebraic subspace
   without being generated sequentially.

4. No compact structure appears.

   Then the layer family may be naturally generated by the original
   kernel expansion rather than a short recurrence.

The distinction matters for the factorization question.

The inverse experiment showed

    top layers -> t -> X -> N -> S -> p,q.

The present experiment asks the complementary question:

    kernel -> layer family.

Only after understanding that construction mechanism can we seriously
ask whether the required layer information could itself be obtained
from N without knowing p and q.
"""
    )

    print()
    print("=" * 78)
    print("FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  involution_decomposition = "
        f"{involution_ok}"
    )

    print(
        f"  fresh_validation = "
        f"{fresh_ok}"
    )

    print(
        f"  all_basic_checks = "
        f"{involution_ok and fresh_ok}"
    )

    print("=" * 78)
    print(
        "EXPERIMENT 74 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

