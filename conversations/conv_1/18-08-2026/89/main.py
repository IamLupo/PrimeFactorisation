#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 91 — EXACT KNOWN-LAYER UNIVERSALITY / STRUCTURE AUDIT
# ==============================================================================

N, X, t, u, j = sp.symbols("N X t u j")


# ==============================================================================
# 1. EXACT BASIC HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr, var):
    return sp.Poly(
        clean(expr),
        var,
        domain=sp.QQ,
    )


def degree(expr, var):
    p = poly(expr, var)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def coeff(expr, var, power):
    return clean(
        poly(expr, var).coeff_monomial(var ** power)
    )


def falling(x, n):
    out = sp.Integer(1)
    for r in range(n):
        out *= x - r
    return clean(out)


# ==============================================================================
# 2. EXACT h_d(t) DATA FROM THE ESTABLISHED (9,16) KERNEL
# ==============================================================================

H = {
    16: -(3*t**2 + 3*t + 1) * (
        3*t**6 + 9*t**5 + 18*t**4 +
        21*t**3 + 15*t**2 + 6*t + 1
    ),

    15: (2*t + 1) * (
        44*t**8 + 176*t**7 + 494*t**6 +
        866*t**5 + 1016*t**4 + 794*t**3 +
        401*t**2 + 119*t + 16
    ),

    14: -(
        276*t**10 + 1380*t**9 + 5460*t**8 +
        13560*t**7 + 23058*t**6 + 27510*t**5 +
        23100*t**4 + 13410*t**3 + 5135*t**2 +
        1169*t + 120
    ),

    13: (2*t + 1) * (
        144*t**10 + 720*t**9 + 4810*t**8 +
        14920*t**7 + 32272*t**6 + 47620*t**5 +
        48955*t**4 + 34510*t**3 + 16020*t**2 +
        4431*t + 560
    ),

    12: -(
        50*t**12 + 300*t**11 + 7436*t**10 +
        34430*t**9 + 117810*t**8 + 267960*t**7 +
        426888*t**6 + 484110*t**5 + 390225*t**4 +
        219010*t**3 + 81510*t**2 + 18109*t + 1820
    ),

    11: 13*(2*t + 1) * (
        50*t**10 + 250*t**9 + 2284*t**8 +
        7636*t**7 + 17434*t**6 + 26626*t**5 +
        28036*t**4 + 20104*t**3 + 9451*t**2 +
        2639*t + 336
    ),

    10: -1001 * (
        10*t**10 + 50*t**9 + 252*t**8 +
        708*t**7 + 1302*t**6 + 1638*t**5 +
        1428*t**4 + 852*t**3 + 333*t**2 +
        77*t + 8
    ),

    9: (2*t + 1) * (
        17875*t**8 + 71500*t**7 + 245960*t**6 +
        487630*t**5 + 622414*t**4 + 515528*t**3 +
        271506*t**2 + 83097*t + 11441
    ),

    8: -(
        71500*t**8 + 286000*t**7 + 755664*t**6 +
        1265992*t**5 + 1375360*t**4 +
        974400*t**3 + 436832*t**2 +
        112964*t + 12879
    ),

    7: 4*(2*t + 1) * (
        11050*t**6 + 33150*t**5 + 63036*t**4 +
        70822*t**3 + 47549*t**2 + 17663*t + 2869
    ),

    6: -4 * (
        17850*t**6 + 53550*t**5 + 86088*t**4 +
        82926*t**3 + 47574*t**2 + 15036*t + 2023
    ),

    5: 6*(2*t + 1) * (
        3230*t**4 + 6460*t**3 + 6688*t**2 +
        3458*t + 749
    ),

    4: -2 * (
        7125*t**4 + 14250*t**3 + 12690*t**2 +
        5565*t + 973
    ),

    3: 14*(2*t + 1) * (
        125*t**2 + 125*t + 46
    ),

    2: -2 * (
        275*t**2 + 275*t + 78
    ),

    1: 25*(2*t + 1),

    0: -2,
}


# ==============================================================================
# 3. EXACT INVOLUTION AUDIT
# ==============================================================================

def involution_audit():
    failures = 0

    print("=" * 78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("=" * 78)

    for d in range(16, -1, -1):

        h = clean(H[d])

        partner = clean(
            h.subs(t, -1 - t)
        )

        invariant = clean(
            partner - h
        ) == 0

        anti = clean(
            partner + h
        ) == 0

        if d % 2 == 0:
            expected = invariant
        else:
            expected = anti

        reconstruction = True

        if d % 2 == 1:
            q, r = sp.div(
                poly(h, t),
                poly(2*t + 1, t),
            )

            reconstruction = (
                r.is_zero
                and clean(
                    (2*t + 1)*q.as_expr() - h
                ) == 0
            )

        ok = expected and reconstruction

        if not ok:
            failures += 1

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"reconstruction={reconstruction}"
        )

    print()
    print(
        "  ALL INVOLUTION CHECKS =",
        failures == 0,
    )

    return failures == 0


# ==============================================================================
# 4. RECOVER A_d(u), B_d(u)
# ==============================================================================

def recover_u(expr):
    """
    Exact reconstruction through u=t(t+1).
    """

    expr = clean(expr)

    deg_t = degree(expr, t)

    if deg_t <= 0:
        return clean(expr)

    deg_u = int(deg_t // 2)

    data = []

    for tv in range(deg_u + 1):

        uv = tv * (tv + 1)

        value = clean(
            expr.subs(t, tv)
        )

        data.append(
            (sp.Integer(uv), value)
        )

    candidate = clean(
        sp.interpolate(
            data,
            u,
        )
    )

    check = clean(
        candidate.subs(
            u,
            t*(t + 1),
        )
        - expr
    )

    if check != 0:
        raise ArithmeticError(
            "u reconstruction failed"
        )

    return candidate


def build_channels():

    A = {}
    B = {}

    for d in range(16, -1, -1):

        h = clean(H[d])

        if d % 2 == 0:

            A[d] = recover_u(h)

        else:

            q, r = sp.div(
                poly(h, t),
                poly(2*t + 1, t),
            )

            if not r.is_zero:
                raise ArithmeticError(
                    f"B channel divisibility failed at d={d}"
                )

            B[d] = recover_u(
                q.as_expr()
            )

    return A, B


# ==============================================================================
# 5. INDEXED CHANNELS
# ==============================================================================

def indexed_channels(A, B):

    Aseq = []
    Bseq = []

    for jj in range(9):
        Aseq.append(
            clean(A[2*jj])
        )

    for jj in range(8):
        Bseq.append(
            clean(B[2*jj + 1])
        )

    return Aseq, Bseq


# ==============================================================================
# 6. EXACT COEFFICIENT POLYNOMIAL P_k(j)
# ==============================================================================

def coefficient_index_polynomial(
    channel,
    k_power,
):

    points = []

    for jj, expr in enumerate(channel):

        value = coeff(
            expr,
            u,
            k_power,
        )

        points.append(
            (
                sp.Integer(jj),
                value,
            )
        )

    return clean(
        sp.interpolate(
            points,
            j,
        )
    )


# ==============================================================================
# 7. EXACT FALLING-FACTORIAL TEST
# ==============================================================================

def falling_test(
    P,
    k_power,
):

    divisor = falling(
        j,
        k_power,
    )

    q, r = sp.div(
        poly(P, j),
        poly(divisor, j),
    )

    return (
        r.is_zero,
        clean(q.as_expr()),
    )


# ==============================================================================
# 8. DIAGONAL C[k,k+s]
# ==============================================================================

def falling_basis_coefficients(P):
    """
    Convert polynomial P(j) to falling-factorial basis.

        P(j) = sum_r c_r j_(r)

    Exact finite difference identity:

        c_r = Delta^r P(0) / r!
    """

    values = [
        clean(P.subs(j, jj))
        for jj in range(
            degree(P, j) + 1
        )
    ]

    out = {}

    current = values[:]

    r = 0

    while current:

        out[r] = clean(
            current[0]
            / sp.factorial(r)
        )

        if len(current) == 1:
            break

        current = [
            clean(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        r += 1

    return out


# ==============================================================================
# 9. EXACT DIAGONAL TABLE
# ==============================================================================

def diagonal_table(channel, max_k):

    C = {}

    for k_power in range(max_k + 1):

        P = coefficient_index_polynomial(
            channel,
            k_power,
        )

        basis = falling_basis_coefficients(
            P
        )

        C[k_power] = basis

    return C


# ==============================================================================
# 10. EXACT RATIO TEST
# ==============================================================================

def exact_ratio(a, b):

    if clean(b) == 0:
        return None

    return clean(
        a / b
    )


def rational_ratio_fit(values):

    """
    Exact low-degree rational search.

    This is deliberately tiny:
      numerator degree <= 2
      denominator degree <= 2

    We reject unless the same rational function satisfies every
    available transition exactly.
    """

    x = sp.symbols("x")

    n0, n1, n2 = sp.symbols(
        "n0 n1 n2"
    )

    d0, d1, d2 = sp.symbols(
        "d0 d1 d2"
    )

    unknowns = [
        n0, n1, n2,
        d0, d1, d2,
    ]

    equations = []

    for idx, (xx, yy) in enumerate(
        values
    ):

        lhs = (
            n0
            + n1*xx
            + n2*xx**2
        )

        rhs = yy * (
            d0
            + d1*xx
            + d2*xx**2
        )

        equations.append(
            clean(lhs - rhs)
        )

    M, vec = sp.linear_eq_to_matrix(
        equations,
        unknowns,
    )

    nullspace = M.nullspace()

    if not nullspace:
        return None

    # A one-dimensional nullspace is the most meaningful case.
    if len(nullspace) != 1:
        return None

    solution = nullspace[0]

    den = clean(
        solution[3]
        + solution[4]*x
        + solution[5]*x**2
    )

    num = clean(
        solution[0]
        + solution[1]*x
        + solution[2]*x**2
    )

    if den == 0:
        return None

    return clean(
        num / den
    )


# ==============================================================================
# 11. MAIN AUDIT
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 91 — EXACT FALLING-FACTORIAL "
        "STRUCTURE STABILITY AUDIT"
    )
    print("=" * 78)

    print()
    print("0. EXACT SETUP")
    print("  kernel instance = G_{9,16}")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    involution_ok = involution_audit()

    A, B = build_channels()

    Aseq, Bseq = indexed_channels(
        A,
        B,
    )

    print()
    print("=" * 78)
    print("2. CHANNEL DEGREE PROFILE")
    print("=" * 78)

    for jj, expr in enumerate(Aseq):
        print(
            f"  A_{jj}: "
            f"degree_u={degree(expr,u)}"
        )

    for jj, expr in enumerate(Bseq):
        print(
            f"  B_{jj}: "
            f"degree_u={degree(expr,u)}"
        )

    print()
    print("=" * 78)
    print("3. FALLING-FACTORIAL STABILITY")
    print("=" * 78)

    A_C = diagonal_table(
        Aseq,
        6,
    )

    B_C = diagonal_table(
        Bseq,
        5,
    )

    A_falling_ok = True
    B_falling_ok = True

    for k_power in range(7):

        P = coefficient_index_polynomial(
            Aseq,
            k_power,
        )

        divides, Q = falling_test(
            P,
            k_power,
        )

        A_falling_ok &= divides

        print(
            f"  A k={k_power}: "
            f"falling_divides={divides} "
            f"quotient_degree={degree(Q,j)}"
        )

    for k_power in range(6):

        P = coefficient_index_polynomial(
            Bseq,
            k_power,
        )

        divides, Q = falling_test(
            P,
            k_power,
        )

        B_falling_ok &= divides

        print(
            f"  B k={k_power}: "
            f"falling_divides={divides} "
            f"quotient_degree={degree(Q,j)}"
        )

    print()
    print("=" * 78)
    print("4. ORIGINAL FALLING-BASIS TRIANGLE")
    print("=" * 78)

    print("  A channel")

    for k_power in range(7):

        row = A_C[k_power]

        print(
            f"    k={k_power}:",
            [
                (
                    r,
                    row.get(
                        r,
                        sp.Integer(0)
                    )
                )
                for r in range(
                    k_power,
                    9,
                )
            ]
        )

    print()
    print("  B channel")

    for k_power in range(6):

        row = B_C[k_power]

        print(
            f"    k={k_power}:",
            [
                (
                    r,
                    row.get(
                        r,
                        sp.Integer(0)
                    )
                )
                for r in range(
                    k_power,
                    8,
                )
            ]
        )

    print()
    print("=" * 78)
    print("5. DIAGONAL / SUPERDIAGONAL RATIO AUDIT")
    print("=" * 78)

    A_ratio_ok = True
    B_ratio_ok = True

    for shift in range(4):

        print()
        print(
            f"  A shift={shift}"
        )

        vals = []

        for k_power in range(
            7 - shift
        ):

            row = A_C[k_power]

            a = row.get(
                k_power + shift,
                sp.Integer(0)
            )

            vals.append(
                (
                    k_power,
                    a,
                )
            )

        print(
            "    values =",
            [v for _, v in vals]
        )

        transition_values = []

        for idx in range(
            len(vals) - 1
        ):

            a = vals[idx][1]
            b = vals[idx + 1][1]

            rr = exact_ratio(
                b,
                a,
            )

            if rr is None:
                continue

            transition_values.append(
                (
                    vals[idx][0],
                    rr,
                )
            )

        print(
            "    successive_ratios =",
            transition_values,
        )

        fit = None

        if len(transition_values) >= 3:
            fit = rational_ratio_fit(
                transition_values
            )

        print(
            "    exact_degree<=2_ratio =",
            fit
        )

        if fit is not None:
            A_ratio_ok = False

    for shift in range(4):

        print()
        print(
            f"  B shift={shift}"
        )

        vals = []

        for k_power in range(
            6 - shift
        ):

            row = B_C[k_power]

            a = row.get(
                k_power + shift,
                sp.Integer(0)
            )

            vals.append(
                (
                    k_power,
                    a,
                )
            )

        print(
            "    values =",
            [v for _, v in vals]
        )

        transition_values = []

        for idx in range(
            len(vals) - 1
        ):

            a = vals[idx][1]
            b = vals[idx + 1][1]

            rr = exact_ratio(
                b,
                a,
            )

            if rr is None:
                continue

            transition_values.append(
                (
                    vals[idx][0],
                    rr,
                )
            )

        print(
            "    successive_ratios =",
            transition_values,
        )

        fit = None

        if len(transition_values) >= 3:
            fit = rational_ratio_fit(
                transition_values
            )

        print(
            "    exact_degree<=2_ratio =",
            fit
        )

        if fit is not None:
            B_ratio_ok = False

    print()
    print("=" * 78)
    print("6. EXACT CROSS-CHANNEL DIAGONAL COMPARISON")
    print("=" * 78)

    for shift in range(4):

        print()
        print(
            f"  shift={shift}"
        )

        max_k = min(
            7 - shift,
            6 - shift,
        )

        for k_power in range(
            max_k
        ):

            a = A_C[k_power].get(
                k_power + shift,
                sp.Integer(0)
            )

            b = B_C[k_power].get(
                k_power + shift,
                sp.Integer(0)
            )

            if b != 0:

                print(
                    f"    k={k_power}: "
                    f"A/B={clean(a/b)}"
                )

    print()
    print("=" * 78)
    print("7. EXACT CONSISTENCY CHECK")
    print("=" * 78)

    consistency = True

    for k_power in range(7):

        P = coefficient_index_polynomial(
            Aseq,
            k_power,
        )

        row = A_C[k_power]

        reconstructed = sp.Integer(0)

        for r, c in row.items():
            reconstructed += (
                c * falling(j, r)
            )

        ok = (
            clean(
                reconstructed - P
            ) == 0
        )

        consistency &= ok

        print(
            f"  A k={k_power}: {ok}"
        )

    for k_power in range(6):

        P = coefficient_index_polynomial(
            Bseq,
            k_power,
        )

        row = B_C[k_power]

        reconstructed = sp.Integer(0)

        for r, c in row.items():
            reconstructed += (
                c * falling(j, r)
            )

        ok = (
            clean(
                reconstructed - P
            ) == 0
        )

        consistency &= ok

        print(
            f"  B k={k_power}: {ok}"
        )

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  This run deliberately does NOT attempt a new large recurrence search.

  The established exact structure is

      P_k(j) = sum_r C[k,r] j_(r),

  with

      j_(r) = j(j-1)...(j-r+1).

  The present audit asks whether the diagonals

      C[k,k+s]

  admit genuinely low-complexity rational transition laws.

  A rational law of degree <= 2 is only accepted if it satisfies every
  available transition exactly.

  This is intentionally conservative: a formula discovered from only
  a handful of points is treated as suspicious rather than promoted
  to a construction law.

  The important established facts remain:

      1. falling-factorial divisibility is exact;

      2. the original coefficient triangle is exactly triangular;

      3. the A/B involution decomposition is exact;

      4. the observed layer family reconstructs exactly.

  A negative result here is useful: it tells us that the diagonal
  coefficients are not obviously hypergeometric at this complexity.

  The next useful direction after this audit is to compare the
  coefficient triangle against the actual combinatorial quantities
  appearing in the original kernel definition, rather than continuing
  unrestricted interpolation searches.
        """
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  involution =",
        involution_ok,
    )

    print(
        "  A falling divisibility =",
        A_falling_ok,
    )

    print(
        "  B falling divisibility =",
        B_falling_ok,
    )

    print(
        "  falling-basis consistency =",
        consistency,
    )

    print(
        "  ALL BASIC CHECKS PASS =",
        (
            involution_ok
            and A_falling_ok
            and B_falling_ok
            and consistency
        ),
    )

    print()
    print(
        "EXPERIMENT 91 COMPLETE"
    )


if __name__ == "__main__":
    main()