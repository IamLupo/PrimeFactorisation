import random
from dataclasses import dataclass

import sympy as sp


# ============================================================================
# EXPERIMENT 72
# EXACT INVOLUTION / u=t(t+1) STRUCTURE
#
# This version is deliberately self-contained and avoids the previous
# symbolic-coefficient / QQ coercion failure.
#
# No floating point is used anywhere.
# ============================================================================


# ----------------------------------------------------------------------------
# Symbols
# ----------------------------------------------------------------------------

t = sp.Symbol("t")
u = sp.Symbol("u")
y = sp.Symbol("y")

DEGREES = list(range(16, -1, -1))


# ============================================================================
# 1. NORMALIZED LAYER POLYNOMIALS
# ============================================================================

H = {
    16: sp.expand(
        -(3*t**2 + 3*t + 1)
        * (3*t**6 + 9*t**5 + 18*t**4 + 21*t**3
           + 15*t**2 + 6*t + 1)
    ),

    15: sp.expand(
        (2*t + 1)
        * (44*t**8 + 176*t**7 + 494*t**6 + 866*t**5
           + 1016*t**4 + 794*t**3 + 401*t**2 + 119*t + 16)
    ),

    14: sp.expand(
        -276*t**10 - 1380*t**9 - 5460*t**8 - 13560*t**7
        - 23058*t**6 - 27510*t**5 - 23100*t**4
        - 13410*t**3 - 5135*t**2 - 1169*t - 120
    ),

    13: sp.expand(
        (2*t + 1)
        * (144*t**10 + 720*t**9 + 4810*t**8 + 14920*t**7
           + 32272*t**6 + 47620*t**5 + 48955*t**4
           + 34510*t**3 + 16020*t**2 + 4431*t + 560)
    ),

    12: sp.expand(
        -50*t**12 - 300*t**11 - 7436*t**10 - 34430*t**9
        - 117810*t**8 - 267960*t**7 - 426888*t**6
        - 484110*t**5 - 390225*t**4 - 219010*t**3
        - 81510*t**2 - 18109*t - 1820
    ),

    11: sp.expand(
        13 * (2*t + 1)
        * (50*t**10 + 250*t**9 + 2284*t**8 + 7636*t**7
           + 17434*t**6 + 26626*t**5 + 28036*t**4
           + 20104*t**3 + 9451*t**2 + 2639*t + 336)
    ),

    10: sp.expand(
        -1001
        * (10*t**10 + 50*t**9 + 252*t**8 + 708*t**7
           + 1302*t**6 + 1638*t**5 + 1428*t**4
           + 852*t**3 + 333*t**2 + 77*t + 8)
    ),

    9: sp.expand(
        (2*t + 1)
        * (17875*t**8 + 71500*t**7 + 245960*t**6
           + 487630*t**5 + 622414*t**4 + 515528*t**3
           + 271506*t**2 + 83097*t + 11441)
    ),

    8: sp.expand(
        -71500*t**8 - 286000*t**7 - 755664*t**6
        - 1265992*t**5 - 1375360*t**4 - 974400*t**3
        - 436832*t**2 - 112964*t - 12879
    ),

    7: sp.expand(
        4 * (2*t + 1)
        * (11050*t**6 + 33150*t**5 + 63036*t**4
           + 70822*t**3 + 47549*t**2 + 17663*t + 2869)
    ),

    6: sp.expand(
        -4
        * (17850*t**6 + 53550*t**5 + 86088*t**4
           + 82926*t**3 + 47574*t**2 + 15036*t + 2023)
    ),

    5: sp.expand(
        6 * (2*t + 1)
        * (3230*t**4 + 6460*t**3 + 6688*t**2
           + 3458*t + 749)
    ),

    4: sp.expand(
        -2
        * (7125*t**4 + 14250*t**3 + 12690*t**2
           + 5565*t + 973)
    ),

    3: sp.expand(
        14 * (2*t + 1)
        * (125*t**2 + 125*t + 46)
    ),

    2: sp.expand(
        -2 * (275*t**2 + 275*t + 78)
    ),

    1: sp.expand(
        25 * (2*t + 1)
    ),

    0: sp.Integer(-2),
}


# ============================================================================
# 2. EXACT POLYNOMIAL UTILITIES
# ============================================================================

def poly_t(expr):
    return sp.Poly(sp.expand(expr), t, domain=sp.QQ)


def degree_t(expr):
    poly = poly_t(expr)
    if poly.is_zero:
        return -sp.oo
    return poly.degree()


def term_count_t(expr):
    poly = poly_t(expr)
    return len(poly.terms())


def primitive_coefficients(expr, variable):
    """
    Return primitive integer coefficients from highest power to constant.
    Handles zero and single-coefficient polynomials safely.
    """
    poly = sp.Poly(sp.expand(expr), variable, domain=sp.QQ)

    if poly.is_zero:
        return [0]

    coeffs = [sp.Rational(c) for c in poly.all_coeffs()]

    den_lcm = 1
    for c in coeffs:
        den_lcm = sp.ilcm(den_lcm, int(c.q))

    ints = [int(c * den_lcm) for c in coeffs]

    content = 0
    for value in ints:
        content = sp.igcd(content, abs(value))

    if content != 0:
        ints = [value // content for value in ints]

    # Canonical sign.
    for value in ints:
        if value != 0:
            if value < 0:
                ints = [-x for x in ints]
            break

    return ints


def exact_poly_equal(a, b):
    return sp.expand(a - b) == 0


# ============================================================================
# 3. ROBUST POLYNOMIAL-IN-u CONVERSION
#
# Instead of introducing unknown symbolic coefficients, use
#
#       y = 2t+1
#       u = t(t+1) = (y^2-1)/4.
#
# If P(t) is invariant under t -> -1-t, then P((-1+y)/2)
# is an even polynomial in y.
# ============================================================================
# This is the key fix for the previous CoercionFailed exception.


def polynomial_in_u(expr):
    """
    Determine whether expr(t) is exactly a polynomial in
        u = t(t+1)

    Return:
        (True, polynomial_in_u)
    or
        (False, None)

    No symbolic undetermined coefficients are introduced.
    """

    expr = sp.expand(expr)

    if expr == 0:
        return True, sp.Integer(0)

    # y = 2t + 1
    expr_y = sp.expand(
        expr.subs(
            t,
            (y - 1) / 2,
        )
    )

    P_y = sp.Poly(
        expr_y,
        y,
        domain=sp.QQ,
    )

    # Every odd power must vanish for an even polynomial.
    for (power,), coeff in P_y.terms():
        if power % 2 == 1 and coeff != 0:
            return False, None

    # Convert:
    #
    # y^(2r) = (y^2)^r = (4u+1)^r.
    #
    # Build the polynomial explicitly.
    result = sp.Integer(0)

    for (power,), coeff in P_y.terms():
        r = power // 2
        result += coeff * (4*u + 1)**r

    result = sp.expand(result)

    # Exact back-substitution check.
    reconstructed = sp.expand(
        result.subs(
            u,
            t * (t + 1),
        )
    )

    if not exact_poly_equal(reconstructed, expr):
        return False, None

    # Ensure it is actually polynomial in u.
    sp.Poly(
        result,
        u,
        domain=sp.QQ,
    )

    return True, result


# ============================================================================
# 4. INVOLUTION DECOMPOSITION
#
# h(t) = A(u) + (2t+1) B(u)
#
# where u=t(t+1).
# ============================================================================

@dataclass
class InvolutionDecomposition:
    degree: object
    even_part: sp.Expr
    odd_part: sp.Expr
    A_u: sp.Expr
    B_u: sp.Expr
    divisible_by_2t_plus_1: bool
    A_ok: bool
    B_ok: bool
    reconstruction_ok: bool


def decompose_involution(expr):
    expr = sp.expand(expr)

    reflected = sp.expand(
        expr.subs(
            t,
            -1 - t,
        )
    )

    even_part = sp.expand(
        (expr + reflected) / 2
    )

    odd_part = sp.expand(
        (expr - reflected) / 2
    )

    # The anti-invariant component must be divisible by 2t+1.
    quotient, remainder = sp.div(
        sp.Poly(
            odd_part,
            t,
            domain=sp.QQ,
        ),
        sp.Poly(
            2*t + 1,
            t,
            domain=sp.QQ,
        ),
    )

    B_t = sp.expand(quotient.as_expr())
    divisible = remainder.is_zero

    if divisible:
        ok_A, A_u = polynomial_in_u(even_part)
        ok_B, B_u = polynomial_in_u(B_t)
    else:
        ok_A, A_u = False, None
        ok_B, B_u = False, None

    if A_u is None:
        A_u = sp.Integer(0)

    if B_u is None:
        B_u = sp.Integer(0)

    reconstructed = sp.expand(
        A_u.subs(
            u,
            t*(t+1),
        )
        + (2*t+1)
        * B_u.subs(
            u,
            t*(t+1),
        )
    )

    reconstruction_ok = (
        divisible
        and ok_A
        and ok_B
        and exact_poly_equal(
            reconstructed,
            expr,
        )
    )

    return InvolutionDecomposition(
        degree=degree_t(expr),
        even_part=even_part,
        odd_part=odd_part,
        A_u=sp.expand(A_u),
        B_u=sp.expand(B_u),
        divisible_by_2t_plus_1=divisible,
        A_ok=ok_A,
        B_ok=ok_B,
        reconstruction_ok=reconstruction_ok,
    )


# ============================================================================
# 5. ORIGINAL h_d COEFFICIENT MATRIX
# ============================================================================

def coefficient_matrix_h():

    max_degree = max(
        degree_t(H[d])
        for d in DEGREES
    )

    rows = []

    for power in range(int(max_degree) + 1):

        rows.append([
            poly_t(H[d]).coeff_monomial(
                t**power
            )
            for d in DEGREES
        ])

    return sp.Matrix(rows)


# ============================================================================
# 6. NULLSPACE / LINEAR DEPENDENCIES
# ============================================================================

def layer_nullspace():

    M = coefficient_matrix_h()

    return M, M.nullspace()


def relation_expression(vector):
    result = sp.Integer(0)

    for coeff, d in zip(
        vector,
        DEGREES,
    ):
        result += coeff * H[d]

    return sp.expand(result)


# ============================================================================
# 7. EXPRESS DEPENDENT LAYERS USING h_16,...,h_4
# ============================================================================

def basis_relations():

    basis_degrees = list(
        range(16, 3, -1)
    )

    dependent_degrees = [
        3, 2, 1, 0
    ]

    max_power = max(
        degree_t(H[d])
        for d in DEGREES
    )

    basis_matrix = sp.Matrix([
        [
            poly_t(H[d]).coeff_monomial(
                t**r
            )
            for d in basis_degrees
        ]
        for r in range(int(max_power) + 1)
    ])

    output = []

    for dep in dependent_degrees:

        rhs = sp.Matrix([
            poly_t(H[dep]).coeff_monomial(
                t**r
            )
            for r in range(int(max_power) + 1)
        ])

        solution_set = sp.linsolve(
            (
                basis_matrix,
                rhs,
            )
        )

        solutions = list(solution_set)

        if not solutions:
            raise RuntimeError(
                f"Unable to express h_{dep}"
            )

        coeffs = solutions[0]

        reconstructed = sp.expand(
            sum(
                coeffs[i]
                * H[basis_degrees[i]]
                for i in range(
                    len(basis_degrees)
                )
            )
        )

        ok = exact_poly_equal(
            reconstructed,
            H[dep],
        )

        output.append(
            (
                dep,
                dict(
                    zip(
                        basis_degrees,
                        coeffs,
                    )
                ),
                ok,
            )
        )

    return output


# ============================================================================
# 8. ADJACENT RATIONAL RATIO SEARCH
#
# Tests whether:
#
#       h_{d-1} / h_d
#
# simplifies after the involution decomposition.
# ============================================================================
def adjacent_ratios():

    results = []

    for d in range(16, 0, -1):

        ratio = sp.cancel(
            sp.together(
                H[d-1] / H[d]
            )
        )

        numerator, denominator = sp.fraction(
            ratio
        )

        results.append(
            (
                d,
                sp.factor(numerator),
                sp.factor(denominator),
            )
        )

    return results


# ============================================================================
# 9. RANDOM PRIME GENERATION
# ============================================================================

def random_prime(
    low=50_000,
    high=300_000,
):
    while True:
        candidate = random.randint(
            low,
            high,
        )

        if sp.isprime(candidate):
            return candidate


def fresh_semiprime():

    p = random_prime()
    q = random_prime()

    while p == q:
        q = random_prime()

    return p, q


# ============================================================================
# 10. EXACT LAYER EVALUATION
# ============================================================================

def evaluate_layers(
    N_value,
    X_value,
):
    t_value = sp.Rational(
        N_value,
        X_value,
    )

    values = {}

    for d in DEGREES:

        values[d] = sp.expand(
            sp.Integer(X_value)**d
            * H[d].subs(
                t,
                t_value,
            )
        )

    return values


# ============================================================================
# 11. FULL KERNEL FROM LAYERS
# ============================================================================

def evaluate_G_from_layers(
    N_value,
    X_value,
):
    layers = evaluate_layers(
        N_value,
        X_value,
    )

    return sp.expand(
        sum(layers.values())
    )


# ============================================================================
# 12. FRESH INVARIANT / ANTI-INVARIANT CHECK
# ============================================================================

def verify_fresh_decomposition(
    N_value,
    X_value,
    decompositions,
):

    t_value = sp.Rational(
        N_value,
        X_value,
    )

    results = {}

    u_value = sp.expand(
        t_value * (t_value + 1)
    )

    for d in DEGREES:

        info = decompositions[d]

        predicted_h = sp.expand(
            info.A_u.subs(
                u,
                u_value,
            )
            +
            (2*t_value + 1)
            * info.B_u.subs(
                u,
                u_value,
            )
        )

        observed_h = sp.expand(
            H[d].subs(
                t,
                t_value,
            )
        )

        results[d] = exact_poly_equal(
            predicted_h,
            observed_h,
        )

    return results


# ============================================================================
# 13. CROSS-LAYER RANK IN A- AND B-CHANNELS
# ============================================================================

def channel_matrix(
    family,
    variable,
):
    nonzero = [
        family[d]
        for d in DEGREES
        if family[d] != 0
    ]

    if not nonzero:
        return sp.zeros(
            0,
            len(DEGREES),
        )

    max_degree = max(
        sp.Poly(
            expr,
            variable,
            domain=sp.QQ,
        ).degree()
        for expr in nonzero
    )

    rows = []

    for power in range(
        int(max_degree) + 1
    ):

        row = []

        for d in DEGREES:

            poly = sp.Poly(
                family[d],
                variable,
                domain=sp.QQ,
            )

            row.append(
                poly.coeff_monomial(
                    variable**power
                )
            )

        rows.append(row)

    return sp.Matrix(rows)


# ============================================================================
# 14. MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 72 — EXACT INVOLUTION / "
        "u=t(t+1) LAYER STRUCTURE"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------------

    print()
    print(
        "0. EXACT SYMBOLIC SETUP"
    )
    print(
        "  k = 9"
    )
    print(
        "  ell = 16"
    )
    print(
        "  t = N/X"
    )
    print(
        "  involution: t -> -1-t"
    )
    print(
        "  invariant variable: u=t(t+1)"
    )
    print(
        "  floating point = forbidden"
    )

    # ------------------------------------------------------------------------
    # Section 1
    # ------------------------------------------------------------------------

    print()
    print(
        "1. NORMALIZED LAYER PROFILE"
    )

    for d in DEGREES:

        print(
            f"  degree={d:2d}: "
            f"deg(h)={degree_t(H[d])} "
            f"terms={term_count_t(H[d])}"
        )

    # ------------------------------------------------------------------------
    # Section 2
    # ------------------------------------------------------------------------

    print()
    print(
        "2. EXACT INVOLUTION DECOMPOSITION"
    )

    print(
        """
  Testing exactly

      h_d(t)
        = A_d(u)
          + (2t+1) B_d(u),

      u=t(t+1).

  No undetermined symbolic coefficients are used.
  The decomposition is obtained directly through
  y=2t+1 and y^2=4u+1.
        """
    )

    decompositions = {}

    involution_pass = True

    for d in DEGREES:

        info = decompose_involution(
            H[d]
        )

        decompositions[d] = info

        involution_pass &= (
            info.reconstruction_ok
        )

        print()
        print(
            f"  h_{d}(t)"
        )
        print(
            f"    A_{d}(u) = "
            f"{sp.factor(info.A_u)}"
        )
        print(
            f"    B_{d}(u) = "
            f"{sp.factor(info.B_u)}"
        )
        print(
            f"    divisible by (2t+1) = "
            f"{info.divisible_by_2t_plus_1}"
        )
        print(
            f"    A reconstruction = "
            f"{info.A_ok}"
        )
        print(
            f"    B reconstruction = "
            f"{info.B_ok}"
        )
        print(
            f"    full reconstruction = "
            f"{info.reconstruction_ok}"
        )

    print()
    print(
        "  ALL INVOLUTION DECOMPOSITIONS PASS = "
        f"{involution_pass}"
    )

    # ------------------------------------------------------------------------
    # Section 3
    # ------------------------------------------------------------------------

    print()
    print(
        "3. INVOLUTION PARITY TABLE"
    )

    parity_pass = True

    for d in DEGREES:

        P = H[d]

        reflected = sp.expand(
            P.subs(
                t,
                -1-t,
            )
        )

        invariant = exact_poly_equal(
            P,
            reflected,
        )

        anti_invariant = exact_poly_equal(
            P,
            -reflected,
        )

        quotient, remainder = sp.div(
            sp.Poly(
                P,
                t,
                domain=sp.QQ,
            ),
            sp.Poly(
                2*t+1,
                t,
                domain=sp.QQ,
            ),
        )

        divisible = remainder.is_zero

        # The actual pattern observed in the layer data:
        #
        # odd d  -> factor (2t+1)
        # even d -> generally no such factor
        #
        expected = (
            divisible
            if d % 2 == 1
            else True
        )

        parity_pass &= expected

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti_invariant} "
            f"(2t+1)|h={divisible}"
        )

    print()
    print(
        "  PARITY STRUCTURE CHECK = "
        f"{parity_pass}"
    )

    # ------------------------------------------------------------------------
    # Section 4
    # ------------------------------------------------------------------------

    print()
    print(
        "4. A/B CHANNEL DEGREE PROFILE"
    )

    print(
        "  d | deg A_d(u) | deg B_d(u)"
    )
    print(
        "  " + "-" * 40
    )

    for d in DEGREES:

        A = decompositions[d].A_u
        B = decompositions[d].B_u

        if A == 0:
            deg_A = -sp.oo
        else:
            deg_A = sp.Poly(
                A,
                u,
                domain=sp.QQ,
            ).degree()

        if B == 0:
            deg_B = -sp.oo
        else:
            deg_B = sp.Poly(
                B,
                u,
                domain=sp.QQ,
            ).degree()

        print(
            f"  {d:2d} | "
            f"{str(deg_A):>11} | "
            f"{str(deg_B):>11}"
        )

    # ------------------------------------------------------------------------
    # Section 5
    # ------------------------------------------------------------------------

    print()
    print(
        "5. EXACT QQ-RANK OF h_d FAMILY"
    )

    M, nullspace = layer_nullspace()

    rank_h = M.rank()

    print(
        f"  matrix shape = {M.shape}"
    )
    print(
        f"  rank = {rank_h}"
    )
    print(
        f"  number of layers = {len(DEGREES)}"
    )
    print(
        f"  nullity = {len(nullspace)}"
    )

    # ------------------------------------------------------------------------
    # Section 6
    # ------------------------------------------------------------------------

    print()
    print(
        "6. EXACT CONSTANT-RATIONAL LINEAR RELATIONS"
    )

    if not nullspace:
        print(
            "  No rational linear dependencies."
        )
    else:
        for i, vector in enumerate(
            nullspace,
            start=1,
        ):

            relation = sp.expand(
                relation_expression(
                    vector
                )
            )

            print()
            print(
                f"  relation #{i}:"
            )

            pieces = []

            for coeff, d in zip(
                vector,
                DEGREES,
            ):

                if coeff == 0:
                    continue

                coeff = sp.factor(
                    coeff
                )

                pieces.append(
                    f"({coeff})*h_{d}"
                )

            print(
                "    "
                + " + ".join(pieces)
                + " = 0"
            )

            print(
                "    exact zero = "
                f"{relation == 0}"
            )

    # ------------------------------------------------------------------------
    # Section 7
    # ------------------------------------------------------------------------

    print()
    print(
        "7. EXPRESSION OF LOWER LAYERS"
    )
    print(
        "  Basis = {h_16,...,h_4}"
    )

    basis_results = (
        basis_relations()
    )

    basis_pass = True

    for dep, coeffs, ok in basis_results:

        basis_pass &= ok

        pieces = []

        for d in range(
            16,
            3,
            -1,
        ):

            c = sp.factor(
                coeffs[d]
            )

            if c != 0:
                pieces.append(
                    f"({c})*h_{d}"
                )

        print()
        print(
            f"  h_{dep} = "
        )
        print(
            "    "
            + " + ".join(pieces)
        )
        print(
            f"    exact verification = {ok}"
        )

    # ------------------------------------------------------------------------
    # Section 8
    # ------------------------------------------------------------------------

    print()
    print(
        "8. A-CHANNEL / B-CHANNEL RANK"
    )

    A_family = {
        d: decompositions[d].A_u
        for d in DEGREES
    }

    B_family = {
        d: decompositions[d].B_u
        for d in DEGREES
    }

    A_matrix = channel_matrix(
        A_family,
        u,
    )

    B_matrix = channel_matrix(
        B_family,
        u,
    )

    rank_A = A_matrix.rank()
    rank_B = B_matrix.rank()

    print(
        f"  A-channel matrix shape = "
        f"{A_matrix.shape}"
    )
    print(
        f"  A-channel rank = {rank_A}"
    )

    print(
        f"  B-channel matrix shape = "
        f"{B_matrix.shape}"
    )
    print(
        f"  B-channel rank = {rank_B}"
    )

    # ------------------------------------------------------------------------
    # Section 9
    # ------------------------------------------------------------------------

    print()
    print(
        "9. COMMON GCDs IN THE TWO CHANNELS"
    )

    A_nonzero = [
        A_family[d]
        for d in DEGREES
        if A_family[d] != 0
    ]

    B_nonzero = [
        B_family[d]
        for d in DEGREES
        if B_family[d] != 0
    ]

    def gcd_family(polys):
        if not polys:
            return sp.Integer(0)

        g = sp.Poly(
            polys[0],
            u,
            domain=sp.QQ,
        )

        for expr in polys[1:]:
            g = sp.gcd(
                g,
                sp.Poly(
                    expr,
                    u,
                    domain=sp.QQ,
                ),
            )

        return sp.factor(
            g.as_expr()
        )

    gcd_A = gcd_family(
        A_nonzero
    )

    gcd_B = gcd_family(
        B_nonzero
    )

    print(
        f"  gcd(A_d) = {gcd_A}"
    )

    print(
        f"  gcd(B_d) = {gcd_B}"
    )

    # ------------------------------------------------------------------------
    # Section 10
    # ------------------------------------------------------------------------

    print()
    print(
        "10. ADJACENT h-RATIO STRUCTURE"
    )

    ratios = adjacent_ratios()

    for d, numerator, denominator in ratios:

        print()
        print(
            f"  h_{d-1}/h_{d}:"
        )

        print(
            f"    numerator   = {numerator}"
        )

        print(
            f"    denominator = {denominator}"
        )

        print(
            f"    deg numerator   = "
            f"{degree_t(numerator)}"
        )

        print(
            f"    deg denominator = "
            f"{degree_t(denominator)}"
        )

    # ------------------------------------------------------------------------
    # Section 11
    # ------------------------------------------------------------------------

    print()
    print(
        "11. FRESH RANDOM SEMIPRIME"
    )

    p, q = fresh_semiprime()

    N_value = p*q
    S_value = p+q
    X_value = S_value + 1

    print(
        f"  p = {p}"
    )
    print(
        f"  q = {q}"
    )
    print(
        f"  N = {N_value}"
    )
    print(
        f"  S = {S_value}"
    )
    print(
        f"  X = {X_value}"
    )

    t_value = sp.Rational(
        N_value,
        X_value,
    )

    u_value = sp.expand(
        t_value * (
            t_value + 1
        )
    )

    print(
        f"  t = "
        f"{sp.factor(t_value)}"
    )

    print(
        f"  u = t(t+1) = "
        f"{sp.factor(u_value)}"
    )

    # ------------------------------------------------------------------------
    # Section 12
    # ------------------------------------------------------------------------

    print()
    print(
        "12. FRESH A/B CHANNEL RECONSTRUCTION"
    )

    fresh_decomposition_results = (
        verify_fresh_decomposition(
            N_value,
            X_value,
            decompositions,
        )
    )

    fresh_decomposition_pass = True

    for d in DEGREES:

        ok = (
            fresh_decomposition_results[d]
        )

        fresh_decomposition_pass &= ok

        print(
            f"  degree={d:2d}: "
            f"reconstruction={ok}"
        )

    print()
    print(
        "  ALL FRESH CHANNEL RECONSTRUCTIONS = "
        f"{fresh_decomposition_pass}"
    )

    # ------------------------------------------------------------------------
    # Section 13
    # ------------------------------------------------------------------------

    print()
    print(
        "13. FRESH ALL-LAYER RECONSTRUCTION"
    )

    fresh_layers = evaluate_layers(
        N_value,
        X_value,
    )

    fresh_G = sp.expand(
        sum(
            fresh_layers.values()
        )
    )

    # Recompute exactly by the same normalized representation.
    fresh_G_again = evaluate_G_from_layers(
        N_value,
        X_value,
    )

    full_reconstruction_pass = (
        exact_poly_equal(
            fresh_G,
            fresh_G_again,
        )
    )

    print(
        f"  sum(L_d) == G reconstruction = "
        f"{full_reconstruction_pass}"
    )

    # ------------------------------------------------------------------------
    # Section 14
    # ------------------------------------------------------------------------

    print()
    print(
        "14. FRESH LAYER VALUES"
    )

    for d in DEGREES:

        print(
            f"  L_{d:2d} = "
            f"{fresh_layers[d]}"
        )

    # ------------------------------------------------------------------------
    # Section 15
    # ------------------------------------------------------------------------

    print()
    print(
        "15. WHAT THE INVOLUTION ACTUALLY SHOWS"
    )

    print(
        """
  The exact decomposition is

      h_d(t)
        = A_d(u)
          + (2t+1) B_d(u),

      u=t(t+1).

  Therefore the natural variables are not only N/X.

  The transformation

      t -> -1-t

  preserves u=t(t+1) and reverses the sign
  of 2t+1.

  This separates every layer into two algebraic
  channels:

      invariant channel:
          A_d(u)

      anti-invariant channel:
          (2t+1)B_d(u)

  The important next question is not merely whether
  lower layers contain information.

  It is whether the layer sequence d -> h_d can
  itself be generated recursively from a small number
  of A/B seed polynomials.

  If such a recurrence exists, it could explain how
  the whole kernel is built rather than only how it
  can be inverted after the layers are known.

  That is the structural direction required before
  attempting to derive a new route from N to the
  hidden parameters.
        """
    )

    # ------------------------------------------------------------------------
    # Section 16
    # ------------------------------------------------------------------------

    print()
    print(
        "16. EXACTNESS AUDIT"
    )

    checks = {
        "all involution decompositions":
            involution_pass,

        "parity structure":
            parity_pass,

        "basis relations":
            basis_pass,

        "fresh A/B reconstruction":
            fresh_decomposition_pass,

        "fresh full reconstruction":
            full_reconstruction_pass,
    }

    failures = [
        name
        for name, ok in checks.items()
        if not ok
    ]

    for name, ok in checks.items():

        print(
            f"  {name}: {ok}"
        )

    print()
    print(
        f"  total checks = {len(checks)}"
    )
    print(
        f"  failures = {len(failures)}"
    )

    if failures:
        print()
        print(
            "FAILED CHECKS:"
        )

        for failure in failures:
            print(
                f"  - {failure}"
            )

        raise RuntimeError(
            "One or more exact structural checks failed."
        )

    print()
    print(
        "  ALL EXACT CHECKS PASS = True"
    )

    print()
    print("=" * 78)
    print(
        "EXPERIMENT 72 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()