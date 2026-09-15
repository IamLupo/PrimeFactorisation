#!/usr/bin/env python3

"""
EXPERIMENT 79 — EXACT LOCAL-OPERATOR INDEX-LAW AUDIT
======================================================================

Purpose
-------
Experiment 78 found exact local coupled operators

    V_j(u) = M_j(u) V_{j-1}(u)

with

    V_j(u) = [ A_{2j}(u)   ]
             [ B_{2j+1}(u) ].

The observed minimal operator degrees were

    1, 2, 3, 4, 5, ...

This experiment investigates whether the local operators themselves
follow a compact law in j.

It performs:

  1. exact involution decomposition;
  2. exact reconstruction of A_d(u), B_d(u);
  3. exact recovery of local 2x2 polynomial operators;
  4. operator degree profile;
  5. determinant and trace profiles;
  6. coefficient sequences in j;
  7. exact finite-difference degrees in j;
  8. exact polynomial-in-j descriptions on the available finite range;
  9. exact validation on an independent semiprime substitution.

No floating point.
No numerical fitting.
No regression.
No approximate linear algebra.

Important:
A polynomial in j reconstructed from the finite available transitions
is reported only as an exact finite-range identity.
It is NOT promoted automatically to a universal formula.
"""


import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")
j = sp.Symbol("j")


# ============================================================================
# NORMALIZED LAYER POLYNOMIALS
# ============================================================================

H = {
    16: (
        -9*t**8
        -36*t**7
        -84*t**6
        -126*t**5
        -126*t**4
        -84*t**3
        -36*t**2
        -9*t
        -1
    ),

    15: (
        88*t**9
        +396*t**8
        +1164*t**7
        +2226*t**6
        +2898*t**5
        +2604*t**4
        +1596*t**3
        +639*t**2
        +151*t
        +16
    ),

    14: (
        -276*t**10
        -1380*t**9
        -5460*t**8
        -13560*t**7
        -23058*t**6
        -27510*t**5
        -23100*t**4
        -13410*t**3
        -5135*t**2
        -1169*t
        -120
    ),

    13: (
        288*t**11
        +1584*t**10
        +10340*t**9
        +34650*t**8
        +79464*t**7
        +127512*t**6
        +145530*t**5
        +117975*t**4
        +66550*t**3
        +24882*t**2
        +5551*t
        +560
    ),

    12: (
        -50*t**12
        -300*t**11
        -7436*t**10
        -34430*t**9
        -117810*t**8
        -267960*t**7
        -426888*t**6
        -484110*t**5
        -390225*t**4
        -219010*t**3
        -81510*t**2
        -18109*t
        -1820
    ),

    11: (
        1300*t**11
        +7150*t**10
        +62634*t**9
        +228228*t**8
        +552552*t**7
        +918918*t**6
        +1075074*t**5
        +887172*t**4
        +507078*t**3
        +191477*t**2
        +43043*t
        +4368
    ),

    10: (
        -10010*t**10
        -50050*t**9
        -252252*t**8
        -708708*t**7
        -1303302*t**6
        -1639638*t**5
        -1429428*t**4
        -852852*t**3
        -333333*t**2
        -77077*t
        -8008
    ),

    9: (
        35750*t**9
        +160875*t**8
        +563420*t**7
        +1221220*t**6
        +1732458*t**5
        +1653470*t**4
        +1058540*t**3
        +437700*t**2
        +105979*t
        +11441
    ),

    8: (
        -71500*t**8
        -286000*t**7
        -755664*t**6
        -1265992*t**5
        -1375360*t**4
        -974400*t**3
        -436832*t**2
        -112964*t
        -12879
    ),

    7: (
        88400*t**7
        +309400*t**6
        +636888*t**5
        +818720*t**4
        +663680*t**3
        +331500*t**2
        +93604*t
        +11476
    ),

    6: (
        -71400*t**6
        -214200*t**5
        -344352*t**4
        -331704*t**3
        -190296*t**2
        -60044*t
        -8092
    ),

    5: (
        38760*t**5
        +96900*t**4
        +119016*t**3
        +81624*t**2
        +29736*t
        +4494
    ),

    4: (
        -14250*t**4
        -28500*t**3
        -25380*t**2
        -11130*t
        -1946
    ),

    3: (
        3500*t**3
        +5250*t**2
        +3038*t
        +644
    ),

    2: (
        -550*t**2
        -550*t
        -156
    ),

    1: (
        50*t
        +25
    ),

    0: sp.Integer(-2),
}


# ============================================================================
# SAFE SYMBOLIC HELPERS
# ============================================================================

def E(expr):
    return sp.expand(sp.sympify(expr))


def is_zero(expr):
    return sp.simplify(E(expr)) == 0


def poly_over_t(expr):
    """
    Polynomial in t over an expression domain.

    This is deliberately NOT QQ because temporary symbolic coefficients
    appear during exact reconstruction.
    """
    return sp.Poly(
        E(expr),
        t,
        domain="EX",
    )


def poly_over_u(expr):
    """
    Polynomial in u with coefficients in QQ after all temporary symbols
    have been eliminated.
    """
    return sp.Poly(
        E(expr),
        u,
        domain=sp.QQ,
    )


def coeff_t(expr, power):
    return sp.expand(
        poly_over_t(expr).coeff_monomial(
            t**power
        )
    )


def coeff_u(expr, power):
    return sp.expand(
        poly_over_u(expr).coeff_monomial(
            u**power
        )
    )


def degree_t(expr):
    p = poly_over_t(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def degree_u(expr):
    if E(expr) == 0:
        return -sp.oo
    return poly_over_u(expr).degree()


def coeff_list_u(expr):
    p = poly_over_u(expr)

    if p.is_zero:
        return [sp.Integer(0)]

    return [
        sp.cancel(
            p.coeff_monomial(
                u**k
            )
        )
        for k in range(
            p.degree(),
            -1,
            -1,
        )
    ]


# ============================================================================
# EXACT RECOVERY OF A(u) FROM P(t(t+1))
# ============================================================================

def recover_u_polynomial(target_t):
    """
    Given an exact polynomial target_t satisfying

        target_t = A(t(t+1))

    recover A(u) exactly.

    Temporary coefficients are handled in EX.
    Only the final solved expression is converted to QQ.
    """

    target_t = E(target_t)

    if target_t == 0:
        return sp.Integer(0)

    deg_t = degree_t(target_t)

    if deg_t % 2 != 0:
        raise RuntimeError(
            "Expected an even t-degree for an invariant polynomial."
        )

    deg_u = deg_t // 2

    unknowns = sp.symbols(
        "a0:" + str(deg_u + 1)
    )

    candidate = E(
        sum(
            unknowns[k] * u**k
            for k in range(
                deg_u + 1
            )
        )
    )

    difference = E(
        candidate.subs(
            u,
            t*(t+1),
        )
        -
        target_t
    )

    # IMPORTANT:
    # use EX here because difference contains symbolic unknowns.
    p_diff = sp.Poly(
        difference,
        t,
        domain="EX",
    )

    equations = []

    for power in range(
        p_diff.degree(),
        -1,
        -1,
    ):
        equations.append(
            sp.Eq(
                p_diff.coeff_monomial(
                    t**power
                ),
                0,
            )
        )

    if not equations:
        return sp.Integer(0)

    solution_set = sp.linsolve(
        equations,
        unknowns,
    )

    if solution_set == sp.EmptySet:
        raise RuntimeError(
            "No exact u-polynomial reconstruction exists."
        )

    tuples = list(solution_set)

    if len(tuples) != 1:
        raise RuntimeError(
            "Unexpectedly non-unique u-polynomial reconstruction."
        )

    solution = tuples[0]

    # Reject unresolved symbolic free parameters.
    free = set()

    for value in solution:
        free |= E(value).free_symbols

    free -= {t, u}

    if free:
        raise RuntimeError(
            "u-polynomial reconstruction left free symbols: "
            f"{sorted(map(str, free))}"
        )

    substitution = dict(
        zip(
            unknowns,
            solution,
        )
    )

    result = sp.cancel(
        E(
            candidate.subs(
                substitution
            )
        )
    )

    # Verify exactly in t.
    reconstructed = E(
        result.subs(
            u,
            t*(t+1),
        )
    )

    if not is_zero(
        reconstructed - target_t
    ):
        raise RuntimeError(
            "Exact u-polynomial reconstruction failed."
        )

    # Convert through QQ only after all temporary
    # coefficient symbols are gone.
    sp.Poly(
        result,
        u,
        domain=sp.QQ,
    )

    return result


# ============================================================================
# INVOLUTION DECOMPOSITION
# ============================================================================

def decompose_layer(d):

    h = E(H[d])

    partner = E(
        h.subs(
            t,
            -1-t,
        )
    )

    even_part = E(
        (h + partner) / 2
    )

    odd_part = E(
        (h - partner) / 2
    )

    A_d = recover_u_polynomial(
        even_part
    )

    if odd_part == 0:
        B_d = sp.Integer(0)
    else:

        odd_poly = sp.Poly(
            odd_part,
            t,
            domain=sp.QQ,
        )

        divisor = sp.Poly(
            2*t+1,
            t,
            domain=sp.QQ,
        )

        quotient, remainder = sp.div(
            odd_poly,
            divisor,
        )

        if not remainder.is_zero:
            raise RuntimeError(
                f"d={d}: odd channel is not divisible by 2t+1."
            )

        B_d = recover_u_polynomial(
            quotient.as_expr()
        )

    reconstructed = E(
        A_d.subs(
            u,
            t*(t+1),
        )
        +
        (2*t+1)
        *
        B_d.subs(
            u,
            t*(t+1),
        )
    )

    if not is_zero(
        reconstructed - h
    ):
        raise RuntimeError(
            f"d={d}: channel reconstruction failed."
        )

    return (
        sp.cancel(A_d),
        sp.cancel(B_d),
    )


# ============================================================================
# DECOMPOSE ALL LAYERS
# ============================================================================

CHANNEL = {}

for d in range(
    16,
    -1,
    -1,
):
    CHANNEL[d] = decompose_layer(d)


A = {
    d: CHANNEL[d][0]
    for d in range(
        0,
        17,
        2,
    )
}

B = {
    d: CHANNEL[d][1]
    for d in range(
        1,
        16,
        2,
    )
}


# ============================================================================
# COUPLED STATES
# ============================================================================

V = []

for j_idx in range(8):

    V.append(
        sp.Matrix([
            E(A[2*j_idx]),
            E(B[2*j_idx+1]),
        ])
    )


# ============================================================================
# OPERATOR CREATION
# ============================================================================

def create_operator(
    degree_bound,
    prefix,
):

    unknowns = []

    entries = []

    for r in range(2):

        row = []

        for c in range(2):

            coeffs = sp.symbols(
                f"{prefix}_{r}_{c}_0:{degree_bound+1}"
            )

            unknowns.extend(
                coeffs
            )

            row.append(
                E(
                    sum(
                        coeffs[k] * u**k
                        for k in range(
                            degree_bound + 1
                        )
                    )
                )
            )

        entries.append(row)

    return (
        sp.Matrix(entries),
        unknowns,
    )


# ============================================================================
# OPERATOR EQUATION GENERATION
# ============================================================================

def operator_equations(
    left,
    right,
    degree_bound,
    prefix,
):

    M, unknowns = create_operator(
        degree_bound,
        prefix,
    )

    predicted = E(
        M[0, 0] * right[0]
        +
        M[0, 1] * right[1]
    )

    predicted_2 = E(
        M[1, 0] * right[0]
        +
        M[1, 1] * right[1]
    )

    differences = [
        E(predicted - left[0]),
        E(predicted_2 - left[1]),
    ]

    equations = []

    for diff in differences:

        if E(diff) == 0:
            continue

        # IMPORTANT:
        # coefficients can contain unknown symbols, therefore EX.
        p = sp.Poly(
            E(diff),
            u,
            domain="EX",
        )

        for monomial_power in range(
            p.degree(),
            -1,
            -1,
        ):

            coefficient = sp.expand(
                p.coeff_monomial(
                    u**monomial_power
                )
            )

            if coefficient != 0:
                equations.append(
                    coefficient
                )

    return (
        M,
        unknowns,
        equations,
    )


# ============================================================================
# EXACT LOCAL OPERATOR SOLVER
# ============================================================================

def solve_operator(
    left,
    right,
    degree_bound,
    prefix,
):

    (
        M,
        unknowns,
        equations,
    ) = operator_equations(
        left,
        right,
        degree_bound,
        prefix,
    )

    if not equations:
        return None

    # The equations are linear in unknown coefficients.
    Aeq, beq = sp.linear_eq_to_matrix(
        equations,
        unknowns,
    )

    rank = Aeq.rank()

    augmented_rank = Aeq.row_join(
        beq
    ).rank()

    if augmented_rank != rank:

        return {
            "exists": False,
            "matrix": None,
            "rank": rank,
            "unknowns": len(unknowns),
            "nullity": len(unknowns)-rank,
            "unique": False,
        }

    solution_set = sp.linsolve(
        (
            Aeq,
            beq,
        ),
        unknowns,
    )

    if solution_set == sp.EmptySet:
        return {
            "exists": False,
            "matrix": None,
            "rank": rank,
            "unknowns": len(unknowns),
            "nullity": len(unknowns)-rank,
            "unique": False,
        }

    tuples = list(
        solution_set
    )

    if len(tuples) != 1:
        raise RuntimeError(
            "Unexpected linsolve result."
        )

    solution = tuples[0]

    free_parameters = set()

    for value in solution:
        free_parameters |= (
            E(value).free_symbols
        )

    free_parameters -= set(
        unknowns
    )

    # Every free symbolic quantity introduced by linsolve
    # represents a genuine non-unique family.
    unique = (
        len(free_parameters) == 0
        and all(
            E(value).free_symbols
            <= set()
            for value in solution
        )
    )

    # For a parameterized system, build one canonical representative
    # by assigning all free symbols to zero.
    substitution = {}

    all_solution_symbols = set()

    for value in solution:
        all_solution_symbols |= (
            E(value).free_symbols
        )

    # Remove genuine mathematical symbols; what remains are
    # linsolve parameters.
    free_linsolve = (
        all_solution_symbols
        &
        set().union(
            *[
                E(value).free_symbols
                for value in solution
            ]
        )
    )

    # Build zero assignment for parameters not among u.
    for symbol in free_linsolve:

        if symbol != u:
            substitution[symbol] = 0

    solved_entries = []

    for r in range(2):

        row = []

        for c in range(2):

            value = E(
                M[r, c].subs(
                    dict(
                        zip(
                            unknowns,
                            solution,
                        )
                    )
                )
            )

            if substitution:
                value = E(
                    value.subs(
                        substitution
                    )
                )

            # The operator may contain u but no other symbols.
            unexpected = (
                value.free_symbols
                -
                {u}
            )

            if unexpected:
                raise RuntimeError(
                    "Operator solution contains unresolved symbols: "
                    + ", ".join(
                        sorted(
                            map(
                                str,
                                unexpected,
                            )
                        )
                    )
                )

            # Final exact QQ coercion.
            value = sp.cancel(value)

            sp.Poly(
                value,
                u,
                domain=sp.QQ,
            )

            row.append(value)

        solved_entries.append(row)

    solved_matrix = sp.Matrix(
        solved_entries
    )

    return {
        "exists": True,
        "matrix": solved_matrix,
        "rank": rank,
        "unknowns": len(unknowns),
        "nullity": len(unknowns)-rank,
        "unique": unique,
        "parameterized": not unique,
    }


# ============================================================================
# RECOVER LOCAL OPERATORS
# ============================================================================

def recover_local_operator(
    transition,
):

    for degree_bound in range(
        0,
        8,
    ):

        result = solve_operator(
            V[transition],
            V[transition-1],
            degree_bound,
            f"M{transition}_{degree_bound}",
        )

        if result is not None and result["exists"]:

            return (
                degree_bound,
                result,
            )

    return None


# ============================================================================
# MATRIX HELPERS
# ============================================================================

def matrix_degree(M):

    finite = []

    for r in range(2):

        for c in range(2):

            d = degree_u(
                M[r, c]
            )

            if d != -sp.oo:
                finite.append(d)

    if not finite:
        return -sp.oo

    return max(finite)


def matrix_coeff(
    M,
    r,
    c,
    power,
):

    return sp.cancel(
        coeff_u(
            M[r, c],
            power,
        )
    )


def determinant(M):

    return sp.factor(
        E(
            M[0, 0]*M[1, 1]
            -
            M[0, 1]*M[1, 0]
        )
    )


def trace(M):

    return sp.factor(
        E(
            M[0, 0]
            +
            M[1, 1]
        )
    )


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def finite_difference(values):

    return [
        sp.cancel(
            values[k+1] - values[k]
        )
        for k in range(
            len(values)-1
        )
    ]


def difference_degree(values):

    current = list(values)
    degree = 0

    while len(current) > 1:

        if all(
            sp.simplify(
                x-current[0]
            ) == 0
            for x in current
        ):
            return degree

        current = finite_difference(
            current
        )

        degree += 1

    return degree


# ============================================================================
# EXACT LAGRANGE INTERPOLATION IN j
# ============================================================================

def interpolate_j(
    values,
    indices,
):

    result = sp.Integer(0)

    for value, ii in zip(
        values,
        indices,
    ):

        basis = sp.Integer(1)

        for jj in indices:

            if jj == ii:
                continue

            basis *= (
                j-jj
            ) / sp.Integer(
                ii-jj
            )

        result += (
            value*basis
        )

    result = sp.cancel(
        result
    )

    sp.Poly(
        result,
        j,
        domain=sp.QQ,
    )

    return result


def verify_j_identity(
    candidate,
    values,
    indices,
):

    for value, ii in zip(
        values,
        indices,
    ):

        check = sp.simplify(
            candidate.subs(
                j,
                ii,
            )
            -
            value
        )

        if check != 0:
            return False

    return True


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print(
        "EXPERIMENT 79 — EXACT LOCAL-OPERATOR INDEX-LAW AUDIT"
    )
    print("="*78)

    # ------------------------------------------------------------------------
    # 0
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("0. EXACT SYMBOLIC SETUP")
    print("="*78)

    print("  t = N/X")
    print("  u = t(t+1)")
    print("  V_j = [A_{2j}, B_{2j+1}]^T")
    print("  M_j: V_j = M_j(u) V_{j-1}")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # ------------------------------------------------------------------------
    # 1
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("="*78)

    involution_ok = True

    for d in range(
        16,
        -1,
        -1,
    ):

        h = E(H[d])

        partner = E(
            h.subs(
                t,
                -1-t,
            )
        )

        A_d, B_d = CHANNEL[d]

        reconstructed = E(
            A_d.subs(
                u,
                t*(t+1),
            )
            +
            (2*t+1)
            *
            B_d.subs(
                u,
                t*(t+1),
            )
        )

        ok = is_zero(
            reconstructed-h
        )

        invariant = is_zero(
            partner-h
        )

        anti = is_zero(
            partner+h
        )

        involution_ok &= ok

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"reconstruction={ok}"
        )

    print()
    print(
        "  ALL INVOLUTION CHECKS =",
        involution_ok,
    )

    # ------------------------------------------------------------------------
    # 2. LOCAL OPERATORS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("2. EXACT LOCAL OPERATORS")
    print("="*78)

    local_ops = {}

    for transition in range(
        1,
        len(V),
    ):

        result = recover_local_operator(
            transition
        )

        if result is None:

            print(
                f"  M_{transition}: NONE"
            )

            continue

        degree_bound, info = result

        local_ops[transition] = (
            degree_bound,
            info["matrix"],
            info,
        )

        M = info["matrix"]

        print()
        print(
            f"  M_{transition}:"
        )

        print(
            f"    V_{transition-1} -> V_{transition}"
        )

        print(
            "    degree =",
            degree_bound,
        )

        print(
            "    rank =",
            info["rank"],
        )

        print(
            "    unknowns =",
            info["unknowns"],
        )

        print(
            "    nullity =",
            info["nullity"],
        )

        print(
            "    unique =",
            info["unique"],
        )

        print(
            "    matrix ="
        )

        print(M)

    # ------------------------------------------------------------------------
    # 3. DEGREE LAW
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("3. LOCAL OPERATOR DEGREE LAW")
    print("="*78)

    degree_data = []

    for transition in sorted(
        local_ops
    ):

        degree_bound = local_ops[
            transition
        ][0]

        degree_data.append(
            (
                transition,
                degree_bound,
            )
        )

        print(
            f"  transition={transition}: "
            f"degree(M)={degree_bound}"
        )

    expected_linear = all(
        degree_bound
        ==
        transition
        for transition, degree_bound
        in degree_data
        if transition >= 1
    )

    print()
    print(
        "  degree(M_j) = j on all available transitions =",
        expected_linear,
    )

    # ------------------------------------------------------------------------
    # 4. DET / TRACE
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("4. DETERMINANT / TRACE PROFILE")
    print("="*78)

    determinants = {}
    traces = {}

    for transition in sorted(
        local_ops
    ):

        M = local_ops[
            transition
        ][1]

        dM = determinant(M)
        trM = trace(M)

        determinants[transition] = dM
        traces[transition] = trM

        print()
        print(
            f"  transition={transition}"
        )

        print(
            "    det degree =",
            degree_u(dM),
        )

        print(
            "    det =",
            dM,
        )

        print(
            "    trace degree =",
            degree_u(trM),
        )

        print(
            "    trace =",
            trM,
        )

    # ------------------------------------------------------------------------
    # 5. MATRIX ENTRY COEFFICIENT DATA
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("5. MATRIX ENTRY COEFFICIENT SEQUENCES")
    print("="*78)

    maximum_operator_degree = max(
        [
            int(local_ops[k][0])
            for k in local_ops
        ]
        or [0]
    )

    for r in range(2):

        for c in range(2):

            print()
            print(
                f"  ENTRY M_j[{r},{c}]"
            )

            for power in range(
                maximum_operator_degree+1
            ):

                indices = []
                values = []

                for transition in sorted(
                    local_ops
                ):

                    M = local_ops[
                        transition
                    ][1]

                    value = matrix_coeff(
                        M,
                        r,
                        c,
                        power,
                    )

                    indices.append(
                        transition
                    )

                    values.append(
                        value
                    )

                if all(
                    v == 0
                    for v in values
                ):
                    continue

                print(
                    f"    u^{power}:"
                )

                print(
                    "      transitions =",
                    indices,
                )

                print(
                    "      values =",
                    values,
                )

    # ------------------------------------------------------------------------
    # 6. FINITE DIFFERENCE DEGREE
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("6. EXACT FINITE-DIFFERENCE DEGREE IN THE INDEX")
    print("="*78)

    for r in range(2):

        for c in range(2):

            print()
            print(
                f"  ENTRY M_j[{r},{c}]"
            )

            for power in range(
                maximum_operator_degree+1
            ):

                values = []

                for transition in sorted(
                    local_ops
                ):

                    M = local_ops[
                        transition
                    ][1]

                    values.append(
                        matrix_coeff(
                            M,
                            r,
                            c,
                            power,
                        )
                    )

                if all(
                    v == 0
                    for v in values
                ):
                    continue

                dd = difference_degree(
                    values
                )

                print(
                    f"    u^{power}: "
                    f"difference_degree={dd}"
                )

    # ------------------------------------------------------------------------
    # 7. EXACT FINITE-RANGE POLYNOMIAL-IN-j DESCRIPTIONS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("7. EXACT FINITE-RANGE POLYNOMIAL-IN-j TEST")
    print("="*78)

    print(
        """
  For each coefficient [u^k] of every operator entry, construct the
  exact polynomial in the transition index j that matches every
  available transition.

  This is an exact finite-range identity only.

  It is NOT interpreted as a universal recurrence beyond the observed
  transition range.
        """
    )

    polynomial_descriptions = []

    for r in range(2):

        for c in range(2):

            for power in range(
                maximum_operator_degree+1
            ):

                indices = []
                values = []

                for transition in sorted(
                    local_ops
                ):

                    value = matrix_coeff(
                        local_ops[
                            transition
                        ][1],
                        r,
                        c,
                        power,
                    )

                    if value != 0:
                        indices.append(
                            transition
                        )
                        values.append(
                            value
                        )
                    else:
                        indices.append(
                            transition
                        )
                        values.append(
                            sp.Integer(0)
                        )

                if all(
                    v == 0
                    for v in values
                ):
                    continue

                candidate = interpolate_j(
                    values,
                    indices,
                )

                verified = verify_j_identity(
                    candidate,
                    values,
                    indices,
                )

                candidate_degree = sp.Poly(
                    candidate,
                    j,
                    domain=sp.QQ,
                ).degree()

                polynomial_descriptions.append(
                    (
                        r,
                        c,
                        power,
                        candidate_degree,
                        candidate,
                    )
                )

                print()
                print(
                    f"  M[{r},{c}] coefficient u^{power}:"
                )

                print(
                    "    j-degree =",
                    candidate_degree,
                )

                print(
                    "    exact on all transitions =",
                    verified,
                )

                print(
                    "    polynomial =",
                    candidate,
                )

    # ------------------------------------------------------------------------
    # 8. DET / TRACE INDEX POLYNOMIALS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("8. DETERMINANT / TRACE AS FUNCTIONS OF THE INDEX")
    print("="*78)

    for label, data in [
        (
            "det",
            determinants,
        ),
        (
            "trace",
            traces,
        ),
    ]:

        indices = sorted(
            data
        )

        values = [
            data[k]
            for k in indices
        ]

        if not values:
            continue

        # Each value is a QQ[u] polynomial.
        # First test coefficientwise in u.
        max_u_degree = max(
            degree_u(v)
            for v in values
        )

        print()
        print(
            f"  {label}:"
        )

        for power in range(
            max_u_degree+1
        ):

            coefficient_values = [
                coeff_u(
                    data[k],
                    power,
                )
                for k in indices
            ]

            if all(
                v == 0
                for v in coefficient_values
            ):
                continue

            candidate = interpolate_j(
                coefficient_values,
                indices,
            )

            verified = verify_j_identity(
                candidate,
                coefficient_values,
                indices,
            )

            print(
                f"    u^{power}: "
                f"j-degree={sp.Poly(candidate, j, domain=sp.QQ).degree()} "
                f"verified={verified}"
            )

            print(
                "      =",
                candidate,
            )

    # ------------------------------------------------------------------------
    # 9. NORMALIZATION OF LOCAL OPERATORS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("9. LOCAL OPERATOR CONTENT / DENOMINATOR AUDIT")
    print("="*78)

    for transition in sorted(
        local_ops
    ):

        M = local_ops[
            transition
        ][1]

        numerators = []
        denominators = []

        for r in range(2):

            for c in range(2):

                d = degree_u(
                    M[r,c]
                )

                if d == -sp.oo:
                    continue

                for power in range(
                    int(d)+1
                ):

                    value = matrix_coeff(
                        M,
                        r,
                        c,
                        power,
                    )

                    value = sp.cancel(
                        value
                    )

                    numerators.append(
                        sp.Integer(
                            sp.numer(value)
                        )
                    )

                    denominators.append(
                        sp.Integer(
                            sp.denom(value)
                        )
                    )

        if numerators:

            content = abs(
                numerators[0]
            )

            for value in numerators[1:]:
                content = sp.igcd(
                    content,
                    abs(value),
                )

            denominator_lcm = abs(
                denominators[0]
            )

            for value in denominators[1:]:
                denominator_lcm = sp.ilcm(
                    denominator_lcm,
                    abs(value),
                )

            print()
            print(
                f"  M_{transition}:"
            )

            print(
                "    numerator content =",
                content,
            )

            print(
                "    denominator LCM =",
                denominator_lcm,
            )

    # ------------------------------------------------------------------------
    # 10. COMPOSITION
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("10. EXACT CONSECUTIVE OPERATOR COMPOSITION")
    print("="*78)

    ordered = sorted(
        local_ops
    )

    for idx in range(
        1,
        len(ordered),
    ):

        first = ordered[idx-1]
        second = ordered[idx]

        M1 = local_ops[first][1]
        M2 = local_ops[second][1]

        composed = sp.Matrix([
            [
                E(
                    M2[r,0]*M1[0,c]
                    +
                    M2[r,1]*M1[1,c]
                )
                for c in range(2)
            ]
            for r in range(2)
        ])

        print()
        print(
            f"  M_{second} * M_{first}:"
        )

        print(
            "    degree =",
            matrix_degree(
                composed
            ),
        )

        print(
            "    determinant degree =",
            degree_u(
                determinant(
                    composed
                )
            ),
        )

        print(
            "    trace degree =",
            degree_u(
                trace(
                    composed
                )
            ),
        )

    # ------------------------------------------------------------------------
    # 11. FRESH EXACT SEMIPRIME
    # ------------------------------------------------------------------------

    p = sp.Integer(
        106621
    )

    q = sp.Integer(
        246473
    )

    if not sp.isprime(p):
        raise RuntimeError(
            "Fresh p is not prime."
        )

    if not sp.isprime(q):
        raise RuntimeError(
            "Fresh q is not prime."
        )

    N = E(p*q)
    S = E(p+q)
    X = E(S+1)

    t_value = sp.cancel(
        sp.Rational(
            int(N),
            int(X),
        )
    )

    u_value = sp.cancel(
        t_value*(t_value+1)
    )

    print()
    print("="*78)
    print("11. FRESH EXACT SEMIPRIME")
    print("="*78)

    print(
        "  p =",
        p,
    )

    print(
        "  q =",
        q,
    )

    print(
        "  N =",
        N,
    )

    print(
        "  S =",
        S,
    )

    print(
        "  X =",
        X,
    )

    print(
        "  t =",
        t_value,
    )

    print(
        "  u =",
        u_value,
    )

    # ------------------------------------------------------------------------
    # 12. FRESH LOCAL OPERATOR VALIDATION
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("12. FRESH LOCAL OPERATOR VALIDATION")
    print("="*78)

    fresh_vectors = []

    for vec in V:

        fresh_vectors.append(
            sp.Matrix([
                E(
                    vec[0].subs(
                        u,
                        u_value,
                    )
                ),
                E(
                    vec[1].subs(
                        u,
                        u_value,
                    )
                ),
            ])
        )

    local_ok = True

    for transition in sorted(
        local_ops
    ):

        M = local_ops[
            transition
        ][1]

        if M.free_symbols - {u}:
            print(
                f"  V_{transition-1} -> "
                f"V_{transition}: parameterized"
            )
            continue

        predicted = (
            M.subs(
                u,
                u_value,
            )
            *
            fresh_vectors[
                transition-1
            ]
        )

        target = fresh_vectors[
            transition
        ]

        ok = all(
            is_zero(
                E(
                    predicted[k]
                    -
                    target[k]
                )
            )
            for k in range(2)
        )

        local_ok &= ok

        print(
            f"  V_{transition-1} -> "
            f"V_{transition}:",
            ok,
        )

    # ------------------------------------------------------------------------
    # 13. FRESH CHANNEL RECONSTRUCTION
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("13. FRESH CHANNEL RECONSTRUCTION")
    print("="*78)

    channel_ok = True

    for d in range(
        16,
        -1,
        -1,
    ):

        A_d, B_d = CHANNEL[d]

        predicted = E(
            A_d.subs(
                u,
                u_value,
            )
            +
            (2*t_value+1)
            *
            B_d.subs(
                u,
                u_value,
            )
        )

        observed = E(
            H[d].subs(
                t,
                t_value,
            )
        )

        ok = is_zero(
            predicted-observed
        )

        channel_ok &= ok

        print(
            f"  d={d:2d}: "
            f"reconstruction={ok}"
        )

    # ------------------------------------------------------------------------
    # 14. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("14. STRUCTURAL INTERPRETATION")
    print("="*78)

    print(
        r"""
  Experiment 78 found exact coupled local operators

      V_j(u) = M_j(u) V_{j-1}(u)

  and the minimal polynomial degrees increased with the transition.

  Experiment 79 now asks whether this apparent degree growth is itself
  governed by a simple index law.

  The principal target is

      M_j(u) = M(j,u).

  This would mean that each of the four matrix entries belongs to a
  low-complexity two-variable rational/polynomial family.

  We separately examine

      det M_j,
      tr  M_j,

  because these scalar invariants can reveal a simpler law even when
  the individual matrix entries are complicated.

  The strongest possible structural outcome would be:

      M_j(u) =
        [ P_00(j,u)  P_01(j,u) ]
        [ P_10(j,u)  P_11(j,u) ]

  with all P_rc having low degree in both j and u.

  That would give an explicit construction mechanism for the coupled
  channel states.

  The important distinction remains:

      INVERSE:
          observed layer values
             -> t
             -> X
             -> N
             -> p,q

      CONSTRUCTION:
          V_0
             -> M_1(j,u)
             -> M_2(j,u)
             -> ...
             -> V_7.

  Even a successful M_j law would still not establish

      N -> layers

  because the initial layer/state and the parameter u remain tied to
  the hidden configuration.

  The experiment therefore targets the algebraic construction mechanism,
  not a factorization claim.
        """
    )

    # ------------------------------------------------------------------------
    # 15. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("="*78)
    print("15. FINAL EXACTNESS")
    print("="*78)

    failures = 0

    if not involution_ok:
        failures += 1

    if not channel_ok:
        failures += 1

    if not local_ok:
        failures += 1

    print(
        "  involution =",
        involution_ok,
    )

    print(
        "  channel_reconstruction =",
        channel_ok,
    )

    print(
        "  local_operator_validation =",
        local_ok,
    )

    print(
        "  local_operator_count =",
        len(local_ops),
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
    print(
        "EXPERIMENT 79 COMPLETE"
    )
    print("="*78)


if __name__ == "__main__":
    main()