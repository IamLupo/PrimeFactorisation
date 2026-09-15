#!/usr/bin/env python3

"""
EXPERIMENT 78 — EXACT COUPLED A/B TRANSFER-OPERATOR SEARCH
======================================================================

Structural setup
----------------
    t = N/X
    u = t(t+1)

The normalized homogeneous layers satisfy

    h_d(t) = A_d(u) + (2t+1) B_d(u)

with

    even d : B_d = 0
    odd  d : A_d = 0.

We group the channels as

    V_j(u) =
        [ A_{2j}(u)   ]
        [ B_{2j+1}(u) ]

and investigate whether

    V_j = M_j(u) V_{j-1}

or, more strongly,

    V_j = M(u) V_{j-1}

for a universal 2x2 polynomial matrix M(u).

IMPORTANT IMPLEMENTATION RULE
------------------------------
SymPy must never be asked to construct

    Poly(expr, u, domain=QQ)

when expr contains unresolved operator coefficients.

That was the source of the previous failures.

This script therefore:

1. extracts u-coefficients directly from symbolic expressions;
2. solves the resulting exact linear systems;
3. distinguishes unique from parameterized solutions;
4. uses a coefficient-based degree routine that works even for
   expressions containing unresolved parameters;
5. only performs QQ-polynomial operations after all parameters have
   been eliminated.

No floating point.
No regression.
No interpolation.
No statistical inference.
"""


import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

t = sp.Symbol("t")
u = sp.Symbol("u")


# ============================================================================
# EXACT NORMALIZED LAYER POLYNOMIALS
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
        -276*t**10 - 1380*t**9 - 5460*t**8
        -13560*t**7 -23058*t**6 -27510*t**5
        -23100*t**4 -13410*t**3 -5135*t**2
        -1169*t -120
    ),

    13: (
        288*t**11 + 1584*t**10 + 10340*t**9 + 34650*t**8
        + 79464*t**7 + 127512*t**6 + 145530*t**5
        + 117975*t**4 + 66550*t**3 + 24882*t**2
        + 5551*t + 560
    ),

    12: (
        -50*t**12 -300*t**11 -7436*t**10 -34430*t**9
        -117810*t**8 -267960*t**7 -426888*t**6
        -484110*t**5 -390225*t**4 -219010*t**3
        -81510*t**2 -18109*t -1820
    ),

    11: (
        1300*t**11 +7150*t**10 +62634*t**9 +228228*t**8
        +552552*t**7 +918918*t**6 +1075074*t**5
        +887172*t**4 +507078*t**3 +191477*t**2
        +43043*t +4368
    ),

    10: (
        -10010*t**10 -50050*t**9 -252252*t**8
        -708708*t**7 -1303302*t**6 -1639638*t**5
        -1429428*t**4 -852852*t**3 -333333*t**2
        -77077*t -8008
    ),

    9: (
        35750*t**9 +160875*t**8 +563420*t**7
        +1221220*t**6 +1732458*t**5 +1653470*t**4
        +1058540*t**3 +437700*t**2 +105979*t +11441
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
        -331704*t**3 -190296*t**2 -60044*t -8092
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
# SAFE SYMBOLIC HELPERS
# ============================================================================

def E(x):
    return sp.expand(sp.sympify(x))


def zero(x):
    return E(x) == 0


def safe_u_degree(x):
    """
    Degree in u without requiring coefficients to belong to QQ.

    This is the crucial repair.

    For an expression such as

        25*c_0_1_0/2 + 78

    the routine simply inspects the powers of u appearing in the expanded
    expression instead of asking SymPy to construct Poly(..., QQ).
    """
    x = E(x)

    if x == 0:
        return -sp.oo

    degree = None

    for term in sp.Add.make_args(x):

        powers = term.as_powers_dict()

        exponent = powers.get(u, sp.Integer(0))

        if exponent.is_integer:
            exponent_int = int(exponent)

            if degree is None or exponent_int > degree:
                degree = exponent_int

    return -sp.oo if degree is None else degree


def safe_u_coeff(x, k):
    """
    Exact coefficient of u^k.

    This works with symbolic parameters.
    """
    x = E(x)
    result = sp.Integer(0)

    for term in sp.Add.make_args(x):

        powers = term.as_powers_dict()

        exponent = powers.get(u, sp.Integer(0))

        if exponent == k:
            result += term / u**k

    return sp.simplify(result)


def all_free_symbols(x):
    return (
        E(x).free_symbols
        - {u, t}
    )


def has_unresolved_parameters(x):
    return bool(all_free_symbols(x))


def matrix_has_parameters(M):
    for i in range(M.rows):
        for j in range(M.cols):
            if has_unresolved_parameters(M[i, j]):
                return True
    return False


def matrix_degree_safe(M):
    """
    Degree of a matrix in u, without QQ coercion.
    """
    values = []

    for i in range(M.rows):
        for j in range(M.cols):
            values.append(
                safe_u_degree(M[i, j])
            )

    finite = [
        x for x in values
        if x != -sp.oo
    ]

    return -sp.oo if not finite else max(finite)


def exact_matrix_zero(M):
    return all(
        zero(M[i, j])
        for i in range(M.rows)
        for j in range(M.cols)
    )


def matrix_substitute(M, variable, value):
    return sp.Matrix([
        [
            E(M[i, j].subs(variable, value))
            for j in range(M.cols)
        ]
        for i in range(M.rows)
    ])


def matrices_equal(M1, M2):
    if M1.shape != M2.shape:
        return False

    return all(
        zero(M1[i, j] - M2[i, j])
        for i in range(M1.rows)
        for j in range(M1.cols)
    )


def vector_equal(v1, v2):
    return (
        v1.shape == v2.shape
        and all(
            zero(v1[i] - v2[i])
            for i in range(v1.rows)
            for _ in [0]
        )
    )


# ============================================================================
# RECOVER A POLYNOMIAL IN u FROM A POLYNOMIAL IN t
# ============================================================================

def recover_u_polynomial(target, prefix):
    """
    Find P(u) such that

        P(t(t+1)) = target(t)

    exactly.
    """

    target = E(target)

    if target == 0:
        return sp.Integer(0)

    target_poly = sp.Poly(
        target,
        t,
        domain=sp.QQ,
    )

    deg_t = target_poly.degree()
    deg_u = deg_t // 2

    coeffs = sp.symbols(
        f"{prefix}_0:{deg_u + 1}"
    )

    candidate = E(
        sum(
            coeffs[k] * u**k
            for k in range(deg_u + 1)
        )
    )

    substituted = E(
        candidate.subs(
            u,
            t*(t+1),
        )
    )

    difference = sp.Poly(
        substituted - target,
        t,
    )

    equations = [
        sp.expand(c)
        for c in difference.all_coeffs()
        if c != 0
    ]

    solution = sp.solve(
        equations,
        coeffs,
        dict=True,
    )

    if not solution:
        raise RuntimeError(
            f"Unable to recover {prefix}(u)"
        )

    solved = E(
        candidate.subs(
            solution[0]
        )
    )

    if not zero(
        solved.subs(u, t*(t+1)) - target
    ):
        raise RuntimeError(
            f"Exact {prefix}(u) reconstruction failed"
        )

    return solved


# ============================================================================
# INVOLUTION DECOMPOSITION
# ============================================================================

def decompose_layer(d):

    h = H[d]

    transformed = E(
        h.subs(
            t,
            -1-t,
        )
    )

    even_part = E(
        (h + transformed) / 2
    )

    odd_part = E(
        (h - transformed) / 2
    )

    A_d = recover_u_polynomial(
        even_part,
        f"A_{d}",
    )

    if odd_part == 0:
        B_d = sp.Integer(0)

    else:

        numerator = sp.Poly(
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
            numerator,
            divisor,
        )

        if not remainder.is_zero:
            raise RuntimeError(
                f"d={d}: odd layer not divisible by 2t+1"
            )

        B_d = recover_u_polynomial(
            quotient.as_expr(),
            f"B_{d}",
        )

    reconstruction = E(
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

    if not zero(
        reconstruction - h
    ):
        raise RuntimeError(
            f"d={d}: reconstruction failed"
        )

    return A_d, B_d


CHANNEL = {
    d: decompose_layer(d)
    for d in range(16, -1, -1)
}


A = {
    d: CHANNEL[d][0]
    for d in range(0, 17, 2)
}

B = {
    d: CHANNEL[d][1]
    for d in range(1, 16, 2)
}


# ============================================================================
# COUPLED STATE VECTORS
# ============================================================================

V = []

for j in range(8):

    V.append(
        sp.Matrix([
            E(A[2*j]),
            E(B[2*j+1]),
        ])
    )


# ============================================================================
# OPERATOR UNKNOWN GENERATION
# ============================================================================

def create_operator(degree_bound, prefix):

    unknowns = []

    M = sp.zeros(2, 2)

    for i in range(2):
        for j in range(2):

            coeffs = sp.symbols(
                f"{prefix}_{i}_{j}_0:{degree_bound+1}"
            )

            unknowns.extend(
                coeffs
            )

            M[i, j] = E(
                sum(
                    coeffs[k] * u**k
                    for k in range(degree_bound + 1)
                )
            )

    return M, unknowns


# ============================================================================
# BUILD EXACT LINEAR SYSTEM
# ============================================================================

def operator_equations(
    left_vectors,
    right_vectors,
    degree_bound,
    prefix,
):

    M, unknowns = create_operator(
        degree_bound,
        prefix,
    )

    equations = []

    for left, right in zip(
        left_vectors,
        right_vectors,
    ):

        predicted = sp.Matrix([
            E(
                M[0, 0]*right[0]
                +
                M[0, 1]*right[1]
            ),
            E(
                M[1, 0]*right[0]
                +
                M[1, 1]*right[1]
            ),
        ])

        for component in range(2):

            diff = E(
                predicted[component]
                -
                left[component]
            )

            # We know all actual channel polynomials have small u-degree.
            # The product can reach degree degree_bound + 6.
            max_degree = degree_bound + 7

            for k in range(max_degree + 1):

                coefficient = safe_u_coeff(
                    diff,
                    k,
                )

                if coefficient != 0:
                    equations.append(
                        coefficient
                    )

    return M, unknowns, equations


# ============================================================================
# SOLVE OPERATOR
# ============================================================================

def solve_operator(
    left_vectors,
    right_vectors,
    degree_bound,
    prefix,
):

    M_symbolic, unknowns, equations = (
        operator_equations(
            left_vectors,
            right_vectors,
            degree_bound,
            prefix,
        )
    )

    if not equations:
        return {
            "exists": False,
            "rank": 0,
            "unknowns": len(unknowns),
            "nullity": len(unknowns),
            "matrix": None,
            "unique": False,
        }

    Aeq, beq = sp.linear_eq_to_matrix(
        equations,
        unknowns,
    )

    rank = Aeq.rank()

    n_unknowns = len(unknowns)

    nullity = n_unknowns - rank

    solutions = sp.linsolve(
        (Aeq, beq),
        unknowns,
    )

    if solutions == sp.EmptySet:

        return {
            "exists": False,
            "rank": rank,
            "unknowns": n_unknowns,
            "nullity": nullity,
            "matrix": None,
            "unique": False,
        }

    tuples = list(solutions)

    if not tuples:

        return {
            "exists": False,
            "rank": rank,
            "unknowns": n_unknowns,
            "nullity": nullity,
            "matrix": None,
            "unique": False,
        }

    solution_tuple = tuples[0]

    substitution = dict(
        zip(
            unknowns,
            solution_tuple,
        )
    )

    solved_matrix = sp.Matrix([
        [
            E(
                M_symbolic[i, j]
                .subs(substitution)
            )
            for j in range(2)
        ]
        for i in range(2)
    ])

    unresolved = set()

    for entry in solution_tuple:

        unresolved |= (
            entry.free_symbols
            & set(unknowns)
        )

    unique = (
        len(unresolved) == 0
    )

    return {
        "exists": True,
        "rank": rank,
        "unknowns": n_unknowns,
        "nullity": nullity,
        "matrix": solved_matrix,
        "unique": unique,
    }


# ============================================================================
# SAFE MATRIX INVARIANTS
# ============================================================================

def safe_matrix_determinant(M):

    return sp.factor(
        E(
            M[0, 0]*M[1, 1]
            -
            M[0, 1]*M[1, 0]
        )
    )


def safe_matrix_trace(M):

    return sp.factor(
        E(
            M[0, 0] + M[1, 1]
        )
    )


def matrix_parameter_count(M):

    symbols = set()

    for i in range(M.rows):
        for j in range(M.cols):
            symbols |= all_free_symbols(
                M[i, j]
            )

    return len(symbols)


# ============================================================================
# FRESH VALIDATION
# ============================================================================

def validate_operator_at_value(
    M,
    left,
    right,
    u_value,
):

    if M is None:
        return None

    if matrix_has_parameters(M):
        return None

    M_eval = matrix_substitute(
        M,
        u,
        u_value,
    )

    right_eval = sp.Matrix([
        E(
            right[i].subs(
                u,
                u_value,
            )
        )
        for i in range(2)
    ])

    left_eval = sp.Matrix([
        E(
            left[i].subs(
                u,
                u_value,
            )
        )
        for i in range(2)
    ])

    predicted = M_eval * right_eval

    return all(
        zero(
            predicted[i]
            -
            left_eval[i]
        )
        for i in range(2)
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("="*78)
    print(
        "EXPERIMENT 78 — EXACT COUPLED A/B TRANSFER-OPERATOR SEARCH"
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
    print("  coupled state V_j = [A_{2j}, B_{2j+1}]^T")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # ----------------------------------------------------------------------
    # 1 INVOLUTION
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("1. EXACT INVOLUTION VALIDATION")
    print("="*78)

    involution_ok = True

    for d in range(16, -1, -1):

        h = H[d]

        partner = E(
            h.subs(
                t,
                -1-t,
            )
        )

        invariant = zero(
            partner - h
        )

        anti = zero(
            partner + h
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

        ok = zero(
            reconstructed - h
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

    # ----------------------------------------------------------------------
    # 2 CHANNEL PROFILE
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("2. COUPLED CHANNEL PROFILE")
    print("="*78)

    for d in range(16, -1, -1):

        A_d, B_d = CHANNEL[d]

        print()
        print(
            f"  d={d}:"
        )

        print(
            "    A_d =",
            sp.factor(A_d),
        )

        print(
            "    B_d =",
            sp.factor(B_d),
        )

    # ----------------------------------------------------------------------
    # 3 COUPLED VECTORS
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("3. COUPLED STATE VECTORS")
    print("="*78)

    for j, vec in enumerate(V):

        print()
        print(
            f"  V_{j} = [A_{2*j}, B_{2*j+1}]^T"
        )

        print(
            "    A =",
            sp.factor(vec[0]),
        )

        print(
            "    B =",
            sp.factor(vec[1]),
        )

    # ----------------------------------------------------------------------
    # 4 UNIVERSAL OPERATOR SEARCH
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("4. UNIVERSAL 2x2 POLYNOMIAL OPERATOR SEARCH")
    print("="*78)

    universal = []

    for degree_bound in range(0, 8):

        result = solve_operator(
            V[1:],
            V[:-1],
            degree_bound,
            f"U{degree_bound}",
        )

        if not result["exists"]:

            print(
                f"  degree={degree_bound}: NONE "
                f"(rank={result['rank']}, "
                f"unknowns={result['unknowns']})"
            )

            continue

        print(
            f"  degree={degree_bound}: EXISTS "
            f"(rank={result['rank']}, "
            f"unknowns={result['unknowns']}, "
            f"nullity={result['nullity']}, "
            f"unique={result['unique']})"
        )

        M = result["matrix"]

        if M is not None:

            print(
                "    matrix degree =",
                matrix_degree_safe(M),
            )

            print(
                "    free parameters =",
                matrix_parameter_count(M),
            )

            print(
                "    M(u) ="
            )

            print(M)

        universal.append(
            (
                degree_bound,
                result,
            )
        )

    # ----------------------------------------------------------------------
    # 5 FIRST UNIQUE UNIVERSAL OPERATOR
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("5. FIRST UNIQUE UNIVERSAL OPERATOR")
    print("="*78)

    found_unique = False

    for degree_bound, result in universal:

        if result["exists"] and result["unique"]:

            M = result["matrix"]

            print(
                f"  UNIQUE operator found at degree={degree_bound}"
            )

            print(
                "  M(u) ="
            )

            print(M)

            print(
                "  det(M) =",
                safe_matrix_determinant(M),
            )

            print(
                "  trace(M) =",
                safe_matrix_trace(M),
            )

            found_unique = True

            break

    if not found_unique:

        print(
            "  No unique universal operator found through degree 7."
        )

    # ----------------------------------------------------------------------
    # 6 LOCAL OPERATOR SEARCH
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("6. LOCAL 2x2 OPERATOR SEARCH")
    print("="*78)

    local = []

    for j in range(1, len(V)):

        print()
        print(
            f"  Transition V_{j-1} -> V_{j}"
        )

        selected = None

        for degree_bound in range(0, 8):

            result = solve_operator(
                [V[j]],
                [V[j-1]],
                degree_bound,
                f"L{j}_{degree_bound}",
            )

            if not result["exists"]:

                print(
                    f"    degree={degree_bound}: NONE"
                )

                continue

            print(
                f"    degree={degree_bound}: EXISTS "
                f"rank={result['rank']} "
                f"unknowns={result['unknowns']} "
                f"nullity={result['nullity']} "
                f"unique={result['unique']}"
            )

            M = result["matrix"]

            if M is not None:

                print(
                    "      degree =",
                    matrix_degree_safe(M),
                )

                print(
                    "      free parameters =",
                    matrix_parameter_count(M),
                )

            selected = (
                degree_bound,
                result,
            )

            # Prefer the first concrete solution.
            break

        if selected is None:

            print(
                "    No operator found in scanned range."
            )

        else:

            degree_bound, result = selected

            M = result["matrix"]

            if M is not None:

                print(
                    "    selected operator:"
                )

                print(M)

                print(
                    "    determinant =",
                    safe_matrix_determinant(M),
                )

                print(
                    "    trace =",
                    safe_matrix_trace(M),
                )

                print(
                    "    upper triangular =",
                    zero(M[1, 0]),
                )

                print(
                    "    lower triangular =",
                    zero(M[0, 1]),
                )

            else:

                print(
                    "    selected solution remains parameterized."
                )

        local.append(
            (
                j,
                selected,
            )
        )

    # ----------------------------------------------------------------------
    # 7 LOCAL OPERATOR PARAMETER AUDIT
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("7. LOCAL OPERATOR PARAMETER AUDIT")
    print("="*78)

    concrete_local = []

    for j, selected in local:

        if selected is None:
            continue

        degree_bound, result = selected

        M = result["matrix"]

        if M is None:
            print(
                f"  V_{j-1} -> V_{j}: "
                "parameterized / no concrete representative"
            )
            continue

        params = matrix_parameter_count(M)

        print(
            f"  V_{j-1} -> V_{j}: "
            f"degree={degree_bound} "
            f"parameters={params}"
        )

        if params == 0:

            concrete_local.append(
                (
                    j,
                    degree_bound,
                    M,
                )
            )

    # ----------------------------------------------------------------------
    # 8 LOCAL OPERATOR SPACE
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("8. CONCRETE LOCAL OPERATOR-SPACE RANK")
    print("="*78)

    if concrete_local:

        max_degree = max(
            matrix_degree_safe(M)
            for _, _, M in concrete_local
        )

        rows = []

        for _, _, M in concrete_local:

            row = []

            for i in range(2):
                for j in range(2):

                    for k in range(
                        max_degree + 1
                    ):

                        row.append(
                            safe_u_coeff(
                                M[i, j],
                                k,
                            )
                        )

            rows.append(row)

        rank = sp.Matrix(rows).rank()

        print(
            "  concrete operators =",
            len(concrete_local),
        )

        print(
            "  coefficient-space dimension =",
            4*(max_degree + 1),
        )

        print(
            "  rank =",
            rank,
        )

        print(
            "  nullity =",
            len(concrete_local) - rank,
        )

    else:

        print(
            "  no fully determined local operators"
        )

    # ----------------------------------------------------------------------
    # 9 CONSECUTIVE OPERATOR COMPARISON
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("9. CONSECUTIVE LOCAL OPERATOR COMPARISON")
    print("="*78)

    if len(concrete_local) >= 2:

        for k in range(
            1,
            len(concrete_local),
        ):

            j0, _, M0 = concrete_local[k-1]
            j1, _, M1 = concrete_local[k]

            D = sp.Matrix([
                [
                    E(
                        M1[i, j]
                        -
                        M0[i, j]
                    )
                    for j in range(2)
                ]
                for i in range(2)
            ])

            print()
            print(
                f"  Delta M: "
                f"transition {j0-1}->{j0} "
                f"to {j1-1}->{j1}"
            )

            print(
                "    matrix ="
            )

            print(D)

            print(
                "    degree =",
                matrix_degree_safe(D),
            )

            print(
                "    determinant =",
                safe_matrix_determinant(D),
            )

            print(
                "    trace =",
                safe_matrix_trace(D),
            )

    else:

        print(
            "  insufficient concrete operators"
        )

    # ----------------------------------------------------------------------
    # 10 FRESH SEMIPRIME
    # ----------------------------------------------------------------------

    p = sp.Integer(106621)
    q = sp.Integer(246473)

    if not sp.isprime(p):
        raise RuntimeError(
            "Fresh p is not prime."
        )

    if not sp.isprime(q):
        raise RuntimeError(
            "Fresh q is not prime."
        )

    N = p*q
    S = p+q
    X = S+1

    t_value = sp.cancel(
        sp.Rational(N, X)
    )

    u_value = sp.cancel(
        t_value*(t_value+1)
    )

    print()
    print("="*78)
    print("10. FRESH EXACT SEMIPRIME")
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

    # ----------------------------------------------------------------------
    # 11 FRESH VALIDATION OF DETERMINED LOCAL OPERATORS
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("11. FRESH LOCAL OPERATOR VALIDATION")
    print("="*78)

    evaluated_vectors = []

    for vec in V:

        evaluated_vectors.append(
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

    local_validation_ok = True

    for j, degree_bound, M in concrete_local:

        result = (
            M.subs(
                u,
                u_value,
            )
            *
            evaluated_vectors[j-1]
        )

        target = evaluated_vectors[j]

        ok = all(
            zero(
                result[i] - target[i]
            )
            for i in range(2)
        )

        local_validation_ok &= ok

        print(
            f"  V_{j-1} -> V_{j}:",
            ok,
        )

    if not concrete_local:

        print(
            "  no determined operators available for fresh validation"
        )

    # ----------------------------------------------------------------------
    # 12 FRESH CHANNEL RECONSTRUCTION
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("12. FRESH CHANNEL RECONSTRUCTION")
    print("="*78)

    channel_ok = True

    for d in range(
        16,
        -1,
        -1,
    ):

        A_d, B_d = CHANNEL[d]

        actual = E(
            H[d].subs(
                t,
                t_value,
            )
        )

        reconstructed = E(
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

        ok = (
            actual == reconstructed
        )

        channel_ok &= ok

        print(
            f"  d={d:2d}: reconstruction={ok}"
        )

    # ----------------------------------------------------------------------
    # 13 STRUCTURAL CONCLUSION
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("13. STRUCTURAL INTERPRETATION")
    print("="*78)

    print(
        """
  The experiment now searches for a genuinely coupled construction law.

  The state is

      V_j(u) =
          [ A_{2j}(u)   ]
          [ B_{2j+1}(u) ].

  We test

      V_j = M(u) V_{j-1}

  both globally and locally.

  There are three distinct outcomes:

    1. UNIQUE:
         the operator is completely determined by exact algebra.

    2. PARAMETERIZED:
         the equations admit a family of operators.
         Such a family is NOT treated as a discovered construction law.

    3. NONE:
         no polynomial operator within the tested degree exists.

  The coupled formulation is stronger than the earlier scalar search
  because the even and odd involution channels are allowed to interact.

  The next structural step, if no small polynomial operator exists,
  is to search for rational-function operators or operators whose
  coefficients depend explicitly on the layer index.

  None of these experiments yet compute the layer values from N alone.
  The established direction remains

      observed layers -> t -> X -> N -> p,q.

  The unresolved direction remains

      N -> layer information -> t -> X -> p,q.
        """
    )

    # ----------------------------------------------------------------------
    # 14 FINAL AUDIT
    # ----------------------------------------------------------------------

    print()
    print("="*78)
    print("14. FINAL EXACTNESS")
    print("="*78)

    failures = 0

    if not involution_ok:
        failures += 1

    if not channel_ok:
        failures += 1

    if concrete_local:
        if not local_validation_ok:
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
        "  concrete_local_operator_validation =",
        (
            local_validation_ok
            if concrete_local
            else "not_applicable"
        ),
    )

    print(
        "  universal_operator_found =",
        found_unique,
    )

    print(
        "  concrete_local_operators =",
        len(concrete_local),
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
    print("EXPERIMENT 78 COMPLETE")
    print("="*78)


if __name__ == "__main__":
    main()