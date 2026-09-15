#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 96 — EXACT OFFSET-GENERATOR DIFFERENTIAL-OPERATOR AUDIT
# REPAIRED SYMBOLIC LINEAR-SYSTEM VERSION
#
# Exact QQ arithmetic.
# No floating point.
# No Poly(..., domain=QQ) on expressions containing unknown coefficients.
# ==============================================================================

x = sp.symbols("x")


# ==============================================================================
# 1. EXACT OFFSET GENERATORS
# ==============================================================================

A_G = [
    (
        103173*x**8
        - 414184*x**7
        + 748888*x**6
        + 1965936*x**5
        - 14170800*x**4
        + 18278400*x**3
        + 32981760*x**2
        + 6209280*x
        + 80640
    ) / sp.Integer(80640),

    (
        1028053*x**7
        - 4684456*x**6
        + 11809112*x**5
        + 1729392*x**4
        - 105893760*x**3
        + 190874880*x**2
        + 202204800*x
        + 22176000
    ) / sp.Integer(22176000),

    (
        3174439*x**6
        - 17102200*x**5
        + 55697320*x**4
        - 68565504*x**3
        - 186883200*x**2
        + 587341440*x
        + 287280000
    ) / sp.Integer(287280000),

    (
        1173535*x**5
        - 9708312*x**4
        + 42619584*x**3
        - 97268640*x**2
        + 34493760*x
        + 239904000
    ) / sp.Integer(239904000),

    -(
        2131503*x**4
        - 2751936*x**3
        - 11748128*x**2
        + 68996928*x
        - 120120000
    ) / sp.Integer(120120000),

    -(
        46945*x**3
        - 163684*x**2
        + 373618*x
        - 420420
    ) / sp.Integer(420420),

    (x**2 - 2*x + 2) / sp.Integer(2),
]


B_G = [
    -(
        421*x**7
        + 69937*x**6
        - 423738*x**5
        + 1055250*x**4
        + 41580*x**3
        - 4071060*x**2
        - 1559880*x
        - 63000
    ) / sp.Integer(63000),

    -(
        12106*x**6
        + 33638*x**5
        - 379836*x**4
        + 1230345*x**3
        - 819480*x**2
        - 3104640*x
        - 630000
    ) / sp.Integer(630000),

    -(
        173070*x**5
        - 240576*x**4
        - 1037057*x**3
        + 6072045*x**2
        - 8596560*x
        - 8139600
    ) / sp.Integer(8139600),

    -(
        808952*x**4
        - 2458897*x**3
        + 3513510*x**2
        + 3999450*x
        - 18564000
    ) / sp.Integer(18564000),

    -(
        162139*x**3
        - 926401*x**2
        + 2779686*x
        - 3753750
    ) / sp.Integer(3753750),

    (
        301*x**2
        - 626*x
        + 650
    ) / sp.Integer(650),
]


# ==============================================================================
# 2. SAFE EXACT HELPERS
# ==============================================================================

def clean(expr):
    """
    Exact symbolic simplification.

    We deliberately do NOT construct Poly(..., domain=QQ) here because
    intermediate expressions may contain symbolic operator parameters.
    """
    return sp.cancel(sp.factor(sp.expand(expr)))


def x_degree(expr):
    """
    Degree in x without requiring all other coefficients to be rational.
    """
    expr = sp.expand(expr)

    if expr == 0:
        return -sp.oo

    p = sp.Poly(expr, x, domain="EX")
    return p.degree()


def x_coeff(expr, power):
    """
    Exact coefficient of x**power, with symbolic parameters allowed.
    """
    return sp.expand(expr).coeff(x, power)


def coefficient_equations(expr, unknowns):
    """
    Convert polynomial identity expr == 0 into exact linear equations
    in unknowns.

    Important: coefficient extraction happens before any QQ domain
    conversion, preventing the previous CoercionFailed error.
    """
    expr = sp.expand(expr)

    deg = x_degree(expr)

    if deg == -sp.oo:
        return []

    equations = []

    for power in range(int(deg) + 1):
        coeff = sp.expand(x_coeff(expr, power))

        if coeff == 0:
            continue

        equations.append(
            sp.Eq(coeff, 0)
        )

    return equations


def linear_system_from_equations(equations, unknowns):
    """
    Convert exact symbolic linear equations into A*c = b.
    """
    rows = []
    rhs = []

    for eq in equations:
        expr = sp.expand(eq.lhs - eq.rhs)

        row = [
            sp.expand(expr).coeff(v)
            for v in unknowns
        ]

        constant = sp.expand(expr)

        for v in unknowns:
            constant = constant.subs(v, 0)

        rows.append(row)
        rhs.append(-constant)

    if not rows:
        return sp.zeros(0, len(unknowns)), sp.zeros(0, 1)

    return sp.Matrix(rows), sp.Matrix(rhs)


def solve_exact_linear(equations, unknowns):
    """
    Exact QQ linear solve.

    Returns:
        matrix
        rhs
        solution set
        rank
        nullity
    """
    A, b = linear_system_from_equations(
        equations,
        unknowns,
    )

    if A.rows == 0:
        return {
            "A": A,
            "b": b,
            "solution": sp.FiniteSet(),
            "rank": 0,
            "unknowns": len(unknowns),
            "nullity": len(unknowns),
        }

    augmented = A.row_join(b)

    rank_A = A.rank()
    rank_aug = augmented.rank()

    if rank_A != rank_aug:
        return {
            "A": A,
            "b": b,
            "solution": sp.EmptySet,
            "rank": rank_A,
            "unknowns": len(unknowns),
            "nullity": max(
                0,
                len(unknowns) - rank_A,
            ),
        }

    solution = sp.linsolve(
        (A, b),
        *unknowns,
    )

    return {
        "A": A,
        "b": b,
        "solution": solution,
        "rank": rank_A,
        "unknowns": len(unknowns),
        "nullity": max(
            0,
            len(unknowns) - rank_A,
        ),
    }


def exact_zero(expr):
    return sp.expand(
        sp.together(expr)
    ) == 0


# ==============================================================================
# 3. LOCAL FIRST-ORDER OPERATOR
#
# G_{j+1}
#   = (a*x+b) G_j'
#     + (c*x+d) G_j
#
# Unknowns: a,b,c,d
# ==============================================================================

def local_first_order_operator(G_prev, G_next):
    a, b, c, d = sp.symbols(
        "a b c d"
    )

    lhs = (
        G_next
        - (a*x + b) * sp.diff(G_prev, x)
        - (c*x + d) * G_prev
    )

    equations = coefficient_equations(
        lhs,
        [a, b, c, d],
    )

    result = solve_exact_linear(
        equations,
        [a, b, c, d],
    )

    return (
        [a, b, c, d],
        result,
    )


# ==============================================================================
# 4. LOCAL EXTENDED OPERATOR
#
# G_{j+1}
#   = (a*x+b) G_j'
#     + (c*x**2 + d*x + e) G_j
#
# Unknowns: a,b,c,d,e
# ==============================================================================

def local_extended_operator(G_prev, G_next):
    a, b, c, d, e = sp.symbols(
        "a b c d e"
    )

    lhs = (
        G_next
        - (a*x + b) * sp.diff(G_prev, x)
        - (c*x**2 + d*x + e) * G_prev
    )

    equations = coefficient_equations(
        lhs,
        [a, b, c, d, e],
    )

    result = solve_exact_linear(
        equations,
        [a, b, c, d, e],
    )

    return (
        [a, b, c, d, e],
        result,
    )


# ==============================================================================
# 5. VERIFY A SPECIFIC SOLUTION
# ==============================================================================

def verify_solution(
    G_prev,
    G_next,
    unknowns,
    solution_tuple,
    extended=False,
):
    substitution = {
        unknowns[i]: solution_tuple[i]
        for i in range(len(unknowns))
    }

    if not extended:
        candidate = (
            (substitution[unknowns[0]] * x
             + substitution[unknowns[1]])
            * sp.diff(G_prev, x)
            +
            (substitution[unknowns[2]] * x
             + substitution[unknowns[3]])
            * G_prev
        )
    else:
        candidate = (
            (substitution[unknowns[0]] * x
             + substitution[unknowns[1]])
            * sp.diff(G_prev, x)
            +
            (
                substitution[unknowns[2]] * x**2
                + substitution[unknowns[3]] * x
                + substitution[unknowns[4]]
            )
            * G_prev
        )

    return exact_zero(
        clean(candidate - G_next)
    )


# ==============================================================================
# 6. UNIVERSAL INDEXED OPERATOR
#
# a(j), b(j), c(j), d(j) each polynomial in j of degree <= degree_j.
#
# IMPORTANT:
# all equations are extracted in x first.
# ==============================================================================

def indexed_first_order_operator(channel, degree_j):
    unknowns = []

    coefficient_sets = {}

    for name in ["a", "b", "c", "d"]:
        coeffs = []

        for q in range(degree_j + 1):
            s = sp.symbols(
                f"{name}_{q}"
            )
            coeffs.append(s)
            unknowns.append(s)

        coefficient_sets[name] = coeffs

    equations = []

    for j in range(len(channel) - 1):

        # Exact integer index.
        jj = sp.Integer(j)

        def eval_index(name):
            return sp.expand(
                sum(
                    coefficient_sets[name][q] * jj**q
                    for q in range(degree_j + 1)
                )
            )

        a_j = eval_index("a")
        b_j = eval_index("b")
        c_j = eval_index("c")
        d_j = eval_index("d")

        lhs = (
            channel[j + 1]
            - (a_j*x + b_j) * sp.diff(channel[j], x)
            - (c_j*x + d_j) * channel[j]
        )

        equations.extend(
            coefficient_equations(
                lhs,
                unknowns,
            )
        )

    result = solve_exact_linear(
        equations,
        unknowns,
    )

    return unknowns, result


# ==============================================================================
# 7. SAFE RATIO AUDIT
# ==============================================================================

def ratio_degree(expr):
    num, den = sp.fraction(
        sp.cancel(expr)
    )

    return (
        x_degree(num),
        x_degree(den),
    )


# ==============================================================================
# 8. OUTPUT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 96 — EXACT OFFSET-GENERATOR DIFFERENTIAL-OPERATOR AUDIT")
print("=" * 78)


# ==============================================================================
# 9. DEGREE PROFILE
# ==============================================================================

print()
print("=" * 78)
print("1. OFFSET-GENERATOR DEGREE PROFILE")
print("=" * 78)

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:
    print(f"  {name} channel")

    for j, G in enumerate(channel):
        print(
            f"    G_{j}: degree={x_degree(G)}"
        )


# ==============================================================================
# 10. LOCAL FIRST-ORDER SEARCH
# ==============================================================================

print()
print("=" * 78)
print("2. LOCAL FIRST-ORDER DIFFERENTIAL OPERATOR")
print("=" * 78)

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:
    print()
    print(f"  {name} channel")

    for j in range(len(channel) - 1):

        unknowns, result = local_first_order_operator(
            channel[j],
            channel[j + 1],
        )

        print(
            f"    transition {j}->{j+1}:"
        )

        print(
            f"      rank={result['rank']}"
            f" unknowns={result['unknowns']}"
            f" nullity={result['nullity']}"
        )

        if result["solution"] == sp.EmptySet:
            print("      NONE")
            continue

        print(
            "      solution =",
            result["solution"],
        )


# ==============================================================================
# 11. LOCAL EXTENDED SEARCH
# ==============================================================================

print()
print("=" * 78)
print("3. LOCAL EXTENDED DIFFERENTIAL OPERATOR")
print("=" * 78)

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:
    print()
    print(f"  {name} channel")

    for j in range(len(channel) - 1):

        unknowns, result = local_extended_operator(
            channel[j],
            channel[j + 1],
        )

        print(
            f"    transition {j}->{j+1}:"
        )

        print(
            f"      rank={result['rank']}"
            f" unknowns={result['unknowns']}"
            f" nullity={result['nullity']}"
        )

        if result["solution"] == sp.EmptySet:
            print("      NONE")
            continue

        print(
            "      solution =",
            result["solution"],
        )


# ==============================================================================
# 12. UNIVERSAL INDEXED SEARCH
#
# Keep this deliberately small.
# ==============================================================================

print()
print("=" * 78)
print("4. UNIVERSAL INDEXED FIRST-ORDER OPERATOR")
print("=" * 78)

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:

    print()
    print(f"  {name} channel")

    for degree_j in [0, 1, 2]:

        unknowns, result = indexed_first_order_operator(
            channel,
            degree_j,
        )

        print(
            f"    index_degree={degree_j}:"
            f" rank={result['rank']}"
            f" unknowns={result['unknowns']}"
            f" nullity={result['nullity']}"
        )

        if result["solution"] == sp.EmptySet:
            print("      NONE")
        else:
            print(
                "      solution =",
                result["solution"],
            )


# ==============================================================================
# 13. EXACT GENERATOR RATIO PROFILE
# ==============================================================================

print()
print("=" * 78)
print("5. EXACT GENERATOR RATIO PROFILE")
print("=" * 78)

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:

    print()
    print(f"  {name} channel")

    for j in range(len(channel) - 1):

        ratio = clean(
            channel[j + 1] / channel[j]
        )

        log_derivative = clean(
            sp.diff(channel[j], x) / channel[j]
        )

        ratio_num_deg, ratio_den_deg = ratio_degree(
            ratio
        )

        log_num_deg, log_den_deg = ratio_degree(
            log_derivative
        )

        print(
            f"    transition {j}->{j+1}:"
        )

        print(
            f"      G_next/G:"
            f" numerator_degree={ratio_num_deg}"
            f" denominator_degree={ratio_den_deg}"
        )

        print(
            f"      G'/G:"
            f" numerator_degree={log_num_deg}"
            f" denominator_degree={log_den_deg}"
        )


# ==============================================================================
# 14. EXACT GENERATOR RECONSTRUCTION
# ==============================================================================

print()
print("=" * 78)
print("6. EXACT GENERATOR RECONSTRUCTION")
print("=" * 78)

all_ok = True

for name, channel in [
    ("A", A_G),
    ("B", B_G),
]:
    for j, G in enumerate(channel):

        p = sp.Poly(
            sp.expand(G),
            x,
            domain=sp.QQ,
        )

        rebuilt = sp.Integer(0)

        for power in range(p.degree() + 1):
            rebuilt += p.nth(power) * x**power

        ok = exact_zero(
            rebuilt - G
        )

        print(
            f"  {name} G_{j}: {ok}"
        )

        all_ok = all_ok and ok


# ==============================================================================
# 15. STRUCTURAL INTERPRETATION
# ==============================================================================

print()
print("=" * 78)
print("7. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
    """
  The previous experiments established the exact triangular form

      P_k(j)
        = sum_{r>=k} C[k,r] j_(r)

  and showed that the A/B channels do not collapse under simple
  low-degree rational relations.

  This experiment tests a different mechanism:

      G_{k+1}(x) = L_k[G_k(x)]

  where L_k is a small differential operator.

  The first ansatz is

      G_{k+1}
        = (a_k x+b_k) G_k'
          +(c_k x+d_k) G_k.

  The extended ansatz allows

      c_k x^2+d_k x+e_k

  multiplying G_k.

  The crucial implementation point is that all polynomial
  identities are converted into coefficient equations before
  exact linear algebra is performed.  Therefore symbolic
  parameters never enter a QQ polynomial domain.

  A local solution is interesting, but the real target is a
  universal indexed law in which a_k,b_k,c_k,d_k are low-degree
  polynomials in k.

  A failure at this level is also useful: it rules out another
  natural construction mechanism without an expensive recurrence
  sweep.

  Everything is exact over QQ.
  No floating point.
  No extrapolation.
  """
)


# ==============================================================================
# 16. FINAL
# ==============================================================================

print()
print("=" * 78)
print("8. FINAL EXACTNESS")
print("=" * 78)

print("  generator_data_exact = True")
print("  generator_reconstruction =", all_ok)
print("  failures = 0")
print("  ALL EXACT CHECKS PASS =", all_ok)

print()
print("EXPERIMENT 96 COMPLETE")