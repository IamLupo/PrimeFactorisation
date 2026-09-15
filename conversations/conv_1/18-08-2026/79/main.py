#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 81 — EXACT SCALAR CHANNEL RECURRENCE / INDEX-FACTOR SEARCH
==============================================================================

Goal
----
Experiment 80 showed that a generic indexed order-1 matrix law

    V_j = M(j,u) V_{j-1}

does not fit within small polynomial degree bounds.

Experiment 81 changes direction.

We study the scalar channel sequences

    A_j = A_{2j}(u)
    B_j = B_{2j+1}(u)

and search for exact order-2 recurrences

    c0(j,u) H_j
  + c1(j,u) H_{j-1}
  + c2(j,u) H_{j-2}
  = 0

but ONLY inside a deliberately small family of structured index forms.

The search considers coefficient templates built from:

    j
    j+1
    2j+1
    2j+3
    u
    u+j
    u+j+1
    u+2j
    u+2j+1

and products of these factors of bounded size.

The experiment also performs:

    * exact constant-coefficient recurrence tests
    * exact polynomial-in-u recurrence tests
    * exact structured-index recurrence tests
    * residual verification
    * recurrence normalization
    * finite-difference diagnostics
    * ratio diagnostics
    * cross-channel comparison
    * fresh numerical-semiprime evaluation

NO FLOATING POINT.
NO FITTING BY NUMERICAL REGRESSION.
NO INTERPOLATION IS ACCEPTED AS A DISCOVERED LAW.

Every candidate recurrence is symbolically verified at every available
transition.

The purpose is to determine whether the channel construction has a compact
hypergeometric / holonomic-style scalar recurrence.

==============================================================================
"""

import sys
import sympy as sp

# ---------------------------------------------------------------------------
# Exact symbols
# ---------------------------------------------------------------------------

t = sp.Symbol("t")
u = sp.Symbol("u")
j = sp.Symbol("j")
z = sp.Symbol("z")

QQ = sp.QQ


# ---------------------------------------------------------------------------
# Canonical polynomial helper
# ---------------------------------------------------------------------------

def S(x):
    return sp.sympify(x)


def poly_u(expr):
    return sp.Poly(sp.expand(S(expr)), u, domain=QQ)


def degree_u(expr):
    p = poly_u(expr)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def content_u(expr):
    p = poly_u(expr)
    if p.is_zero:
        return sp.Integer(0)

    vals = [sp.Rational(c) for c in p.all_coeffs()]
    den = sp.ilcm(*[c.q for c in vals]) if len(vals) >= 2 else vals[0].q
    ints = [int(c * den) for c in vals]

    g = 0
    for v in ints:
        g = sp.igcd(g, abs(v))

    if g == 0:
        return sp.Integer(0)

    return sp.Rational(g, den)


def primitive_coefficients(expr):
    p = poly_u(expr)
    if p.is_zero:
        return []

    coeffs = [sp.Rational(c) for c in p.all_coeffs()]
    den = sp.ilcm(*[c.q for c in coeffs]) if len(coeffs) >= 2 else coeffs[0].q

    ints = [int(c * den) for c in coeffs]

    g = 0
    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return ints

    ints = [x // g for x in ints]

    # Normalize sign so leading coefficient is positive.
    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break

    return ints


# ---------------------------------------------------------------------------
# h_d(t) from the exact Experiment 72 structure
# ---------------------------------------------------------------------------

H = {
    16: (
        -(3 * t**2 + 3 * t + 1)
        * (3 * t**6 + 9 * t**5 + 18 * t**4 + 21 * t**3
           + 15 * t**2 + 6 * t + 1)
    ),

    15: (
        (2 * t + 1)
        * (
            44 * t**8 + 176 * t**7 + 494 * t**6 + 866 * t**5
            + 1016 * t**4 + 794 * t**3 + 401 * t**2 + 119 * t + 16
        )
    ),

    14: (
        -276 * t**10 - 1380 * t**9 - 5460 * t**8 - 13560 * t**7
        - 23058 * t**6 - 27510 * t**5 - 23100 * t**4
        - 13410 * t**3 - 5135 * t**2 - 1169 * t - 120
    ),

    13: (
        (2 * t + 1)
        * (
            144 * t**10 + 720 * t**9 + 4810 * t**8 + 14920 * t**7
            + 32272 * t**6 + 47620 * t**5 + 48955 * t**4
            + 34510 * t**3 + 16020 * t**2 + 4431 * t + 560
        )
    ),

    12: (
        -50 * t**12 - 300 * t**11 - 7436 * t**10 - 34430 * t**9
        - 117810 * t**8 - 267960 * t**7 - 426888 * t**6
        - 484110 * t**5 - 390225 * t**4 - 219010 * t**3
        - 81510 * t**2 - 18109 * t - 1820
    ),

    11: (
        13 * (2 * t + 1)
        * (
            50 * t**10 + 250 * t**9 + 2284 * t**8 + 7636 * t**7
            + 17434 * t**6 + 26626 * t**5 + 28036 * t**4
            + 20104 * t**3 + 9451 * t**2 + 2639 * t + 336
        )
    ),

    10: (
        -1001
        * (
            10 * t**10 + 50 * t**9 + 252 * t**8 + 708 * t**7
            + 1302 * t**6 + 1638 * t**5 + 1428 * t**4
            + 852 * t**3 + 333 * t**2 + 77 * t + 8
        )
    ),

    9: (
        (2 * t + 1)
        * (
            17875 * t**8 + 71500 * t**7 + 245960 * t**6
            + 487630 * t**5 + 622414 * t**4 + 515528 * t**3
            + 271506 * t**2 + 83097 * t + 11441
        )
    ),

    8: (
        -71500 * t**8 - 286000 * t**7 - 755664 * t**6
        - 1265992 * t**5 - 1375360 * t**4 - 974400 * t**3
        - 436832 * t**2 - 112964 * t - 12879
    ),

    7: (
        4 * (2 * t + 1)
        * (
            11050 * t**6 + 33150 * t**5 + 63036 * t**4
            + 70822 * t**3 + 47549 * t**2 + 17663 * t + 2869
        )
    ),

    6: (
        -4
        * (
            17850 * t**6 + 53550 * t**5 + 86088 * t**4
            + 82926 * t**3 + 47574 * t**2 + 15036 * t + 2023
        )
    ),

    5: (
        6 * (2 * t + 1)
        * (3230 * t**4 + 6460 * t**3 + 6688 * t**2 + 3458 * t + 749)
    ),

    4: (
        -2
        * (7125 * t**4 + 14250 * t**3 + 12690 * t**2 + 5565 * t + 973)
    ),

    3: (
        14 * (2 * t + 1) * (125 * t**2 + 125 * t + 46)
    ),

    2: (
        -2 * (275 * t**2 + 275 * t + 78)
    ),

    1: 25 * (2 * t + 1),

    0: -2,
}


# ---------------------------------------------------------------------------
# Convert h_d(t) into A_d(u), B_d(u)
# ---------------------------------------------------------------------------

def involution_decompose(expr):
    """
    Exact decomposition

        h(t) = A(u) + (2t+1)B(u)
        u = t(t+1).

    Uses y = 2t+1, hence u=(y^2-1)/4.

    We first symmetrize:

        even = (h(t)+h(-1-t))/2
        odd  = (h(t)-h(-1-t))/2

    Then recover A and B from the resulting invariant pieces.
    """

    f = sp.expand(expr)
    partner = sp.expand(f.subs(t, -1 - t))

    even = sp.cancel((f + partner) / 2)
    odd = sp.cancel((f - partner) / 2)

    # A(u): replace t^2+t by u using polynomial division/rewrite.
    # Work by introducing y=2t+1.
    y = sp.Symbol("_y")

    even_y = sp.expand(even.subs(t, (y - 1) / 2))
    odd_y = sp.expand(odd.subs(t, (y - 1) / 2))

    # even should contain only even powers of y.
    A = 0
    for (pow_y,), coeff in sp.Poly(even_y, y, domain=QQ).terms():
        if pow_y % 2 != 0:
            raise ValueError("Even channel unexpectedly contains odd y-power.")
        A += coeff * (4 * u + 1) ** (pow_y // 2)

    # For odd channel divide by y first.
    odd_div = sp.cancel(odd_y / y)

    B = 0
    for (pow_y,), coeff in sp.Poly(odd_div, y, domain=QQ).terms():
        if pow_y % 2 != 0:
            raise ValueError("Odd channel quotient unexpectedly contains odd y-power.")
        B += coeff * (4 * u + 1) ** (pow_y // 2)

    # The substitution y^2=4u+1 above produces the correct expression,
    # but powers are expressed in expanded form.
    A = sp.expand(A)
    B = sp.expand(B)

    # Verify exactly.
    reconstructed = sp.expand(
        A.subs(u, t * (t + 1))
        + (2 * t + 1) * B.subs(u, t * (t + 1))
    )

    if sp.simplify(reconstructed - f) != 0:
        raise AssertionError("Involution decomposition failed.")

    return A, B


# ---------------------------------------------------------------------------
# Build channel sequences
# ---------------------------------------------------------------------------

A = {}
B = {}

for d, expr in H.items():
    A[d], B[d] = involution_decompose(expr)


A_sequence = [A[d] for d in range(0, 17, 2)]
B_sequence = [B[d] for d in range(1, 17, 2)]

# A_sequence:
# A_0, A_2, ..., A_16
#
# B_sequence:
# B_1, B_3, ..., B_15


# ---------------------------------------------------------------------------
# Exact polynomial over QQ(u) utilities
# ---------------------------------------------------------------------------

def rat_simplify(x):
    return sp.cancel(sp.factor(sp.sympify(x)))


def exact_zero(x):
    return sp.simplify(sp.cancel(x)) == 0


def coeff(expr, k):
    return sp.expand(expr).coeff(u, k)


# ---------------------------------------------------------------------------
# Test constant-in-j order-2 recurrence
#
# H_j = c1(u)H_{j-1}+c2(u)H_{j-2}
#
# This is deliberately tested first because it is cheap.
# ---------------------------------------------------------------------------

def solve_constant_order2(seq):
    """
    Search for

        H_j = c1(u) H_{j-1} + c2(u) H_{j-2}

    where c1,c2 are rational functions in u.

    Since the sequence is polynomial, we solve the first usable
    equations over QQ(u).

    Returns:
        (found, c1, c2)
    """

    if len(seq) < 3:
        return False, None, None

    # Use j=2 first:
    # H2 = c1 H1 + c2 H0
    #
    # One polynomial identity gives many equations in u, so this determines
    # only a one-dimensional family. Add the remaining transitions.

    c1, c2 = sp.symbols("_c1 _c2")

    equations = []

    # Rather than creating symbolic unknown polynomials,
    # solve pointwise in the rational function field QQ(u).
    # Use the first two transitions as a 2x2 linear system:
    #
    # [H1 H0] [c1] = H2
    # [H2 H1] [c2]   H3
    #
    # interpreted over QQ(u).

    if len(seq) < 4:
        return False, None, None

    M = sp.Matrix([
        [seq[1], seq[0]],
        [seq[2], seq[1]],
    ])

    rhs = sp.Matrix([seq[2], seq[3]])

    det = sp.simplify(M.det())

    if exact_zero(det):
        return False, None, None

    sol = M.inv() * rhs
    c1_val = rat_simplify(sol[0])
    c2_val = rat_simplify(sol[1])

    for idx in range(2, len(seq)):
        residual = rat_simplify(
            seq[idx]
            - c1_val * seq[idx - 1]
            - c2_val * seq[idx - 2]
        )
        if not exact_zero(residual):
            return False, None, None

    return True, c1_val, c2_val


# ---------------------------------------------------------------------------
# Structured index factors
# ---------------------------------------------------------------------------

def factor_library(index):
    return {
        "1": sp.Integer(1),
        "j": index,
        "j+1": index + 1,
        "2j+1": 2 * index + 1,
        "2j+3": 2 * index + 3,
        "u": u,
        "u+j": u + index,
        "u+j+1": u + index + 1,
        "u+2j": u + 2 * index,
        "u+2j+1": u + 2 * index + 1,
    }


# ---------------------------------------------------------------------------
# Structured first-order recurrence
#
#    a(j,u) H_j = b(j,u) H_{j-1}
#
# with a,b selected from products of at most two primitive factors.
# ---------------------------------------------------------------------------

def generate_structured_coefficients(index, max_factor_count=2):
    lib = factor_library(index)
    names = list(lib.keys())

    expressions = []

    # Single factors
    for name in names:
        expressions.append((name, sp.expand(lib[name])))

    # Products of two factors
    if max_factor_count >= 2:
        for i, name1 in enumerate(names):
            for name2 in names[i:]:
                expressions.append(
                    (
                        f"({name1})*({name2})",
                        sp.expand(lib[name1] * lib[name2]),
                    )
                )

    # Remove exact duplicates.
    unique = []
    seen = set()

    for name, expr in expressions:
        key = sp.srepr(sp.expand(expr))
        if key not in seen:
            seen.add(key)
            unique.append((name, expr))

    return unique


def search_structured_order1(seq, seq_name):
    """
    Search exact relations of the form

        a(j,u) H_j = c b(j,u) H_{j-1}

    using the same structured factor family for all transitions.

    The scalar multiplier c is allowed to be rational.

    This is a deliberately small search, not a generic polynomial fit.
    """

    results = []

    # Sequence index:
    # seq[0] = H_0
    # seq[1] = H_1
    # ...
    #
    # transition r means H_r / H_{r-1}

    coeff_pairs_by_r = {}

    for r in range(1, len(seq)):
        coeff_pairs_by_r[r] = generate_structured_coefficients(sp.Integer(r))

    # Compare ratio structures.
    for r0 in range(1, len(seq)):
        for a_name, a_expr in coeff_pairs_by_r[r0]:
            lhs = rat_simplify(a_expr)

            for b_name, b_expr in coeff_pairs_by_r[r0]:
                rhs = rat_simplify(b_expr)

                ratio = rat_simplify(lhs / rhs)

                # Test whether H_r/H_{r-1} / ratio is constant across
                # every transition.
                base = rat_simplify(
                    seq[r0] / seq[r0 - 1] / ratio
                )

                if any(
                    exact_zero(seq[k - 1]) for k in range(1, len(seq))
                ):
                    continue

                ok = True

                for r in range(1, len(seq)):
                    coeffs = factor_library(sp.Integer(r))

                    a_r = rat_simplify(a_expr.xreplace({
                        sp.Symbol("_dummy"): 0
                    }))

                    # Rebuild structurally for the actual index.
                    # Find expressions by string names.
                    actual_a = dict(
                        generate_structured_coefficients(sp.Integer(r))
                    ).get(a_name)

                    actual_b = dict(
                        generate_structured_coefficients(sp.Integer(r))
                    ).get(b_name)

                    if actual_a is None or actual_b is None:
                        ok = False
                        break

                    candidate = rat_simplify(
                        seq[r] / seq[r - 1]
                        * actual_b / actual_a
                    )

                    if not exact_zero(candidate - base):
                        ok = False
                        break

                if ok:
                    results.append((a_name, b_name, base))

    return results


# ---------------------------------------------------------------------------
# Structured order-2 recurrence
#
#    a(j,u) H_j + b(j,u) H_{j-1} + c(j,u) H_{j-2} = 0
#
# We search only coefficient shapes from a small factor library.
# ---------------------------------------------------------------------------

def search_structured_order2(seq, seq_name, max_factor_count=2):
    """
    Exact search for a compact order-2 relation.

    We choose a coefficient template a(j,u), b(j,u), c(j,u) from the
    structured library and test whether the three-term combination vanishes
    identically for every transition.

    To keep the search finite, the coefficient expressions are normalized
    by fixing the leading coefficient template to 1.
    """

    if len(seq) < 3:
        return []

    results = []

    # Use a relatively small family. This is intentional.
    # A large unrestricted search becomes interpolation rather than
    # structural discovery.
    base_names = [
        "1",
        "j",
        "j+1",
        "2j+1",
        "2j+3",
        "u",
        "u+j",
        "u+j+1",
        "u+2j",
        "u+2j+1",
    ]

    candidates = []

    for r in range(1, len(seq)):
        items = dict(
            generate_structured_coefficients(
                sp.Integer(r),
                max_factor_count=max_factor_count,
            )
        )
        candidates.append(items)

    # To control runtime, use only single factors plus selected products.
    candidate_names = base_names + [
        "u+j+1",
        "u+2j+1",
    ]

    for a_name in candidate_names:
        for b_name in candidate_names:
            for c_name in candidate_names:

                # At each transition reconstruct the actual coefficient.
                ok = True
                scale = None

                for r in range(2, len(seq)):
                    items = candidates[r - 1]

                    if a_name not in items or b_name not in items or c_name not in items:
                        ok = False
                        break

                    a_r = items[a_name]
                    b_r = items[b_name]
                    c_r = items[c_name]

                    residual = rat_simplify(
                        a_r * seq[r]
                        + b_r * seq[r - 1]
                        + c_r * seq[r - 2]
                    )

                    if exact_zero(residual):
                        continue

                    # No constant scale can repair an identity that already
                    # fails for the exact polynomial at one transition.
                    ok = False
                    break

                if ok:
                    results.append(
                        (a_name, b_name, c_name)
                    )

    return results


# ---------------------------------------------------------------------------
# Exact ratio diagnostics
# ---------------------------------------------------------------------------

def ratio_profile(seq):
    out = []
    for r in range(1, len(seq)):
        out.append(
            (
                r,
                rat_simplify(seq[r] / seq[r - 1]),
            )
        )
    return out


# ---------------------------------------------------------------------------
# Finite-difference diagnostics with respect to channel index
# ---------------------------------------------------------------------------

def finite_difference(values):
    result = list(values)
    levels = [result]

    while len(result) > 1:
        result = [
            rat_simplify(result[k + 1] - result[k])
            for k in range(len(result) - 1)
        ]
        levels.append(result)

    return levels


def polynomial_index_degree(values):
    """
    Since the values are polynomials in u, determine the degree in j of each
    coefficient sequence using finite differences.

    For each coefficient u^k we only have a finite exact sample. We report
    whether the finite difference becomes identically zero before the final
    step.
    """

    max_u = 0
    for expr in values:
        max_u = max(max_u, degree_u(expr))

    profile = {}

    for k in range(max_u + 1):
        coeff_values = [sp.expand(expr).coeff(u, k) for expr in values]
        diffs = finite_difference(coeff_values)

        degree = None
        for d, level in enumerate(diffs):
            if all(x == 0 for x in level):
                degree = d
                break

        if degree is None:
            degree = len(values) - 1

        profile[k] = degree

    return profile


# ---------------------------------------------------------------------------
# Fresh semiprime check
# ---------------------------------------------------------------------------

def is_prime(n):
    return bool(sp.isprime(int(n)))


def fresh_check():
    p = sp.Integer(106621)
    q = sp.Integer(246473)

    N = p * q
    S0 = p + q
    X = S0 + 1
    t0 = sp.Rational(N, X)
    u0 = sp.factor(t0 * (t0 + 1))

    print("===============================================================================")
    print("FRESH EXACT SEMIPRIME")
    print("===============================================================================")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {N}")
    print(f"  S = {S0}")
    print(f"  X = {X}")
    print(f"  t = {t0}")
    print(f"  u = {u0}")
    print()

    checks = []

    for d in range(17):
        f = sp.expand(H[d])
        a = A[d]
        b = B[d]

        expected = sp.expand(
            a.subs(u, u0)
            + (2 * t0 + 1) * b.subs(u, u0)
        )
        observed = sp.expand(f.subs(t, t0))

        ok = sp.simplify(observed - expected) == 0
        checks.append(ok)

        print(f"  d={d:2d}: channel evaluation={ok}")

    print()
    print(f"  all fresh evaluations = {all(checks)}")
    print()

    return all(checks)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    print("==============================================================================")
    print("EXPERIMENT 81 — EXACT SCALAR CHANNEL RECURRENCE / INDEX-FACTOR SEARCH")
    print("==============================================================================")
    print()
    print("0. EXACT SYMBOLIC SETUP")
    print("  t = N/X")
    print("  u = t(t+1)")
    print("  A_j = A_{2j}(u)")
    print("  B_j = B_{2j+1}(u)")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")
    print()

    # -----------------------------------------------------------------------
    # Involution validation
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("1. EXACT INVOLUTION VALIDATION")
    print("==============================================================================")

    involution_ok = True

    for d in range(16, -1, -1):

        f = sp.expand(H[d])
        partner = sp.expand(f.subs(t, -1 - t))

        even = sp.expand((f + partner) / 2)
        odd = sp.expand((f - partner) / 2)

        a = A[d]
        b = B[d]

        reconstructed = sp.expand(
            a.subs(u, t * (t + 1))
            + (2 * t + 1) * b.subs(u, t * (t + 1))
        )

        ok = sp.simplify(reconstructed - f) == 0

        invariant = sp.simplify(partner - f) == 0
        anti = sp.simplify(partner + f) == 0
        divisible = sp.rem(
            sp.Poly(odd, t, domain=QQ),
            sp.Poly(2 * t + 1, t, domain=QQ),
        ).is_zero

        involution_ok = involution_ok and ok

        print(
            f"  d={d:2d}: "
            f"invariant={invariant} "
            f"anti_invariant={anti} "
            f"divisible={divisible} "
            f"reconstruction={ok}"
        )

    print()
    print(f"  ALL INVOLUTION CHECKS = {involution_ok}")
    print()

    # -----------------------------------------------------------------------
    # Print channel profiles
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("2. SCALAR CHANNEL SEQUENCES")
    print("==============================================================================")

    print("  A-sequence:")
    for idx, expr in enumerate(A_sequence):
        d = 2 * idx
        print(
            f"    j={idx}: A_{d} "
            f"degree_u={degree_u(expr)} "
            f"primitive={primitive_coefficients(expr)}"
        )

    print()
    print("  B-sequence:")
    for idx, expr in enumerate(B_sequence):
        d = 2 * idx + 1
        print(
            f"    j={idx}: B_{d} "
            f"degree_u={degree_u(expr)} "
            f"primitive={primitive_coefficients(expr)}"
        )

    print()

    # -----------------------------------------------------------------------
    # Constant-in-j order 2
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("3. EXACT CONSTANT-IN-j ORDER-2 RECURRENCE")
    print("==============================================================================")

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:
        found, c1, c2 = solve_constant_order2(seq)

        print(f"  {name}-channel:")
        print(f"    recurrence found = {found}")

        if found:
            print(f"    c1(u) = {c1}")
            print(f"    c2(u) = {c2}")

        print()

    # -----------------------------------------------------------------------
    # Ratio profile
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("4. EXACT ADJACENT RATIO PROFILE")
    print("==============================================================================")

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:
        print(f"  {name}-channel")

        ratios = ratio_profile(seq)

        for r, ratio in ratios:
            print(f"    H_{r}/H_{r-1} = {ratio}")

        print()

    # -----------------------------------------------------------------------
    # Finite-difference degree diagnostics
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("5. INDEX FINITE-DIFFERENCE DIAGNOSTIC")
    print("==============================================================================")

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:
        profile = polynomial_index_degree(seq)

        print(f"  {name}-channel coefficient sequences")

        for k, deg in profile.items():
            print(f"    [u^{k}] : finite-difference degree = {deg}")

        print()

    # -----------------------------------------------------------------------
    # Structured order-1 search
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("6. STRUCTURED ORDER-1 INDEX-FACTOR SEARCH")
    print("==============================================================================")

    structured_results = {}

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:

        print(f"  {name}-channel")

        found = search_structured_order1(seq, name)
        structured_results[name] = found

        if not found:
            print("    NONE")
        else:
            for a_name, b_name, scale in found:
                print(
                    f"    EXISTS: "
                    f"{a_name}*H_j = "
                    f"({scale})*{b_name}*H_(j-1)"
                )

        print()

    # -----------------------------------------------------------------------
    # Structured order-2 search
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("7. STRUCTURED ORDER-2 INDEX-FACTOR SEARCH")
    print("==============================================================================")

    order2_results = {}

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:

        print(f"  {name}-channel")

        results = search_structured_order2(
            seq,
            name,
            max_factor_count=1,
        )

        order2_results[name] = results

        if not results:
            print("    NONE")
        else:
            for a_name, b_name, c_name in results:
                print(
                    f"    EXISTS: "
                    f"({a_name}) H_j "
                    f"+ ({b_name}) H_(j-1) "
                    f"+ ({c_name}) H_(j-2) = 0"
                )

        print()

    # -----------------------------------------------------------------------
    # Explicit small-template search using the most plausible factors.
    # This is the most important structural search in this experiment.
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("8. TARGETED INDEX-FACTOR TEMPLATES")
    print("==============================================================================")

    templates = {
        "j": j,
        "j+1": j + 1,
        "2j+1": 2 * j + 1,
        "2j+3": 2 * j + 3,
        "u+j": u + j,
        "u+j+1": u + j + 1,
        "u+2j": u + 2 * j,
        "u+2j+1": u + 2 * j + 1,
    }

    # Search relations of the form
    #
    # a(j,u) H_j + b(j,u) H_(j-1) + c(j,u) H_(j-2) = 0
    #
    # where each coefficient is one of the basic templates above.
    #
    # Signs are independently considered.

    targeted_hits = []

    for name, seq in [
        ("A", A_sequence),
        ("B", B_sequence),
    ]:

        for an, aexpr in templates.items():
            for bn, bexpr in templates.items():
                for cn, cexpr in templates.items():

                    for sa in [1, -1]:
                        for sb in [1, -1]:
                            for sc in [1, -1]:

                                ok = True

                                for r in range(2, len(seq)):
                                    ar = sp.expand(
                                        sa * aexpr.subs(j, r)
                                    )
                                    br = sp.expand(
                                        sb * bexpr.subs(j, r)
                                    )
                                    cr = sp.expand(
                                        sc * cexpr.subs(j, r)
                                    )

                                    residual = sp.cancel(
                                        ar * seq[r]
                                        + br * seq[r - 1]
                                        + cr * seq[r - 2]
                                    )

                                    if residual != 0:
                                        ok = False
                                        break

                                if ok:
                                    targeted_hits.append(
                                        (
                                            name,
                                            sa,
                                            an,
                                            sb,
                                            bn,
                                            sc,
                                            cn,
                                        )
                                    )

    if targeted_hits:
        for hit in targeted_hits:
            name, sa, an, sb, bn, sc, cn = hit
            print(
                f"  FOUND {name}: "
                f"{sa}*({an})H_j + "
                f"{sb}*({bn})H_(j-1) + "
                f"{sc}*({cn})H_(j-2) = 0"
            )
    else:
        print("  NO targeted order-2 recurrence found.")

    print()

    # -----------------------------------------------------------------------
    # Cross-channel relation diagnostics
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("9. CROSS-CHANNEL RATIO DIAGNOSTIC")
    print("==============================================================================")

    max_idx = min(len(A_sequence), len(B_sequence))

    for r in range(max_idx):
        a_expr = A_sequence[r]
        b_expr = B_sequence[r]

        print(
            f"  j={r}: "
            f"deg(A_j)={degree_u(a_expr)} "
            f"deg(B_j)={degree_u(b_expr)}"
        )

        if not exact_zero(a_expr) and not exact_zero(b_expr):
            print(f"       A_j/B_j = {rat_simplify(a_expr / b_expr)}")

    print()

    # -----------------------------------------------------------------------
    # Fresh semiprime validation
    # -----------------------------------------------------------------------

    fresh_ok = fresh_check()

    # -----------------------------------------------------------------------
    # Final summary
    # -----------------------------------------------------------------------

    print("==============================================================================")
    print("10. FINAL DIAGNOSTIC")
    print("==============================================================================")

    print(
        """
  Experiment 80 showed that a generic small-degree indexed matrix
  recurrence was already inconsistent.

  Experiment 81 therefore restricts the search to structured scalar
  recurrences motivated by the actual algebraic variables appearing
  in the layer system.

  The hierarchy is:

      constant in j
          ↓
      simple index factors
          ↓
      simple order-2 index factors
          ↓
      only then consider a larger rational-function ansatz.

  A recurrence is accepted only when the symbolic residual is exactly
  zero for every available layer transition.

  Importantly, a recurrence over the finite observed range is not by
  itself evidence for a universal law.  The next stage, when a candidate
  is found, must test it on independently generated kernels / parameter
  values.

  The desired structural outcome remains:

      kernel
        -> channel seed
        -> explicit index operator
        -> all layers

  The factorization question remains separate:

      N
        -> computable layer information
        -> t
        -> X
        -> p,q
  """
    )

    print("==============================================================================")
    print("FINAL EXACTNESS")
    print("==============================================================================")

    all_basic = involution_ok and fresh_ok

    print(f"  involution = {involution_ok}")
    print(f"  fresh_validation = {fresh_ok}")
    print(f"  targeted_order2_hits = {len(targeted_hits)}")
    print(f"  ALL BASIC CHECKS PASS = {all_basic}")

    print()
    print("==============================================================================")
    print("EXPERIMENT 81 COMPLETE")
    print("==============================================================================")


if __name__ == "__main__":
    main()

