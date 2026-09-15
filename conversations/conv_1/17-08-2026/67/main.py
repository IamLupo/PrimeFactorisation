#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 132R
WEIGHTED HOMOGENEOUS NEWTON KERNEL
DIRECT TOP-EDGE DERIVATION

TARGET:
    j = ell-k-d-1

WEIGHTING:
    deg(S) = 1
    deg(N) = 2

SUBSTITUTION:
    S = t*x
    N = t^2

NEWTON LIMIT:
    Pi_0 = 2
    Pi_1 = x
    Pi_j = x Pi_{j-1} - Pi_{j-2}

GOAL
----
Determine whether the leading N-coefficient of

    c_{ell-k-d-1}(N)

comes entirely from the weighted homogeneous top of the exact quotient.

NO FULL NEWTON TENSOR
NO C/D TENSOR
NO POLYNOMIAL FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")
t, x = sp.symbols("t x")


# ============================================================================
# Exact finite kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient
# ============================================================================

def quotient_Q(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(PF, PD)

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric reduction to S,N
# ============================================================================

def symmetric_to_NS(expr):
    expr_sub = sp.expand(
        expr.subs(q, S - p)
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        expr_sub,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    p_coeff = sp.expand(
        rem.coeff(p, 1)
    )

    if p_coeff != 0:
        raise ArithmeticError(
            "Symmetric reduction retained p-dependence: "
            f"{sp.factor(p_coeff)}"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Safe exponent helper
# ============================================================================

def safe_integer_exponent(value):
    """
    Normalize SymPy/native integer exponents safely.

    The previous version assumed every exponent had an .is_Integer
    attribute. Native Python ints do not.
    """

    value = sp.sympify(value)

    if not value.is_Integer:
        raise ArithmeticError(
            f"Non-integer exponent encountered: {value}"
        )

    return int(value)


# ============================================================================
# Weighted monomial degree
# ============================================================================

def weighted_monomial_degree(term):
    """
    Weight:
        deg(S) = 1
        deg(N) = 2

    Constants therefore have degree 0.
    """

    term = sp.sympify(term)

    powers = term.as_powers_dict()

    s_pow = safe_integer_exponent(
        powers.get(S, 0)
    )

    n_pow = safe_integer_exponent(
        powers.get(N, 0)
    )

    return s_pow + 2 * n_pow


# ============================================================================
# Weighted top
# ============================================================================

def weighted_top(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return 0, sp.Integer(0)

    terms = sp.Add.make_args(expr)

    degree_pairs = []

    for term in terms:
        degree_pairs.append(
            (
                weighted_monomial_degree(term),
                term,
            )
        )

    top_degree = max(
        degree
        for degree, _ in degree_pairs
    )

    top_expr = sp.expand(
        sum(
            term
            for degree, term in degree_pairs
            if degree == top_degree
        )
    )

    return top_degree, top_expr


# ============================================================================
# Weighted substitution
# ============================================================================

def weighted_image(expr):
    return sp.expand(
        expr.subs(
            {
                S: t*x,
                N: t**2,
            }
        )
    )


def leading_t_form(expr):
    scaled = weighted_image(expr)

    poly = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x),
    )

    if poly.is_zero:
        return None, sp.Integer(0)

    degree = int(
        poly.degree()
    )

    leading = sp.expand(
        poly.LC()
    )

    return degree, leading


# ============================================================================
# Homogeneous Newton basis
# ============================================================================

def build_pi(max_j):
    """
    Pi_0 = 2
    Pi_1 = x
    Pi_j = x Pi_{j-1} - Pi_{j-2}
    """

    Pi = [sp.Integer(2)]

    if max_j >= 1:
        Pi.append(x)

    for j in range(
        2,
        max_j + 1
    ):
        Pi.append(
            sp.expand(
                x * Pi[j - 1]
                - Pi[j - 2]
            )
        )

    return Pi


# ============================================================================
# Homogeneous Newton decomposition
# ============================================================================

def homogeneous_newton_coefficients(H, max_j):
    """
    Decompose:

        H(x) = sum_j a_j Pi_j(x).

    The Pi_j are monic in x for j>=1.
    """

    Pi = build_pi(
        max_j
    )

    work = sp.Poly(
        sp.expand(H),
        x,
        domain=sp.QQ,
    )

    coeffs = {}

    for j in range(
        max_j,
        0,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(
                x**j
            )
        )

        coeffs[j] = c

        if c != 0:
            work -= sp.Poly(
                sp.expand(
                    c * Pi[j]
                ),
                x,
                domain=sp.QQ,
            )

    residual = sp.expand(
        work.as_expr()
    )

    if residual != 0:

        if residual.is_Rational:
            coeffs[0] = sp.simplify(
                residual / 2
            )
            residual = sp.Integer(0)

        else:
            raise ArithmeticError(
                "Nonconstant residual after homogeneous "
                "Newton decomposition: "
                f"{sp.factor(residual)}"
            )

    return coeffs


# ============================================================================
# Established formulas
# ============================================================================

def expected_degree(k, d):
    return k + (d - 1) // 2


def expected_lead(k, ell, d):
    if d % 2 == 1:
        return sp.expand(
            (-1) ** ((d - 1) // 2)
            * (k + ell)
        )

    return sp.expand(
        (-1) ** (d // 2)
        * sp.Rational(k + ell, 2)
        * (k - ell + d + 1)
    )


# ============================================================================
# Exact Newton basis
# ============================================================================

def build_exact_newton_basis(max_j):
    P = [sp.Integer(2)]

    if max_j >= 1:
        P.append(S)

    for j in range(
        2,
        max_j + 1
    ):
        P.append(
            sp.expand(
                S * P[j - 1]
                - N * P[j - 2]
            )
        )

    return P


# ============================================================================
# Exact target extraction
# ============================================================================

def extract_exact_target(
    Q_NS,
    ell,
    target_j,
):
    P = build_exact_newton_basis(
        ell - 1
    )

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    for j in range(
        ell - 1,
        target_j - 1,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(
                S**j
            )
        )

        if j == target_j:
            return c

        if c != 0:
            work -= sp.Poly(
                sp.expand(
                    c * P[j]
                ),
                S,
                domain=sp.QQ.frac_field(N),
            )

    raise ArithmeticError(
        f"Target j={target_j} not found."
    )


# ============================================================================
# Analyze one case
# ============================================================================

def analyze_case(k, ell, d):

    j = ell - k - d - 1

    if j < 1:
        return None

    Qpq = quotient_Q(
        k,
        ell
    )

    QNS = symmetric_to_NS(
        Qpq
    )

    # Independent exact coefficient.
    c = extract_exact_target(
        QNS,
        ell,
        j
    )

    poly = sp.Poly(
        sp.expand(c),
        N,
        domain=sp.QQ,
    )

    exact_degree = int(
        poly.degree()
    )

    exact_lead = sp.expand(
        poly.LC()
    )

    expected_D = expected_degree(
        k,
        d
    )

    expected_L = expected_lead(
        k,
        ell,
        d
    )

    # Weighted top of full quotient.
    W, H = weighted_top(
        QNS
    )

    t_degree, Hx = leading_t_form(
        QNS
    )

    H_check = sp.expand(
        H.subs(
            {
                S: x,
                N: 1,
            }
        )
    )

    if sp.expand(
        Hx - H_check
    ) != 0:
        raise ArithmeticError(
            "Weighted top substitution mismatch."
        )

    coeffs = homogeneous_newton_coefficients(
        Hx,
        ell - 1
    )

    homogeneous_target = sp.expand(
        coeffs.get(
            j,
            0
        )
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "c": c,
        "exact_degree": exact_degree,
        "exact_lead": exact_lead,
        "weighted_degree": W,
        "weighted_t_degree": t_degree,
        "H": H,
        "H_x": Hx,
        "homogeneous_target": homogeneous_target,
        "expected_D": expected_D,
        "expected_L": expected_L,
        "degree_ok": (
            exact_degree == expected_D
        ),
        "direct_lead_ok": (
            sp.expand(
                exact_lead
                - expected_L
            ) == 0
        ),
        "homogeneous_ok": (
            sp.expand(
                homogeneous_target
                - expected_L
            ) == 0
        ),
    }


# ============================================================================
# Grid
# ============================================================================

def build_grid(max_ell):
    rows = []

    for ell in range(
        3,
        max_ell + 1,
        2
    ):

        for k in range(
            1,
            ell,
            2
        ):

            for d in range(
                1,
                ell - k - 1
            ):

                row = analyze_case(
                    k,
                    ell,
                    d
                )

                if row is not None:
                    rows.append(
                        row
                    )

    return rows


# ============================================================================
# Validation
# ============================================================================

def validate(rows):

    degree_fail = [
        row for row in rows
        if not row["degree_ok"]
    ]

    direct_fail = [
        row for row in rows
        if not row["direct_lead_ok"]
    ]

    homogeneous_fail = [
        row for row in rows
        if not row["homogeneous_ok"]
    ]

    print("=" * 78)
    print("1. WEIGHTED TOP-KERNEL CERTIFICATE")
    print("=" * 78)

    print(
        f"rows = {len(rows)}"
    )

    print(
        f"degree failures = "
        f"{len(degree_fail)}"
    )

    print(
        f"direct leading failures = "
        f"{len(direct_fail)}"
    )

    print(
        f"homogeneous leading failures = "
        f"{len(homogeneous_fail)}"
    )

    if homogeneous_fail:

        print()
        print(
            "FIRST HOMOGENEOUS FAILURES"
        )
        print("-" * 78)

        for row in homogeneous_fail[:12]:

            print(
                f"({row['k']},{row['ell']}) "
                f"d={row['d']} "
                f"j={row['j']}"
            )

            print(
                f"  H(x) = "
                f"{sp.factor(row['H_x'])}"
            )

            print(
                f"  homogeneous coefficient = "
                f"{row['homogeneous_target']}"
            )

            print(
                f"  expected = "
                f"{row['expected_L']}"
            )

            print()

    return (
        degree_fail,
        direct_fail,
        homogeneous_fail
    )


# ============================================================================
# Examples
# ============================================================================

def print_examples(rows):

    wanted = {
        (1, 5, 1),
        (3, 7, 1),
        (1, 7, 2),
        (5, 9, 1),
        (1, 9, 3),
        (3, 11, 4),
        (1, 13, 6),
        (5, 13, 5),
        (3, 17, 8),
    }

    print("=" * 78)
    print(
        "2. HOMOGENEOUS TOP-KERNEL EXAMPLES"
    )
    print("=" * 78)

    for row in rows:

        key = (
            row["k"],
            row["ell"],
            row["d"]
        )

        if key not in wanted:
            continue

        print(
            f"({row['k']},{row['ell']}) "
            f"d={row['d']} "
            f"j={row['j']}"
        )

        print(
            f"  weighted degree = "
            f"{row['weighted_degree']}"
        )

        print(
            f"  H(x) = "
            f"{sp.factor(row['H_x'])}"
        )

        print(
            f"  Pi_j coefficient = "
            f"{sp.factor(row['homogeneous_target'])}"
        )

        print(
            f"  expected LC = "
            f"{row['expected_L']}"
        )

        print(
            f"  exact c_j(N) = "
            f"{sp.factor(row['c'])}"
        )

        print()


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ells):

    rows = []

    for ell in ells:

        for k in range(
            1,
            ell,
            2
        ):

            for d in range(
                1,
                ell - k - 1
            ):

                row = analyze_case(
                    k,
                    ell,
                    d
                )

                if row is not None:
                    rows.append(
                        row
                    )

    degree_fail = sum(
        not row["degree_ok"]
        for row in rows
    )

    direct_fail = sum(
        not row["direct_lead_ok"]
        for row in rows
    )

    homogeneous_fail = sum(
        not row["homogeneous_ok"]
        for row in rows
    )

    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    print(
        f"forward rows = {len(rows)}"
    )

    print(
        f"degree failures = "
        f"{degree_fail}"
    )

    print(
        f"direct leading failures = "
        f"{direct_fail}"
    )

    print(
        f"homogeneous leading failures = "
        f"{homogeneous_fail}"
    )

    return rows


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 132R")
    print("WEIGHTED HOMOGENEOUS NEWTON KERNEL")
    print("DIRECT TOP-EDGE DERIVATION")
    print("TARGET j = ell-k-d-1")
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 13
    FORWARD_ELL = [15]

    # ------------------------------------------------------------------
    # Sanity
    # ------------------------------------------------------------------

    print("=" * 78)
    print("0. QUOTIENT SANITY")
    print("=" * 78)

    sanity = [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]

    failures = 0

    for k, ell in sanity:

        try:
            quotient_Q(
                k,
                ell
            )

            print(
                f"({k},{ell}) quotient = PASS"
            )

        except Exception as exc:

            failures += 1

            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = "
        f"{failures}/{len(sanity)}"
    )

    if failures:
        raise ArithmeticError(
            "Quotient sanity failed."
        )

    print()

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    print("=" * 78)
    print("BUILDING WEIGHTED TOP-KERNEL DATASET")
    print("=" * 78)

    rows = build_grid(
        TRAIN_MAX_ELL
    )

    print(
        f"training rows = {len(rows)}"
    )

    print()

    degree_fail, direct_fail, homogeneous_fail = validate(
        rows
    )

    print()

    print_examples(
        rows
    )

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    forward = forward_holdout(
        FORWARD_ELL
    )

    forward_degree_fail = sum(
        not row["degree_ok"]
        for row in forward
    )

    forward_direct_fail = sum(
        not row["direct_lead_ok"]
        for row in forward
    )

    forward_homogeneous_fail = sum(
        not row["homogeneous_ok"]
        for row in forward
    )

    # ------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    total_failures = (
        len(degree_fail)
        + len(direct_fail)
        + len(homogeneous_fail)
        + forward_degree_fail
        + forward_direct_fail
        + forward_homogeneous_fail
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The highest-N coefficient of"
        )
        print(
            "    c_{ell-k-d-1}(N)"
        )
        print(
            "is recovered exactly from the weighted homogeneous"
        )
        print(
            "top of the full quotient."
        )
        print()
        print(
            "The homogeneous Newton recurrence"
        )
        print(
            "    Pi_0=2"
        )
        print(
            "    Pi_1=x"
        )
        print(
            "    Pi_j=x Pi_{j-1}-Pi_{j-2}"
        )
        print(
            "produces the same leading coefficient."
        )
        print()
        print(
            "This is the computational bridge from the exact"
        )
        print(
            "finite kernel to the observed post-boundary edge law."
        )

    else:

        print("STATUS = FAIL")

        print(
            f"total failures = "
            f"{total_failures}"
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print()

        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)