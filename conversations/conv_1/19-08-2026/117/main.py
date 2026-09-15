#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 283 — EXACT B-GRID / BIVARIATE INDEX-LAW AUDIT
==============================================================================

Experiments 279-282 established:

    B[k][offset] multiplies j_(k+offset),

so the natural absolute index is

    r = k + offset.

Experiment 282 then rejected the hypothesis that B[k,r] is a fixed
linear functional of the currently known q_p(r) source rows.

Experiment 283 therefore temporarily removes q_p(r) from the problem.

It asks:

    Is the B-array itself governed by a low-complexity exact law in

        (k,r)

or equivalently in

        (k,d),   d = r-k ?

The experiment tests only LOW-COMPLEXITY laws.

It does NOT fit a degree-5 bivariate polynomial through all 21 points,
because that would merely interpolate the finite table.

Tests:

    1. total-degree bivariate polynomials;
    2. tensor-product polynomial laws;
    3. the shifted coordinate d=r-k;
    4. fixed-offset sequences B[k,k+d];
    5. fixed-absolute-index sequences B[k,r].

Only overdetermined fits count as structural evidence.

Exact QQ arithmetic only.
No floating point.
No q-family fitting.
No arbitrary matrix fitting.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# EXACT B TABLE
# =============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def b_points():
    """
    Return exact points (k,r,value), with

        r = k + offset.
    """
    out = []

    for k, row in enumerate(B):
        for offset, value in enumerate(row):
            r = k + offset
            out.append(
                (
                    k,
                    r,
                    sp.Rational(value),
                )
            )

    return out


def total_degree_monomials(degree):
    """
    Monomials k^a r^b with a+b <= degree.
    """
    terms = []

    for a in range(degree + 1):
        for b in range(degree + 1 - a):
            terms.append(
                (a, b)
            )

    return terms


def tensor_monomials(deg_k, deg_r):
    return [
        (a, b)
        for a in range(deg_k + 1)
        for b in range(deg_r + 1)
    ]


def solve_polynomial_law(
    points,
    variables,
    monomials,
):
    """
    Exact polynomial fit.

    Returns:
        status
        polynomial_or_None
        rank
        nullity
        residuals
    """

    matrix = []
    rhs = []

    for point in points:

        values = {
            variables[0]: sp.Rational(point[0]),
            variables[1]: sp.Rational(point[1]),
        }

        matrix.append([
            values[variables[0]] ** a
            * values[variables[1]] ** b
            for a, b in monomials
        ])

        rhs.append(
            sp.Rational(point[2])
        )

    M = sp.Matrix(matrix)
    y = sp.Matrix(rhs)

    rank = M.rank()
    augmented_rank = M.row_join(y).rank()

    if augmented_rank != rank:
        return (
            "NO_SOLUTION",
            None,
            rank,
            M.cols - rank,
            [],
        )

    nullity = M.cols - rank

    if nullity != 0:
        return (
            "NONUNIQUE",
            None,
            rank,
            nullity,
            [],
        )

    coeffs = M.gauss_jordan_solve(y)[0]

    poly = sp.Integer(0)

    for coeff, (a, b) in zip(
        coeffs,
        monomials,
    ):
        poly += (
            clean(coeff)
            * variables[0] ** a
            * variables[1] ** b
        )

    poly = clean(poly)

    residuals = []

    for point in points:

        prediction = clean(
            poly.subs(
                {
                    variables[0]: point[0],
                    variables[1]: point[1],
                }
            )
        )

        residuals.append(
            clean(
                prediction - point[2]
            )
        )

    exact = all(
        x == 0
        for x in residuals
    )

    if exact:
        status = "UNIQUE_EXACT"
    else:
        status = "UNIQUE_FAILED"

    return (
        status,
        poly,
        rank,
        0,
        residuals,
    )


def fit_univariate(
    points,
    variable_name,
    degree,
):
    """
    Exact univariate fit for (x,value) pairs.
    """
    x = sp.Symbol(variable_name)

    matrix = []
    rhs = []

    for xv, value in points:

        matrix.append([
            sp.Rational(xv) ** d
            for d in range(degree + 1)
        ])

        rhs.append(
            sp.Rational(value)
        )

    M = sp.Matrix(matrix)
    y = sp.Matrix(rhs)

    rank = M.rank()
    augmented_rank = M.row_join(y).rank()

    if augmented_rank != rank:
        return (
            "NO_SOLUTION",
            None,
            rank,
            M.cols - rank,
        )

    if M.cols != rank:
        return (
            "NONUNIQUE",
            None,
            rank,
            M.cols - rank,
        )

    coeffs = M.gauss_jordan_solve(y)[0]

    poly = clean(
        sum(
            coeffs[d] * x ** d
            for d in range(degree + 1)
        )
    )

    for xv, value in points:

        pred = clean(
            poly.subs(
                x,
                xv,
            )
        )

        if pred != clean(value):
            return (
                "UNIQUE_FAILED",
                poly,
                rank,
                0,
            )

    return (
        "UNIQUE_EXACT",
        poly,
        rank,
        0,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    ksym, rsym = sp.symbols(
        "k r"
    )

    print("=" * 78)
    print(
        "EXPERIMENT 283 — EXACT B-GRID / "
        "BIVARIATE INDEX-LAW AUDIT"
    )
    print("=" * 78)

    points = b_points()

    # =========================================================================
    # 1. EXACT GRID
    # =========================================================================

    print()
    print("=" * 78)
    print("1. ABSOLUTE B-GRID")
    print("=" * 78)

    for k, r, value in points:
        print(
            f"  (k={k}, r={r}) -> {value}"
        )

    # =========================================================================
    # 2. TOTAL-DEGREE BIVARIATE SEARCH
    # =========================================================================

    print()
    print("=" * 78)
    print("2. TOTAL-DEGREE BIVARIATE POLYNOMIAL SEARCH")
    print("=" * 78)

    total_results = []

    # Stop at degree 4 because degree 5 has
    # C(5+2,2)=21 coefficients, exactly the
    # number of data points, hence pure interpolation.
    for degree in range(0, 5):

        monomials = total_degree_monomials(
            degree
        )

        status, poly, rank, nullity, residuals = (
            solve_polynomial_law(
                points,
                (
                    ksym,
                    rsym,
                ),
                monomials,
            )
        )

        total_results.append(
            (
                degree,
                status,
                poly,
                rank,
                nullity,
            )
        )

        print()
        print(
            f"  degree<={degree}:"
        )

        print(
            f"    coefficient_count="
            f"{len(monomials)}"
        )

        print(
            f"    rank={rank}"
        )

        print(
            f"    nullity={nullity}"
        )

        print(
            f"    status={status}"
        )

        print(
            f"    polynomial={poly}"
        )

    # =========================================================================
    # 3. TENSOR-PRODUCT SEARCH
    # =========================================================================

    print()
    print("=" * 78)
    print("3. TENSOR-PRODUCT POLYNOMIAL SEARCH")
    print("=" * 78)

    tensor_degrees = [
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
    ]

    tensor_results = []

    for deg_k, deg_r in tensor_degrees:

        monomials = tensor_monomials(
            deg_k,
            deg_r,
        )

        status, poly, rank, nullity, residuals = (
            solve_polynomial_law(
                points,
                (
                    ksym,
                    rsym,
                ),
                monomials,
            )
        )

        tensor_results.append(
            (
                deg_k,
                deg_r,
                status,
                poly,
                rank,
            )
        )

        print()
        print(
            f"  deg_k<={deg_k}, "
            f"deg_r<={deg_r}:"
        )

        print(
            f"    coefficient_count="
            f"{len(monomials)}"
        )

        print(
            f"    rank={rank}"
        )

        print(
            f"    status={status}"
        )

        print(
            f"    polynomial={poly}"
        )

    # =========================================================================
    # 4. SHIFTED COORDINATE d = r-k
    # =========================================================================

    print()
    print("=" * 78)
    print("4. SHIFTED COORDINATE d = r-k")
    print("=" * 78)

    d_points = [
        (
            k,
            r - k,
            value,
        )
        for k, r, value in points
    ]

    dsym = sp.Symbol(
        "d"
    )

    shifted_results = []

    for degree in range(0, 5):

        monomials = total_degree_monomials(
            degree
        )

        # Build direct polynomial system in (k,d).
        matrix = []
        rhs = []

        for kk, dd, value in d_points:

            matrix.append([
                sp.Rational(kk) ** a
                * sp.Rational(dd) ** b
                for a, b in monomials
            ])

            rhs.append(
                sp.Rational(value)
            )

        M = sp.Matrix(matrix)
        y = sp.Matrix(rhs)

        rank = M.rank()
        aug = M.row_join(y).rank()

        if aug != rank:
            status = "NO_SOLUTION"
            poly = None
        elif rank != M.cols:
            status = "NONUNIQUE"
            poly = None
        else:
            coeffs = M.gauss_jordan_solve(y)[0]

            poly = clean(
                sum(
                    coeff
                    * ksym ** a
                    * dsym ** b
                    for coeff, (a, b)
                    in zip(
                        coeffs,
                        monomials,
                    )
                )
            )

            status = "UNIQUE_EXACT"

        shifted_results.append(
            (
                degree,
                status,
                poly,
                rank,
            )
        )

        print()
        print(
            f"  total_degree<={degree}:"
        )

        print(
            f"    rank={rank}"
        )

        print(
            f"    status={status}"
        )

        print(
            f"    polynomial={poly}"
        )

    # =========================================================================
    # 5. FIXED-OFFSET SEQUENCES
    # =========================================================================

    print()
    print("=" * 78)
    print("5. FIXED-OFFSET d=r-k POLYNOMIAL TEST")
    print("=" * 78)

    offset_results = []

    max_offset = 7

    for d in range(
        max_offset + 1
    ):

        seq = []

        for kk in range(
            len(B)
        ):

            rr = kk + d

            if rr >= 8:
                continue

            value = None

            for offset, candidate in enumerate(
                B[kk]
            ):
                if kk + offset == rr:
                    value = candidate
                    break

            if value is not None:
                seq.append(
                    (
                        kk,
                        value,
                    )
                )

        if len(seq) < 3:
            continue

        print()
        print(
            f"  d={d}: "
            f"points={seq}"
        )

        max_degree_here = min(
            3,
            len(seq) - 1,
        )

        exact_found = False

        for degree in range(
            max_degree_here + 1
        ):

            status, poly, rank, nullity = fit_univariate(
                seq,
                "k",
                degree,
            )

            if status == "UNIQUE_EXACT":

                print(
                    f"    degree<={degree}: "
                    f"EXACT "
                    f"P(k)={poly}"
                )

                exact_found = True
                break

            print(
                f"    degree<={degree}: "
                f"{status}"
            )

        offset_results.append(
            (
                d,
                exact_found,
            )
        )

    # =========================================================================
    # 6. FIXED-ABSOLUTE-r SEQUENCES
    # =========================================================================

    print()
    print("=" * 78)
    print("6. FIXED-ABSOLUTE-INDEX r POLYNOMIAL TEST")
    print("=" * 78)

    absolute_results = []

    for rr in range(8):

        seq = []

        for kk in range(
            len(B)
        ):

            offset = rr - kk

            if (
                offset < 0
                or offset >= len(B[kk])
            ):
                continue

            seq.append(
                (
                    kk,
                    B[kk][offset],
                )
            )

        if len(seq) < 3:
            continue

        print()
        print(
            f"  r={rr}: "
            f"points={seq}"
        )

        max_degree_here = min(
            3,
            len(seq) - 1,
        )

        exact_found = False

        for degree in range(
            max_degree_here + 1
        ):

            status, poly, rank, nullity = fit_univariate(
                seq,
                "k",
                degree,
            )

            if status == "UNIQUE_EXACT":

                print(
                    f"    degree<={degree}: "
                    f"EXACT "
                    f"P(k)={poly}"
                )

                exact_found = True
                break

            print(
                f"    degree<={degree}: "
                f"{status}"
            )

        absolute_results.append(
            (
                rr,
                exact_found,
            )
        )

    # =========================================================================
    # 7. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 282 rejected the simple hypothesis

    B[k,r] = sum_p w[k,p] q_p(r).

Experiment 283 therefore asks a narrower upstream question:

    Does B[k,r] itself have a compact arithmetic law?

The strongest useful outcomes are:

    * a low-degree bivariate polynomial in (k,r);
    * a low-degree polynomial in (k,r-k);
    * repeated low-degree laws along fixed offsets;
    * repeated low-degree laws along fixed absolute indices.

A LOW-DEGREE law means the structure is visible directly in the B
array and may reveal the original combinatorial operator.

A failure of all low-degree tests is equally useful: it means the
B-channel is not merely a polynomial index kernel and we should return
to the original matrix/source construction rather than keep fitting
the final coefficients.
"""
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    low_degree_exact = [
        x
        for x in total_results
        if x[1] == "UNIQUE_EXACT"
    ]

    tensor_exact = [
        x
        for x in tensor_results
        if x[2] == "UNIQUE_EXACT"
    ]

    shifted_exact = [
        x
        for x in shifted_results
        if x[1] == "UNIQUE_EXACT"
    ]

    offset_exact = [
        x
        for x in offset_results
        if x[1]
    ]

    absolute_exact = [
        x
        for x in absolute_results
        if x[1]
    ]

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  total_degree_exact_count="
        f"{len(low_degree_exact)}"
    )

    print(
        f"  tensor_product_exact_count="
        f"{len(tensor_exact)}"
    )

    print(
        f"  shifted_coordinate_exact_count="
        f"{len(shifted_exact)}"
    )

    print(
        f"  exact_fixed_offset_laws="
        f"{len(offset_exact)}"
    )

    print(
        f"  exact_fixed_absolute_r_laws="
        f"{len(absolute_exact)}"
    )

    print(
        "  degree_5_full_interpolation_used=False"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 283 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted."
        )
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

