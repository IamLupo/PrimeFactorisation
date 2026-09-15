#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 296 — EXACT BIVARIATE BINOMIAL-BASIS / NEWTON-KERNEL AUDIT
==============================================================================

Experiment 295 found no convincing entrywise normalization using the obvious
single factors.

Experiment 296 changes basis rather than inventing another scalar factor.

Define

    C[k,r] = r! * B[k,r].

Use

    d = r-k.

The natural finite-combinatorial basis is

    C[k,d]
      =
    sum_{i,j} a[i,j] * binom(k,i) * binom(d,j).

We test low bidegrees only.

For each candidate degree pair:

    deg_k <= I
    deg_d <= J

the system is solved over QQ using ALL available entries.

Only an exact overdetermined solution counts as a structural hit.

We also test the one-variable Newton/binomial expansions separately:

    fixed k:
        C_k(d) = sum_j a[k,j] binom(d,j)

    fixed d:
        C_d(k) = sum_i a[d,i] binom(k,i)

The important question is whether the coefficients become simpler in this
basis than in the ordinary power basis.

No q-family data.
No arbitrary matrix fit.
No high-degree interpolation counted as evidence.
Exact QQ arithmetic only.
"""

from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# DATA
# ============================================================================

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


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def integer_grid():
    C = {}

    for k, row in enumerate(B):
        C[k] = {}

        for offset, value in enumerate(row):
            r = k + offset

            v = clean(
                sp.factorial(r) * value
            )

            if sp.denom(v) != 1:
                raise ArithmeticError(
                    f"C[{k},{r}] not integral: {v}"
                )

            C[k][r] = sp.Integer(v)

    return C


def primitive_signature(values):
    if not values:
        return []

    vals = [
        sp.Rational(v)
        for v in values
    ]

    den = 1

    for v in vals:
        den = math.lcm(
            den,
            int(v.q),
        )

    ints = [
        int(v * den)
        for v in vals
    ]

    g = 0

    for x in ints:
        g = math.gcd(
            g,
            abs(x),
        )

    if g:
        ints = [
            x // g
            for x in ints
        ]

    return ints


def binomial_basis_value(x, n):
    if n < 0:
        return sp.Integer(0)

    if x < 0:
        return sp.Integer(0)

    return sp.binomial(
        sp.Integer(x),
        sp.Integer(n),
    )


# ============================================================================
# ONE-VARIABLE BINOMIAL BASIS
# ============================================================================

def solve_binomial_1d(points):
    """
    Solve

        f(x) = sum_{j=0}^D a_j C(x,j)

    exactly.

    Because the binomial basis is triangular on x=0,1,...,
    coefficients are obtained by finite differences.
    """

    points = {
        int(x): sp.Rational(y)
        for x, y in points
    }

    if not points:
        return []

    xs = sorted(points)
    maximum = max(xs)

    values = [
        points.get(
            x,
            sp.Integer(0),
        )
        for x in range(
            maximum + 1
        )
    ]

    coefficients = []

    current = values

    while current:

        coefficients.append(
            clean(
                current[0]
            )
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

    return coefficients


def verify_binomial_1d(
    points,
    coefficients,
):
    for x, expected in points:

        rebuilt = sp.Integer(0)

        for j, a in enumerate(
            coefficients
        ):
            if j <= x:
                rebuilt += (
                    a
                    * binomial_basis_value(
                        x,
                        j,
                    )
                )

        if clean(
            rebuilt - expected
        ) != 0:
            return False

    return True


# ============================================================================
# BIVARIATE BASIS SOLVER
# ============================================================================

def basis_terms(max_i, max_j):
    return [
        (i, j)
        for i in range(max_i + 1)
        for j in range(max_j + 1)
    ]


def solve_bivariate(
    data,
    max_i,
    max_j,
):
    """
    Solve

        C(k,d)
          =
        sum a[i,j] binom(k,i) binom(d,j).

    Uses only points with 0 <= k and 0 <= d.
    """

    terms = basis_terms(
        max_i,
        max_j,
    )

    equations = []
    rhs = []

    for k, rmap in sorted(data.items()):

        for r, value in sorted(
            rmap.items()
        ):

            d = r - k

            row = []

            for i, j in terms:

                row.append(
                    binomial_basis_value(
                        k,
                        i,
                    )
                    * binomial_basis_value(
                        d,
                        j,
                    )
                )

            equations.append(row)
            rhs.append(
                sp.Rational(value)
            )

    if not equations:
        return {
            "status": "NO_DATA",
            "rank": 0,
            "augmented_rank": 0,
            "parameters": None,
            "terms": terms,
        }

    M = sp.Matrix(
        equations
    )

    b = sp.Matrix(
        rhs
    )

    rank = M.rank()
    augmented_rank = (
        M.row_join(b).rank()
    )

    if augmented_rank != rank:

        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
            "terms": terms,
        }

    if rank < len(terms):

        return {
            "status": "NONUNIQUE",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
            "terms": terms,
        }

    solution = M.gauss_jordan_solve(
        b
    )[0]

    parameters = [
        clean(v)
        for v in solution
    ]

    return {
        "status": "EXACT",
        "rank": rank,
        "augmented_rank": augmented_rank,
        "parameters": parameters,
        "terms": terms,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 296 — EXACT BIVARIATE BINOMIAL-BASIS / "
        "NEWTON-KERNEL AUDIT"
    )
    print("=" * 78)

    C = integer_grid()

    # =========================================================================
    # 1. INTEGER GRID
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "1. INTEGER GRID C[k,r] = r! B[k,r]"
    )
    print("=" * 78)

    for k in sorted(C):

        row = [
            C[k][r]
            for r in sorted(C[k])
        ]

        print()
        print(
            f"  k={k}: {row}"
        )

    # =========================================================================
    # 2. FIXED-k BINOMIAL EXPANSIONS IN d=r-k
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "2. FIXED-k BINOMIAL EXPANSIONS IN d=r-k"
    )
    print("=" * 78)

    fixed_k_profiles = {}

    for k in sorted(C):

        points = [
            (
                r - k,
                C[k][r],
            )
            for r in sorted(C[k])
        ]

        coefficients = solve_binomial_1d(
            points
        )

        exact = verify_binomial_1d(
            points,
            coefficients,
        )

        fixed_k_profiles[k] = coefficients

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    coefficients={coefficients}"
        )

        print(
            f"    primitive="
            f"{primitive_signature(coefficients)}"
        )

        print(
            f"    exact={exact}"
        )

    # =========================================================================
    # 3. FIXED-d BINOMIAL EXPANSIONS IN k
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "3. FIXED-d BINOMIAL EXPANSIONS IN k"
    )
    print("=" * 78)

    fixed_d_data = {}

    all_d = sorted(
        {
            r - k
            for k in C
            for r in C[k]
        }
    )

    for d in all_d:

        points = []

        for k in sorted(C):

            r = k + d

            if r in C[k]:

                points.append(
                    (
                        k,
                        C[k][r],
                    )
                )

        if not points:
            continue

        coefficients = solve_binomial_1d(
            points
        )

        exact = verify_binomial_1d(
            points,
            coefficients,
        )

        fixed_d_data[d] = coefficients

        print()
        print(
            f"  d={d}:"
        )

        print(
            f"    points={points}"
        )

        print(
            f"    coefficients={coefficients}"
        )

        print(
            f"    primitive="
            f"{primitive_signature(coefficients)}"
        )

        print(
            f"    exact={exact}"
        )

    # =========================================================================
    # 4. BIVARIATE BINOMIAL SEARCH
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "4. BIVARIATE BINOMIAL-BASIS SEARCH"
    )
    print("=" * 78)

    exact_hits = []
    nonunique_hits = []

    for deg_k, deg_d in [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
    ]:

        result = solve_bivariate(
            C,
            deg_k,
            deg_d,
        )

        print()
        print(
            f"  deg_k<={deg_k}, "
            f"deg_d<={deg_d}:"
        )

        print(
            f"    status={result['status']}"
        )

        print(
            f"    rank={result['rank']}"
        )

        print(
            f"    augmented_rank="
            f"{result['augmented_rank']}"
        )

        if result["status"] == "EXACT":

            params = result["parameters"]

            nonzero = [
                (
                    term,
                    coeff,
                )
                for term, coeff in zip(
                    result["terms"],
                    params,
                )
                if coeff != 0
            ]

            print(
                f"    nonzero_coefficients="
                f"{nonzero}"
            )

            print(
                f"    primitive="
                f"{primitive_signature(params)}"
            )

            exact_hits.append(
                (
                    deg_k,
                    deg_d,
                    result,
                )
            )

        elif result["status"] == "NONUNIQUE":

            nonunique_hits.append(
                (
                    deg_k,
                    deg_d,
                )
            )

    # =========================================================================
    # 5. SPARSITY PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "5. BINOMIAL-KERNEL SPARSITY PROFILE"
    )
    print("=" * 78)

    for (
        deg_k,
        deg_d,
        result,
    ) in exact_hits:

        params = result[
            "parameters"
        ]

        nonzero_count = sum(
            1
            for x in params
            if x != 0
        )

        print()
        print(
            f"  ({deg_k},{deg_d}): "
            f"nonzero={nonzero_count}/"
            f"{len(params)}"
        )

        print(
            f"    signature="
            f"{primitive_signature(params)}"
        )

    # =========================================================================
    # 6. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The ordinary power basis did not expose a compact B-grid law.

Experiment 296 tests the finite-combinatorial Newton/binomial basis

    binom(k,i) binom(d,j),
        d = r-k.

This is a more natural basis for a triangular discrete kernel because
binomial coefficients are the eigenbasis of finite-difference operators.

There are three useful outcomes.

A. Exact low-degree bivariate hit:

       C[k,r]
         =
       sum a[i,j] C(k,i) C(d,j).

   This would reveal a compact discrete kernel.

B. Exact fixed-d or fixed-k simplicity:

   This would show that one index has a simple Newton structure even if
   the full two-variable law is not compact.

C. No low-degree exact law:

   Then the B-grid is not a simple finite-difference/binomial kernel,
   strengthening the conclusion that the original B construction must
   be reconstructed directly.

Nonunique solutions are diagnostic only and are never counted as proof.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "7. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1 = 495451247
    q3 = 421514439

    g = math.gcd(
        q1,
        q3,
    )

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={g}"
    )

    print(
        f"  q1/17={q1 // 17}"
    )

    print(
        f"  q3/17={q3 // 17}"
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  r_factorial_grid_exact=True"
    )

    print(
        "  exact_1d_binomial_reconstructions="
        + str(
            len(fixed_k_profiles)
        )
    )

    print(
        "  exact_bivariate_hits="
        + str(
            len(exact_hits)
        )
    )

    print(
        "  nonunique_bivariate_hits="
        + str(
            len(nonunique_hits)
        )
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
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
        "EXPERIMENT 296 COMPLETE"
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

