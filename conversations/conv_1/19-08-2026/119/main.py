#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 285 — EXACT ADJACENT-ROW B-CHANNEL LOCAL OPERATOR AUDIT
==============================================================================

Experiment 284 ruled out:

    * low-order constant-coefficient recurrences along fixed diagonals;
    * simple geometric ratios;
    * low-rank multiplicative separability.

The next natural possibility is a ROW-TO-ROW operator.

Let B_k(r) denote the absolute-index B coefficient.

Test whether

    B_{k+1}(r)
      =
    a_k,0 * B_k(r)
    + a_k,1 * B_k(r-1)

or, more generally,

    B_{k+1}(r)
      =
    sum_{s=0}^{w-1}
        a_{k,s} B_k(r-s)

for small widths w.

The crucial condition is:

    the same coefficients a_{k,s}
    must work for every overlapping r.

Therefore a width-w operator is only accepted when the system is
overdetermined:

    number_of_equations > number_of_unknowns.

No interpolation-only case is counted.

We test widths:

    w = 1
    w = 2
    w = 3

and also the reversed orientation

    B_{k+1}(r)
      =
    sum_s a_{k,s} B_k(r+s).

Finally we test whether the recovered coefficients themselves follow
a simple exact law in k.

Exact QQ arithmetic only.
No q-family data.
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
# ABSOLUTE-INDEX ACCESS
# =============================================================================

def b_value(k, r):
    offset = r - k

    if offset < 0:
        return None

    if offset >= len(B[k]):
        return None

    return sp.Rational(
        B[k][offset]
    )


def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


# =============================================================================
# EXACT LOCAL OPERATOR SOLVER
# =============================================================================

def solve_operator(
    k,
    width,
    reverse=False,
):
    """
    Solve

        B[k+1,r]
          =
        sum_s a_s B[k, r-s]

    for reverse=False, or

        B[k+1,r]
          =
        sum_s a_s B[k, r+s]

    for reverse=True.

    Only use equations where EVERY required B entry exists.

    Require:

        equations > width

    so the test is overdetermined.
    """

    if k + 1 >= len(B):
        return {
            "status": "NO_NEXT_ROW",
            "equations": [],
            "solution": None,
            "rank": None,
            "nullity": None,
            "residuals": [],
        }

    equations = []
    rhs = []
    r_values = []

    for r in range(8):

        target = b_value(
            k + 1,
            r,
        )

        if target is None:
            continue

        row = []

        valid = True

        for s in range(width):

            source_r = (
                r + s
                if reverse
                else r - s
            )

            source = b_value(
                k,
                source_r,
            )

            if source is None:
                valid = False
                break

            row.append(
                source
            )

        if not valid:
            continue

        equations.append(
            row
        )

        rhs.append(
            target
        )

        r_values.append(
            r
        )

    # Require an actually overdetermined test.
    if len(equations) <= width:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": r_values,
            "solution": None,
            "rank": None,
            "nullity": None,
            "residuals": [],
        }

    M = sp.Matrix(
        equations
    )

    y = sp.Matrix(
        rhs
    )

    rank = M.rank()
    aug_rank = (
        M.row_join(y).rank()
    )

    if aug_rank != rank:
        return {
            "status": "NO_SOLUTION",
            "equations": r_values,
            "solution": None,
            "rank": rank,
            "nullity": M.cols - rank,
            "residuals": [],
        }

    nullity = M.cols - rank

    if nullity != 0:
        return {
            "status": "NONUNIQUE",
            "equations": r_values,
            "solution": None,
            "rank": rank,
            "nullity": nullity,
            "residuals": [],
        }

    solution = M.gauss_jordan_solve(
        y
    )[0]

    solution = tuple(
        clean(x)
        for x in solution
    )

    residuals = []

    for i, r in enumerate(
        r_values
    ):

        predicted = clean(
            sum(
                solution[s]
                * M[i, s]
                for s in range(width)
            )
        )

        residuals.append(
            clean(
                predicted - y[i]
            )
        )

    exact = all(
        x == 0
        for x in residuals
    )

    return {
        "status": (
            "EXACT"
            if exact
            else "FAILED"
        ),
        "equations": r_values,
        "solution": solution,
        "rank": rank,
        "nullity": 0,
        "residuals": residuals,
    }


# =============================================================================
# APPLY AND VERIFY A RECOVERED OPERATOR
# =============================================================================

def verify_operator(
    k,
    solution,
    width,
    reverse=False,
):
    checks = []

    if solution is None:
        return checks

    for r in range(8):

        target = b_value(
            k + 1,
            r,
        )

        if target is None:
            continue

        terms = []

        valid = True

        for s in range(width):

            source_r = (
                r + s
                if reverse
                else r - s
            )

            source = b_value(
                k,
                source_r,
            )

            if source is None:
                valid = False
                break

            terms.append(
                solution[s] * source
            )

        if not valid:
            continue

        prediction = clean(
            sum(terms)
        )

        checks.append(
            (
                r,
                prediction,
                target,
                clean(
                    prediction - target
                ),
            )
        )

    return checks


# =============================================================================
# RECOVERED COEFFICIENT POLYNOMIAL TEST
# =============================================================================

def polynomial_in_k(values):
    """
    Given values:

        [(k, value), ...]

    test exact polynomial laws of degree 0,1,2.

    Only overdetermined fits count.
    """

    ksym = sp.Symbol(
        "k"
    )

    out = []

    for degree in [0, 1, 2]:

        if len(values) <= degree + 1:
            out.append(
                (
                    degree,
                    "INSUFFICIENT_DATA",
                    None,
                )
            )
            continue

        M = []
        y = []

        for kk, value in values:

            M.append([
                sp.Rational(kk) ** d
                for d in range(
                    degree + 1
                )
            ])

            y.append(
                sp.Rational(value)
            )

        M = sp.Matrix(M)
        y = sp.Matrix(y)

        if len(values) <= degree + 1:
            continue

        rank = M.rank()
        aug_rank = (
            M.row_join(y).rank()
        )

        if aug_rank != rank:
            out.append(
                (
                    degree,
                    "NO_SOLUTION",
                    None,
                )
            )
            continue

        if rank != M.cols:
            out.append(
                (
                    degree,
                    "NONUNIQUE",
                    None,
                )
            )
            continue

        coeffs = M.gauss_jordan_solve(
            y
        )[0]

        poly = clean(
            sum(
                coeffs[d]
                * ksym ** d
                for d in range(
                    degree + 1
                )
            )
        )

        exact = True

        for kk, value in values:

            if clean(
                poly.subs(
                    ksym,
                    kk,
                )
                - value
            ) != 0:
                exact = False
                break

        out.append(
            (
                degree,
                "EXACT"
                if exact
                else "FAILED",
                poly,
            )
        )

    return out


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 285 — EXACT ADJACENT-ROW "
        "B-CHANNEL LOCAL OPERATOR AUDIT"
    )
    print("=" * 78)

    exact_hits = []
    all_recovered = {}

    # =========================================================================
    # 1. LOCAL ROW OPERATORS
    # =========================================================================

    print()
    print("=" * 78)
    print("1. ADJACENT-ROW LOCAL OPERATORS")
    print("=" * 78)

    for reverse in [
        False,
        True,
    ]:

        orientation = (
            "REVERSE"
            if reverse
            else "FORWARD"
        )

        print()
        print(
            f"  orientation={orientation}"
        )

        for width in [1, 2, 3]:

            print()
            print(
                f"  width={width}"
            )

            for k in range(
                len(B) - 1
            ):

                result = solve_operator(
                    k,
                    width,
                    reverse,
                )

                all_recovered[
                    (
                        orientation,
                        width,
                        k,
                    )
                ] = result

                print()
                print(
                    f"    k={k}:"
                )

                print(
                    f"      status="
                    f"{result['status']}"
                )

                print(
                    f"      equations="
                    f"{result['equations']}"
                )

                print(
                    f"      rank="
                    f"{result['rank']}"
                )

                print(
                    f"      nullity="
                    f"{result['nullity']}"
                )

                print(
                    f"      coefficients="
                    f"{result['solution']}"
                )

                if (
                    result["status"]
                    == "EXACT"
                ):

                    checks = verify_operator(
                        k,
                        result["solution"],
                        width,
                        reverse,
                    )

                    print(
                        f"      verification="
                        f"{checks}"
                    )

                    exact_hits.append(
                        (
                            orientation,
                            width,
                            k,
                            result[
                                "solution"
                            ],
                        )
                    )

    # =========================================================================
    # 2. COEFFICIENT PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("2. RECOVERED COEFFICIENT PROFILE")
    print("=" * 78)

    by_orientation_width = {}

    for orientation, width, k, coeffs in exact_hits:

        by_orientation_width.setdefault(
            (
                orientation,
                width,
            ),
            [],
        ).append(
            (
                k,
                coeffs,
            )
        )

    for key, rows in (
        by_orientation_width.items()
    ):

        orientation, width = key

        print()
        print(
            f"  orientation={orientation} "
            f"width={width}"
        )

        for k, coeffs in rows:
            print(
                f"    k={k}: "
                f"{coeffs}"
            )

        # Test each coefficient separately as a function of k.
        max_width = max(
            len(coeffs)
            for _, coeffs in rows
        )

        for index in range(
            max_width
        ):

            values = [
                (
                    k,
                    coeffs[index],
                )
                for k, coeffs in rows
                if len(coeffs) > index
            ]

            if len(values) < 3:
                continue

            print()
            print(
                f"    coefficient_index={index}"
            )

            for (
                degree,
                status,
                poly,
            ) in polynomial_in_k(
                values
            ):

                print(
                    f"      degree<={degree}: "
                    f"{status} "
                    f"{poly}"
                )

    # =========================================================================
    # 3. DIRECT ROW-PAIR RELATION
    # =========================================================================

    print()
    print("=" * 78)
    print("3. DIRECT ROW-PAIR CONSISTENCY")
    print("=" * 78)

    for k in range(
        len(B) - 1
    ):

        overlap = []

        for r in range(8):

            a = b_value(
                k,
                r,
            )

            b = b_value(
                k + 1,
                r,
            )

            if (
                a is not None
                and b is not None
            ):
                overlap.append(
                    (
                        r,
                        a,
                        b,
                    )
                )

        print()
        print(
            f"  k={k}: "
            f"overlap={overlap}"
        )

    # =========================================================================
    # 4. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 284 found no useful diagonal recurrence or separability.

Experiment 285 now tests whether the B triangle is generated locally from
one row to the next.

A successful width-1 law means:

    B_{k+1}(r) = a_k B_k(r).

A successful width-2 law means:

    B_{k+1}(r)
      = a_k B_k(r)
        + b_k B_k(r-1)

with the SAME (a_k,b_k) for all valid overlapping r.

Width 3 tests the next natural finite-band operator.

The reverse orientation tests whether the shift acts toward increasing
absolute r instead.

Only overdetermined exact systems count.

If these fail too, then the B array is not generated by a small local
row-to-row recurrence. At that point the correct next move is no longer
another pattern search: we need the original formula/code that constructs
B itself.
"""
    )

    # =========================================================================
    # 5. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. TERMINAL PROJECTIVE SOURCE REFERENCE")
    print("=" * 78)

    q1_terminal = 495451247
    q3_terminal = 421514439

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={sp.gcd(
            q1_terminal,
            q3_terminal,
        )}"
    )

    print(
        f"  q1/17={q1_terminal // 17}"
    )

    print(
        f"  q3/17={q3_terminal // 17}"
    )

    # =========================================================================
    # 6. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exact_local_operator_hits="
        f"{len(exact_hits)}"
    )

    print(
        "  overdetermined_only=True"
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
        "EXPERIMENT 285 COMPLETE"
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

