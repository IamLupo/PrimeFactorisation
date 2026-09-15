#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 90 — EXACT HYPERGEOMETRIC / BINOMIAL DIAGONAL LAW SEARCH
# ==============================================================================

k = sp.Symbol("k")


def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def fact(n):
    return sp.factorial(n)


def factor_ratio(expr):
    return sp.factor(sp.together(clean(expr)))


def integer_signature(values):
    vals = [sp.Rational(v) for v in values]

    lcm = sp.Integer(1)
    for q in vals:
        lcm = sp.ilcm(lcm, int(q.q))

    ints = [int(q * lcm) for q in vals]

    g = 0
    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return ints

    return [x // g for x in ints]


# ------------------------------------------------------------------------------
# Exact low-degree rational-ratio fitting
# ------------------------------------------------------------------------------

def ratio_fit(values, max_num_degree=2, max_den_degree=2):
    """
    Search for

        D(k+1)/D(k) = P(k)/Q(k)

    with small polynomial degrees.

    The fit is exact over QQ and is only claimed on the observed finite set.
    """

    ratios = []

    for i in range(len(values) - 1):
        x0, y0 = values[i]
        x1, y1 = values[i + 1]

        if y0 == 0:
            return None

        ratios.append((sp.Integer(x0), clean(y1 / y0)))

    if not ratios:
        return None

    for deg_num in range(max_num_degree + 1):
        for deg_den in range(max_den_degree + 1):

            # Q is normalized to constant coefficient 1.
            num_vars = deg_num + 1
            den_vars = deg_den

            unknown_count = num_vars + den_vars

            if len(ratios) < unknown_count:
                continue

            a = sp.symbols(
                f"a0:{num_vars}"
            )
            b = sp.symbols(
                f"b0:{den_vars}"
            )

            equations = []

            for x, y in ratios:
                numerator = sum(
                    a[r] * x**r
                    for r in range(num_vars)
                )

                denominator = sp.Integer(1)

                for r in range(den_vars):
                    denominator += b[r] * x**(r + 1)

                equations.append(
                    sp.expand(
                        numerator
                        - y * denominator
                    )
                )

            matrix, rhs = sp.linear_eq_to_matrix(
                equations,
                list(a) + list(b),
            )

            try:
                solution_set = sp.linsolve(
                    (matrix, rhs)
                )
            except Exception:
                continue

            if solution_set == sp.EmptySet:
                continue

            sols = list(solution_set)
            if not sols:
                continue

            sol_tuple = sols[0]

            # Reject underdetermined families.
            free_symbols = set()
            for expr in sol_tuple:
                free_symbols |= expr.free_symbols

            if free_symbols & set(
                list(a) + list(b)
            ):
                continue

            substitution = dict(
                zip(
                    list(a) + list(b),
                    sol_tuple,
                )
            )

            x_dummy = k

            numerator = clean(
                sum(
                    a[r] * x_dummy**r
                    for r in range(num_vars)
                ).subs(substitution)
            )

            denominator = clean(
                (
                    1
                    + sum(
                        b[r] * x_dummy**(r + 1)
                        for r in range(den_vars)
                    )
                ).subs(substitution)
            )

            if denominator == 0:
                continue

            candidate = clean(
                numerator / denominator
            )

            ok = True

            for x, y in ratios:
                lhs = clean(
                    candidate.subs(k, x)
                )

                if lhs != y:
                    ok = False
                    break

            if ok:
                return candidate

    return None


# ------------------------------------------------------------------------------
# Exact diagonal data
# ------------------------------------------------------------------------------

A_C = {
    0: [
        -2,
        -154,
        -818,
        sp.Rational(-1360, 3),
        sp.Rational(8435, 24),
        sp.Rational(-5851, 120),
        sp.Rational(-13373, 720),
        sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    1: [
        0,
        -550,
        -5015,
        -4734,
        sp.Rational(7879, 3),
        sp.Rational(-5147, 120),
        sp.Rational(-210877, 720),
        sp.Rational(83651, 720),
        sp.Rational(-1028053, 40320),
    ],
    2: [
        0,
        0,
        -7125,
        -14567,
        4635,
        sp.Rational(25508, 15),
        sp.Rational(-198919, 144),
        sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    3: [
        0,
        0,
        0,
        -11900,
        -1711,
        sp.Rational(28949, 6),
        sp.Rational(-31711, 15),
        sp.Rational(404513, 840),
        sp.Rational(-234707, 4032),
    ],
    4: [
        0,
        0,
        0,
        0,
        sp.Rational(-17875, 6),
        sp.Rational(51337, 30),
        sp.Rational(-52447, 180),
        sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    5: [
        0,
        0,
        0,
        0,
        0,
        sp.Rational(-1001, 12),
        sp.Rational(26687, 360),
        sp.Rational(-40921, 1260),
        sp.Rational(9389, 1008),
    ],
    6: [
        0,
        0,
        0,
        0,
        0,
        0,
        sp.Rational(-5, 72),
        sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
}


B_C = {
    0: [
        25,
        619,
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    1: [
        0,
        1750,
        8624,
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    2: [
        0,
        0,
        9690,
        10234,
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    3: [
        0,
        0,
        0,
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    4: [
        0,
        0,
        0,
        0,
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    5: [
        0,
        0,
        0,
        0,
        0,
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
}


# ------------------------------------------------------------------------------
# Diagonal extraction
# ------------------------------------------------------------------------------

def diagonal(matrix, shift):
    result = []

    for row in sorted(matrix):
        col = row + shift

        if col >= len(matrix[row]):
            continue

        value = sp.Rational(matrix[row][col])

        if value != 0:
            result.append(
                (row, clean(value))
            )

    return result


def print_diagonal(name, values):
    print(f"  {name}")
    print(
        "    values =",
        [v for _, v in values],
    )
    print(
        "    primitive_signature =",
        integer_signature(
            [v for _, v in values]
        ),
    )


# ------------------------------------------------------------------------------
# Diagonal analysis
# ------------------------------------------------------------------------------

def analyze_channel(name, matrix, max_shift=3):
    print()
    print("-" * 78)
    print(f"{name} CHANNEL")
    print("-" * 78)

    diagonals = {}

    for shift in range(max_shift + 1):
        values = diagonal(matrix, shift)
        diagonals[shift] = values

        print()
        print(f"{name} shift={shift}")
        print_diagonal("diagonal", values)

        if len(values) < 2:
            continue

        print("    successive_ratios:")

        for i in range(len(values) - 1):
            kk, v0 = values[i]
            _, v1 = values[i + 1]

            if v0 == 0:
                ratio = sp.nan
            else:
                ratio = clean(v1 / v0)

            print(
                f"      k={kk}: {factor_ratio(ratio)}"
            )

        fitted = ratio_fit(
            values,
            max_num_degree=2,
            max_den_degree=2,
        )

        print(
            "    low_degree_ratio =",
            (
                factor_ratio(fitted)
                if fitted is not None
                else "NONE"
            ),
        )

    return diagonals


# ------------------------------------------------------------------------------
# Corrected normalization functions
# ------------------------------------------------------------------------------

def norm_raw(value, kk):
    return clean(value)


def norm_times_k_factorial(value, kk):
    return clean(
        value * fact(kk)
    )


def norm_divide_k_factorial(value, kk):
    return clean(
        value / fact(kk)
    )


def norm_times_k_shift_factorial(value, kk, shift):
    return clean(
        value * fact(kk + shift)
    )


def factorial_normalization(name, diagonals):
    print()
    print("=" * 78)
    print(f"{name} FACTORIAL NORMALIZATION")
    print("=" * 78)

    for shift, values in diagonals.items():

        if not values:
            continue

        print()
        print(f"  shift={shift}")

        transforms = [
            (
                "raw",
                norm_raw,
            ),
            (
                "times_k_factorial",
                norm_times_k_factorial,
            ),
            (
                "divide_k_factorial",
                norm_divide_k_factorial,
            ),
            (
                "times_(k+s)_factorial",
                lambda value, kk, s=shift:
                    norm_times_k_shift_factorial(
                        value,
                        kk,
                        s,
                    ),
            ),
        ]

        for label, fn in transforms:

            transformed = []

            for kk, value in values:
                transformed.append(
                    (
                        kk,
                        fn(value, kk),
                    )
                )

            vals = [
                value
                for _, value in transformed
            ]

            print(
                f"    {label}:"
            )
            print(
                "      values =",
                vals,
            )

            fitted = ratio_fit(
                transformed,
                max_num_degree=2,
                max_den_degree=2,
            )

            print(
                "      ratio_fit =",
                (
                    factor_ratio(fitted)
                    if fitted is not None
                    else "NONE"
                ),
            )


# ------------------------------------------------------------------------------
# Cross-channel comparison
# ------------------------------------------------------------------------------

def cross_channel(A_diag, B_diag):
    print()
    print("=" * 78)
    print("CROSS-CHANNEL HYPERGEOMETRIC AUDIT")
    print("=" * 78)

    for shift in sorted(
        set(A_diag)
        & set(B_diag)
    ):

        A = dict(A_diag[shift])
        B = dict(B_diag[shift])

        common = sorted(
            set(A)
            & set(B)
        )

        values = []

        for kk in common:
            if B[kk] != 0:
                values.append(
                    (
                        kk,
                        clean(
                            A[kk] / B[kk]
                        ),
                    )
                )

        if len(values) < 2:
            continue

        print()
        print(
            f"  shift={shift}"
        )

        print(
            "    A/B =",
            [v for _, v in values],
        )

        fitted = ratio_fit(
            values,
            max_num_degree=2,
            max_den_degree=2,
        )

        print(
            "    A/B successive-ratio law =",
            (
                factor_ratio(fitted)
                if fitted is not None
                else "NONE"
            ),
        )


# ------------------------------------------------------------------------------
# Product reconstruction
# ------------------------------------------------------------------------------

def test_product_reconstruction(
    name,
    values,
    ratio_expr,
):
    if ratio_expr is None:
        return

    print()
    print(name)

    predicted = {
        values[0][0]: clean(
            values[0][1]
        )
    }

    for i in range(len(values) - 1):

        kk = values[i][0]

        predicted[kk + 1] = clean(
            predicted[kk]
            * ratio_expr.subs(
                k,
                kk,
            )
        )

    ok = True

    for kk, actual in values:
        result = clean(
            predicted[kk]
        )

        match = (
            result == actual
        )

        print(
            f"  k={kk}: exact={match}"
        )

        if not match:
            ok = False

    print(
        "  product reconstruction =",
        ok,
    )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 90 — EXACT HYPERGEOMETRIC / "
        "BINOMIAL DIAGONAL LAW SEARCH"
    )
    print("=" * 78)

    print()
    print("0. EXACT SETUP")
    print("  P_k(j) = sum_r C[k,r] j_under_r")
    print("  D_s(k) = C[k,k+s]")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    A_diag = analyze_channel(
        "A",
        A_C,
        max_shift=3,
    )

    B_diag = analyze_channel(
        "B",
        B_C,
        max_shift=3,
    )

    factorial_normalization(
        "A",
        A_diag,
    )

    factorial_normalization(
        "B",
        B_diag,
    )

    cross_channel(
        A_diag,
        B_diag,
    )

    print()
    print("=" * 78)
    print("EXACT PRODUCT-FORM RECONSTRUCTION")
    print("=" * 78)

    for name, diagonals in [
        ("A", A_diag),
        ("B", B_diag),
    ]:

        for shift, values in diagonals.items():

            if len(values) < 2:
                continue

            fitted = ratio_fit(
                values,
                max_num_degree=2,
                max_den_degree=2,
            )

            if fitted is not None:
                test_product_reconstruction(
                    f"{name} shift={shift}",
                    values,
                    fitted,
                )

    print()
    print("=" * 78)
    print("STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 88 established

      P_k(j)
        = sum_{r>=k} C[k,r] j_under_r.

  Experiment 89 isolated the first few diagonals

      D_s(k) = C[k,k+s].

  The present experiment tests a natural criterion for
  combinatorial coefficient families:

      D_s(k+1) / D_s(k).

  A low-degree rational expression in k would mean that
  the diagonal is hypergeometric and therefore admits an
  exact product representation.

  We also test factorial-normalized versions because a
  hidden construction may differ from a hypergeometric
  sequence by simple k! or (k+s)! factors.

  Positive result:
      exact product structure.

  Negative result:
      no small hypergeometric law at the tested complexity.

  The computation is finite and exact.
  No extrapolation is claimed.
  No factorization algorithm is inferred.
        """
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 90 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()