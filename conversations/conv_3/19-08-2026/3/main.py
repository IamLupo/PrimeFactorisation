import sympy as sp
from fractions import Fraction
from math import isqrt


# ============================================================
# EXACT SYMBOLIC SETUP
# ============================================================

t = sp.symbols("t")


def h16_sym(t):
    return (
        -9 * t**8
        -36 * t**7
        -84 * t**6
        -126 * t**5
        -126 * t**4
        -84 * t**3
        -36 * t**2
        -9 * t
        -1
    )


def h15_sym(t):
    return (
        88 * t**9
        +396 * t**8
        +1164 * t**7
        +2226 * t**6
        +2898 * t**5
        +2604 * t**4
        +1596 * t**3
        +639 * t**2
        +151 * t
        +16
    )


def h14_sym(t):
    return (
        -276 * t**10
        -1380 * t**9
        -5460 * t**8
        -13560 * t**7
        -23058 * t**6
        -27510 * t**5
        -23100 * t**4
        -13410 * t**3
        -5135 * t**2
        -1169 * t
        -120
    )


H16 = sp.expand(h16_sym(t))
H15 = sp.expand(h15_sym(t))
H14 = sp.expand(h14_sym(t))


# ============================================================
# R(t) AND ITS DERIVATIVE
# ============================================================

R_EXPR = sp.cancel((H15**2) / (H16 * H14))
R_DERIV = sp.cancel(sp.diff(R_EXPR, t))

R_NUM, R_DEN = sp.fraction(R_DERIV)

R_NUM = sp.Poly(sp.expand(R_NUM), t)
R_DEN = sp.Poly(sp.expand(R_DEN), t)


# ============================================================
# REAL ROOTS
# ============================================================

def real_roots_of_poly(poly):
    roots = sp.nroots(poly.as_expr(), n=40, maxsteps=500)

    real_roots = []

    for r in roots:
        real_part = sp.re(r)
        imag_part = sp.im(r)

        if abs(float(imag_part)) < 1e-25:
            real_roots.append(float(real_part))

    real_roots.sort()

    return real_roots


# ============================================================
# EXACT NUMERIC R(t) USING FRACTION
# ============================================================

def h16_frac(x):
    return (
        -9 * x**8
        -36 * x**7
        -84 * x**6
        -126 * x**5
        -126 * x**4
        -84 * x**3
        -36 * x**2
        -9 * x
        -1
    )


def h15_frac(x):
    return (
        88 * x**9
        +396 * x**8
        +1164 * x**7
        +2226 * x**6
        +2898 * x**5
        +2604 * x**4
        +1596 * x**3
        +639 * x**2
        +151 * x
        +16
    )


def h14_frac(x):
    return (
        -276 * x**10
        -1380 * x**9
        -5460 * x**8
        -13560 * x**7
        -23058 * x**6
        -27510 * x**5
        -23100 * x**4
        -13410 * x**3
        -5135 * x**2
        -1169 * x
        -120
    )


def layer_invariant(n, x):
    tx = Fraction(n, x)

    a = h16_frac(tx)
    b = h15_frac(tx)
    c = h14_frac(tx)

    if a == 0 or c == 0:
        raise ZeroDivisionError("Layer denominator vanished.")

    return (b * b) / (a * c)


# ============================================================
# KAPPA
# ============================================================

def kappa_from_x(n, x):
    s = x - 1

    return (
        -s * s
        + (n + 1) * s
        - n * n
        + n
    )


# ============================================================
# DISCRIMINANT
# ============================================================

def discriminant(n, x):
    s = x - 1
    return s * s - 4 * n


def is_square(n):
    if n < 0:
        return False

    r = isqrt(n)

    return r * r == n


# ============================================================
# TEST PRIME PAIRS
# ============================================================

TEST_PAIRS = (
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
)


# ============================================================
# RELATIVE X PERTURBATIONS
#
# Much wider than +/-1000.
#
# For each instance we test:
#
#   X * factor
#
# approximately below and above the true value.
# ============================================================

RELATIVE_FACTORS = (
    Fraction(1, 100),
    Fraction(1, 20),
    Fraction(1, 10),
    Fraction(1, 5),
    Fraction(1, 2),
    Fraction(3, 4),
    Fraction(9, 10),
    Fraction(19, 20),
    Fraction(99, 100),
    Fraction(101, 100),
    Fraction(21, 20),
    Fraction(11, 10),
    Fraction(5, 4),
    Fraction(3, 2),
    Fraction(2, 1),
    Fraction(5, 2),
    Fraction(5, 1),
)


# ============================================================
# SIGN
# ============================================================

def sign(x):
    if x < 0:
        return -1
    if x > 0:
        return 1
    return 0


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 110)
    print("EXPERIMENT: SYMBOLIC LAYER-MONOTONICITY / UNIVERSAL X ORACLE")
    print("=" * 110)

    print()
    print("OBJECTIVE")
    print("  Determine whether the layer invariant")
    print()
    print("      R(X) = R_layer(N / X)")
    print()
    print("  is monotone in X on the full relevant positive domain.")
    print()
    print("  The exact derivative is analyzed symbolically first.")
    print("  Then the result is tested against many independent semiprimes.")
    print()

    # --------------------------------------------------------
    # SYMBOLIC DERIVATIVE INFORMATION
    # --------------------------------------------------------

    print("=" * 110)
    print("1. EXACT R(t)")
    print("=" * 110)

    print()
    print("h16(t) =")
    print(H16)

    print()
    print("h15(t) =")
    print(H15)

    print()
    print("h14(t) =")
    print(H14)

    print()
    print("R(t) =")
    print(R_EXPR)

    print()
    print("=" * 110)
    print("2. EXACT DERIVATIVE")
    print("=" * 110)

    print()
    print("R'(t) numerator degree =", R_NUM.degree())
    print("R'(t) denominator degree =", R_DEN.degree())

    print()
    print("FACTORED NUMERATOR OF R'(t):")
    print(sp.factor(R_NUM.as_expr()))

    print()
    print("FACTORED DENOMINATOR OF R'(t):")
    print(sp.factor(R_DEN.as_expr()))

    # --------------------------------------------------------
    # REAL ROOTS OF DERIVATIVE NUMERATOR
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("3. REAL ROOTS OF DERIVATIVE NUMERATOR")
    print("=" * 110)

    derivative_roots = real_roots_of_poly(R_NUM)

    if not derivative_roots:
        print()
        print("  No real roots detected.")
    else:
        print()
        for i, root in enumerate(derivative_roots, start=1):
            print(f"  root[{i}] = {root:.40g}")

    # --------------------------------------------------------
    # REAL ROOTS OF DENOMINATOR
    # --------------------------------------------------------

    denominator_roots = real_roots_of_poly(R_DEN)

    print()
    print("=" * 110)
    print("4. REAL ROOTS OF R'(t) DENOMINATOR")
    print("=" * 110)

    if not denominator_roots:
        print()
        print("  No real denominator roots detected.")
    else:
        print()
        for i, root in enumerate(denominator_roots, start=1):
            print(f"  denominator_root[{i}] = {root:.40g}")

    # --------------------------------------------------------
    # SIGN SAMPLING BETWEEN ROOTS
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("5. SIGN OF R'(t) ON POSITIVE AXIS")
    print("=" * 110)

    positive_breaks = [0.0]

    for root in derivative_roots:
        if root > 0:
            positive_breaks.append(root)

    for root in denominator_roots:
        if root > 0:
            positive_breaks.append(root)

    positive_breaks = sorted(set(positive_breaks))

    # Add large endpoint.
    positive_breaks.append(float("inf"))

    intervals = []

    for i in range(len(positive_breaks) - 1):
        left = positive_breaks[i]
        right = positive_breaks[i + 1]

        if right == float("inf"):
            if left == 0:
                sample = 1.0
            else:
                sample = left * 2.0
        else:
            if left == 0:
                sample = right / 2.0
            else:
                sample = (left + right) / 2.0

        try:
            deriv_value = float(
                sp.N(
                    R_DERIV.subs(t, sample),
                    30
                )
            )

            deriv_sign = sign(deriv_value)

        except Exception:
            deriv_value = float("nan")
            deriv_sign = 0

        intervals.append(
            (left, right, sample, deriv_sign, deriv_value)
        )

    print()

    for i, (left, right, sample, sgn, value) in enumerate(
        intervals,
        start=1
    ):

        left_text = "0" if left == 0 else f"{left:.12g}"
        right_text = (
            "+inf"
            if right == float("inf")
            else f"{right:.12g}"
        )

        print(
            f"  interval {i:2d}: "
            f"({left_text}, {right_text})"
        )

        print(
            f"      sample = {sample:.12g}"
        )

        print(
            f"      R'(t) sign = {sgn:+d}"
        )

        print(
            f"      R'(t) value = {value:.12g}"
        )

    # --------------------------------------------------------
    # CHECK WHETHER DERIVATIVE IS POSITIVE FOR ALL t > 0
    # --------------------------------------------------------

    derivative_positive_all_positive_t = True

    for left, right, sample, sgn, value in intervals:
        if sgn <= 0:
            derivative_positive_all_positive_t = False

    print()
    print("=" * 110)
    print("6. GLOBAL POSITIVE-T DERIVATIVE TEST")
    print("=" * 110)

    print()
    print(
        "  R'(t) appears positive on every sampled positive interval:",
        derivative_positive_all_positive_t
    )

    # --------------------------------------------------------
    # SEMIPRIME TESTS
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("7. SEMIPRIME MONOTONICITY TEST")
    print("=" * 110)

    total_directional_tests = 0
    directional_failures = 0

    total_relative_tests = 0
    relative_failures = 0

    derivative_domain_failures = 0

    for instance_id, (p, q) in enumerate(TEST_PAIRS, start=1):

        n = p * q
        s_true = p + q
        x_true = s_true + 1

        t_true = Fraction(n, x_true)

        r_true = layer_invariant(n, x_true)
        k_true = kappa_from_x(n, x_true)

        print()
        print("-" * 110)
        print(f"INSTANCE {instance_id}")
        print("-" * 110)

        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {n}")
        print(f"  X_TRUE  = {x_true}")

        print(
            f"  t_TRUE  = {t_true.numerator}/{t_true.denominator}"
        )

        t_true_float = float(t_true)

        print(
            f"  t_TRUE(float) = {t_true_float:.15g}"
        )

        # ----------------------------------------------------
        # Determine derivative sign at TRUE t.
        # ----------------------------------------------------

        derivative_at_true = sp.N(
            R_DERIV.subs(
                t,
                sp.Rational(t_true.numerator, t_true.denominator)
            ),
            40
        )

        derivative_sign_true = sign(float(derivative_at_true))

        print(
            f"  sign R'(t_TRUE) = {derivative_sign_true:+d}"
        )

        if derivative_sign_true <= 0:
            derivative_domain_failures += 1

        # ----------------------------------------------------
        # Relative directional tests.
        # ----------------------------------------------------

        print()
        print(
            f"  {'factor':>10} "
            f"{'X_guess':>14} "
            f"{'X_side':>8} "
            f"{'layer_sign':>12} "
            f"{'kappa_sign':>12} "
            f"{'disc_sq':>10}"
        )

        print("  " + "-" * 85)

        for factor in RELATIVE_FACTORS:

            # Below true X:
            x_below = int(
                Fraction(x_true) * factor
            )

            # Above true X:
            # reciprocal-style symmetric factor
            if factor == 0:
                continue

            x_above = int(
                Fraction(x_true) * (Fraction(1, factor))
            )

            # Avoid accidental equality.
            guesses = (
                ("BELOW", x_below),
                ("ABOVE", x_above),
            )

            for side, x_guess in guesses:

                if x_guess <= 0 or x_guess == x_true:
                    continue

                total_relative_tests += 1

                try:
                    r_guess = layer_invariant(n, x_guess)
                    dr = r_guess - r_true
                except ZeroDivisionError:
                    relative_failures += 1
                    continue

                kr = kappa_from_x(n, x_guess) - k_true

                actual_side = (
                    "BELOW"
                    if x_guess < x_true
                    else "ABOVE"
                )

                expected_layer_sign = (
                    1
                    if x_guess < x_true
                    else -1
                )

                expected_kappa_sign = (
                    -1
                    if x_guess < x_true
                    else 1
                )

                if sign(dr) != expected_layer_sign:
                    relative_failures += 1

                if sign(kr) != expected_kappa_sign:
                    relative_failures += 1

                total_directional_tests += 1

                square = is_square(
                    discriminant(n, x_guess)
                )

                print(
                    f"  {str(factor):>10} "
                    f"{x_guess:14d} "
                    f"{actual_side:>8} "
                    f"{sign(dr):12d} "
                    f"{sign(kr):12d} "
                    f"{str(square):>10}"
                )

        # ----------------------------------------------------
        # Local finite-difference monotonicity test.
        # Test every integer X in a relative window around TRUE.
        # ----------------------------------------------------

        local_radius = max(
            100,
            min(5000, x_true // 20)
        )

        print()
        print(
            f"  Integer-window test: +/- {local_radius}"
        )

        local_failures = 0

        for dx in range(
            -local_radius,
            local_radius + 1
        ):

            if dx == 0:
                continue

            x_guess = x_true + dx

            r_guess = layer_invariant(
                n,
                x_guess
            )

            dr = r_guess - r_true

            expected = (
                1 if dx < 0 else -1
            )

            if sign(dr) != expected:
                local_failures += 1

        print(
            f"  Local layer sign failures = {local_failures}"
        )

        if local_failures != 0:
            directional_failures += local_failures

    # --------------------------------------------------------
    # GLOBAL SUMMARY
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("8. GLOBAL SUMMARY")
    print("=" * 110)

    print()
    print(
        "  Total directional tests =",
        total_directional_tests
    )

    print(
        "  Relative tests          =",
        total_relative_tests
    )

    print(
        "  Relative-sign failures  =",
        relative_failures
    )

    print(
        "  Local-window failures   =",
        directional_failures
    )

    print(
        "  Derivative-at-true failures =",
        derivative_domain_failures
    )

    print()
    print(
        "  DERIVATIVE POSITIVE ON ALL SAMPLED POSITIVE INTERVALS :",
        derivative_positive_all_positive_t
    )

    print(
        "  ALL RELATIVE SIGN TESTS PASS :",
        relative_failures == 0
    )

    print(
        "  ALL LOCAL WINDOW TESTS PASS  :",
        directional_failures == 0
    )

    print(
        "  ALL TRUE-T DERIVATIVE TESTS PASS :",
        derivative_domain_failures == 0
    )

    # --------------------------------------------------------
    # FINAL INTERPRETATION
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("9. INTERPRETATION")
    print("=" * 110)

    if derivative_positive_all_positive_t:
        print()
        print("  RESULT:")
        print("    The exact derivative appears strictly positive over")
        print("    every sampled positive t-interval.")
        print()
        print("    Since t = N/X and dt/dX < 0, this is consistent with")
        print("    R_layer(X) being strictly decreasing in X.")
    else:
        print()
        print("  RESULT:")
        print("    R'(t) changes sign or reaches a non-positive interval.")
        print("    A universal monotonicity theorem cannot be inferred")
        print("    from the current symbolic sign analysis.")

    print()
    print("  IMPORTANT:")
    print("    This experiment still assumes access to the reference")
    print("    invariant R_TRUE generated from the hidden X_TRUE.")
    print("    It tests the mathematical directional property, not")
    print("    independent accessibility of the invariant from N.")

    print()
    print("=" * 110)
    print("END")
    print("=" * 110)


if __name__ == "__main__":
    main()
