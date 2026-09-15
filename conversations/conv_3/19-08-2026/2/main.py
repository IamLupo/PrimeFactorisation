from fractions import Fraction
from math import isqrt


# ============================================================
# EXACT HOMOGENEOUS-LAYER POLYNOMIALS
# ============================================================

def h16(t: Fraction) -> Fraction:
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


def h15(t: Fraction) -> Fraction:
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


def h14(t: Fraction) -> Fraction:
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


def layer_invariant(n: int, x: int) -> Fraction:
    t = Fraction(n, x)

    a = h16(t)
    b = h15(t)
    c = h14(t)

    if a == 0 or c == 0:
        raise ZeroDivisionError("Layer invariant denominator vanished.")

    return (b * b) / (a * c)


# ============================================================
# KAPPA FOR TWO FACTORS, EXPRESSED ONLY THROUGH N AND X
# X = p + q + 1
# ============================================================

def kappa_from_x(n: int, x: int) -> int:
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

def discriminant(n: int, x: int) -> int:
    s = x - 1
    return s * s - 4 * n


def is_square(n: int) -> bool:
    if n < 0:
        return False

    r = isqrt(n)
    return r * r == n


# ============================================================
# SIMPLE PRIMALITY TEST
# Enough for the test values below.
# ============================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # More than sufficient for the test integers used here.
    bases = (2, 3, 5, 7, 11, 13, 17)

    for a in bases:
        if a >= n:
            continue

        y = pow(a, d, n)

        if y == 1 or y == n - 1:
            continue

        for _ in range(s - 1):
            y = (y * y) % n

            if y == n - 1:
                break
        else:
            return False

    return True


# ============================================================
# FIXED PRIME TEST SET
#
# Deliberately varied sizes and balance.
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
# DELTA SET
#
# We test both:
#   small local perturbations
#   larger perturbations
#
# ============================================================

DELTAS = (
    1,
    2,
    5,
    10,
    25,
    50,
    100,
    250,
    500,
    1000,
)


# ============================================================
# SIGN
# ============================================================

def sign(x) -> int:
    if x < 0:
        return -1

    if x > 0:
        return 1

    return 0


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():
    print("=" * 110)
    print("EXPERIMENT: UNIVERSAL X-DIRECTIONAL ORACLE")
    print("=" * 110)

    print()
    print("GOAL")
    print("  For many independent semiprimes N=p*q, test whether:")
    print()
    print("    X_guess < X_true  =>  kappa residual < 0")
    print("    X_guess > X_true  =>  kappa residual > 0")
    print()
    print("  and independently:")
    print()
    print("    X_guess < X_true  =>  layer residual > 0")
    print("    X_guess > X_true  =>  layer residual < 0")
    print()
    print("  We also test whether the true X is the only nearby value")
    print("  producing a square discriminant.")
    print()

    total_cases = 0

    kappa_failures = 0
    layer_failures = 0
    discriminant_false_positives = 0

    all_instances_pass = True

    # --------------------------------------------------------
    # INSTANCE LOOP
    # --------------------------------------------------------

    for instance_id, (p, q) in enumerate(TEST_PAIRS, start=1):

        if not is_prime(p) or not is_prime(q):
            print()
            print(f"WARNING: test pair {p}, {q} is not prime.")
            continue

        n = p * q
        s_true = p + q
        x_true = s_true + 1

        k_true = kappa_from_x(n, x_true)
        r_true = layer_invariant(n, x_true)

        instance_kappa_ok = True
        instance_layer_ok = True
        instance_square_ok = True

        print()
        print("=" * 110)
        print(f"INSTANCE {instance_id}")
        print("=" * 110)

        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {n}")
        print(f"  S       = {s_true}")
        print(f"  X_TRUE  = {x_true}")

        print()
        print("  LOCAL SIGN TABLE")
        print(
            f"  {'delta':>8} "
            f"{'kappa':>8} "
            f"{'layer':>8} "
            f"{'disc_sq':>10}"
        )
        print("  " + "-" * 42)

        for delta in DELTAS:

            # -------------------------
            # Left guess
            # -------------------------

            x_left = x_true - delta

            k_left = kappa_from_x(n, x_left)
            r_left = layer_invariant(n, x_left)

            dk_left = k_left - k_true
            dr_left = r_left - r_true

            left_kappa_expected = dk_left < 0
            left_layer_expected = dr_left > 0

            if not left_kappa_expected:
                kappa_failures += 1
                instance_kappa_ok = False

            if not left_layer_expected:
                layer_failures += 1
                instance_layer_ok = False

            # -------------------------
            # Right guess
            # -------------------------

            x_right = x_true + delta

            k_right = kappa_from_x(n, x_right)
            r_right = layer_invariant(n, x_right)

            dk_right = k_right - k_true
            dr_right = r_right - r_true

            right_kappa_expected = dk_right > 0
            right_layer_expected = dr_right < 0

            if not right_kappa_expected:
                kappa_failures += 1
                instance_kappa_ok = False

            if not right_layer_expected:
                layer_failures += 1
                instance_layer_ok = False

            # -------------------------
            # Discriminant
            # -------------------------

            left_disc = discriminant(n, x_left)
            right_disc = discriminant(n, x_right)

            left_square = is_square(left_disc)
            right_square = is_square(right_disc)

            if left_square:
                discriminant_false_positives += 1
                instance_square_ok = False

            if right_square:
                discriminant_false_positives += 1
                instance_square_ok = False

            total_cases += 2

            print(
                f"  {-delta:8d} "
                f"{sign(dk_left):8d} "
                f"{sign(dr_left):8d} "
                f"{str(left_square):>10}"
            )

            print(
                f"  {delta:8d} "
                f"{sign(dk_right):8d} "
                f"{sign(dr_right):8d} "
                f"{str(right_square):>10}"
            )

        # ----------------------------------------------------
        # Exact local derivative check for kappa
        # ----------------------------------------------------

        # K(X+delta)-K(X) =
        # delta*(N+3-2X) - delta^2
        #
        # This gives an exact interval of monotonicity.

        turning_point = Fraction(n + 3, 2)

        kappa_monotone_right = x_true < turning_point

        print()
        print("  INSTANCE RESULT")
        print(f"    kappa orientation : {instance_kappa_ok}")
        print(f"    layer orientation : {instance_layer_ok}")
        print(f"    discr. uniqueness : {instance_square_ok}")
        print(
            f"    kappa increasing at X_TRUE : "
            f"{kappa_monotone_right}"
        )

        if not (
            instance_kappa_ok
            and instance_layer_ok
            and instance_square_ok
        ):
            all_instances_pass = False

    # --------------------------------------------------------
    # GLOBAL RESULT
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("GLOBAL RESULT")
    print("=" * 110)

    print(f"  Total directional tests       = {total_cases}")
    print(f"  Kappa sign failures           = {kappa_failures}")
    print(f"  Layer sign failures           = {layer_failures}")
    print(
        "  Nearby square-discriminant "
        f"false positives = {discriminant_false_positives}"
    )

    print()
    print(f"  KAPPA ORACLE PASS : {kappa_failures == 0}")
    print(f"  LAYER ORACLE PASS : {layer_failures == 0}")
    print(
        "  DISCRIMINANT PASS : "
        f"{discriminant_false_positives == 0}"
    )

    print()
    print(f"  ALL TESTS PASS    : {all_instances_pass}")

    # --------------------------------------------------------
    # FINAL INTERPRETATION
    # --------------------------------------------------------

    print()
    print("=" * 110)
    print("INTERPRETATION")
    print("=" * 110)

    if kappa_failures == 0:
        print("  KAPPA:")
        print("    The tested semiprimes all obey the same directional")
        print("    kappa sign rule around X_TRUE.")
    else:
        print("  KAPPA:")
        print("    The directional rule has counterexamples.")

    if layer_failures == 0:
        print("  LAYER:")
        print("    The tested semiprimes all obey the same directional")
        print("    homogeneous-layer sign rule.")
    else:
        print("  LAYER:")
        print("    The layer directional rule has counterexamples.")

    if discriminant_false_positives == 0:
        print("  DISCRIMINANT:")
        print("    No tested X_TRUE +/- delta produced another square")
        print("    discriminant.")
    else:
        print("  DISCRIMINANT:")
        print("    Additional square discriminants occurred.")

    print()
    print("  IMPORTANT:")
    print("    This experiment tests the directional property only.")
    print("    It does not test whether the reference invariants can")
    print("    be computed from N without knowing p and q.")

    print()
    print("=" * 110)
    print("END")
    print("=" * 110)


if __name__ == "__main__":
    main()
