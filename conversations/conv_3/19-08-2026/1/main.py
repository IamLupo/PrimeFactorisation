from fractions import Fraction
from math import isqrt


# ============================================================
# Known test instance
# ============================================================

P = 50387
Q = 282589

N = P * Q
S_TRUE = P + Q
X_TRUE = S_TRUE + 1


# ============================================================
# KAPPA
#
# For N = p*q and X = p+q+1:
#
#   kappa_2 =
#       -S^2 + (N+1)S - N^2 + N
#
# with S = X-1.
# ============================================================

def kappa_from_x(n: int, x: int) -> int:
    s = x - 1
    return (
        -s * s
        + (n + 1) * s
        - n * n
        + n
    )


KAPPA_TRUE = kappa_from_x(N, X_TRUE)


# ============================================================
# Homogeneous-layer invariant
#
# R(t) = h15(t)^2 / (h16(t) * h14(t))
#
# t = N / X
#
# The h-polynomials are taken from the exact reconstructed
# G_{9,16} homogeneous-layer formulas.
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


def layer_invariant_from_x(n: int, x: int) -> Fraction:
    t = Fraction(n, x)

    a = h16(t)
    b = h15(t)
    c = h14(t)

    if a == 0 or c == 0:
        raise ZeroDivisionError("Layer invariant denominator vanished.")

    return (b * b) / (a * c)


R_TRUE = layer_invariant_from_x(N, X_TRUE)


# ============================================================
# Residuals
#
# residual < 0 : guess is below/above according to orientation
# residual = 0 : exact match
# residual > 0 : opposite side
#
# For kappa, we define:
#
#   residual = guessed_value - true_value
#
# For layer invariant we do the same.
# ============================================================

def kappa_residual(x_guess: int) -> int:
    return kappa_from_x(N, x_guess) - KAPPA_TRUE


def layer_residual(x_guess: int) -> Fraction:
    return layer_invariant_from_x(N, x_guess) - R_TRUE


# ============================================================
# Discriminant feasibility
#
# A valid X must give:
#
#   S = X-1
#   Delta = S^2 - 4N
#
# and Delta must be a non-negative perfect square.
# ============================================================

def discriminant(x_guess: int) -> int:
    s = x_guess - 1
    return s * s - 4 * N


def square_discriminant(x_guess: int):
    delta = discriminant(x_guess)

    if delta < 0:
        return False, None

    root = isqrt(delta)

    if root * root != delta:
        return False, root

    return True, root


# ============================================================
# Candidate prime reconstruction
# ============================================================

def reconstruct_factors(x_guess: int):
    valid_square, root = square_discriminant(x_guess)

    if not valid_square:
        return None, None

    s = x_guess - 1

    if (s + root) % 2 != 0 or (s - root) % 2 != 0:
        return None, None

    p = (s + root) // 2
    q = (s - root) // 2

    return p, q


# ============================================================
# Simple deterministic primality test
#
# This is only used for the small test instance / candidate check.
# ============================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic bases for the integer range used here.
    bases = (2, 3, 5, 7, 11, 13, 17)

    for a in bases:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


# ============================================================
# Print helpers
# ============================================================

def sign_of(value) -> str:
    if value < 0:
        return "-"
    if value > 0:
        return "+"
    return "0"


def compact_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)

    return f"{value.numerator}/{value.denominator}"


# ============================================================
# Main experiment
# ============================================================

def main():
    print("=" * 100)
    print("X-DIRECTIONAL ORACLE EXPERIMENT")
    print("=" * 100)

    print()
    print("KNOWN INSTANCE")
    print(f"  p       = {P}")
    print(f"  q       = {Q}")
    print(f"  N       = {N}")
    print(f"  S       = {S_TRUE}")
    print(f"  X_TRUE  = {X_TRUE}")

    print()
    print("TRUE INVARIANTS")
    print(f"  kappa_2 = {KAPPA_TRUE}")
    print(f"  R_layer = {compact_fraction(R_TRUE)}")

    # --------------------------------------------------------
    # Deltas around the true X.
    # --------------------------------------------------------

    deltas = (
        -10000,
        -5000,
        -2000,
        -1000,
        -500,
        -250,
        -100,
        -50,
        -25,
        -10,
        -5,
        -2,
        -1,
        0,
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
        2000,
        5000,
        10000,
    )

    print()
    print("=" * 100)
    print("DIRECTION TEST")
    print("=" * 100)

    print(
        f"{'delta':>8} "
        f"{'X_guess':>10} "
        f"{'kappa_sign':>12} "
        f"{'layer_sign':>12} "
        f"{'delta_discr':>16} "
        f"{'square?':>8} "
        f"{'factors':>28}"
    )

    print("-" * 100)

    for delta in deltas:
        x_guess = X_TRUE + delta

        k_res = kappa_residual(x_guess)
        r_res = layer_residual(x_guess)

        delta_guess = discriminant(x_guess)
        square, root = square_discriminant(x_guess)

        factors = "-"

        if square:
            p_guess, q_guess = reconstruct_factors(x_guess)

            if p_guess is not None:
                prime_status = (
                    is_prime(p_guess)
                    and is_prime(q_guess)
                )

                factors = (
                    f"{p_guess} * {q_guess}"
                    + (" [PRIME]" if prime_status else " [NOT PRIME]")
                )

        print(
            f"{delta:8d} "
            f"{x_guess:10d} "
            f"{sign_of(k_res):>12} "
            f"{sign_of(r_res):>12} "
            f"{delta_guess:16d} "
            f"{str(square):>8} "
            f"{factors:>28}"
        )

    # --------------------------------------------------------
    # Exact residual values near X.
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("EXACT RESIDUAL VALUES NEAR X_TRUE")
    print("=" * 100)

    for delta in (-5, -2, -1, 0, 1, 2, 5):
        x_guess = X_TRUE + delta

        k_res = kappa_residual(x_guess)
        r_res = layer_residual(x_guess)

        print()
        print(f"delta      = {delta}")
        print(f"X_guess    = {x_guess}")
        print(f"kappa_diff = {k_res}")
        print(f"layer_diff = {compact_fraction(r_res)}")
        print(f"disc       = {discriminant(x_guess)}")

    # --------------------------------------------------------
    # Local monotonicity test.
    #
    # This is particularly important:
    # determine whether the sign of the residual consistently
    # identifies whether X_guess is below or above X_TRUE.
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("LOCAL MONOTONICITY TEST")
    print("=" * 100)

    left_kappa = True
    right_kappa = True
    left_layer = True
    right_layer = True

    for delta in range(1, 101):
        xl = X_TRUE - delta
        xr = X_TRUE + delta

        if not (kappa_residual(xl) < 0):
            left_kappa = False

        if not (kappa_residual(xr) > 0):
            right_kappa = False

        if not (layer_residual(xl) * layer_residual(X_TRUE - 1000) > 0):
            pass

        # For the layer invariant, we record the actual sign separately.
        if layer_residual(xl) == 0:
            left_layer = False

        if layer_residual(xr) == 0:
            right_layer = False

    print(f"kappa: X_guess < X_TRUE gives negative residual: {left_kappa}")
    print(f"kappa: X_guess > X_TRUE gives positive residual: {right_kappa}")

    print()
    print("Layer invariant local sign samples:")

    for delta in (1, 2, 5, 10, 25, 50, 100):
        left = layer_residual(X_TRUE - delta)
        right = layer_residual(X_TRUE + delta)

        print(
            f"  delta={delta:3d}  "
            f"left_sign={sign_of(left):>2}  "
            f"right_sign={sign_of(right):>2}"
        )

    # --------------------------------------------------------
    # Feasibility around the true value.
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("DISCRIMINANT FEASIBILITY")
    print("=" * 100)

    feasible = []

    for delta in range(-1000, 1001):
        x_guess = X_TRUE + delta

        square, root = square_discriminant(x_guess)

        if square:
            feasible.append(x_guess)

    print(f"Number of X values with square discriminant in +/-1000: {len(feasible)}")

    if feasible:
        print("Feasible X values:")
        print("  " + ", ".join(str(x) for x in feasible))

    print()
    print("=" * 100)
    print("END")
    print("=" * 100)


if __name__ == "__main__":
    main()
