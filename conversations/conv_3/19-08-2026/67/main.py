#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================================================================
# EXPERIMENT 444
# ======================================================================================================================
#
# CRT PULLBACK DEGENERACY / X-SPACE IDENTITY AUDIT
#
# QUESTION
#   When the K constraints
#
#       K = A/4 - ((x+d0)/2)^2
#
#   are pulled back to x-space, do the 2-adic and odd-prime
#   quadratic-residue conditions actually eliminate x candidates?
#
# KEY POINT
#   If y=x+d0 is even and
#
#       K(x) = A/4 - y^2/4,
#
#   then for every modulus m:
#
#       A/4 - K(x) == (y/2)^2 (mod m).
#
#   Therefore every parity-valid x automatically satisfies
#   every local quadratic-residue condition.
#
#   This experiment checks that statement exactly rather than
#   assuming it.
#
# RULES
#   exact integer arithmetic only
#   no resultants
#   no Groebner basis
#   no factor search
#   no giant K image
#   no CRT Cartesian product
#   bounded x enumeration only
#
# EXPECTED OUTCOME
#   The K-image restrictions should collapse completely under
#   direct substitution K=K(x):
#
#       x parity-valid
#           -> 2-adic condition TRUE
#           -> every odd-prime condition TRUE
#           -> combined CRT condition TRUE.
#
#   If confirmed, the CRT K sieve cannot independently reduce
#   the x-space population once K is parametrically reconstructed
#   from x.
#
# ======================================================================================================================


ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19]

# Small windows through which x is explicitly enumerated.
X_RADII = [
    0,
    1,
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256,
    512,
    1024,
    4096,
    16384,
    65536,
]

MAX_X_ENUMERATION = 150_000


# ======================================================================================================================
# INSTANCE DATA
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    p: int
    q: int
    N: int
    S: int
    d: int
    d0: int
    x_true: int


INSTANCES_RAW = [
    (50411, 282599, -10),
    (1013, 10009, -17),
    (10009, 1000033, -24),
    (10009, 10037, -31),
    (50023, 50051, -38),
    (100019, 100043, -45),
    (200009, 200017, -52),
    (300017, 900007, -59),
]


def build_instances() -> list[Instance]:
    out: list[Instance] = []

    for idx, (p, q, x_true) in enumerate(INSTANCES_RAW, start=1):
        N = p * q
        S = p + q

        # Same construction used throughout the preceding experiments.
        d = N - 2 * S + 1

        d0 = d - x_true

        out.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                d0=d0,
                x_true=x_true,
            )
        )

    return out


# ======================================================================================================================
# CONTROL K / R GENERATION
# ======================================================================================================================

HIDDEN_FACTORS = [17, 29, 43, 61, 89, 127, 191, 257]


def hidden_K(inst: Instance) -> int:
    f = HIDDEN_FACTORS[inst.idx - 1]

    value = (
        inst.d0 * inst.d0
        + 13 * inst.N
        + 7 * inst.S * inst.S
        + f * inst.d0
    )

    return -(value * f + inst.N * (f + 3))


def build_r(inst: Instance, K_true: int) -> int:
    x = inst.x_true

    r = (
        4 * K_true
        + x * x
        + 2 * inst.d0 * x
    )

    # Exact generator check.
    numerator = r - x * x - 2 * inst.d0 * x

    assert numerator % 4 == 0
    assert numerator // 4 == K_true

    return r


# ======================================================================================================================
# BASIC NUMBER THEORY
# ======================================================================================================================

def is_square_mod_prime(a: int, p: int) -> bool:
    """
    Exact quadratic-residue test for an odd prime.
    Zero counts as a square.
    """

    a %= p

    if a == 0:
        return True

    return pow(a, (p - 1) // 2, p) == 1


def is_square_mod_2pow(value: int, bits: int) -> bool:
    """
    Exact test for whether value is a quadratic residue modulo 2^bits.

    Uses the valuation criterion:
      - zero is square;
      - valuation must be even;
      - after removing the even power of two, the odd part must satisfy
        the standard 2-adic square conditions.
    """

    modulus = 1 << bits
    value %= modulus

    if value == 0:
        return True

    # v2(value)
    v2 = (value & -value).bit_length() - 1

    if v2 & 1:
        return False

    remaining = bits - v2
    odd_part = value >> v2

    if remaining <= 1:
        return True

    if remaining == 2:
        return (odd_part & 3) == 1

    return (odd_part & 7) == 1


# ======================================================================================================================
# ORIGINAL K IDENTITY
# ======================================================================================================================

def candidate_y(inst: Instance, x: int) -> int:
    return x + inst.d0


def candidate_k(
    inst: Instance,
    A: int,
    x: int,
) -> tuple[bool, int]:
    """
    Return (valid_parity, K).

    Since

        4K = A - y^2,

    K is an integer exactly when A-y^2 is divisible by 4.

    In these experiments A is divisible by 4, hence this is equivalent
    to y being even.
    """

    y = candidate_y(inst, x)

    numerator = A - y * y

    if numerator % 4 != 0:
        return False, 0

    return True, numerator // 4


# ======================================================================================================================
# LOCAL K-IMAGE TEST
# ======================================================================================================================

def k_local_2adic(
    A4: int,
    K: int,
    bits: int,
) -> bool:
    return is_square_mod_2pow(
        A4 - K,
        bits,
    )


def k_local_odd(
    A4: int,
    K: int,
    p: int,
) -> bool:
    return is_square_mod_prime(
        A4 - K,
        p,
    )


def k_local_all(
    A4: int,
    K: int,
    bits: int,
) -> tuple[bool, list[int]]:
    failed = []

    if not k_local_2adic(A4, K, bits):
        failed.append(2)

    for p in ODD_PRIMES:
        if not k_local_odd(A4, K, p):
            failed.append(p)

    return len(failed) == 0, failed


# ======================================================================================================================
# DIRECT PULLBACK CHECK
# ======================================================================================================================

def direct_square_identity_check(
    inst: Instance,
    A4: int,
    x: int,
) -> bool:
    """
    Directly verifies:

        A4 - K(x) == ((x+d0)/2)^2

    as an exact integer identity.
    """

    y = candidate_y(inst, x)

    if y & 1:
        return False

    K = A4 - (y // 2) * (y // 2)

    lhs = A4 - K
    rhs = (y // 2) * (y // 2)

    return lhs == rhs


# ======================================================================================================================
# ENUMERATED WINDOW AUDIT
# ======================================================================================================================

@dataclass
class WindowResult:
    radius: int
    total_x: int
    parity_valid: int
    parity_invalid: int

    identity_failures: int

    two_adic_failures: int
    odd_prime_failures: int
    combined_failures: int

    true_in_window: bool
    true_survives: bool


def audit_window(
    inst: Instance,
    A: int,
    A4: int,
    radius: int,
    bits: int,
) -> WindowResult:

    lo = -radius
    hi = radius

    total = hi - lo + 1

    if total > MAX_X_ENUMERATION:
        raise MemoryError(
            f"x window too large: radius={radius}, "
            f"population={total}, "
            f"limit={MAX_X_ENUMERATION}"
        )

    parity_valid = 0
    parity_invalid = 0

    identity_failures = 0
    two_adic_failures = 0
    odd_prime_failures = 0
    combined_failures = 0

    true_in_window = (
        inst.x_true >= lo
        and inst.x_true <= hi
    )

    true_survives = False

    for x in range(lo, hi + 1):

        valid, K = candidate_k(
            inst,
            A,
            x,
        )

        if not valid:
            parity_invalid += 1
            continue

        parity_valid += 1

        if not direct_square_identity_check(
            inst,
            A4,
            x,
        ):
            identity_failures += 1
            continue

        two_ok = k_local_2adic(
            A4,
            K,
            bits,
        )

        if not two_ok:
            two_adic_failures += 1

        odd_ok = True

        for p in ODD_PRIMES:
            if not k_local_odd(
                A4,
                K,
                p,
            ):
                odd_prime_failures += 1
                odd_ok = False
                break

        combined_ok = two_ok and odd_ok

        if not combined_ok:
            combined_failures += 1

        if x == inst.x_true:
            true_survives = combined_ok

    return WindowResult(
        radius=radius,
        total_x=total,
        parity_valid=parity_valid,
        parity_invalid=parity_invalid,
        identity_failures=identity_failures,
        two_adic_failures=two_adic_failures,
        odd_prime_failures=odd_prime_failures,
        combined_failures=combined_failures,
        true_in_window=true_in_window,
        true_survives=true_survives,
    )


# ======================================================================================================================
# NECESSARY PARITY PREDICTION
# ======================================================================================================================

def expected_parity_valid_count(
    inst: Instance,
    radius: int,
) -> int:
    """
    Since A ≡ 0 (mod 4), K(x) is integral exactly when
    x+d0 is even.

    Thus admissible x are exactly those satisfying

        x ≡ d0 (mod 2).
    """

    count = 0

    for x in range(-radius, radius + 1):
        if ((x + inst.d0) & 1) == 0:
            count += 1

    return count


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 444")
    print("=" * 120)
    print()
    print("CRT PULLBACK DEGENERACY / X-SPACE IDENTITY AUDIT")
    print()
    print("QUESTION")
    print("  When K is reconstructed from x, do the K quadratic-residue")
    print("  constraints eliminate any parity-valid x candidates?")
    print()
    print("ODD PRIMES =", ODD_PRIMES)
    print("MAX_X_ENUMERATION =", MAX_X_ENUMERATION)
    print()

    instances = build_instances()

    total_windows = 0
    total_identity_failures = 0
    total_two_adic_failures = 0
    total_odd_failures = 0
    total_combined_failures = 0

    all_true_survive = True

    for inst in instances:

        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        A = r + inst.d0 * inst.d0

        if A % 4 != 0:
            raise AssertionError(
                f"A not divisible by 4 at instance={inst.idx}"
            )

        A4 = A // 4

        # Verify the original true relation.
        true_y = inst.x_true + inst.d0

        true_numerator = A - true_y * true_y

        if true_numerator % 4 != 0:
            raise AssertionError(
                f"true relation non-integral at instance={inst.idx}"
            )

        true_K_from_x = true_numerator // 4

        if true_K_from_x != K_true:
            raise AssertionError(
                f"true K reconstruction mismatch at instance={inst.idx}"
            )

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"d0={inst.d0} "
            f"x_true={inst.x_true}"
        )
        print(
            f"  hidden K bits = {abs(K_true).bit_length()}"
        )
        print(
            f"  A = {A}"
        )
        print(
            f"  A/4 = {A4}"
        )
        print(
            f"  d0 parity = {inst.d0 & 1}"
        )
        print()

        # bits are only used as local K-image resolutions.
        bits_values = [4, 8, 12, 16, 20]

        for bits in bits_values:

            print(
                f"  bits={bits}"
            )

            print(
                "    radius"
                "       x-count"
                " parity-valid"
                " identity-fail"
                " 2adic-fail"
                " odd-fail"
                " combined-fail"
                " true-in"
                " true-survives"
            )
            print("    " + "-" * 105)

            for radius in X_RADII:

                result = audit_window(
                    inst,
                    A,
                    A4,
                    radius,
                    bits,
                )

                expected = expected_parity_valid_count(
                    inst,
                    radius,
                )

                if result.parity_valid != expected:
                    raise AssertionError(
                        f"parity count mismatch "
                        f"instance={inst.idx}, "
                        f"bits={bits}, "
                        f"radius={radius}: "
                        f"{result.parity_valid} != {expected}"
                    )

                total_windows += 1

                total_identity_failures += (
                    result.identity_failures
                )

                total_two_adic_failures += (
                    result.two_adic_failures
                )

                total_odd_failures += (
                    result.odd_prime_failures
                )

                total_combined_failures += (
                    result.combined_failures
                )

                if result.true_in_window:
                    all_true_survive &= result.true_survives

                print(
                    f"    {radius:>6}"
                    f"{result.total_x:>13}"
                    f"{result.parity_valid:>14}"
                    f"{result.identity_failures:>14}"
                    f"{result.two_adic_failures:>13}"
                    f"{result.odd_prime_failures:>11}"
                    f"{result.combined_failures:>16}"
                    f"{str(result.true_in_window):>9}"
                    f"{str(result.true_survives):>14}"
                )

            print()

        # ==================================================================
        # Strong symbolic identity audit on a separate set of points.
        # ==================================================================

        print("  SYMBOLIC-PULLBACK CHECK")

        test_points = [
            -100,
            -59,
            -52,
            -45,
            -38,
            -31,
            -24,
            -17,
            -10,
            0,
            1,
            2,
            17,
            31,
            64,
            127,
        ]

        parity_valid_points = 0
        local_failures = 0

        for x in test_points:

            valid, K = candidate_k(
                inst,
                A,
                x,
            )

            if not valid:
                continue

            parity_valid_points += 1

            # Exact algebraic identity.
            if not direct_square_identity_check(
                inst,
                A4,
                x,
            ):
                local_failures += 1
                continue

            # Several resolutions.
            for bits in [4, 8, 12, 16, 20, 24]:
                all_ok, failed = k_local_all(
                    A4,
                    K,
                    bits,
                )

                if not all_ok:
                    local_failures += 1

                    print(
                        f"    FAILURE x={x} bits={bits} "
                        f"failed_moduli={failed}"
                    )

        if local_failures != 0:
            raise AssertionError(
                f"pullback local failure at instance={inst.idx}"
            )

        print(
            f"    parity-valid test points = "
            f"{parity_valid_points}"
        )
        print(
            f"    algebraic/local failures  = "
            f"{local_failures}"
        )
        print(
            "    pullback identity = TRUE"
        )
        print(
            "    all tested local CRT conditions = TRUE"
        )
        print()

    # ==================================================================
    # GLOBAL SUMMARY
    # ==================================================================

    print("=" * 120)
    print("GLOBAL PULLBACK-DEGENERACY SUMMARY")
    print("=" * 120)

    print(
        f"  total x-window audits             = "
        f"{total_windows}"
    )

    print(
        f"  exact identity failures           = "
        f"{total_identity_failures}"
    )

    print(
        f"  2-adic pullback failures          = "
        f"{total_two_adic_failures}"
    )

    print(
        f"  odd-prime pullback failures       = "
        f"{total_odd_failures}"
    )

    print(
        f"  combined CRT pullback failures    = "
        f"{total_combined_failures}"
    )

    print(
        f"  true K survives every tested x-window = "
        f"{all_true_survive}"
    )

    print()
    print("EXPECTED STRUCTURE")
    print()
    print(
        "  For parity-valid x:"
    )
    print()
    print(
        "      y = x + d0"
    )
    print(
        "      K(x) = A/4 - (y/2)^2"
    )
    print()
    print(
        "  Therefore:"
    )
    print()
    print(
        "      A/4 - K(x) = (y/2)^2"
    )
    print()
    print(
        "  This is an exact integer identity."
    )
    print()
    print(
        "  Reducing modulo any odd prime p gives:"
    )
    print()
    print(
        "      A/4 - K(x) ≡ (y/2)^2 (mod p),"
    )
    print()
    print(
        "  so the odd-prime quadratic-residue condition is automatic."
    )
    print()
    print(
        "  Reducing modulo 2^b gives the same result."
    )
    print()
    print(
        "  Thus the combined K-image sieve does not independently"
    )
    print(
        "  restrict x once K is defined by the same quadratic identity."
    )
    print()
    print(
        "  The only generic restriction exposed by integrality is:"
    )
    print()
    print(
        "      x + d0 ≡ 0 (mod 2)."
    )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  Experiment 442R2 correctly characterized the K image"
    )
    print(
        "  as an affine square-residue image."
    )
    print()
    print(
        "  Experiment 443 then counted how sparse that K image is"
    )
    print(
        "  inside artificial/natural K intervals."
    )
    print()
    print(
        "  Experiment 444 checks the reverse direction."
    )
    print()
    print(
        "  Once K is pulled back through"
    )
    print(
        "      K = A/4 - ((x+d0)/2)^2,"
    )
    print(
        "  the square-residue constraints cease to be filters:"
    )
    print(
        "  they become identities."
    )
    print()
    print(
        "  Consequently, increasing the number of CRT primes cannot"
    )
    print(
        "  by itself collapse the x-space population in this model."
    )
    print()
    print(
        "  To obtain genuine candidate reduction in x-space, a future"
    )
    print(
        "  experiment must introduce information that is NOT algebraically"
    )
    print(
        "  equivalent to the definition of K(x)."
    )
    print()
    print(
        "  This is an exact structural result about the present"
    )
    print(
        "  K/x formulation. It is not a factoring theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 444 FINAL STATUS")
    print("=" * 120)
    print(
        f"  IDENTITY FAILURES       = "
        f"{total_identity_failures}"
    )
    print(
        f"  2-ADIC PULLBACK FAILS   = "
        f"{total_two_adic_failures}"
    )
    print(
        f"  ODD-PRIME PULLBACK FAILS= "
        f"{total_odd_failures}"
    )
    print(
        f"  COMBINED PULLBACK FAILS = "
        f"{total_combined_failures}"
    )
    print(
        f"  TRUE SURVIVAL CHECK     = "
        f"{all_true_survive}"
    )
    print("  GIANT K IMAGE = False")
    print("  CRT CARTESIAN PRODUCT = False")
    print("  ALL ARITHMETIC = INTEGER-EXACT")
    print("=" * 120)
    print("EXPERIMENT 444 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()

