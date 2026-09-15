#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 213 — EXACT 7-ADIC LOG/EXP INVERSE / PRINCIPAL-COORDINATE
                       UNIQUENESS AUDIT
==============================================================================

Experiment 212 established:

    rho = q1/q3,

    rho/6 = 1 + O(7^2),

    729 = 1 + O(7),

and independently verified

    m = log(rho/6) / log(729)

modulo 7^11,

where

    m = (k-1)/6.

Experiment 213 tests the inverse direction.

It computes finite 7-adic exponentials and verifies:

    exp(log(rho/6)) == rho/6,

and

    exp(m*log(729)) == rho/6,

to the tested 7-adic precision.

It also checks that if

    m' = m + 7^N,

then

    6*729^m' == 6*729^m mod 7^(N+1),

but that distinct residues modulo 7^N do not collapse.

This separates:

    logarithmic reconstruction,
    exponential reconstruction,
    uniqueness of the principal exponent coordinate.

No floating point.
No SymPy.
No giant powers of 3.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

P = 7
BASE = 2
MULT = 3

END_E = 12
PRECISION = END_E + 2


# ============================================================================
# BASIC 7-ADIC HELPERS
# ============================================================================

def valuation_7(x: int) -> int | None:
    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % P == 0:
        x //= P
        v += 1

    return v


def fraction_valuation_7(x: Fraction) -> int | None:
    if x == 0:
        return None

    vn = valuation_7(x.numerator)
    vd = valuation_7(x.denominator)

    vn = 0 if vn is None else vn
    vd = 0 if vd is None else vd

    return vn - vd


def inverse_mod(a: int, m: int) -> int:
    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r, r = a, m
    old_s, s = 1, 0

    while r:
        q = old_r // r

        old_r, r = (
            r,
            old_r - q * r,
        )

        old_s, s = (
            s,
            old_s - q * s,
        )

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def fraction_mod(
    x: Fraction,
    modulus: int,
) -> int:
    """
    Reduce a 7-adic integral rational modulo modulus.
    The denominator must be coprime to 7.
    """

    numerator = x.numerator
    denominator = x.denominator

    if denominator % P == 0:
        raise ArithmeticError(
            "fraction_mod received a denominator divisible by 7."
        )

    return (
        (numerator % modulus)
        * inverse_mod(
            denominator % modulus,
            modulus,
        )
    ) % modulus


def ratio_mod(modulus: int) -> int:
    return (
        Q1
        * inverse_mod(Q3, modulus)
    ) % modulus


# ============================================================================
# 7-ADIC LOGARITHM
# ============================================================================

def log_one_plus_x(
    x: Fraction,
    precision: int,
) -> Fraction:
    """
    Exact finite truncation of

        log(1+x)
          = sum_{n>=1} (-1)^(n+1) x^n/n.

    Terms are added until their 7-adic valuation is safely beyond
    the requested precision.
    """

    vx = fraction_valuation_7(x)

    if vx is None or vx <= 0:
        raise ArithmeticError(
            "log_one_plus_x requires v7(x) > 0."
        )

    total = Fraction(0)

    n = 1

    while True:

        term = (
            x ** n
        ) / n

        term_v = fraction_valuation_7(
            term
        )

        if term_v is None:
            break

        sign = (
            1
            if n % 2 == 1
            else -1
        )

        total += (
            sign * term
        )

        # Exact lower bound:
        #
        # v7(x^n/n) = n*v7(x) - v7(n).
        #
        # Once the valuation is safely beyond precision and future
        # valuations are monotone, the tail is invisible.
        if term_v >= precision:

            safe = True

            for j in range(
                n + 1,
                n + 8,
            ):

                vx_j = (
                    j * vx
                    - (
                        valuation_7(j)
                        or 0
                    )
                )

                if vx_j < precision:
                    safe = False
                    break

            if safe:
                break

        n += 1

        if n > 4 * precision + 50:
            raise ArithmeticError(
                "log truncation bound exceeded."
            )

    return total


# ============================================================================
# 7-ADIC EXPONENTIAL
# ============================================================================

def factorial_valuation_7(n: int) -> int:
    total = 0
    d = P

    while d <= n:
        total += n // d
        d *= P

    return total


def exp_7adic(
    y: Fraction,
    precision: int,
) -> Fraction:
    """
    Exact finite truncation of

        exp(y) = sum y^n/n!

    for v7(y) >= 1.

    The truncation rule uses

        v7(y^n/n!)
          = n*v7(y) - v7(n!).
    """

    vy = fraction_valuation_7(y)

    if vy is None:
        return Fraction(1)

    if vy <= 0:
        raise ArithmeticError(
            "exp_7adic requires v7(y) > 0."
        )

    total = Fraction(1)

    n = 1

    while True:

        term = (
            y ** n
        ) / factorial(n)

        term_v = fraction_valuation_7(
            term
        )

        if term_v is None:
            break

        total += term

        if term_v >= precision:

            safe = True

            for j in range(
                n + 1,
                n + 8,
            ):

                future_v = (
                    j * vy
                    - factorial_valuation_7(j)
                )

                if future_v < precision:
                    safe = False
                    break

            if safe:
                break

        n += 1

        if n > 8 * precision + 100:
            raise ArithmeticError(
                "exp truncation bound exceeded."
            )

    return total


def factorial(n: int) -> int:
    result = 1

    for i in range(2, n + 1):
        result *= i

    return result


# ============================================================================
# MODULAR PRINCIPAL EXPONENT
# ============================================================================

def orbit_mod(
    k: int,
    modulus: int,
) -> int:
    return (
        BASE
        * pow(
            MULT,
            k,
            modulus,
        )
    ) % modulus


def principal_orbit_mod(
    m: int,
    modulus: int,
) -> int:
    return (
        6
        * pow(
            729,
            m,
            modulus,
        )
    ) % modulus


def base_exponent_mod7() -> int:

    rho = ratio_mod(7)

    for k in range(6):

        if orbit_mod(k, 7) == rho:
            return k

    raise ArithmeticError(
        "No base exponent found."
    )


def defect_digit(
    e: int,
    m: int,
) -> int:

    modulus_e = P ** e
    modulus_next = P ** (
        e + 1
    )

    rho = ratio_mod(
        modulus_next
    )

    orbit = principal_orbit_mod(
        m,
        modulus_next,
    )

    diff = (
        orbit
        - rho
    ) % modulus_next

    if diff % modulus_e != 0:
        raise ArithmeticError(
            f"m={m} is not valid at level e={e}."
        )

    return (
        diff // modulus_e
    ) % P


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 213 — EXACT 7-ADIC LOG/EXP INVERSE / "
        "PRINCIPAL-COORDINATE UNIQUENESS AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    q1_v7 = valuation_7(Q1)
    q3_v7 = valuation_7(Q3)

    rho = Fraction(
        Q1,
        Q3,
    )

    rho_over_6 = (
        rho
        / 6
    )

    x_source = (
        rho_over_6
        - 1
    )

    x_base = (
        Fraction(729, 1)
        - 1
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  v7(q1)={q1_v7}"
    )

    print(
        f"  v7(q3)={q3_v7}"
    )

    print(
        f"  v7(rho/6-1)="
        f"{fraction_valuation_7(x_source)}"
    )

    print(
        f"  v7(729-1)="
        f"{fraction_valuation_7(x_base)}"
    )

    # ------------------------------------------------------------------
    # 2. INDEPENDENT LOGARITHMS
    # ------------------------------------------------------------------

    log_source = log_one_plus_x(
        x_source,
        PRECISION,
    )

    log_base = log_one_plus_x(
        x_base,
        PRECISION,
    )

    print()
    print("=" * 78)
    print("2. LOGARITHMIC DATA")
    print("=" * 78)

    log_source_v7 = fraction_valuation_7(
        log_source
    )

    log_base_v7 = fraction_valuation_7(
        log_base
    )

    print(
        f"  precision={PRECISION}"
    )

    print(
        f"  v7(log(rho/6))="
        f"{log_source_v7}"
    )

    print(
        f"  v7(log(729))="
        f"{log_base_v7}"
    )

    # ------------------------------------------------------------------
    # 3. HENSEL PRINCIPAL COORDINATE
    # ------------------------------------------------------------------

    k1 = base_exponent_mod7()

    current_m = (
        k1 - 1
    ) // 6

    m_sequence = [
        current_m
    ]

    print()
    print("=" * 78)
    print("3. INDEPENDENT HENSEL PRINCIPAL COORDINATE")
    print("=" * 78)

    print(
        f"  k1={k1}"
    )

    print(
        f"  m1={current_m}"
    )

    for e in range(
        1,
        END_E,
    ):

        delta = defect_digit(
            e,
            current_m,
        )

        t = (
            -delta
        ) % P

        current_m += (
            t
            * (P ** (e - 1))
        )

        m_sequence.append(
            current_m
        )

        print(
            f"  e={e}: "
            f"delta={delta} "
            f"t={t} "
            f"m_next={current_m}"
        )

    # ------------------------------------------------------------------
    # 4. LOGARITHMIC QUOTIENT
    # ------------------------------------------------------------------

    quotient = (
        log_source
        / log_base
    )

    quotient_v7 = fraction_valuation_7(
        quotient
    )

    print()
    print("=" * 78)
    print("4. LOGARITHMIC QUOTIENT")
    print("=" * 78)

    print(
        f"  v7(log_source/log_base)="
        f"{quotient_v7}"
    )

    # ------------------------------------------------------------------
    # 5. LOG/EXP RECONSTRUCTION
    # ------------------------------------------------------------------

    exp_log_source = exp_7adic(
        log_source,
        PRECISION,
    )

    reconstructed_ratio = exp_log_source

    print()
    print("=" * 78)
    print("5. EXP(LOG(rho/6)) RECONSTRUCTION")
    print("=" * 78)

    compare_modulus = P ** (
        END_E
    )

    rho_fraction_mod = fraction_mod(
        rho_over_6,
        compare_modulus,
    )

    exp_log_mod = fraction_mod(
        reconstructed_ratio,
        compare_modulus,
    )

    exp_log_exact = (
        rho_fraction_mod
        == exp_log_mod
    )

    print(
        f"  modulus={compare_modulus}"
    )

    print(
        f"  rho_over_6_mod={rho_fraction_mod}"
    )

    print(
        f"  exp_log_mod={exp_log_mod}"
    )

    print(
        f"  exact={exp_log_exact}"
    )

    # ------------------------------------------------------------------
    # 6. exp(m log 729)
    # ------------------------------------------------------------------

    # The final Hensel m is known independently.
    final_m = m_sequence[-1]

    m_log_base = (
        final_m
        * log_base
    )

    exp_m_log_base = exp_7adic(
        m_log_base,
        PRECISION,
    )

    exp_m_mod = fraction_mod(
        exp_m_log_base,
        compare_modulus,
    )

    exp_m_exact = (
        exp_m_mod
        == rho_fraction_mod
    )

    print()
    print("=" * 78)
    print("6. EXP(m * LOG(729)) RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  final_m={final_m}"
    )

    print(
        f"  rho_over_6_mod={rho_fraction_mod}"
    )

    print(
        f"  exp(m*log729)_mod={exp_m_mod}"
    )

    print(
        f"  exact={exp_m_exact}"
    )

    # ------------------------------------------------------------------
    # 7. LOGARITHMIC m COMPARISON
    # ------------------------------------------------------------------

    # Compare quotient modulo 7^(END_E-1).
    target_modulus = P ** (
        END_E - 1
    )

    quotient_mod = fraction_mod(
        quotient,
        target_modulus,
    )

    final_m_mod = (
        final_m
        % target_modulus
    )

    quotient_match = (
        quotient_mod
        == final_m_mod
    )

    print()
    print("=" * 78)
    print("7. LOGARITHMIC m COMPARISON")
    print("=" * 78)

    print(
        f"  modulus={target_modulus}"
    )

    print(
        f"  Hensel_m_mod={final_m_mod}"
    )

    print(
        f"  log_quotient_m_mod={quotient_mod}"
    )

    print(
        f"  exact={quotient_match}"
    )

    # ------------------------------------------------------------------
    # 8. UNIQUENESS TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. PRINCIPAL-COORDINATE UNIQUENESS")
    print("=" * 78)

    # Compare nearby principal coordinates modulo 7^(END_E-1).
    candidates = [
        final_m,
        final_m + target_modulus,
        final_m - target_modulus,
        final_m + P ** (END_E - 2),
        final_m - P ** (END_E - 2),
    ]

    uniqueness_rows = []

    for candidate in candidates:

        candidate_orbit = (
            principal_orbit_mod(
                candidate,
                compare_modulus,
            )
        )

        rho_mod = ratio_mod(
            compare_modulus
        )

        exact = (
            candidate_orbit
            == rho_mod
        )

        uniqueness_rows.append(
            {
                "candidate": candidate,
                "exact": exact,
            }
        )

        print(
            f"  candidate={candidate} "
            f"matches={exact}"
        )

    exact_candidates = [
        row["candidate"]
        for row in uniqueness_rows
        if row["exact"]
    ]

    # ------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiment established the closed logarithmic candidate

    m = log(rho/6) / log(729).

Experiment 213 checks the inverse map.

The two key identities are

    exp(log(rho/6)) = rho/6,

and

    exp(m log(729)) = rho/6.

If both hold to the tested precision, the logarithmic coordinate is not
just a numerical quotient: it is the inverse of the exponential
parameterization

    m -> 6*729^m.

The uniqueness test then checks whether nearby principal coordinates
produce the same residue at the tested modulus.

Together these give a much stronger basis for a theorem that the
canonical Hensel coordinate is the 7-adic logarithmic coordinate.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    uniqueness_ok = (
        exact_candidates == [final_m]
        or (
            final_m in exact_candidates
            and len(exact_candidates) == 1
        )
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and log_source_v7 == 2
        and log_base_v7 == 1
        and exp_log_exact
        and exp_m_exact
        and quotient_match
        and uniqueness_ok
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  log_source_v7_is_2="
        f"{log_source_v7 == 2}"
    )

    print(
        f"  log_base_v7_is_1="
        f"{log_base_v7 == 1}"
    )

    print(
        f"  exp_log_reconstructs_source="
        f"{exp_log_exact}"
    )

    print(
        f"  exp_m_log_reconstructs_source="
        f"{exp_m_exact}"
    )

    print(
        f"  log_quotient_matches_hensel_m="
        f"{quotient_match}"
    )

    print(
        f"  exact_candidates="
        f"{exact_candidates}"
    )

    print(
        f"  principal_coordinate_unique="
        f"{uniqueness_ok}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 213 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:
        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

