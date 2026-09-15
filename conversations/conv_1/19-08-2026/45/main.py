#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 212 — EXACT 7-ADIC LOGARITHMIC PRINCIPAL-COORDINATE AUDIT
==============================================================================

Goal
----

Experiment 211 explained the persistent slope-one Hensel law by writing

    k = 1 + 6m,

so that

    rho = q1/q3 = 6 * 729^m,

with

    729 = 1 + 728,
    v7(729-1) = 1.

The natural closed-form candidate is therefore

    m = log(rho/6) / log(729).

Experiment 212 tests this candidate p-adically.

Important details:

    rho/6 - 1 has v7 = 2 for the present data;
    729 - 1 has v7 = 1.

Hence

    v7(log(rho/6)) = 2,
    v7(log(729))   = 1,

and the quotient should have v7(m) = 1.

The computation is entirely finite:

    * truncated 7-adic logarithms are computed as rational numbers;
    * truncation stops once the remaining terms are beyond the tested
      7-adic precision;
    * the quotient is reduced modulo powers of 7;
    * the result is compared with the independently reconstructed
      principal coordinate m=(k-1)/6.

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

# Independently reconstructed through Experiment 211.
END_E = 12


# ============================================================================
# BASIC HELPERS
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


def ratio_mod(modulus: int) -> int:
    return (
        Q1
        * inverse_mod(Q3, modulus)
    ) % modulus


def orbit_mod(k: int, modulus: int) -> int:
    return (
        2
        * pow(3, k, modulus)
    ) % modulus


def order_3_mod_7e(e: int) -> int:
    return 6 * (P ** (e - 1))


def base_exponent_mod7() -> int:
    rho = ratio_mod(7)

    matches = []

    for k in range(6):

        if orbit_mod(k, 7) == rho:
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique base exponent, got {matches}."
        )

    return matches[0]


# ============================================================================
# MODULAR RATIONAL REDUCTION
# ============================================================================

def fraction_mod(frac: Fraction, modulus: int) -> int:
    """
    Reduce a rational number modulo modulus, provided its denominator is
    coprime to 7.

    This helper is used after extracting any 7-adic denominator valuation.
    """
    numerator = frac.numerator
    denominator = frac.denominator

    v_num = valuation_7(numerator)
    v_den = valuation_7(denominator)

    v_num = 0 if v_num is None else v_num
    v_den = 0 if v_den is None else v_den

    if v_den != 0:
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


def fraction_valuation_7(frac: Fraction) -> int | None:
    vn = valuation_7(frac.numerator)
    vd = valuation_7(frac.denominator)

    vn = 0 if vn is None else vn
    vd = 0 if vd is None else vd

    if frac.numerator == 0:
        return None

    return vn - vd


def normalized_unit_mod(
    frac: Fraction,
    modulus: int,
    target_valuation: int,
) -> int:
    """
    If v7(frac)=target_valuation, return

        frac / 7^target_valuation mod modulus.
    """

    if fraction_valuation_7(frac) != target_valuation:
        raise ArithmeticError(
            "Unexpected 7-adic valuation."
        )

    num = frac.numerator
    den = frac.denominator

    if target_valuation >= 0:
        num //= P ** target_valuation
    else:
        den *= P ** (-target_valuation)

    if den % P == 0:
        raise ArithmeticError(
            "Normalized denominator still divisible by 7."
        )

    return (
        (num % modulus)
        * inverse_mod(
            den % modulus,
            modulus,
        )
    ) % modulus


# ============================================================================
# 7-ADIC LOGARITHM
# ============================================================================

def log_one_plus_x(
    x: Fraction,
    precision: int,
) -> Fraction:
    """
    Compute the truncated 7-adic logarithm

        log(1+x)
          = x - x^2/2 + x^3/3 - ...

    sufficiently far that omitted terms have valuation >= precision.

    The function assumes v7(x) > 0.

    Since v7(x^n/n) = n*v7(x)-v7(n), the stopping rule is exact.
    """

    vx = fraction_valuation_7(x)

    if vx is None or vx <= 0:
        raise ArithmeticError(
            "log_one_plus_x requires v7(x) > 0."
        )

    total = Fraction(0)
    n = 1

    # Safe finite stopping bound. The valuation grows at least linearly.
    while n <= 4 * precision + 20:

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

        # Once terms are beyond precision and increasing thereafter,
        # the remaining tail is invisible at the requested precision.
        if (
            term_v >= precision
            and n * vx - valuation_7(n) >= precision
        ):
            # Check a few subsequent terms to make the monotonic
            # cutoff conservative.
            safe = True

            for j in range(
                n + 1,
                n + 6,
            ):
                if (
                    j * vx
                    - (
                        valuation_7(j)
                        or 0
                    )
                ) < precision:
                    safe = False
                    break

            if safe:
                break

        n += 1

    else:
        raise ArithmeticError(
            "Logarithm truncation bound was insufficient."
        )

    return total


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 212 — EXACT 7-ADIC LOGARITHMIC "
        "PRINCIPAL-COORDINATE AUDIT"
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

    # ------------------------------------------------------------------
    # 2. BASE EXPONENT
    # ------------------------------------------------------------------

    k1 = base_exponent_mod7()

    print()
    print("=" * 78)
    print("2. BASE EXPONENT CLASS")
    print("=" * 78)

    print(
        f"  rho_mod7={ratio_mod(7)}"
    )

    print(
        f"  k1={k1}"
    )

    # ------------------------------------------------------------------
    # 3. INDEPENDENT HENSEL RECONSTRUCTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. INDEPENDENT HENSEL RECONSTRUCTION")
    print("=" * 78)

    current_k = k1

    exponent_sequence = [
        k1
    ]

    transition_digits = []

    for e in range(
        1,
        END_E,
    ):

        modulus_e = P ** e
        modulus_next = P ** (
            e + 1
        )

        rho_next = ratio_mod(
            modulus_next
        )

        orbit_current = orbit_mod(
            current_k,
            modulus_next
        )

        diff = (
            orbit_current
            - rho_next
        ) % modulus_next

        if diff % modulus_e != 0:
            raise ArithmeticError(
                f"Current k does not solve level e={e}."
            )

        delta = (
            diff // modulus_e
        ) % P

        t = (
            -delta
        ) % P

        order = order_3_mod_7e(e)

        next_k = (
            current_k
            + t * order
        )

        if (
            orbit_mod(
                next_k,
                modulus_next,
            )
            != rho_next
        ):
            raise ArithmeticError(
                f"Lift failed at e={e}."
            )

        transition_digits.append(t)
        exponent_sequence.append(next_k)

        print(
            f"  e={e}: "
            f"k={current_k} "
            f"delta={delta} "
            f"t={t} "
            f"next_k={next_k}"
        )

        current_k = next_k

    # ------------------------------------------------------------------
    # 4. PRINCIPAL COORDINATE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PRINCIPAL COORDINATE m=(k-1)/6")
    print("=" * 78)

    m_sequence = []

    for E, k in enumerate(
        exponent_sequence,
        start=1,
    ):

        if (k - 1) % 6 != 0:
            raise ArithmeticError(
                f"k={k} is not 1 mod 6."
            )

        m = (
            k - 1
        ) // 6

        m_sequence.append(m)

        print(
            f"  E={E}: "
            f"k={k} "
            f"m={m}"
        )

    # ------------------------------------------------------------------
    # 5. LOG ARGUMENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT 7-ADIC LOG ARGUMENTS")
    print("=" * 78)

    # rho/6 - 1
    rho = Fraction(
        Q1,
        Q3,
    )

    x_source = (
        rho
        / 6
        - 1
    )

    x_base = Fraction(
        729,
        1,
    ) - 1

    print(
        f"  rho/6 - 1={x_source}"
    )

    print(
        f"  v7(rho/6-1)="
        f"{fraction_valuation_7(x_source)}"
    )

    print(
        f"  729-1={x_base}"
    )

    print(
        f"  v7(729-1)="
        f"{fraction_valuation_7(x_base)}"
    )

    # ------------------------------------------------------------------
    # 6. TRUNCATED 7-ADIC LOGS
    # ------------------------------------------------------------------

    # Need one extra digit because log(729) has valuation 1 and will
    # be divided by 7.
    log_precision = END_E + 2

    log_source = log_one_plus_x(
        x_source,
        log_precision,
    )

    log_base = log_one_plus_x(
        x_base,
        log_precision,
    )

    source_log_v7 = fraction_valuation_7(
        log_source
    )

    base_log_v7 = fraction_valuation_7(
        log_base
    )

    print()
    print("=" * 78)
    print("6. TRUNCATED 7-ADIC LOGARITHMS")
    print("=" * 78)

    print(
        f"  precision={log_precision}"
    )

    print(
        f"  v7(log(rho/6))={source_log_v7}"
    )

    print(
        f"  v7(log(729))={base_log_v7}"
    )

    print(
        f"  expected_source_v7=2"
    )

    print(
        f"  expected_base_v7=1"
    )

    # ------------------------------------------------------------------
    # 7. LOGARITHMIC QUOTIENT
    # ------------------------------------------------------------------

    if source_log_v7 is None or base_log_v7 is None:
        raise ArithmeticError(
            "Unexpected zero logarithm."
        )

    if base_log_v7 != 1:
        raise ArithmeticError(
            "Expected v7(log(729))=1."
        )

    if source_log_v7 != 2:
        raise ArithmeticError(
            "Expected v7(log(rho/6))=2."
        )

    # Normalize both logs.
    source_unit = (
        log_source
        / P ** 2
    )

    base_unit = (
        log_base
        / P
    )

    if fraction_valuation_7(
        source_unit
    ) != 0:
        raise ArithmeticError(
            "Source logarithm unit is not a 7-adic unit."
        )

    if fraction_valuation_7(
        base_unit
    ) != 0:
        raise ArithmeticError(
            "Base logarithm unit is not a 7-adic unit."
        )

    # The quotient is
    #
    #   log(rho/6)/log(729)
    #     = 7 * source_unit/base_unit.
    #
    # To compare modulo 7^(E-1), compare the unit quotient modulo
    # 7^(E-2).
    #
    # The largest required principal-coordinate precision is END_E-1.
    quotient_modulus = P ** (
        END_E - 2
    )

    source_unit_residue = fraction_mod(
        source_unit,
        quotient_modulus,
    )

    base_unit_residue = fraction_mod(
        base_unit,
        quotient_modulus,
    )

    unit_quotient = (
        source_unit_residue
        * inverse_mod(
            base_unit_residue,
            quotient_modulus,
        )
    ) % quotient_modulus

    logarithmic_m = (
        P
        * unit_quotient
    )

    print()
    print("=" * 78)
    print("7. LOGARITHMIC PRINCIPAL COORDINATE")
    print("=" * 78)

    print(
        f"  quotient_unit_mod={unit_quotient}"
    )

    print(
        f"  logarithmic_m_mod_7^{END_E - 1}="
        f"{logarithmic_m}"
    )

    # ------------------------------------------------------------------
    # 8. PRECISION-BY-PRECISION COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. PRECISION-BY-PRECISION LOG / HENSEL COMPARISON")
    print("=" * 78)

    log_match_flags = []

    for E, m in enumerate(
        m_sequence,
        start=1,
    ):

        if E == 1:
            modulus = 1
            log_residue = 0
        elif E == 2:
            modulus = P
            log_residue = (
                logarithmic_m
                % modulus
            )
        else:
            modulus = P ** (
                E - 1
            )
            log_residue = (
                logarithmic_m
                % modulus
            )

        exact = (
            m % modulus
            == log_residue
        )

        log_match_flags.append(
            exact
        )

        print(
            f"  E={E}: "
            f"modulus={modulus} "
            f"m={m % modulus} "
            f"log_m={log_residue} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 9. HIGH-PRECISION FINAL COMPARISON
    # ------------------------------------------------------------------

    final_modulus = P ** (
        END_E - 1
    )

    final_m = (
        m_sequence[-1]
        % final_modulus
    )

    log_m_final = (
        logarithmic_m
        % final_modulus
    )

    print()
    print("=" * 78)
    print("9. FINAL LOGARITHMIC COMPARISON")
    print("=" * 78)

    print(
        f"  modulus={final_modulus}"
    )

    print(
        f"  Hensel_m={final_m}"
    )

    print(
        f"  logarithmic_m={log_m_final}"
    )

    print(
        f"  exact={final_m == log_m_final}"
    )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The principal-coordinate identity is

    rho = 6 * 729^m.

Since

    rho/6 = 1 + O(7^2),

and

    729 = 1 + O(7),

the 7-adic logarithm is defined and gives formally

    log(rho/6)
       = m * log(729).

Thus

    m
      = log(rho/6) / log(729).

The experiment computes both sides with finite 7-adic precision.

The important point is that this is an independent reconstruction:
the logarithmic quotient is not built from the Hensel digits.

Agreement would convert the recursive Hensel description into a closed
7-adic logarithmic formula for the principal exponent coordinate.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    all_log_matches = all(
        log_match_flags
    )

    final_exact = (
        final_m
        == log_m_final
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and source_log_v7 == 2
        and base_log_v7 == 1
        and all_log_matches
        and final_exact
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_log_v7_is_2="
        f"{source_log_v7 == 2}"
    )

    print(
        f"  base_log_v7_is_1="
        f"{base_log_v7 == 1}"
    )

    print(
        f"  all_precision_log_matches="
        f"{all_log_matches}"
    )

    print(
        f"  final_log_match="
        f"{final_exact}"
    )

    print(
        f"  Hensel_m={final_m}"
    )

    print(
        f"  logarithmic_m={log_m_final}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 212 COMPLETE")


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

