#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 214RR — EXACT LOGARITHMIC-COORDINATE UNIQUENESS MOD 7^N AUDIT
==============================================================================

Corrected Experiment 214R.

The previous script defined q1_v7/q3_v7 only implicitly for printing but
then referenced them in the final check without storing them.

This version computes those values explicitly at the beginning of main()
and reuses them consistently.

It also keeps all modular exponents strictly integral.

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
PRINCIPAL_SCALE = 6
PRINCIPAL_BASE = 729

END_E = 12
PRECISION = END_E + 2


# ============================================================================
# BASIC HELPERS
# ============================================================================

def valuation_7(x: int) -> int | None:
    if x == 0:
        return None

    x = abs(int(x))
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
    m = int(m)
    a = int(a) % m

    if m <= 1:
        return 0

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

    return int(old_s % m)


def fraction_mod(
    x: Fraction,
    modulus: int,
) -> int:

    modulus = int(modulus)

    if x.denominator % P == 0:
        raise ArithmeticError(
            "Denominator divisible by 7."
        )

    numerator = int(x.numerator)
    denominator = int(x.denominator)

    return int(
        (
            (numerator % modulus)
            * inverse_mod(
                denominator % modulus,
                modulus,
            )
        )
        % modulus
    )


def ratio_mod(modulus: int) -> int:
    modulus = int(modulus)

    return int(
        (
            Q1
            * inverse_mod(
                Q3,
                modulus,
            )
        )
        % modulus
    )


def principal_orbit_mod(
    m: int,
    modulus: int,
) -> int:

    m = int(m)
    modulus = int(modulus)

    return int(
        (
            PRINCIPAL_SCALE
            * pow(
                PRINCIPAL_BASE,
                m,
                modulus,
            )
        )
        % modulus
    )


# ============================================================================
# 7-ADIC LOG / EXP
# ============================================================================

def factorial(n: int) -> int:
    result = 1

    for i in range(2, int(n) + 1):
        result *= i

    return result


def factorial_valuation_7(n: int) -> int:
    n = int(n)

    total = 0
    d = P

    while d <= n:
        total += n // d
        d *= P

    return total


def log_one_plus_x(
    x: Fraction,
    precision: int,
) -> Fraction:

    vx = fraction_valuation_7(x)

    if vx is None or vx <= 0:
        raise ArithmeticError(
            "log requires v7(x)>0."
        )

    total = Fraction(0)
    n = 1

    while True:

        term = (
            x ** n
        ) / n

        tv = fraction_valuation_7(term)

        sign = (
            1
            if n % 2 == 1
            else -1
        )

        total += (
            sign * term
        )

        if tv is not None and tv >= precision:

            safe = True

            for j in range(
                n + 1,
                n + 8,
            ):

                future_v = (
                    j * vx
                    - (
                        valuation_7(j)
                        or 0
                    )
                )

                if future_v < precision:
                    safe = False
                    break

            if safe:
                break

        n += 1

        if n > 4 * precision + 100:
            raise ArithmeticError(
                "log truncation failed."
            )

    return total


def exp_7adic(
    y: Fraction,
    precision: int,
) -> Fraction:

    vy = fraction_valuation_7(y)

    if vy is None:
        return Fraction(1)

    if vy <= 0:
        raise ArithmeticError(
            "exp requires v7(y)>0."
        )

    total = Fraction(1)
    n = 1

    while True:

        term = (
            y ** n
        ) / factorial(n)

        tv = fraction_valuation_7(term)

        total += term

        if tv is not None and tv >= precision:

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
                "exp truncation failed."
            )

    return total


# ============================================================================
# HENSEL PRINCIPAL COORDINATE
# ============================================================================

def base_exponent_mod7() -> int:

    rho = ratio_mod(7)
    matches = []

    for k in range(6):

        k = int(k)

        value = int(
            (
                2
                * pow(
                    3,
                    k,
                    7,
                )
            )
            % 7
        )

        if value == rho:
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique base exponent, got {matches}."
        )

    return int(matches[0])


def defect_digit(
    e: int,
    m: int,
) -> int:

    e = int(e)
    m = int(m)

    modulus_e = int(
        P ** e
    )

    modulus_next = int(
        P ** (e + 1)
    )

    rho_next = int(
        ratio_mod(
            modulus_next
        )
    )

    orbit_next = int(
        principal_orbit_mod(
            m,
            modulus_next,
        )
    )

    diff = int(
        (
            orbit_next
            - rho_next
        )
        % modulus_next
    )

    if diff % modulus_e != 0:
        raise ArithmeticError(
            f"m={m} is not valid at level e={e}."
        )

    return int(
        (
            diff // modulus_e
        )
        % P
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 214RR — EXACT LOGARITHMIC-COORDINATE "
        "UNIQUENESS MOD 7^N AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

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
        Fraction(
            PRINCIPAL_BASE,
            1,
        )
        - 1
    )

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

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
    # 2. LOGARITHMIC COORDINATE
    # ------------------------------------------------------------------

    log_source = log_one_plus_x(
        x_source,
        PRECISION,
    )

    log_base = log_one_plus_x(
        x_base,
        PRECISION,
    )

    logarithmic_m = (
        log_source
        / log_base
    )

    log_source_v7 = fraction_valuation_7(
        log_source
    )

    log_base_v7 = fraction_valuation_7(
        log_base
    )

    print()
    print("=" * 78)
    print("2. LOGARITHMIC PRINCIPAL COORDINATE")
    print("=" * 78)

    print(
        f"  v7(log_source)={log_source_v7}"
    )

    print(
        f"  v7(log_base)={log_base_v7}"
    )

    print(
        f"  v7(logarithmic_m)="
        f"{fraction_valuation_7(logarithmic_m)}"
    )

    # ------------------------------------------------------------------
    # 3. HENSEL PRINCIPAL COORDINATE
    # ------------------------------------------------------------------

    k1 = int(
        base_exponent_mod7()
    )

    current_m = int(
        (k1 - 1) // 6
    )

    m_sequence = [
        int(current_m)
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

        e = int(e)

        delta = int(
            defect_digit(
                e,
                current_m,
            )
        )

        t = int(
            (-delta) % P
        )

        increment = int(
            t
            * (
                P ** (
                    e - 1
                )
            )
        )

        current_m = int(
            current_m
            + increment
        )

        m_sequence.append(
            int(current_m)
        )

        print(
            f"  e={e}: "
            f"delta={delta} "
            f"t={t} "
            f"increment={increment} "
            f"m_next={current_m}"
        )

    if not all(
        type(m) is int
        for m in m_sequence
    ):
        raise TypeError(
            "Hensel m_sequence contains a non-int value."
        )

    # ------------------------------------------------------------------
    # 4. PRECISION CHECKS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. PRECISION-BY-PRECISION "
        "LOG/HENSEL/UNIQUENESS"
    )
    print("=" * 78)

    precision_flags = []

    for E in range(
        1,
        END_E + 1,
    ):

        E = int(E)

        orbit_modulus = int(
            P ** E
        )

        principal_modulus = int(
            P ** (
                E - 1
            )
        )

        if principal_modulus > 1:
            m_hensel = int(
                m_sequence[E - 1]
                % principal_modulus
            )

            m_log = int(
                fraction_mod(
                    logarithmic_m,
                    principal_modulus,
                )
            )
        else:
            m_hensel = 0
            m_log = 0

        log_match = (
            m_hensel
            == m_log
        )

        orbit_match = (
            principal_orbit_mod(
                int(m_hensel),
                int(orbit_modulus),
            )
            ==
            ratio_mod(
                int(orbit_modulus),
            )
        )

        separated = True

        if E >= 2:

            smaller_step = int(
                P ** (
                    E - 2
                )
            )

            plus_candidate = int(
                m_hensel
                + smaller_step
            )

            minus_candidate = int(
                m_hensel
                - smaller_step
            )

            target = int(
                ratio_mod(
                    orbit_modulus
                )
            )

            plus_match = (
                principal_orbit_mod(
                    plus_candidate,
                    orbit_modulus,
                )
                == target
            )

            minus_match = (
                principal_orbit_mod(
                    minus_candidate,
                    orbit_modulus,
                )
                == target
            )

            separated = (
                not plus_match
                and not minus_match
            )

        precision_ok = (
            log_match
            and orbit_match
            and separated
        )

        precision_flags.append(
            precision_ok
        )

        print(
            f"  E={E}: "
            f"modulus={orbit_modulus} "
            f"m_hensel={m_hensel} "
            f"m_log={m_log} "
            f"log_match={log_match} "
            f"orbit_match={orbit_match} "
            f"separated={separated}"
        )

    # ------------------------------------------------------------------
    # 5. FINAL CLASS REPRESENTATIVES
    # ------------------------------------------------------------------

    principal_modulus = int(
        P ** (
            END_E - 1
        )
    )

    final_hensel = int(
        m_sequence[-1]
        % principal_modulus
    )

    final_log = int(
        fraction_mod(
            logarithmic_m,
            principal_modulus,
        )
    )

    orbit_modulus = int(
        P ** END_E
    )

    target = int(
        ratio_mod(
            orbit_modulus
        )
    )

    class_representatives = [
        int(final_hensel),
        int(
            final_hensel
            + principal_modulus
        ),
        int(
            final_hensel
            - principal_modulus
        ),
    ]

    class_match_flags = []

    for candidate in class_representatives:

        candidate = int(candidate)

        match = (
            principal_orbit_mod(
                candidate,
                orbit_modulus,
            )
            == target
        )

        class_match_flags.append(
            match
        )

    smallest_distinguishing_step = int(
        P ** (
            END_E - 2
        )
    )

    neighbor_candidates = [
        int(
            final_hensel
            + smallest_distinguishing_step
        ),
        int(
            final_hensel
            - smallest_distinguishing_step
        ),
    ]

    neighbor_flags = []

    for candidate in neighbor_candidates:

        candidate = int(candidate)

        match = (
            principal_orbit_mod(
                candidate,
                orbit_modulus,
            )
            == target
        )

        neighbor_flags.append(
            match
        )

    print()
    print("=" * 78)
    print("5. FINAL CLASS UNIQUENESS")
    print("=" * 78)

    print(
        f"  principal_modulus={principal_modulus}"
    )

    print(
        f"  Hensel_residue={final_hensel}"
    )

    print(
        f"  log_residue={final_log}"
    )

    print(
        f"  Hensel_equals_log="
        f"{final_hensel == final_log}"
    )

    print(
        f"  class_representatives="
        f"{class_representatives}"
    )

    print(
        f"  class_match_flags="
        f"{class_match_flags}"
    )

    print(
        f"  distinguishing_step="
        f"{smallest_distinguishing_step}"
    )

    print(
        f"  neighbor_candidates="
        f"{neighbor_candidates}"
    )

    print(
        f"  neighbor_match_flags="
        f"{neighbor_flags}"
    )

    # ------------------------------------------------------------------
    # 6. EXP/LOG CHECK
    # ------------------------------------------------------------------

    exp_log_source = exp_7adic(
        log_source,
        PRECISION,
    )

    final_m_int = int(
        m_sequence[-1]
    )

    m_log_base = (
        final_m_int
        * log_base
    )

    exp_m_log = exp_7adic(
        m_log_base,
        PRECISION,
    )

    compare_modulus = int(
        P ** END_E
    )

    source_mod = int(
        fraction_mod(
            rho_over_6,
            compare_modulus,
        )
    )

    exp_log_mod = int(
        fraction_mod(
            exp_log_source,
            compare_modulus,
        )
    )

    exp_m_mod = int(
        fraction_mod(
            exp_m_log,
            compare_modulus,
        )
    )

    exp_log_ok = (
        source_mod
        == exp_log_mod
    )

    exp_m_ok = (
        source_mod
        == exp_m_mod
    )

    print()
    print("=" * 78)
    print("6. EXP/LOG INVERSE CHECKS")
    print("=" * 78)

    print(
        f"  modulus={compare_modulus}"
    )

    print(
        f"  rho_over_6={source_mod}"
    )

    print(
        f"  exp_log={exp_log_mod}"
    )

    print(
        f"  exp_log_exact={exp_log_ok}"
    )

    print(
        f"  exp_m_log={exp_m_mod}"
    )

    print(
        f"  exp_m_log_exact={exp_m_ok}"
    )

    # ------------------------------------------------------------------
    # 7. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The uniqueness statement is about a residue class modulo 7^(E-1).

Thus

    m,
    m + 7^(E-1),
    m - 7^(E-1)

must represent the same principal-coordinate class and therefore give
the same orbit residue modulo 7^E.

In contrast,

    m +/- 7^(E-2)

must change the orbit at the tested precision.

The logarithmic quotient gives the same class:

    m = log(rho/6) / log(729).

The exp/log tests then verify the inverse parameterization

    m -> 6*729^m.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    all_precision = all(
        precision_flags
    )

    same_class_exact = all(
        class_match_flags
    )

    neighbors_separated = all(
        not flag
        for flag in neighbor_flags
    )

    log_hensel_exact = (
        final_hensel
        == final_log
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and all_precision
        and same_class_exact
        and neighbors_separated
        and exp_log_ok
        and exp_m_ok
        and log_hensel_exact
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_v7_is_0={q1_v7 == 0}"
    )

    print(
        f"  q3_v7_is_0={q3_v7 == 0}"
    )

    print(
        f"  all_precision_checks_exact="
        f"{all_precision}"
    )

    print(
        f"  final_log_hensel_class_exact="
        f"{log_hensel_exact}"
    )

    print(
        f"  all_same_class_representatives_match="
        f"{same_class_exact}"
    )

    print(
        f"  neighboring_classes_separated="
        f"{neighbors_separated}"
    )

    print(
        f"  exp_log_inverse_exact="
        f"{exp_log_ok}"
    )

    print(
        f"  exp_m_log_inverse_exact="
        f"{exp_m_ok}"
    )

    print(
        f"  final_m={final_m_int}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 214RR COMPLETE")


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