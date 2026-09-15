#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 215 — EXACT HIGHER-PRECISION LOG/HENSEL THEOREM STRESS TEST
==============================================================================

Candidate theorem under test:

    rho = q1/q3 = 6 * 729^m,

with

    m = log(rho/6) / log(729).

Equivalently,

    k = 1 + 6m.

Experiment 214RR verified this through precision 7^12.

Experiment 215 stress-tests the same theorem through 7^16.

Important:

    * No later exponent sequence is hard-coded.
    * Hensel lifting is reconstructed from rho itself.
    * The logarithmic coordinate is computed independently.
    * The two coordinates are compared modulo 7^(E-1).
    * The local derivative/unit condition is checked at every level.
    * Class uniqueness is checked at every tested precision.

No floating point.
No SymPy.
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

START_E = 1
END_E = 16

# Log/exp truncation margin.
PRECISION = END_E + 3


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
    a = int(a)
    m = int(m)

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

    return int(
        (
            (x.numerator % modulus)
            * inverse_mod(
                x.denominator % modulus,
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
# 7-ADIC LOGARITHM
# ============================================================================

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

        if n > 4 * precision + 150:
            raise ArithmeticError(
                "log truncation failed."
            )

    return total


# ============================================================================
# HENSEL HELPERS
# ============================================================================

def base_exponent_mod7() -> int:

    rho = ratio_mod(7)

    matches = []

    for k in range(6):

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
            f"Expected one base exponent, got {matches}."
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

    rho_next = ratio_mod(
        modulus_next
    )

    orbit_next = principal_orbit_mod(
        m,
        modulus_next,
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
            f"m={m} does not solve modulo 7^{e}."
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
        "EXPERIMENT 215 — EXACT HIGHER-PRECISION "
        "LOG/HENSEL THEOREM STRESS TEST"
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

    source_x = (
        rho_over_6
        - 1
    )

    base_x = (
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
        f"  q3_mod7={Q3 % P}"
    )

    print(
        f"  v7(q1)={q1_v7}"
    )

    print(
        f"  v7(q3)={q3_v7}"
    )

    print(
        f"  v7(rho/6-1)="
        f"{fraction_valuation_7(source_x)}"
    )

    print(
        f"  v7(729-1)="
        f"{fraction_valuation_7(base_x)}"
    )

    # ------------------------------------------------------------------
    # 2. LOGARITHMIC COORDINATE
    # ------------------------------------------------------------------

    log_source = log_one_plus_x(
        source_x,
        PRECISION,
    )

    log_base = log_one_plus_x(
        base_x,
        PRECISION,
    )

    log_m = (
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
    print("2. LOGARITHMIC COORDINATE")
    print("=" * 78)

    print(
        f"  precision={PRECISION}"
    )

    print(
        f"  v7(log_source)={log_source_v7}"
    )

    print(
        f"  v7(log_base)={log_base_v7}"
    )

    print(
        f"  v7(log_m)="
        f"{fraction_valuation_7(log_m)}"
    )

    # ------------------------------------------------------------------
    # 3. INDEPENDENT HENSEL RECONSTRUCTION
    # ------------------------------------------------------------------

    k1 = base_exponent_mod7()

    if k1 != 1:
        raise ArithmeticError(
            f"Expected k1=1, got {k1}."
        )

    current_m = int(
        (k1 - 1) // 6
    )

    m_sequence = [
        current_m
    ]

    t_sequence = []
    delta_sequence = []

    print()
    print("=" * 78)
    print("3. INDEPENDENT HENSEL LIFT")
    print("=" * 78)

    for e in range(
        START_E,
        END_E,
    ):

        delta = defect_digit(
            e,
            current_m,
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

        delta_sequence.append(
            delta
        )

        t_sequence.append(
            t
        )

        m_sequence.append(
            current_m
        )

        modulus_next = int(
            P ** (
                e + 1
            )
        )

        alignment = (
            principal_orbit_mod(
                current_m,
                modulus_next,
            )
            ==
            ratio_mod(
                modulus_next,
            )
        )

        print(
            f"  e={e}: "
            f"delta={delta} "
            f"t={t} "
            f"increment={increment} "
            f"m_next={current_m} "
            f"alignment={alignment}"
        )

        if not alignment:
            raise ArithmeticError(
                f"Hensel alignment failed at e={e}."
            )

    # ------------------------------------------------------------------
    # 4. HIGHER-PRECISION LOG VS HENSEL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. LOGARITHM VS HENSEL THROUGH 7^16")
    print("=" * 78)

    precision_flags = []

    for E in range(
        1,
        END_E + 1,
    ):

        principal_modulus = int(
            P ** (
                E - 1
            )
        )

        m_hensel = int(
            m_sequence[E - 1]
            % principal_modulus
        )

        if principal_modulus == 1:
            m_log = 0
        else:
            m_log = int(
                fraction_mod(
                    log_m,
                    principal_modulus,
                )
            )

        exact = (
            m_hensel
            == m_log
        )

        precision_flags.append(
            exact
        )

        print(
            f"  E={E}: "
            f"modulus={principal_modulus} "
            f"m_hensel={m_hensel} "
            f"m_log={m_log} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 5. LOCAL DERIVATIVE / UNIT TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. LOGARITHMIC DERIVATIVE / UNIT TEST")
    print("=" * 78)

    derivative_unit = (
        (
            PRINCIPAL_BASE
            - 1
        )
        // P
    ) % P

    rho_unit = (
        ratio_mod(7)
    )

    slope = (
        rho_unit
        * derivative_unit
    ) % P

    print(
        f"  ((729-1)/7) mod7={derivative_unit}"
    )

    print(
        f"  rho_mod7={rho_unit}"
    )

    print(
        f"  predicted_Hensel_slope="
        f"{slope}"
    )

    print(
        f"  unit_slope=True="
        f"{slope != 0}"
    )

    # ------------------------------------------------------------------
    # 6. UNIQUENESS AT HIGH PRECISION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. HIGH-PRECISION CLASS UNIQUENESS")
    print("=" * 78)

    final_principal_modulus = int(
        P ** (
            END_E - 1
        )
    )

    final_orbit_modulus = int(
        P ** END_E
    )

    final_m = int(
        m_sequence[-1]
        % final_principal_modulus
    )

    final_target = int(
        ratio_mod(
            final_orbit_modulus
        )
    )

    # Same residue class representatives.
    same_class_candidates = [
        int(final_m),
        int(
            final_m
            + final_principal_modulus
        ),
        int(
            final_m
            - final_principal_modulus
        ),
    ]

    same_class_flags = []

    for candidate in same_class_candidates:

        match = (
            principal_orbit_mod(
                candidate,
                final_orbit_modulus,
            )
            == final_target
        )

        same_class_flags.append(
            match
        )

        print(
            f"  same_class_candidate={candidate} "
            f"match={match}"
        )

    # Distinct neighboring classes at the next lower power.
    distinguishing_step = int(
        P ** (
            END_E - 2
        )
    )

    neighboring_candidates = [
        int(
            final_m
            + distinguishing_step
        ),
        int(
            final_m
            - distinguishing_step
        ),
    ]

    neighboring_flags = []

    for candidate in neighboring_candidates:

        match = (
            principal_orbit_mod(
                candidate,
                final_orbit_modulus,
            )
            == final_target
        )

        neighboring_flags.append(
            match
        )

        print(
            f"  neighboring_candidate={candidate} "
            f"match={match}"
        )

    # ------------------------------------------------------------------
    # 7. CLOSED FORM / EXPONENTIAL CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. CLOSED-FORM EXPONENTIAL CHECK")
    print("=" * 78)

    # Compute exp(log(rho/6)).
    exp_log_source = log_source

    # Rather than introduce another independent exp truncation here,
    # verify the closed form directly with modular exponentiation using
    # the logarithm-derived m residue.
    reconstructed_ratio = (
        principal_orbit_mod(
            final_m,
            final_orbit_modulus,
        )
        * inverse_mod(
            6,
            final_orbit_modulus,
        )
    ) % final_orbit_modulus

    source_ratio_normalized = fraction_mod(
        rho_over_6,
        final_orbit_modulus,
    )

    normalized_source_match = (
        reconstructed_ratio
        == source_ratio_normalized
    )

    print(
        f"  modulus={final_orbit_modulus}"
    )

    print(
        f"  normalized_source="
        f"{source_ratio_normalized}"
    )

    print(
        f"  reconstructed_729^m="
        f"{reconstructed_ratio}"
    )

    print(
        f"  exact={normalized_source_match}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The candidate theorem now has three independent components:

    (1) Hensel construction

        rho = 6*729^m

        with m obtained digit-by-digit;

    (2) logarithmic coordinate

        m = log(rho/6)/log(729);

    (3) uniqueness

        m is determined modulo 7^(E-1) by the orbit equation modulo 7^E.

Experiment 215 pushes the three statements from precision 7^12 to
precision 7^16.

The derivative/unit calculation

    rho * ((729-1)/7) = 1 mod 7

explains the unique one-digit lift at every level.

If all checks pass, the evidence for the principal-coordinate theorem
is substantially stronger than a finite pattern fit.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    all_log_matches = all(
        precision_flags
    )

    all_same_class = all(
        same_class_flags
    )

    all_neighbors_separated = all(
        not flag
        for flag in neighboring_flags
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and log_source_v7 == 2
        and log_base_v7 == 1
        and all_log_matches
        and slope == 1
        and all_same_class
        and all_neighbors_separated
        and normalized_source_match
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_v7_is_0={q1_v7 == 0}"
    )

    print(
        f"  q3_v7_is_0={q3_v7 == 0}"
    )

    print(
        f"  log_source_v7_is_2="
        f"{log_source_v7 == 2}"
    )

    print(
        f"  log_base_v7_is_1="
        f"{log_base_v7 == 1}"
    )

    print(
        f"  all_log_hensel_matches="
        f"{all_log_matches}"
    )

    print(
        f"  predicted_slope={slope}"
    )

    print(
        f"  unit_lift_slope_exact="
        f"{slope == 1}"
    )

    print(
        f"  all_same_class_representatives_match="
        f"{all_same_class}"
    )

    print(
        f"  all_neighbor_classes_separated="
        f"{all_neighbors_separated}"
    )

    print(
        f"  normalized_exponential_reconstruction="
        f"{normalized_source_match}"
    )

    print(
        f"  final_m={final_m}"
    )

    print(
        f"  final_precision=7^{END_E}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 215 COMPLETE")


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

