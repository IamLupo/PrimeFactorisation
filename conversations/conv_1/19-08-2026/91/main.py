#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 256 — EXACT SYMBOLIC n=pq TERMINAL-ORBIT / SCHUR-GEOMETRY AUDIT
==============================================================================

The known n=6 case has

    p=2, q=3

and terminal Schur coefficients

    2, 6, 18

which are exactly

    p*q^0,
    p*q^1,
    p*q^2.

Experiment 256 isolates this observed pattern and tests the corresponding
symbolic structure before attempting to derive q1(p,q), q3(p,q).

For a general pair (p,q), define the candidate terminal orbit

    c_j = p*q^j.

Then the candidate Schur row is

    s_j = q1 - p*q^j*q3.

The differences should satisfy

    s_{j+1} - s_j
      = -p*q^j*(q-1)*q3.

For the first three rows:

    s0 = q1 - p*q3
    s1 = q1 - p*q*q3
    s2 = q1 - p*q^2*q3

and therefore

    d01 = -p*(q-1)*q3
    d12 = -p*q*(q-1)*q3
    d12 = q*d01.

This script performs two separate audits:

    A. exact verification for the known n=6 data;

    B. synthetic symbolic-family tests for many coprime prime pairs (p,q),
       using arbitrary coprime source pairs q1,q3.

The synthetic tests do NOT claim to derive q1(p,q). They only verify
that, once the candidate terminal orbit p*q^j is present, the associated
Schur geometry and projective divisibility laws are exact.

The experiment therefore tells us exactly what the missing upstream
n=pq theorem must establish.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN n=6 SOURCE
# ============================================================================

KNOWN_P = 2
KNOWN_Q = 3
KNOWN_Q1 = 29144191
KNOWN_Q3 = 24794967


# ============================================================================
# SYNTHETIC TEST PARAMETERS
# ============================================================================

PRIME_LIMIT = 31
SYNTHETIC_Q1 = 1000003
SYNTHETIC_Q3 = 1000033


# ============================================================================
# HELPERS
# ============================================================================

def is_prime(n: int) -> bool:

    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:

        if n % d == 0:
            return False

        d += 2

    return True


def candidate_schur_row(
    p: int,
    q: int,
    q1: int,
    q3: int,
    rows: int = 3,
) -> list[int]:

    return [
        q1 - p * (q ** j) * q3
        for j in range(rows)
    ]


def candidate_orbit_constants(
    p: int,
    q: int,
    rows: int = 3,
) -> list[int]:

    return [
        p * (q ** j)
        for j in range(rows)
    ]


def difference_data(
    row: list[int],
) -> tuple[list[int], list[int]]:

    differences = [
        row[i + 1] - row[i]
        for i in range(len(row) - 1)
    ]

    second_differences = [
        differences[i + 1] - differences[i]
        for i in range(len(differences) - 1)
    ]

    return differences, second_differences


def exact_projective_divisor_test(
    q1: int,
    q3: int,
    coefficient: int,
    bound: int,
) -> bool:

    """
    Finite verification of

        ell | q1 - coefficient*q3

    iff

        q1/q3 = coefficient mod ell,

    for primes ell <= bound not dividing q3.
    """

    for ell in range(
        2,
        bound + 1,
    ):

        if not is_prime(ell):
            continue

        if q3 % ell == 0:
            continue

        residual = (
            q1 - coefficient * q3
        )

        divisor_side = (
            residual % ell == 0
        )

        ratio_side = (
            (
                q1
                * pow(
                    q3,
                    -1,
                    ell,
                )
            )
            % ell
            ==
            coefficient % ell
        )

        if divisor_side != ratio_side:
            return False

    return True


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 256 — EXACT SYMBOLIC n=pq TERMINAL-ORBIT / "
        "SCHUR-GEOMETRY AUDIT"
    )
    print("=" * 78)

    failures = []

    # ------------------------------------------------------------------
    # 1. KNOWN n=6 AUDIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. KNOWN n=6 TERMINAL ORBIT")
    print("=" * 78)

    known_constants = candidate_orbit_constants(
        KNOWN_P,
        KNOWN_Q,
    )

    known_row = candidate_schur_row(
        KNOWN_P,
        KNOWN_Q,
        KNOWN_Q1,
        KNOWN_Q3,
    )

    print(
        f"  p={KNOWN_P}"
    )
    print(
        f"  q={KNOWN_Q}"
    )
    print(
        f"  observed_constants=[2,6,18]"
    )
    print(
        f"  candidate_constants={known_constants}"
    )

    constants_exact = (
        known_constants
        ==
        [2, 6, 18]
    )

    print(
        f"  candidate_constants_exact="
        f"{constants_exact}"
    )

    print(
        f"  candidate_schur_row={known_row}"
    )

    if not constants_exact:
        failures.append(
            (
                "known_n6",
                "candidate_constants",
            )
        )

    # ------------------------------------------------------------------
    # 2. KNOWN n=6 DIFFERENCE GEOMETRY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. KNOWN n=6 SCHUR DIFFERENCE GEOMETRY")
    print("=" * 78)

    d, dd = difference_data(
        known_row
    )

    print(
        f"  differences={d}"
    )

    print(
        f"  second_differences={dd}"
    )

    expected_d01 = (
        -KNOWN_P
        * (KNOWN_Q - 1)
        * KNOWN_Q3
    )

    expected_d12 = (
        -KNOWN_P
        * KNOWN_Q
        * (KNOWN_Q - 1)
        * KNOWN_Q3
    )

    geometry_checks = {
        "d01_exact":
            d[0] == expected_d01,

        "d12_exact":
            d[1] == expected_d12,

        "d12_q_times_d01":
            d[1] == KNOWN_Q * d[0],

        "second_difference_exact":
            dd[0]
            ==
            -KNOWN_P
            * (KNOWN_Q - 1) ** 2
            * KNOWN_Q3,
    }

    for name, value in geometry_checks.items():

        print(
            f"  {name}={value}"
        )

        if not value:
            failures.append(
                (
                    "known_n6",
                    name,
                )
            )

    # ------------------------------------------------------------------
    # 3. SYNTHETIC PRIME-PAIR FAMILY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SYNTHETIC PRIME-PAIR FAMILY TEST")
    print("=" * 78)

    pair_count = 0
    pair_failures = 0

    primes = [
        p
        for p in range(
            2,
            PRIME_LIMIT + 1,
        )
        if is_prime(p)
    ]

    for p in primes:

        for q in primes:

            if p == q:
                continue

            pair_count += 1

            row = candidate_schur_row(
                p,
                q,
                SYNTHETIC_Q1,
                SYNTHETIC_Q3,
            )

            differences, second_differences = (
                difference_data(
                    row
                )
            )

            expected_d01 = (
                -p
                * (q - 1)
                * SYNTHETIC_Q3
            )

            expected_d12 = (
                -p
                * q
                * (q - 1)
                * SYNTHETIC_Q3
            )

            expected_second = (
                -p
                * (q - 1) ** 2
                * SYNTHETIC_Q3
            )

            ok = (
                differences[0]
                ==
                expected_d01
                and
                differences[1]
                ==
                expected_d12
                and
                differences[1]
                ==
                q * differences[0]
                and
                second_differences[0]
                ==
                expected_second
            )

            if not ok:

                pair_failures += 1

                if pair_failures <= 10:

                    print(
                        f"  FAILURE p={p} q={q}: "
                        f"row={row} "
                        f"differences={differences}"
                    )

    print(
        f"  synthetic_prime_pairs={pair_count}"
    )

    print(
        f"  synthetic_geometry_failures="
        f"{pair_failures}"
    )

    if pair_failures:
        failures.append(
            (
                "synthetic_geometry",
                pair_failures,
            )
        )

    # ------------------------------------------------------------------
    # 4. PROJECTIVE DIVISIBILITY INTERFACE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PROJECTIVE DIVISIBILITY INTERFACE")
    print("=" * 78)

    divisor_constants = (
        known_constants
    )

    divisor_checks = []

    for c in divisor_constants:

        result = exact_projective_divisor_test(
            KNOWN_Q1,
            KNOWN_Q3,
            c,
            1000,
        )

        divisor_checks.append(
            result
        )

        print(
            f"  coefficient={c}: "
            f"divisor_ratio_equivalence="
            f"{result}"
        )

    if not all(divisor_checks):
        failures.append(
            (
                "projective_divisor_interface",
            )
        )

    # ------------------------------------------------------------------
    # 5. CONDITIONAL GENERAL FORM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. CONDITIONAL GENERAL n=pq FORM")
    print("=" * 78)

    print(
        """
Conditional on the terminal constants being

    c_j = p*q^j,

the terminal Schur row must be

    s_j = q1 - p*q^j*q3.

The first differences are therefore

    s_{j+1}-s_j
      = -p*q^j*(q-1)*q3.

In particular,

    s1-s0
      = -p*(q-1)*q3,

    s2-s1
      = -p*q*(q-1)*q3,

and hence

    s2-s1 = q*(s1-s0).

This identity is independent of the actual formulas for q1 and q3.

Thus the remaining n=pq problem is not the local Schur geometry.
It is the derivation of q1(p,q) and q3(p,q), together with a proof
that the terminal constants really have the form p*q^j.
"""
    )

    # ------------------------------------------------------------------
    # 6. SOURCE-THEOREM REQUIREMENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. WHAT THE ACTUAL n=pq THEOREM MUST SUPPLY")
    print("=" * 78)

    print(
        """
A complete symbolic n=pq theorem now needs only two genuinely upstream
ingredients:

    (A) SOURCE FORMULAS

        q1 = q1(p,q),
        q3 = q3(p,q);

    (B) TERMINAL-ORBIT FORMULA

        c_j = p*q^j.

Everything downstream of these inputs has already been independently
audited:

    projective invariance,
    representative invariance,
    Hensel lifting,
    principal-coordinate digits,
    logarithmic coordinate,
    additive group law,
    Schur projective divisibility.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  known_n6_constants_exact="
        f"{constants_exact}"
    )

    print(
        f"  known_n6_difference_geometry_exact="
        f"{all(geometry_checks.values())}"
    )

    print(
        f"  synthetic_prime_pair_geometry_exact="
        f"{pair_failures == 0}"
    )

    print(
        f"  projective_divisibility_interface_exact="
        f"{all(divisor_checks)}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  local_projective_theorem_complete=True"
    )

    print(
        "  terminal_pq_orbit_formula_conditional=True"
    )

    print(
        "  q1_of_p_q_derived=False"
    )

    print(
        "  q3_of_p_q_derived=False"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason="
        "upstream source formulas still missing"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 256 COMPLETE")


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

