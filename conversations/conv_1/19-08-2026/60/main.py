#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 227 — EXACT SCHUR-PRIME / PROJECTIVE-RATIO CONGRUENCE AUDIT
==============================================================================

Known case:

    p = 2
    q = 3
    n = 6

    q1 = 29144191
    q3 = 24794967

Terminal Schur coefficients:

    s0 = q1 - 2*q3
    s1 = q1 - 6*q3
    s2 = q1 - 18*q3

Experiment 226 found:

    source primes that cancel:
        2, 3, 839, 9851, 229583, 971473, 29144191

    Schur-created primes:
        7, 53, 73, 137, 631, 1187, 70289, 149239.

Experiment 227 asks the structural question:

    Why do these new primes appear?

For every prime ell not dividing q3,

    ell | (q1 - c*q3)

is equivalent to

    q1/q3 = c (mod ell),

where

    c in {2,6,18}.

Thus every prime divisor of a Schur coefficient should correspond
to an exact projective-ratio congruence.

The experiment checks both directions:

    divisor -> ratio congruence

and, for a finite prime range,

    ratio congruence -> divisor.

It also checks the exceptional case ell | q3 separately, because
division by q3 is then impossible modulo ell.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN INSTANCE
# ============================================================================

P = 2
Q = 3
N = P * Q

Q1 = 29144191
Q3 = 24794967

COEFFICIENTS = {
    "s0": 2,
    "s1": 6,
    "s2": 18,
}

SEARCH_LIMIT = 20000


# ============================================================================
# EXACT HELPERS
# ============================================================================

def factor_integer(
    x: int,
) -> dict[int, int]:

    x = abs(int(x))

    if x < 2:
        return {}

    factors: dict[int, int] = {}

    while x % 2 == 0:

        factors[2] = factors.get(2, 0) + 1
        x //= 2

    d = 3

    while d * d <= x:

        while x % d == 0:

            factors[d] = factors.get(d, 0) + 1
            x //= d

        d += 2

    if x > 1:
        factors[x] = factors.get(x, 0) + 1

    return factors


def is_prime(
    n: int,
) -> bool:

    n = int(n)

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


def primes_up_to(
    limit: int,
) -> list[int]:

    primes = []

    for n in range(
        2,
        int(limit) + 1,
    ):

        if is_prime(n):
            primes.append(n)

    return primes


def inverse_mod(
    a: int,
    m: int,
) -> int:

    a = int(a)
    m = int(m)

    if m <= 1:
        return 0

    a %= m

    if gcd(
        a,
        m,
    ) != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return pow(
        a,
        -1,
        m,
    )


def ratio_mod(
    q1: int,
    q3: int,
    ell: int,
) -> int | None:

    ell = int(ell)

    if gcd(
        int(q3),
        ell,
    ) != 1:
        return None

    return int(
        (
            (int(q1) % ell)
            * inverse_mod(
                int(q3),
                ell,
            )
        )
        % ell
    )


def valuation(
    x: int,
    ell: int,
) -> int | None:

    x = int(x)
    ell = int(ell)

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % ell == 0:

        x //= ell
        v += 1

    return v


# ============================================================================
# SCHUR VALUES
# ============================================================================

S0 = int(
    Q1 - COEFFICIENTS["s0"] * Q3
)

S1 = int(
    Q1 - COEFFICIENTS["s1"] * Q3
)

S2 = int(
    Q1 - COEFFICIENTS["s2"] * Q3
)

SCHUR = {
    "s0": S0,
    "s1": S1,
    "s2": S2,
}


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 227 — EXACT SCHUR-PRIME / "
        "PROJECTIVE-RATIO CONGRUENCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(
        f"  p={P}"
    )

    print(
        f"  q={Q}"
    )

    print(
        f"  n={N}"
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. EXACT SCHUR VALUES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT SCHUR COEFFICIENTS")
    print("=" * 78)

    for name in (
        "s0",
        "s1",
        "s2",
    ):

        value = SCHUR[name]

        print(
            f"  {name}="
            f"{value}"
        )

        print(
            f"    coefficient="
            f"{COEFFICIENTS[name]}"
        )

        print(
            f"    factorization="
            f"{factor_integer(value)}"
        )

    # ------------------------------------------------------------------
    # 3. PRIME DIVISORS OF SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    schur_prime_map = {}

    for name, value in SCHUR.items():

        schur_prime_map[name] = sorted(
            factor_integer(value).keys()
        )

    print()
    print("=" * 78)
    print("3. SCHUR PRIME DIVISORS")
    print("=" * 78)

    for name in (
        "s0",
        "s1",
        "s2",
    ):

        print(
            f"  {name}: "
            f"{schur_prime_map[name]}"
        )

    # ------------------------------------------------------------------
    # 4. DIVISOR -> PROJECTIVE CONGRUENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DIVISOR -> PROJECTIVE RATIO CONGRUENCE")
    print("=" * 78)

    divisor_forward_failures = []

    for name, coeff in COEFFICIENTS.items():

        value = SCHUR[name]

        for ell in schur_prime_map[name]:

            rho = ratio_mod(
                Q1,
                Q3,
                ell,
            )

            expected = int(
                coeff % ell
            )

            actual_divisibility = (
                value % ell == 0
            )

            ratio_match = (
                rho is not None
                and rho == expected
            )

            print(
                f"  {name}: ell={ell} "
                f"q3_invertible="
                f"{rho is not None} "
                f"rho={rho} "
                f"expected={expected} "
                f"divides={actual_divisibility} "
                f"ratio_match={ratio_match}"
            )

            if not ratio_match:

                divisor_forward_failures.append(
                    (
                        name,
                        ell,
                        rho,
                        expected,
                    )
                )

    # ------------------------------------------------------------------
    # 5. EXCEPTIONAL q3 PRIMES
    # ------------------------------------------------------------------

    q3_primes = sorted(
        factor_integer(Q3).keys()
    )

    print()
    print("=" * 78)
    print("5. PRIMES DIVIDING q3")
    print("=" * 78)

    print(
        f"  q3_prime_divisors={q3_primes}"
    )

    for ell in q3_primes:

        rho = ratio_mod(
            Q1,
            Q3,
            ell,
        )

        print(
            f"  ell={ell}: "
            f"q3_divisible=True "
            f"rho={rho}"
        )

        for name, coeff in COEFFICIENTS.items():

            print(
                f"    {name}: "
                f"valuation={valuation(SCHUR[name],ell)}"
            )

    # ------------------------------------------------------------------
    # 6. CREATED PRIME CLASSIFICATION
    # ------------------------------------------------------------------

    source_primes = set(
        factor_integer(Q1).keys()
    )

    source_primes.update(
        factor_integer(Q3).keys()
    )

    schur_primes = set()

    for plist in schur_prime_map.values():
        schur_primes.update(
            plist
        )

    created_primes = sorted(
        schur_primes
        - source_primes
    )

    canceled_primes = sorted(
        source_primes
        - schur_primes
    )

    print()
    print("=" * 78)
    print("6. CREATED / CANCELED PRIME CLASSES")
    print("=" * 78)

    print(
        f"  source_primes="
        f"{sorted(source_primes)}"
    )

    print(
        f"  schur_primes="
        f"{sorted(schur_primes)}"
    )

    print(
        f"  created_primes="
        f"{created_primes}"
    )

    print(
        f"  canceled_primes="
        f"{canceled_primes}"
    )

    # ------------------------------------------------------------------
    # 7. RATIO VALUES AT CREATED PRIMES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PROJECTIVE VALUES AT CREATED SCHUR PRIMES")
    print("=" * 78)

    created_match_failures = []

    for ell in created_primes:

        rho = ratio_mod(
            Q1,
            Q3,
            ell,
        )

        print(
            f"  ell={ell}: "
            f"rho={rho}"
        )

        matched = []

        for name, coeff in COEFFICIENTS.items():

            if rho == coeff % ell:

                matched.append(
                    name
                )

        print(
            f"    matched_schur_rows="
            f"{matched}"
        )

        if len(matched) == 0:

            created_match_failures.append(
                (ell, rho)
            )

    # ------------------------------------------------------------------
    # 8. FINITE CONVERSE TEST
    # ------------------------------------------------------------------
    #
    # Test all primes up to SEARCH_LIMIT for:
    #
    #    rho == c mod ell
    #
    # versus
    #
    #    ell | s_c.
    #
    # This gives a finite exact bidirectional test.
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "8. FINITE CONVERSE TEST: "
        "PROJECTIVE CONGRUENCE -> SCHUR DIVISIBILITY"
    )
    print("=" * 78)

    converse_failures = []

    checked = 0

    for ell in primes_up_to(
        SEARCH_LIMIT
    ):

        rho = ratio_mod(
            Q1,
            Q3,
            ell,
        )

        if rho is None:
            continue

        checked += 1

        for name, coeff in COEFFICIENTS.items():

            condition_ratio = (
                rho == coeff % ell
            )

            condition_divides = (
                SCHUR[name] % ell == 0
            )

            if condition_ratio != condition_divides:

                converse_failures.append(
                    (
                        ell,
                        name,
                        rho,
                        coeff % ell,
                        condition_ratio,
                        condition_divides,
                    )
                )

    print(
        f"  primes_checked={checked}"
    )

    print(
        f"  converse_failures="
        f"{converse_failures}"
    )

    # ------------------------------------------------------------------
    # 9. PRIME-POWER VALUATION COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. PRIME-POWER VALUATION COMPARISON")
    print("=" * 78)

    for name, coeff in COEFFICIENTS.items():

        print(
            f"  {name}: coefficient={coeff}"
        )

        for ell in schur_prime_map[name]:

            v_s = valuation(
                SCHUR[name],
                ell,
            )

            print(
                f"    ell={ell}: "
                f"v_ell(s)={v_s}"
            )

            if rho_mod := ratio_mod(
                Q1,
                Q3,
                ell,
            ):

                defect = int(
                    rho_mod
                    - (coeff % ell)
                ) % ell

                print(
                    f"      rho_mod={rho_mod} "
                    f"coefficient_mod={coeff % ell} "
                    f"first_defect_mod_ell={defect}"
                )

    # ------------------------------------------------------------------
    # 10. EXACT DIFFERENCE CONNECTION
    # ------------------------------------------------------------------

    d01 = S1 - S0
    d12 = S2 - S1
    d02 = S2 - S0

    print()
    print("=" * 78)
    print("10. DIFFERENCE CONNECTION")
    print("=" * 78)

    print(
        f"  s1-s0={d01}"
    )

    print(
        f"  expected=-4*q3="
        f"{-4*Q3}"
    )

    print(
        f"  s2-s1={d12}"
    )

    print(
        f"  expected=-12*q3="
        f"{-12*Q3}"
    )

    print(
        f"  s2-s0={d02}"
    )

    print(
        f"  expected=-16*q3="
        f"{-16*Q3}"
    )

    # ------------------------------------------------------------------
    # 11. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The Schur coefficients are

    s_k = q1 - c_k q3,

with

    c_k in {2,6,18}.

Whenever ell does not divide q3,

    ell | s_k

is equivalent to

    q1/q3 = c_k (mod ell).

Therefore the prime factorization of the Schur row is not arbitrary:
its newly created primes are precisely candidates for primes at which
the source projective ratio lands on one of the Schur constants.

This experiment tests that statement directly.

The finite converse test is especially important. It checks that, among
all primes up to the chosen bound, no prime satisfies the projective
congruence without actually dividing the corresponding Schur coefficient.

If both directions pass, the prime spectrum of the terminal Schur row
has an exact projective interpretation.

This is stronger and more structural than simply factoring the large
integers.

The statement remains elementary and exact; it does not yet depend on
the special n=pq construction.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    forward_exact = (
        len(divisor_forward_failures) == 0
    )

    created_exact = (
        len(created_match_failures) == 0
    )

    converse_exact = (
        len(converse_failures) == 0
    )

    difference_exact = (
        d01 == -4 * Q3
        and d12 == -12 * Q3
        and d02 == -16 * Q3
    )

    final_ok = (
        forward_exact
        and created_exact
        and converse_exact
        and difference_exact
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  divisor_to_ratio_exact="
        f"{forward_exact}"
    )

    print(
        f"  created_prime_projective_explanation_exact="
        f"{created_exact}"
    )

    print(
        f"  finite_converse_exact="
        f"{converse_exact}"
    )

    print(
        f"  Schur_difference_identities_exact="
        f"{difference_exact}"
    )

    print(
        f"  created_primes="
        f"{created_primes}"
    )

    print(
        f"  canceled_primes="
        f"{canceled_primes}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 227 COMPLETE")


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

