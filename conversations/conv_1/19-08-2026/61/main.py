#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 228 — EXACT SCHUR PRIME-POWER / PROJECTIVE-DEPTH AUDIT
==============================================================================

Known case:

    q1 = 29144191
    q3 = 24794967

Schur coefficients:

    s0 = q1 - 2*q3
    s1 = q1 - 6*q3
    s2 = q1 - 18*q3

Experiment 227 established, for all tested primes ell <= 20000,

    ell | s_c
        iff
    q1/q3 = c (mod ell),

provided ell does not divide q3.

Experiment 228 upgrades this to prime powers.

For each relevant prime ell and each Schur coefficient c, test:

    v_ell(s_c) >= e

against

    q1/q3 == c (mod ell^e).

Then determine:

    * exact congruence depth;
    * exact Schur valuation;
    * first nonzero projective defect digit;
    * whether the valuation equals the projective distance.

This is the finite-prime analogue of the 7-adic projective analysis.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN INSTANCE
# ============================================================================

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

            factors[d] = (
                factors.get(d, 0) + 1
            )

            x //= d

        d += 2

    if x > 1:

        factors[x] = (
            factors.get(x, 0) + 1
        )

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

    return [
        n
        for n in range(
            2,
            int(limit) + 1,
        )
        if is_prime(n)
    ]


def inverse_mod(
    a: int,
    m: int,
) -> int:

    a = int(a)
    m = int(m)

    if gcd(
        a,
        m,
    ) != 1:

        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return pow(
        a % m,
        -1,
        m,
    )


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int | None:

    modulus = int(modulus)

    if gcd(
        int(q3),
        modulus,
    ) != 1:

        return None

    return int(
        (
            (int(q1) % modulus)
            * inverse_mod(
                q3,
                modulus,
            )
        )
        % modulus
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 228 — EXACT SCHUR PRIME-POWER / "
        "PROJECTIVE-DEPTH AUDIT"
    )
    print("=" * 78)

    schur_values = {
        "s0": int(
            Q1 - 2 * Q3
        ),
        "s1": int(
            Q1 - 6 * Q3
        ),
        "s2": int(
            Q1 - 18 * Q3
        ),
    }

    # ------------------------------------------------------------------
    # 1. SOURCE DATA
    # ------------------------------------------------------------------

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
        f"  gcd(q1,q3)="
        f"{gcd(Q1,Q3)}"
    )

    print(
        f"  q3_factorization="
        f"{factor_integer(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. SCHUR DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SCHUR DATA")
    print("=" * 78)

    for name in (
        "s0",
        "s1",
        "s2",
    ):

        value = schur_values[name]

        print(
            f"  {name}="
            f"{value} "
            f"c={COEFFICIENTS[name]} "
            f"factorization="
            f"{factor_integer(value)}"
        )

    # ------------------------------------------------------------------
    # 3. PRIME-POWER DEPTH TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT PRIME-POWER PROJECTIVE DEPTH")
    print("=" * 78)

    relevant_primes = set()

    for value in schur_values.values():

        relevant_primes.update(
            factor_integer(value).keys()
        )

    relevant_primes = sorted(
        relevant_primes
    )

    all_prime_power_checks = True
    depth_records = []

    for name, coeff in COEFFICIENTS.items():

        s = schur_values[name]

        print()
        print(
            f"  {name}: coefficient={coeff}"
        )

        for ell in relevant_primes:

            if gcd(
                Q3,
                ell,
            ) != 1:

                continue

            v = valuation(
                s,
                ell,
            )

            if v is None or v == 0:
                continue

            print(
                f"    ell={ell}: "
                f"v_ell(s)={v}"
            )

            for e in range(
                1,
                v + 2,
            ):

                modulus = int(
                    ell ** e
                )

                rho = ratio_mod(
                    Q1,
                    Q3,
                    modulus,
                )

                congruent = (
                    rho
                    ==
                    (
                        coeff
                        % modulus
                    )
                )

                expected = (
                    e <= v
                )

                exact = (
                    congruent
                    == expected
                )

                print(
                    f"      e={e} "
                    f"mod={modulus} "
                    f"rho={rho} "
                    f"c_mod={coeff % modulus} "
                    f"match={congruent} "
                    f"expected={expected} "
                    f"exact={exact}"
                )

                if not exact:
                    all_prime_power_checks = False

            depth_records.append(
                (
                    name,
                    ell,
                    v,
                )
            )

    # ------------------------------------------------------------------
    # 4. EXACT MAXIMUM PROJECTIVE DEPTH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. MAXIMUM PROJECTIVE DEPTH")
    print("=" * 78)

    for name, ell, v in depth_records:

        print(
            f"  {name}: "
            f"ell={ell} "
            f"depth={v}"
        )

    maximum_depth = 0
    maximum_records = []

    for name, ell, v in depth_records:

        if v > maximum_depth:

            maximum_depth = v
            maximum_records = [
                (
                    name,
                    ell,
                )
            ]

        elif v == maximum_depth:

            maximum_records.append(
                (
                    name,
                    ell,
                )
            )

    print(
        f"  maximum_depth="
        f"{maximum_depth}"
    )

    print(
        f"  maximizers="
        f"{maximum_records}"
    )

    # ------------------------------------------------------------------
    # 5. FIRST NONZERO PROJECTIVE DEFECT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FIRST NONZERO PROJECTIVE DEFECT")
    print("=" * 78)

    defect_checks = True

    for name, ell, v in depth_records:

        coeff = COEFFICIENTS[name]

        modulus = int(
            ell ** (v + 1)
        )

        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        if rho is None:
            defect_checks = False
            continue

        defect = int(
            (
                rho
                - (
                    coeff
                    % modulus
                )
            )
            % modulus
        )

        required_divisor = int(
            ell ** v
        )

        divisible = (
            defect
            % required_divisor
            == 0
        )

        normalized = int(
            (
                defect
                // required_divisor
            )
            % ell
        )

        nonzero = (
            normalized != 0
        )

        exact = (
            divisible
            and
            nonzero
        )

        print(
            f"  {name}: "
            f"ell={ell} "
            f"v={v} "
            f"rho_mod_ell^(v+1)={rho} "
            f"c_mod={coeff % modulus} "
            f"defect={defect} "
            f"normalized_digit={normalized} "
            f"exact={exact}"
        )

        if not exact:
            defect_checks = False

    # ------------------------------------------------------------------
    # 6. SPECIAL 7-ADIC COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. 7-ADIC SPECIAL CASE")
    print("=" * 78)

    s1 = schur_values["s1"]

    v7_s1 = valuation(
        s1,
        7,
    )

    rho7 = ratio_mod(
        Q1,
        Q3,
        7 ** 3,
    )

    print(
        f"  s1="
        f"{s1}"
    )

    print(
        f"  v7(s1)="
        f"{v7_s1}"
    )

    print(
        f"  rho_mod_343="
        f"{rho7}"
    )

    print(
        f"  coefficient_6_mod_343="
        f"{6 % 343}"
    )

    if v7_s1 is not None:

        print(
            f"  7-adic_depth="
            f"{v7_s1}"
        )

    # ------------------------------------------------------------------
    # 7. PRIME-POWER VS ORDINARY PRIME FORMULATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PRIME-POWER / PROJECTIVE-EQUIVALENCE STATEMENT")
    print("=" * 78)

    equivalence_checks = []

    for name, coeff in COEFFICIENTS.items():

        s = schur_values[name]

        for ell in relevant_primes:

            if gcd(
                ell,
                Q3,
            ) != 1:

                continue

            for e in range(
                1,
                4,
            ):

                modulus = int(
                    ell ** e
                )

                rho = ratio_mod(
                    Q1,
                    Q3,
                    modulus,
                )

                if rho is None:
                    continue

                left = (
                    s
                    % modulus
                    == 0
                )

                right = (
                    rho
                    == (
                        coeff
                        % modulus
                    )
                )

                exact = (
                    left == right
                )

                equivalence_checks.append(
                    exact
                )

                if not exact:

                    print(
                        f"  FAILURE: "
                        f"{name} ell={ell} e={e}"
                    )

    equivalence_exact = all(
        equivalence_checks
    )

    print(
        f"  checks="
        f"{len(equivalence_checks)}"
    )

    print(
        f"  equivalence_exact="
        f"{equivalence_exact}"
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
For every prime ell with ell not dividing q3,

    s_c = q1 - c*q3

can be rewritten modulo ell^e as

    s_c = 0
    iff
    q1/q3 = c.

Experiment 228 checks this not merely at e=1, but through the complete
valuation depth of every observed Schur prime.

Thus a factor

    ell^v || s_c

is equivalent to a projective congruence of exact depth v:

    q1/q3 = c (mod ell^v),

but

    q1/q3 != c (mod ell^(v+1)).

This is the finite prime analogue of the 7-adic distance calculation.

If the checks pass, the Schur factorization is naturally interpreted
as a collection of local projective approximation events.

The primes 137 and 149239 are therefore not merely factors of s0:
they are primes for which the source ratio is exactly congruent to 2
to one prime-adic digit beyond the trivial level.

Likewise the exceptional factor 7^2 in s1 is a depth-two projective
approximation to the constant 6.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        all_prime_power_checks
        and defect_checks
        and equivalence_exact
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_prime_power_checks_exact="
        f"{all_prime_power_checks}"
    )

    print(
        f"  first_defect_checks_exact="
        f"{defect_checks}"
    )

    print(
        f"  prime_power_equivalence_exact="
        f"{equivalence_exact}"
    )

    print(
        f"  maximum_depth="
        f"{maximum_depth}"
    )

    print(
        f"  maximum_depth_events="
        f"{maximum_records}"
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
    print("EXPERIMENT 228 COMPLETE")


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

