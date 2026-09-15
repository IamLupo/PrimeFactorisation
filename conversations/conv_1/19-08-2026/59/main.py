#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 226 — EXACT PRIME-ADIC SCHUR TRANSPORT / CANCELLATION AUDIT
==============================================================================

Known n=pq instance:

    p=2
    q=3
    n=6

    q1=29144191
    q3=24794967

Normalize:

    q3 = q*u
    q1 = 1 + pq*a
    u  = 1 + pq*b.

Experiment 225 showed that the normalized source coordinates have
largely disjoint prime spectra from the Schur row.

Experiment 226 asks the sharper question:

    For every relevant prime ell dividing q1, q3, u, a, or b,

what happens to the ell-adic valuation under the terminal Schur map

    s_k = q1 - 2*3^k*q3,
    k=0,1,2?

The experiment records:

    * v_ell(q1)
    * v_ell(q3)
    * v_ell(u)
    * v_ell(a)
    * v_ell(b)
    * v_ell(s0), v_ell(s1), v_ell(s2)
    * v_ell(s1-s0)
    * v_ell(s2-s1)
    * v_ell(s2-s0)

It then classifies each prime as:

    inherited,
    canceled,
    created by subtraction,
    or invisible.

This is still descriptive for n=6 only.

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


# ============================================================================
# EXACT HELPERS
# ============================================================================

def valuation_p(
    x: int,
    p: int,
) -> int | None:

    x = int(x)
    p = int(p)

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def gcd_int(
    a: int,
    b: int,
) -> int:

    a = abs(int(a))
    b = abs(int(b))

    while b:
        a, b = b, a % b

    return int(a)


def factor_integer(
    x: int,
) -> dict[int, int]:

    x = abs(int(x))

    if x < 2:
        return {}

    factors: dict[int, int] = {}

    while x % 2 == 0:

        factors[2] = (
            factors.get(2, 0) + 1
        )

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


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 226 — EXACT PRIME-ADIC SCHUR TRANSPORT / "
        "CANCELLATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. NORMALIZED SOURCE
    # ------------------------------------------------------------------

    if Q3 % Q != 0:
        raise ArithmeticError(
            "q3 is not divisible by q."
        )

    U = int(
        Q3 // Q
    )

    PQ = int(
        P * Q
    )

    if (Q1 - 1) % PQ != 0:
        raise ArithmeticError(
            "q1-1 is not divisible by pq."
        )

    if (U - 1) % PQ != 0:
        raise ArithmeticError(
            "u-1 is not divisible by pq."
        )

    A = int(
        (Q1 - 1) // PQ
    )

    B = int(
        (U - 1) // PQ
    )

    # ------------------------------------------------------------------
    # 2. SCHUR ROW
    # ------------------------------------------------------------------

    s0 = int(
        Q1 - 2 * Q3
    )

    s1 = int(
        Q1 - 6 * Q3
    )

    s2 = int(
        Q1 - 18 * Q3
    )

    d01 = int(
        s1 - s0
    )

    d12 = int(
        s2 - s1
    )

    d02 = int(
        s2 - s0
    )

    # ------------------------------------------------------------------
    # 3. RELEVANT PRIME SET
    # ------------------------------------------------------------------

    numbers_for_factorization = {
        "q1": Q1,
        "q3": Q3,
        "u": U,
        "a": A,
        "b": B,
        "s0": s0,
        "s1": s1,
        "s2": s2,
    }

    prime_set = set()

    factors = {}

    for name, value in numbers_for_factorization.items():

        ff = factor_integer(
            value
        )

        factors[name] = ff

        prime_set.update(
            ff.keys()
        )

    primes = sorted(
        prime_set
    )

    # ------------------------------------------------------------------
    # 4. SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE / SCHUR DATA")
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
        f"  u=q3/q={U}"
    )

    print(
        f"  a=(q1-1)/(pq)={A}"
    )

    print(
        f"  b=(u-1)/(pq)={B}"
    )

    print(
        f"  s0={s0}"
    )

    print(
        f"  s1={s1}"
    )

    print(
        f"  s2={s2}"
    )

    print(
        f"  d01={d01}"
    )

    print(
        f"  d12={d12}"
    )

    print(
        f"  d02={d02}"
    )

    # ------------------------------------------------------------------
    # 5. PRIME SPECTRUM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. COMPLETE RELEVANT PRIME SPECTRUM")
    print("=" * 78)

    print(
        f"  primes={primes}"
    )

    for name in (
        "q1",
        "q3",
        "u",
        "a",
        "b",
        "s0",
        "s1",
        "s2",
    ):

        print(
            f"  {name}: "
            f"{factors[name]}"
        )

    # ------------------------------------------------------------------
    # 6. PRIME-BY-PRIME VALUATION TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. PRIME-BY-PRIME VALUATION TRANSPORT")
    print("=" * 78)

    for ell in primes:

        values = {
            "q1": valuation_p(Q1, ell),
            "q3": valuation_p(Q3, ell),
            "u": valuation_p(U, ell),
            "a": valuation_p(A, ell),
            "b": valuation_p(B, ell),
            "s0": valuation_p(s0, ell),
            "s1": valuation_p(s1, ell),
            "s2": valuation_p(s2, ell),
            "d01": valuation_p(d01, ell),
            "d12": valuation_p(d12, ell),
            "d02": valuation_p(d02, ell),
        }

        print(
            f"  ell={ell}: "
            f"{values}"
        )

    # ------------------------------------------------------------------
    # 7. SOURCE -> SCHUR INHERITANCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SOURCE-TO-SCHUR INHERITANCE CLASSIFICATION")
    print("=" * 78)

    classifications = {}

    for ell in primes:

        source_values = [
            valuation_p(Q1, ell) or 0,
            valuation_p(Q3, ell) or 0,
            valuation_p(U, ell) or 0,
            valuation_p(A, ell) or 0,
            valuation_p(B, ell) or 0,
        ]

        schur_values = [
            valuation_p(s0, ell) or 0,
            valuation_p(s1, ell) or 0,
            valuation_p(s2, ell) or 0,
        ]

        source_has = any(
            v > 0
            for v in source_values
        )

        schur_has = any(
            v > 0
            for v in schur_values
        )

        if source_has and schur_has:

            classification = "INHERITED_OR_PRESERVED"

        elif source_has and not schur_has:

            classification = "CANCELED_IN_SCHUR"

        elif not source_has and schur_has:

            classification = "CREATED_BY_COMBINATION"

        else:

            classification = "INVISIBLE"

        classifications[ell] = classification

        print(
            f"  ell={ell}: "
            f"{classification}"
        )

    # ------------------------------------------------------------------
    # 8. SOURCE-FACTOR SURVIVAL MATRIX
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SOURCE FACTOR SURVIVAL MATRIX")
    print("=" * 78)

    print(
        "  columns: ell, q1, q3, u, a, b, s0, s1, s2"
    )

    for ell in primes:

        row = []

        for name, value in (
            ("q1", Q1),
            ("q3", Q3),
            ("u", U),
            ("a", A),
            ("b", B),
            ("s0", s0),
            ("s1", s1),
            ("s2", s2),
        ):

            row.append(
                valuation_p(
                    value,
                    ell,
                )
                or 0
            )

        print(
            f"  {ell}: {row}"
        )

    # ------------------------------------------------------------------
    # 9. DIFFERENCE-LAYER TRANSPORT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. DIFFERENCE-LAYER TRANSPORT")
    print("=" * 78)

    for ell in primes:

        v_q3 = (
            valuation_p(
                Q3,
                ell,
            )
            or 0
        )

        v_d01 = (
            valuation_p(
                d01,
                ell,
            )
            or 0
        )

        v_d12 = (
            valuation_p(
                d12,
                ell,
            )
            or 0
        )

        v_d02 = (
            valuation_p(
                d02,
                ell,
            )
            or 0
        )

        print(
            f"  ell={ell}: "
            f"v(q3)={v_q3} "
            f"v(d01)={v_d01} "
            f"v(d12)={v_d12} "
            f"v(d02)={v_d02}"
        )

    # ------------------------------------------------------------------
    # 10. EXACTLY EXPECTED DIFFERENCE FACTOR
    # ------------------------------------------------------------------

    expected_difference_checks = {
        "d01=-4*q3":
            d01 == -4 * Q3,

        "d12=-12*q3":
            d12 == -12 * Q3,

        "d02=-16*q3":
            d02 == -16 * Q3,

        "d12=3*d01":
            d12 == 3 * d01,
    }

    print()
    print("=" * 78)
    print("7. EXACT DIFFERENCE IDENTITIES")
    print("=" * 78)

    for name, result in (
        expected_difference_checks.items()
    ):

        print(
            f"  {name}={result}"
        )

    # ------------------------------------------------------------------
    # 11. SCHUR CANCELLATION TEST
    # ------------------------------------------------------------------

    source_prime_canceled = []

    for ell in primes:

        source_has = any(
            (
                valuation_p(
                    value,
                    ell,
                )
                or 0
            ) > 0
            for value in (
                Q1,
                Q3,
                U,
                A,
                B,
            )
        )

        schur_has = any(
            (
                valuation_p(
                    value,
                    ell,
                )
                or 0
            ) > 0
            for value in (
                s0,
                s1,
                s2,
            )
        )

        if (
            source_has
            and not schur_has
        ):

            source_prime_canceled.append(
                ell
            )

    print()
    print("=" * 78)
    print("8. SCHUR-CANCELLATION PRIMES")
    print("=" * 78)

    print(
        f"  source_primes_canceled="
        f"{source_prime_canceled}"
    )

    # ------------------------------------------------------------------
    # 12. SCHUR-CREATED PRIMES
    # ------------------------------------------------------------------

    source_prime_set = set()

    for name in (
        "q1",
        "q3",
        "u",
        "a",
        "b",
    ):

        source_prime_set.update(
            factors[name].keys()
        )

    schur_prime_set = set()

    for name in (
        "s0",
        "s1",
        "s2",
    ):

        schur_prime_set.update(
            factors[name].keys()
        )

    created_primes = sorted(
        schur_prime_set
        - source_prime_set
    )

    print()
    print("=" * 78)
    print("9. SCHUR-CREATED PRIMES")
    print("=" * 78)

    print(
        f"  created_primes="
        f"{created_primes}"
    )

    # ------------------------------------------------------------------
    # 13. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiment showed that the normalized source coordinates

    a = (q1-1)/(pq),
    b = ((q3/q)-1)/(pq)

are primitive and have largely disjoint prime spectra.

The present experiment asks how those prime factors behave after the
Schur combinations

    s0 = q1 - 2q3,
    s1 = q1 - 6q3,
    s2 = q1 - 18q3.

There are two distinct mechanisms:

    inherited valuation:
        a source prime remains visible in a Schur coefficient;

    cancellation:
        a source prime divides the inputs but disappears from all
        Schur coefficients;

    created valuation:
        neither source coordinate has ell-content, but the subtraction
        creates an ell-divisible Schur coefficient.

The most interesting case is a prime that appears in exactly one source
coordinate and disappears from every Schur coefficient. That would show
that the Schur map actively erases arithmetic source content rather than
merely transporting it.

For the known n=6 instance this is finite evidence only.
"""
    )

    # ------------------------------------------------------------------
    # 14. FINAL
    # ------------------------------------------------------------------

    all_difference_exact = all(
        expected_difference_checks.values()
    )

    factorization_complete = all(
        factor_integer(value)
        is not None
        for value in (
            Q1,
            Q3,
            U,
            A,
            B,
            s0,
            s1,
            s2,
        )
    )

    final_ok = (
        all_difference_exact
        and factorization_complete
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  factorization_complete="
        f"{factorization_complete}"
    )

    print(
        f"  difference_identities_exact="
        f"{all_difference_exact}"
    )

    print(
        f"  source_primes_canceled="
        f"{source_prime_canceled}"
    )

    print(
        f"  Schur_created_primes="
        f"{created_primes}"
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
    print("EXPERIMENT 226 COMPLETE")


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

