#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 230 — EXACT PRIME-LOCAL ORBIT PERIOD / PROJECTIVE CLASS AUDIT
==============================================================================

Known source:

    q1 = 29144191
    q3 = 24794967

Orbit:

    c_k = 2*3^k

Residual:

    R_k = q1 - c_k*q3.

Experiment 229 found repeated prime divisors, e.g.

    ell = 7:
        k = 1,7,13

    ell = 5:
        k = 2,6,10,14

Experiment 230 asks whether these repetitions are explained exactly by
the multiplicative order of 3 modulo ell.

For every prime ell appearing in the factorizations of R_k, k=0,...,KMAX:

    * collect all observed k with ell | R_k;
    * compute the residue class of k modulo ord_ell(3) whenever feasible;
    * test whether all observed hits lie in one exponent class;
    * compare gaps between repeated hits with ord_ell(3);
    * record the maximum ell-adic valuation seen;
    * identify primes with repeated projective alignment;
    * compare prime-power depth against repeated orbit returns.

This is the prime-local analogue of the 7-adic exponent-lift analysis.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

KMAX = 20


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

    return int(v)


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


def multiplicative_order_3(
    ell: int,
    cap: int = 1000000,
) -> int | None:
    """
    Compute ord_ell(3) by direct iteration.

    This is deliberately conservative: for very large prime factors,
    the order may exceed the cap. In that case return None rather than
    performing an expensive factorization of ell-1.
    """

    ell = int(ell)

    if ell <= 2:
        return None

    if ell == 3:
        return None

    value = 3 % ell

    if value == 1:
        return 1

    current = value

    for k in range(
        1,
        min(
            ell,
            cap,
        ) + 1,
    ):

        if current == 1:
            return int(k)

        current = (
            current * 3
        ) % ell

    return None


def projective_ratio_mod(
    ell: int,
) -> int | None:

    ell = int(ell)

    if gcd(Q3, ell) != 1:
        return None

    return int(
        (
            (Q1 % ell)
            * pow(
                Q3 % ell,
                -1,
                ell,
            )
        )
        % ell
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 230 — EXACT PRIME-LOCAL ORBIT PERIOD / "
        "PROJECTIVE CLASS AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE PROJECTIVE DATA")
    print("=" * 78)

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
    # 2. ORBIT RESIDUALS
    # ------------------------------------------------------------------

    residuals = {}

    prime_factorizations = {}

    print()
    print("=" * 78)
    print("2. ORBIT RESIDUALS")
    print("=" * 78)

    for k in range(
        KMAX + 1
    ):

        c = int(
            2 * pow(
                3,
                k,
            )
        )

        residual = int(
            Q1 - c * Q3
        )

        residuals[k] = residual

        ff = factor_integer(
            residual
        )

        prime_factorizations[k] = ff

        print(
            f"  k={k}: "
            f"c_k={c} "
            f"residual={residual}"
        )

        print(
            f"    factors={ff}"
        )

    # ------------------------------------------------------------------
    # 3. COLLECT PRIME EVENTS
    # ------------------------------------------------------------------

    prime_events: dict[int, list[int]] = {}
    prime_depths: dict[int, list[tuple[int,int]]] = {}

    for k, factors in prime_factorizations.items():

        for ell, exponent in factors.items():

            prime_events.setdefault(
                ell,
                [],
            ).append(
                k
            )

            prime_depths.setdefault(
                ell,
                [],
            ).append(
                (
                    k,
                    exponent,
                )
            )

    for ell in prime_events:

        prime_events[ell] = sorted(
            prime_events[ell]
        )

        prime_depths[ell] = sorted(
            prime_depths[ell]
        )

    primes = sorted(
        prime_events
    )

    # ------------------------------------------------------------------
    # 4. PRIME EVENT SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. PRIME EVENT SUMMARY")
    print("=" * 78)

    for ell in primes:

        print(
            f"  ell={ell}: "
            f"k_hits={prime_events[ell]} "
            f"depths={prime_depths[ell]}"
        )

    # ------------------------------------------------------------------
    # 5. MULTIPLICATIVE ORDER / EXPONENT CLASS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. MULTIPLICATIVE ORDER / EXPONENT CLASS")
    print("=" * 78)

    order_results = {}

    for ell in primes:

        if ell in (
            2,
            3,
        ):

            order = None

        else:

            order = multiplicative_order_3(
                ell
            )

        order_results[ell] = order

        rho = projective_ratio_mod(
            ell
        )

        hits = prime_events[ell]

        if order is not None:

            classes = sorted(
                {
                    k % order
                    for k in hits
                }
            )

        else:

            classes = None

        print(
            f"  ell={ell}: "
            f"ord_ell(3)={order} "
            f"rho_mod_ell={rho} "
            f"k_classes={classes}"
        )

    # ------------------------------------------------------------------
    # 6. GAP ANALYSIS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. REPEATED-HIT GAP ANALYSIS")
    print("=" * 78)

    gap_failures = []

    repeated_primes = []

    for ell in primes:

        hits = prime_events[ell]

        if len(hits) < 2:
            continue

        repeated_primes.append(
            ell
        )

        gaps = [
            hits[i + 1] - hits[i]
            for i in range(
                len(hits) - 1
            )
        ]

        order = order_results[ell]

        print(
            f"  ell={ell}: "
            f"hits={hits} "
            f"gaps={gaps} "
            f"order={order}"
        )

        if order is not None:

            exact_gap = all(
                gap == order
                for gap in gaps
            )

            print(
                f"    gaps_equal_order="
                f"{exact_gap}"
            )

            if not exact_gap:

                gap_failures.append(
                    (
                        ell,
                        hits,
                        gaps,
                        order,
                    )
                )

    # ------------------------------------------------------------------
    # 7. SINGLE PROJECTIVE EXPONENT CLASS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SINGLE PROJECTIVE EXPONENT CLASS TEST")
    print("=" * 78)

    class_failures = []

    for ell in repeated_primes:

        order = order_results[ell]

        if order is None:
            continue

        classes = sorted(
            {
                k % order
                for k in prime_events[ell]
            }
        )

        single_class = (
            len(classes) == 1
        )

        print(
            f"  ell={ell}: "
            f"order={order} "
            f"classes={classes} "
            f"single_class={single_class}"
        )

        if not single_class:

            class_failures.append(
                (
                    ell,
                    classes,
                )
            )

    # ------------------------------------------------------------------
    # 8. PROJECTIVE RATIO MATCH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PROJECTIVE RATIO MATCH AT REPEATED PRIMES")
    print("=" * 78)

    ratio_failures = []

    for ell in repeated_primes:

        rho = projective_ratio_mod(
            ell
        )

        if rho is None:
            continue

        hits = prime_events[ell]

        print(
            f"  ell={ell}: "
            f"rho={rho}"
        )

        for k in hits:

            c = int(
                2 * pow(
                    3,
                    k,
                    ell,
                )
            )

            match = (
                c == rho
            )

            print(
                f"    k={k}: "
                f"c_k_mod_ell={c} "
                f"match={match}"
            )

            if not match:

                ratio_failures.append(
                    (
                        ell,
                        k,
                        c,
                        rho,
                    )
                )

    # ------------------------------------------------------------------
    # 9. PRIME-POWER DEPTH REPEAT PATTERNS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. REPEATED PRIME / DEPTH PATTERN")
    print("=" * 78)

    for ell in repeated_primes:

        print(
            f"  ell={ell}:"
        )

        for k, depth in prime_depths[ell]:

            print(
                f"    k={k}: "
                f"depth={depth}"
            )

    # ------------------------------------------------------------------
    # 10. SPECIAL PRIMES 5 AND 7
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. SPECIAL LOCAL ORBITS: ell=5 AND ell=7")
    print("=" * 78)

    for ell in (
        5,
        7,
    ):

        if ell not in prime_events:
            continue

        order = order_results[ell]
        rho = projective_ratio_mod(
            ell
        )

        print(
            f"  ell={ell}: "
            f"order={order} "
            f"rho={rho} "
            f"hits={prime_events[ell]}"
        )

        if order is not None:

            print(
                f"    hit_classes="
                f"{sorted({k % order for k in prime_events[ell]})}"
            )

    # ------------------------------------------------------------------
    # 11. GLOBAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 229 showed that the same prime can divide several orbit
residuals.

Experiment 230 asks whether those repeated appearances are simply the
periodic return structure of

    c_k = 2*3^k mod ell.

If

    ell | q1 - c_k q3

and ell does not divide q3, then

    c_k = q1/q3 mod ell.

Hence all exponent hits for a fixed ell should lie in the discrete-log
class determined by the source ratio.

When ord_ell(3) is known, repeated hits should recur with that period.

The two especially informative examples here are:

    ell=7:
        k=1,7,13,...

    ell=5:
        k=2,6,10,14,...

The different periods reflect different multiplicative orders.

Prime-power depth is then a second layer:
the same orbit class can return modulo ell while occasionally giving a
higher-order coincidence modulo ell^2.

This is the exact finite-prime analogue of the 7-adic lifting picture.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    all_ratio_exact = (
        len(ratio_failures) == 0
    )

    all_class_exact = (
        len(class_failures) == 0
    )

    all_gap_exact = (
        len(gap_failures) == 0
    )

    repeated_case_exists = (
        len(repeated_primes) > 0
    )

    final_ok = (
        all_ratio_exact
        and all_class_exact
        and all_gap_exact
        and repeated_case_exists
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  repeated_primes="
        f"{repeated_primes}"
    )

    print(
        f"  all_repeated_hits_match_ratio="
        f"{all_ratio_exact}"
    )

    print(
        f"  all_repeated_hits_single_exponent_class="
        f"{all_class_exact}"
    )

    print(
        f"  all_repeated_gaps_equal_order="
        f"{all_gap_exact}"
    )

    print(
        f"  special_ell7_hits="
        f"{prime_events.get(7, [])}"
    )

    print(
        f"  special_ell5_hits="
        f"{prime_events.get(5, [])}"
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
    print("EXPERIMENT 230 COMPLETE")


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

