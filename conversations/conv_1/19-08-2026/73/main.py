#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 241 — EXACT THIRD-LEVEL PRIME-LOCAL DERIVATIVE / HENSEL AUDIT
==============================================================================

Experiment 239:
    unique lifts ell -> ell^2 -> ell^3

Experiment 240:
    exact affine residual law and derivative formula through ell^3

Experiment 241:
    continue one more level:

        ell^3 -> ell^4

For each admissible prime ell <= 200:

    k1  solves modulo ell
    k2  solves modulo ell^2
    k3  solves modulo ell^3
    k4  solves modulo ell^4

At the third lift level define

    L3 = ord_{ell^3}(3)

and

    k(t) = k3 + t*L3.

The normalized residual is

    R3(t)
      = F(k(t))/ell^3 mod ell

with

    F(k) = q1 - 2*3^k*q3.

The predicted slope is

    B3 = -q3*rho*D3 mod ell,

where

    D3 = (3^L3 - 1)/ell^3 mod ell.

The experiment verifies:

    * unique first lift;
    * unique second lift;
    * unique third lift;
    * affine residual law at level ell^3;
    * derivative formula at level ell^3;
    * algebraic recovery of the third lift digit;
    * exact k4 reconstruction modulo ell^4.

The order at ell^3 is computed efficiently:

    ord_{ell^3}(3)

is either

    ord_{ell^2}(3)

or

    ell * ord_{ell^2}(3),

so no long brute-force order search is needed.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# SOURCE
# ============================================================================

Q1 = 29144191
Q3 = 24794967

PRIME_LIMIT = 200


# ============================================================================
# BASIC HELPERS
# ============================================================================

def is_prime(n: int) -> bool:
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


def ratio_mod(modulus: int) -> int | None:
    modulus = int(modulus)

    if gcd(Q3, modulus) != 1:
        return None

    return int(
        (
            (Q1 % modulus)
            * pow(
                Q3 % modulus,
                -1,
                modulus,
            )
        )
        % modulus
    )


def orbit_mod(k: int, modulus: int) -> int:
    return int(
        (
            2
            * pow(
                3,
                int(k),
                int(modulus),
            )
        )
        % modulus
    )


def residual_mod(k: int, modulus: int) -> int:
    return int(
        (
            Q1
            -
            (
                2
                * pow(
                    3,
                    int(k),
                    int(modulus),
                )
                * Q3
            )
        )
        % modulus
    )


def multiplicative_order_small(
    base: int,
    prime: int,
) -> int:
    """
    Exact ord_prime(base) by direct iteration.
    Only used modulo a small prime <= 200.
    """

    current = 1

    for e in range(
        1,
        prime + 1,
    ):
        current = (
            current
            * base
        ) % prime

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Could not determine ord_{prime}({base})."
    )


def lift_order_once(
    base: int,
    prime: int,
    previous_order: int,
    previous_power: int,
) -> int:
    """
    Given ord_{p^r}(base) = previous_order, determine ord_{p^(r+1)}.

    It is either previous_order or p*previous_order.
    """

    next_modulus = prime ** (
        previous_power + 1
    )

    if (
        pow(
            base,
            previous_order,
            next_modulus,
        )
        == 1
    ):
        return int(previous_order)

    return int(
        previous_order
        * prime
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 241 — EXACT THIRD-LEVEL PRIME-LOCAL "
        "DERIVATIVE / HENSEL AUDIT"
    )
    print("=" * 78)

    failures = []

    admissible = []

    # ------------------------------------------------------------------
    # 1. FIND ADMISSIBLE CHARTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ADMISSIBLE PRIME CHARTS")
    print("=" * 78)

    for ell in range(
        2,
        PRIME_LIMIT + 1,
    ):

        if not is_prime(ell):
            continue

        if ell in (2, 3):
            continue

        if Q3 % ell == 0:
            continue

        rho = ratio_mod(
            ell
        )

        order1 = multiplicative_order_small(
            3,
            ell,
        )

        base_matches = [
            k
            for k in range(
                order1
            )
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        if len(base_matches) != 1:
            continue

        k1 = base_matches[0]

        # --------------------------------------------------------------
        # ell -> ell^2
        # --------------------------------------------------------------

        mod2 = ell ** 2
        rho2 = ratio_mod(mod2)

        lift1 = [
            (
                t,
                k1 + t * order1,
            )
            for t in range(ell)
            if orbit_mod(
                k1 + t * order1,
                mod2,
            )
            == rho2
        ]

        if len(lift1) != 1:
            failures.append(
                (
                    ell,
                    "lift1",
                    lift1,
                )
            )
            continue

        t1, k2 = lift1[0]

        # --------------------------------------------------------------
        # ord_{ell^2}(3)
        # --------------------------------------------------------------

        order2 = lift_order_once(
            3,
            ell,
            order1,
            1,
        )

        # --------------------------------------------------------------
        # ell^2 -> ell^3
        # --------------------------------------------------------------

        mod3 = ell ** 3
        rho3 = ratio_mod(mod3)

        lift2 = [
            (
                t,
                k2 + t * order2,
            )
            for t in range(ell)
            if orbit_mod(
                k2 + t * order2,
                mod3,
            )
            == rho3
        ]

        if len(lift2) != 1:
            failures.append(
                (
                    ell,
                    "lift2",
                    lift2,
                )
            )
            continue

        t2, k3 = lift2[0]

        # --------------------------------------------------------------
        # ord_{ell^3}(3)
        # --------------------------------------------------------------

        order3 = lift_order_once(
            3,
            ell,
            order2,
            2,
        )

        admissible.append(
            (
                ell,
                rho,
                k1,
                order1,
                t1,
                k2,
                order2,
                t2,
                k3,
                order3,
            )
        )

        print(
            f"  ell={ell}: "
            f"k1={k1} "
            f"t1={t1} "
            f"k2={k2} "
            f"t2={t2} "
            f"k3={k3} "
            f"ord3={order3}"
        )

    print(
        f"  admissible_count={len(admissible)}"
    )

    # ------------------------------------------------------------------
    # 2. THIRD-LEVEL AFFINE RESIDUAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. THIRD-LEVEL AFFINE RESIDUAL")
    print("=" * 78)

    for (
        ell,
        rho,
        k1,
        order1,
        t1,
        k2,
        order2,
        t2,
        k3,
        order3,
    ) in admissible:

        mod3 = ell ** 3
        mod4 = ell ** 4

        rho4 = ratio_mod(
            mod4
        )

        # --------------------------------------------------------------
        # Third-level multiplicative increment
        # --------------------------------------------------------------

        pow3L3 = pow(
            3,
            order3,
            mod4,
        )

        increment = (
            pow3L3 - 1
        )

        if (
            increment % mod3
            != 0
        ):
            failures.append(
                (
                    ell,
                    "D3_not_integral",
                )
            )
            continue

        D3 = int(
            (
                increment
                // mod3
            )
            % ell
        )

        q3_mod = Q3 % ell

        predicted_B3 = int(
            (
                -q3_mod
                * rho
                * D3
            )
            % ell
        )

        # --------------------------------------------------------------
        # All ell possible third-level digits
        # --------------------------------------------------------------

        normalized = []

        division_ok = True

        for t3 in range(ell):

            candidate_k4 = int(
                k3
                + t3 * order3
            )

            F_mod = residual_mod(
                candidate_k4,
                mod4,
            )

            if (
                F_mod % mod3
                != 0
            ):
                division_ok = False
                break

            R3 = int(
                (
                    F_mod
                    // mod3
                )
                % ell
            )

            normalized.append(
                (
                    t3,
                    R3,
                )
            )

        if not division_ok:

            failures.append(
                (
                    ell,
                    "third_normalization",
                )
            )

            continue

        A3 = normalized[0][1]

        slopes = [
            int(
                (
                    normalized[i][1]
                    - normalized[i - 1][1]
                )
                % ell
            )
            for i in range(
                1,
                ell
            )
        ]

        slope_set = sorted(
            set(slopes)
        )

        affine_exact = (
            len(slope_set) == 1
        )

        measured_B3 = (
            slope_set[0]
            if affine_exact
            else None
        )

        slope_formula_exact = (
            affine_exact
            and measured_B3 == predicted_B3
        )

        full_affine_exact = False

        if measured_B3 is not None:

            full_affine_exact = all(
                normalized[t][1]
                ==
                (
                    A3
                    + t * measured_B3
                ) % ell
                for t in range(ell)
            )

        # --------------------------------------------------------------
        # Recover t3
        # --------------------------------------------------------------

        if (
            measured_B3 is not None
            and gcd(
                measured_B3,
                ell,
            ) == 1
        ):

            recovered_t3 = int(
                (
                    (-A3)
                    * pow(
                        measured_B3,
                        -1,
                        ell,
                    )
                )
                % ell
            )

        else:

            recovered_t3 = None

        # --------------------------------------------------------------
        # Direct validation
        # --------------------------------------------------------------

        if recovered_t3 is not None:

            recovered_k4 = int(
                k3
                + recovered_t3 * order3
            )

            recovered_orbit = orbit_mod(
                recovered_k4,
                mod4,
            )

            direct_match = (
                recovered_orbit
                == rho4
            )

            direct_lifts = [
                (
                    t,
                    k3 + t * order3,
                )
                for t in range(ell)
                if orbit_mod(
                    k3 + t * order3,
                    mod4,
                )
                == rho4
            ]

            direct_unique = (
                len(direct_lifts) == 1
            )

            affine_matches_direct = (
                direct_unique
                and
                direct_lifts[0][0]
                == recovered_t3
            )

        else:

            recovered_k4 = None
            direct_match = False
            direct_lifts = []
            direct_unique = False
            affine_matches_direct = False

        prime_ok = (
            division_ok
            and affine_exact
            and slope_formula_exact
            and full_affine_exact
            and recovered_t3 is not None
            and direct_match
            and direct_unique
            and affine_matches_direct
        )

        if not prime_ok:

            failures.append(
                (
                    ell,
                    {
                        "affine_exact": affine_exact,
                        "slope_formula_exact": slope_formula_exact,
                        "full_affine_exact": full_affine_exact,
                        "recovered_t3": recovered_t3,
                        "direct_match": direct_match,
                        "direct_unique": direct_unique,
                    },
                )
            )

        print(
            f"  ell={ell}: "
            f"k3={k3} "
            f"ord3={order3} "
            f"D3={D3} "
            f"A3={A3} "
            f"B_pred={predicted_B3} "
            f"B_meas={measured_B3} "
            f"t3={recovered_t3} "
            f"k4={recovered_k4} "
            f"OK={prime_ok}"
        )

    # ------------------------------------------------------------------
    # 3. EXPLICIT EXAMPLES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXPLICIT THIRD-LEVEL EXAMPLES")
    print("=" * 78)

    for target in (
        5,
        7,
        17,
        73,
    ):

        rows = [
            row
            for row in admissible
            if row[0] == target
        ]

        if not rows:
            continue

        (
            ell,
            rho,
            k1,
            order1,
            t1,
            k2,
            order2,
            t2,
            k3,
            order3,
        ) = rows[0]

        mod4 = ell ** 4
        rho4 = ratio_mod(
            mod4
        )

        lift3 = [
            (
                t,
                k3 + t * order3,
            )
            for t in range(ell)
            if orbit_mod(
                k3 + t * order3,
                mod4,
            )
            == rho4
        ]

        if len(lift3) == 1:

            t3, k4 = lift3[0]

            print(
                f"  ell={ell}: "
                f"k1={k1} "
                f"t1={t1} "
                f"k2={k2} "
                f"t2={t2} "
                f"k3={k3} "
                f"t3={t3} "
                f"k4={k4}"
            )

    # ------------------------------------------------------------------
    # 4. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 240 verified the affine derivative law through

    ell -> ell^2 -> ell^3.

Experiment 241 continues to

    ell^3 -> ell^4.

At this level

    L3 = ord_{ell^3}(3),

    k(t) = k3 + t*L3,

and

    3^L3
      = 1 + ell^3*D3
      (mod ell^4).

Therefore

    R3(t)
      = F(k(t))/ell^3
      (mod ell)

should satisfy

    R3(t)
      = A3 + B3*t
      (mod ell),

with

    B3
      = -q3*rho*D3
      (mod ell).

If the complete admissible sample passes again, then the affine
local derivative law is stable through three successive prime-power
lifts:

    ell,
    ell^2,
    ell^3,
    ell^4.

That would be a substantially stronger finite-prime analogue of the
7-adic lifting mechanism.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(admissible) > 0
        and len(failures) == 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{len(admissible)}"
    )

    print(
        f"  all_third_level_affine_laws_exact="
        f"{final_ok}"
    )

    print(
        f"  all_third_level_derivative_formulas_exact="
        f"{final_ok}"
    )

    print(
        f"  all_third_level_lift_digits_recovered="
        f"{final_ok}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 241 COMPLETE")


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

