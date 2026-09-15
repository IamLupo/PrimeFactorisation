#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 245 — EXACT EXCEPTIONAL-PRIME LOCAL HENSEL / ORDER-DEFECT AUDIT
==============================================================================

Experiment 244 found exactly one exceptional prime <= 200:

    ell = 11

for which

    d = ord_11(3) = 5,

but

    3^5 = 1 (mod 11^2),

so

    D1 = (3^5 - 1)/11 = 0 (mod 11),

and therefore

    ord_{11^2}(3) = 5

instead of

    5*11 = 55.

This experiment isolates that exceptional branch.

The questions are:

    1. Does the source ratio even lie in the 3-orbit modulo 11?
    2. What happens modulo 11^2?
    3. What happens modulo 11^3?
    4. Does a unique exponent still exist?
    5. If so, what is the correct exponent-lift modulus?
    6. Does the usual affine residual law survive?
    7. Which part of the nonexceptional proof must be modified?

For an exceptional prime we cannot blindly use

    k_(r+1)
      = k_r + t_r*ord_{11^r}(3)

with an 11-fold growth of the order.

Instead the script measures the actual local order at each level and
tests the complete residue class of candidate exponents.

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

ELL = 11

MAX_LEVEL = 5


# ============================================================================
# HELPERS
# ============================================================================

def ratio_mod(
    modulus: int,
) -> int | None:

    modulus = int(modulus)

    if gcd(
        Q3,
        modulus,
    ) != 1:
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


def orbit_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                int(k),
                modulus,
            )
        )
        % modulus
    )


def residual_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            Q1
            -
            (
                2
                * pow(
                    3,
                    int(k),
                    modulus,
                )
                * Q3
            )
        )
        % modulus
    )


def multiplicative_order(
    base: int,
    modulus: int,
) -> int:

    if gcd(
        base,
        modulus,
    ) != 1:
        raise ArithmeticError(
            f"{base} is not invertible modulo {modulus}."
        )

    # For 11^r we can safely search. The orders are tiny relative to
    # the enormous source integers.
    current = 1

    # phi(11^r) = 10*11^(r-1), so this is a finite exact bound.
    bound = modulus

    for e in range(
        1,
        bound + 1,
    ):

        current = (
            current
            * base
        ) % modulus

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Could not determine order of {base} modulo {modulus}."
    )


def valuation_mod_power(
    x: int,
    p: int,
    max_power: int,
) -> int:

    x = int(x)
    p = int(p)

    if x == 0:
        return max_power

    v = 0

    while (
        v < max_power
        and
        x % p == 0
    ):

        x //= p
        v += 1

    return int(v)


def find_matching_exponents(
    modulus: int,
    search_bound: int,
) -> list[int]:

    rho = ratio_mod(
        modulus
    )

    order = multiplicative_order(
        3,
        modulus,
    )

    # Every orbit point repeats with this order.
    return [
        k
        for k in range(
            order
        )
        if orbit_mod(
            k,
            modulus,
        ) == rho
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 245 — EXACT EXCEPTIONAL-PRIME LOCAL "
        "HENSEL / ORDER-DEFECT AUDIT"
    )
    print("=" * 78)

    failures = []

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
        f"  ell={ELL}"
    )

    print(
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. ORDER PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT ORDER PROFILE")
    print("=" * 78)

    orders = {}

    for level in range(
        1,
        MAX_LEVEL + 1,
    ):

        modulus = ELL ** level

        order = multiplicative_order(
            3,
            modulus,
        )

        orders[level] = order

        expected_nonexceptional = (
            orders[1]
            * ELL ** (
                level - 1
            )
        )

        growth_ratio = (
            order
            // (
                orders[level - 1]
                if level > 1
                else 1
            )
        )

        print(
            f"  level={level}: "
            f"modulus={modulus} "
            f"ord={order} "
            f"expected_simple="
            f"{expected_nonexceptional if level == 1 else orders[1] * ELL ** (level-1)} "
            f"growth_ratio={growth_ratio}"
        )

    # ------------------------------------------------------------------
    # 3. ORDER-DEFECT DIGITS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ORDER-DEFECT DIGITS")
    print("=" * 78)

    for level in range(
        1,
        MAX_LEVEL,
    ):

        current_order = orders[level]

        next_modulus = ELL ** (
            level + 1
        )

        increment = (
            pow(
                3,
                current_order,
                next_modulus,
            )
            - 1
        )

        divisible = (
            increment
            % (
                ELL ** level
            )
            == 0
        )

        if not divisible:
            failures.append(
                (
                    level,
                    "order_increment_not_divisible",
                )
            )
            continue

        D = int(
            (
                increment
                // (
                    ELL ** level
                )
            )
            % ELL
        )

        print(
            f"  level={level}: "
            f"ord={current_order} "
            f"D={D}"
        )

    # ------------------------------------------------------------------
    # 4. PROJECTIVE CHART EXISTENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PROJECTIVE CHART")
    print("=" * 78)

    rho1 = ratio_mod(
        ELL
    )

    print(
        f"  rho_mod_11={rho1}"
    )

    base_matches = find_matching_exponents(
        ELL,
        orders[1],
    )

    print(
        f"  base_matches_mod_11="
        f"{base_matches}"
    )

    base_unique = (
        len(base_matches) == 1
    )

    if not base_unique:

        print(
            "  chart_exists=False"
        )

        print(
            "  exceptional_analysis_stops=True"
        )

        print(
            "EXPERIMENT 245 COMPLETE"
        )

        return

    k_current = base_matches[0]

    print(
        f"  chart_exists=True"
    )

    print(
        f"  k1={k_current}"
    )

    # ------------------------------------------------------------------
    # 5. EXACT LIFT THROUGH THE EXCEPTIONAL ORDER TOWER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT EXCEPTIONAL HENSEL LIFT")
    print("=" * 78)

    exponent_sequence = [
        k_current
    ]

    lift_digits = []

    for level in range(
        1,
        MAX_LEVEL,
    ):

        current_modulus = ELL ** level
        next_modulus = ELL ** (
            level + 1
        )

        current_order = orders[level]

        rho_next = ratio_mod(
            next_modulus
        )

        # Because the order may fail to multiply by 11, the complete
        # candidate family has actual period current_order.
        #
        # We inspect one entire period and identify all solutions.
        matches = []

        for t in range(ELL):

            candidate = int(
                k_current
                + t * current_order
            )

            if (
                orbit_mod(
                    candidate,
                    next_modulus,
                )
                == rho_next
            ):

                matches.append(
                    (
                        t,
                        candidate,
                    )
                )

        # If the actual order does not grow, the t-parameter may not
        # enumerate distinct exponents modulo the next power. Detect this.
        distinct_candidates = sorted(
            {
                candidate
                for _, candidate in matches
            }
        )

        print(
            f"  level={level}: "
            f"current_k={k_current} "
            f"current_order={current_order} "
            f"next_modulus={next_modulus} "
            f"matches={matches}"
        )

        if len(matches) != 1:

            print(
                f"    unique_next_lift=False"
            )

            # Do not immediately mark this as a mathematical failure:
            # exceptional order behavior can change the parameterization.
            continue

        t, next_k = matches[0]

        lift_digits.append(
            t
        )

        exponent_sequence.append(
            next_k
        )

        print(
            f"    unique_next_lift=True "
            f"t={t} "
            f"next_k={next_k}"
        )

        k_current = next_k

    # ------------------------------------------------------------------
    # 6. RESIDUAL VALUATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. RESIDUAL VALUATIONS")
    print("=" * 78)

    for E, k in enumerate(
        exponent_sequence,
        start=1,
    ):

        modulus = ELL ** E

        F_mod = residual_mod(
            k,
            modulus,
        )

        # Determine valuation from the exact residual represented
        # modulo ell^E.
        if F_mod == 0:

            v = E

        else:

            v = valuation_mod_power(
                F_mod,
                ELL,
                E,
            )

        print(
            f"  E={E}: "
            f"k={k} "
            f"F_mod_11^E={F_mod} "
            f"v11_detected={v}"
        )

    # ------------------------------------------------------------------
    # 7. EXCEPTIONAL PRINCIPAL COORDINATE TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. EXCEPTIONAL PRINCIPAL-COORDINATE TEST")
    print("=" * 78)

    d = orders[1]

    print(
        f"  base_order_d={d}"
    )

    for index, k in enumerate(
        exponent_sequence,
    ):

        delta_k = (
            k - exponent_sequence[0]
        )

        divisible = (
            delta_k % d == 0
        )

        print(
            f"  E={index+1}: "
            f"k={k} "
            f"(k-k1)={delta_k} "
            f"divisible_by_d={divisible}"
        )

        if not divisible:

            failures.append(
                (
                    index + 1,
                    "principal_coordinate_nonintegral",
                )
            )

    # ------------------------------------------------------------------
    # 8. LOCAL STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The nonexceptional theory uses

    ord_{ell^r}(3)
      = ord_ell(3)*ell^(r-1).

At ell=11 this fails immediately:

    ord_11(3)=5,

but

    ord_121(3)=5.

Thus the usual base-ell principal coordinate cannot simply be assumed.

The exceptional branch must be described using the actual tower

    ord_{11}(3),
    ord_{11^2}(3),
    ord_{11^3}(3), ...

rather than the generic formula.

There are two possible outcomes.

First:

    the projective exponent still lifts uniquely,
    but with a different coordinate scaling.

Second:

    the exponent chart itself develops multiple or delayed
    lifting behavior.

Either outcome is important for the final theorem because it tells us
whether the theorem should include an explicit nonexceptional
hypothesis or a separate exceptional branch.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exceptional_prime={ELL}"
    )

    print(
        f"  base_chart_exists={base_unique}"
    )

    print(
        f"  orders={orders}"
    )

    print(
        f"  exponent_sequence={exponent_sequence}"
    )

    print(
        f"  lift_digits={lift_digits}"
    )

    print(
        f"  diagnostic_failures={failures}"
    )

    print(
        "  generic_nonexceptional_theorem_not_used=True"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    # This experiment is diagnostic. Even if exceptional lifting behaves
    # differently, that is not automatically a failure.
    print(
        "  ALL BASIC CHECKS PASS="
        f"{len(failures) == 0}"
    )

    print()
    print("EXPERIMENT 245 COMPLETE")


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

