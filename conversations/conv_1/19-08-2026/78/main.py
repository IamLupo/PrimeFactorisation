#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 246 — EXACT EXCEPTIONAL-PRIME SYNTHETIC ORBIT / DELAYED HENSEL AUDIT
==============================================================================

Experiment 245 found:

    ord_11(3)   = 5
    ord_121(3)  = 5
    ord_1331(3) = 55
    ord_14641(3)= 605

So ell=11 is exceptional for order growth, but the original source ratio
rho=q1/q3 satisfies

    rho = 1 (mod 11),

which is outside the orbit

    2*<3> mod 11.

Therefore Experiment 245 could not test the exceptional Hensel mechanism
inside an actual admissible chart.

Experiment 246 removes that ambiguity by constructing the five admissible
projective orbit classes synthetically.

For ell=11 define

    rho_k = 2*3^k mod 11^4,

for k=0,...,4.

Use the synthetic source pair

    q3 = 1,
    q1 = rho_k.

Then the source ratio is exactly rho_k.

For each of the five classes, test:

    mod 11
    mod 11^2
    mod 11^3
    mod 11^4.

The experiment determines:

    * whether the base class lifts to 11^2;
    * whether a digit is visible at the exceptional 11 -> 121 step;
    * how the exponent period changes at 1331;
    * whether the first nontrivial exponent correction occurs at the
      delayed order-growth level;
    * whether the affine residual law can be formulated with the ACTUAL
      local order rather than assuming order multiplication by 11.

This is a diagnostic experiment for the exceptional branch.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


ELL = 11
MAX_LEVEL = 4


# ============================================================================
# HELPERS
# ============================================================================

def multiplicative_order(
    base: int,
    modulus: int,
) -> int:

    if gcd(base, modulus) != 1:
        raise ArithmeticError(
            f"{base} is not invertible modulo {modulus}."
        )

    current = 1

    for e in range(
        1,
        modulus + 1,
    ):

        current = (
            current * base
        ) % modulus

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Order not found for {base} modulo {modulus}."
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
                int(modulus),
            )
        )
        % modulus
    )


def synthetic_residual(
    q1: int,
    q3: int,
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            q1
            -
            (
                2
                * pow(
                    3,
                    int(k),
                    modulus,
                )
                * q3
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
        "EXPERIMENT 246 — EXACT EXCEPTIONAL-PRIME SYNTHETIC "
        "ORBIT / DELAYED HENSEL AUDIT"
    )
    print("=" * 78)

    failures = []

    print()
    print("=" * 78)
    print("1. ORDER PROFILE")
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

        print(
            f"  level={level}: "
            f"modulus={modulus} "
            f"ord={order}"
        )

    # ------------------------------------------------------------------
    # Synthetic admissible projective classes.
    # ------------------------------------------------------------------

    base_order = orders[1]

    base_classes = [
        orbit_mod(
            k,
            ELL,
        )
        for k in range(
            base_order
        )
    ]

    print()
    print("=" * 78)
    print("2. SYNTHETIC ADMISSIBLE CLASSES")
    print("=" * 78)

    print(
        f"  base_order={base_order}"
    )

    print(
        f"  rho_mod_11_classes={base_classes}"
    )

    # ------------------------------------------------------------------
    # Each class is represented by the exact orbit value modulo 11^4.
    # q3=1, q1=rho.
    # ------------------------------------------------------------------

    for base_k in range(
        base_order
    ):

        q3 = 1

        full_modulus = ELL ** MAX_LEVEL

        q1 = orbit_mod(
            base_k,
            full_modulus,
        )

        print()
        print("=" * 78)
        print(
            f"3. EXCEPTIONAL CHART base_k={base_k}"
        )
        print("=" * 78)

        print(
            f"  synthetic_q1={q1}"
        )

        print(
            f"  synthetic_q3={q3}"
        )

        print(
            f"  gcd={gcd(q1,q3)}"
        )

        # --------------------------------------------------------------
        # Verify source ratio at each precision.
        # --------------------------------------------------------------

        print(
            "  orbit_alignment:"
        )

        alignment = []

        for level in range(
            1,
            MAX_LEVEL + 1,
        ):

            modulus = ELL ** level

            expected = q1 % modulus

            actual = orbit_mod(
                base_k,
                modulus,
            )

            ok = (
                actual == expected
            )

            alignment.append(ok)

            print(
                f"    level={level}: "
                f"modulus={modulus} "
                f"orbit={actual} "
                f"rho={expected} "
                f"exact={ok}"
            )

        if not all(alignment):

            failures.append(
                (
                    base_k,
                    "base_orbit_alignment",
                )
            )

        # --------------------------------------------------------------
        # Analyze each lifting transition using the actual order.
        # --------------------------------------------------------------

        current_k = base_k

        exponent_sequence = [
            current_k
        ]

        print(
            "  actual-order lift transitions:"
        )

        for level in range(
            1,
            MAX_LEVEL
        ):

            current_order = orders[level]

            next_modulus = ELL ** (
                level + 1
            )

            target_rho = q1 % next_modulus

            matches = [
                (
                    t,
                    current_k
                    + t * current_order,
                )
                for t in range(
                    ELL
                )
                if orbit_mod(
                    current_k
                    + t * current_order,
                    next_modulus,
                )
                == target_rho
            ]

            print(
                f"    level={level}: "
                f"current_k={current_k} "
                f"current_order={current_order} "
                f"next_modulus={next_modulus} "
                f"matches={matches}"
            )

            if len(matches) == 0:

                print(
                    "      lift_status=NO_LIFT_IN_STANDARD_11-CANDIDATE_FAMILY"
                )

                continue

            if len(matches) > 1:

                print(
                    "      lift_status=MULTIPLE_LIFTS"
                )

            else:

                t, next_k = matches[0]

                print(
                    f"      lift_status=UNIQUE "
                    f"t={t} "
                    f"next_k={next_k}"
                )

                exponent_sequence.append(
                    next_k
                )

                current_k = next_k

        print(
            f"  exponent_sequence={exponent_sequence}"
        )

        # --------------------------------------------------------------
        # Compare naive and actual-order candidate families.
        # --------------------------------------------------------------

        print(
            "  candidate-family comparison:"
        )

        for level in range(
            1,
            MAX_LEVEL
        ):

            actual_order = orders[level]

            naive_order = (
                orders[1]
                * ELL ** (
                    level - 1
                )
            )

            actual_candidates = [
                current_k
                + t * actual_order
                for t in range(ELL)
            ]

            naive_candidates = [
                current_k
                + t * naive_order
                for t in range(ELL)
            ]

            print(
                f"    level={level}: "
                f"actual_order={actual_order} "
                f"naive_order={naive_order}"
            )

            if actual_order != naive_order:

                print(
                    "      order_exception=True"
                )

                print(
                    f"      actual_candidates="
                    f"{actual_candidates}"
                )

                print(
                    f"      naive_candidates="
                    f"{naive_candidates}"
                )

        # --------------------------------------------------------------
        # Residual valuation at the known exact orbit exponent.
        # --------------------------------------------------------------

        print(
            "  residual_valuation_profile:"
        )

        for E, k in enumerate(
            exponent_sequence,
            start=1,
        ):

            modulus = ELL ** E

            Fmod = synthetic_residual(
                q1,
                q3,
                k,
                modulus,
            )

            if Fmod == 0:

                detected = E

            else:

                detected = 0
                x = Fmod

                while (
                    detected < E
                    and x % ELL == 0
                ):
                    x //= ELL
                    detected += 1

            print(
                f"    E={E}: "
                f"k={k} "
                f"F_mod={Fmod} "
                f"v11_detected={detected}"
            )

    # ------------------------------------------------------------------
    # 4. ORDER-DELAY INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXCEPTIONAL-BRANCH INTERPRETATION")
    print("=" * 78)

    print(
        """
For ell=11 the generic order law fails at the first step:

    ord_11(3) = 5,
    ord_121(3) = 5.

Therefore an ordinary digit t_1 in

    k2 = k1 + t_1*ord_11(3)

cannot create eleven distinct classes modulo 121.

The first genuine order enlargement occurs only at

    ord_1331(3) = 55.

Experiment 246 therefore tests whether the exceptional orbit should be
described by a delayed digit system:

    level 1:
        no 11-fold exponent splitting;

    level 2:
        still the same order;

    level 3:
        the order expands by 11 and a new exponent digit becomes
        available.

This determines whether the final theorem needs:

    (a) a nonexceptional hypothesis,

or

    (b) an explicit delayed-growth exceptional branch.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exceptional_prime={ELL}"
    )

    print(
        f"  order_profile={orders}"
    )

    print(
        f"  synthetic_classes_tested={base_order}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        f"  failure_details={failures}"
    )

    print(
        "  generic_nonexceptional_formula_not_assumed=True"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  ALL BASIC CHECKS PASS={len(failures) == 0}"
    )

    print()
    print("EXPERIMENT 246 COMPLETE")


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

