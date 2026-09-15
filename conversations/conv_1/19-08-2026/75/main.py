#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 243 — EXACT PRIME-LOCAL PRINCIPAL-COORDINATE DIGIT AUDIT
==============================================================================

Experiment 242 reconstructed, for every admissible prime ell <= 200,

    k1 -> k2 -> k3 -> k4

from the local affine residual digits.

Experiment 243 asks whether those Hensel digits are literally the
canonical base-ell digits of the principal exponent coordinate.

Let

    d = ord_ell(3).

Write

    k_r = k1 + d*m_r.

When the order grows by ell at every prime-power step,

    ord_{ell^r}(3) = d * ell^(r-1),

and the recurrence

    k_(r+1)
      = k_r + t_r * ord_{ell^r}(3)

becomes

    m_(r+1)
      = m_r + t_r * ell^(r-1).

Therefore

    m_E
      = t_1
        + t_2*ell
        + t_3*ell^2
        + ...

Experiment 243 checks this identity directly.

For every admissible prime:

    1. reconstruct k1,k2,k3,k4;
    2. compute the Hensel digits t1,t2,t3;
    3. compute the actual principal coordinates
           m_r=(k_r-k1)/d;
    4. compute ordinary integer base-ell digits of m_r;
    5. compare them with the Hensel digits;
    6. verify the order-growth condition;
    7. independently verify every exponent against rho modulo ell^r.

This is the ordinary-prime analogue of the canonical base-7 digit
reconstruction already obtained.

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


def multiplicative_order_prime(
    base: int,
    prime: int,
) -> int:

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
        f"Could not determine order modulo {prime}."
    )


def lift_order_once(
    base: int,
    prime: int,
    previous_order: int,
    current_power: int,
) -> int:

    next_modulus = prime ** (
        current_power + 1
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
        previous_order * prime
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


def integer_base_digits(
    value: int,
    base: int,
) -> list[int]:

    value = int(value)
    base = int(base)

    if value < 0:
        raise ValueError(
            "Principal coordinate must be nonnegative."
        )

    if value == 0:
        return [0]

    digits = []

    while value > 0:

        digits.append(
            value % base
        )

        value //= base

    return digits


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 243 — EXACT PRIME-LOCAL PRINCIPAL-COORDINATE "
        "DIGIT AUDIT"
    )
    print("=" * 78)

    failures = []
    admissible_count = 0

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
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )

    # ------------------------------------------------------------------
    # PRIME LOOP
    # ------------------------------------------------------------------

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

        d = multiplicative_order_prime(
            3,
            ell,
        )

        # --------------------------------------------------------------
        # Base exponent k1
        # --------------------------------------------------------------

        base_matches = [
            k
            for k in range(d)
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        if len(base_matches) != 1:
            continue

        admissible_count += 1

        k1 = int(
            base_matches[0]
        )

        # --------------------------------------------------------------
        # Reconstruct k2, k3, k4 from affine digits.
        # --------------------------------------------------------------

        exponents = [k1]
        digits = []

        orders = [d]

        k_current = k1
        order_current = d

        local_ok = True

        intercepts = []
        slopes = []

        for level in (
            1,
            2,
            3,
        ):

            modulus_current = ell ** level
            modulus_next = ell ** (
                level + 1
            )

            normalized = []

            for t in range(ell):

                candidate = int(
                    k_current
                    + t * order_current
                )

                F = residual_mod(
                    candidate,
                    modulus_next,
                )

                if (
                    F % modulus_current
                    != 0
                ):

                    local_ok = False
                    break

                R = int(
                    (
                        F
                        // modulus_current
                    )
                    % ell
                )

                normalized.append(R)

            if not local_ok:
                break

            A = normalized[0]

            slope_set = sorted(
                {
                    int(
                        (
                            normalized[i]
                            - normalized[i - 1]
                        )
                        % ell
                    )
                    for i in range(
                        1,
                        ell
                    )
                }
            )

            if len(slope_set) != 1:
                local_ok = False
                break

            B = slope_set[0]

            if gcd(B, ell) != 1:
                local_ok = False
                break

            affine_exact = all(
                normalized[t]
                ==
                (
                    A
                    + t * B
                ) % ell
                for t in range(ell)
            )

            if not affine_exact:
                local_ok = False
                break

            t = int(
                (
                    (-A)
                    * pow(
                        B,
                        -1,
                        ell,
                    )
                )
                % ell
            )

            k_next = int(
                k_current
                + t * order_current
            )

            next_rho = ratio_mod(
                modulus_next
            )

            alignment = (
                orbit_mod(
                    k_next,
                    modulus_next,
                )
                == next_rho
            )

            if not alignment:
                local_ok = False
                break

            intercepts.append(A)
            slopes.append(B)
            digits.append(t)
            exponents.append(k_next)

            # Advance order.
            order_next = lift_order_once(
                3,
                ell,
                order_current,
                level,
            )

            orders.append(
                order_next
            )

            k_current = k_next
            order_current = order_next

        if not local_ok:

            failures.append(
                (
                    ell,
                    "local_reconstruction_failure",
                )
            )

            continue

        # --------------------------------------------------------------
        # Principal coordinates
        # --------------------------------------------------------------

        principal = []

        principal_exact = True

        for k in exponents:

            numerator = (
                k - k1
            )

            if numerator % d != 0:

                principal_exact = False
                break

            principal.append(
                numerator // d
            )

        if not principal_exact:

            failures.append(
                (
                    ell,
                    "principal_coordinate_not_integral",
                )
            )

            continue

        # --------------------------------------------------------------
        # Compare Hensel digits with integer base-ell digits.
        # --------------------------------------------------------------

        digit_matches = []

        for index, m in enumerate(
            principal
        ):

            actual_digits = integer_base_digits(
                m,
                ell,
            )

            expected_digits = digits[:index]

            # For m=0 the canonical representation may be [0], while
            # the empty truncated Hensel list is the natural expansion.
            if index == 0 and m == 0:

                match = (
                    expected_digits == []
                    or expected_digits == [0]
                )

            else:

                # Remove only harmless high zero padding.
                while (
                    len(expected_digits) > 1
                    and
                    expected_digits[-1] == 0
                ):
                    expected_digits.pop()

                while (
                    len(actual_digits) > 1
                    and
                    actual_digits[-1] == 0
                ):
                    actual_digits.pop()

                match = (
                    expected_digits
                    ==
                    actual_digits
                )

            digit_matches.append(
                match
            )

        all_digit_matches = all(
            digit_matches
        )

        # --------------------------------------------------------------
        # Order growth test.
        # --------------------------------------------------------------

        order_growth = [
            orders[i]
            ==
            d * (ell ** i)
            for i in range(
                len(orders)
            )
        ]

        all_order_growth = all(
            order_growth
        )

        # --------------------------------------------------------------
        # Closed-form reconstruction
        # --------------------------------------------------------------

        closed_form_matches = []

        for r in range(
            len(exponents)
        ):

            if r == 0:

                m_closed = 0

            else:

                m_closed = sum(
                    digits[j]
                    * (
                        ell ** j
                    )
                    for j in range(r)
                )

            k_closed = int(
                k1
                + d * m_closed
            )

            closed_form_matches.append(
                k_closed
                ==
                exponents[r]
            )

        all_closed = all(
            closed_form_matches
        )

        # --------------------------------------------------------------
        # Independent orbit checks
        # --------------------------------------------------------------

        orbit_checks = []

        for E, k in enumerate(
            exponents,
            start=1,
        ):

            modulus = ell ** E

            rho_E = ratio_mod(
                modulus
            )

            orbit_E = orbit_mod(
                k,
                modulus,
            )

            orbit_checks.append(
                orbit_E == rho_E
            )

        all_orbit_checks = all(
            orbit_checks
        )

        prime_ok = (
            all_digit_matches
            and
            all_order_growth
            and
            all_closed
            and
            all_orbit_checks
        )

        if not prime_ok:

            failures.append(
                (
                    ell,
                    {
                        "digit_matches": digit_matches,
                        "order_growth": order_growth,
                        "closed": closed_form_matches,
                        "orbit": orbit_checks,
                    },
                )
            )

        print()
        print(
            f"  ell={ell}:"
        )

        print(
            f"    k_sequence={exponents}"
        )

        print(
            f"    hensel_digits={digits}"
        )

        print(
            f"    principal_m={principal}"
        )

        print(
            f"    integer_base_ell_matches="
            f"{digit_matches}"
        )

        print(
            f"    order_growth={orders}"
        )

        print(
            f"    order_growth_expected="
            f"{[d * (ell ** i) for i in range(len(orders))]}"
        )

        print(
            f"    closed_form_matches="
            f"{closed_form_matches}"
        )

        print(
            f"    orbit_checks="
            f"{orbit_checks}"
        )

        print(
            f"    OK={prime_ok}"
        )

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. GLOBAL SUMMARY")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{admissible_count}"
    )

    print(
        f"  local_failures="
        f"{len(failures)}"
    )

    print(
        f"  failures_detail="
        f"{failures}"
    )

    # ------------------------------------------------------------------
    # SPECIAL CASES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SPECIAL PRIME DIGIT TABLE")
    print("=" * 78)

    for target in (
        5,
        7,
        17,
        73,
    ):

        if not is_prime(target):
            continue

        if Q3 % target == 0:
            continue

        rho = ratio_mod(
            target
        )

        d = multiplicative_order_prime(
            3,
            target,
        )

        base = [
            k
            for k in range(d)
            if orbit_mod(k, target) == rho
        ]

        if len(base) != 1:
            continue

        k1 = base[0]

        # Reconstruct four levels.
        k_values = [k1]
        t_values = []
        current_k = k1
        current_order = d

        for level in (
            1,
            2,
            3,
        ):

            modulus_current = target ** level
            modulus_next = target ** (
                level + 1
            )

            normalized = []

            for t in range(target):

                candidate = (
                    current_k
                    + t * current_order
                )

                F = residual_mod(
                    candidate,
                    modulus_next,
                )

                R = (
                    F
                    // modulus_current
                ) % target

                normalized.append(R)

            A = normalized[0]

            B = (
                normalized[1]
                - normalized[0]
            ) % target

            t = (
                (-A)
                * pow(
                    B,
                    -1,
                    target,
                )
            ) % target

            current_k = (
                current_k
                + t * current_order
            )

            current_order = lift_order_once(
                3,
                target,
                current_order,
                level,
            )

            t_values.append(t)
            k_values.append(current_k)

        m_values = [
            (
                k - k1
            ) // d
            for k in k_values
        ]

        print(
            f"  ell={target}: "
            f"k={k_values} "
            f"digits={t_values} "
            f"m={m_values}"
        )

    # ------------------------------------------------------------------
    # STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 242 reconstructed the local exponent digit-by-digit.

Experiment 243 asks whether those digits have a canonical arithmetic
meaning.

Let

    d = ord_ell(3),

and write

    k_r = k1 + d*m_r.

If

    ord_{ell^r}(3)
      = d*ell^(r-1),

then the Hensel recurrence becomes

    m_(r+1)
      = m_r + t_r*ell^(r-1).

Therefore

    m_E
      = t_1
        + t_2*ell
        + t_3*ell^2
        + ...

The Hensel digits are then literally the ordinary base-ell digits of
the principal exponent coordinate m_E.

This is the exact finite-prime analogue of the canonical base-7
principal-coordinate expansion found earlier.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        admissible_count > 0
        and
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{admissible_count}"
    )

    print(
        f"  all_hensel_digits_match_base_ell_digits="
        f"{final_ok}"
    )

    print(
        f"  all_order_growth_relations_exact="
        f"{final_ok}"
    )

    print(
        f"  all_closed_form_reconstructions_exact="
        f"{final_ok}"
    )

    print(
        f"  all_orbit_checks_exact="
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
    print("EXPERIMENT 243 COMPLETE")


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

