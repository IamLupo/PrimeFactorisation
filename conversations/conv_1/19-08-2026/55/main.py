#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 222 — EXACT TERMINAL SOURCE CONGRUENCE /
                 NORMALIZED PRIME-UNIT AUDIT
==============================================================================

Known instance:

    n = p*q = 2*3 = 6

    q1 = 29144191
    q3 = 24794967

Experiment 221R showed:

    q1 == 1 (mod p),
    q1 == 1 (mod q),

    q3 == 1 (mod p),
    q3 == 0 (mod q),

and

    q3/q = 8264989.

Experiment 222 isolates this signature.

The purpose is to determine exactly:

    * q-adic order of q3;
    * p-adic unit residues;
    * normalized source pair;
    * source behavior modulo p*q;
    * source behavior modulo powers of q;
    * source behavior modulo powers of p;
    * whether simple first-order congruence patterns are visible.

This does NOT claim a general p,q theorem.
It identifies the exact local source signature that must be tested
on the next independent n=pq instance.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN n=pq CASE
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
        a, b = (
            b,
            a % b,
        )

    return a


def residues(
    x: int,
    p: int,
    levels: int,
) -> list[tuple[int, int]]:

    out = []

    for e in range(
        1,
        levels + 1,
    ):

        modulus = int(
            p ** e
        )

        out.append(
            (
                modulus,
                int(x % modulus),
            )
        )

    return out


def exact_division(
    x: int,
    d: int,
) -> tuple[bool, int | None]:

    if d == 0:
        return (
            False,
            None,
        )

    if x % d != 0:
        return (
            False,
            None,
        )

    return (
        True,
        int(x // d),
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 222 — EXACT TERMINAL SOURCE CONGRUENCE / "
        "NORMALIZED PRIME-UNIT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. PARAMETERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. n=pq PARAMETERS")
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

    # ------------------------------------------------------------------
    # 2. p-ADIC / q-ADIC ORDERS
    # ------------------------------------------------------------------

    v_p_q1 = valuation_p(
        Q1,
        P,
    )

    v_p_q3 = valuation_p(
        Q3,
        P,
    )

    v_q_q1 = valuation_p(
        Q1,
        Q,
    )

    v_q_q3 = valuation_p(
        Q3,
        Q,
    )

    print()
    print("=" * 78)
    print("2. EXACT PRIME-ADIC ORDERS")
    print("=" * 78)

    print(
        f"  v_p(q1)={v_p_q1}"
    )

    print(
        f"  v_p(q3)={v_p_q3}"
    )

    print(
        f"  v_q(q1)={v_q_q1}"
    )

    print(
        f"  v_q(q3)={v_q_q3}"
    )

    # ------------------------------------------------------------------
    # 3. FIRST CONGRUENCE SIGNATURE
    # ------------------------------------------------------------------

    q1_mod_p = Q1 % P
    q3_mod_p = Q3 % P

    q1_mod_q = Q1 % Q
    q3_mod_q = Q3 % Q

    print()
    print("=" * 78)
    print("3. FIRST CONGRUENCE SIGNATURE")
    print("=" * 78)

    print(
        f"  q1 mod p={q1_mod_p}"
    )

    print(
        f"  q3 mod p={q3_mod_p}"
    )

    print(
        f"  q1 mod q={q1_mod_q}"
    )

    print(
        f"  q3 mod q={q3_mod_q}"
    )

    print(
        f"  q1_equals_1_mod_p="
        f"{q1_mod_p == 1 % P}"
    )

    print(
        f"  q1_equals_1_mod_q="
        f"{q1_mod_q == 1 % Q}"
    )

    print(
        f"  q3_equals_1_mod_p="
        f"{q3_mod_p == 1 % P}"
    )

    print(
        f"  q3_equals_0_mod_q="
        f"{q3_mod_q == 0}"
    )

    # ------------------------------------------------------------------
    # 4. NORMALIZED q-SOURCE
    # ------------------------------------------------------------------

    q3_divisible_by_q, q3_unit_q = exact_division(
        Q3,
        Q,
    )

    print()
    print("=" * 78)
    print("4. q-NORMALIZED SOURCE")
    print("=" * 78)

    print(
        f"  q3_divisible_by_q="
        f"{q3_divisible_by_q}"
    )

    print(
        f"  q3/q={q3_unit_q}"
    )

    if q3_unit_q is not None:

        print(
            f"  (q3/q) mod p="
            f"{q3_unit_q % P}"
        )

        print(
            f"  (q3/q) mod q="
            f"{q3_unit_q % Q}"
        )

        print(
            f"  v_q(q3/q)="
            f"{valuation_p(q3_unit_q,Q)}"
        )

        print(
            f"  v_p(q3/q)="
            f"{valuation_p(q3_unit_q,P)}"
        )

    # ------------------------------------------------------------------
    # 5. p-NORMALIZED CONTROLS
    # ------------------------------------------------------------------

    q1_minus_1 = (
        Q1 - 1
    )

    q3_minus_1 = (
        Q3 - 1
    )

    p_div_q1_minus_1, q1_p_quotient = exact_division(
        q1_minus_1,
        P,
    )

    p_div_q3_minus_1, q3_p_quotient = exact_division(
        q3_minus_1,
        P,
    )

    print()
    print("=" * 78)
    print("5. p-NORMALIZED FIRST-ORDER CONTROLS")
    print("=" * 78)

    print(
        f"  q1-1={q1_minus_1}"
    )

    print(
        f"  q3-1={q3_minus_1}"
    )

    print(
        f"  p | (q1-1)={p_div_q1_minus_1}"
    )

    print(
        f"  p | (q3-1)={p_div_q3_minus_1}"
    )

    print(
        f"  (q1-1)/p={q1_p_quotient}"
    )

    print(
        f"  (q3-1)/p={q3_p_quotient}"
    )

    # ------------------------------------------------------------------
    # 6. HIGHER q-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. HIGHER q-ADIC SOURCE PROFILE")
    print("=" * 78)

    print(
        "  q1:"
    )

    for modulus, value in residues(
        Q1,
        Q,
        8,
    ):

        print(
            f"    mod {modulus}: {value}"
        )

    print(
        "  q3:"
    )

    for modulus, value in residues(
        Q3,
        Q,
        8,
    ):

        print(
            f"    mod {modulus}: {value}"
        )

    # ------------------------------------------------------------------
    # 7. HIGHER p-ADIC PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. HIGHER p-ADIC SOURCE PROFILE")
    print("=" * 78)

    print(
        "  q1:"
    )

    for modulus, value in residues(
        Q1,
        P,
        8,
    ):

        print(
            f"    mod {modulus}: {value}"
        )

    print(
        "  q3:"
    )

    for modulus, value in residues(
        Q3,
        P,
        8,
    ):

        print(
            f"    mod {modulus}: {value}"
        )

    # ------------------------------------------------------------------
    # 8. COMPOSITE MODULUS PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. MODULO n=pq PROFILE")
    print("=" * 78)

    print(
        f"  q1 mod n={Q1 % N}"
    )

    print(
        f"  q3 mod n={Q3 % N}"
    )

    print(
        f"  q1-q3 mod n="
        f"{(Q1-Q3) % N}"
    )

    print(
        f"  q1+q3 mod n="
        f"{(Q1+Q3) % N}"
    )

    # ------------------------------------------------------------------
    # 9. SOURCE DIFFERENCE STRUCTURE
    # ------------------------------------------------------------------

    delta = (
        Q1 - Q3
    )

    print()
    print("=" * 78)
    print("9. SOURCE DIFFERENCE STRUCTURE")
    print("=" * 78)

    print(
        f"  q1-q3={delta}"
    )

    print(
        f"  v_p(q1-q3)="
        f"{valuation_p(delta,P)}"
    )

    print(
        f"  v_q(q1-q3)="
        f"{valuation_p(delta,Q)}"
    )

    print(
        f"  gcd(q1-q3,q1)="
        f"{gcd_int(delta,Q1)}"
    )

    print(
        f"  gcd(q1-q3,q3)="
        f"{gcd_int(delta,Q3)}"
    )

    # ------------------------------------------------------------------
    # 10. PRIME-PARAMETER RECONSTRUCTION CANDIDATES
    # ------------------------------------------------------------------
    #
    # These are intentionally simple candidate forms. A True result is
    # interesting, but one False result does not disprove a more complex
    # source formula.
    # ------------------------------------------------------------------

    candidates = {
        "q3 = q * unit":
            q3_divisible_by_q,

        "q1 = 1 mod q":
            q1_mod_q == 1,

        "q1 = 1 mod p":
            q1_mod_p == 1 % P,

        "q3 = 1 mod p":
            q3_mod_p == 1 % P,

        "q1 and q3 both p-units":
            (
                v_p_q1 == 0
                and v_p_q3 == 0
            ),

        "q1 q-unit":
            v_q_q1 == 0,

        "q3 exact q-order 1":
            v_q_q3 == 1,
    }

    print()
    print("=" * 78)
    print("10. CANDIDATE SOURCE SIGNATURE")
    print("=" * 78)

    for name, value in candidates.items():

        print(
            f"  {name}={value}"
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
The strongest exact observation from the current n=6 source layer is
not a closed formula for q1 and q3.

It is a congruence pattern:

    q1 ≡ 1 (mod p),
    q3 ≡ 1 (mod p),

while

    q1 ≡ 1 (mod q),
    q3 ≡ 0 (mod q),

with

    v_q(q3)=1.

Equivalently,

    q3 = q*u_q,

where u_q is a q-adic unit.

For the current instance,

    u_q = q3/q = 8264989.

The next genuine n=pq case should be tested specifically against this
signature.

If the same pattern occurs for n=10, n=15, n=21, etc., that would be
much stronger evidence for a universal terminal-source congruence law.

This experiment therefore prepares the exact checklist for the next
independent n=pq instance.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    core_signature = (
        q1_mod_p == 1 % P
        and q3_mod_p == 1 % P
        and q1_mod_q == 1 % Q
        and q3_mod_q == 0
        and v_q_q3 == 1
        and q3_divisible_by_q
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  q1_equals_1_mod_p="
        f"{q1_mod_p == 1 % P}"
    )

    print(
        f"  q3_equals_1_mod_p="
        f"{q3_mod_p == 1 % P}"
    )

    print(
        f"  q1_equals_1_mod_q="
        f"{q1_mod_q == 1}"
    )

    print(
        f"  q3_equals_0_mod_q="
        f"{q3_mod_q == 0}"
    )

    print(
        f"  q3_exact_q_adic_order_1="
        f"{v_q_q3 == 1}"
    )

    print(
        f"  normalized_q_source_exists="
        f"{q3_divisible_by_q}"
    )

    print(
        f"  observed_terminal_source_signature="
        f"{core_signature}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=only_n=6_has_been_supplied"
    )

    print(
        f"  failures={0 if core_signature else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={core_signature}"
    )

    print()
    print("EXPERIMENT 222 COMPLETE")


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

