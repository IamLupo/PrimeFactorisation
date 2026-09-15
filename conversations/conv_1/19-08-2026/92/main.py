#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 257 — EXACT n=pq SOURCE-FORMULA CANDIDATE / INVARIANT SEARCH
==============================================================================

Experiment 256 isolated the remaining upstream problem:

    q1 = q1(p,q),
    q3 = q3(p,q).

With only one genuine n=pq instance available, no universal formula can
be proved. The correct next step is therefore a disciplined candidate
search, not an invented closed form.

Experiment 257 searches a finite, explicitly stated family of low-complexity
candidate expressions built from p and q and tests them against every exact
source-side identity currently known for n=6.

The search includes:

    * monomials p^a q^b;
    * signed monomials;
    * short linear combinations;
    * combinations involving n=pq;
    * normalized quantities

          q3/q,
          (q1-1)/(pq),
          ((q3/q)-1)/(pq).

The script does NOT claim that a surviving candidate is a theorem.
It merely identifies the lowest-complexity formulas compatible with the
known n=6 data and the exact congruence interface.

The output is intended to tell us what symbolic n=pq data should be
collected next.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd


# ============================================================================
# KNOWN n=6 SOURCE
# ============================================================================

P = 2
Q = 3
N = P * Q

Q1 = 29144191
Q3 = 24794967

U = Q3 // Q
A = (Q1 - 1) // (P * Q)
B = (U - 1) // (P * Q)

# Exact downstream identities.
S0 = Q1 - P * Q3
S1 = Q1 - P * Q * Q3
S2 = Q1 - P * Q * Q * Q3


# ============================================================================
# CANDIDATE EXPRESSION ENGINE
# ============================================================================

def monomial(
    p: int,
    q: int,
    ap: int,
    aq: int,
) -> int:

    return (
        p ** ap
        * q ** aq
    )


def candidate_library(
    p: int,
    q: int,
) -> dict[str, Fraction]:

    n = p * q

    values: dict[str, Fraction] = {}

    # --------------------------------------------------------------
    # Basic quantities
    # --------------------------------------------------------------

    values["1"] = Fraction(1)
    values["p"] = Fraction(p)
    values["q"] = Fraction(q)
    values["n=pq"] = Fraction(n)

    # --------------------------------------------------------------
    # Monomials
    # --------------------------------------------------------------

    for ap in range(0, 4):

        for aq in range(0, 5):

            name = (
                f"p^{ap}q^{aq}"
            )

            values[name] = Fraction(
                monomial(
                    p,
                    q,
                    ap,
                    aq,
                )
            )

    # --------------------------------------------------------------
    # Simple combinations
    # --------------------------------------------------------------

    values["p+q"] = Fraction(
        p + q
    )

    values["q-p"] = Fraction(
        q - p
    )

    values["p-q"] = Fraction(
        p - q
    )

    values["pq+1"] = Fraction(
        p * q + 1
    )

    values["pq-1"] = Fraction(
        p * q - 1
    )

    values["p+q+1"] = Fraction(
        p + q + 1
    )

    values["q^2-p"] = Fraction(
        q * q - p
    )

    values["p*q+q"] = Fraction(
        p * q + q
    )

    values["p*q-p"] = Fraction(
        p * q - p
    )

    return values


# ============================================================================
# FACTORIZED RELATION SEARCH
# ============================================================================

def relation_value(
    name1: str,
    value1: Fraction,
    name2: str,
    value2: Fraction,
) -> Fraction:

    return value1 - value2


def exact_integer_relation(
    target: int,
    values: dict[str, Fraction],
    max_coeff: int = 10,
) -> list[tuple[int, str, int, str]]:

    """
    Search

        target = a*X + b*Y

    with small integer coefficients.
    """

    found = []

    items = list(
        values.items()
    )

    for i in range(
        len(items)
    ):

        name_x, x = items[i]

        for j in range(
            i,
            len(items)
        ):

            name_y, y = items[j]

            for a in range(
                -max_coeff,
                max_coeff + 1,
            ):

                for b in range(
                    -max_coeff,
                    max_coeff + 1,
                ):

                    if a == 0 and b == 0:
                        continue

                    candidate = (
                        a * x
                        + b * y
                    )

                    if candidate == target:

                        found.append(
                            (
                                a,
                                name_x,
                                b,
                                name_y,
                            )
                        )

    return found


# ============================================================================
# NORMALIZED SOURCE RELATIONS
# ============================================================================

def normalized_identities() -> dict[str, Fraction]:

    return {
        "q3/q": Fraction(U),
        "(q1-1)/(pq)": Fraction(A),
        "((q3/q)-1)/(pq)": Fraction(B),
        "a-b": Fraction(A - B),
        "a+b": Fraction(A + B),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 257 — EXACT n=pq SOURCE-FORMULA CANDIDATE / "
        "INVARIANT SEARCH"
    )
    print("=" * 78)

    failures = []

    values = candidate_library(
        P,
        Q,
    )

    normalized = normalized_identities()

    # ------------------------------------------------------------------
    # 1. RAW SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. KNOWN n=6 SOURCE TARGETS")
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
        f"  q3/q={U}"
    )
    print(
        f"  a=(q1-1)/(pq)={A}"
    )
    print(
        f"  b=((q3/q)-1)/(pq)={B}"
    )

    # ------------------------------------------------------------------
    # 2. DIRECT MONOMIAL MATCHES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. DIRECT LOW-COMPLEXITY MATCHES")
    print("=" * 78)

    targets = {
        "q1": Q1,
        "q3": Q3,
        "u=q3/q": U,
        "a": A,
        "b": B,
        "a-b": A - B,
        "a+b": A + B,
    }

    direct_hits = {}

    for target_name, target in targets.items():

        hits = []

        for name, value in values.items():

            if value == target:

                hits.append(
                    name
                )

        direct_hits[target_name] = hits

        print(
            f"  {target_name}: "
            f"{hits}"
        )

    # ------------------------------------------------------------------
    # 3. LOW-COMPLEXITY LINEAR RELATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TWO-TERM LINEAR RELATION SEARCH")
    print("=" * 78)

    linear_hits = {}

    for target_name, target in targets.items():

        hits = exact_integer_relation(
            target,
            values,
            max_coeff=5,
        )

        linear_hits[target_name] = hits

        print()
        print(
            f"  target={target_name}"
        )

        if not hits:

            print(
                "    none"
            )

        else:

            for a, x, b, y in hits[:20]:

                print(
                    f"    {a}*{x} "
                    f"+ {b}*{y} = {target}"
                )

            if len(hits) > 20:

                print(
                    f"    ... "
                    f"{len(hits)-20} more"
                )

    # ------------------------------------------------------------------
    # 4. NORMALIZED SOURCE SEARCH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. NORMALIZED SOURCE RELATIONS")
    print("=" * 78)

    for name, value in normalized.items():

        hits = [
            expr
            for expr, candidate in values.items()
            if candidate == value
        ]

        print(
            f"  {name}={value}: "
            f"direct_matches={hits}"
        )

    # ------------------------------------------------------------------
    # 5. SIMPLE RATIO TESTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SIMPLE RATIO TESTS")
    print("=" * 78)

    ratio_targets = {
        "q1/q3": Fraction(Q1, Q3),
        "a/b": Fraction(A, B),
        "q3/q": Fraction(U),
    }

    for target_name, target in ratio_targets.items():

        matches = []

        item_list = list(
            values.items()
        )

        for i, (name_x, x) in enumerate(
            item_list
        ):

            if x == 0:
                continue

            for name_y, y in item_list:

                if y == 0:
                    continue

                if (
                    x / y
                    == target
                ):

                    matches.append(
                        (
                            name_x,
                            name_y,
                        )
                    )

        print(
            f"  {target_name}: "
            f"{matches[:20]}"
        )

    # ------------------------------------------------------------------
    # 6. PRIME PARAMETER SIGNATURES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. PRIME-PARAMETER SIGNATURES")
    print("=" * 78)

    signature_checks = {
        "q1 mod p":
            Q1 % P,

        "q3 mod p":
            Q3 % P,

        "q1 mod q":
            Q1 % Q,

        "q3 mod q":
            Q3 % Q,

        "q1 mod pq":
            Q1 % N,

        "q3 mod pq":
            Q3 % N,

        "a mod p":
            A % P,

        "a mod q":
            A % Q,

        "b mod p":
            B % P,

        "b mod q":
            B % Q,
    }

    for name, value in signature_checks.items():

        print(
            f"  {name}={value}"
        )

    # ------------------------------------------------------------------
    # 7. EXACT KNOWN IDENTITIES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. KNOWN EXACT IDENTITIES")
    print("=" * 78)

    identities = {
        "q3=q*u":
            Q3 == Q * U,

        "q1=1+pq*a":
            Q1 == 1 + P * Q * A,

        "u=1+pq*b":
            U == 1 + P * Q * B,

        "s0=q1-p*q3":
            S0 == Q1 - P * Q3,

        "s1=q1-p*q*q3":
            S1 == Q1 - P * Q * Q3,

        "s2=q1-p*q^2*q3":
            S2 == Q1 - P * Q * Q * Q3,
    }

    for name, result in identities.items():

        print(
            f"  {name}={result}"
        )

        if not result:
            failures.append(
                (
                    "identity",
                    name,
                )
            )

    # ------------------------------------------------------------------
    # 8. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This experiment intentionally makes no family claim.

With only one genuine n=pq source pair, infinitely many formulas can
fit the data. Therefore a surviving expression is only a candidate.

The useful output is instead the complexity profile:

    * direct monomial matches;
    * short linear combinations;
    * simple ratios;
    * normalized-coordinate identities;
    * p/q congruence signatures.

The next genuinely decisive step is to obtain a second actual n=pq
source pair, for example another coprime prime pair (p,q).

Then every candidate produced here can be tested cross-case.

In particular, the most important targets are:

    q1(p,q),
    q3(p,q),
    q3/q,
    (q1-1)/(pq),
    ((q3/q)-1)/(pq),

and the hypothesis

    c_j = p*q^j.

"""
    )

    # ------------------------------------------------------------------
    # 9. PRIORITY OUTPUT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. PRIORITY CANDIDATE SUMMARY")
    print("=" * 78)

    for target_name in [
        "q1",
        "q3",
        "u=q3/q",
        "a",
        "b",
    ]:

        print(
            f"  {target_name}: "
            f"direct={direct_hits[target_name]} "
            f"linear_relation_count="
            f"{len(linear_hits[target_name])}"
        )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_known_source_identities_exact="
        f"{final_ok}"
    )

    print(
        f"  candidate_search_completed=True"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  q1_of_p_q_proved=False"
    )

    print(
        "  q3_of_p_q_proved=False"
    )

    print(
        "  terminal_c_j=p*q^j_proved=False"
    )

    print(
        "  next_required_step="
        "second_genuine_n_pq_source_case"
    )

    print(
        "  ALL BASIC CHECKS PASS="
        f"{final_ok}"
    )

    print()
    print("EXPERIMENT 257 COMPLETE")


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
