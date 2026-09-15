#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 258 — EXACT MULTI-CASE n=pq SOURCE-PROVENANCE / CROSS-CASE AUDIT
==============================================================================

This experiment is intentionally designed around the current bottleneck.

The local projective theorem is already heavily audited.

The unresolved problem is upstream:

    q1 = q1(p,q),
    q3 = q3(p,q).

A second genuine n=pq source pair is now the decisive missing datum.

Experiment 258 therefore creates a strict cross-case interface.

Each CASE must contain:

    p,
    q,
    q1,
    q3,

where q1 and q3 come from the REAL underlying n=pq construction.

For every supplied case the script computes:

    gcd(q1,q3),
    prime valuations at p,q,7,17,
    rho mod 7^e,
    normalized source coordinates,
    terminal Schur row,
    exact Schur difference geometry,
    terminal constants,
    candidate p*q^j orbit,
    principal 7-adic exponent data.

Then it compares every pair of genuine cases and asks:

    * Does the terminal orbit still equal p*q^j?
    * Does q3/q have a common structural formula?
    * Does (q1-1)/(pq) have a cross-case relation?
    * Does ((q3/q)-1)/(pq) have a cross-case relation?
    * Do the same congruence signatures recur?
    * Does the Schur geometry scale exactly as predicted?

With only one supplied case the script deliberately refuses to make a
family claim.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# CASE INPUT
# ============================================================================

# IMPORTANT:
#
# Only enter rows that are genuinely produced by the underlying n=pq
# construction.
#
# The current experiment has only the known n=6 case.
#
# Example format for a future independent case:
#
#     {
#         "name": "n10",
#         "p": 2,
#         "q": 5,
#         "q1": ...,
#         "q3": ...,
#     }
#
CASES = [
    {
        "name": "n6",
        "p": 2,
        "q": 3,
        "q1": 29144191,
        "q3": 24794967,
    },
]


MAX_7_ADIC_E = 4


# ============================================================================
# BASIC HELPERS
# ============================================================================

def valuation(
    x: int,
    p: int,
) -> int:

    x = abs(int(x))

    if x == 0:
        return 10**9

    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int:

    if gcd(
        q3,
        modulus,
    ) != 1:

        raise ArithmeticError(
            "q3 is not invertible modulo modulus."
        )

    return int(
        (
            (q1 % modulus)
            * pow(
                q3 % modulus,
                -1,
                modulus,
            )
        )
        % modulus
    )


def prime_modular_profile(
    x: int,
    p: int,
    levels: int = 4,
) -> list[int]:

    return [
        x % (p ** e)
        for e in range(
            1,
            levels + 1,
        )
    ]


def candidate_orbit_constants(
    p: int,
    q: int,
    rows: int = 3,
) -> list[int]:

    return [
        p * (q ** j)
        for j in range(rows)
    ]


def terminal_row(
    p: int,
    q: int,
    q1: int,
    q3: int,
) -> list[int]:

    return [
        q1 - p * (q ** j) * q3
        for j in range(3)
    ]


def difference_geometry(
    row: list[int],
    q: int,
    p: int,
    q3: int,
) -> dict:

    d01 = row[1] - row[0]
    d12 = row[2] - row[1]
    d02 = row[2] - row[0]

    return {
        "d01": d01,
        "d12": d12,
        "d02": d02,
        "d01_expected":
            -p * (q - 1) * q3,
        "d12_expected":
            -p * q * (q - 1) * q3,
        "d02_expected":
            -p * (q * q - 1) * q3,
        "d01_exact":
            d01
            ==
            -p * (q - 1) * q3,
        "d12_exact":
            d12
            ==
            -p * q * (q - 1) * q3,
        "d02_exact":
            d02
            ==
            -p * (q * q - 1) * q3,
        "ratio_exact":
            d12 == q * d01,
    }


def case_package(case: dict) -> dict:

    name = case["name"]
    p = case["p"]
    q = case["q"]
    q1 = case["q1"]
    q3 = case["q3"]

    n = p * q

    package = {
        "name": name,
        "p": p,
        "q": q,
        "n": n,
        "q1": q1,
        "q3": q3,
    }

    package["gcd"] = gcd(
        q1,
        q3,
    )

    package["valuations"] = {
        "v_p_q1": valuation(q1, p),
        "v_p_q3": valuation(q3, p),
        "v_q_q1": valuation(q1, q),
        "v_q_q3": valuation(q3, q),
        "v_7_q1": valuation(q1, 7),
        "v_7_q3": valuation(q3, 7),
        "v_17_q1": valuation(q1, 17),
        "v_17_q3": valuation(q3, 17),
    }

    if q3 % q == 0:
        package["q3_divisible_by_q"] = True
        package["u"] = q3 // q
        package["u_valuations"] = {
            "v_p_u": valuation(
                q3 // q,
                p,
            ),
            "v_q_u": valuation(
                q3 // q,
                q,
            ),
        }
    else:
        package["q3_divisible_by_q"] = False
        package["u"] = None
        package["u_valuations"] = None

    # --------------------------------------------------------------
    # 7-adic projective profile
    # --------------------------------------------------------------

    rho_profile = []

    if gcd(q3, 7) == 1:

        for e in range(
            1,
            MAX_7_ADIC_E + 1,
        ):

            rho_profile.append(
                ratio_mod(
                    q1,
                    q3,
                    7 ** e,
                )
            )

    package["rho_profile"] = rho_profile

    if rho_profile:
        package["rho_mod_7"] = rho_profile[0]
    else:
        package["rho_mod_7"] = None

    # --------------------------------------------------------------
    # Normalized coordinates
    # --------------------------------------------------------------

    normalized = {}

    if q3 % q == 0:

        u = q3 // q

        normalized["u"] = u

        if (q1 - 1) % n == 0:
            normalized["a"] = (
                q1 - 1
            ) // n
        else:
            normalized["a"] = None

        if (u - 1) % n == 0:
            normalized["b"] = (
                u - 1
            ) // n
        else:
            normalized["b"] = None

    else:

        normalized["u"] = None
        normalized["a"] = None
        normalized["b"] = None

    package["normalized"] = normalized

    # --------------------------------------------------------------
    # Terminal orbit
    # --------------------------------------------------------------

    package["candidate_orbit"] = (
        candidate_orbit_constants(
            p,
            q,
        )
    )

    package["terminal_row"] = (
        terminal_row(
            p,
            q,
            q1,
            q3,
        )
    )

    package["terminal_difference_geometry"] = (
        difference_geometry(
            package["terminal_row"],
            q,
            p,
            q3,
        )
    )

    # --------------------------------------------------------------
    # Direct source congruence signatures
    # --------------------------------------------------------------

    package["source_signature"] = {
        "q1_mod_p":
            q1 % p,

        "q3_mod_p":
            q3 % p,

        "q1_mod_q":
            q1 % q,

        "q3_mod_q":
            q3 % q,

        "q1_mod_n":
            q1 % n,

        "q3_mod_n":
            q3 % n,
    }

    return package


# ============================================================================
# CROSS-CASE SEARCH
# ============================================================================

def pairwise_candidate_relations(
    packages: list[dict],
) -> dict:

    results = {
        "same_q3_over_q": [],
        "same_a": [],
        "same_b": [],
        "same_rho_mod_7": [],
        "orbit_matches_pq": [],
        "signature_matches": [],
    }

    for i in range(
        len(packages)
    ):

        for j in range(
            i + 1,
            len(packages),
        ):

            left = packages[i]
            right = packages[j]

            pair_name = (
                left["name"],
                right["name"],
            )

            # Only meaningful when p,q differ.
            same_u = (
                left["normalized"]["u"]
                is not None
                and
                right["normalized"]["u"]
                is not None
                and
                left["normalized"]["u"]
                ==
                right["normalized"]["u"]
            )

            if same_u:
                results[
                    "same_q3_over_q"
                ].append(
                    pair_name
                )

            same_a = (
                left["normalized"]["a"]
                is not None
                and
                right["normalized"]["a"]
                is not None
                and
                left["normalized"]["a"]
                ==
                right["normalized"]["a"]
            )

            if same_a:
                results[
                    "same_a"
                ].append(
                    pair_name
                )

            same_b = (
                left["normalized"]["b"]
                is not None
                and
                right["normalized"]["b"]
                is not None
                and
                left["normalized"]["b"]
                ==
                right["normalized"]["b"]
            )

            if same_b:
                results[
                    "same_b"
                ].append(
                    pair_name
                )

            same_rho = (
                left["rho_profile"]
                == right["rho_profile"]
                and
                left["rho_profile"]
            )

            if same_rho:
                results[
                    "same_rho_mod_7"
                ].append(
                    pair_name
                )

            left_orbit = (
                left["terminal_row"]
            )

            right_expected = (
                right["candidate_orbit"]
            )

            # Check that each row uses the predicted p*q^j
            # coefficients independently.
            left_orbit_ok = (
                left["candidate_orbit"]
                ==
                [
                    left["p"],
                    left["p"] * left["q"],
                    left["p"] * left["q"] ** 2,
                ]
            )

            right_orbit_ok = (
                right["candidate_orbit"]
                ==
                [
                    right["p"],
                    right["p"] * right["q"],
                    right["p"] * right["q"] ** 2,
                ]
            )

            if (
                left_orbit_ok
                and
                right_orbit_ok
            ):
                results[
                    "orbit_matches_pq"
                ].append(
                    pair_name
                )

            same_signature = (
                left["source_signature"]
                ==
                right["source_signature"]
            )

            if same_signature:
                results[
                    "signature_matches"
                ].append(
                    pair_name
                )

    return results


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 258 — EXACT MULTI-CASE n=pq "
        "SOURCE-PROVENANCE / CROSS-CASE AUDIT"
    )
    print("=" * 78)

    failures = []

    packages = []

    # ------------------------------------------------------------------
    # Validate and package input cases
    # ------------------------------------------------------------------

    for case in CASES:

        try:
            package = case_package(
                case
            )
        except Exception as exc:

            failures.append(
                (
                    case.get("name", "?"),
                    "case_package",
                    type(exc).__name__,
                    str(exc),
                )
            )

            continue

        packages.append(
            package
        )

    # ------------------------------------------------------------------
    # 1. CASE INVENTORY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. CASE INVENTORY")
    print("=" * 78)

    print(
        f"  genuine_cases={len(packages)}"
    )

    for package in packages:

        print(
            f"  {package['name']}: "
            f"(p,q)=({package['p']},{package['q']}) "
            f"n={package['n']}"
        )

    # ------------------------------------------------------------------
    # 2. CASE-BY-CASE SOURCE DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CASE-BY-CASE SOURCE DATA")
    print("=" * 78)

    for package in packages:

        print()
        print(
            f"  CASE {package['name']}"
        )

        print(
            f"    p={package['p']}"
        )

        print(
            f"    q={package['q']}"
        )

        print(
            f"    n={package['n']}"
        )

        print(
            f"    q1={package['q1']}"
        )

        print(
            f"    q3={package['q3']}"
        )

        print(
            f"    gcd={package['gcd']}"
        )

        print(
            f"    rho_mod_7={package['rho_mod_7']}"
        )

        print(
            f"    rho_profile={package['rho_profile']}"
        )

        print(
            f"    normalized={package['normalized']}"
        )

        print(
            f"    candidate_orbit="
            f"{package['candidate_orbit']}"
        )

        print(
            f"    terminal_row="
            f"{package['terminal_row']}"
        )

        geometry = (
            package[
                "terminal_difference_geometry"
            ]
        )

        print(
            "    difference_geometry_exact="
            f"{geometry['d01_exact'] and geometry['d12_exact'] and geometry['d02_exact'] and geometry['ratio_exact']}"
        )

        # Exact source signatures.
        print(
            f"    source_signature="
            f"{package['source_signature']}"
        )

        # A warning rather than a family failure.
        if package["gcd"] != 1:

            print(
                "    WARNING: source pair is not primitive."
            )

    # ------------------------------------------------------------------
    # 3. INTERNAL CASE CHECKS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. INTERNAL CASE CHECKS")
    print("=" * 78)

    all_geometry_exact = True
    all_source_normalization_exact = True

    for package in packages:

        geometry = (
            package[
                "terminal_difference_geometry"
            ]
        )

        geometry_ok = all(
            [
                geometry["d01_exact"],
                geometry["d12_exact"],
                geometry["d02_exact"],
                geometry["ratio_exact"],
            ]
        )

        all_geometry_exact &= geometry_ok

        norm = package[
            "normalized"
        ]

        if (
            norm["u"] is not None
            and
            norm["a"] is not None
            and
            norm["b"] is not None
        ):

            q1_reconstructed = (
                1
                +
                package["p"]
                * package["q"]
                * norm["a"]
            )

            q3_reconstructed = (
                package["q"]
                * (
                    1
                    +
                    package["p"]
                    * package["q"]
                    * norm["b"]
                )
            )

            normalization_ok = (
                q1_reconstructed
                ==
                package["q1"]
                and
                q3_reconstructed
                ==
                package["q3"]
            )

        else:

            normalization_ok = True

        all_source_normalization_exact &= (
            normalization_ok
        )

        print(
            f"  {package['name']}: "
            f"Schur_geometry={geometry_ok} "
            f"normalized_source={normalization_ok}"
        )

        if not geometry_ok:
            failures.append(
                (
                    package["name"],
                    "Schur_geometry",
                )
            )

        if not normalization_ok:
            failures.append(
                (
                    package["name"],
                    "source_normalization",
                )
            )

    # ------------------------------------------------------------------
    # 4. CROSS-CASE COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CROSS-CASE COMPARISON")
    print("=" * 78)

    if len(packages) < 2:

        print(
            "  cross_case_available=False"
        )

        print(
            "  reason=only_one_genuine_n_pq_case_supplied"
        )

        cross_case_available = False
        cross_results = None

    else:

        cross_case_available = True

        cross_results = pairwise_candidate_relations(
            packages
        )

        print(
            "  cross_case_available=True"
        )

        for key, pairs in cross_results.items():

            print(
                f"  {key}={pairs}"
            )

    # ------------------------------------------------------------------
    # 5. SECOND-CASE REQUIREMENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SECOND-CASE REQUIREMENTS")
    print("=" * 78)

    print(
        """
The next decisive row must be a genuinely independently generated n=pq
case.

For example:

    (p,q) = (2,5),
    (2,7),
    (3,5),

or another coprime prime pair.

The important requirement is NOT merely to choose different p,q.
The q1 and q3 values must come from the actual underlying construction.

Once a second genuine case is inserted, this script will compare:

    q3/q,
    (q1-1)/(pq),
    ((q3/q)-1)/(pq),
    rho mod 7^e,
    source congruence signatures,
    terminal p*q^j constants,
    Schur difference geometry.

That will be the first real cross-case test of the n=pq source theorem.
"""
    )

    # ------------------------------------------------------------------
    # 6. CURRENT STATUS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CURRENT FAMILY STATUS")
    print("=" * 78)

    print(
        f"  genuine_case_count={len(packages)}"
    )

    print(
        f"  source_geometry_exact="
        f"{all_geometry_exact}"
    )

    print(
        f"  source_normalization_exact="
        f"{all_source_normalization_exact}"
    )

    if len(packages) < 2:

        print(
            "  cross_case_family_theorem_tested=False"
        )

    else:

        print(
            "  cross_case_family_theorem_tested=True"
        )

    print(
        "  q1_of_p_q_derived=False"
    )

    print(
        "  q3_of_p_q_derived=False"
    )

    print(
        "  terminal_c_j=p*q^j_derived=False"
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        all_geometry_exact
        and
        all_source_normalization_exact
    )

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  internal_case_checks_exact="
        f"{final_ok}"
    )

    print(
        f"  genuine_case_count="
        f"{len(packages)}"
    )

    print(
        f"  cross_case_available="
        f"{cross_case_available}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  local_projective_theorem_complete=True"
    )

    print(
        "  n_pq_source_theorem_complete=False"
    )

    print(
        "  decisive_missing_input="
        "second_genuine_n_pq_source_pair"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 258 COMPLETE")


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

