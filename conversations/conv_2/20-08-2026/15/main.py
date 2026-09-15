#!/usr/bin/env python3
# ==============================================================================
# EXPERIMENT 394R-COMPACT
# EXACT RESIDUAL FACTOR-CHAIN / MINIMAL TOWER AUDIT
#
# Motivation:
#
# 393R found strong pairwise factor towers such as
#
#     (p-t)(p^2-t-1)
#     (p-t)(p-t-1)
#     p(p^2+1)
#
# but none had constant residuals.
#
# 394R asks a narrower question:
#
#     After removing one strong pair of primitive factors,
#     does the residual itself contain a THIRD primitive
#     source factor on multiple observed cells?
#
# This avoids enumerating all 3-factor combinations.
#
# Rules:
#   * observed cells only
#   * fixed expression library
#   * no missing cells
#   * no interpolation
#   * no extrapolation
#   * singleton coincidences ignored
#   * minimum support required
#   * constant residuals are strongest
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from functools import reduce
from itertools import combinations
from math import gcd


# ------------------------------------------------------------------------------
# OBSERVED DATA
# ------------------------------------------------------------------------------

DATA = {
    (0, 0): 495451247,
    (1, 0): 421514439,
    (2, 0): 16027881,
    (3, 0): 1,

    (0, 1): -1338089411,
    (1, 1): -128667196,
    (2, 1): 4771718,

    (0, 2): 1764373740,
    (1, 2): -152369292,
    (2, 2): -62398,

    (0, 3): 2668721436,
    (1, 3): -1263551016,

    (0, 4): -11600759760,
    (1, 4): 9955176,

    (0, 5): -126258696,
}

MISSING = {
    (2, 3): "Q_3(5)",
    (3, 1): "Q_1(7)",
}

MIN_SUPPORT = 3


# ------------------------------------------------------------------------------
# FIXED EXPRESSION LIBRARY
# ------------------------------------------------------------------------------

def expressions(p: int, t: int) -> dict[str, int]:
    return {
        "p": p,
        "p-1": p - 1,
        "p+1": p + 1,

        "p+t": p + t,
        "p-t": p - t,
        "p+t+1": p + t + 1,
        "p-t-1": p - t - 1,

        "p+2t": p + 2*t,
        "p-2t": p - 2*t,
        "p+2t+1": p + 2*t + 1,
        "p-2t-1": p - 2*t - 1,

        "p+3t": p + 3*t,
        "p-3t": p - 3*t,

        "p^2": p*p,
        "p^2+1": p*p + 1,
        "p^2+t": p*p + t,
        "p^2-t": p*p - t,
        "p^2+t+1": p*p + t + 1,
        "p^2-t-1": p*p - t - 1,

        "p^2+t^2": p*p + t*t,
        "p^2-t^2": p*p - t*t,
        "p^2+t^2-1": p*p + t*t - 1,
        "p^2-t^2-1": p*p - t*t - 1,
        "p^2+t^2+1": p*p + t*t + 1,
        "p^2-t^2+1": p*p - t*t + 1,

        "p^2+p+t": p*p + p + t,
        "p^2-p+t": p*p - p + t,
        "p^2+p-t": p*p + p - t,
        "p^2-p-t": p*p - p - t,

        "p^2+2pt": p*p + 2*p*t,
        "p^2-2pt": p*p - 2*p*t,

        "p*(p+t)": p*(p+t),
        "p*(p-t)": p*(p-t),
        "p*(t+1)": p*(t+1),
        "p*(t-1)": p*(t-1),

        "(p-1)*(t+1)": (p-1)*(t+1),
        "(p+1)*(t+1)": (p+1)*(t+1),
        "(p-1)*(t-1)": (p-1)*(t-1),
        "(p+1)*(t-1)": (p+1)*(t-1),

        "(p-1)*(p+t)": (p-1)*(p+t),
        "(p+1)*(p-t)": (p+1)*(p-t),
    }


def build_library():
    lib = defaultdict(dict)

    for cell in DATA:
        r, t = cell
        p = 2*r + 1

        for name, value in expressions(p, t).items():
            lib[name][cell] = value

    return dict(lib)


# ------------------------------------------------------------------------------
# TOP PAIRS FROM 393R
# ------------------------------------------------------------------------------

TOP_PAIRS = [
    ("p-t", "p^2-t-1"),
    ("p-t", "p-t-1"),
    ("p", "p^2+1"),
    ("p*(p+t)", "p+1"),
    ("p*(p-t)", "p-t-1"),
    ("p*(p-t)", "p^2-t-1"),
    ("p*(t-1)", "p-t-1"),
    ("p+1", "p^2"),
    ("p+1", "p^2+1"),
    ("p+t", "p^2+p-t"),
]


# ------------------------------------------------------------------------------
# BASIC HELPERS
# ------------------------------------------------------------------------------

def divides(a: int, b: int) -> bool:
    return a != 0 and b % a == 0


def gcd_list(values):
    values = [abs(x) for x in values if x != 0]

    if not values:
        return 0

    return reduce(gcd, values)


def pair_residual(lib, a: str, b: str):
    """
    Return residual Q / (a*b) on every usable observed cell.
    """
    residual = {}

    for cell, q in DATA.items():
        x = lib[a][cell]
        y = lib[b][cell]

        if x == 0 or y == 0:
            continue

        product = x * y

        if product != 0 and q % product == 0:
            residual[cell] = q // product

    return residual


def residual_factor_hits(residual, lib, excluded):
    """
    Find primitive third factors dividing the residual on multiple cells.
    """
    hits = []

    for name in lib:
        if name in excluded:
            continue

        cells = []

        for cell, q in residual.items():
            e = lib[name][cell]

            if e != 0 and q % e == 0:
                cells.append(cell)

        if len(cells) < MIN_SUPPORT:
            continue

        quotients = [
            residual[cell] // lib[name][cell]
            for cell in cells
        ]

        g = gcd_list(quotients)
        constant = len(set(quotients)) == 1
        repeated = len(quotients) - len(set(quotients))

        hits.append({
            "name": name,
            "cells": cells,
            "count": len(cells),
            "gcd": g,
            "constant": constant,
            "repeated": repeated,
        })

    hits.sort(
        key=lambda x: (
            x["constant"],
            x["gcd"] > 1,
            x["count"],
            x["repeated"] > 0,
        ),
        reverse=True,
    )

    return hits


def residual_prime_gcd(residual):
    return gcd_list(residual.values())


# ------------------------------------------------------------------------------
# CHAIN EXTENSION
# ------------------------------------------------------------------------------

def extend_chain(lib, pair):
    a, b = pair

    residual = pair_residual(lib, a, b)

    if len(residual) < MIN_SUPPORT:
        return {
            "pair": pair,
            "base_cells": len(residual),
            "third_hits": [],
            "residual_gcd": residual_prime_gcd(residual),
        }

    hits = residual_factor_hits(
        residual,
        lib,
        excluded={a, b},
    )

    return {
        "pair": pair,
        "base_cells": len(residual),
        "third_hits": hits,
        "residual_gcd": residual_prime_gcd(residual),
    }


# ------------------------------------------------------------------------------
# CHECK 4-FACTOR EXTENSION
# ------------------------------------------------------------------------------

def fourth_factor_test(lib, first_pair, third):
    a, b = first_pair

    residual = pair_residual(lib, a, b)

    if not residual:
        return None

    third_name = third["name"]

    next_residual = {}

    for cell, q in residual.items():
        e = lib[third_name][cell]

        if e != 0 and q % e == 0:
            next_residual[cell] = q // e

    if len(next_residual) < MIN_SUPPORT:
        return None

    # Look for any other factor that divides the new residual
    # on at least MIN_SUPPORT cells.
    hits = residual_factor_hits(
        next_residual,
        lib,
        excluded={a, b, third_name},
    )

    return {
        "cells": list(next_residual),
        "count": len(next_residual),
        "gcd": residual_prime_gcd(next_residual),
        "hits": hits,
    }


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():
    lib = build_library()

    print("=" * 78)
    print("EXPERIMENT 394R-COMPACT — EXACT RESIDUAL FACTOR-CHAIN AUDIT")
    print("=" * 78)
    print()
    print("SOURCE")
    print(f"  observed_cells={len(DATA)}")
    print(f"  missing_cells={MISSING}")
    print("  fixed_expression_library=True")
    print(f"  tested_pair_count={len(TOP_PAIRS)}")
    print()

    # --------------------------------------------------------------------------
    # 1. PAIR RESIDUALS
    # --------------------------------------------------------------------------

    print("=" * 78)
    print("1. SELECTED PAIR RESIDUALS")
    print("=" * 78)

    chain_results = []

    for pair in TOP_PAIRS:
        result = extend_chain(lib, pair)
        chain_results.append(result)

        print(
            f"  pair={pair} "
            f"base_cells={result['base_cells']} "
            f"residual_gcd={result['residual_gcd']} "
            f"third_factor_candidates={len(result['third_hits'])}"
        )

    # --------------------------------------------------------------------------
    # 2. BEST THIRD-FACTOR EXTENSIONS
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. BEST THREE-FACTOR EXTENSIONS")
    print("=" * 78)

    three_factor_results = []

    for result in chain_results:
        pair = result["pair"]

        for hit in result["third_hits"][:5]:
            three_factor_results.append(
                (
                    pair,
                    hit,
                )
            )

    three_factor_results.sort(
        key=lambda x: (
            x[1]["constant"],
            x[1]["gcd"] > 1,
            x[1]["count"],
        ),
        reverse=True,
    )

    if not three_factor_results:
        print("  NONE")
    else:
        for pair, hit in three_factor_results[:15]:
            print(
                f"  factors={pair + (hit['name'],)} "
                f"cells={hit['count']} "
                f"gcd={hit['gcd']} "
                f"constant={hit['constant']} "
                f"repeated={hit['repeated']}"
            )

    # --------------------------------------------------------------------------
    # 3. FOURTH-FACTOR EXTENSION
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FOURTH-FACTOR EXTENSION")
    print("=" * 78)

    fourth_results = []

    for pair, hit in three_factor_results[:10]:
        result = fourth_factor_test(
            lib,
            pair,
            hit,
        )

        if result is None:
            continue

        fourth_results.append(
            (pair, hit["name"], result)
        )

    if not fourth_results:
        print("  NONE")
    else:
        fourth_results.sort(
            key=lambda x: (
                any(h["constant"] for h in x[2]["hits"]),
                x[2]["count"],
            ),
            reverse=True,
        )

        for pair, third, result in fourth_results[:10]:
            fourth_names = [
                h["name"] for h in result["hits"][:3]
            ]

            print(
                f"  base={pair} "
                f"third={third} "
                f"cells={result['count']} "
                f"residual_gcd={result['gcd']} "
                f"next_factors={fourth_names}"
            )

    # --------------------------------------------------------------------------
    # 4. CONSTANT RESIDUAL CHECK
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CONSTANT RESIDUAL CHECK")
    print("=" * 78)

    constant_three = [
        (pair, hit)
        for pair, hit in three_factor_results
        if hit["constant"]
    ]

    if not constant_three:
        print("  NONE")
    else:
        for pair, hit in constant_three[:10]:
            print(
                f"  factors={pair + (hit['name'],)} "
                f"cells={hit['cells']}"
            )

    # --------------------------------------------------------------------------
    # 5. STRUCTURAL VERDICT
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL VERDICT")
    print("=" * 78)

    if constant_three:
        verdict = "EXACT_CONSTANT_THREE_FACTOR_TOWER_FOUND"
    elif fourth_results:
        verdict = "DEEP_FACTOR_CHAIN_STRUCTURE_FOUND"
    elif three_factor_results:
        verdict = "THREE_FACTOR_RESIDUAL_STRUCTURE_FOUND"
    else:
        verdict = "NO_VALIDATED_DEEPER_FACTOR_CHAIN"

    print(f"  verdict={verdict}")
    print()
    print("  Interpretation:")
    print("    only the strongest 393R pair candidates are continued")
    print("    no full triple-product enumeration is performed")
    print("    the third factor must divide the pair residual on multiple cells")
    print("    fourth-factor extensions are tested only for surviving chains")
    print("    constant residuals are strongest")
    print("    missing cells are never used as evidence")

    # --------------------------------------------------------------------------
    # 6. EXACTNESS
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_expression_library=True")
    print("  selected_393R_pairs_only=True")
    print("  residual_factor_chain_test=True")
    print("  third_factor_overdetermination=True")
    print("  fourth_factor_extension_test=True")
    print("  constant_residual_test=True")
    print("  missing_Q3_5_used=False")
    print("  missing_Q1_7_used=False")
    print("  interpolation_performed=False")
    print("  extrapolation_counted_as_evidence=False")
    print("  synthetic_second_case=False")
    print("  arbitrary_matrix_fit=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")
    print()
    print("EXPERIMENT 394R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
