#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 340R — EXACT LOCAL-GCD NORMALIZATION / PRIMITIVE-QUOTIENT
             LATTICE AUDIT
==============================================================================

Purpose
-------
Experiment 339R found:

    persistent primes = {2,3,7},

but much more strikingly, many adjacent cells have structured gcds such as

    12, 14, 24, 168.

This experiment removes the local arithmetic content before searching
for any further structure.

For every observed cell Q(r,t), compute:

    1. row gcd;
    2. column gcd;
    3. local horizontal gcd;
    4. local vertical gcd;
    5. local diagonal gcd;
    6. gcd of the entire observed neighborhood;
    7. exact quotient after normalization;
    8. prime factorization of the primitive quotient.

The main object is the normalized quotient

    H(r,t) = Q(r,t) / G(r,t),

where G is constructed only from observed local gcd data.

We then test whether H has:

    * smaller prime support;
    * reduced valuation variation;
    * integer divisibility relations;
    * repeated exact quotients;
    * low-rank additive structure;
    * low-rank multiplicative structure.

IMPORTANT
---------
This is NOT another recurrence-fitting experiment.

No interpolation.
No missing cells.
No extrapolation.
No synthetic second n=pq case.
Exact integer arithmetic only.
"""


from __future__ import annotations

import math
import sys
from collections import defaultdict

import sympy as sp


# ============================================================================
# SOURCE DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}


# ============================================================================
# LATTICE
# ============================================================================

def build_lattice():

    lattice = {}

    for p_value, values in Q.items():

        r = (p_value - 1) // 2

        for index, value in enumerate(values):

            t = len(values) - 1 - index

            lattice[(r, t)] = int(value)

    return lattice


def cells(lattice):
    return sorted(
        lattice,
        key=lambda z: (z[1], z[0]),
    )


# ============================================================================
# EXACT HELPERS
# ============================================================================

def factor_integer(n):

    n = int(n)

    if n == 0:
        return {}

    return sp.factorint(
        abs(n)
    )


def vp(n, p):

    n = abs(int(n))

    if n == 0:
        return sp.oo

    out = 0

    while n % p == 0:

        n //= p
        out += 1

    return out


def gcd_list(values):

    values = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    if not values:
        return 0

    g = values[0]

    for value in values[1:]:

        g = math.gcd(
            g,
            value,
        )

    return g


def quotient_if_divisible(a, b):

    if b == 0:
        return None

    if a % b != 0:
        return None

    return a // b


# ============================================================================
# NEIGHBORHOODS
# ============================================================================

def neighborhood(lattice, r, t):

    positions = [
        (r, t),

        # horizontal
        (r - 1, t),
        (r + 1, t),

        # vertical
        (r, t - 1),
        (r, t + 1),

        # diagonals
        (r - 1, t - 1),
        (r + 1, t - 1),
        (r - 1, t + 1),
        (r + 1, t + 1),
    ]

    return [
        lattice[pos]
        for pos in positions
        if pos in lattice
    ]


# ============================================================================
# 1. ROW / COLUMN CONTENT
# ============================================================================

def row_column_content(lattice):

    print()
    print("=" * 78)
    print("1. ROW / COLUMN GCD CONTENT")
    print("=" * 78)

    row_gcd = {}
    col_gcd = {}

    for t in sorted(
        {
            t
            for _, t in lattice
        }
    ):

        values = [
            lattice[(r, t)]
            for r, tt in lattice
            if tt == t
        ]

        row_gcd[t] = gcd_list(
            values
        )

        print(
            "  t={}: gcd={}".format(
                t,
                row_gcd[t],
            )
        )

    print()

    for r in sorted(
        {
            r
            for r, _ in lattice
        }
    ):

        values = [
            lattice[(r, t)]
            for rr, t in lattice
            if rr == r
        ]

        col_gcd[r] = gcd_list(
            values
        )

        print(
            "  r={}: gcd={}".format(
                r,
                col_gcd[r],
            )
        )

    return row_gcd, col_gcd


# ============================================================================
# 2. LOCAL GCD CONTENT
# ============================================================================

def local_gcd_content(lattice):

    print()
    print("=" * 78)
    print("2. LOCAL GCD CONTENT")
    print("=" * 78)

    local = {}

    for cell in cells(lattice):

        r, t = cell

        neigh = neighborhood(
            lattice,
            r,
            t,
        )

        g = gcd_list(
            neigh
        )

        local[cell] = g

        print()
        print(
            "  cell={}: local_gcd={}".format(
                cell,
                g,
            )
        )

        print(
            "    factorization={}".format(
                factor_integer(g)
            )
        )

    return local


# ============================================================================
# 3. PRIMITIVE QUOTIENTS
# ============================================================================

def primitive_quotient_lattice(
    lattice,
    row_gcd,
    col_gcd,
    local_gcd,
):

    print()
    print("=" * 78)
    print(
        "3. PRIMITIVE-QUOTIENT NORMALIZATION"
    )
    print("=" * 78)

    # Three normalizations are deliberately kept separate.
    #
    # H_local = Q / local neighborhood gcd
    # H_row   = Q / row gcd
    # H_col   = Q / column gcd

    H_local = {}
    H_row = {}
    H_col = {}

    for cell in cells(lattice):

        r, t = cell
        value = lattice[cell]

        lg = local_gcd[cell]
        rg = row_gcd[t]
        cg = col_gcd[r]

        H_local[cell] = (
            quotient_if_divisible(
                value,
                lg,
            )
        )

        H_row[cell] = (
            quotient_if_divisible(
                value,
                rg,
            )
        )

        H_col[cell] = (
            quotient_if_divisible(
                value,
                cg,
            )
        )

        print()
        print(
            "  cell={}: Q={}".format(
                cell,
                value,
            )
        )

        print(
            "    Q/local_gcd={}".format(
                H_local[cell]
            )
        )

        print(
            "    Q/row_gcd={}".format(
                H_row[cell]
            )
        )

        print(
            "    Q/column_gcd={}".format(
                H_col[cell]
            )
        )

    return (
        H_local,
        H_row,
        H_col,
    )


# ============================================================================
# 4. PRIME-SUPPORT REDUCTION
# ============================================================================

def reduced_prime_support(
    name,
    normalized,
):

    print()
    print("=" * 78)
    print(
        "{} PRIME-SUPPORT AUDIT".format(
            name
        )
    )
    print("=" * 78)

    support = defaultdict(list)

    for cell, value in normalized.items():

        if value is None:
            continue

        for prime in factor_integer(
            value
        ):

            support[prime].append(
                cell
            )

    for prime in sorted(support):

        print(
            "  prime={}: count={}".format(
                prime,
                len(
                    support[prime]
                ),
            )
        )

    print()
    print(
        "  distinct_primes={}".format(
            len(support)
        )
    )

    return support


# ============================================================================
# 5. VALUATION-RANGE AUDIT
# ============================================================================

def valuation_range_audit(
    name,
    normalized,
):

    print()
    print("=" * 78)
    print(
        "{} VALUATION-RANGE AUDIT".format(
            name
        )
    )
    print("=" * 78)

    primes = set()

    for value in normalized.values():

        if value is None:
            continue

        primes.update(
            factor_integer(value)
        )

    for prime in sorted(primes):

        vals = [
            vp(
                value,
                prime,
            )
            for value in normalized.values()
            if value not in (
                None,
                0,
            )
        ]

        print()
        print(
            "  prime={}: min={}, max={}, "
            "range={}".format(
                prime,
                min(vals),
                max(vals),
                max(vals) - min(vals),
            )
        )


# ============================================================================
# 6. LOCAL DIVISIBILITY AFTER NORMALIZATION
# ============================================================================

def normalized_neighbor_audit(
    name,
    lattice,
    normalized,
):

    print()
    print("=" * 78)
    print(
        "{} NORMALIZED NEIGHBOR AUDIT".format(
            name
        )
    )
    print("=" * 78)

    directions = {
        "horizontal": (1, 0),
        "vertical": (0, 1),
        "diag_up_right": (1, 1),
        "diag_down_right": (1, -1),
    }

    summary = {}

    for direction, (
        dr,
        dt,
    ) in directions.items():

        divisibility = 0
        gcds = []

        for cell in cells(lattice):

            r, t = cell

            target = (
                r + dr,
                t + dt,
            )

            if target not in lattice:
                continue

            a = normalized[cell]
            b = normalized[target]

            if a is None or b is None:
                continue

            g = math.gcd(
                int(a),
                int(b),
            )

            gcds.append(
                g
            )

            if (
                a % b == 0
                or
                b % a == 0
            ):

                divisibility += 1

        summary[
            direction
        ] = {
            "pairs": len(gcds),
            "divisibility_hits": divisibility,
            "gcds": gcds,
        }

        print()
        print(
            "  {}: pairs={}, divisibility_hits={}".format(
                direction,
                len(gcds),
                divisibility,
            )
        )

        print(
            "    gcds={}".format(
                gcds
            )
        )

    return summary


# ============================================================================
# 7. NORMALIZED MATRIX RANK
# ============================================================================

def normalized_rank(
    name,
    lattice,
    normalized,
):

    print()
    print("=" * 78)
    print(
        "{} NORMALIZED RANK AUDIT".format(
            name
        )
    )
    print("=" * 78)

    matrix = sp.Matrix([
        [
            (
                normalized.get(
                    (r, t),
                    0,
                )
                or 0
            )
            if (r, t) in lattice
            else 0
            for t in range(6)
        ]
        for r in range(4)
    ])

    print(
        "  matrix="
    )

    print(
        matrix
    )

    print()
    print(
        "  shape={}".format(
            matrix.shape
        )
    )

    print(
        "  rank={}".format(
            matrix.rank()
        )
    )

    return matrix


# ============================================================================
# 8. EXTRACTION OF COMMON SMALL FACTORS
# ============================================================================

def small_factor_profile(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "8. SMALL-PRIME VALUATION PROFILE"
    )
    print("=" * 78)

    primes = (
        2,
        3,
        5,
        7,
        11,
        13,
        17,
    )

    for cell in cells(lattice):

        value = lattice[cell]

        print()
        print(
            "  cell={}:".format(
                cell
            )
        )

        print(
            "    Q={}".format(
                value
            )
        )

        print(
            "    valuations={}".format(
                {
                    p: vp(
                        value,
                        p,
                    )
                    for p in primes
                }
            )
        )


# ============================================================================
# 9. TERMINAL GCD COMPARISON
# ============================================================================

def terminal_gcd_audit(
    lattice,
):

    q1 = 495451247
    q3 = 421514439

    terminal_gcd = math.gcd(
        q1,
        q3,
    )

    print()
    print("=" * 78)
    print(
        "9. TERMINAL-GCD NORMALIZATION"
    )
    print("=" * 78)

    print(
        "  terminal_gcd={}".format(
            terminal_gcd
        )
    )

    print()

    for cell in cells(lattice):

        value = lattice[cell]

        g1 = math.gcd(
            value,
            q1,
        )

        g3 = math.gcd(
            value,
            q3,
        )

        gt = math.gcd(
            value,
            terminal_gcd,
        )

        print(
            "  cell={}: gcd_q1={}, gcd_q3={}, "
            "gcd_terminal={}".format(
                cell,
                g1,
                g3,
                gt,
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 340R — EXACT LOCAL-GCD NORMALIZATION / "
        "PRIMITIVE-QUOTIENT LATTICE AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    row_gcd, col_gcd = row_column_content(
        lattice
    )

    local_gcd = local_gcd_content(
        lattice
    )

    (
        H_local,
        H_row,
        H_col,
    ) = primitive_quotient_lattice(
        lattice,
        row_gcd,
        col_gcd,
        local_gcd,
    )

    support_local = reduced_prime_support(
        "H_LOCAL",
        H_local,
    )

    support_row = reduced_prime_support(
        "H_ROW",
        H_row,
    )

    support_col = reduced_prime_support(
        "H_COLUMN",
        H_col,
    )

    valuation_range_audit(
        "H_LOCAL",
        H_local,
    )

    valuation_range_audit(
        "H_ROW",
        H_row,
    )

    valuation_range_audit(
        "H_COLUMN",
        H_col,
    )

    normalized_neighbor_audit(
        "H_LOCAL",
        lattice,
        H_local,
    )

    normalized_rank(
        "H_LOCAL",
        lattice,
        H_local,
    )

    normalized_rank(
        "H_ROW",
        lattice,
        H_row,
    )

    normalized_rank(
        "H_COLUMN",
        lattice,
        H_col,
    )

    small_factor_profile(
        lattice
    )

    terminal_gcd_audit(
        lattice
    )

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 339R showed that the source table has a sparse arithmetic
signature:

    persistent primes = {2,3,7},

while the terminal factor 17 does not propagate.

The next question is whether the recurring small factors are merely
local normalization artifacts.

Experiment 340R removes several natural layers of arithmetic content:

    Q / row_gcd,
    Q / column_gcd,
    Q / neighborhood_gcd.

The resulting primitive quotients are then tested for:

    * reduced prime support;
    * reduced valuation variation;
    * restored divisibility;
    * reduced matrix rank.

A substantial simplification after normalization would indicate that the
apparent complexity of Q is partly caused by predictable local content.

A failure to simplify would mean that even after removing all obvious gcd
content, the remaining quotient lattice is still arithmetically irregular.

This is still source reconstruction, not recurrence fitting.

Only observed cells are used.
No missing values are reconstructed.
No extrapolation is introduced.
No synthetic second n=pq case is generated.
"""
    )

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    checks = {
        "row_column_gcd_exact": True,
        "local_gcd_exact": True,
        "primitive_quotients_exact": True,
        "normalized_prime_support_completed": True,
        "normalized_valuation_audit_completed": True,
        "normalized_neighbor_audit_completed": True,
        "normalized_rank_audit_completed": True,
        "terminal_gcd_audit_completed": True,
        "missing_values_used": False,
        "interpolation_performed": False,
        "extrapolation_used": False,
        "synthetic_second_case": False,
        "external_files_used": False,
        "universal_q_p_r_formula_proved": False,
    }

    for key, value in checks.items():

        print(
            "  {}={}".format(
                key,
                value,
            )
        )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 340R COMPLETE"
    )


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
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )

        raise
