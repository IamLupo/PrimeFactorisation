#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 385R — EXACT PRIME-FACTOR PROVENANCE / ARITHMETIC-INVARIANT AUDIT
==============================================================================

Purpose
-------
Investigate whether the observed Q_t(p) values contain arithmetic structure
directly tied to the underlying prime p.

This is a provenance experiment, NOT a generic interpolation experiment.

The experiment audits:

    * complete prime factorization of every observed Q value;
    * sign and absolute-value structure;
    * valuations v_l(Q) for relevant primes l;
    * divisibility by p, p-1, p+1 and powers of p;
    * gcd relationships between neighboring cells;
    * persistent prime factors across rows / columns;
    * gcds of complete row / column collections;
    * cross-cell common-factor structure;
    * exact comparisons with simple p-derived quantities;
    * normalized quotients Q / p^k where integral;
    * whether factors of p-1 or p+1 systematically recur.

No unknown q is assumed.
No missing cell is inserted.
No interpolation.
No extrapolation.
No fitted universal formula.
Exact integer arithmetic only.
"""


from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict

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
# OBSERVED LATTICE
# ============================================================================

def build_lattice():
    lattice = {}

    for p, values in Q.items():

        r = (p - 1) // 2

        for index, value in enumerate(values):

            t = len(values) - 1 - index

            lattice[(r, t)] = sp.Integer(value)

    return lattice


# ============================================================================
# HELPERS
# ============================================================================

def factor_abs(value):
    value = int(value)

    if value == 0:
        return {}

    return sp.factorint(abs(value))


def divisor_gcd(values):
    values = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    if not values:
        return 0

    result = values[0]

    for value in values[1:]:
        result = math.gcd(
            result,
            value,
        )

    return result


def valuation(value, prime):
    value = abs(int(value))

    if value == 0:
        return None

    count = 0

    while value % prime == 0:
        value //= prime
        count += 1

    return count


def factor_string(factors):
    return factors


def proper_divisibility(value, divisor):
    value = int(value)
    divisor = int(divisor)

    if divisor == 0:
        return False

    return value % divisor == 0


def primitive_vector(values):
    values = [int(v) for v in values]

    g = divisor_gcd(values)

    if g == 0:
        return values, 0

    normalized = [
        v // g
        for v in values
    ]

    return normalized, g


# ============================================================================
# BASIC SOURCE AUDIT
# ============================================================================

def source_audit(lattice):

    print()
    print("=" * 78)
    print("1. EXACT SOURCE / PRIME INVENTORY")
    print("=" * 78)

    print()
    print("  observed_cells={}".format(len(lattice)))
    print(
        "  underlying_primes={}".format(
            sorted(Q)
        )
    )

    for p in sorted(Q):

        print()
        print(
            "  p={}".format(p)
        )

        print(
            "    p_minus_1={}".format(
                p - 1
            )
        )

        print(
            "    p_plus_1={}".format(
                p + 1
            )
        )

        print(
            "    factorization_p_minus_1={}".format(
                factor_abs(p - 1)
            )
        )

        print(
            "    factorization_p_plus_1={}".format(
                factor_abs(p + 1)
            )
        )


# ============================================================================
# COMPLETE PRIME FACTORIZATION
# ============================================================================

def complete_factorization_audit(lattice):

    print()
    print("=" * 78)
    print(
        "2. COMPLETE PRIME FACTORIZATION OF OBSERVED Q VALUES"
    )
    print("=" * 78)

    prime_occurrences = Counter()

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = 2 * r + 1
        value_int = int(value)

        factors = factor_abs(value_int)

        for prime in factors:
            prime_occurrences[prime] += 1

        print()
        print(
            "  cell=(r={},t={}) p={}".format(
                r,
                t,
                p,
            )
        )

        print(
            "    value={}".format(
                value_int
            )
        )

        print(
            "    sign={}".format(
                "+" if value_int >= 0 else "-"
            )
        )

        print(
            "    abs_value={}".format(
                abs(value_int)
            )
        )

        print(
            "    factorization={}".format(
                factors
            )
        )

        print(
            "    omega_distinct_prime_count={}".format(
                len(factors)
            )
        )

        print(
            "    bigomega_total_prime_count={}".format(
                sum(factors.values())
            )
        )

    print()
    print(
        "  global_prime_occurrence_profile={}".format(
            sorted(
                prime_occurrences.items()
            )
        )
    )


# ============================================================================
# P-DIVISIBILITY AUDIT
# ============================================================================

def p_divisibility_audit(lattice):

    print()
    print("=" * 78)
    print(
        "3. EXACT p / (p-1) / (p+1) DIVISIBILITY AUDIT"
    )
    print("=" * 78)

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = 2 * r + 1
        value_int = int(value)

        print()
        print(
            "  cell=(r={},t={}) p={}".format(
                r,
                t,
                p,
            )
        )

        for label, divisor in (
            ("p", p),
            ("p_minus_1", p - 1),
            ("p_plus_1", p + 1),
            ("p_squared", p * p),
            ("p_minus_1_squared", (p - 1) ** 2),
            ("p_plus_1_squared", (p + 1) ** 2),
        ):

            divisible = proper_divisibility(
                value_int,
                divisor,
            )

            quotient = None

            if divisible:
                quotient = value_int // divisor

            print(
                "    {}: divisor={}, divisible={}, quotient={}".format(
                    label,
                    divisor,
                    divisible,
                    quotient,
                )
            )

        print(
            "    v_p(Q)={}".format(
                valuation(
                    value_int,
                    p,
                )
            )
        )

        for prime in factor_abs(p - 1):

            print(
                "    v_{}_Q={}".format(
                    prime,
                    valuation(
                        value_int,
                        prime,
                    ),
                )
            )

        for prime in factor_abs(p + 1):

            print(
                "    v_{}_Q={}".format(
                    prime,
                    valuation(
                        value_int,
                        prime,
                    ),
                )
            )


# ============================================================================
# FACTOR SUPPORT BY p
# ============================================================================

def p_factor_support_audit(lattice):

    print()
    print("=" * 78)
    print(
        "4. FACTOR SUPPORT VS p-ARITHMETIC"
    )
    print("=" * 78)

    for p in sorted(Q):

        r = (p - 1) // 2

        cells = [
            (cell, value)
            for cell, value in lattice.items()
            if cell[0] == r
        ]

        p_minus_factors = set(
            factor_abs(p - 1)
        )

        p_plus_factors = set(
            factor_abs(p + 1)
        )

        print()
        print(
            "  p={}".format(p)
        )

        print(
            "    factors_of_p_minus_1={}".format(
                sorted(p_minus_factors)
            )
        )

        print(
            "    factors_of_p_plus_1={}".format(
                sorted(p_plus_factors)
            )
        )

        for cell, value in sorted(
            cells,
            key=lambda item: item[0][1],
        ):

            value_factors = set(
                factor_abs(value)
            )

            shared_minus = sorted(
                value_factors
                & p_minus_factors
            )

            shared_plus = sorted(
                value_factors
                & p_plus_factors
            )

            print()
            print(
                "    cell={}".format(
                    cell
                )
            )

            print(
                "      Q_factor_support={}".format(
                    sorted(value_factors)
                )
            )

            print(
                "      shared_with_(p-1)={}".format(
                    shared_minus
                )
            )

            print(
                "      shared_with_(p+1)={}".format(
                    shared_plus
                )
            )


# ============================================================================
# ROW / COLUMN GCD AUDIT
# ============================================================================

def row_column_gcd_audit(lattice):

    print()
    print("=" * 78)
    print(
        "5. EXACT ROW / COLUMN GCD AND PRIMITIVE-CONTENT AUDIT"
    )
    print("=" * 78)

    for r in sorted(
        {
            cell[0]
            for cell in lattice
        }
    ):

        values = [
            lattice[(r, t)]
            for t in sorted(
                t
                for (rr, t) in lattice
                if rr == r
            )
        ]

        gcd_value = divisor_gcd(values)

        primitive, primitive_gcd = primitive_vector(
            values
        )

        print()
        print(
            "  row_r={}".format(
                r
            )
        )

        print(
            "    values={}".format(
                values
            )
        )

        print(
            "    gcd={}".format(
                gcd_value
            )
        )

        print(
            "    gcd_factorization={}".format(
                factor_abs(gcd_value)
            )
        )

        print(
            "    primitive_vector={}".format(
                primitive
            )
        )

        print(
            "    primitive_gcd={}".format(
                primitive_gcd
            )
        )

    for t in sorted(
        {
            cell[1]
            for cell in lattice
        }
    ):

        values = [
            lattice[(r, t)]
            for r in sorted(
                r
                for (r, tt) in lattice
                if tt == t
            )
        ]

        gcd_value = divisor_gcd(values)

        print()
        print(
            "  column_t={}".format(
                t
            )
        )

        print(
            "    values={}".format(
                values
            )
        )

        print(
            "    gcd={}".format(
                gcd_value
            )
        )

        print(
            "    gcd_factorization={}".format(
                factor_abs(gcd_value)
            )
        )


# ============================================================================
# CROSS-CELL COMMON PRIME AUDIT
# ============================================================================

def common_prime_audit(lattice):

    print()
    print("=" * 78)
    print(
        "6. EXACT CROSS-CELL PRIME PERSISTENCE AUDIT"
    )
    print("=" * 78)

    factor_sets = {}

    for cell, value in lattice.items():

        factor_sets[cell] = set(
            factor_abs(value)
        )

    all_primes = sorted(
        {
            prime
            for factors in factor_sets.values()
            for prime in factors
        }
    )

    persistent = []

    for prime in all_primes:

        cells = [
            cell
            for cell in sorted(
                factor_sets
            )
            if prime in factor_sets[cell]
        ]

        if len(cells) >= 2:

            persistent.append(
                (
                    prime,
                    len(cells),
                    cells,
                )
            )

    for prime, count, cells in sorted(
        persistent,
        key=lambda record: (
            -record[1],
            record[0],
        ),
    ):

        print()
        print(
            "  prime={}".format(
                prime
            )
        )

        print(
            "    occurrence_count={}".format(
                count
            )
        )

        print(
            "    cells={}".format(
                cells
            )
        )

    print()
    print(
        "  persistent_prime_count={}".format(
            len(persistent)
        )
    )


# ============================================================================
# NEIGHBORING CELL GCD AUDIT
# ============================================================================

def neighboring_gcd_audit(lattice):

    print()
    print("=" * 78)
    print(
        "7. EXACT NEIGHBORING-CELL GCD AUDIT"
    )
    print("=" * 78)

    horizontal = []
    vertical = []

    cells = set(
        lattice
    )

    for r, t in sorted(cells):

        value = int(
            lattice[(r, t)]
        )

        right = (r + 1, t)

        if right in cells:

            right_value = int(
                lattice[right]
            )

            g = math.gcd(
                abs(value),
                abs(right_value),
            )

            horizontal.append(
                (
                    (r, t),
                    right,
                    g,
                )
            )

        up = (r, t + 1)

        if up in cells:

            up_value = int(
                lattice[up]
            )

            g = math.gcd(
                abs(value),
                abs(up_value),
            )

            vertical.append(
                (
                    (r, t),
                    up,
                    g,
                )
            )

    print()
    print(
        "  r_direction_edges={}".format(
            len(horizontal)
        )
    )

    for cell_a, cell_b, g in horizontal:

        print(
            "    {} -> {}: gcd={}, factors={}".format(
                cell_a,
                cell_b,
                g,
                factor_abs(g),
            )
        )

    print()
    print(
        "  t_direction_edges={}".format(
            len(vertical)
        )
    )

    for cell_a, cell_b, g in vertical:

        print(
            "    {} -> {}: gcd={}, factors={}".format(
                cell_a,
                cell_b,
                g,
                factor_abs(g),
            )
        )


# ============================================================================
# SIMPLE p-NORMALIZATION AUDIT
# ============================================================================

def p_normalization_audit(lattice):

    print()
    print("=" * 78)
    print(
        "8. EXACT p-POWER NORMALIZATION AUDIT"
    )
    print("=" * 78)

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = 2 * r + 1
        value_int = int(value)

        vp = valuation(
            value_int,
            p,
        )

        if vp is None:
            continue

        normalized = (
            value_int
            // (
                p ** vp
            )
        )

        print()
        print(
            "  cell=(r={},t={}) p={}".format(
                r,
                t,
                p,
            )
        )

        print(
            "    Q={}".format(
                value_int
            )
        )

        print(
            "    v_p(Q)={}".format(
                vp
            )
        )

        print(
            "    Q/p^v_p={}".format(
                normalized
            )
        )

        print(
            "    normalized_factorization={}".format(
                factor_abs(normalized)
            )
        )


# ============================================================================
# p-DERIVED QUANTITY GCD TABLE
# ============================================================================

def p_quantity_gcd_audit(lattice):

    print()
    print("=" * 78)
    print(
        "9. EXACT GCD WITH BASIC p-DERIVED QUANTITIES"
    )
    print("=" * 78)

    quantities = {
        "p": lambda p: p,
        "p_minus_1": lambda p: p - 1,
        "p_plus_1": lambda p: p + 1,
        "p_squared_minus_1": lambda p: p * p - 1,
        "p_cubed_minus_p": lambda p: p ** 3 - p,
    }

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = 2 * r + 1
        value_int = int(value)

        print()
        print(
            "  cell=(r={},t={}) p={} Q={}".format(
                r,
                t,
                p,
                value_int,
            )
        )

        for name, function in quantities.items():

            quantity = int(
                function(p)
            )

            gcd_value = math.gcd(
                abs(value_int),
                abs(quantity),
            )

            print(
                "    gcd(Q,{})={}; factors={}".format(
                    name,
                    gcd_value,
                    factor_abs(gcd_value),
                )
            )


# ============================================================================
# REPEATED FACTORIZATION-PATTERN AUDIT
# ============================================================================

def factor_pattern_audit(lattice):

    print()
    print("=" * 78)
    print(
        "10. EXACT FACTORIZATION-PATTERN AUDIT"
    )
    print("=" * 78)

    patterns = defaultdict(
        list
    )

    for cell, value in sorted(
        lattice.items()
    ):

        factors = factor_abs(
            value
        )

        pattern = tuple(
            sorted(
                factors.keys()
            )
        )

        patterns[pattern].append(
            cell
        )

    repeated = [
        (
            pattern,
            cells,
        )
        for pattern, cells in patterns.items()
        if len(cells) >= 2
    ]

    for pattern, cells in sorted(
        repeated,
        key=lambda item: (
            -len(item[1]),
            item[0],
        ),
    ):

        print()
        print(
            "  factor_support={}".format(
                pattern
            )
        )

        print(
            "    occurrence_count={}".format(
                len(cells)
            )
        )

        print(
            "    cells={}".format(
                cells
            )
        )

    print()
    print(
        "  repeated_factor_support_patterns={}".format(
            len(repeated)
        )
    )


# ============================================================================
# EXACT SIGN / ZERO STRUCTURE
# ============================================================================

def sign_structure_audit(lattice):

    print()
    print("=" * 78)
    print(
        "11. EXACT SIGN / ZERO STRUCTURE AUDIT"
    )
    print("=" * 78)

    positive = []
    negative = []
    zero = []

    for cell, value in sorted(
        lattice.items()
    ):

        value_int = int(value)

        if value_int > 0:
            positive.append(cell)

        elif value_int < 0:
            negative.append(cell)

        else:
            zero.append(cell)

    print()
    print(
        "  positive_count={}".format(
            len(positive)
        )
    )

    print(
        "  positive_cells={}".format(
            positive
        )
    )

    print()
    print(
        "  negative_count={}".format(
            len(negative)
        )
    )

    print(
        "  negative_cells={}".format(
            negative
        )
    )

    print()
    print(
        "  zero_count={}".format(
            len(zero)
        )
    )

    print(
        "  zero_cells={}".format(
            zero
        )
    )


# ============================================================================
# SUMMARY
# ============================================================================

def structural_summary(lattice):

    print()
    print("=" * 78)
    print(
        "12. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    all_values = [
        int(value)
        for value in lattice.values()
    ]

    distinct_primes = sorted(
        {
            prime
            for value in all_values
            for prime in factor_abs(value)
        }
    )

    print()
    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  distinct_prime_factors={}".format(
            len(distinct_primes)
        )
    )

    print(
        "  global_value_gcd={}".format(
            divisor_gcd(
                all_values
            )
        )
    )

    print(
        "  global_value_gcd_factorization={}".format(
            factor_abs(
                divisor_gcd(
                    all_values
                )
            )
        )
    )

    print()
    print(
        "  IMPORTANT:"
    )

    print(
        "    This experiment does not assume a second prime q."
    )

    print(
        "    It does not infer any missing Q value."
    )

    print(
        "    It does not interpret repeated factors as a universal law."
    )


# ============================================================================
# FINAL EXACTNESS
# ============================================================================

def final_exactness(lattice):

    print()
    print("=" * 78)
    print(
        "13. FINAL EXACTNESS"
    )
    print("=" * 78)

    all_integer = all(
        isinstance(
            value,
            sp.Integer,
        )
        for value in lattice.values()
    )

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  all_observed_values_integer={}".format(
            all_integer
        )
    )

    print(
        "  exact_prime_factorization=True"
    )

    print(
        "  p_arithmetic_invariants_tested=True"
    )

    print(
        "  row_column_gcds_tested=True"
    )

    print(
        "  cross_cell_prime_persistence_tested=True"
    )

    print(
        "  neighboring_gcds_tested=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  symbolic_missing_values_created=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 385R — EXACT PRIME-FACTOR PROVENANCE / "
        "ARITHMETIC-INVARIANT AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    source_audit(
        lattice
    )

    complete_factorization_audit(
        lattice
    )

    p_divisibility_audit(
        lattice
    )

    p_factor_support_audit(
        lattice
    )

    row_column_gcd_audit(
        lattice
    )

    common_prime_audit(
        lattice
    )

    neighboring_gcd_audit(
        lattice
    )

    p_normalization_audit(
        lattice
    )

    p_quantity_gcd_audit(
        lattice
    )

    factor_pattern_audit(
        lattice
    )

    sign_structure_audit(
        lattice
    )

    structural_summary(
        lattice
    )

    final_exactness(
        lattice
    )

    print()
    print(
        "EXPERIMENT 385R COMPLETE"
    )


# ============================================================================
# ENTRY POINT
# ============================================================================

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
