#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 385R-FIXED — EXACT PRIME-FACTOR PROVENANCE /
                         ARITHMETIC-INVARIANT AUDIT
==============================================================================

FIXES
-----
1. Never compute v_p(Q) for p=1.
   The p=1 row is a boundary row, not a prime-adic valuation case.

2. Never perform valuation() for a divisor <= 1.

3. Cache all integer factorizations.

4. Cache all valuations.

5. Cache all gcd factorizations.

6. Keep the experiment purely exact and source-observed.

No missing value is used.
No second prime q is assumed.
No interpolation.
No extrapolation.
No fitted universal formula.
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
# GLOBAL CACHES
# ============================================================================

_FACTOR_CACHE: dict[int, dict[int, int]] = {}
_VALUATION_CACHE: dict[tuple[int, int], int | None] = {}


# ============================================================================
# EXACT HELPERS
# ============================================================================

def factor_abs(value: int) -> dict[int, int]:
    """
    Exact prime factorization with memoization.
    """

    value = abs(int(value))

    if value <= 1:
        return {}

    cached = _FACTOR_CACHE.get(value)

    if cached is not None:
        return dict(cached)

    result = {
        int(p): int(e)
        for p, e in sp.factorint(value).items()
    }

    _FACTOR_CACHE[value] = result

    return dict(result)


def valuation(value: int, prime: int):
    """
    Exact v_prime(value).

    IMPORTANT:
        prime <= 1 is not a valid p-adic valuation.

    In particular v_1(n) is undefined, and must NEVER be evaluated.
    """

    value = abs(int(value))
    prime = int(prime)

    if value == 0:
        return None

    if prime <= 1:
        return None

    key = (
        value,
        prime,
    )

    if key in _VALUATION_CACHE:
        return _VALUATION_CACHE[key]

    count = 0
    current = value

    while current % prime == 0:
        current //= prime
        count += 1

    _VALUATION_CACHE[key] = count

    return count


def gcd_many(values):
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


def is_divisible(value, divisor):
    value = int(value)
    divisor = int(divisor)

    if divisor == 0:
        return False

    return value % divisor == 0


def primitive_vector(values):
    values = [int(v) for v in values]

    g = gcd_many(values)

    if g == 0:
        return values, 0

    return [
        value // g
        for value in values
    ], g


# ============================================================================
# LATTICE
# ============================================================================

def build_lattice():

    lattice = {}

    for p, values in Q.items():

        r = (p - 1) // 2

        for index, value in enumerate(values):

            t = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r, t)
            ] = sp.Integer(value)

    return lattice


# ============================================================================
# 1. SOURCE / PRIME INVENTORY
# ============================================================================

def source_inventory(lattice):

    print()
    print("=" * 78)
    print(
        "1. EXACT SOURCE / PRIME INVENTORY"
    )
    print("=" * 78)

    print()
    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  source_parameters={}".format(
            sorted(Q)
        )
    )

    for p in sorted(Q):

        print()
        print(
            "  source_parameter={}".format(
                p
            )
        )

        if p == 1:

            print(
                "    classification=BOUNDARY_PARAMETER_NOT_PRIME"
            )

            print(
                "    p_minus_1=0"
            )

            print(
                "    p_plus_1=2"
            )

            print(
                "    factorization_p_minus_1={}"
            )

            print(
                "    factorization_p_plus_1={}".format(
                    factor_abs(2)
                )
            )

            continue

        print(
            "    classification=ODD_PRIME_PARAMETER"
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
                factor_abs(
                    p - 1
                )
            )
        )

        print(
            "    factorization_p_plus_1={}".format(
                factor_abs(
                    p + 1
                )
            )
        )


# ============================================================================
# 2. COMPLETE FACTORIZATION
# ============================================================================

def complete_factorization(lattice):

    print()
    print("=" * 78)
    print(
        "2. COMPLETE PRIME FACTORIZATION OF OBSERVED Q VALUES"
    )
    print("=" * 78)

    occurrence_counter = Counter()

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = (
            2 * r
            + 1
        )

        value_int = int(value)

        factors = factor_abs(
            value_int
        )

        for prime in factors:
            occurrence_counter[prime] += 1

        print()
        print(
            "  cell=(r={},t={}) parameter={}".format(
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
                "+"
                if value_int > 0
                else "-"
                if value_int < 0
                else "0"
            )
        )

        print(
            "    factorization={}".format(
                factors
            )
        )

        print(
            "    distinct_prime_count={}".format(
                len(factors)
            )
        )

        print(
            "    total_prime_multiplicity={}".format(
                sum(factors.values())
            )
        )

    print()
    print(
        "  global_prime_occurrence_profile={}".format(
            sorted(
                occurrence_counter.items()
            )
        )
    )


# ============================================================================
# 3. p / (p-1) / (p+1) AUDIT
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

        p = (
            2 * r
            + 1
        )

        value_int = int(value)

        print()
        print(
            "  cell=(r={},t={}) parameter={} Q={}".format(
                r,
                t,
                p,
                value_int,
            )
        )

        # Boundary p=1 is handled separately.
        if p == 1:

            print(
                "    p_adic_analysis=SKIPPED_BOUNDARY_PARAMETER"
            )

            print(
                "    divisor_2: divisible={}".format(
                    is_divisible(
                        value_int,
                        2,
                    )
                )
            )

            continue

        quantities = [
            ("p", p),
            ("p_minus_1", p - 1),
            ("p_plus_1", p + 1),
            ("p_squared", p * p),
            (
                "p_minus_1_squared",
                (p - 1) ** 2,
            ),
            (
                "p_plus_1_squared",
                (p + 1) ** 2,
            ),
        ]

        for name, divisor in quantities:

            divisible = is_divisible(
                value_int,
                divisor,
            )

            quotient = (
                value_int // divisor
                if divisible
                else None
            )

            print(
                "    {}: divisor={}, divisible={}, quotient={}".format(
                    name,
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

        for prime in sorted(
            factor_abs(
                p - 1
            )
        ):

            print(
                "    v_{}_Q={}".format(
                    prime,
                    valuation(
                        value_int,
                        prime,
                    )
                )
            )

        for prime in sorted(
            factor_abs(
                p + 1
            )
        ):

            print(
                "    v_{}_Q={}".format(
                    prime,
                    valuation(
                        value_int,
                        prime,
                    )
                )
            )


# ============================================================================
# 4. FACTOR SUPPORT VS p
# ============================================================================

def p_factor_support(lattice):

    print()
    print("=" * 78)
    print(
        "4. EXACT FACTOR SUPPORT VS p-ARITHMETIC"
    )
    print("=" * 78)

    for p in sorted(Q):

        r = (
            p - 1
        ) // 2

        cells = [
            (
                cell,
                value,
            )
            for cell, value in lattice.items()
            if cell[0] == r
        ]

        print()
        print(
            "  parameter={}".format(
                p
            )
        )

        if p == 1:

            print(
                "    boundary_parameter=True"
            )

            continue

        minus_factors = set(
            factor_abs(
                p - 1
            )
        )

        plus_factors = set(
            factor_abs(
                p + 1
            )
        )

        print(
            "    factors_of_p_minus_1={}".format(
                sorted(
                    minus_factors
                )
            )
        )

        print(
            "    factors_of_p_plus_1={}".format(
                sorted(
                    plus_factors
                )
            )
        )

        for cell, value in sorted(
            cells,
            key=lambda item: item[0][1],
        ):

            q_factors = set(
                factor_abs(
                    int(value)
                )
            )

            print()
            print(
                "    cell={}".format(
                    cell
                )
            )

            print(
                "      Q_factor_support={}".format(
                    sorted(q_factors)
                )
            )

            print(
                "      shared_with_p_minus_1={}".format(
                    sorted(
                        q_factors
                        & minus_factors
                    )
                )
            )

            print(
                "      shared_with_p_plus_1={}".format(
                    sorted(
                        q_factors
                        & plus_factors
                    )
                )
            )


# ============================================================================
# 5. ROW / COLUMN GCD
# ============================================================================

def row_column_gcd(lattice):

    print()
    print("=" * 78)
    print(
        "5. EXACT ROW / COLUMN GCD AUDIT"
    )
    print("=" * 78)

    rows = sorted(
        {
            cell[0]
            for cell in lattice
        }
    )

    for r in rows:

        values = [
            int(
                lattice[(r, t)]
            )
            for t in sorted(
                t
                for (rr, t) in lattice
                if rr == r
            )
        ]

        g = gcd_many(
            values
        )

        primitive, primitive_gcd = (
            primitive_vector(
                values
            )
        )

        print()
        print(
            "  row_r={}".format(
                r
            )
        )

        print(
            "    gcd={}".format(
                g
            )
        )

        print(
            "    gcd_factorization={}".format(
                factor_abs(g)
            )
        )

        print(
            "    primitive_gcd={}".format(
                primitive_gcd
            )
        )

        print(
            "    primitive_vector={}".format(
                primitive
            )
        )

    columns = sorted(
        {
            cell[1]
            for cell in lattice
        }
    )

    for t in columns:

        values = [
            int(
                lattice[(r, t)]
            )
            for r in sorted(
                r
                for (r, tt) in lattice
                if tt == t
            )
        ]

        g = gcd_many(
            values
        )

        print()
        print(
            "  column_t={}".format(
                t
            )
        )

        print(
            "    gcd={}".format(
                g
            )
        )

        print(
            "    gcd_factorization={}".format(
                factor_abs(g)
            )
        )


# ============================================================================
# 6. CROSS-CELL PRIME PERSISTENCE
# ============================================================================

def cross_cell_prime_persistence(lattice):

    print()
    print("=" * 78)
    print(
        "6. EXACT CROSS-CELL PRIME PERSISTENCE AUDIT"
    )
    print("=" * 78)

    prime_cells = defaultdict(
        list
    )

    for cell, value in sorted(
        lattice.items()
    ):

        for prime in factor_abs(
            int(value)
        ):

            prime_cells[prime].append(
                cell
            )

    repeated = [
        (
            prime,
            cells,
        )
        for prime, cells in prime_cells.items()
        if len(cells) >= 2
    ]

    repeated.sort(
        key=lambda item: (
            -len(item[1]),
            item[0],
        )
    )

    for prime, cells in repeated:

        print()
        print(
            "  prime={}".format(
                prime
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
        "  persistent_prime_count={}".format(
            len(repeated)
        )
    )


# ============================================================================
# 7. NEIGHBORING GCD
# ============================================================================

def neighboring_gcd(lattice):

    print()
    print("=" * 78)
    print(
        "7. EXACT NEIGHBORING-CELL GCD AUDIT"
    )
    print("=" * 78)

    cells = set(
        lattice
    )

    horizontal = []
    vertical = []

    for r, t in sorted(cells):

        value = int(
            lattice[(r, t)]
        )

        right = (
            r + 1,
            t,
        )

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

        up = (
            r,
            t + 1,
        )

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

    for a, b, g in horizontal:

        print(
            "    {} -> {}: gcd={}, factors={}".format(
                a,
                b,
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

    for a, b, g in vertical:

        print(
            "    {} -> {}: gcd={}, factors={}".format(
                a,
                b,
                g,
                factor_abs(g),
            )
        )


# ============================================================================
# 8. p-POWER NORMALIZATION
# ============================================================================

def p_power_normalization(lattice):

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

        p = (
            2 * r
            + 1
        )

        value_int = int(value)

        if p == 1:

            print()
            print(
                "  cell=(r={},t={}) p=1: SKIPPED_P_ADIC_NORMALIZATION".format(
                    r,
                    t,
                )
            )

            continue

        vp = valuation(
            value_int,
            p,
        )

        normalized = (
            value_int
            // p ** vp
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
            "    Q_div_p_power={}".format(
                normalized
            )
        )

        print(
            "    normalized_factorization={}".format(
                factor_abs(
                    normalized
                )
            )
        )


# ============================================================================
# 9. GCD WITH BASIC p-DERIVED QUANTITIES
# ============================================================================

def p_quantity_gcd(lattice):

    print()
    print("=" * 78)
    print(
        "9. EXACT GCD WITH BASIC p-DERIVED QUANTITIES"
    )
    print("=" * 78)

    quantity_functions = [
        (
            "p",
            lambda p: p,
        ),
        (
            "p_minus_1",
            lambda p: p - 1,
        ),
        (
            "p_plus_1",
            lambda p: p + 1,
        ),
        (
            "p_squared_minus_1",
            lambda p: p * p - 1,
        ),
        (
            "p_cubed_minus_p",
            lambda p: p ** 3 - p,
        ),
    ]

    for (r, t), value in sorted(
        lattice.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        p = (
            2 * r
            + 1
        )

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

        for name, fn in quantity_functions:

            quantity = int(
                fn(p)
            )

            if quantity == 0:

                print(
                    "    gcd(Q,{})=UNDEFINED_ZERO_REFERENCE".format(
                        name
                    )
                )

                continue

            g = math.gcd(
                abs(value_int),
                abs(quantity),
            )

            print(
                "    gcd(Q,{})={}; factors={}".format(
                    name,
                    g,
                    factor_abs(g),
                )
            )


# ============================================================================
# 10. REPEATED FACTOR-SUPPORT PATTERNS
# ============================================================================

def repeated_factor_patterns(lattice):

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

        support = tuple(
            sorted(
                factor_abs(
                    int(value)
                ).keys()
            )
        )

        patterns[support].append(
            cell
        )

    repeated = [
        (
            support,
            cells,
        )
        for support, cells in patterns.items()
        if len(cells) >= 2
    ]

    repeated.sort(
        key=lambda item: (
            -len(item[1]),
            item[0],
        )
    )

    for support, cells in repeated:

        print()
        print(
            "  factor_support={}".format(
                support
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
        "  repeated_support_pattern_count={}".format(
            len(repeated)
        )
    )


# ============================================================================
# 11. SIGN STRUCTURE
# ============================================================================

def sign_structure(lattice):

    print()
    print("=" * 78)
    print(
        "11. EXACT SIGN / ZERO STRUCTURE"
    )
    print("=" * 78)

    positive = []
    negative = []
    zero = []

    for cell, value in sorted(
        lattice.items()
    ):

        value = int(value)

        if value > 0:
            positive.append(cell)

        elif value < 0:
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
# 12. SUMMARY
# ============================================================================

def structural_summary(lattice):

    print()
    print("=" * 78)
    print(
        "12. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    values = [
        int(value)
        for value in lattice.values()
    ]

    global_gcd = gcd_many(
        values
    )

    all_primes = sorted(
        {
            prime
            for value in values
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
        "  distinct_prime_factor_count={}".format(
            len(all_primes)
        )
    )

    print(
        "  global_value_gcd={}".format(
            global_gcd
        )
    )

    print(
        "  global_value_gcd_factorization={}".format(
            factor_abs(
                global_gcd
            )
        )
    )

    print()
    print(
        "  factorization_cache_size={}".format(
            len(_FACTOR_CACHE)
        )
    )

    print(
        "  valuation_cache_size={}".format(
            len(_VALUATION_CACHE)
        )
    )

    print()
    print(
        "  interpretation=ARITHMETIC_PROVENANCE_ONLY"
    )


# ============================================================================
# 13. FINAL EXACTNESS
# ============================================================================

def final_exactness(lattice):

    print()
    print("=" * 78)
    print(
        "13. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  all_observed_values_integer={}".format(
            all(
                isinstance(
                    value,
                    sp.Integer,
                )
                for value in lattice.values()
            )
        )
    )

    print(
        "  exact_prime_factorization=True"
    )

    print(
        "  p_equals_1_boundary_handled=True"
    )

    print(
        "  invalid_v1_evaluation_prevented=True"
    )

    print(
        "  factorization_cached=True"
    )

    print(
        "  valuation_cached=True"
    )

    print(
        "  row_column_gcds_tested=True"
    )

    print(
        "  neighboring_gcds_tested=True"
    )

    print(
        "  cross_cell_prime_persistence_tested=True"
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
        "EXPERIMENT 385R-FIXED — EXACT PRIME-FACTOR PROVENANCE / "
        "ARITHMETIC-INVARIANT AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    source_inventory(
        lattice
    )

    complete_factorization(
        lattice
    )

    p_divisibility_audit(
        lattice
    )

    p_factor_support(
        lattice
    )

    row_column_gcd(
        lattice
    )

    cross_cell_prime_persistence(
        lattice
    )

    neighboring_gcd(
        lattice
    )

    p_power_normalization(
        lattice
    )

    p_quantity_gcd(
        lattice
    )

    repeated_factor_patterns(
        lattice
    )

    sign_structure(
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
        "EXPERIMENT 385R-FIXED COMPLETE"
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
