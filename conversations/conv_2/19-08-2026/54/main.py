#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 369R-FIXED — EXACT DETERMINANTAL-DIVISOR / SMITH-INVARIANT /
                         MODULAR-EXCEPTION AUDIT
==============================================================================

Purpose
-------
Experiment 368R found exceptional-prime modular annihilators.

369R determines whether those exceptional primes are simply explained by
integer-matrix determinantal divisors.

For every support:

    * build the exact observed window matrix M over Z;
    * compute rational rank;
    * distinguish rational nullspace from full-column-rank cases;
    * compute determinantal divisors;
    * compute invariant factors;
    * identify exact possible bad primes;
    * independently compute modular ranks;
    * compare with the 368R exceptional-prime profile.

IMPORTANT
---------
If rank_Q(M) < number_of_columns, then there is already a rational
annihilator. In that case there is NO finite "exceptional-prime set"
coming from maximal minors of full column rank.

Such cases are reported as:

    RATIONAL_NULLSPACE

with bad_primes=[]

rather than bad_primes=None.

Only fully observed windows are used.

No missing cells.
No symbolic Z.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact integer arithmetic.
"""


from __future__ import annotations

import itertools
import math
import sys

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
# SUPPORTS
# ============================================================================

SUPPORTS = {
    "horizontal_width3": (
        (0, 0),
        (1, 0),
        (2, 0),
    ),

    "vertical_width3": (
        (0, 0),
        (0, 1),
        (0, 2),
    ),

    "rectangle_2x2": (
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
    ),

    "transport_3": (
        (0, 0),
        (0, 1),
        (1, 1),
    ),

    "lower_transport_3": (
        (0, 0),
        (1, 0),
        (1, 1),
    ),

    "diamond_5": (
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (0, 2),
    ),

    "rectangle_2x3": (
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (0, 2),
        (1, 2),
    ),

    "triangle_6": (
        (0, 0),
        (1, 0),
        (2, 0),
        (0, 1),
        (1, 1),
        (0, 2),
    ),
}


TEST_PRIMES = [
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
]


# ============================================================================
# HELPERS
# ============================================================================

def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r_value, t_value)
            ] = int(value)

    return lattice


def canonical_support(support):
    return tuple(
        sorted(
            support,
            key=lambda z: (
                z[1],
                z[0],
            ),
        )
    )


def support_windows(
    lattice,
    support,
):
    support = canonical_support(
        support
    )

    origins = set()

    for r, t in lattice:

        for dr, dt in support:

            origins.add(
                (
                    r - dr,
                    t - dt,
                )
            )

    windows = []

    for origin in sorted(
        origins,
        key=lambda z: (
            z[1],
            z[0],
        ),
    ):

        r0, t0 = origin

        cells = [
            (
                r0 + dr,
                t0 + dt,
            )
            for dr, dt in support
        ]

        if all(
            cell in lattice
            for cell in cells
        ):
            windows.append(
                origin
            )

    return windows


def build_window_matrix(
    lattice,
    support,
    windows,
):
    support = canonical_support(
        support
    )

    return sp.Matrix([
        [
            sp.Integer(
                lattice[
                    (
                        r0 + dr,
                        t0 + dt,
                    )
                ]
            )
            for dr, dt in support
        ]
        for r0, t0 in windows
    ])


def gcd_list(values):
    result = 0

    for value in values:
        result = math.gcd(
            result,
            abs(int(value)),
        )

    return abs(result)


def factor_integer(n):
    n = int(n)

    if n == 0:
        return {}

    return sp.factorint(
        abs(n)
    )


# ============================================================================
# EXACT MODULAR RANK
# ============================================================================

def modular_rank(matrix, prime):
    """
    Exact rank over F_p by ordinary Gaussian elimination.

    Avoids SymPy DomainMatrix implementation differences.
    """

    p = int(prime)

    data = [
        [
            int(matrix[i, j]) % p
            for j in range(matrix.cols)
        ]
        for i in range(matrix.rows)
    ]

    rows = len(data)
    cols = (
        len(data[0])
        if rows
        else 0
    )

    rank = 0
    pivot_col = 0

    while (
        rank < rows
        and pivot_col < cols
    ):

        pivot = None

        for r in range(
            rank,
            rows,
        ):

            if data[r][pivot_col] % p != 0:
                pivot = r
                break

        if pivot is None:
            pivot_col += 1
            continue

        if pivot != rank:
            data[rank], data[pivot] = (
                data[pivot],
                data[rank],
            )

        inverse = pow(
            data[rank][pivot_col],
            -1,
            p,
        )

        for c in range(
            pivot_col,
            cols,
        ):
            data[rank][c] = (
                data[rank][c]
                * inverse
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = (
                data[r][pivot_col]
                % p
            )

            if factor == 0:
                continue

            for c in range(
                pivot_col,
                cols,
            ):

                data[r][c] = (
                    data[r][c]
                    - factor
                    * data[rank][c]
                ) % p

        rank += 1
        pivot_col += 1

    return rank


# ============================================================================
# MINORS
# ============================================================================

def all_k_minors(
    matrix,
    k,
):
    """
    Return all nonzero k x k minors.

    Each record is:

        (row_indices, column_indices, determinant)
    """

    rows = matrix.rows
    cols = matrix.cols

    if k <= 0 or k > rows or k > cols:
        return []

    records = []

    for row_indices in itertools.combinations(
        range(rows),
        k,
    ):

        for col_indices in itertools.combinations(
            range(cols),
            k,
        ):

            submatrix = matrix.extract(
                row_indices,
                col_indices,
            )

            determinant = int(
                submatrix.det()
            )

            if determinant != 0:

                records.append(
                    (
                        row_indices,
                        col_indices,
                        determinant,
                    )
                )

    return records


# ============================================================================
# DETERMINANTAL DIVISORS / BAD PRIMES
# ============================================================================

def exact_rank_drop_primes(matrix):
    """
    Correctly distinguishes:

        1. RATIONAL_NULLSPACE
           rank_Q(M) < columns.

        2. FULL_COLUMN_RANK
           maximal minors have a nonzero gcd.

    For full column rank, a prime can cause rank collapse modulo p only
    if it divides every maximal minor, i.e. divides the determinantal gcd.

    For rational-nullspace cases, no "bad-prime set" is assigned because
    the relation already exists over Q.
    """

    rank = matrix.rank()
    cols = matrix.cols

    if rank < cols:

        return {
            "mode": "RATIONAL_NULLSPACE",
            "rank": rank,
            "maximal_size": rank,
            "gcd_maximal_minors": None,
            "bad_primes": [],
            "maximal_minors": [],
        }

    maximal_minors = all_k_minors(
        matrix,
        cols,
    )

    gcd_maximal = gcd_list(
        determinant
        for (
            _,
            _,
            determinant,
        ) in maximal_minors
    )

    bad_primes = sorted(
        factor_integer(
            gcd_maximal
        )
    )

    return {
        "mode": "FULL_COLUMN_RANK",
        "rank": rank,
        "maximal_size": cols,
        "gcd_maximal_minors": gcd_maximal,
        "bad_primes": bad_primes,
        "maximal_minors": maximal_minors,
    }


# ============================================================================
# SMITH INVARIANTS
# ============================================================================

def smith_invariants(matrix):
    """
    Compute Smith invariant factors through determinantal divisors.

        Δ_k = gcd(all k x k minors)

    and

        s_1 = Δ_1
        s_k = Δ_k / Δ_{k-1}.

    This is sufficient for the small matrices in this experiment.
    """

    rank = matrix.rank()

    if rank == 0:
        return {
            "rank": 0,
            "determinantal_divisors": [],
            "invariant_factors": [],
        }

    divisors = []

    for k in range(
        1,
        rank + 1,
    ):

        minors = all_k_minors(
            matrix,
            k,
        )

        divisor = gcd_list(
            determinant
            for (
                _,
                _,
                determinant,
            ) in minors
        )

        divisors.append(
            divisor
        )

    invariant_factors = []

    previous = 1

    for divisor in divisors:

        if previous == 0:

            invariant = 0

        else:

            invariant = (
                divisor
                // previous
            )

        invariant_factors.append(
            invariant
        )

        previous = divisor

    return {
        "rank": rank,
        "determinantal_divisors": divisors,
        "invariant_factors": invariant_factors,
    }


# ============================================================================
# SECTION 1
# ============================================================================

def inventory_audit(lattice):

    print()
    print("=" * 78)
    print(
        "1. EXACT SUPPORT / WINDOW INVENTORY"
    )
    print("=" * 78)

    inventory = {}

    for (
        support_name,
        support,
    ) in SUPPORTS.items():

        windows = support_windows(
            lattice,
            support,
        )

        matrix = build_window_matrix(
            lattice,
            support,
            windows,
        )

        print()
        print(
            "  {}:".format(
                support_name
            )
        )

        print(
            "    support={}".format(
                canonical_support(
                    support
                )
            )
        )

        print(
            "    windows={}".format(
                windows
            )
        )

        print(
            "    shape={}".format(
                matrix.shape
            )
        )

        print(
            "    rational_rank={}".format(
                matrix.rank()
            )
        )

        inventory[
            support_name
        ] = (
            canonical_support(
                support
            ),
            windows,
            matrix,
        )

    return inventory


# ============================================================================
# SECTION 2
# ============================================================================

def determinantal_audit(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "2. EXACT DETERMINANTAL-DIVISOR AUDIT"
    )
    print("=" * 78)

    results = {}

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        result = exact_rank_drop_primes(
            matrix
        )

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    shape={}".format(
                matrix.shape
            )
        )

        print(
            "    rational_rank={}".format(
                result["rank"]
            )
        )

        print(
            "    mode={}".format(
                result["mode"]
            )
        )

        if (
            result["mode"]
            == "FULL_COLUMN_RANK"
        ):

            print(
                "    maximal_minor_size={}".format(
                    result[
                        "maximal_size"
                    ]
                )
            )

            print(
                "    gcd_of_maximal_minors={}".format(
                    result[
                        "gcd_maximal_minors"
                    ]
                )
            )

            print(
                "    gcd_factorization={}".format(
                    factor_integer(
                        result[
                            "gcd_maximal_minors"
                        ]
                    )
                )
            )

            print(
                "    exact_exceptional_primes={}".format(
                    result[
                        "bad_primes"
                    ]
                )
            )

        else:

            print(
                "    rational_nullspace=True"
            )

            print(
                "    exact_exceptional_primes=[]"
            )

        results[
            support_name
        ] = result

    return results


# ============================================================================
# SECTION 3
# ============================================================================

def smith_audit(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "3. EXACT SMITH-INVARIANT AUDIT"
    )
    print("=" * 78)

    results = {}

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        result = smith_invariants(
            matrix
        )

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    rank={}".format(
                result["rank"]
            )
        )

        print(
            "    determinantal_divisors={}".format(
                result[
                    "determinantal_divisors"
                ]
            )
        )

        print(
            "    invariant_factors={}".format(
                result[
                    "invariant_factors"
                ]
            )
        )

        print(
            "    invariant_factor_factorizations={}".format(
                [
                    factor_integer(
                        value
                    )
                    for value in result[
                        "invariant_factors"
                    ]
                    if value != 0
                ]
            )
        )

        results[
            support_name
        ] = result

    return results


# ============================================================================
# SECTION 4
# ============================================================================

def modular_exception_audit(
    inventory,
    determinant_results,
):

    print()
    print("=" * 78)
    print(
        "4. EXACT MODULAR EXCEPTION VERIFICATION"
    )
    print("=" * 78)

    comparison = {}

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        rational_rank = matrix.rank()

        result = determinant_results[
            support_name
        ]

        exact_exceptional = set(
            result.get(
                "bad_primes",
                [],
            )
            or []
        )

        observed_exceptions = []

        for prime in TEST_PRIMES:

            rank_mod = modular_rank(
                matrix,
                prime,
            )

            if rank_mod < rational_rank:

                observed_exceptions.append(
                    prime
                )

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    rational_rank={}".format(
                rational_rank
            )
        )

        print(
            "    matrix_columns={}".format(
                matrix.cols
            )
        )

        print(
            "    determinant_mode={}".format(
                result["mode"]
            )
        )

        print(
            "    exact_exceptional_primes={}".format(
                sorted(
                    exact_exceptional
                )
            )
        )

        print(
            "    observed_modular_rank_drop_primes={}".format(
                observed_exceptions
            )
        )

        if (
            result["mode"]
            == "RATIONAL_NULLSPACE"
        ):

            print(
                "    interpretation="
                "rational relation already exists; modular rank drops "
                "are not classified as exceptional-prime-only phenomena"
            )

            agreement = True

        else:

            agreement = (
                sorted(
                    exact_exceptional
                )
                ==
                sorted(
                    observed_exceptions
                )
            )

            print(
                "    exact_match_on_tested_primes={}".format(
                    agreement
                )

            )

        comparison[
            support_name
        ] = {
            "exact_exceptional": sorted(
                exact_exceptional
            ),
            "observed_exceptional": observed_exceptions,
            "match": agreement,
            "mode": result["mode"],
        }

    return comparison


# ============================================================================
# SECTION 5
# ============================================================================

def modular_nullspace_explanation(
    inventory,
    determinant_results,
):

    print()
    print("=" * 78)
    print(
        "5. MODULAR NULLSPACE EXPLANATION"
    )
    print("=" * 78)

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        result = determinant_results[
            support_name
        ]

        if (
            result["mode"]
            != "FULL_COLUMN_RANK"
        ):
            continue

        bad_primes = result[
            "bad_primes"
        ]

        if not bad_primes:
            continue

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    rational_rank={}".format(
                matrix.rank()
            )
        )

        print(
            "    exact_bad_primes={}".format(
                bad_primes
            )
        )

        for prime in TEST_PRIMES:

            if prime not in bad_primes:
                continue

            rank_mod = modular_rank(
                matrix,
                prime,
            )

            print(
                "    prime={}: rank_Q={} -> rank_mod={} "
                "exceptional_rank_collapse=True".format(
                    prime,
                    matrix.rank(),
                    rank_mod,
                )
            )


# ============================================================================
# SECTION 6
# ============================================================================

def minor_prime_profile(
    inventory,
    determinant_results,
):

    print()
    print("=" * 78)
    print(
        "6. MAXIMAL-MINOR PRIME PROFILE"
    )
    print("=" * 78)

    for (
        support_name,
        result,
    ) in determinant_results.items():

        if (
            result["mode"]
            != "FULL_COLUMN_RANK"
        ):
            continue

        gcd_value = result[
            "gcd_maximal_minors"
        ]

        factors = factor_integer(
            gcd_value
        )

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    gcd_of_maximal_minors={}".format(
                gcd_value
            )
        )

        print(
            "    factorization={}".format(
                factors
            )
        )

        matrix = inventory[
            support_name
        ][2]

        k = matrix.cols

        minors = all_k_minors(
            matrix,
            k,
        )

        if not factors:

            print(
                "    prime_profile={}"
            )

            continue

        for prime, exponent in factors.items():

            divisible_count = sum(
                1
                for (
                    _,
                    _,
                    determinant,
                ) in minors
                if (
                    determinant
                    % prime
                    == 0
                )
            )

            print(
                "    prime={}: exponent_in_gcd={}, "
                "divisible_maximal_minors={}/{}".format(
                    prime,
                    exponent,
                    divisible_count,
                    len(minors),
                )
            )


# ============================================================================
# SECTION 7
# ============================================================================

def rational_relation_audit(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "7. RATIONAL VS MODULAR NULLSPACE AUDIT"
    )
    print("=" * 78)

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        rational_rank = matrix.rank()

        rational_nullity = (
            matrix.cols
            - rational_rank
        )

        print()
        print(
            "  {}:".format(
                support_name
            )
        )

        print(
            "    rational_rank={}".format(
                rational_rank
            )
        )

        print(
            "    column_count={}".format(
                matrix.cols
            )
        )

        print(
            "    rational_nullity={}".format(
                rational_nullity
            )
        )

        print(
            "    rational_relation_exists={}".format(
                rational_nullity > 0
            )
        )

        if rational_nullity > 0:

            basis = matrix.nullspace()

            print(
                "    rational_nullspace_basis={}".format(
                    basis
                )
            )


# ============================================================================
# SECTION 8
# ============================================================================

def previous_exception_crosscheck(
    comparison,
):

    print()
    print("=" * 78)
    print(
        "8. 368R CROSS-CHECK"
    )
    print("=" * 78)

    expected = {
        "horizontal_width3": [2],
        "vertical_width3": [2],
        "rectangle_2x2": [2, 3],
        "transport_3": [2],
        "lower_transport_3": [2],
        "diamond_5": [],
        "rectangle_2x3": [],
        "triangle_6": [],
    }

    failures = []

    for (
        support_name,
        expected_primes,
    ) in expected.items():

        actual = comparison[
            support_name
        ][
            "observed_exceptional"
        ]

        match = (
            sorted(actual)
            ==
            sorted(expected_primes)
        )

        print()
        print(
            "  {}:".format(
                support_name
            )
        )

        print(
            "    expected_368R={}".format(
                expected_primes
            )
        )

        print(
            "    observed_now={}".format(
                actual
            )
        )

        print(
            "    match={}".format(
                match
            )
        )

        if not match:
            failures.append(
                support_name
            )

    print()
    print(
        "  crosscheck_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# SECTION 9
# ============================================================================

def structural_summary(
    determinant_results,
    comparison,
    crosscheck_failures,
):

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    exceptional_union = set()

    rational_nullspace_supports = []

    exact_exceptional_supports = []

    for (
        support_name,
        result,
    ) in determinant_results.items():

        if (
            result["mode"]
            == "RATIONAL_NULLSPACE"
        ):

            rational_nullspace_supports.append(
                support_name
            )

        for prime in result.get(
            "bad_primes",
            [],
        ) or []:

            exceptional_union.add(
                prime
            )

            exact_exceptional_supports.append(
                support_name
            )

    print(
        "  exact_exceptional_prime_union={}".format(
            sorted(
                exceptional_union
            )
        )
    )

    print(
        "  rational_nullspace_supports={}".format(
            rational_nullspace_supports
        )
    )

    print(
        "  full_column_rank_exceptional_supports={}".format(
            sorted(
                set(
                    exact_exceptional_supports
                )
            )
        )
    )

    print(
        "  supports_matching_exact_exception_profile={}".format(
            [
                support_name
                for support_name, record
                in comparison.items()
                if record["match"]
            ]
        )
    )

    print(
        "  368R_crosscheck_failures={}".format(
            crosscheck_failures
        )
    )

    if (
        crosscheck_failures
    ):

        interpretation = (
            "MODULAR_PROFILE_MISMATCH"
        )

    elif (
        exceptional_union
        and not rational_nullspace_supports
    ):

        interpretation = (
            "EXCEPTIONAL_PRIME_RANK_COLLAPSE_EXACTLY_IDENTIFIED"
        )

    else:

        interpretation = (
            "MIXED_RATIONAL_AND_EXCEPTIONAL_MODULAR_STRUCTURE"
        )

    print(
        "  interpretation={}".format(
            interpretation
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 369R-FIXED — EXACT DETERMINANTAL-DIVISOR / "
        "SMITH-INVARIANT / MODULAR-EXCEPTION AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    inventory = inventory_audit(
        lattice
    )

    determinant_results = determinantal_audit(
        inventory
    )

    smith_results = smith_audit(
        inventory
    )

    comparison = modular_exception_audit(
        inventory,
        determinant_results,
    )

    modular_nullspace_explanation(
        inventory,
        determinant_results,
    )

    minor_prime_profile(
        inventory,
        determinant_results,
    )

    rational_relation_audit(
        inventory
    )

    crosscheck_failures = (
        previous_exception_crosscheck(
            comparison
        )
    )

    structural_summary(
        determinant_results,
        comparison,
        crosscheck_failures,
    )

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  support_count={}".format(
            len(SUPPORTS)
        )
    )

    print(
        "  exact_determinantal_divisors=True"
    )

    print(
        "  smith_invariants_computed=True"
    )

    print(
        "  exact_exceptional_prime_sets_computed=True"
    )

    print(
        "  modular_368R_crosschecked=True"
    )

    print(
        "  rational_nullspace_distinguished_from_modular=True"
    )

    print(
        "  none_bad_prime_bug_fixed=True"
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
        "  failures={}".format(
            len(crosscheck_failures)
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            len(crosscheck_failures) == 0
        )
    )

    print()
    print(
        "EXPERIMENT 369R-FIXED COMPLETE"
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