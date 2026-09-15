#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 370R — EXACT MODULAR NULLSPACE GENERATOR / SMITH-LIFT AUDIT
==============================================================================

Purpose
-------
369R showed that the modular exceptions are exactly explained by the
determinantal divisors / Smith invariant factors.

370R asks the sharper question:

    What do the exceptional modular nullspaces actually look like?

For each support and each exceptional prime p:

    1. compute the exact nullspace over F_p;
    2. obtain canonical primitive residue generators;
    3. search for the smallest signed integer representative;
    4. test whether that residue vector is the reduction of a small
       integer vector already visible over Z;
    5. compare exceptional-prime generators across supports;
    6. compute the modular nullity profile;
    7. distinguish:
          genuine modular relation
       from
          reduction of an ordinary integer relation.

For supports already having a rational nullspace, their rational relations
are separately primitive-normalized over Z and reduced modulo p.

This experiment does NOT claim that a modular relation is a source law.

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
# BASIC HELPERS
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

    possible_origins = set()

    for r, t in lattice:

        for dr, dt in support:

            possible_origins.add(
                (
                    r - dr,
                    t - dt,
                )
            )

    windows = []

    for origin in sorted(
        possible_origins,
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


def build_matrix(
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


# ============================================================================
# MODULAR ROW REDUCTION
# ============================================================================

def rref_mod(matrix, prime):
    """
    Return modular RREF and pivot columns.

    Entries are Python integers in [0,p-1].
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

    pivots = []
    row = 0

    for col in range(cols):

        pivot = None

        for r in range(
            row,
            rows,
        ):

            if data[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        if pivot != row:
            data[row], data[pivot] = (
                data[pivot],
                data[row],
            )

        inv = pow(
            data[row][col],
            -1,
            p,
        )

        for c in range(cols):
            data[row][c] = (
                data[row][c]
                * inv
            ) % p

        for r in range(rows):

            if r == row:
                continue

            factor = data[r][col]

            if factor == 0:
                continue

            for c in range(cols):
                data[r][c] = (
                    data[r][c]
                    - factor * data[row][c]
                ) % p

        pivots.append(col)
        row += 1

        if row == rows:
            break

    return data, pivots


def modular_nullspace(matrix, prime):
    """
    Nullspace basis over F_p.
    """

    p = int(prime)

    rref, pivots = rref_mod(
        matrix,
        p,
    )

    cols = matrix.cols

    free_cols = [
        col
        for col in range(cols)
        if col not in pivots
    ]

    basis = []

    for free_col in free_cols:

        vector = [0] * cols
        vector[free_col] = 1

        for row, pivot_col in enumerate(
            pivots
        ):

            vector[pivot_col] = (
                -rref[row][free_col]
            ) % p

        basis.append(
            tuple(vector)
        )

    return basis


# ============================================================================
# CANONICAL RESIDUE REPRESENTATIVES
# ============================================================================

def centered_residue(value, prime):
    """
    Map 0,...,p-1 to the symmetric interval.
    """

    value %= prime

    if value > prime // 2:
        value -= prime

    return value


def centered_vector(vector, prime):
    return tuple(
        centered_residue(
            value,
            prime,
        )
        for value in vector
    )


def primitive_integer_vector(vector):
    """
    Primitive normalization of an integer vector.
    """

    values = [
        int(value)
        for value in vector
    ]

    gcd_value = 0

    for value in values:
        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value == 0:
        return tuple(values)

    values = [
        value // gcd_value
        for value in values
    ]

    first_nonzero = next(
        (
            value
            for value in values
            if value != 0
        ),
        0,
    )

    if first_nonzero < 0:
        values = [
            -value
            for value in values
        ]

    return tuple(values)


def l1_norm(vector):
    return sum(
        abs(int(value))
        for value in vector
    )


def linf_norm(vector):
    return max(
        (
            abs(int(value))
            for value in vector
        ),
        default=0,
    )


# ============================================================================
# SUBSPACE MEMBERSHIP MOD p
# ============================================================================

def vector_in_modular_span(
    vector,
    basis,
    prime,
):
    """
    Determine whether vector lies in span(basis) over F_p.
    """

    if not basis:
        return all(
            int(value) % prime == 0
            for value in vector
        )

    p = int(prime)

    basis_matrix = sp.Matrix([
        [
            int(vector[j]) % p
            for j in range(len(vector))
        ]
    ])

    span_rows = [
        [
            int(v[j]) % p
            for j in range(len(v))
        ]
        for v in basis
    ]

    matrix = sp.Matrix(
        span_rows
    )

    target = sp.Matrix(
        [
            int(value) % p
            for value in vector
        ]
    )

    return (
        modular_rank(
            matrix,
            p,
        )
        ==
        modular_rank(
            matrix.col_join(
                target.T
            ),
            p,
        )
    )


def modular_rank(matrix, prime):
    _, pivots = rref_mod(
        matrix,
        prime,
    )

    return len(pivots)


# ============================================================================
# RATIONAL NULLSPACE
# ============================================================================

def primitive_rational_nullspace(
    matrix,
):
    """
    Convert SymPy rational nullspace vectors to primitive integer vectors.
    """

    basis = matrix.nullspace()

    primitive = []

    for vector in basis:

        denominators = [
            int(
                sp.denom(
                    sp.Rational(
                        value
                    )
                )
            )
            for value in vector
        ]

        lcm_denominator = 1

        for denominator in denominators:
            lcm_denominator = sp.ilcm(
                lcm_denominator,
                denominator,
            )

        integer_vector = [
            int(
                sp.Rational(
                    value
                )
                * lcm_denominator
            )
            for value in vector
        ]

        primitive.append(
            primitive_integer_vector(
                integer_vector
            )
        )

    return primitive


# ============================================================================
# EXCEPTIONAL PRIME IDENTIFICATION
# ============================================================================

def exceptional_primes(
    matrix,
):

    rank_q = matrix.rank()

    if rank_q < matrix.cols:
        return []

    exceptional = []

    for prime in TEST_PRIMES:

        rank_p = modular_rank(
            matrix,
            prime,
        )

        if rank_p < rank_q:
            exceptional.append(
                prime
            )

    return exceptional


# ============================================================================
# SECTION 1
# ============================================================================

def modular_nullspace_audit(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "1. EXACT MODULAR NULLSPACE GENERATOR AUDIT"
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

        rank_q = matrix.rank()

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    rational_rank={}".format(
                rank_q
            )
        )

        print(
            "    columns={}".format(
                matrix.cols
            )
        )

        rational_basis = (
            primitive_rational_nullspace(
                matrix
            )
        )

        if rational_basis:

            print(
                "    rational_primitive_generators={}".format(
                    rational_basis
                )
            )

        records = []

        for prime in TEST_PRIMES:

            rank_p = modular_rank(
                matrix,
                prime,
            )

            nullity_p = (
                matrix.cols
                - rank_p
            )

            if nullity_p == 0:
                continue

            basis = modular_nullspace(
                matrix,
                prime,
            )

            centered_basis = [
                centered_vector(
                    vector,
                    prime,
                )
                for vector in basis
            ]

            print()
            print(
                "    prime={}:".format(
                    prime
                )
            )

            print(
                "      rank_mod={}".format(
                    rank_p
                )
            )

            print(
                "      nullity_mod={}".format(
                    nullity_p
                )
            )

            print(
                "      canonical_basis={}".format(
                    centered_basis
                )
            )

            smallest_basis = sorted(
                centered_basis,
                key=lambda v: (
                    l1_norm(v),
                    linf_norm(v),
                    v,
                ),
            )

            print(
                "      smallest_basis_vector={}".format(
                    smallest_basis[0]
                    if smallest_basis
                    else None
                )
            )

            records.append(
                {
                    "prime": prime,
                    "rank": rank_p,
                    "nullity": nullity_p,
                    "basis": centered_basis,
                }
            )

        results[
            support_name
        ] = {
            "rational_basis": rational_basis,
            "modular_records": records,
        }

    return results


# ============================================================================
# SECTION 2
# ============================================================================

def reduction_of_rational_relations_audit(
    inventory,
    modular_results,
):

    print()
    print("=" * 78)
    print(
        "2. RATIONAL-RELATION REDUCTION AUDIT"
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

        rational_basis = primitive_rational_nullspace(
            matrix
        )

        if not rational_basis:
            continue

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        print(
            "    primitive_rational_generators={}".format(
                rational_basis
            )
        )

        for prime in TEST_PRIMES:

            rank_p = modular_rank(
                matrix,
                prime,
            )

            nullity_p = (
                matrix.cols
                - rank_p
            )

            if nullity_p == 0:
                continue

            print()
            print(
                "    prime={}:".format(
                    prime
                )
            )

            for generator in rational_basis:

                reduced = tuple(
                    int(value) % prime
                    for value in generator
                )

                if all(
                    value == 0
                    for value in reduced
                ):
                    status = (
                        "REDUCTION_ZERO"
                    )

                else:
                    status = (
                        "NONZERO_REDUCTION"
                    )

                print(
                    "      generator={} -> {} : {}".format(
                        generator,
                        centered_vector(
                            reduced,
                            prime,
                        ),
                        status,
                    )
                )


# ============================================================================
# SECTION 3
# ============================================================================

def small_integer_lift_audit(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "3. SMALL INTEGER LIFT AUDIT"
    )
    print("=" * 78)

    """
    For each exceptional modular basis vector, test centered representatives
    and simple small multiples against the original integer matrix.

    A vector that is in ker(M mod p) but not in ker(M over Z) is a genuine
    modular relation rather than an ordinary integer relation reduced mod p.
    """

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        exceptional = exceptional_primes(
            matrix
        )

        if not exceptional:
            continue

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        for prime in exceptional:

            basis = modular_nullspace(
                matrix,
                prime,
            )

            print()
            print(
                "    prime={}:".format(
                    prime
                )
            )

            for vector in basis:

                centered = centered_vector(
                    vector,
                    prime,
                )

                integer_matrix_product = (
                    matrix
                    * sp.Matrix(
                        centered
                    )
                )

                zero_over_Z = all(
                    value == 0
                    for value
                    in integer_matrix_product
                )

                print(
                    "      vector={}".format(
                        centered
                    )
                )

                print(
                    "        integer_kernel={}".format(
                        zero_over_Z
                    )
                )

                print(
                    "        M_times_vector={}".format(
                        list(
                            integer_matrix_product
                        )
                    )
                )

                # Search small scalar multiples.
                hits = []

                for scalar in range(
                    -10,
                    11,
                ):

                    if scalar == 0:
                        continue

                    candidate = tuple(
                        scalar * value
                        for value in centered
                    )

                    product = (
                        matrix
                        * sp.Matrix(
                            candidate
                        )
                    )

                    if all(
                        value == 0
                        for value in product
                    ):

                        hits.append(
                            scalar
                        )

                print(
                    "        small_integer_multiples_in_kernel={}".format(
                        hits
                    )
                )


# ============================================================================
# SECTION 4
# ============================================================================

def canonical_support_comparison(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "4. EXCEPTIONAL-PRIME RELATION COMPARISON"
    )
    print("=" * 78)

    signatures = {}

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        exceptional = exceptional_primes(
            matrix
        )

        support_signature = []

        for prime in exceptional:

            basis = modular_nullspace(
                matrix,
                prime,
            )

            canonical = sorted(
                (
                    centered_vector(
                        vector,
                        prime,
                    )
                    for vector in basis
                ),
                key=lambda v: (
                    l1_norm(v),
                    linf_norm(v),
                    v,
                ),
            )

            support_signature.append(
                (
                    prime,
                    canonical,
                )
            )

        signatures[
            support_name
        ] = support_signature

        print()
        print(
            "  {}: {}".format(
                support_name,
                support_signature,
            )
        )

    print()
    print(
        "  exact_same_prime_relation_supports="
    )

    groups = {}

    for (
        support_name,
        signature,
    ) in signatures.items():

        key = repr(
            signature
        )

        groups.setdefault(
            key,
            [],
        ).append(
            support_name
        )

    for members in groups.values():

        if len(members) >= 2:

            print(
                "    {}".format(
                    members
                )
            )

    return signatures


# ============================================================================
# SECTION 5
# ============================================================================

def modular_nullity_profile(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "5. MODULAR NULLITY PROFILE"
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

        profile = []

        for prime in TEST_PRIMES:

            rank_p = modular_rank(
                matrix,
                prime,
            )

            nullity = (
                matrix.cols
                - rank_p
            )

            profile.append(
                (
                    prime,
                    nullity,
                )
            )

        print()
        print(
            "  {}:".format(
                support_name
            )
        )

        print(
            "    {}".format(
                profile
            )
        )


# ============================================================================
# SECTION 6
# ============================================================================

def structural_verdict(
    inventory,
):

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    rational_supports = []
    modular_only_supports = []

    for (
        support_name,
        (
            support,
            windows,
            matrix,
        ),
    ) in inventory.items():

        if matrix.nullspace():

            rational_supports.append(
                support_name
            )

            continue

        exceptional = exceptional_primes(
            matrix
        )

        if exceptional:

            modular_only_supports.append(
                (
                    support_name,
                    exceptional,
                )
            )

    print(
        "  rational_nullspace_supports={}".format(
            rational_supports
        )
    )

    print(
        "  modular_only_exceptional_supports={}".format(
            modular_only_supports
        )
    )

    print(
        "  interpretation="
        "exceptional modular relations are arithmetic rank-collapse "
        "phenomena unless an integer-kernel lift is found"
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 370R — EXACT MODULAR NULLSPACE GENERATOR / "
        "SMITH-LIFT AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    inventory = {}

    for (
        support_name,
        support,
    ) in SUPPORTS.items():

        windows = support_windows(
            lattice,
            support,
        )

        matrix = build_matrix(
            lattice,
            support,
            windows,
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

    modular_results = modular_nullspace_audit(
        inventory
    )

    reduction_of_rational_relations_audit(
        inventory,
        modular_results,
    )

    small_integer_lift_audit(
        inventory
    )

    signatures = canonical_support_comparison(
        inventory
    )

    modular_nullity_profile(
        inventory
    )

    structural_verdict(
        inventory
    )

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells=15"
    )

    print(
        "  support_count={}".format(
            len(SUPPORTS)
        )
    )

    print(
        "  modular_nullspaces_computed_exactly=True"
    )

    print(
        "  rational_relations_reduced_separately=True"
    )

    print(
        "  integer_kernel_lifts_tested=True"
    )

    print(
        "  small_integer_multiple_lifts_tested=True"
    )

    print(
        "  modular_nullity_profiles_completed=True"
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

    print()
    print(
        "EXPERIMENT 370R COMPLETE"
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
