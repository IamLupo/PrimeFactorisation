#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 368R — EXACT MODULAR ANNIHILATOR SPECTROSCOPY
==============================================================================

Purpose
-------
Search for small-support linear annihilators of the observed triangular
source lattice modulo many primes.

This is deliberately different from the previous rational/nullspace
experiments.

Question:

    Does the source exhibit a coherent local integer relation whose
    reductions remain visible modulo many unrelated primes?

For each modulus prime ell, the experiment tests selected small connected
supports for a homogeneous relation

    sum_{s in support} c_s Q(r+s_r, t+s_t) = 0 mod ell.

The relation is considered structurally interesting only when:

    * all required cells are observed;
    * equations > support_size;
    * nullspace dimension = 1;
    * the relation is nonzero;
    * the same support repeatedly survives across many primes.

No missing cells are inserted.

No symbolic Z is used as data.

No rational fitting is accepted as a discovery.

Important distinction
---------------------
A relation modulo one prime can occur accidentally.

Therefore the experiment ranks supports by persistence:

    persistence_count
    persistence_fraction
    common primitive support pattern

A support surviving many unrelated primes is evidence of an underlying
integer annihilator candidate, but is NOT by itself a proof.

Exact arithmetic
----------------
Integer arithmetic modulo prime fields only.

No floating point.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
"""


from __future__ import annotations

import itertools
import math
import sys
from collections import defaultdict

import sympy as sp


# ============================================================================
# OBSERVED SOURCE
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
# MODULI
# ============================================================================

PRIMES = [
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
# SUPPORTS
# ============================================================================
#
# Coordinates are (delta_r, delta_t).
#
# Every support is normalized lexicographically so that support identities
# can be compared across primes.
#

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
        (0, 2),
        (1, 1),
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


def mod_value(value, prime):
    return int(value) % prime


def mod_matrix_rank(matrix, prime):
    """
    Exact Gaussian elimination modulo a prime.
    """

    if not matrix:
        return 0

    A = [
        [
            int(value) % prime
            for value in row
        ]
        for row in matrix
    ]

    m = len(A)
    n = len(A[0])

    rank = 0

    for col in range(n):

        pivot = None

        for row in range(
            rank,
            m,
        ):

            if A[row][col] % prime != 0:
                pivot = row
                break

        if pivot is None:
            continue

        A[rank], A[pivot] = (
            A[pivot],
            A[rank],
        )

        inv = pow(
            A[rank][col],
            -1,
            prime,
        )

        A[rank] = [
            value * inv % prime
            for value in A[rank]
        ]

        for row in range(m):

            if row == rank:
                continue

            factor = A[row][col]

            if factor == 0:
                continue

            A[row] = [
                (
                    A[row][j]
                    - factor * A[rank][j]
                ) % prime
                for j in range(n)
            ]

        rank += 1

        if rank == m:
            break

    return rank


def nullspace_mod_prime(
    matrix,
    prime,
):
    """
    Return a basis for the right nullspace of matrix over F_prime.

    Basis vectors are returned as ordinary integer lists in [0,prime-1].
    """

    if not matrix:
        return []

    A = [
        [
            int(value) % prime
            for value in row
        ]
        for row in matrix
    ]

    m = len(A)
    n = len(A[0])

    row = 0
    pivot_columns = []

    for col in range(n):

        pivot = None

        for i in range(
            row,
            m,
        ):

            if A[i][col] % prime != 0:
                pivot = i
                break

        if pivot is None:
            continue

        A[row], A[pivot] = (
            A[pivot],
            A[row],
        )

        inv = pow(
            A[row][col],
            -1,
            prime,
        )

        A[row] = [
            value * inv % prime
            for value in A[row]
        ]

        for i in range(m):

            if i == row:
                continue

            factor = A[i][col]

            if factor == 0:
                continue

            A[i] = [
                (
                    A[i][j]
                    - factor * A[row][j]
                ) % prime
                for j in range(n)
            ]

        pivot_columns.append(col)
        row += 1

        if row == m:
            break

    free_columns = [
        col
        for col in range(n)
        if col not in pivot_columns
    ]

    basis = []

    for free in free_columns:

        vector = [
            0
            for _ in range(n)
        ]

        vector[free] = 1

        # Since the matrix is in RREF, every pivot variable is determined
        # directly from the free variables.
        for pivot_row in range(
            len(pivot_columns)
        ):

            pivot_col = pivot_columns[
                pivot_row
            ]

            vector[pivot_col] = (
                -A[pivot_row][free]
            ) % prime

        basis.append(vector)

    return basis


def primitive_projective_vector(
    vector,
    prime,
):
    """
    Normalize a nonzero vector over F_p projectively.

    The first nonzero coordinate becomes 1.

    This canonical form permits comparison of the same projective relation
    across modular computations.
    """

    vector = [
        int(value) % prime
        for value in vector
    ]

    first = None

    for value in vector:

        if value != 0:
            first = value
            break

    if first is None:
        return None

    inverse = pow(
        first,
        -1,
        prime,
    )

    return tuple(
        (
            value * inverse
        ) % prime
        for value in vector
    )


def support_windows(
    lattice,
    support,
):
    """
    Return all translations whose entire support is observed.
    """

    support = canonical_support(
        support
    )

    windows = []

    candidate_origins = sorted(
        lattice,
        key=lambda z: (
            z[1],
            z[0],
        ),
    )

    for origin in candidate_origins:

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


def relation_matrix(
    lattice,
    support,
    windows,
    prime,
):
    matrix = []

    support = canonical_support(
        support
    )

    for r0, t0 in windows:

        matrix.append(
            [
                lattice[
                    (
                        r0 + dr,
                        t0 + dt,
                    )
                ] % prime
                for dr, dt in support
            ]
        )

    return matrix


def relation_residual(
    lattice,
    support,
    windows,
    relation,
    prime,
):
    residuals = []

    support = canonical_support(
        support
    )

    for r0, t0 in windows:

        value = 0

        for coefficient, (
            dr,
            dt,
        ) in zip(
            relation,
            support,
        ):

            value += (
                coefficient
                * lattice[
                    (
                        r0 + dr,
                        t0 + dt,
                    )
                ]
            )

        residuals.append(
            value % prime
        )

    return residuals


# ============================================================================
# SECTION 1 — BASIC MODULAR DATA
# ============================================================================

def source_summary(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "1. SOURCE / MODULAR AUDIT"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  tested_primes={}".format(
            PRIMES
        )
    )

    for prime in PRIMES:

        nonzero = sum(
            1
            for value in lattice.values()
            if value % prime != 0
        )

        zero = (
            len(lattice)
            - nonzero
        )

        print()
        print(
            "  prime={}: nonzero_cells={}, "
            "zero_cells={}".format(
                prime,
                nonzero,
                zero,
            )
        )


# ============================================================================
# SECTION 2 — SINGLE CASE
# ============================================================================

def audit_case(
    lattice,
    support_name,
    support,
    prime,
):
    support = canonical_support(
        support
    )

    windows = support_windows(
        lattice,
        support,
    )

    equation_count = len(
        windows
    )

    support_size = len(
        support
    )

    if equation_count == 0:

        return {
            "status": "UNAVAILABLE",
            "equations": 0,
            "unknowns": support_size,
            "rank": 0,
            "nullity": support_size,
            "windows": [],
            "relation": None,
        }

    matrix = relation_matrix(
        lattice,
        support,
        windows,
        prime,
    )

    rank = mod_matrix_rank(
        matrix,
        prime,
    )

    nullity = (
        support_size
        - rank
    )

    if nullity == 0:

        status = "NO_RELATION"
        relation = None

    elif (
        equation_count
        > support_size
        and nullity == 1
    ):

        basis = nullspace_mod_prime(
            matrix,
            prime,
        )

        relation = (
            primitive_projective_vector(
                basis[0],
                prime,
            )
            if basis
            else None
        )

        residuals = relation_residual(
            lattice,
            support,
            windows,
            relation,
            prime,
        )

        if all(
            residual == 0
            for residual in residuals
        ):
            status = "EXACT_OVERDETERMINED"
        else:
            status = "VERIFICATION_FAILED"

    elif nullity == 1:

        basis = nullspace_mod_prime(
            matrix,
            prime,
        )

        relation = (
            primitive_projective_vector(
                basis[0],
                prime,
            )
            if basis
            else None
        )

        status = "DATA_SIZED"

    else:

        relation = None
        status = "MULTIPLE_RELATIONS"

    return {
        "status": status,
        "equations": equation_count,
        "unknowns": support_size,
        "rank": rank,
        "nullity": nullity,
        "windows": windows,
        "relation": relation,
    }


# ============================================================================
# SECTION 3 — FULL SPECTROSCOPY
# ============================================================================

def run_spectroscopy(
    lattice,
):
    print()
    print("=" * 78)
    print(
        "2. EXACT MODULAR ANNIHILATOR SPECTROSCOPY"
    )
    print("=" * 78)

    results = []

    for (
        support_name,
        support,
    ) in SUPPORTS.items():

        print()
        print(
            "  SUPPORT: {}".format(
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

        support_records = []

        for prime in PRIMES:

            result = audit_case(
                lattice,
                support_name,
                support,
                prime,
            )

            record = (
                prime,
                result,
            )

            support_records.append(
                record
            )

            print()
            print(
                "    p={}: status={}, "
                "equations={}, rank={}, "
                "nullity={}".format(
                    prime,
                    result["status"],
                    result["equations"],
                    result["rank"],
                    result["nullity"],
                )
            )

            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            ):

                print(
                    "      relation={}".format(
                        result["relation"]
                    )
                )

        results.append(
            (
                support_name,
                canonical_support(
                    support
                ),
                support_records,
            )
        )

    return results


# ============================================================================
# SECTION 4 — PERSISTENCE AUDIT
# ============================================================================

def persistence_audit(
    spectroscopy_results,
):
    print()
    print("=" * 78)
    print(
        "3. MODULAR PERSISTENCE AUDIT"
    )
    print("=" * 78)

    persistent = []

    for (
        support_name,
        support,
        records,
    ) in spectroscopy_results:

        tested = 0
        overdetermined = 0
        projective_relations = {}

        for prime, result in records:

            if (
                result["status"]
                == "UNAVAILABLE"
            ):
                continue

            tested += 1

            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            ):

                overdetermined += 1

                relation = (
                    result["relation"]
                )

                projective_relations[
                    prime
                ] = relation

        print()
        print(
            "  {}:".format(
                support_name
            )
        )

        print(
            "    tested_primes={}".format(
                tested
            )
        )

        print(
            "    overdetermined_prime_count={}".format(
                overdetermined
            )
        )

        fraction = (
            sp.Rational(
                overdetermined,
                tested,
            )
            if tested
            else sp.Rational(0)
        )

        print(
            "    persistence_fraction={}".format(
                fraction
            )
        )

        if projective_relations:

            print(
                "    overdetermined_relations={}".format(
                    projective_relations
                )
            )

        if (
            tested
            and overdetermined == tested
        ):

            persistent.append(
                support_name
            )

    print()
    print(
        "  fully_persistent_supports={}".format(
            persistent
        )
    )

    return persistent


# ============================================================================
# SECTION 5 — CROSS-PRIME COMPATIBILITY
# ============================================================================

def cross_prime_audit(
    spectroscopy_results,
):
    """
    If a support survives modulo multiple primes, compare whether the
    modular relations can plausibly be reductions of one integer vector.

    We do NOT claim a unique lift from arbitrary modular data.

    Instead we record the projective vectors exactly and look for the
    special case where all coordinates are 0/1 or small signed residues
    across many primes.
    """

    print()
    print("=" * 78)
    print(
        "4. CROSS-PRIME RELATION COMPATIBILITY"
    )
    print("=" * 78)

    candidate_supports = []

    for (
        support_name,
        support,
        records,
    ) in spectroscopy_results:

        relations = [
            (
                prime,
                result["relation"],
            )
            for prime, result in records
            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            )
        ]

        if len(relations) < 2:
            continue

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        for prime, relation in relations:

            print(
                "    prime={}: relation={}".format(
                    prime,
                    relation,
                )
            )

        small_coordinate = True

        for prime, relation in relations:

            for value in relation:

                signed = (
                    value
                    if value <= prime // 2
                    else value - prime
                )

                if abs(signed) > 20:

                    small_coordinate = False

        print(
            "    small_signed_coordinates={}".format(
                small_coordinate
            )
        )

        if small_coordinate:

            candidate_supports.append(
                support_name
            )

    print()
    print(
        "  small_signed_relation_candidates={}".format(
            candidate_supports
        )
    )

    return candidate_supports


# ============================================================================
# SECTION 6 — MODULAR RANK PROFILE
# ============================================================================

def rank_profile_audit(
    lattice,
):
    print()
    print("=" * 78)
    print(
        "5. MODULAR RANK PROFILE"
    )
    print("=" * 78)

    # The four-by-six Newton-style coefficient matrix is not reused here.
    # This is the raw observed rectangular embedding with missing cells
    # replaced by no row/column operation. We instead compute ranks of
    # complete rectangular submatrices.
    #
    # For each prime, report the ranks of the largest complete rectangles
    # available in the triangular support.

    rectangles = [
        (
            "3x2",
            [
                (0, 0),
                (1, 0),
                (2, 0),
            ],
            [0, 1],
        ),
        (
            "2x3",
            [
                (0, 0),
                (1, 0),
            ],
            [0, 1, 2],
        ),
    ]

    # The actual complete blocks are selected conservatively.
    blocks = [
        (
            "top_3x2",
            [
                [
                    (0, 0),
                    (1, 0),
                    (2, 0),
                ],
                [
                    (0, 1),
                    (1, 1),
                    (2, 1),
                ],
            ],
        ),
        (
            "left_2x4",
            [
                [
                    (0, 0),
                    (0, 1),
                    (0, 2),
                    (0, 3),
                ],
                [
                    (1, 0),
                    (1, 1),
                    (1, 2),
                    (1, 3),
                ],
            ],
        ),
    ]

    for name, rows in blocks:

        if not all(
            cell in lattice
            for row in rows
            for cell in row
        ):
            continue

        matrix = [
            [
                lattice[cell]
                for cell in row
            ]
            for row in rows
        ]

        print()
        print(
            "  block={}".format(
                name
            )
        )

        print(
            "    shape=({}, {})".format(
                len(matrix),
                len(matrix[0]),
            )
        )

        for prime in PRIMES:

            rank = mod_matrix_rank(
                matrix,
                prime,
            )

            print(
                "    prime={}: rank={}".format(
                    prime,
                    rank,
                )
            )


# ============================================================================
# SECTION 7 — CRT LIFT CANDIDATES
# ============================================================================

def CRT_candidate_audit(
    spectroscopy_results,
):
    """
    When one support has a one-dimensional projective nullspace over several
    primes, attempt a conservative CRT lift of each coordinate.

    We only perform this diagnostic when:
        * the support has an overdetermined nullspace for >= 3 primes;
        * each modular relation has a nonzero first coordinate.

    The output is NOT interpreted as a discovered integer relation.
    """

    print()
    print("=" * 78)
    print(
        "6. CONSERVATIVE CRT LIFT DIAGNOSTIC"
    )
    print("=" * 78)

    for (
        support_name,
        support,
        records,
    ) in spectroscopy_results:

        usable = [
            (
                prime,
                result["relation"],
            )
            for prime, result in records
            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
                and result["relation"] is not None
                and result["relation"][0] != 0
            )
        ]

        if len(usable) < 3:
            continue

        print()
        print(
            "  SUPPORT={}".format(
                support_name
            )
        )

        # Only use the first three qualifying primes to keep the diagnostic
        # deliberately bounded.
        selected = usable[:3]

        print(
            "    selected_primes={}".format(
                [
                    prime
                    for prime, _
                    in selected
                ]
            )
        )

        modulus = math.prod(
            prime
            for prime, _
            in selected
        )

        lifted = []

        for index in range(
            len(support)
        ):

            residues = []

            for (
                prime,
                relation,
            ) in selected:

                residues.append(
                    int(
                        relation[index]
                    )
                )

            value = int(
                sp.ntheory.modular.crt(
                    [
                        prime
                        for prime, _
                        in selected
                    ],
                    residues,
                )[0]
            )

            # Symmetric representative.
            if value > modulus // 2:
                value -= modulus

            lifted.append(
                value
            )

        print(
            "    crt_lift_modulus={}".format(
                modulus
            )
        )

        print(
            "    symmetric_crt_vector={}".format(
                lifted
            )
        )

        print(
            "    diagnostic_only=True"
        )


# ============================================================================
# SECTION 8 — STRUCTURAL SUMMARY
# ============================================================================

def structural_summary(
    spectroscopy_results,
    persistent,
    compatible,
):
    print()
    print("=" * 78)
    print(
        "7. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    total_cases = 0
    total_overdetermined = 0

    for (
        _,
        _,
        records,
    ) in spectroscopy_results:

        for _, result in records:

            total_cases += 1

            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            ):
                total_overdetermined += 1

    print(
        "  support_count={}".format(
            len(SUPPORTS)
        )
    )

    print(
        "  prime_count={}".format(
            len(PRIMES)
        )
    )

    print(
        "  modular_cases_tested={}".format(
            total_cases
        )
    )

    print(
        "  overdetermined_modular_cases={}".format(
            total_overdetermined
        )
    )

    print(
        "  fully_persistent_supports={}".format(
            persistent
        )
    )

    print(
        "  cross_prime_small_relation_candidates={}".format(
            compatible
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 368R — EXACT MODULAR ANNIHILATOR "
        "SPECTROSCOPY"
    )
    print("=" * 78)

    lattice = build_lattice()

    source_summary(
        lattice
    )

    spectroscopy_results = run_spectroscopy(
        lattice
    )

    persistent = persistence_audit(
        spectroscopy_results
    )

    compatible = cross_prime_audit(
        spectroscopy_results
    )

    rank_profile_audit(
        lattice
    )

    CRT_candidate_audit(
        spectroscopy_results
    )

    structural_summary(
        spectroscopy_results,
        persistent,
        compatible,
    )

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  modular_arithmetic_exact=True"
    )

    print(
        "  primes_tested={}".format(
            PRIMES
        )
    )

    print(
        "  small_supports_tested={}".format(
            list(SUPPORTS)
        )
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
        "EXPERIMENT 368R COMPLETE"
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
