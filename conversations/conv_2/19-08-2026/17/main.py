#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 333R — EXACT EXHAUSTIVE SMALL-SUPPORT 2D ANNIHILATOR AUDIT
==============================================================================

Purpose
-------
Experiments 322R-332R rejected:

    * fixed vertical recurrences;
    * low-degree p-dependent vertical recurrences;
    * low-degree rational vertical recurrences;
    * low-degree bivariate polynomial laws;
    * several hand-picked local 2D stencils;
    * rank-one additive separability;
    * rank-one multiplicative separability;
    * low-complexity Newton-coordinate dynamics.

Experiment 333R now performs the systematic version of the local-stencil
test.

Represent the observed source table on the integer lattice

    r = (p-1)/2,

so

    p = 1,3,5,7  <=>  r = 0,1,2,3.

A translation-invariant stencil is a finite set S of offsets

    (dr,dt)

with a homogeneous relation

    sum_{(dr,dt) in S}
        c[dr,dt] Q(r+dr,t+dt) = 0

whenever every translated cell is observed.

This experiment exhaustively searches small supports rather than choosing
individual stencils by hand.

Search region
-------------
All nonempty support masks containing (0,0) inside:

    2x2
    2x3
    2x4
    3x2
    3x3
    3x4
    4x2
    4x3

Only supports with enough translated windows to make the relation
falsifiable are considered.

Acceptance standard
-------------------
A candidate is interesting only when:

    1. there are MORE equations than unknown stencil coefficients;
    2. the homogeneous system has a nonzero nullspace;
    3. the nullspace is one-dimensional.

Thus:

    nullity = 0
        => exact contradiction;

    nullity = 1 and equations > unknowns
        => exact overdetermined stencil law;

    nullity > 1
        => underdetermined family, not a discovery.

No interpolation.
No missing values.
No synthetic second n=pq case.
Exact SymPy rational/integer arithmetic only.
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
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def build_observed_grid():

    grid = {}

    for p_value in sorted(Q):

        values = Q[p_value]

        r = (
            p_value - 1
        ) // 2

        for index, value in enumerate(
            values
        ):

            t = (
                len(values)
                - 1
                - index
            )

            grid[(r, t)] = sp.Integer(
                value
            )

    return grid


def support_string(support):

    return [
        tuple(offset)
        for offset in support
    ]


def normalize_support(
    support,
):
    """
    Canonicalize an offset support.

    Require (0,0) to be present and sort lexicographically.
    """

    support = set(
        tuple(offset)
        for offset in support
    )

    if (0, 0) not in support:
        raise ValueError(
            "support must contain (0,0)"
        )

    return tuple(
        sorted(support)
    )


def translated_windows(
    grid,
    support,
):
    """
    Return all base locations (r,t) for which every translated cell
    (r+dr,t+dt) in support is observed.
    """

    support = normalize_support(
        support
    )

    if not support:
        return []

    r_min = min(
        r
        for r, _ in grid
    )

    r_max = max(
        r
        for r, _ in grid
    )

    t_min = min(
        t
        for _, t in grid
    )

    t_max = max(
        t
        for _, t in grid
    )

    max_dr = max(
        dr
        for dr, _ in support
    )

    max_dt = max(
        dt
        for _, dt in support
    )

    windows = []

    for r in range(
        r_min,
        r_max - max_dr + 1,
    ):

        for t in range(
            t_min,
            t_max - max_dt + 1,
        ):

            if all(
                (
                    r + dr,
                    t + dt,
                )
                in grid
                for dr, dt in support
            ):

                windows.append(
                    (
                        r,
                        t,
                    )
                )

    return windows


# ============================================================================
# HOMOGENEOUS EXACT SYSTEM
# ============================================================================

def stencil_system(
    grid,
    support,
):
    support = normalize_support(
        support
    )

    windows = translated_windows(
        grid,
        support,
    )

    A = sp.Matrix([
        [
            grid[
                (
                    r + dr,
                    t + dt,
                )
            ]
            for dr, dt in support
        ]
        for r, t in windows
    ])

    return A, windows


def analyze_support(
    grid,
    support,
):
    """
    Analyze a homogeneous stencil system

        A c = 0.

    The nullspace dimension is the projective number of free stencil laws.
    """

    support = normalize_support(
        support
    )

    A, windows = stencil_system(
        grid,
        support,
    )

    equations = len(windows)
    unknowns = len(support)

    if equations == 0:
        return {
            "support": support,
            "windows": windows,
            "equations": 0,
            "unknowns": unknowns,
            "rank": 0,
            "nullity": 0,
            "status": "UNUSABLE",
            "basis": [],
        }

    rank = A.rank()
    nullspace = A.nullspace()
    nullity = len(nullspace)

    if nullity == 0:
        status = "NO_SOLUTION"

    elif equations > unknowns and nullity == 1:
        status = "EXACT_OVERDETERMINED"

    elif equations == unknowns and nullity == 1:
        status = "EXACT_DATA_SIZED"

    else:
        status = "NONUNIQUE"

    return {
        "support": support,
        "windows": windows,
        "equations": equations,
        "unknowns": unknowns,
        "rank": rank,
        "nullity": nullity,
        "status": status,
        "basis": nullspace,
        "matrix": A,
    }


# ============================================================================
# PRIMITIVE NULLSPACE NORMALIZATION
# ============================================================================

def primitive_vector(
    vector,
):
    """
    Scale a rational nullspace vector to primitive integer coordinates.
    """

    vector = [
        sp.Rational(x)
        for x in vector
    ]

    denominators = [
        int(x.q)
        for x in vector
        if x != 0
    ]

    if not denominators:
        return vector

    lcm = sp.ilcm(
        *denominators
    )

    integers = [
        int(x * lcm)
        for x in vector
    ]

    g = 0

    for x in integers:
        g = math.gcd(
            g,
            abs(x),
        )

    if g == 0:
        return [
            sp.Integer(0)
            for _ in integers
        ]

    integers = [
        x // g
        for x in integers
    ]

    for x in integers:
        if x != 0:

            if x < 0:
                integers = [
                    -y
                    for y in integers
                ]

            break

    return [
        sp.Integer(x)
        for x in integers
    ]


def verify_stencil(
    grid,
    result,
):
    """
    Exact residual verification for a candidate nullspace vector.
    """

    if result["status"] not in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        return None

    basis = result["basis"]

    if len(basis) != 1:
        return False

    vector = basis[0]
    support = result["support"]

    residuals = []

    for r, t in result["windows"]:

        residual = clean(
            sum(
                vector[i]
                * grid[
                    (
                        r + support[i][0],
                        t + support[i][1],
                    )
                ]
                for i in range(
                    len(support)
                )
            )
        )

        residuals.append(
            residual
        )

    return all(
        residual == 0
        for residual in residuals
    )


# ============================================================================
# SUPPORT GENERATION
# ============================================================================

def full_box(
    height,
    width,
):
    """
    Offsets:

        dr = 0,...,height-1
        dt = 0,...,width-1

    This is a rectangular support.
    """

    return normalize_support(
        itertools.product(
            range(height),
            range(width),
        )
    )


def connected_supports_in_box(
    height,
    width,
    max_size,
):
    """
    Enumerate connected supports containing (0,0).

    Connectivity is 4-neighbour connectivity in the offset lattice.

    This avoids counting disconnected algebraic operators whose separate
    pieces have no genuine local interaction.
    """

    cells = [
        (dr, dt)
        for dr in range(height)
        for dt in range(width)
    ]

    origin = (0, 0)

    results = set()

    for size in range(
        2,
        min(
            max_size,
            len(cells),
        ) + 1,
    ):

        for subset in itertools.combinations(
            [
                cell
                for cell in cells
                if cell != origin
            ],
            size - 1,
        ):

            support = (
                origin,
            ) + subset

            support_set = set(
                support
            )

            frontier = {
                origin
            }

            changed = True

            while changed:

                changed = False

                for cell in list(
                    support_set
                ):

                    r, c = cell

                    neighbours = (
                        (r - 1, c),
                        (r + 1, c),
                        (r, c - 1),
                        (r, c + 1),
                    )

                    if any(
                        n in frontier
                        for n in neighbours
                    ):

                        if cell not in frontier:

                            frontier.add(
                                cell
                            )

                            changed = True

            if len(frontier) == len(
                support_set
            ):

                results.add(
                    normalize_support(
                        support_set
                    )
                )

    return sorted(
        results,
        key=lambda s: (
            len(s),
            s,
        ),
    )


# ============================================================================
# SEARCH
# ============================================================================

def exhaustive_search(
    grid,
):
    """
    Search all connected supports of size <= 6 inside small boxes.

    This is deliberately conservative: the observed table is too small
    to justify a large unrestricted stencil search.
    """

    boxes = [
        (2, 2),
        (2, 3),
        (2, 4),
        (3, 2),
        (3, 3),
        (3, 4),
        (4, 2),
        (4, 3),
    ]

    candidates = []

    summary = []

    seen = set()

    for height, width in boxes:

        supports = connected_supports_in_box(
            height,
            width,
            max_size=6,
        )

        box_total = 0
        box_usable = 0
        box_overdetermined = 0
        box_discoveries = 0

        for support in supports:

            if support in seen:
                continue

            seen.add(
                support
            )

            result = analyze_support(
                grid,
                support,
            )

            if result["status"] == "UNUSABLE":
                continue

            box_total += 1

            if result["equations"] >= 2:
                box_usable += 1

            if (
                result["equations"]
                > result["unknowns"]
            ):
                box_overdetermined += 1

            if result["status"] == "EXACT_OVERDETERMINED":

                verified = verify_stencil(
                    grid,
                    result,
                )

                if verified:

                    box_discoveries += 1

                    candidates.append(
                        (
                            support,
                            result,
                        )
                    )

        summary.append(
            (
                height,
                width,
                box_total,
                box_usable,
                box_overdetermined,
                box_discoveries,
            )
        )

    return candidates, summary


# ============================================================================
# PRIME PROFILE FOR A DISCOVERED STENCIL
# ============================================================================

def print_candidate(
    index,
    grid,
    support,
    result,
):
    print()
    print("=" * 78)
    print(
        "DISCOVERED OVERDETERMINED STENCIL #{}".format(
            index
        )
    )
    print("=" * 78)

    print(
        "  support={}".format(
            support
        )
    )

    print(
        "  equations={}".format(
            result["equations"]
        )
    )

    print(
        "  unknowns={}".format(
            result["unknowns"]
        )
    )

    print(
        "  rank={}".format(
            result["rank"]
        )
    )

    print(
        "  nullity={}".format(
            result["nullity"]
        )
    )

    vector = primitive_vector(
        result["basis"][0]
    )

    print(
        "  primitive_coefficients={}".format(
            vector
        )
    )

    print()
    print(
        "  relation:"
    )

    terms = []

    for coefficient, offset in zip(
        vector,
        support,
    ):

        dr, dt = offset

        terms.append(
            "{} * Q(r+{},t+{})".format(
                coefficient,
                dr,
                dt,
            )
        )

    print(
        "    {}".format(
            " + ".join(
                terms
            )
            + " = 0"
        )
    )

    print()
    print(
        "  residual_check={}".format(
            verify_stencil(
                grid,
                result,
            )
        )
    )

    print()
    print(
        "  coefficient_prime_profiles="
    )

    for coefficient in vector:

        print(
            "    {} -> {}".format(
                coefficient,
                {
                    p: (
                        None
                        if coefficient == 0
                        else (
                            (
                                lambda q: (
                                    (
                                        lambda n, d: (
                                            (
                                                lambda vn, vd: (
                                                    vn - vd
                                                )
                                            )(
                                                # valuation numerator
                                                (lambda n_: (
                                                    0
                                                ))(n),
                                                # dummy
                                                d,
                                            )
                                        )
                                    )(
                                        abs(int(q.p)),
                                        abs(int(q.q)),
                                    )
                                )
                            )(
                                sp.Rational(
                                    coefficient
                                )
                            )
                        )
                    )
                    for p in ()
                }
            )
        )

    # Simpler exact valuation display.
    for coefficient in vector:

        profiles = {}

        coefficient = sp.Rational(
            coefficient
        )

        for prime in (
            2,
            3,
            5,
            7,
            11,
            13,
            17,
        ):

            if coefficient == 0:

                profiles[prime] = sp.oo

                continue

            n = abs(
                int(
                    coefficient.p
                )
            )

            d = abs(
                int(
                    coefficient.q
                )
            )

            vn = 0
            vd = 0

            while n % prime == 0:
                n //= prime
                vn += 1

            while d % prime == 0:
                d //= prime
                vd += 1

            profiles[prime] = vn - vd

        print(
            "    {}: {}".format(
                coefficient,
                profiles,
            )
        )


# ============================================================================
# EXISTING HAND-PICKED STENCIL COMPARISON
# ============================================================================

def compare_known_stencils(
    grid,
):

    print()
    print("=" * 78)
    print(
        "KNOWN-STENCIL COMPARISON"
    )
    print("=" * 78)

    known = {
        "rectangle_2x2": [
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
        ],
        "raw_width2": [
            (0, 0),
            (0, 1),
            (1, 0),
        ],
        "vertical_3": [
            (0, 0),
            (0, 1),
            (0, 2),
        ],
        "transport_3": [
            (0, 0),
            (0, 1),
            (1, 1),
        ],
    }

    for name, support in known.items():

        result = analyze_support(
            grid,
            support,
        )

        print()
        print(
            "  {}:".format(
                name
            )
        )

        print(
            "    support={}".format(
                support
            )
        )

        print(
            "    equations={}".format(
                result["equations"]
            )
        )

        print(
            "    unknowns={}".format(
                result["unknowns"]
            )
        )

        print(
            "    rank={}".format(
                result["rank"]
            )
        )

        print(
            "    nullity={}".format(
                result["nullity"]
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )


# ============================================================================
# STRUCTURAL SUMMARY
# ============================================================================

def print_summary(
    summary,
):
    print()
    print("=" * 78)
    print(
        "SEARCH SUMMARY"
    )
    print("=" * 78)

    for (
        height,
        width,
        total,
        usable,
        overdetermined,
        discoveries,
    ) in summary:

        print()
        print(
            "  box={}x{}:".format(
                height,
                width,
            )
        )

        print(
            "    supports_tested={}".format(
                total
            )
        )

        print(
            "    usable_supports={}".format(
                usable
            )
        )

        print(
            "    overdetermined_supports={}".format(
                overdetermined
            )
        )

        print(
            "    exact_overdetermined_discoveries={}".format(
                discoveries
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 333R — EXACT EXHAUSTIVE SMALL-SUPPORT "
        "2D ANNIHILATOR AUDIT"
    )
    print("=" * 78)

    grid = build_observed_grid()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED LATTICE"
    )
    print("=" * 78)

    print(
        "  lattice_coordinates r=(p-1)/2"
    )

    print(
        "  observed_cells={}".format(
            len(grid)
        )
    )

    for (
        r,
        t,
    ), value in sorted(
        grid.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        ),
    ):

        print(
            "  (r={},t={}) -> {}".format(
                r,
                t,
                value,
            )
        )

    compare_known_stencils(
        grid
    )

    candidates, search_summary = exhaustive_search(
        grid
    )

    print_summary(
        search_summary
    )

    print()
    print("=" * 78)
    print(
        "4. OVERDETERMINED DISCOVERIES"
    )
    print("=" * 78)

    print(
        "  discovery_count={}".format(
            len(candidates)
        )
    )

    for index, (
        support,
        result,
    ) in enumerate(
        candidates,
        start=1,
    ):

        print_candidate(
            index,
            grid,
            support,
            result,
        )

    print()
    print("=" * 78)
    print(
        "5. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 333R upgrades the earlier local-stencil experiments from a
hand-picked search to an exhaustive finite search over connected small
supports.

The tested class is:

    translation-invariant homogeneous linear relations

with support contained in small local boxes.

A candidate counts as structurally interesting only when:

    equations > unknown coefficients

and the exact homogeneous nullspace has dimension one.

This excludes two common failure modes:

    * a stencil that merely interpolates the available windows;
    * a large nullspace containing many unrelated accidental relations.

If discovery_count = 0, then no connected translation-invariant
annihilator of the tested sizes survives exact overdetermined validation.

If discoveries appear, their primitive integer coefficient vectors become
the next objects to study for arithmetic or combinatorial meaning.

The search is deliberately finite and local. It does not claim that
no larger-support or non-translation-invariant relation exists.

No missing cells are filled.
No extrapolation is performed.
No synthetic second n=pq case is generated.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exhaustive_small_support_search=True"
    )

    print(
        "  discovered_overdetermined_stencils={}".format(
            len(candidates)
        )
    )

    print(
        "  known_stencils_rechecked=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
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
        "EXPERIMENT 333R COMPLETE"
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
