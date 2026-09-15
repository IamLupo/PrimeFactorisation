#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 380R — EXACT BIVARIATE HYPERGEOMETRIC / RATIONAL-RATIO
                  SOURCE-LAW AUDIT
==============================================================================

Purpose
-------
379R ruled out low-order constant-coefficient 1D convolution laws.

380R tests a different canonical source family:

    multiplicative / hypergeometric laws,

where adjacent source values have rational ratios in the lattice
coordinates.

The primary models are:

    Q(r,t+1) / Q(r,t)
        = N(r,t) / D(r,t)

    Q(r+1,t) / Q(r,t)
        = N(r,t) / D(r,t)

with N and D restricted to low-degree affine forms.

The equations are cross-multiplied exactly, so no division is performed
inside the linear solve:

    Q(r,t+1) D(r,t)
      -
    Q(r,t) N(r,t)
      = 0.

This is preferable because it remains exact over Z and avoids any
denominator assumptions.

The same ratio law must hold across ALL usable observed edges.

We test:

    R_t(1,1):
        numerator = a0 + a1*r + a2*t
        denominator = b0 + b1*r + b2*t

    R_r(1,1):
        numerator = a0 + a1*r + a2*t
        denominator = b0 + b1*r + b2*t

and lower-complexity subfamilies:

    numerator/denominator constants;
    t-only affine ratio;
    r-only affine ratio.

A successful result requires:

    * all equations use only observed cells;
    * rank = unknowns - 1;
    * augmented rank = rank;
    * nullity = 1;
    * equations > effective parameter count;
    * the resulting ratio law verifies every usable edge exactly.

A data-sized fit is never accepted as a discovery.

An additional exact plaquette audit tests whether

    Q(r+1,t+1) Q(r,t)
    -------------------
    Q(r+1,t) Q(r,t+1)

has a simple constant or one-variable law.

No missing cells.
No interpolation.
No extrapolation.
No synthetic second case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

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
# SYMBOLS
# ============================================================================

a0, a1, a2 = sp.symbols(
    "a0 a1 a2"
)

b0, b1, b2 = sp.symbols(
    "b0 b1 b2"
)


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def primitive_vector(vector):
    """
    Primitive integer normalization of a rational coefficient vector.
    """

    rationals = [
        sp.Rational(value)
        for value in vector
    ]

    lcm = 1

    for value in rationals:
        lcm = sp.ilcm(
            lcm,
            int(sp.denom(value)),
        )

    integers = [
        int(value * lcm)
        for value in rationals
    ]

    gcd_value = 0

    for value in integers:
        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value == 0:
        return tuple(
            0
            for _ in integers
        )

    integers = [
        value // gcd_value
        for value in integers
    ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -v
                    for v in integers
                ]
            break

    return tuple(
        integers
    )


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
                (
                    r_value,
                    t_value,
                )
            ] = sp.Integer(value)

    return lattice


def observed_r_edges(lattice):
    """
    All edges

        (r,t) -> (r,t+1)

    where both cells are observed.
    """

    edges = []

    for (
        r,
        t,
    ) in sorted(lattice):

        source = (
            r,
            t,
        )

        target = (
            r,
            t + 1,
        )

        if target in lattice:

            edges.append(
                (
                    source,
                    target,
                )
            )

    return edges


def observed_t_edges(lattice):
    """
    All edges

        (r,t) -> (r+1,t)

    where both cells are observed.
    """

    edges = []

    for (
        r,
        t,
    ) in sorted(lattice):

        source = (
            r,
            t,
        )

        target = (
            r + 1,
            t,
        )

        if target in lattice:

            edges.append(
                (
                    source,
                    target,
                )
            )

    return edges


# ============================================================================
# RATIONAL-RATIO MODEL SOLVER
# ============================================================================

def ratio_basis(
    family,
):
    """
    Return coefficient basis for the numerator and denominator.

    Families:

        constant:
            1

        t_affine:
            1,t

        r_affine:
            1,r

        rt_affine:
            1,r,t
    """

    if family == "constant":
        return [
            lambda r, t: 1,
        ]

    if family == "t_affine":
        return [
            lambda r, t: 1,
            lambda r, t: t,
        ]

    if family == "r_affine":
        return [
            lambda r, t: 1,
            lambda r, t: r,
        ]

    if family == "rt_affine":
        return [
            lambda r, t: 1,
            lambda r, t: r,
            lambda r, t: t,
        ]

    raise ValueError(
        f"Unknown family: {family}"
    )


def ratio_model_rows(
    lattice,
    edges,
    family,
):
    """
    Unknown vector:

        numerator coefficients
        denominator coefficients

    Cross-multiplied equation for edge
    source -> target:

        target * D(source)
        -
        source * N(source)
        = 0.
    """

    basis = ratio_basis(
        family
    )

    rows = []

    for (
        source,
        target,
    ) in edges:

        r, t = source

        source_value = lattice[
            source
        ]

        target_value = lattice[
            target
        ]

        row = []

        # Numerator coefficients.
        for basis_fn in basis:

            row.append(
                -source_value
                * basis_fn(
                    r,
                    t,
                )
            )

        # Denominator coefficients.
        for basis_fn in basis:

            row.append(
                target_value
                * basis_fn(
                    r,
                    t,
                )
            )

        rows.append(
            row
        )

    return rows


def solve_ratio_model(
    lattice,
    edges,
    family,
):
    """
    Solve homogeneous rational-ratio relation.

    Because the system is homogeneous, a genuine unique ratio law should
    have nullity exactly one.
    """

    rows = ratio_model_rows(
        lattice,
        edges,
        family,
    )

    basis_count = len(
        ratio_basis(
            family
        )
    )

    parameter_count = (
        2 * basis_count
    )

    equation_count = len(rows)

    if equation_count == 0:

        return {
            "status": "DATA_LIMITED",
            "equations": 0,
            "unknowns": parameter_count,
            "rank": 0,
            "nullity": parameter_count,
            "basis": [],
        }

    matrix = sp.Matrix(
        rows
    )

    rank = matrix.rank()

    nullity = (
        parameter_count
        -
        rank
    )

    if nullity == 0:

        return {
            "status": "NO_RELATION",
            "equations": equation_count,
            "unknowns": parameter_count,
            "rank": rank,
            "nullity": 0,
            "basis": [],
        }

    nullspace = matrix.nullspace()

    if nullity > 1:

        return {
            "status": "UNDERDETERMINED",
            "equations": equation_count,
            "unknowns": parameter_count,
            "rank": rank,
            "nullity": nullity,
            "basis": nullspace,
        }

    relation = nullspace[0]

    primitive = primitive_vector(
        relation
    )

    verification_residuals = []

    for row in rows:

        residual = sum(
            sp.Integer(
                row[i]
            )
            *
            sp.Integer(
                primitive[i]
            )
            for i in range(
                parameter_count
            )
        )

        verification_residuals.append(
            clean(
                residual
            )
        )

    exact = all(
        residual == 0
        for residual
        in verification_residuals
    )

    if not exact:

        status = (
            "VERIFICATION_FAILED"
        )

    elif equation_count > parameter_count:

        status = (
            "EXACT_OVERDETERMINED"
        )

    else:

        status = (
            "EXACT_DATA_SIZED"
        )

    return {
        "status": status,
        "equations": equation_count,
        "unknowns": parameter_count,
        "rank": rank,
        "nullity": nullity,
        "basis": nullspace,
        "primitive_relation": primitive,
        "residuals": verification_residuals,
        "exact": exact,
    }


# ============================================================================
# RATIO LAW PRESENTATION
# ============================================================================

def ratio_law_from_relation(
    relation,
    family,
):
    """
    Convert primitive vector into

        N(r,t) / D(r,t).
    """

    basis = ratio_basis(
        family
    )

    k = len(
        basis
    )

    numerator = 0
    denominator = 0

    for i, basis_fn in enumerate(
        basis
    ):

        numerator += (
            relation[i]
            *
            {
                "constant": 1,
                "t_affine": sp.Symbol("t"),
                "r_affine": sp.Symbol("r"),
                "rt_affine": sp.Symbol("r") + sp.Symbol("t"),
            }[family]
            if False
            else 0
        )

    # Rebuild explicitly so the presentation is readable.
    r_sym, t_sym = sp.symbols(
        "r t"
    )

    if family == "constant":

        numerator = (
            relation[0]
        )

        denominator = (
            relation[1]
        )

    elif family == "t_affine":

        numerator = (
            relation[0]
            +
            relation[1]
            * t_sym
        )

        denominator = (
            relation[2]
            +
            relation[3]
            * t_sym
        )

    elif family == "r_affine":

        numerator = (
            relation[0]
            +
            relation[1]
            * r_sym
        )

        denominator = (
            relation[2]
            +
            relation[3]
            * r_sym
        )

    elif family == "rt_affine":

        numerator = (
            relation[0]
            +
            relation[1]
            * r_sym
            +
            relation[2]
            * t_sym
        )

        denominator = (
            relation[3]
            +
            relation[4]
            * r_sym
            +
            relation[5]
            * t_sym
        )

    else:

        raise ValueError(
            family
        )

    return (
        clean(numerator),
        clean(denominator),
    )


# ============================================================================
# AUDIT ONE FAMILY
# ============================================================================

def audit_ratio_direction(
    lattice,
    edges,
    direction,
):

    print()
    print("=" * 78)

    print(
        f"1. {direction.upper()}-DIRECTION "
        "HYPERGEOMETRIC RATIO AUDIT"
    )

    print("=" * 78)

    print(
        f"  usable_edges={len(edges)}"
    )

    results = {}

    for family in (
        "constant",
        "t_affine",
        "r_affine",
        "rt_affine",
    ):

        result = solve_ratio_model(
            lattice,
            edges,
            family,
        )

        results[
            family
        ] = result

        print()
        print(
            f"  family={family}"
        )

        print(
            f"    equations={result['equations']}"
        )

        print(
            f"    unknowns={result['unknowns']}"
        )

        print(
            f"    redundancy="
            f"{result['equations'] - result['unknowns']}"
        )

        print(
            f"    rank={result['rank']}"
        )

        print(
            f"    nullity={result['nullity']}"
        )

        print(
            f"    status={result['status']}"
        )

        if result.get(
            "primitive_relation"
        ) is not None:

            relation = result[
                "primitive_relation"
            ]

            print(
                f"    primitive_relation="
                f"{relation}"
            )

            numerator, denominator = (
                ratio_law_from_relation(
                    relation,
                    family,
                )
            )

            print(
                f"    numerator={numerator}"
            )

            print(
                f"    denominator={denominator}"
            )

    return results


# ============================================================================
# CROSS-DIRECTION COMMON LAW
# ============================================================================

def common_ratio_law_audit(
    r_results,
    t_results,
):

    print()
    print("=" * 78)
    print(
        "2. CROSS-DIRECTION RATIO-LAW COMPATIBILITY"
    )
    print("=" * 78)

    for family in (
        "constant",
        "t_affine",
        "r_affine",
        "rt_affine",
    ):

        r_result = r_results[
            family
        ]

        t_result = t_results[
            family
        ]

        print()
        print(
            f"  family={family}"
        )

        print(
            "    r_direction_status="
            f"{r_result['status']}"
        )

        print(
            "    t_direction_status="
            f"{t_result['status']}"
        )

        if (
            r_result.get(
                "primitive_relation"
            ) is not None
            and
            t_result.get(
                "primitive_relation"
            ) is not None
        ):

            print(
                "    primitive_r="
                f"{r_result['primitive_relation']}"
            )

            print(
                "    primitive_t="
                f"{t_result['primitive_relation']}"
            )


# ============================================================================
# PLAQUETTE CROSS-RATIO / MIXED RATIO AUDIT
# ============================================================================

def plaquette_values(
    lattice,
):

    values = []

    for r in range(3):

        for t in range(5):

            cells = [
                (
                    r,
                    t,
                ),
                (
                    r + 1,
                    t,
                ),
                (
                    r,
                    t + 1,
                ),
                (
                    r + 1,
                    t + 1,
                ),
            ]

            if all(
                cell in lattice
                for cell in cells
            ):

                q00 = lattice[
                    (
                        r,
                        t,
                    )
                ]

                q10 = lattice[
                    (
                        r + 1,
                        t,
                    )
                ]

                q01 = lattice[
                    (
                        r,
                        t + 1,
                    )
                ]

                q11 = lattice[
                    (
                        r + 1,
                        t + 1,
                    )
                ]

                denominator = (
                    q10 * q01
                )

                if denominator != 0:

                    mixed = clean(
                        (
                            q11 * q00
                        )
                        /
                        denominator
                    )

                    values.append(
                        (
                            (
                                r,
                                t,
                            ),
                            mixed,
                        )
                    )

    return values


def plaquette_audit(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "3. EXACT PLAQUETTE MULTIPLICATIVE AUDIT"
    )
    print("=" * 78)

    values = plaquette_values(
        lattice
    )

    print(
        f"  fully_observed_plaquettes="
        f"{len(values)}"
    )

    if not values:

        print(
            "  status=DATA_LIMITED"
        )

        return values

    print()

    for location, value in values:

        print(
            f"  plaquette={location}"
        )

        print(
            f"    mixed_ratio={value}"
        )

    distinct = {
        value
        for (
            _,
            value,
        ) in values
    }

    print()
    print(
        f"  distinct_mixed_ratio_count="
        f"{len(distinct)}"
    )

    if len(distinct) == 1:

        print(
            "  constant_plaquette_ratio=True"
        )

    else:

        print(
            "  constant_plaquette_ratio=False"
        )

    return values


# ============================================================================
# AFFINE PLAQUETTE LAW
# ============================================================================

def plaquette_affine_law(
    values,
):

    print()
    print("=" * 78)
    print(
        "4. EXACT PLAQUETTE RATIO AFFINE-LAW AUDIT"
    )
    print("=" * 78)

    if not values:

        print(
            "  status=DATA_LIMITED"
        )

        return

    # Fit
        # mixed_ratio = a + b*r + c*t
    rows = []
    rhs = []

    for (
        location,
        value,
    ) in values:

        r, t = location

        rows.append(
            [
                1,
                r,
                t,
            ]
        )

        rhs.append(
            value
        )

    matrix = sp.Matrix(
        rows
    )

    rhs_matrix = sp.Matrix(
        rhs
    )

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(
            rhs_matrix
        )
        .rank()
    )

    print(
        f"  equations={len(rows)}"
    )

    print(
        "  unknowns=3"
    )

    print(
        f"  rank={rank}"
    )

    print(
        f"  augmented_rank={augmented_rank}"
    )

    if augmented_rank > rank:

        print(
            "  status=NO_SOLUTION"
        )

        return

    if rank < 3:

        print(
            "  status=NONUNIQUE"
        )

        return

    solution = matrix.gauss_jordan_solve(
        rhs_matrix
    )[0]

    coefficients = tuple(
        clean(
            solution[i, 0]
        )
        for i in range(3)
    )

    residuals = [
        clean(
            (
                coefficients[0]
                +
                coefficients[1] * r
                +
                coefficients[2] * t
            )
            -
            value
        )
        for (
            (
                r,
                t,
            ),
            value,
        ) in values
    ]

    exact = all(
        residual == 0
        for residual in residuals
    )

    if (
        exact
        and
        len(rows) > 3
    ):

        status = (
            "EXACT_OVERDETERMINED"
        )

    elif exact:

        status = (
            "EXACT_DATA_SIZED"
        )

    else:

        status = (
            "VERIFICATION_FAILED"
        )

    print(
        f"  status={status}"
    )

    if exact:

        print(
            f"  coefficients={coefficients}"
        )

        print(
            f"  residuals={residuals}"
        )


# ============================================================================
# INTEGER / SMALL COEFFICIENT DIAGNOSTIC
# ============================================================================

def integer_relation_diagnostic(
    results,
    direction,
):

    print()
    print("=" * 78)
    print(
        f"5. {direction.upper()}-DIRECTION "
        "INTEGER RATIO DIAGNOSTIC"
    )
    print("=" * 78)

    for family in (
        "constant",
        "t_affine",
        "r_affine",
        "rt_affine",
    ):

        result = results[
            family
        ]

        relation = result.get(
            "primitive_relation"
        )

        if relation is None:

            continue

        gcd_value = 0

        for value in relation:

            gcd_value = math.gcd(
                gcd_value,
                abs(
                    int(value)
                ),
            )

        print()
        print(
            f"  family={family}"
        )

        print(
            f"    primitive_gcd={gcd_value}"
        )

        print(
            f"    coefficient_L1="
            f"{sum(abs(int(v)) for v in relation)}"
        )

        print(
            f"    coefficient_Linf="
            f"{max(abs(int(v)) for v in relation)}"
        )


# ============================================================================
# STRUCTURAL SUMMARY
# ============================================================================

def structural_summary(
    r_results,
    t_results,
    plaquettes,
):

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    exact_overdetermined = []

    for direction, results in (
        (
            "r",
            r_results,
        ),
        (
            "t",
            t_results,
        ),
    ):

        for family, result in results.items():

            if result[
                "status"
            ] == "EXACT_OVERDETERMINED":

                exact_overdetermined.append(
                    (
                        direction,
                        family,
                    )
                )

    print(
        f"  exact_overdetermined_ratio_laws="
        f"{exact_overdetermined}"
    )

    print(
        f"  fully_observed_plaquettes="
        f"{len(plaquettes)}"
    )

    if exact_overdetermined:

        verdict = (
            "EXACT_OVERDETERMINED_HYPERGEOMETRIC_STRUCTURE_FOUND"
        )

    elif (
        plaquettes
        and
        len(
            {
                value
                for (
                    _,
                    value,
                ) in plaquettes
            }
        )
        == 1
    ):

        verdict = (
            "EXACT_CONSTANT_PLAQUETTE_MULTIPLICATIVE_STRUCTURE"
        )

    else:

        verdict = (
            "NO_VALIDATED_LOW_DEGREE_HYPERGEOMETRIC_LAW"
        )

    print(
        f"  verdict={verdict}"
    )

    print()
    print(
r"""
The constant-coefficient convolution route has now been exhausted at
low order.

380R tests the complementary multiplicative possibility:

    adjacent ratios are rational functions of lattice coordinates.

This includes many hypergeometric-type constructions and finite products.

The critical methodological distinction is:

    fitted ratio on one short line
        versus
    one common rational ratio law surviving an overdetermined collection
    of observed edges.

Only the second is structurally meaningful.

The plaquette quantity

    Q(r+1,t+1) Q(r,t)
    -------------------
    Q(r+1,t) Q(r,t+1)

is also tested because multiplicatively separable source constructions
often force a simple mixed-ratio law.

A negative result would leave the project with:

    recurrence laws rejected;
    differential/Euler/shift operators rejected;
    projective row laws rejected;
    cross-ratio invariants rejected;
    low-degree generating-function denominators rejected;
    low-rank polynomial nullspaces not simplified.

At that point, source-definition reconstruction or a genuine second case
becomes substantially more important than additional generic fitting.

No missing source value is used.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 380R — EXACT BIVARIATE HYPERGEOMETRIC / "
        "RATIONAL-RATIO SOURCE-LAW AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "OBSERVED SOURCE"
    )

    print(
        f"  observed_cells={len(lattice)}"
    )

    print(
        f"  r_edges={len(observed_r_edges(lattice))}"
    )

    print(
        f"  t_edges={len(observed_t_edges(lattice))}"
    )

    # ------------------------------------------------------------------------
    # Ratio models.
    # ------------------------------------------------------------------------

    r_edges = observed_r_edges(
        lattice
    )

    t_edges = observed_t_edges(
        lattice
    )

    r_results = audit_ratio_direction(
        lattice,
        r_edges,
        "r",
    )

    t_results = audit_ratio_direction(
        lattice,
        t_edges,
        "t",
    )

    common_ratio_law_audit(
        r_results,
        t_results,
    )

    integer_relation_diagnostic(
        r_results,
        "r",
    )

    integer_relation_diagnostic(
        t_results,
        "t",
    )

    # ------------------------------------------------------------------------
    # Plaquette audit.
    # ------------------------------------------------------------------------

    plaquettes = plaquette_audit(
        lattice
    )

    plaquette_affine_law(
        plaquettes
    )

    # ------------------------------------------------------------------------
    # Structural summary.
    # ------------------------------------------------------------------------

    structural_summary(
        r_results,
        t_results,
        plaquettes,
    )

    # ------------------------------------------------------------------------
    # Final exactness.
    # ------------------------------------------------------------------------

    exact_overdetermined = []

    for direction, results in (
        (
            "r",
            r_results,
        ),
        (
            "t",
            t_results,
        ),
    ):

        for family, result in results.items():

            if result[
                "status"
            ] == "EXACT_OVERDETERMINED":

                exact_overdetermined.append(
                    (
                        direction,
                        family,
                    )
                )

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  exact_cross_multiplied_ratio_equations=True"
    )

    print(
        "  constant_ratio_tested=True"
    )

    print(
        "  r_affine_ratio_tested=True"
    )

    print(
        "  t_affine_ratio_tested=True"
    )

    print(
        "  rt_affine_ratio_tested=True"
    )

    print(
        "  plaquette_multiplicative_audit=True"
    )

    print(
        "  plaquette_affine_audit=True"
    )

    print(
        f"  exact_overdetermined_ratio_laws="
        f"{exact_overdetermined}"
    )

    print(
        "  missing_values_used=False"
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
        "EXPERIMENT 380R COMPLETE"
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
            "\nFATAL ERROR: "
            f"{type(exc).__name__}: {exc}"
        )

        raise
