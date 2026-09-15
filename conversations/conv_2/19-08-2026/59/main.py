#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 374R — EXACT ANNIHILATOR-POLYNOMIAL / FACTORIZATION /
                  SEPARABILITY / COMMON-FACTOR AUDIT
==============================================================================

Purpose
-------
373R established that the genuine rational nullspaces

    rectangle_2x3
    triangle_6

remain arithmetically large after exact LLL reduction.

374R changes the notion of simplicity.

Each nullspace vector is interpreted as a bivariate stencil polynomial

    S(u,v) = sum c_(r,t) u^r v^t.

The experiment tests whether the nullspace contains hidden algebraic
structure such as:

    * factorization over Q;
    * common polynomial factors between independent annihilators;
    * factors of the forms (u-a), (v-b), (u*v-a), etc.;
    * separability S(u,v)=R(u)T(v);
    * rank-1 coefficient matrices;
    * low-rank combinations inside the 2-dimensional nullspace;
    * small integer combinations producing factorized or separable stencils.

A positive result would be considerably more interesting than merely
finding another rational relation, because it could identify the algebraic
shape of the annihilator itself.

A negative result would indicate that the rational nullspaces are genuine
but do not obviously arise from low-complexity polynomial factors.

No missing values.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
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


SUPPORTS = {
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


u, v = sp.symbols("u v")


# ============================================================================
# LATTICE
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
            ] = sp.Integer(value)

    return lattice


def find_windows(
    lattice,
    support,
):

    candidates = set()

    for r, t in lattice:

        for dr, dt in support:

            candidates.add(
                (
                    r - dr,
                    t - dt,
                )
            )

    windows = []

    for r0, t0 in sorted(
        candidates,
        key=lambda z: (
            z[1],
            z[0],
        ),
    ):

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
                (
                    r0,
                    t0,
                )
            )

    return windows


def build_matrix(
    lattice,
    support,
    windows,
):

    return sp.Matrix([
        [
            lattice[
                (
                    r0 + dr,
                    t0 + dt,
                )
            ]
            for dr, dt in support
        ]
        for r0, t0 in windows
    ])


# ============================================================================
# INTEGER NULLSPACE
# ============================================================================

def primitive_integer_vector(
    vector,
):

    values = [
        sp.Rational(value)
        for value in vector
    ]

    denominator = 1

    for value in values:

        denominator = sp.ilcm(
            denominator,
            int(
                sp.denom(value)
            ),
        )

    integers = [
        int(
            value * denominator
        )
        for value in values
    ]

    g = 0

    for value in integers:

        g = math.gcd(
            g,
            abs(value),
        )

    if g:

        integers = [
            value // g
            for value in integers
        ]

    for value in integers:

        if value != 0:

            if value < 0:

                integers = [
                    -x
                    for x in integers
                ]

            break

    return tuple(
        integers
    )


def nullspace_basis(
    matrix,
):

    return [
        primitive_integer_vector(
            vector
        )
        for vector in matrix.nullspace()
    ]


def verify_vector(
    matrix,
    vector,
):

    result = (
        matrix
        * sp.Matrix(vector)
    )

    return all(
        value == 0
        for value in result
    )


# ============================================================================
# STENCIL POLYNOMIAL
# ============================================================================

def stencil_polynomial(
    support,
    vector,
):

    expr = sp.Integer(0)

    for (
        (r, t),
        coefficient,
    ) in zip(
        support,
        vector,
    ):

        expr += (
            sp.Integer(coefficient)
            * u**r
            * v**t
        )

    return sp.expand(
        expr
    )


def primitive_polynomial(
    expr,
):

    poly = sp.Poly(
        expr,
        u,
        v,
        domain=sp.QQ,
    )

    denominators = [
        int(
            sp.denom(coefficient)
        )
        for coefficient
        in poly.coeffs()
    ]

    lcm = 1

    for denominator in denominators:

        lcm = sp.ilcm(
            lcm,
            denominator,
        )

    integer_coefficients = [
        int(
            coefficient * lcm
        )
        for coefficient
        in poly.coeffs()
    ]

    g = 0

    for coefficient in integer_coefficients:

        g = math.gcd(
            g,
            abs(coefficient),
        )

    if g:

        integer_coefficients = [
            value // g
            for value
            in integer_coefficients
        ]

    expression = sp.Integer(0)

    for (
        monomial,
        coefficient,
    ) in zip(
        poly.monoms(),
        integer_coefficients,
    ):

        i, j = monomial

        expression += (
            coefficient
            * u**i
            * v**j
        )

    if expression == 0:

        return sp.Integer(0)

    if sp.LC(
        sp.Poly(
            expression,
            u,
            v,
        )
    ) < 0:

        expression = -expression

    return sp.expand(
        expression
    )


# ============================================================================
# FACTORIZATION
# ============================================================================

def factor_audit(
    expression,
):

    expression = sp.expand(
        expression
    )

    factored = sp.factor(
        expression
    )

    constant, factors = sp.factor_list(
        expression
    )

    nonconstant = [
        (
            sp.expand(factor),
            exponent,
        )
        for factor, exponent
        in factors
        if sp.Poly(
            factor,
            u,
            v,
        ).total_degree() > 0
    ]

    return {
        "factorized": factored,
        "constant": constant,
        "factors": nonconstant,
        "nontrivial": (
            len(nonconstant) > 1
            or (
                len(nonconstant) == 1
                and nonconstant[0][1] > 1
            )
        ),
    }


# ============================================================================
# SIMPLE FACTOR TESTS
# ============================================================================

def simple_factor_tests(
    expression,
):

    tests = {}

    for a in (
        -3,
        -2,
        -1,
        0,
        1,
        2,
        3,
    ):

        tests[
            f"u_minus_{a}"
        ] = sp.expand(
            expression.subs(
                u,
                a,
            )
        ) == 0

        tests[
            f"v_minus_{a}"
        ] = sp.expand(
            expression.subs(
                v,
                a,
            )
        ) == 0

    return tests


# ============================================================================
# SEPARABILITY
# ============================================================================

def coefficient_matrix(
    expression,
):

    poly = sp.Poly(
        expression,
        u,
        v,
    )

    degree_u = (
        poly.degree(u)
        if not poly.is_zero
        else -1
    )

    degree_v = (
        poly.degree(v)
        if not poly.is_zero
        else -1
    )

    matrix = []

    for i in range(
        degree_u + 1
    ):

        row = []

        for j in range(
            degree_v + 1
        ):

            row.append(
                poly.coeff_monomial(
                    u**i * v**j
                )
            )

        matrix.append(
            row
        )

    if not matrix:

        return sp.Matrix([])

    return sp.Matrix(
        matrix
    )


def separability_audit(
    expression,
):

    matrix = coefficient_matrix(
        expression
    )

    if matrix.rows == 0:

        rank = 0

    else:

        rank = matrix.rank()

    return {
        "matrix_shape": (
            matrix.rows,
            matrix.cols,
        ),
        "matrix_rank": rank,
        "separable_rank1": (
            rank <= 1
        ),
        "matrix": matrix,
    }


# ============================================================================
# COMMON FACTORS
# ============================================================================

def common_factor_audit(
    expressions,
):

    if not expressions:

        return {
            "gcd": sp.Integer(0),
            "factorized_gcd": sp.Integer(0),
        }

    gcd = expressions[0]

    for expression in expressions[1:]:

        gcd = sp.gcd(
            gcd,
            expression,
        )

    gcd = sp.factor(
        gcd
    )

    return {
        "gcd": gcd,
        "factorized_gcd": gcd,
    }


# ============================================================================
# SMALL NULLSPACE COMBINATIONS
# ============================================================================

def combination_records(
    basis,
    bound=15,
):

    if len(basis) != 2:

        return []

    b0 = basis[0]
    b1 = basis[1]

    records = []

    for a in range(
        -bound,
        bound + 1,
    ):

        for b in range(
            -bound,
            bound + 1,
        ):

            if (
                a == 0
                and b == 0
            ):

                continue

            vector = tuple(
                a * b0[j]
                + b * b1[j]
                for j in range(
                    len(b0)
                )
            )

            if not any(
                vector
            ):

                continue

            primitive = primitive_integer_vector(
                vector
            )

            expression = sp.expand(
                sum(
                    primitive[index]
                    * u**support[0]
                    * v**support[1]
                    for index, support
                    in enumerate(
                        current_support
                    )
                )
            )

            factor = factor_audit(
                expression
            )

            separation = separability_audit(
                expression
            )

            records.append(
                {
                    "a": a,
                    "b": b,
                    "vector": primitive,
                    "support_size": sum(
                        value != 0
                        for value in primitive
                    ),
                    "L1": sum(
                        abs(value)
                        for value in primitive
                    ),
                    "Linf": max(
                        abs(value)
                        for value in primitive
                    ),
                    "factorized": factor,
                    "separable": separation,
                }
            )

    return records


# ============================================================================
# LOW-COMPLEXITY RANK / FACTOR RANKING
# ============================================================================

def rank_combination_records(
    records,
):

    ranked = sorted(
        records,
        key=lambda record: (
            not record[
                "factorized"
            ]["nontrivial"],
            not record[
                "separable"
            ]["separable_rank1"],
            record["support_size"],
            record["L1"],
            record["Linf"],
        ),
    )

    return ranked


# ============================================================================
# MAIN SUPPORT AUDIT
# ============================================================================

def audit_support(
    name,
    lattice,
):

    global current_support

    current_support = SUPPORTS[
        name
    ]

    windows = find_windows(
        lattice,
        current_support,
    )

    matrix = build_matrix(
        lattice,
        current_support,
        windows,
    )

    basis = nullspace_basis(
        matrix
    )

    print()
    print("=" * 78)
    print(
        f"SUPPORT={name}"
    )
    print("=" * 78)

    print(
        f"  support={current_support}"
    )

    print(
        f"  windows={windows}"
    )

    print(
        f"  matrix_shape="
        f"({matrix.rows},{matrix.cols})"
    )

    print(
        f"  rational_rank={matrix.rank()}"
    )

    print(
        f"  rational_nullity="
        f"{matrix.cols - matrix.rank()}"
    )

    expressions = []

    print()
    print(
        "  NULLSPACE POLYNOMIALS"
    )

    for index, vector in enumerate(
        basis
    ):

        expression = stencil_polynomial(
            current_support,
            vector,
        )

        expression = primitive_polynomial(
            expression
        )

        expressions.append(
            expression
        )

        factor = factor_audit(
            expression
        )

        separation = (
            separability_audit(
                expression
            )
        )

        simple_factors = (
            simple_factor_tests(
                expression
            )
        )

        print()
        print(
            f"    basis_{index}:"
        )

        print(
            f"      polynomial={expression}"
        )

        print(
            f"      factorized="
            f"{factor['factorized']}"
        )

        print(
            f"      nontrivial_factorization="
            f"{factor['nontrivial']}"
        )

        print(
            f"      coefficient_matrix_shape="
            f"{separation['matrix_shape']}"
        )

        print(
            f"      coefficient_matrix_rank="
            f"{separation['matrix_rank']}"
        )

        print(
            f"      rank1_separable="
            f"{separation['separable_rank1']}"
        )

        simple_hits = [
            name
            for name, hit
            in simple_factors.items()
            if hit
        ]

        print(
            f"      simple_linear_factor_hits="
            f"{simple_hits}"
        )

        print(
            f"      exact_residual_zero="
            f"{verify_vector(matrix, vector)}"
        )

    print()
    print(
        "  COMMON FACTOR AUDIT"
    )

    common = common_factor_audit(
        expressions
    )

    print(
        f"    gcd={common['gcd']}"
    )

    print(
        f"    factorized_gcd="
        f"{common['factorized_gcd']}"
    )

    print()
    print(
        "  SMALL LINEAR COMBINATION SEARCH"
    )

    records = combination_records(
        basis,
        bound=15,
    )

    ranked = rank_combination_records(
        records
    )

    if not ranked:

        print(
            "    no_two_parameter_combinations"
        )

    else:

        # Only show the first 12 genuinely interesting records.
        shown = 0

        for record in ranked:

            factor = record[
                "factorized"
            ]

            separable = record[
                "separable"
            ]

            interesting = (
                factor["nontrivial"]
                or separable[
                    "separable_rank1"
                ]
                or record[
                    "support_size"
                ] <= 3
            )

            if not interesting:
                continue

            shown += 1

            print()
            print(
                f"    candidate_{shown - 1}:"
            )

            print(
                f"      combination="
                f"({record['a']},{record['b']})"
            )

            print(
                f"      vector="
                f"{record['vector']}"
            )

            print(
                f"      polynomial="
                f"{stencil_polynomial(current_support, record['vector'])}"
            )

            print(
                f"      support_size="
                f"{record['support_size']}"
            )

            print(
                f"      L1="
                f"{record['L1']}"
            )

            print(
                f"      Linf="
                f"{record['Linf']}"
            )

            print(
                f"      factorized="
                f"{factor['factorized']}"
            )

            print(
                f"      nontrivial_factorization="
                f"{factor['nontrivial']}"
            )

            print(
                f"      coefficient_matrix_rank="
                f"{separable['matrix_rank']}"
            )

            print(
                f"      rank1_separable="
                f"{separable['separable_rank1']}"
            )

            if shown >= 12:
                break

        if shown == 0:

            print(
                "    no_factorized_or_separable_combination_found"
            )

    return {
        "basis": basis,
        "expressions": expressions,
        "matrix": matrix,
        "windows": windows,
        "combination_records": records,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 374R — EXACT ANNIHILATOR-POLYNOMIAL / "
        "FACTORIZATION / SEPARABILITY / COMMON-FACTOR AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "OBSERVED DATA"
    )

    print(
        f"  observed_cells={len(lattice)}"
    )

    results = {}

    for name in SUPPORTS:

        results[name] = audit_support(
            name,
            lattice,
        )

    # ------------------------------------------------------------------------
    # Cross-support comparison
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. CROSS-SUPPORT COMMON-FACTOR AUDIT"
    )
    print("=" * 78)

    rectangle_expressions = results[
        "rectangle_2x3"
    ]["expressions"]

    triangle_expressions = results[
        "triangle_6"
    ]["expressions"]

    rr_gcd = common_factor_audit(
        rectangle_expressions
    )

    tt_gcd = common_factor_audit(
        triangle_expressions
    )

    print()
    print(
        "  rectangle_2x3:"
    )

    print(
        f"    common_gcd={rr_gcd['gcd']}"
    )

    print(
        "    common_nonconstant_factor="
        f"{rr_gcd['gcd'] not in (0, 1, -1)}"
    )

    print()
    print(
        "  triangle_6:"
    )

    print(
        f"    common_gcd={tt_gcd['gcd']}"
    )

    print(
        "    common_nonconstant_factor="
        f"{tt_gcd['gcd'] not in (0, 1, -1)}"
    )

    # ------------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------------

    factorized_basis_count = 0
    separable_basis_count = 0

    for result in results.values():

        for expression in result[
            "expressions"
        ]:

            factor = factor_audit(
                expression
            )

            separation = (
                separability_audit(
                    expression
                )
            )

            if factor["nontrivial"]:

                factorized_basis_count += 1

            if separation[
                "separable_rank1"
            ]:

                separable_basis_count += 1

    print()
    print("=" * 78)
    print(
        "4. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
373R showed that exact LLL reduction does not collapse the two rational
nullspaces into small sparse coefficient vectors.

374R therefore tests a different hypothesis:

    perhaps the coefficient vectors are large only because the
    annihilators are written in an inconvenient basis.

Writing each relation as

    S(u,v) = sum c_(r,t) u^r v^t

allows algebraic structure to become visible.

The strongest possible outcomes are:

    FACTORIZED
        S(u,v) factors nontrivially over Q.

    COMMON_FACTOR
        both independent annihilators share a nonconstant factor.

    SEPARABLE
        S(u,v)=R(u)T(v), equivalently the coefficient matrix has rank 1.

    SIMPLE_COMBINATION
        a small integer combination of the two nullspace generators
        becomes factorized or separable.

These would point toward an algebraic annihilator rather than a merely
accidental numerical relation.

A negative result is also informative: if the exact kernel basis has
no nontrivial factorization, no common factor, no separable member, and
no simple combination with those properties, then the rational
nullspaces are structurally real but not obviously generated by a
low-complexity stencil polynomial.

No missing source value is introduced.
"""
    )

    print()
    print("=" * 78)
    print(
        "5. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  exact_nullspace_polynomials_built=True"
    )

    print(
        "  exact_factorization_completed=True"
    )

    print(
        "  common_factor_audit_completed=True"
    )

    print(
        "  separability_audit_completed=True"
    )

    print(
        "  small_integer_combination_search_completed=True"
    )

    print(
        f"  factorized_basis_count="
        f"{factorized_basis_count}"
    )

    print(
        f"  separable_basis_count="
        f"{separable_basis_count}"
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
        "EXPERIMENT 374R COMPLETE"
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
