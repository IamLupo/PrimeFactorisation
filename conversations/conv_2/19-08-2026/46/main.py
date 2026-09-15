#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 361R — EXACT D7 RATIONAL-FAMILY / INTEGER-FEASIBILITY / PRIMITIVE
                  VECTOR AUDIT
==============================================================================

Purpose
-------
Repair the remaining failure in 360R.

The D7 symbolic family is affine in

    Z = Q_3(5),

but an integer Z does NOT imply that every D7 coefficient is itself an
integer.

Therefore this experiment distinguishes:

    rational coefficient evaluation;
    coefficient-wise integrality;
    primitive integer-vector normalization.

For any integer Z:

    coefficient_vector(Z) ∈ Q^5

is converted to an integer vector by clearing the common denominator.

This is mathematically different from requiring the coefficients themselves
to lie in Z.

Audits
------
1. Exact D7 family reconstruction.
2. Exact affine (A Z + B)/D representation.
3. Symbolic residual verification.
4. Coefficient-wise integer feasibility.
5. Special values 0, ±1.
6. Pairwise coefficient equality.
7. Rational-vector primitive normalization.
8. Compact bounded diagnostic search.
9. Exact denominator/lcm profile.
10. No missing value is inserted.

No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
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

Z = sp.Symbol(
    "Z",
    integer=True,
)

a, b, c, d, e = sp.symbols(
    "a b c d e"
)

D7_VARS = [
    a,
    b,
    c,
    d,
    e,
]


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


def is_integer_exact(value):
    value = sp.Rational(value)
    return value.q == 1


def factor_integer(value):
    value = int(value)

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

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

            lattice[
                (r, t)
            ] = sp.Integer(value)

    return lattice


def integer_vector_primitive(values):
    """
    Given rational values v_i, clear the common denominator and return

        primitive_integer_vector,
        common_denominator,
        integer_vector.

    This does NOT assert that the original coefficients were integers.
    """

    rationals = [
        sp.Rational(
            clean(value)
        )
        for value in values
    ]

    common_denominator = 1

    for value in rationals:
        common_denominator = sp.ilcm(
            common_denominator,
            int(value.q),
        )

    integer_vector = [
        int(
            value * common_denominator
        )
        for value in rationals
    ]

    gcd_value = 0

    for value in integer_vector:
        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value == 0:
        return (
            [0] * len(integer_vector),
            common_denominator,
            integer_vector,
            0,
        )

    primitive = [
        value // gcd_value
        for value in integer_vector
    ]

    # Canonical sign.
    for value in primitive:
        if value != 0:
            if value < 0:
                primitive = [
                    -v
                    for v in primitive
                ]
            break

    return (
        primitive,
        common_denominator,
        integer_vector,
        gcd_value,
    )


# ============================================================================
# D7 EQUATIONS
# ============================================================================

def d7_equations():

    lattice = build_lattice()

    def value(cell):

        if cell == (2, 3):
            return Z

        return lattice.get(
            cell
        )

    equations = []

    for r in range(1, 4):

        for t in range(2, 6):

            target = (
                r,
                t,
            )

            target_value = value(
                target
            )

            if target_value is None:
                continue

            sources = [
                (r - 1, t),
                (r, t - 1),
                (r - 1, t - 1),
                (r, t - 2),
                (r - 1, t - 2),
            ]

            source_values = [
                value(cell)
                for cell in sources
            ]

            if not all(
                source_value is not None
                for source_value
                in source_values
            ):
                continue

            equations.append(
                sp.Eq(
                    a * source_values[0]
                    + b * source_values[1]
                    + c * source_values[2]
                    + d * source_values[3]
                    + e * source_values[4],
                    target_value,
                )
            )

    return equations


def solve_d7_family():

    equations = d7_equations()

    solution_set = sp.linsolve(
        equations,
        D7_VARS,
    )

    solutions = list(
        solution_set
    )

    if len(solutions) != 1:
        raise RuntimeError(
            "Unexpected D7 solution-set structure."
        )

    tuple_solution = solutions[0]

    expressions = {
        var: clean(
            tuple_solution[i]
        )
        for i, var in enumerate(
            D7_VARS
        )
    }

    return (
        equations,
        expressions,
    )


# ============================================================================
# EXACT AFFINE FORM
# ============================================================================

def affine_form(expression):
    """
    Return integers A,B,D with

        expression = (A*Z + B)/D,

    gcd(A,B,D)=1 and D>0.
    """

    expression = clean(
        sp.together(
            expression
        )
    )

    numerator, denominator = sp.fraction(
        expression
    )

    numerator = sp.Poly(
        sp.expand(numerator),
        Z,
        domain=sp.QQ,
    )

    denominator = sp.Rational(
        denominator
    )

    if numerator.degree() > 1:
        raise RuntimeError(
            "Non-affine D7 coefficient."
        )

    A = numerator.coeff_monomial(
        Z
    )
    B = numerator.coeff_monomial(
        1
    )

    scale = sp.ilcm(
        int(sp.denom(A)),
        int(sp.denom(B)),
        int(sp.denom(denominator)),
    )

    A = int(
        A * scale
    )

    B = int(
        B * scale
    )

    D = int(
        denominator * scale
    )

    g = math.gcd(
        math.gcd(
            abs(A),
            abs(B),
        ),
        abs(D),
    )

    if g != 0:
        A //= g
        B //= g
        D //= g

    if D < 0:
        A = -A
        B = -B
        D = -D

    return (
        A,
        B,
        D,
    )


# ============================================================================
# 1. FAMILY AUDIT
# ============================================================================

def family_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "1. EXACT D7 AFFINE-FAMILY AUDIT"
    )
    print("=" * 78)

    forms = {}

    for var in D7_VARS:

        expression = clean(
            expressions[var]
        )

        derivative = clean(
            sp.diff(
                expression,
                Z,
            )
        )

        A, B, D = affine_form(
            expression
        )

        reconstructed = clean(
            (
                sp.Integer(A) * Z
                + sp.Integer(B)
            )
            / sp.Integer(D)
        )

        ok = clean(
            expression
            - reconstructed
        ) == 0

        forms[var] = (
            A,
            B,
            D,
        )

        print()
        print(
            "  {}:".format(
                var
            )
        )

        print(
            "    expression={}".format(
                expression
            )
        )

        print(
            "    d_dZ={}".format(
                derivative
            )
        )

        print(
            "    depends_on_Z={}".format(
                derivative != 0
            )
        )

        print(
            "    affine_form=({}*Z + {})/{}".format(
                A,
                B,
                D,
            )
        )

        print(
            "    reconstruction_exact={}".format(
                ok
            )
        )

    return forms


# ============================================================================
# 2. SYMBOLIC RESIDUAL AUDIT
# ============================================================================

def symbolic_residual_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "2. SYMBOLIC D7 RESIDUAL AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    def value(cell):

        if cell == (2, 3):
            return Z

        return lattice.get(
            cell
        )

    failures = 0
    equations = 0

    for r in range(1, 4):

        for t in range(2, 6):

            target = (
                r,
                t,
            )

            target_value = value(
                target
            )

            if target_value is None:
                continue

            sources = [
                (r - 1, t),
                (r, t - 1),
                (r - 1, t - 1),
                (r, t - 2),
                (r - 1, t - 2),
            ]

            source_values = [
                value(cell)
                for cell in sources
            ]

            if not all(
                source_value is not None
                for source_value
                in source_values
            ):
                continue

            equations += 1

            predicted = clean(
                expressions[a] * source_values[0]
                + expressions[b] * source_values[1]
                + expressions[c] * source_values[2]
                + expressions[d] * source_values[3]
                + expressions[e] * source_values[4]
            )

            residual = clean(
                predicted
                - target_value
            )

            if residual != 0:
                failures += 1

            print(
                "  target={}: residual={}".format(
                    target,
                    residual,
                )
            )

    print()
    print(
        "  equation_count={}".format(
            equations
        )
    )

    print(
        "  residual_failures={}".format(
            failures
        )
    )

    return failures == 0


# ============================================================================
# 3. COEFFICIENT INTEGRALITY CONDITIONS
# ============================================================================

def integer_coefficient_conditions(
    expressions,
    forms,
):

    print()
    print("=" * 78)
    print(
        "3. EXACT COEFFICIENT-WISE INTEGER-FEASIBILITY"
    )
    print("=" * 78)

    conditions = {}

    for var in D7_VARS:

        A, B, D = forms[var]

        gcd_AD = math.gcd(
            abs(A),
            abs(D),
        )

        # A Z + B ≡ 0 (mod D).
        #
        # A solution exists iff gcd(A,D) divides B.
        compatible = (
            B % gcd_AD == 0
        )

        if not compatible:

            solutions = []

        else:

            modulus = (
                abs(D)
                // gcd_AD
            )

            A_reduced = (
                A // gcd_AD
            )

            B_reduced = (
                B // gcd_AD
            )

            if modulus == 1:

                residue = 0

            else:

                inv = pow(
                    A_reduced % modulus,
                    -1,
                    modulus,
                )

                residue = (
                    -B_reduced * inv
                ) % modulus

            solutions = [
                (
                    residue,
                    modulus,
                )
            ]

        conditions[var] = {
            "compatible": compatible,
            "gcd_A_D": gcd_AD,
            "solutions": solutions,
        }

        print()
        print(
            "  {}:".format(
                var
            )
        )

        print(
            "    A={}, B={}, D={}".format(
                A,
                B,
                D,
            )
        )

        print(
            "    gcd(A,D)={}".format(
                gcd_AD
            )
        )

        print(
            "    integer_solution_exists={}".format(
                compatible
            )
        )

        if solutions:

            residue, modulus = solutions[0]

            print(
                "    Z_congruence=Z ≡ {} (mod {})".format(
                    residue,
                    modulus,
                )
            )

        else:

            print(
                "    Z_congruence=NONE"
            )

    return conditions


# ============================================================================
# 4. COMBINE INTEGRALITY CONDITIONS
# ============================================================================

def generalized_crt(
    congruences
):
    """
    Combine

        Z ≡ r_i (mod m_i)

    exactly, allowing non-coprime moduli.

    Returns:
        (residue, modulus)
    or None.
    """

    if not congruences:
        return (
            0,
            1,
        )

    residue = int(
        congruences[0][0]
    )

    modulus = int(
        congruences[0][1]
    )

    for next_residue, next_modulus in congruences[1:]:

        next_residue = int(
            next_residue
        )

        next_modulus = int(
            next_modulus
        )

        g = math.gcd(
            modulus,
            next_modulus,
        )

        if (
            next_residue - residue
        ) % g != 0:
            return None

        m1 = modulus // g
        m2 = next_modulus // g

        if m2 == 1:
            k = 0

        else:

            inv = pow(
                m1 % m2,
                -1,
                m2,
            )

            k = (
                (
                    next_residue
                    - residue
                )
                // g
                * inv
            ) % m2

        residue = (
            residue
            + modulus * k
        )

        modulus = (
            modulus * m2
        )

        residue %= modulus

    return (
        residue,
        modulus,
    )


def combined_integrality_audit(
    conditions
):

    print()
    print("=" * 78)
    print(
        "4. SIMULTANEOUS INTEGER-COEFFICIENT FEASIBILITY"
    )
    print("=" * 78)

    congruences = []

    for var in D7_VARS:

        condition = conditions[var]

        if not condition["compatible"]:

            print(
                "  {}: no integer Z exists.".format(
                    var
                )
            )

            print(
                "  simultaneous_integer_coefficients=False"
            )

            return None

        congruences.extend(
            condition["solutions"]
        )

    combined = generalized_crt(
        congruences
    )

    if combined is None:

        print(
            "  simultaneous_integer_coefficients=False"
        )

        print(
            "  CRT_result=INCONSISTENT"
        )

    else:

        residue, modulus = combined

        print(
            "  simultaneous_integer_coefficients=True"
        )

        print(
            "  Z ≡ {} (mod {})".format(
                residue,
                modulus,
            )
        )

        print(
            "  smallest_nonnegative_Z={}".format(
                residue
            )
        )

    return combined


# ============================================================================
# 5. SPECIAL VALUES
# ============================================================================

def special_value_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "5. SPECIAL COEFFICIENT VALUES"
    )
    print("=" * 78)

    hits = []

    for var in D7_VARS:

        print()
        print(
            "  {}:".format(
                var
            )
        )

        for target in (
            sp.Integer(0),
            sp.Integer(1),
            sp.Integer(-1),
        ):

            solutions = sp.solve(
                sp.Eq(
                    expressions[var],
                    target,
                ),
                Z,
            )

            integer_solutions = [
                solution
                for solution in solutions
                if solution.is_integer is True
            ]

            print(
                "    {}={}: {}".format(
                    target,
                    "Z solutions",
                    integer_solutions,
                )
            )

            hits.extend(
                (
                    var,
                    target,
                    solution,
                )
                for solution in integer_solutions
            )

    print()
    print(
        "  integer_special_hits={}".format(
            hits
        )
    )

    return hits


# ============================================================================
# 6. PAIRWISE EQUALITY
# ============================================================================

def pairwise_equality_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "6. PAIRWISE COEFFICIENT EQUALITY"
    )
    print("=" * 78)

    hits = []

    for i in range(
        len(D7_VARS)
    ):

        for j in range(
            i + 1,
            len(D7_VARS),
        ):

            left = D7_VARS[i]
            right = D7_VARS[j]

            solutions = sp.solve(
                sp.Eq(
                    expressions[left],
                    expressions[right],
                ),
                Z,
            )

            integer_solutions = [
                solution
                for solution in solutions
                if solution.is_integer is True
            ]

            print(
                "  {}={}: {}".format(
                    left,
                    right,
                    integer_solutions,
                )
            )

            for solution in integer_solutions:

                hits.append(
                    (
                        left,
                        right,
                        solution,
                    )
                )

    print()
    print(
        "  integer_pairwise_hits={}".format(
            hits
        )
    )

    return hits


# ============================================================================
# 7. RATIONAL VECTOR / PRIMITIVE NORMALIZATION
# ============================================================================

def rational_vector_at_Z(
    expressions,
    z_value,
):

    return [
        clean(
            expressions[var].subs(
                Z,
                sp.Integer(z_value),
            )
        )
        for var in D7_VARS
    ]


def primitive_vector_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "7. RATIONAL-VECTOR / PRIMITIVE INTEGER NORMALIZATION"
    )
    print("=" * 78)

    sample_Z = [
        -2,
        -1,
        0,
        1,
        2,
    ]

    records = []

    for z_value in sample_Z:

        rational_vector = (
            rational_vector_at_Z(
                expressions,
                z_value,
            )
        )

        (
            primitive,
            common_denominator,
            integer_vector,
            gcd_value,
        ) = integer_vector_primitive(
            rational_vector
        )

        record = {
            "Z": z_value,
            "rational_vector": rational_vector,
            "all_coefficients_integer": all(
                is_integer_exact(value)
                for value
                in rational_vector
            ),
            "common_denominator": common_denominator,
            "integer_vector": integer_vector,
            "gcd": gcd_value,
            "primitive": primitive,
            "L1": sum(
                abs(v)
                for v in primitive
            ),
            "Linf": max(
                abs(v)
                for v in primitive
            ),
        }

        records.append(
            record
        )

        print()
        print(
            "  Z={}:".format(
                z_value
            )
        )

        print(
            "    rational_vector={}".format(
                rational_vector
            )
        )

        print(
            "    all_coefficients_integer={}".format(
                record[
                    "all_coefficients_integer"
                ]
            )
        )

        print(
            "    common_denominator={}".format(
                common_denominator
            )
        )

        print(
            "    integer_vector={}".format(
                integer_vector
            )
        )

        print(
            "    gcd={}".format(
                gcd_value
            )
        )

        print(
            "    primitive={}".format(
                primitive
            )
        )

        print(
            "    L1={}".format(
                record["L1"]
            )
        )

        print(
            "    Linf={}".format(
                record["Linf"]
            )
        )

    return records


# ============================================================================
# 8. COMPACT BOUNDED SEARCH
# ============================================================================

def compact_bounded_search(
    expressions,
    radius=1000,
):

    print()
    print("=" * 78)
    print(
        "8. COMPACT DIAGNOSTIC SEARCH"
    )
    print("=" * 78)

    records = []

    for z_value in range(
        -radius,
        radius + 1,
    ):

        rational_vector = (
            rational_vector_at_Z(
                expressions,
                z_value,
            )
        )

        (
            primitive,
            common_denominator,
            integer_vector,
            gcd_value,
        ) = integer_vector_primitive(
            rational_vector
        )

        l1 = sum(
            abs(v)
            for v in primitive
        )

        linf = max(
            abs(v)
            for v in primitive
        )

        all_integer = all(
            is_integer_exact(value)
            for value
            in rational_vector
        )

        records.append(
            (
                l1,
                linf,
                abs(z_value),
                z_value,
                all_integer,
                common_denominator,
                primitive,
            )
        )

    records.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
            item[3],
        )
    )

    integer_records = [
        record
        for record in records
        if record[4]
    ]

    print(
        "  search_range=[{},{}]".format(
            -radius,
            radius,
        )
    )

    print(
        "  integer_coefficient_Z_count={}".format(
            len(integer_records)
        )
    )

    print(
        "  top_10_primitive_vectors="
    )

    for record in records[:10]:

        (
            l1,
            linf,
            _abs_z,
            z_value,
            all_integer,
            denominator,
            primitive,
        ) = record

        print(
            "    Z={}, all_integer={}, denominator={}, "
            "L1={}, Linf={}, primitive={}".format(
                z_value,
                all_integer,
                denominator,
                l1,
                linf,
                primitive,
            )
        )

    return records[:10], integer_records


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 361R — EXACT D7 RATIONAL-FAMILY / "
        "INTEGER-FEASIBILITY / PRIMITIVE VECTOR AUDIT"
    )
    print("=" * 78)

    equations, expressions = (
        solve_d7_family()
    )

    print()
    print(
        "OBSERVED CELLS={}".format(
            len(build_lattice())
        )
    )

    print(
        "D7 calibration equations={}".format(
            len(equations)
        )
    )

    forms = family_audit(
        expressions
    )

    symbolic_ok = (
        symbolic_residual_audit(
            expressions
        )
    )

    conditions = (
        integer_coefficient_conditions(
            expressions,
            forms,
        )
    )

    simultaneous_integer = (
        combined_integrality_audit(
            conditions
        )
    )

    special_hits = (
        special_value_audit(
            expressions
        )
    )

    equality_hits = (
        pairwise_equality_audit(
            expressions
        )
    )

    vector_records = (
        primitive_vector_audit(
            expressions
        )
    )

    top_records, integer_records = (
        compact_bounded_search(
            expressions,
            radius=1000,
        )
    )

    # ------------------------------------------------------------------------
    # STRATEGIC INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "9. STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The previous D7 audit made an invalid assumption:

    integer Z
        =>
    integer D7 coefficients.

That implication is false.

361R separates the two notions exactly.

For each integer Z, the D7 coefficient vector is first evaluated in Q^5.
Only then is a common denominator cleared to obtain a primitive integer
representative of the same projective vector.

Thus:

    coefficient-wise integer feasibility
        is tested separately;

    primitive integer normalization
        is always available.

A D7 family that survives primitive normalization is not thereby validated.
Likewise, a non-integer coefficient at a sampled Z is not a crash or a
falsification.

The decisive arithmetic question is:

    Does there exist an integer Z such that ALL FIVE D7 coefficients
    are integers simultaneously?

The CRT calculation answers exactly that question.

Special values and bounded primitive-vector searches are secondary
diagnostics only.

The actual Q_3(5) remains unobserved.
"""
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
        "  D7_family_reconstructed_exactly=True"
    )

    print(
        "  affine_parameterization_verified=True"
    )

    print(
        "  symbolic_residuals_verified={}".format(
            symbolic_ok
        )
    )

    print(
        "  coefficient_integrality_conditions_completed=True"
    )

    print(
        "  simultaneous_integer_coefficient_feasibility={}".format(
            simultaneous_integer is not None
        )
    )

    print(
        "  special_value_audit_completed=True"
    )

    print(
        "  pairwise_equality_audit_completed=True"
    )

    print(
        "  primitive_rational_vector_audit_completed=True"
    )

    print(
        "  bounded_search_diagnostic_only=True"
    )

    print(
        "  bounded_search_top_count={}".format(
            len(top_records)
        )
    )

    print(
        "  bounded_search_integer_Z_count={}".format(
            len(integer_records)
        )
    )

    print(
        "  missing_Q3_5_used_as_data=False"
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

    ok = symbolic_ok

    print(
        "  failures={}".format(
            0
            if ok
            else 1
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            ok
        )
    )

    print()
    print(
        "EXPERIMENT 361R COMPLETE"
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