#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 359R — EXACT D7 PRIMITIVE-NORMALIZATION / PARAMETER-SELECTION AUDIT
==============================================================================

Purpose
-------
358R established:

    D5 -> forced non-integer Z
    D6 -> forced non-integer Z
    D7 -> integer coefficients possible for every integer Z

Therefore integrality alone cannot distinguish the D7 family.

359R asks whether the D7 coefficient family contains a distinguished
integer parameter value for any exact, non-arbitrary normalization such as:

    * primitive coefficient vector;
    * coefficient gcd;
    * vanishing of one coefficient;
    * coefficient = ±1;
    * repeated coefficient;
    * exact sign-pattern transitions;
    * common divisibility of coefficient differences;
    * exceptionally small primitive coefficient vectors for a bounded
      diagnostic search.

The bounded search is explicitly diagnostic only.

No value of Z is promoted to observed data.

No missing cell is inserted numerically.
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


# ============================================================================
# SYMBOLS
# ============================================================================

Z = sp.symbols("Z")

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


def coefficient_vector_at(
    expressions,
    z_value,
):
    return [
        sp.Integer(
            sp.Rational(
                clean(
                    expressions[var].subs(
                        Z,
                        z_value,
                    )
                )
            )
        )
        for var in D7_VARS
    ]


def integer_gcd(values):
    g = 0

    for value in values:

        value = int(value)

        g = math.gcd(
            g,
            abs(value),
        )

    return abs(g)


def primitive_vector(values):
    values = [
        int(value)
        for value in values
    ]

    g = integer_gcd(
        values
    )

    if g == 0:
        return values, 0

    primitive = [
        value // g
        for value in values
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

    return primitive, g


def l1_norm(values):
    return sum(
        abs(int(value))
        for value in values
    )


def linf_norm(values):
    return max(
        abs(int(value))
        for value in values
    )


# ============================================================================
# D7 FAMILY
# ============================================================================

def d7_family():

    equations = []

    lattice = build_lattice()

    def value(cell):
        if cell == (2, 3):
            return Z

        return lattice.get(
            cell,
            None,
        )

    def sources(r, t):
        return [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
            (r - 1, t - 2),
        ]

    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            target_value = value(
                target
            )

            if target_value is None:
                continue

            if r < 1:
                continue

            source_cells = sources(
                r,
                t,
            )

            source_values = [
                value(cell)
                for cell in source_cells
            ]

            if not all(
                source_value is not None
                for source_value
                in source_values
            ):
                continue

            equation = sp.Eq(
                (
                    a * source_values[0]
                    + b * source_values[1]
                    + c * source_values[2]
                    + d * source_values[3]
                    + e * source_values[4]
                ),
                target_value,
            )

            equations.append(
                equation
            )

    solution = sp.linsolve(
        equations,
        D7_VARS,
    )

    tuples = list(
        solution
    )

    if len(tuples) != 1:
        raise RuntimeError(
            "D7 solution family is not uniquely parameterized over Q(Z)."
        )

    tuple_solution = tuples[0]

    expressions = {
        var: clean(
            tuple_solution[i]
        )
        for i, var in enumerate(
            D7_VARS
        )
    }

    return equations, expressions


# ============================================================================
# EXACT AFFINE FORM
# ============================================================================

def affine_parts(
    expression
):
    expression = clean(
        expression
    )

    numerator, denominator = (
        sp.fraction(
            sp.together(
                expression
            )
        )
    )

    numerator = sp.Poly(
        numerator,
        Z,
        domain=sp.QQ,
    )

    if numerator.degree() > 1:
        raise RuntimeError(
            "D7 coefficient is not affine in Z."
        )

    A = sp.Integer(
        numerator.coeff_monomial(Z)
    )

    B = sp.Integer(
        numerator.coeff_monomial(1)
    )

    D = sp.Integer(
        denominator
    )

    g = math.gcd(
        math.gcd(
            abs(int(A)),
            abs(int(B)),
        ),
        abs(int(D)),
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
        int(A),
        int(B),
        int(D),
    )


# ============================================================================
# SECTION 1 — FAMILY
# ============================================================================

def print_family(expressions):

    print()
    print("=" * 78)
    print(
        "1. EXACT D7 INTEGER-AFFINE FAMILY"
    )
    print("=" * 78)

    for var in D7_VARS:

        expression = expressions[
            var
        ]

        A, B, D = affine_parts(
            expression
        )

        print()
        print(
            "  {}(Z)={}".format(
                var,
                expression,
            )
        )

        print(
            "    normalized_affine_form="
            "({}*Z + {})/{}".format(
                A,
                B,
                D,
            )
        )

        print(
            "    denominator={}".format(
                D
            )
        )

        print(
            "    slope={}".format(
                sp.Rational(
                    A,
                    D,
                )
            )
        )

        print(
            "    intercept={}".format(
                sp.Rational(
                    B,
                    D,
                )
            )
        )


# ============================================================================
# SECTION 2 — EXACT ZERO / ±1 CONDITIONS
# ============================================================================

def special_value_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "2. EXACT SPECIAL-PARAMETER AUDIT"
    )
    print("=" * 78)

    for var in D7_VARS:

        expression = expressions[
            var
        ]

        print()
        print(
            "  {}:".format(
                var
            )
        )

        for target in (
            0,
            1,
            -1,
        ):

            equation = sp.Eq(
                expression,
                target,
            )

            solutions = sp.solve(
                equation,
                Z,
            )

            integer_solutions = [
                clean(solution)
                for solution in solutions
                if solution.is_integer
            ]

            print(
                "    {}={} : rational_solutions={}, "
                "integer_solutions={}".format(
                    var,
                    target,
                    solutions,
                    integer_solutions,
                )
            )


# ============================================================================
# SECTION 3 — PAIRWISE EQUALITY
# ============================================================================

def pairwise_equalities(
    expressions
):

    print()
    print("=" * 78)
    print(
        "3. EXACT COEFFICIENT-EQUALITY AUDIT"
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

            v1 = D7_VARS[i]
            v2 = D7_VARS[j]

            equation = sp.Eq(
                expressions[v1],
                expressions[v2],
            )

            solutions = sp.solve(
                equation,
                Z,
            )

            integer_solutions = [
                clean(solution)
                for solution in solutions
                if solution.is_integer
            ]

            if integer_solutions:

                hits.append(
                    (
                        v1,
                        v2,
                        integer_solutions,
                    )
                )

            print()
            print(
                "  {}={}:".format(
                    v1,
                    v2,
                )
            )

            print(
                "    rational_solutions={}".format(
                    solutions
                )
            )

            print(
                "    integer_solutions={}".format(
                    integer_solutions
                )
            )

    print()
    print(
        "  pairwise_integer_equality_hits={}".format(
            hits
        )
    )

    return hits


# ============================================================================
# SECTION 4 — PRIMITIVE VECTORS
# ============================================================================

def primitive_vector_audit(
    expressions
):

    print()
    print("=" * 78)
    print(
        "4. PRIMITIVE INTEGER COEFFICIENT VECTOR AUDIT"
    )
    print("=" * 78)

    #
    # Because 358R found each coefficient integer for integer Z,
    # evaluate exact integer coefficient vectors.
    #
    sample_values = [
        -10,
        -1,
        0,
        1,
        2,
        10,
    ]

    records = []

    for z_value in sample_values:

        vector = coefficient_vector_at(
            expressions,
            z_value,
        )

        primitive, gcd_value = (
            primitive_vector(
                vector
            )
        )

        record = (
            z_value,
            vector,
            gcd_value,
            primitive,
            l1_norm(
                primitive
            ),
            linf_norm(
                primitive
            ),
        )

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
            "    coefficient_vector={}".format(
                vector
            )
        )

        print(
            "    coefficient_gcd={}".format(
                gcd_value
            )
        )

        print(
            "    primitive_vector={}".format(
                primitive
            )
        )

        print(
            "    primitive_L1={}".format(
                l1_norm(
                    primitive
                )
            )
        )

        print(
            "    primitive_Linf={}".format(
                linf_norm(
                    primitive
                )
            )
        )

    return records


# ============================================================================
# SECTION 5 — BOUNDED DIAGNOSTIC SEARCH
# ============================================================================

def bounded_search(
    expressions,
    radius=1000,
):

    print()
    print("=" * 78)
    print(
        "5. BOUNDED INTEGER-Z PRIMITIVE-NORM SEARCH"
    )
    print("=" * 78)

    print(
        "  search_range=[{},{}]".format(
            -radius,
            radius,
        )
    )

    best_l1 = None
    best_linf = None

    best_l1_records = []
    best_linf_records = []

    for z_value in range(
        -radius,
        radius + 1,
    ):

        vector = coefficient_vector_at(
            expressions,
            z_value,
        )

        primitive, gcd_value = (
            primitive_vector(
                vector
            )
        )

        l1 = l1_norm(
            primitive
        )

        linf = linf_norm(
            primitive
        )

        record = {
            "Z": z_value,
            "primitive": primitive,
            "gcd": gcd_value,
            "L1": l1,
            "Linf": linf,
        }

        if (
            best_l1 is None
            or l1 < best_l1
        ):

            best_l1 = l1
            best_l1_records = [
                record
            ]

        elif l1 == best_l1:

            best_l1_records.append(
                record
            )

        if (
            best_linf is None
            or linf < best_linf
        ):

            best_linf = linf
            best_linf_records = [
                record
            ]

        elif linf == best_linf:

            best_linf_records.append(
                record
            )

    print()
    print(
        "  best_L1={}".format(
            best_l1
        )
    )

    print(
        "  best_L1_records={}".format(
            best_l1_records
        )
    )

    print()
    print(
        "  best_Linf={}".format(
            best_linf
        )
    )

    print(
        "  best_Linf_records={}".format(
            best_linf_records
        )
    )

    return (
        best_l1_records,
        best_linf_records,
    )


# ============================================================================
# SECTION 6 — GCD STRUCTURE
# ============================================================================

def gcd_structure(
    expressions,
    radius=100,
):

    print()
    print("=" * 78)
    print(
        "6. COEFFICIENT-GCD PARAMETER AUDIT"
    )
    print("=" * 78)

    gcd_values = {}

    for z_value in range(
        -radius,
        radius + 1,
    ):

        vector = coefficient_vector_at(
            expressions,
            z_value,
        )

        gcd_values[
            z_value
        ] = integer_gcd(
            vector
        )

    distinct_gcds = sorted(
        set(
            gcd_values.values()
        )
    )

    print(
        "  search_range=[{},{}]".format(
            -radius,
            radius,
        )
    )

    print(
        "  distinct_gcd_values={}".format(
            distinct_gcds
        )
    )

    gcd_one_values = [
        z_value
        for (
            z_value,
            gcd_value,
        ) in gcd_values.items()
        if gcd_value == 1
    ]

    print(
        "  gcd_one_Z_values={}".format(
            gcd_one_values
        )
    )

    return gcd_values


# ============================================================================
# SECTION 7 — PARAMETERIZED RESIDUAL CHECK
# ============================================================================

def symbolic_verification(
    expressions,
):

    print()
    print("=" * 78)
    print(
        "7. SYMBOLIC D7 FAMILY RESIDUAL AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    equations = []

    def value(cell):

        if cell == (2, 3):
            return Z

        return lattice.get(
            cell,
            None,
        )

    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            target_value = value(
                target
            )

            if target_value is None:
                continue

            if r < 1:
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
                source is not None
                for source in source_values
            ):
                continue

            lhs = (
                expressions[a]
                * source_values[0]
                + expressions[b]
                * source_values[1]
                + expressions[c]
                * source_values[2]
                + expressions[d]
                * source_values[3]
                + expressions[e]
                * source_values[4]
            )

            equations.append(
                (
                    target,
                    clean(
                        lhs
                        - target_value
                    ),
                )
            )

    failures = [
        item
        for item in equations
        if item[1] != 0
    ]

    for target, residual in equations:

        print(
            "  target={}: residual={}".format(
                target,
                residual,
            )
        )

    print()
    print(
        "  all_symbolic_residuals_zero={}".format(
            len(failures) == 0
        )
    )

    return len(failures) == 0


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 359R — EXACT D7 PRIMITIVE-NORMALIZATION / "
        "PARAMETER-SELECTION AUDIT"
    )
    print("=" * 78)

    equations, expressions = d7_family()

    print_family(
        expressions
    )

    special_value_audit(
        expressions
    )

    equality_hits = pairwise_equalities(
        expressions
    )

    primitive_records = primitive_vector_audit(
        expressions
    )

    best_l1, best_linf = bounded_search(
        expressions,
        radius=1000,
    )

    gcd_values = gcd_structure(
        expressions,
        radius=100,
    )

    symbolic_ok = symbolic_verification(
        expressions
    )

    # ------------------------------------------------------------------------
    # Strategic interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "8. STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
358R showed that D7 survives integer-coefficient compatibility for every
integer value of Z.

Therefore there is no reason to repeat a generic integrality audit.

359R asks whether D7 nevertheless contains a canonical parameter value.

The exact tests are:

    * coefficient = 0 or ±1;
    * pairwise coefficient equality;
    * primitive integer normalization;
    * coefficient gcd;
    * bounded small-norm search.

The bounded search is deliberately labeled diagnostic.

It does NOT prove that the smallest-norm Z is mathematically preferred.

The strongest positive result would be a parameter-independent exact
normalization condition selecting a unique integer Z.

The strongest negative result would be that the family remains generic:
no zero/±1 coefficient condition, no equality condition, and no isolated
primitive-minimum mechanism.

In that case D7 is simply another flexible family fitted to the current
five symbolic equations.

The actual Q_3(5) remains the decisive independent datum.

No missing value is used as observed data.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "9. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  D7_family_built_exactly=True"
    )

    print(
        "  symbolic_residuals_verified={}".format(
            symbolic_ok
        )
    )

    print(
        "  special_value_audit_completed=True"
    )

    print(
        "  pairwise_equality_audit_completed=True"
    )

    print(
        "  primitive_normalization_audit_completed=True"
    )

    print(
        "  coefficient_gcd_audit_completed=True"
    )

    print(
        "  bounded_search_is_diagnostic_only=True"
    )

    print(
        "  missing_values_inserted_numerically=False"
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
            0
            if symbolic_ok
            else 1
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            symbolic_ok
        )
    )

    print()
    print(
        "EXPERIMENT 359R COMPLETE"
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
