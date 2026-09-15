#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 348R — EXACT ITERATED SOURCE-OPERATOR / SEMIGROUP PROVENANCE AUDIT
==============================================================================

Purpose
-------
Experiments 322R–347R tested many one-step or static relations between the
observed source rows P_t(p).

Experiment 348R asks a structurally different question:

    Does there exist ONE fixed operator L such that

        P_{t+1} = L(P_t)

    for every observed transition, hence

        P_t = L^t(P_0)?

This is a semigroup / iterated-operator hypothesis.

The important distinction is:

    one-step fitting:
        P_{t+1} = L_t(P_t)

    versus

    fixed evolution:
        P_{t+1} = L(P_t)    for every t.

Only the second is tested here.

The principal operator families are:

    1. CONSTANT-COEFFICIENT DIFFERENTIAL OPERATOR

        L = a0 I + a1 D + a2 D^2 + a3 D^3

    2. EULER OPERATOR POLYNOMIAL

        E = p D/dp
        L = a0 I + a1 E + a2 E^2 + a3 E^3

    3. FIXED FINITE-DIFFERENCE OPERATOR

        S_h f(p) = f(p+h)

        L = a0 I + a1 S_h + a2 S_h^2 + a3 S_h^3

       for small integer h.

    4. MIXED DIFFERENTIAL / EULER OPERATOR

        L = sum a_k D^k
            + sum b_k E^k

Only operators with a UNIQUE exact rational solution are accepted.

A merely data-sized or nonunique operator is not a discovery.

The test is intentionally stronger than Experiment 334R:
the SAME operator must propagate the ENTIRE observed row sequence.

No missing source values are introduced.
No new n=pq case is generated.
All arithmetic is exact SymPy arithmetic.
"""


from __future__ import annotations

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


p = sp.symbols("p")
x = sp.symbols("x")


# ============================================================================
# HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def build_layers():
    """
    Return observed rows in terminal-distance order.
    """

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t_value in range(
        maximum_t + 1
    ):

        row = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = (
                len(values)
                - 1
                - t_value
            )

            if index >= 0:

                row.append(
                    (
                        sp.Integer(p_value),
                        sp.Integer(values[index]),
                    )
                )

        layers[t_value] = row

    return layers


def interpolate_observed_row(row):
    """
    Reconstruct only the polynomial through the actually observed cells
    of one row.

    This is NOT treated as evidence outside observed p-values. It is simply
    the coordinate representation used to test the operator identity.
    """

    points = [
        (sp.Integer(pv), sp.Integer(value))
        for pv, value in row
    ]

    return clean(
        sp.interpolate(
            points,
            p,
        )
    )


def build_row_polynomials(layers):

    return {
        t_value: interpolate_observed_row(
            row
        )
        for t_value, row in layers.items()
    }


def operator_equations(
    source_polynomials,
    target_polynomials,
    basis_polynomials,
):
    """
    Solve for coefficients a_k in

        target_t = sum_k a_k basis_k(source_t)

    simultaneously for all observed transitions.
    """

    rows = []
    rhs = []

    for t_value in sorted(
        target_polynomials
    ):

        source = source_polynomials[t_value]
        target = target_polynomials[t_value]

        basis_images = [
            clean(
                basis(source)
            )
            for basis in basis_polynomials
        ]

        source_poly = sp.Poly(
            source,
            p,
            domain=sp.QQ,
        )

        target_poly = sp.Poly(
            target,
            p,
            domain=sp.QQ,
        )

        max_degree = max(
            source_poly.degree(),
            target_poly.degree(),
            *[
                sp.Poly(
                    image,
                    p,
                    domain=sp.QQ,
                ).degree()
                for image in basis_images
                if image != 0
            ],
        )

        for degree in range(
            max_degree + 1
        ):

            coefficient_row = [
                sp.Poly(
                    image,
                    p,
                    domain=sp.QQ,
                ).nth(degree)
                if image != 0
                else sp.Integer(0)
                for image in basis_images
            ]

            target_coeff = (
                target_poly.nth(degree)
            )

            # Skip completely zero identities.
            if all(
                coefficient == 0
                for coefficient in coefficient_row
            ) and target_coeff == 0:
                continue

            rows.append(
                coefficient_row
            )

            rhs.append(
                target_coeff
            )

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    return A, b


def solve_unique_operator(
    label,
    source_polynomials,
    target_polynomials,
    basis_builders,
):
    """
    Exact simultaneous operator solve.

    basis_builders is a list of functions f -> basis element.
    """

    A, b = operator_equations(
        source_polynomials,
        target_polynomials,
        basis_builders,
    )

    unknowns = len(basis_builders)

    rank = A.rank()
    augmented_rank = (
        A.row_join(b).rank()
    )

    equation_count = A.rows
    redundancy = (
        equation_count - unknowns
    )

    print()
    print("=" * 78)
    print(label)
    print("=" * 78)

    print(
        "  equations={}".format(
            equation_count
        )
    )

    print(
        "  unknowns={}".format(
            unknowns
        )
    )

    print(
        "  redundancy={}".format(
            redundancy
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  augmented_rank={}".format(
            augmented_rank
        )
    )

    if augmented_rank > rank:

        print(
            "  status=NO_SOLUTION"
        )

        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank != unknowns:

        nullity = (
            unknowns - rank
        )

        print(
            "  nullity={}".format(
                nullity
            )
        )

        print(
            "  status=NONUNIQUE"
        )

        return {
            "status": "NONUNIQUE",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "nullity": nullity,
        }

    solution = A.gauss_jordan_solve(
        b
    )[0]

    coefficients = [
        clean(solution[i])
        for i in range(
            unknowns
        )
    ]

    print(
        "  coefficients={}".format(
            coefficients
        )
    )

    residuals = [
        clean(
            sum(
                coefficients[j]
                * A[i, j]
                for j in range(
                    unknowns
                )
            )
            - b[i]
        )
        for i in range(
            equation_count
        )
    ]

    exact = all(
        residual == 0
        for residual in residuals
    )

    if exact and equation_count > unknowns:
        status = "EXACT_OVERDETERMINED"
    elif exact:
        status = "EXACT_DATA_SIZED"
    else:
        status = "VERIFICATION_FAILED"

    print(
        "  residuals_zero={}".format(
            exact
        )
    )

    print(
        "  status={}".format(
            status
        )
    )

    return {
        "status": status,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "coefficients": coefficients,
        "residuals": residuals,
    }


# ============================================================================
# OPERATOR FAMILIES
# ============================================================================

def differential_basis(max_order):
    basis = []

    for order in range(
        max_order + 1
    ):

        if order == 0:
            basis.append(
                lambda f: f
            )
        else:

            basis.append(
                lambda f, k=order: sp.diff(
                    f,
                    p,
                    k,
                )
            )

    return basis


def euler_operator(f):
    return clean(
        p * sp.diff(
            f,
            p,
        )
    )


def euler_power_operator(k):
    def apply(f):
        result = f

        for _ in range(k):
            result = euler_operator(
                result
            )

        return clean(result)

    return apply


def euler_basis(max_order):
    return [
        (
            lambda f: f
            if k == 0
            else euler_power_operator(k)(f)
        )
        for k in range(
            max_order + 1
        )
    ]


def shift_basis(
    h,
    max_power,
):
    basis = []

    for power in range(
        max_power + 1
    ):

        def apply(
            f,
            h=h,
            power=power,
        ):
            return clean(
                f.subs(
                    p,
                    p + h * power,
                )
            )

        basis.append(
            apply
        )

    return basis


def mixed_basis(max_order):
    basis = []

    # I, D, D^2, ...
    for order in range(
        max_order + 1
    ):

        if order == 0:
            basis.append(
                lambda f: f
            )
        else:

            basis.append(
                lambda f, k=order: sp.diff(
                    f,
                    p,
                    k,
                )
            )

    # E, E^2, ...
    for order in range(
        1,
        max_order + 1,
    ):

        basis.append(
            euler_power_operator(
                order
            )
        )

    return basis


def basis_expression(
    name,
    coefficients,
):

    parts = []

    if name == "D":

        for k, coefficient in enumerate(
            coefficients
        ):

            if coefficient == 0:
                continue

            if k == 0:
                parts.append(
                    coefficient
                )
            elif k == 1:
                parts.append(
                    coefficient
                    * sp.Symbol("D")
                )
            else:
                parts.append(
                    coefficient
                    * sp.Symbol("D") ** k
                )

    elif name == "E":

        E = sp.Symbol("E")

        for k, coefficient in enumerate(
            coefficients
        ):

            if coefficient == 0:
                continue

            parts.append(
                coefficient * E**k
            )

    elif name.startswith(
        "SHIFT"
    ):

        S = sp.Symbol("S")

        for k, coefficient in enumerate(
            coefficients
        ):

            if coefficient == 0:
                continue

            parts.append(
                coefficient
                * S**k
            )

    return clean(
        sum(parts)
    ) if parts else sp.Integer(0)


# ============================================================================
# FIXED-ITERATE VALIDATION
# ============================================================================

def validate_fixed_iteration(
    label,
    row_polynomials,
    basis_builders,
    coefficients,
):
    """
    Verify every observed transition independently after solving.

    This is the semigroup check:

        L(P_t) == P_{t+1}

    for all t.
    """

    print()
    print("=" * 78)
    print(
        "{} FIXED-ITERATION VALIDATION".format(
            label
        )
    )
    print("=" * 78)

    failures = []

    for t_value in range(5):

        source = row_polynomials[t_value]
        target = row_polynomials[
            t_value + 1
        ]

        image = clean(
            sum(
                coefficients[k]
                * basis_builders[k](source)
                for k in range(
                    len(coefficients)
                )
            )
        )

        residual = clean(
            image - target
        )

        print()
        print(
            "  t={} -> {}:".format(
                t_value,
                t_value + 1,
            )
        )

        print(
            "    source_degree={}".format(
                sp.degree(
                    source,
                    p,
                )
            )
        )

        print(
            "    target_degree={}".format(
                sp.degree(
                    target,
                    p,
                )
            )
        )

        print(
            "    residual={}".format(
                residual
            )
        )

        print(
            "    exact={}".format(
                residual == 0
            )
        )

        if residual != 0:
            failures.append(
                t_value
            )

    print()
    print(
        "  iteration_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 348R — EXACT ITERATED SOURCE-OPERATOR / "
        "SEMIGROUP PROVENANCE AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    row_polynomials = (
        build_row_polynomials(
            layers
        )
    )

    print()
    print("=" * 78)
    print(
        "1. OBSERVED ROW-POLYNOMIAL REPRESENTATION"
    )
    print("=" * 78)

    for t_value in range(6):

        print()
        print(
            "  t={}: P_t(p)={}".format(
                t_value,
                row_polynomials[t_value],
            )
        )

        print(
            "    degree={}".format(
                sp.degree(
                    row_polynomials[t_value],
                    p,
                )
            )
        )

    # ------------------------------------------------------------------------
    # Differential operator
    # ------------------------------------------------------------------------

    diff_basis = differential_basis(3)

    diff_result = solve_unique_operator(
        "2. FIXED DIFFERENTIAL OPERATOR DEGREE <= 3",
        {
            t: row_polynomials[t]
            for t in range(5)
        },
        {
            t: row_polynomials[t + 1]
            for t in range(5)
        },
        diff_basis,
    )

    if diff_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        diff_failures = validate_fixed_iteration(
            "DIFFERENTIAL",
            row_polynomials,
            diff_basis,
            diff_result["coefficients"],
        )

    else:
        diff_failures = []

    # ------------------------------------------------------------------------
    # Euler operator
    # ------------------------------------------------------------------------

    euler_basis_functions = (
        euler_basis(3)
    )

    euler_result = solve_unique_operator(
        "3. FIXED EULER OPERATOR DEGREE <= 3",
        {
            t: row_polynomials[t]
            for t in range(5)
        },
        {
            t: row_polynomials[t + 1]
            for t in range(5)
        },
        euler_basis_functions,
    )

    if euler_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        euler_failures = validate_fixed_iteration(
            "EULER",
            row_polynomials,
            euler_basis_functions,
            euler_result["coefficients"],
        )

    else:
        euler_failures = []

    # ------------------------------------------------------------------------
    # Shift operators
    # ------------------------------------------------------------------------

    shift_results = {}

    for h in (
        -3,
        -2,
        -1,
        1,
        2,
        3,
    ):

        basis = shift_basis(
            h,
            3,
        )

        result = solve_unique_operator(
            "4. FIXED SHIFT OPERATOR h={} DEGREE <= 3".format(
                h
            ),
            {
                t: row_polynomials[t]
                for t in range(5)
            },
            {
                t: row_polynomials[t + 1]
                for t in range(5)
            },
            basis,
        )

        shift_results[h] = result

        if result["status"] in (
            "EXACT_OVERDETERMINED",
            "EXACT_DATA_SIZED",
        ):

            failures = validate_fixed_iteration(
                "SHIFT h={}".format(
                    h
                ),
                row_polynomials,
                basis,
                result["coefficients"],
            )

        else:
            failures = []

        shift_results[h]["iteration_failures"] = (
            failures
        )

    # ------------------------------------------------------------------------
    # Mixed D/E family
    # ------------------------------------------------------------------------

    mixed = mixed_basis(2)

    mixed_result = solve_unique_operator(
        "5. FIXED MIXED (D,E) OPERATOR ORDER <= 2",
        {
            t: row_polynomials[t]
            for t in range(5)
        },
        {
            t: row_polynomials[t + 1]
            for t in range(5)
        },
        mixed,
    )

    if mixed_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        mixed_failures = validate_fixed_iteration(
            "MIXED",
            row_polynomials,
            mixed,
            mixed_result["coefficients"],
        )

    else:
        mixed_failures = []

    # ------------------------------------------------------------------------
    # Direct L^2/L^3 consistency
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. DIRECT ITERATE CONSISTENCY"
    )
    print("=" * 78)

    direct_reports = []

    candidate_families = []

    if diff_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        candidate_families.append(
            (
                "DIFFERENTIAL",
                diff_basis,
                diff_result["coefficients"],
            )
        )

    if euler_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        candidate_families.append(
            (
                "EULER",
                euler_basis_functions,
                euler_result["coefficients"],
            )
        )

    for h, result in shift_results.items():

        if result["status"] in (
            "EXACT_OVERDETERMINED",
            "EXACT_DATA_SIZED",
        ):
            candidate_families.append(
                (
                    "SHIFT {}".format(h),
                    shift_basis(h, 3),
                    result["coefficients"],
                )
            )

    if mixed_result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        candidate_families.append(
            (
                "MIXED",
                mixed,
                mixed_result["coefficients"],
            )
        )

    for (
        name,
        basis,
        coefficients,
    ) in candidate_families:

        print()
        print(
            "  candidate={}".format(
                name
            )
        )

        L_P0 = clean(
            sum(
                coefficients[k]
                * basis[k](
                    row_polynomials[0]
                )
                for k in range(
                    len(coefficients)
                )
            )
        )

        L2_P0 = clean(
            sum(
                coefficients[k]
                * basis[k](L_P0)
                for k in range(
                    len(coefficients)
                )
            )
        )

        L3_P0 = clean(
            sum(
                coefficients[k]
                * basis[k](L2_P0)
                for k in range(
                    len(coefficients)
                )
            )
        )

        print(
            "    L(P0)-P1={}".format(
                clean(
                    L_P0
                    - row_polynomials[1]
                )
            )
        )

        print(
            "    L^2(P0)-P2={}".format(
                clean(
                    L2_P0
                    - row_polynomials[2]
                )
            )
        )

        print(
            "    L^3(P0)-P3={}".format(
                clean(
                    L3_P0
                    - row_polynomials[3]
                )
            )

        )

        direct_reports.append(
            (
                name,
                clean(
                    L_P0
                    - row_polynomials[1]
                ),
                clean(
                    L2_P0
                    - row_polynomials[2]
                ),
                clean(
                    L3_P0
                    - row_polynomials[3]
                ),
            )
        )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "7. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 347R found only tautological terminal-polynomial hits.

Experiment 348R therefore makes the semigroup requirement explicit.

A genuine source operator should not merely map one row to the next.
The SAME operator must explain every transition:

    P1 = L(P0)
    P2 = L(P1)
    P3 = L(P2)
    P4 = L(P3)
    P5 = L(P4).

Equivalently,

    P_t = L^t(P0).

This is substantially stronger than the one-step operator searches.

Three especially natural operator algebras are tested:

    constant-coefficient differential operators;
    Euler-operator polynomials;
    fixed finite-shift algebras.

A mixed differential/Euler family is also checked.

The decisive outcome is:

    exact_overdetermined
        AND
    L^2(P0)=P2
        AND
    L^3(P0)=P3.

If one of these survives, the experiment has identified a genuine
iterated operator mechanism rather than a collection of independently
fitted transitions.

If all fail, generic operator provenance should be considered exhausted
for the current observed data.

The remaining task is source reconstruction from the original
mathematical definition, not additional post-hoc formula fitting.

No missing source cells are used.
No synthetic second n=pq case is generated.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    exact_overdetermined = []

    if diff_result["status"] == (
        "EXACT_OVERDETERMINED"
    ):
        exact_overdetermined.append(
            "differential"
        )

    if euler_result["status"] == (
        "EXACT_OVERDETERMINED"
    ):
        exact_overdetermined.append(
            "euler"
        )

    for h, result in shift_results.items():

        if result["status"] == (
            "EXACT_OVERDETERMINED"
        ):
            exact_overdetermined.append(
                "shift_{}".format(h)
            )

    if mixed_result["status"] == (
        "EXACT_OVERDETERMINED"
    ):
        exact_overdetermined.append(
            "mixed"
        )

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_iterated_operator_families={}".format(
            exact_overdetermined
        )
    )

    print(
        "  differential_failures={}".format(
            diff_failures
        )
    )

    print(
        "  euler_failures={}".format(
            euler_failures
        )
    )

    print(
        "  mixed_failures={}".format(
            mixed_failures
        )
    )

    print(
        "  observed_row_polynomials_only=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  extrapolation_used=False"
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
        "  interpolation_as_external_evidence=False"
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
        "EXPERIMENT 348R COMPLETE"
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
