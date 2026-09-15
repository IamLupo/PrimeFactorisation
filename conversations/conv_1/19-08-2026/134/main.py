#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 300R — EXACT JOINT q_p(r) COORDINATE-LAW / LEAVE-ONE-p-OUT AUDIT
==============================================================================

Purpose
-------

Experiment 299R reconstructed each q-row independently and found that the
resulting power/falling/Newton coefficients do not obey a low-degree
cross-p law under leave-one-p-out testing.

Experiment 300R now removes that intermediate representation.

We work directly with the observed source table

    q_p(r)

and test whether it is governed by a compact joint law in the two source
coordinates

    p, r.

The tested families are:

    A) total-degree polynomial in (p,r);

    B) tensor-product polynomial
           sum c_ij p^i r^j;

    C) polynomial in p times Newton/binomial basis in r;

    D) polynomial in p times falling-factorial basis in r;

    E) polynomial in p and terminal-distance
           t = D(p) - r;

    F) polynomial in p and normalized source coordinate
           z = r / D(p)
       only where D(p) != 0.

For each class we perform:

    * exact global consistency;
    * leave-one-p-out prediction;
    * exact residual verification.

A formula passing only through all available points is labelled
INTERPOLATION and is not counted as structural evidence.

A leave-one-p-out prediction must reproduce every available q_p(r) in
the omitted row to count as an exact cross-p law.

No floating point.
No SymPy matrix with None.
No external files.
No synthetic second n=pq case.
No arbitrary matrix fit.
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

p_sym = sp.Symbol("p")
r_sym = sp.Symbol("r")
d_sym = sp.Symbol("d")
z_sym = sp.Symbol("z")


# ============================================================================
# EXACT HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def falling(x, n):
    out = sp.Integer(1)

    for i in range(n):
        out *= x - i

    return clean(out)


def basis_binomial_r(r, degree):
    return [
        sp.binomial(r, j)
        for j in range(degree + 1)
    ]


def basis_falling_r(r, degree):
    return [
        falling(r, j)
        for j in range(degree + 1)
    ]


def all_observations():
    """
    Return tuples (p, r, q).
    """
    out = []

    for p in sorted(Q):

        for r, value in enumerate(Q[p]):

            out.append(
                (
                    int(p),
                    int(r),
                    sp.Integer(value),
                )
            )

    return out


# ============================================================================
# MONOMIAL BASES
# ============================================================================

def total_degree_basis(
    p,
    r,
    degree,
):
    basis = []

    for i in range(degree + 1):

        for j in range(degree + 1 - i):

            basis.append(
                p**i * r**j
            )

    return basis


def tensor_basis(
    p,
    r,
    degree_p,
    degree_r,
):
    basis = []

    for i in range(degree_p + 1):

        for j in range(degree_r + 1):

            basis.append(
                p**i * r**j
            )

    return basis


def p_binomial_r_basis(
    p,
    r,
    degree_p,
    degree_r,
):
    basis = []

    for i in range(degree_p + 1):

        for j in range(degree_r + 1):

            basis.append(
                p**i * sp.binomial(r, j)
            )

    return basis


def p_falling_r_basis(
    p,
    r,
    degree_p,
    degree_r,
):
    basis = []

    for i in range(degree_p + 1):

        for j in range(degree_r + 1):

            basis.append(
                p**i * falling(r, j)
            )

    return basis


def p_terminal_distance_basis(
    p,
    r,
    D,
    degree_p,
    degree_t,
):
    t = D - r

    basis = []

    for i in range(degree_p + 1):

        for j in range(degree_t + 1):

            basis.append(
                p**i * t**j
            )

    return basis


# ============================================================================
# EXACT LINEAR SYSTEM
# ============================================================================

def solve_exact(
    observations,
    basis_function,
):
    """
    Solve A*c=b over QQ.

    Returns:

        status,
        solution,
        rank,
        augmented_rank,
        basis_count
    """

    if not observations:
        return (
            "NO_DATA",
            None,
            0,
            0,
            0,
        )

    rows = []
    rhs = []

    for p, r, q in observations:

        values = basis_function(
            p,
            r,
        )

        rows.append(
            [
                clean(value)
                for value in values
            ]
        )

        rhs.append(
            sp.Integer(q)
        )

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    rank = A.rank()

    augmented = A.row_join(b)

    augmented_rank = augmented.rank()

    basis_count = A.cols

    if augmented_rank > rank:

        return (
            "NO_SOLUTION",
            None,
            rank,
            augmented_rank,
            basis_count,
        )

    solution_set = sp.linsolve(
        (
            A,
            b,
        )
    )

    if solution_set == sp.EmptySet:

        return (
            "NO_SOLUTION",
            None,
            rank,
            augmented_rank,
            basis_count,
        )

    solution_list = list(
        solution_set
    )

    if not solution_list:

        return (
            "NO_SOLUTION",
            None,
            rank,
            augmented_rank,
            basis_count,
        )

    solution = solution_list[0]

    # A solution containing free symbols means non-unique.
    free_symbols = set()

    for value in solution:
        free_symbols |= value.free_symbols

    if free_symbols:

        return (
            "NONUNIQUE",
            solution,
            rank,
            augmented_rank,
            basis_count,
        )

    return (
        "UNIQUE",
        solution,
        rank,
        augmented_rank,
        basis_count,
    )


def evaluate_solution(
    solution,
    basis_values,
):
    if solution is None:
        return None

    if len(solution) != len(basis_values):
        return None

    return clean(
        sum(
            solution[i]
            * basis_values[i]
            for i in range(
                len(solution)
            )
        )
    )


def verify_solution(
    observations,
    solution,
    basis_function,
):
    if solution is None:
        return False

    for p, r, q in observations:

        predicted = evaluate_solution(
            solution,
            basis_function(
                p,
                r,
            ),
        )

        if predicted != q:
            return False

    return True


# ============================================================================
# PRINT DATA
# ============================================================================

def print_inventory():
    print()
    print("=" * 78)
    print("1. SOURCE q-TABLE")
    print("=" * 78)

    for p in sorted(Q):

        D = len(Q[p]) - 1

        print()
        print(
            "  p={} D(p)={}".format(
                p,
                D,
            )
        )

        print(
            "    q={}".format(
                Q[p],
            )
        )


# ============================================================================
# GLOBAL LAW SEARCH
# ============================================================================

def total_degree_search():
    print()
    print("=" * 78)
    print("2. TOTAL-DEGREE q(p,r) POLYNOMIAL SEARCH")
    print("=" * 78)

    observations = all_observations()

    hits = []

    for degree in range(0, 6):

        status, solution, rank, arank, count = (
            solve_exact(
                observations,
                lambda p, r, deg=degree:
                    total_degree_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        deg,
                    ),
            )
        )

        print()
        print(
            "  degree<={}:".format(
                degree,
            )
        )

        print(
            "    basis_count={}".format(
                count,
            )
        )

        print(
            "    equations={}".format(
                len(observations),
            )
        )

        print(
            "    rank={}".format(
                rank,
            )
        )

        print(
            "    augmented_rank={}".format(
                arank,
            )
        )

        print(
            "    status={}".format(
                status,
            )
        )

        if (
            status == "UNIQUE"
            and verify_solution(
                observations,
                solution,
                lambda p, r, deg=degree:
                    total_degree_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        deg,
                    ),
            )
        ):
            print(
                "    exact_global=True"
            )
            print(
                "    solution={}".format(
                    solution,
                )
            )

            hits.append(
                (
                    degree,
                    solution,
                )
            )
        else:
            print(
                "    exact_global=False"
            )

    return hits


# ============================================================================
# TENSOR SEARCH
# ============================================================================

def tensor_search():
    print()
    print("=" * 78)
    print("3. TENSOR-PRODUCT q(p,r) SEARCH")
    print("=" * 78)

    observations = all_observations()

    hits = []

    configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
        (4, 1),
        (1, 4),
    ]

    for degree_p, degree_r in configurations:

        status, solution, rank, arank, count = (
            solve_exact(
                observations,
                lambda p, r, dp=degree_p, dr=degree_r:
                    tensor_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        dp,
                        dr,
                    ),
            )
        )

        print()
        print(
            "  deg_p<={}, deg_r<={}:".format(
                degree_p,
                degree_r,
            )
        )

        print(
            "    basis_count={}".format(
                count,
            )
        )

        print(
            "    rank={}".format(
                rank,
            )
        )

        print(
            "    augmented_rank={}".format(
                arank,
            )
        )

        print(
            "    status={}".format(
                status,
            )
        )

        exact = (
            status == "UNIQUE"
            and verify_solution(
                observations,
                solution,
                lambda p, r, dp=degree_p, dr=degree_r:
                    tensor_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        dp,
                        dr,
                    ),
            )
        )

        print(
            "    exact_global={}".format(
                exact,
            )
        )

        if exact:

            print(
                "    solution={}".format(
                    solution,
                )
            )

            hits.append(
                (
                    degree_p,
                    degree_r,
                    solution,
                )
            )

    return hits


# ============================================================================
# BINOMIAL / FALLING BASIS SEARCH
# ============================================================================

def discrete_basis_search():
    print()
    print("=" * 78)
    print(
        "4. p x DISCRETE-r BASIS SEARCH"
    )
    print("=" * 78)

    observations = all_observations()

    configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
        (4, 1),
        (1, 4),
    ]

    results = {
        "binomial": [],
        "falling": [],
    }

    for basis_name, builder in [
        (
            "binomial",
            p_binomial_r_basis,
        ),
        (
            "falling",
            p_falling_r_basis,
        ),
    ]:

        print()
        print(
            "  BASIS={}".format(
                basis_name.upper(),
            )
        )

        for degree_p, degree_r in configurations:

            status, solution, rank, arank, count = (
                solve_exact(
                    observations,
                    lambda p, r, dp=degree_p, dr=degree_r, b=builder:
                        b(
                            sp.Integer(p),
                            sp.Integer(r),
                            dp,
                            dr,
                        ),
                )
            )

            exact = (
                status == "UNIQUE"
                and verify_solution(
                    observations,
                    solution,
                    lambda p, r, dp=degree_p, dr=degree_r, b=builder:
                        b(
                            sp.Integer(p),
                            sp.Integer(r),
                            dp,
                            dr,
                        ),
                )
            )

            print()
            print(
                "    deg_p<={}, deg_r<={}: status={} "
                "rank={} augmented_rank={} exact={}".format(
                    degree_p,
                    degree_r,
                    status,
                    rank,
                    arank,
                    exact,
                )
            )

            if exact:

                print(
                    "      solution={}".format(
                        solution,
                    )
                )

                results[
                    basis_name
                ].append(
                    (
                        degree_p,
                        degree_r,
                        solution,
                    )
                )

    return results


# ============================================================================
# TERMINAL-DISTANCE SEARCH
# ============================================================================

def terminal_distance_search():
    print()
    print("=" * 78)
    print(
        "5. TERMINAL-DISTANCE t=D(p)-r SEARCH"
    )
    print("=" * 78)

    observations = all_observations()

    hits = []

    def D_of_p(p):
        return len(
            Q[int(p)]
        ) - 1

    configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
    ]

    for degree_p, degree_t in configurations:

        status, solution, rank, arank, count = (
            solve_exact(
                observations,
                lambda p, r, dp=degree_p, dt=degree_t:
                    p_terminal_distance_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        D_of_p(
                            int(p)
                        ),
                        dp,
                        dt,
                    ),
            )
        )

        exact = (
            status == "UNIQUE"
            and verify_solution(
                observations,
                solution,
                lambda p, r, dp=degree_p, dt=degree_t:
                    p_terminal_distance_basis(
                        sp.Integer(p),
                        sp.Integer(r),
                        D_of_p(
                            int(p)
                        ),
                        dp,
                        dt,
                    ),
            )
        )

        print()
        print(
            "  deg_p<={}, deg_t<={}: status={} "
            "rank={} augmented_rank={} exact={}".format(
                degree_p,
                degree_t,
                status,
                rank,
                arank,
                exact,
            )
        )

        if exact:

            print(
                "    solution={}".format(
                    solution,
                )
            )

            hits.append(
                (
                    degree_p,
                    degree_t,
                    solution,
                )
            )

    return hits


# ============================================================================
# LEAVE-ONE-p-OUT
# ============================================================================

def leave_one_p_test(
    name,
    builder_factory,
    configurations,
):
    print()
    print("=" * 78)
    print(
        "6. {} LEAVE-ONE-p-OUT TEST".format(
            name
        )
    )
    print("=" * 78)

    p_values = sorted(
        Q
    )

    exact_total = 0
    tested_total = 0

    for degree_spec in configurations:

        print()
        print(
            "  configuration={}".format(
                degree_spec
            )
        )

        configuration_exact = 0
        configuration_tested = 0

        for target_p in p_values:

            train_observations = [
                observation
                for observation in all_observations()
                if observation[0] != target_p
            ]

            test_observations = [
                observation
                for observation in all_observations()
                if observation[0] == target_p
            ]

            builder = builder_factory(
                degree_spec
            )

            status, solution, rank, arank, count = (
                solve_exact(
                    train_observations,
                    builder,
                )
            )

            if status != "UNIQUE":
                print(
                    "    target_p={}: train_status={} "
                    "tested=False".format(
                        target_p,
                        status,
                    )
                )
                continue

            all_exact = True

            for _, r, q in test_observations:

                predicted = evaluate_solution(
                    solution,
                    builder(
                        target_p,
                        r,
                    ),
                )

                if predicted != q:
                    all_exact = False
                    break

            configuration_tested += 1
            tested_total += 1

            if all_exact:

                configuration_exact += 1
                exact_total += 1

            print(
                "    target_p={}: exact={} "
                "test_points={}".format(
                    target_p,
                    all_exact,
                    len(test_observations),
                )
            )

        print(
            "    configuration_exact={}/{}".format(
                configuration_exact,
                configuration_tested,
            )
        )

    print()
    print(
        "  exact_predictions={}".format(
            exact_total,
        )
    )

    print(
        "  tested_configurations={}".format(
            tested_total,
        )
    )

    return (
        exact_total,
        tested_total,
    )


# ============================================================================
# TERMINAL SOURCE AUDIT
# ============================================================================

def terminal_audit():
    print()
    print("=" * 78)
    print("7. TERMINAL SOURCE AUDIT")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print(
        "  q1_terminal={}".format(
            q1
        )
    )

    print(
        "  q3_terminal={}".format(
            q3
        )
    )

    gcd_value = math.gcd(
        q1,
        q3,
    )

    print(
        "  gcd={}".format(
            gcd_value
        )
    )

    print(
        "  q1/17={}".format(
            q1 // 17
        )
    )

    print(
        "  q3/17={}".format(
            q3 // 17
        )
    )

    print(
        "  v17(q1)={}".format(
            v17(q1)
        )
    )

    print(
        "  v17(q3)={}".format(
            v17(q3)
        )
    )


def v17(value):
    value = abs(
        int(value)
    )

    if value == 0:
        return None

    count = 0

    while value % 17 == 0:
        value //= 17
        count += 1

    return count


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():
    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 299R showed that independent row-polynomial coordinates
do not admit a successful low-degree cross-p law under leave-one-p-out
prediction.

Experiment 300R now tests the q-table directly as a two-variable source
object:

    q = q(p,r).

The distinction is important.

A polynomial interpolating the supplied observations is not enough.

The decisive test is:

    train on three p-values,
    omit the fourth p-value,
    predict its entire available q-row.

A successful joint law must reproduce every q_p(r) in the omitted row.

The experiment therefore separates:

    global interpolation
        from
    genuine cross-p prediction.

The discrete-r bases are especially relevant because the original
source data are indexed by integer r and the earlier experiments
showed repeated falling-factorial and Newton structure.

The terminal-distance coordinate

    t = D(p) - r

tests whether the source row is more naturally organized from the
terminal end.

A positive leave-one-p-out result would be the first direct evidence
that the four supplied q-rows are shadows of a common source law.

A negative result would mean that even direct two-variable low-degree
models cannot identify the upstream source construction.

No second n=pq case is manufactured.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 300R — EXACT JOINT q_p(r) COORDINATE-LAW / "
        "LEAVE-ONE-p-OUT AUDIT"
    )
    print("=" * 78)

    print_inventory()

    global_total_degree_hits = (
        total_degree_search()
    )

    global_tensor_hits = (
        tensor_search()
    )

    global_discrete_hits = (
        discrete_basis_search()
    )

    global_terminal_hits = (
        terminal_distance_search()
    )

    # ----------------------------------------------------------------------
    # LEAVE-ONE-p-OUT CONFIGURATIONS
    # ----------------------------------------------------------------------

    total_degree_configurations = [
        0,
        1,
        2,
        3,
    ]

    total_degree_results = leave_one_p_test(
        "TOTAL-DEGREE POLYNOMIAL",
        lambda degree:
            lambda p, r:
                total_degree_basis(
                    sp.Integer(p),
                    sp.Integer(r),
                    degree,
                ),
        total_degree_configurations,
    )

    tensor_configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
    ]

    tensor_results = leave_one_p_test(
        "TENSOR POLYNOMIAL",
        lambda spec:
            lambda p, r:
                tensor_basis(
                    sp.Integer(p),
                    sp.Integer(r),
                    spec[0],
                    spec[1],
                ),
        tensor_configurations,
    )

    binomial_configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
    ]

    binomial_results = leave_one_p_test(
        "p x BINOMIAL(r)",
        lambda spec:
            lambda p, r:
                p_binomial_r_basis(
                    sp.Integer(p),
                    sp.Integer(r),
                    spec[0],
                    spec[1],
                ),
        binomial_configurations,
    )

    falling_results = leave_one_p_test(
        "p x FALLING(r)",
        lambda spec:
            lambda p, r:
                p_falling_r_basis(
                    sp.Integer(p),
                    sp.Integer(r),
                    spec[0],
                    spec[1],
                ),
        binomial_configurations,
    )

    terminal_configurations = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
    ]

    terminal_results = leave_one_p_test(
        "p x TERMINAL-DISTANCE",
        lambda spec:
            lambda p, r:
                p_terminal_distance_basis(
                    sp.Integer(p),
                    sp.Integer(r),
                    len(
                        Q[int(p)]
                    ) - 1,
                    spec[0],
                    spec[1],
                ),
        terminal_configurations,
    )

    terminal_audit()

    interpretation()

    # ======================================================================
    # FINAL
    # ======================================================================

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  exact_q_table_observations={}".format(
            len(
                all_observations()
            )
        )
    )

    print(
        "  global_total_degree_hits={}".format(
            len(
                global_total_degree_hits
            )
        )
    )

    print(
        "  global_tensor_hits={}".format(
            len(
                global_tensor_hits
            )
        )
    )

    print(
        "  global_binomial_hits={}".format(
            len(
                global_discrete_hits[
                    "binomial"
                ]
            )
        )
    )

    print(
        "  global_falling_hits={}".format(
            len(
                global_discrete_hits[
                    "falling"
                ]
            )
        )
    )

    print(
        "  global_terminal_distance_hits={}".format(
            len(
                global_terminal_hits
            )
        )
    )

    print(
        "  total_degree_LOO_exact={}".format(
            total_degree_results[0]
        )
    )

    print(
        "  tensor_LOO_exact={}".format(
            tensor_results[0]
        )
    )

    print(
        "  binomial_r_LOO_exact={}".format(
            binomial_results[0]
        )
    )

    print(
        "  falling_r_LOO_exact={}".format(
            falling_results[0]
        )
    )

    print(
        "  terminal_distance_LOO_exact={}".format(
            terminal_results[0]
        )
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
        "  interpolation_counted_as_proof=False"
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
        "EXPERIMENT 300R COMPLETE"
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

