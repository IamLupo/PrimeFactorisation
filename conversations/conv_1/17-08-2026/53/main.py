"""
==============================================================================
KAPPA EXPERIMENT 120R-P
C/D RECURRENCE KERNEL
ROBUST QUADRATIC-BRANCH DECOMPOSITION
MODULO S^2-S-X
PRINT-SAFE DEGREE PROFILE
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sympy as sp
from dataclasses import dataclass


# ============================================================================
# SYMBOLS
# ============================================================================

N, S, X = sp.symbols("N S X")


# ============================================================================
# BASIC HELPERS
# ============================================================================

def canon(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def degree_string(expr, var):
    """
    Print-safe polynomial degree.

    SymPy returns -oo for the zero polynomial. Never apply a numeric
    format specification directly to SymPy's NegativeInfinity.
    """
    expr = sp.expand(expr)

    if expr == 0:
        return "-oo"

    poly = sp.Poly(
        expr,
        var,
    )

    if poly.is_zero:
        return "-oo"

    return str(poly.degree())


def degree_int(expr, var):
    """
    Numeric degree helper for comparisons.
    Returns None for the zero polynomial.
    """
    expr = sp.expand(expr)

    if expr == 0:
        return None

    poly = sp.Poly(
        expr,
        var,
    )

    if poly.is_zero:
        return None

    return int(poly.degree())


def term_count(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return 0

    poly = sp.Poly(
        expr,
        N,
        X,
    )

    return len(poly.terms())


# ============================================================================
# DIRECT DETECTOR
# ============================================================================

def detector_raw(k, ell):
    p, q = sp.symbols("p q")

    F = sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ArithmeticError(
            f"Symmetrization remainder for ({k},{ell}): "
            f"{canon(remainder)}"
        )

    if len(mapping) != 2:
        raise ArithmeticError(
            f"Unexpected symmetric mapping for ({k},{ell}): {mapping}"
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    expr = canon(
        sym_expr.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    check = sp.expand(
        expr.subs(
            {
                S: p + q,
                N: p * q,
            }
        ) - F
    )

    if check != 0:
        raise ArithmeticError(
            f"Symmetric rewrite failed for ({k},{ell}): "
            f"{canon(check)}"
        )

    return expr


# ============================================================================
# EXACT QUOTIENT
# ============================================================================

def quotient_exact(k, ell):
    F = detector_raw(k, ell)

    P_F = sp.Poly(
        F,
        S,
        domain=sp.QQ.frac_field(N),
    )

    P_D = sp.Poly(
        S + 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    Q_poly, R_poly = sp.div(
        P_F,
        P_D,
    )

    remainder = canon(
        R_poly.as_expr()
    )

    if remainder != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{remainder}"
        )

    Q = canon(
        Q_poly.as_expr()
    )

    residual = canon(
        F - (S + 1) * Q
    )

    if residual != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed ({k},{ell}): "
            f"{residual}"
        )

    return Q


# ============================================================================
# REDUCTION MODULO S^2-S-X
# ============================================================================

def reduce_quadratic_relation(expr):
    """
    Reduce modulo

        S^2 - S - X = 0.

    Returns

        A(N,X)*S + B(N,X).
    """

    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N, X),
    )

    modulus = sp.Poly(
        S**2 - S - X,
        S,
        domain=sp.QQ.frac_field(N, X),
    )

    _, remainder = sp.div(
        poly,
        modulus,
    )

    rem_expr = sp.expand(
        remainder.as_expr()
    )

    rem = sp.Poly(
        rem_expr,
        S,
        domain=sp.QQ.frac_field(N, X),
    )

    A = sp.expand(
        rem.coeff_monomial(S)
    )

    B = sp.expand(
        rem.coeff_monomial(1)
    )

    reconstructed = sp.expand(
        A*S + B
    )

    residual = canon(
        rem_expr - reconstructed
    )

    if residual != 0:
        raise ArithmeticError(
            "Quadratic reduction returned unexpected remainder."
        )

    return canon(A), canon(B)


# ============================================================================
# C/D DECOMPOSITION
# ============================================================================

def branch_decompose(Q):
    """
    Exact decomposition

        Q(N,S) = C(N,X) + (2S-1)D(N,X),

        X=S(S-1).
    """

    Q = sp.expand(Q)

    Q_reflected = sp.expand(
        Q.subs(
            S,
            1 - S,
        )
    )

    Q_plus = sp.expand(
        (Q + Q_reflected) / 2
    )

    Q_minus = sp.expand(
        (Q - Q_reflected) / 2
    )

    # ----------------------------------------------------------------------
    # C
    # ----------------------------------------------------------------------

    A_C, B_C = reduce_quadratic_relation(
        Q_plus
    )

    if canon(A_C) != 0:
        raise ArithmeticError(
            "Invariant component unexpectedly contains S:\n"
            f"{canon(A_C)}"
        )

    C = canon(B_C)

    # ----------------------------------------------------------------------
    # D
    # ----------------------------------------------------------------------

    numerator = sp.Poly(
        Q_minus,
        S,
        domain=sp.QQ.frac_field(N),
    )

    divisor = sp.Poly(
        2*S - 1,
        S,
        domain=sp.QQ.frac_field(N),
    )

    D_poly, D_rem = sp.div(
        numerator,
        divisor,
    )

    D_remainder = canon(
        D_rem.as_expr()
    )

    if D_remainder != 0:
        raise ArithmeticError(
            "Odd component not divisible by 2S-1:\n"
            f"{D_remainder}"
        )

    D_raw = sp.expand(
        D_poly.as_expr()
    )

    A_D, B_D = reduce_quadratic_relation(
        D_raw
    )

    if canon(A_D) != 0:
        raise ArithmeticError(
            "D component unexpectedly contains S:\n"
            f"{canon(A_D)}"
        )

    D = canon(B_D)

    # ----------------------------------------------------------------------
    # Exact reconstruction
    # ----------------------------------------------------------------------

    reconstruction = sp.expand(
        C.subs(
            X,
            S * (S - 1),
        )
        +
        (2*S - 1)
        *
        D.subs(
            X,
            S * (S - 1),
        )
    )

    residual = canon(
        Q - reconstruction
    )

    if residual != 0:
        raise ArithmeticError(
            "C/D reconstruction failed:\n"
            f"{residual}"
        )

    return C, D


# ============================================================================
# RECORD
# ============================================================================

@dataclass(frozen=True)
class CDRecord:
    k: int
    ell: int
    r: int
    m: int
    C: sp.Expr
    D: sp.Expr


# ============================================================================
# DATASET
# ============================================================================

def odd_weights(max_ell):
    return list(
        range(
            3,
            max_ell + 1,
            2,
        )
    )


def build_cd_dataset(max_ell):
    records = []

    for ell in odd_weights(max_ell):

        for k in range(
            1,
            ell,
            2,
        ):

            Q = quotient_exact(
                k,
                ell,
            )

            C, D = branch_decompose(
                Q
            )

            records.append(
                CDRecord(
                    k=k,
                    ell=ell,
                    r=(k - 1) // 2,
                    m=(ell - 1) // 2,
                    C=C,
                    D=D,
                )
            )

    return records


# ============================================================================
# VALIDATION
# ============================================================================

def validate_dataset(records):
    failures = []

    for rec in records:

        Q = quotient_exact(
            rec.k,
            rec.ell,
        )

        reconstructed = canon(
            rec.C.subs(
                X,
                S * (S - 1),
            )
            +
            (2*S - 1)
            *
            rec.D.subs(
                X,
                S * (S - 1),
            )
        )

        if canon(
            Q - reconstructed
        ) != 0:
            failures.append(
                (
                    rec.k,
                    rec.ell,
                )
            )

    return failures


# ============================================================================
# INVOLUTION VALIDATION
# ============================================================================

def involution_validation(records):
    failures = []

    for rec in records:

        Q = quotient_exact(
            rec.k,
            rec.ell,
        )

        C_S = canon(
            rec.C.subs(
                X,
                S * (S - 1),
            )
        )

        D_S = canon(
            rec.D.subs(
                X,
                S * (S - 1),
            )
        )

        q_plus = canon(
            C_S + (2*S - 1) * D_S
        )

        q_minus = canon(
            C_S - (2*S - 1) * D_S
        )

        if canon(
            Q - q_plus
        ) != 0:
            failures.append(
                (rec.k, rec.ell, "plus")
            )

        reflected = canon(
            Q.subs(
                S,
                1 - S,
            )
        )

        if canon(
            reflected - q_minus
        ) != 0:
            failures.append(
                (rec.k, rec.ell, "minus")
            )

    return failures


# ============================================================================
# COMPLEXITY PROFILE
# ============================================================================

def complexity_profile(records):

    print("=" * 78)
    print("C/D DEGREE PROFILE")
    print("=" * 78)

    for rec in records:

        c_dx = degree_string(
            rec.C,
            X,
        )

        c_dn = degree_string(
            rec.C,
            N,
        )

        d_dx = degree_string(
            rec.D,
            X,
        )

        d_dn = degree_string(
            rec.D,
            N,
        )

        print(
            f"({rec.k},{rec.ell}) "
            f"m={rec.m:2d} "
            f"r={rec.r:2d} | "
            f"C: degX={c_dx:>3} "
            f"degN={c_dn:>3} "
            f"terms={term_count(rec.C):>3} | "
            f"D: degX={d_dx:>3} "
            f"degN={d_dn:>3} "
            f"terms={term_count(rec.D):>3}"
        )


# ============================================================================
# DEGREE PROFILE SUMMARY
# ============================================================================

def profile_summary(records):

    print("\n" + "=" * 78)
    print("PROFILE SUMMARY")
    print("=" * 78)

    by_k = {}

    for rec in records:
        by_k.setdefault(
            rec.k,
            []
        ).append(rec)

    for k in sorted(by_k):

        rows = sorted(
            by_k[k],
            key=lambda rec: rec.ell,
        )

        print(f"\nk={k}")

        print(
            "  ell | "
            "C_degX C_degN | "
            "D_degX D_degN"
        )

        for rec in rows:

            print(
                f"  {rec.ell:3d} | "
                f"{degree_string(rec.C, X):>6} "
                f"{degree_string(rec.C, N):>6} | "
                f"{degree_string(rec.D, X):>6} "
                f"{degree_string(rec.D, N):>6}"
            )


# ============================================================================
# COMMON FACTOR CHECKS
# ============================================================================

def common_factor_report(records):

    print("\n" + "=" * 78)
    print("C/D COMMON-FACTOR REPORT")
    print("=" * 78)

    for rec in records:

        if rec.D == 0:
            print(
                f"({rec.k},{rec.ell}) D=0"
            )
            continue

        c_num = sp.together(
            rec.C
        ).as_numer_denom()[0]

        d_num = sp.together(
            rec.D
        ).as_numer_denom()[0]

        gcd = canon(
            sp.gcd(
                sp.Poly(
                    c_num,
                    N,
                    X,
                ),
                sp.Poly(
                    d_num,
                    N,
                    X,
                ),
            ).as_expr()
        )

        print(
            f"({rec.k},{rec.ell}) "
            f"gcd(C,D)={gcd}"
        )


# ============================================================================
# BASIS REPRESENTATION
# ============================================================================

def coefficient_matrix(expr):
    """
    Return coefficients of polynomial in N,X.
    """

    P = sp.Poly(
        sp.expand(expr),
        N,
        X,
    )

    return {
        monomial: coeff
        for monomial, coeff
        in P.terms()
    }


# ============================================================================
# FINITE DIFFERENCE IN m
# ============================================================================

def finite_difference_sequence(records, field, fixed_k):

    seq = {}

    for rec in records:

        if rec.k != fixed_k:
            continue

        expr = (
            rec.C
            if field == "C"
            else rec.D
        )

        seq[rec.m] = expr

    return seq


def difference_operator(seq, order):

    current = seq

    for _ in range(order):

        keys = sorted(current)

        nxt = {}

        for a, b in zip(
            keys[:-1],
            keys[1:],
        ):

            if b != a + 1:
                continue

            nxt[a] = canon(
                current[b]
                - current[a]
            )

        current = nxt

    return current


def finite_difference_report(
    records,
    field,
    max_order=6,
):

    print("\n" + "=" * 78)
    print(
        f"FINITE-DIFFERENCE REPORT — {field}"
    )
    print("=" * 78)

    ks = sorted(
        {
            rec.k
            for rec in records
        }
    )

    for k in ks:

        seq = finite_difference_sequence(
            records,
            field,
            k,
        )

        print(
            f"\nk={k}"
        )

        for order in range(
            1,
            max_order + 1,
        ):

            diff = difference_operator(
                seq,
                order,
            )

            if not diff:
                break

            complexity = sorted(
                {
                    (
                        degree_int(expr, N),
                        degree_int(expr, X),
                    )
                    for expr in diff.values()
                }
            )

            print(
                f"  Δ^{order}: "
                f"rows={len(diff)} "
                f"degree_pairs={complexity}"
            )


# ============================================================================
# SIMPLE RECURRENCE SEARCH
# ============================================================================

def recurrence_matrix(
    sequence,
    order,
):
    """
    Build a simple scalar recurrence test:

        f_m - f_(m-1) ...
    """

    keys = sorted(sequence)

    equations = []

    for idx in range(
        order,
        len(keys),
    ):

        m = keys[idx]

        if any(
            keys[idx-j] != m-j
            for j in range(
                1,
                order + 1,
            )
        ):
            continue

        # Polynomial identity is tested numerically as symbolic
        # zero only for constant coefficients  — intentionally conservative.
        equations.append(
            (
                m,
                [
                    sequence[
                        keys[idx-j]
                    ]
                    for j in range(
                        0,
                        order + 1,
                    )
                ],
            )
        )

    return equations


def constant_recurrence_test(
    sequence,
    max_order=5,
):
    """
    Search for a constant-coefficient recurrence

        f_m = a1 f_(m-1) + ... + at f_(m-t).

    This deliberately tests only rational constants.
    """

    keys = sorted(sequence)

    for order in range(
        1,
        max_order + 1,
    ):

        if len(keys) <= order:
            continue

        symbols = sp.symbols(
            f"u0:{order}"
        )

        equations = []

        for idx in range(
            order,
            len(keys),
        ):

            m = keys[idx]

            if any(
                keys[idx-j] != m-j
                for j in range(
                    1,
                    order + 1,
                )
            ):
                continue

            lhs = sequence[m]

            rhs = sum(
                symbols[j-1]
                * sequence[
                    m-j
                ]
                for j in range(
                    1,
                    order + 1,
                )
            )

            diff = sp.Poly(
                sp.expand(
                    lhs-rhs
                ),
                N,
                X,
            )

            equations.extend(
                coeff
                for _, coeff
                in diff.terms()
            )

        if not equations:
            continue

        sol = sp.linsolve(
            equations,
            symbols,
        )

        if sol is sp.EmptySet:
            continue

        solutions = list(sol)

        if not solutions:
            continue

        candidate = solutions[0]

        if any(
            v.free_symbols
            for v in candidate
        ):
            continue

        # Verify.
        valid = True

        for idx in range(
            order,
            len(keys),
        ):

            m = keys[idx]

            if any(
                keys[idx-j] != m-j
                for j in range(
                    1,
                    order + 1,
                )
            ):
                continue

            predicted = sum(
                candidate[j-1]
                * sequence[m-j]
                for j in range(
                    1,
                    order + 1,
                )
            )

            if canon(
                sequence[m]-predicted
            ) != 0:
                valid = False
                break

        if valid:
            return order, candidate

    return None


# ============================================================================
# HOLDOUT STRUCTURAL TEST
# ============================================================================

def holdout_build(
    train_max_ell,
    test_ell,
):
    train = build_cd_dataset(
        train_max_ell
    )

    test = build_cd_dataset(
        test_ell
    )

    return train, test


def holdout_degree_report(
    train,
    test,
):

    print("\n" + "=" * 78)
    print("HOLDOUT DEGREE PROFILE")
    print("=" * 78)

    train_map = {
        (rec.k, rec.m): rec
        for rec in train
    }

    test_map = {
        (rec.k, rec.m): rec
        for rec in test
    }

    failures = []

    for key, rec_test in sorted(
        test_map.items()
    ):

        rec_train = train_map.get(
            key
        )

        if rec_train is None:
            continue

        fields = (
            ("C", rec_train.C, rec_test.C),
            ("D", rec_train.D, rec_test.D),
        )

        for field, train_expr, test_expr in fields:

            if degree_int(
                train_expr,
                X,
            ) != degree_int(
                test_expr,
                X,
            ):

                failures.append(
                    (
                        key,
                        field,
                        "degX",
                    )
                )

            if degree_int(
                train_expr,
                N,
            ) != degree_int(
                test_expr,
                N,
            ):

                failures.append(
                    (
                        key,
                        field,
                        "degN",
                    )
                )

    print(
        f"degree-profile failures = "
        f"{len(failures)}"
    )

    return failures


# ============================================================================
# SAMPLE EXPRESSIONS
# ============================================================================

def sample_report(records):

    wanted = {
        (1, 3),
        (1, 5),
        (1, 7),
        (3, 5),
        (3, 7),
        (5, 7),
    }

    print("\n" + "=" * 78)
    print("SAMPLE C/D DECOMPOSITIONS")
    print("=" * 78)

    for rec in records:

        if (
            rec.k,
            rec.ell,
        ) not in wanted:
            continue

        print(
            f"\n({rec.k},{rec.ell})"
        )
        print(
            f"  C(N,X) = {rec.C}"
        )
        print(
            f"  D(N,X) = {rec.D}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    MAX_ELL = 23

    print("=" * 78)
    print("KAPPA EXPERIMENT 120R-P")
    print("C/D RECURRENCE KERNEL")
    print("PRINT-SAFE DEGREE ANALYSIS")
    print("ROBUST QUADRATIC-BRANCH DECOMPOSITION")
    print("MODULO S^2-S-X")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1
    # ----------------------------------------------------------------------

    print(
        "\n1. BUILDING EXACT C/D DATASET"
    )

    records = build_cd_dataset(
        MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    # ----------------------------------------------------------------------
    # 2
    # ----------------------------------------------------------------------

    print(
        "\n2. FULL C/D RECONSTRUCTION"
    )

    reconstruction_failures = validate_dataset(
        records
    )

    print(
        f"reconstruction failures = "
        f"{len(reconstruction_failures)}/"
        f"{len(records)}"
    )

    if reconstruction_failures:
        for item in reconstruction_failures[:20]:
            print(
                f"  failure = {item}"
            )

        raise ArithmeticError(
            "C/D reconstruction validation failed."
        )

    print("STATUS = PASS")

    # ----------------------------------------------------------------------
    # 3
    # ----------------------------------------------------------------------

    print(
        "\n3. INVOLUTION VALIDATION"
    )

    involution_failures = involution_validation(
        records
    )

    print(
        f"involution/reconstruction failures = "
        f"{len(involution_failures)}"
    )

    if involution_failures:
        for item in involution_failures[:20]:
            print(
                f"  failure = {item}"
            )

        raise ArithmeticError(
            "Involution validation failed."
        )

    print("STATUS = PASS")

    # ----------------------------------------------------------------------
    # 4
    # ----------------------------------------------------------------------

    print(
        "\n4. DEGREE PROFILE"
    )

    complexity_profile(
        records
    )

    # ----------------------------------------------------------------------
    # 5
    # ----------------------------------------------------------------------

    profile_summary(
        records
    )

    # ----------------------------------------------------------------------
    # 6
    # ----------------------------------------------------------------------

    sample_report(
        records
    )

    # ----------------------------------------------------------------------
    # 7
    # ----------------------------------------------------------------------

    common_factor_report(
        records
    )

    # ----------------------------------------------------------------------
    # 8
    # ----------------------------------------------------------------------

    print(
        "\n8. FINITE DIFFERENCE TESTS"
    )

    finite_difference_report(
        records,
        "C",
        max_order=5,
    )

    finite_difference_report(
        records,
        "D",
        max_order=5,
    )

    # ----------------------------------------------------------------------
    # 9
    # ----------------------------------------------------------------------

    print(
        "\n9. CONSTANT-COEFFICIENT RECURRENCE SEARCH"
    )

    for field in ("C", "D"):

        print(
            f"\nFIELD = {field}"
        )

        for k in sorted(
            {
                rec.k
                for rec in records
            }
        ):

            seq = finite_difference_sequence(
                records,
                field,
                k,
            )

            result = constant_recurrence_test(
                seq,
                max_order=5,
            )

            if result is None:

                print(
                    f"  k={k}: "
                    f"no constant recurrence"
                )

            else:

                order, coeffs = result

                print(
                    f"  k={k}: "
                    f"order={order} "
                    f"coeffs={coeffs}"
                )

    # ----------------------------------------------------------------------
    # 10
    # ----------------------------------------------------------------------

    print(
        "\n10. HOLDOUT STRUCTURE"
    )

    train, test = holdout_build(
        train_max_ell=19,
        test_ell=21,
    )

    holdout_failures = holdout_degree_report(
        train,
        test,
    )

    if holdout_failures:
        print(
            "STATUS = DEGREE PATTERN CHANGED"
        )
    else:
        print(
            "STATUS = DEGREE PROFILE STABLE"
        )

    # ----------------------------------------------------------------------
    # 11
    # ----------------------------------------------------------------------

    print(
        "\n11. DIRECT ALGEBRAIC BRIDGE"
    )

    bridge_failures = 0

    for rec in records:

        Q = quotient_exact(
            rec.k,
            rec.ell,
        )

        C_S = rec.C.subs(
            X,
            S * (S - 1),
        )

        D_S = rec.D.subs(
            X,
            S * (S - 1),
        )

        lhs = canon(
            Q - C_S
        )

        rhs = canon(
            (2*S - 1)*D_S
        )

        if canon(
            lhs-rhs
        ) != 0:
            bridge_failures += 1

        discriminant = canon(
            (2*S-1)**2
            -
            (4*S*(S-1)+1)
        )

        if discriminant != 0:
            bridge_failures += 1

    print(
        f"bridge failures = "
        f"{bridge_failures}"
    )

    if bridge_failures:
        raise ArithmeticError(
            "Quadratic branch bridge failed."
        )

    print("STATUS = PASS")

    # ----------------------------------------------------------------------
    # 12
    # ----------------------------------------------------------------------

    print(
        "\n12. FINAL DIAGNOSTIC"
    )
    print("=" * 78)

    print(
        """
The previous execution stopped for a formatting bug only:

    SymPy returns -oo for the zero polynomial,
    and NegativeInfinity does not support numeric formatting
    such as >3.

This patch converts symbolic degree values to strings before
formatting. The underlying decomposition is unchanged.

The mathematical object remains

    Q(N,S) = C(N,X) + (2S-1)D(N,X)

with

    X = S(S-1)

and

    (2S-1)^2 = 4X+1.

The experiment now investigates whether C and D themselves
display finite-difference or recurrence structure in the two
indices.

The important results are:

    * exact C/D reconstruction;
    * involution validation;
    * degree stability;
    * finite-difference structure;
    * constant-coefficient recurrence candidates;
    * ell holdout degree stability.

No numerical fitting is used.
No factor-pair search is used.
No N-only factoring claim is made.
"""
    )

    print(
        "\n" + "=" * 78
    )
    print(
        "EXPERIMENT 120R-P COMPLETE"
    )
    print(
        "=" * 78
    )


if __name__ == "__main__":
    main()