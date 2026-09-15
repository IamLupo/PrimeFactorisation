#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 129
CORRECT POST-BOUNDARY SINGLE-COEFFICIENT LEADING LAW
DIRECT INDEX:
    j = ell - k - d - 1
FINITE BINOMIAL KERNEL -> ONE NEWTON COEFFICIENT

Purpose
-------
Experiment 128 established exactly where the Experiment-126 coefficients live
inside the Newton tensor:

    j = ell - k - d - 1.

This experiment uses that index directly.

For each admissible (k, ell, d):

    r = k + d
    j = ell - 1 - r

we extract ONLY c_j(N), then verify:

    1. exact Newton coefficient,
    2. actual N-degree,
    3. actual leading coefficient,
    4. Experiment-126 leading law,
    5. leave-one-ell-out consistency,
    6. forward ell holdout.

No full Newton tensor is retained after extraction.
No C/D tensor.
No polynomial interpolation of the leading coefficient.
No factor-pair search.
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")


# ============================================================================
# Exact kernel and quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(PF, PD)

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric conversion
# ============================================================================

def symmetric_to_NS(expr):
    """
    q = S-p, reduced modulo p^2-Sp+N.
    """

    expr_sub = sp.expand(expr.subs(q, S - p))

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        expr_sub,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    p_coeff = sp.expand(rem.coeff(p, 1))

    if p_coeff != 0:
        raise ArithmeticError(
            "Symmetric reduction retained p-dependence: "
            + str(sp.factor(p_coeff))
        )

    return sp.expand(rem.coeff(p, 0))


# ============================================================================
# Newton basis
# ============================================================================

def newton_basis_up_to(max_j):
    """
    P0 = 2
    P1 = S
    Pj = S*P(j-1) - N*P(j-2)
    """

    P = [sp.Integer(2)]

    if max_j >= 1:
        P.append(S)

    for j in range(2, max_j + 1):
        P.append(
            sp.expand(
                S * P[j - 1]
                - N * P[j - 2]
            )
        )

    return P


# ============================================================================
# Single-target extraction
# ============================================================================

def extract_target_coefficient(Q_NS, ell, target_j):
    """
    Extract c_target_j in

        Q = sum_j c_j(N) P_j.

    Important:
        target_j is supplied explicitly.

    We descend only from ell-1 to target_j.
    """

    if target_j < 1:
        raise ValueError(
            "target_j must be >= 1"
        )

    if target_j >= ell:
        raise ValueError(
            f"target_j={target_j} invalid for ell={ell}"
        )

    P = newton_basis_up_to(ell - 1)

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    target = None

    for j in range(ell - 1, target_j - 1, -1):

        c = sp.expand(
            work.coeff_monomial(S**j)
        )

        if j == target_j:
            target = c
            break

        if c != 0:
            subtract = sp.Poly(
                sp.expand(c * P[j]),
                S,
                domain=sp.QQ.frac_field(N),
            )

            work = work - subtract

    if target is None:
        raise ArithmeticError(
            f"Could not extract target j={target_j}"
        )

    return sp.expand(target)


# ============================================================================
# Experiment-126 leading law
# ============================================================================

def expected_leading(k, ell, d):
    """
    Closed law established empirically by Experiment 126.

    a=(k-1)/2
    b=(ell-1)/2

    d odd:
        2*(-1)^((d-1)/2)*(a+b+1)

    d even:
        (-1)^(d/2)*(a+b+1)*(2a-2b+d+1)
    """

    a = sp.Rational(k - 1, 2)
    b = sp.Rational(ell - 1, 2)

    if d % 2 == 1:
        return sp.expand(
            2 * (-1)**((d - 1) // 2)
            * (a + b + 1)
        )

    return sp.expand(
        (-1)**(d // 2)
        * (a + b + 1)
        * (2*a - 2*b + d + 1)
    )


# ============================================================================
# Experiment-126 target degree
# ============================================================================

def expected_degree(k, ell, d):
    """
    From the verified post-boundary examples:

        r = k+d
        j = ell-1-r

    and the observed N-degree is

        floor((ell-1)/2) - r
        +
        parity-dependent offset.

    Rather than impose a guessed degree law, this experiment verifies only
    the leading coefficient law. Degree is reported diagnostically.

    Returns None intentionally.
    """
    return None


# ============================================================================
# One exact row
# ============================================================================

def row_record(k, ell, d):
    r = k + d
    j = ell - 1 - r

    if j < 1:
        return None

    Qpq = quotient_Q(k, ell)
    QNS = symmetric_to_NS(Qpq)

    c_j = extract_target_coefficient(
        QNS,
        ell,
        j,
    )

    PN = sp.Poly(
        c_j,
        N,
        domain=sp.QQ,
    )

    if PN.is_zero:
        degree = None
        lead = sp.Integer(0)
    else:
        degree = int(PN.degree())
        lead = sp.expand(PN.LC())

    expected = sp.expand(
        expected_leading(k, ell, d)
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "r": r,
        "j": j,
        "coefficient": c_j,
        "degree": degree,
        "lead": lead,
        "expected": expected,
        "match": (
            sp.simplify(
                lead - expected
            ) == 0
        ),
    }


# ============================================================================
# Admissible rows
# ============================================================================

def admissible_rows(max_ell):
    rows = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):

            max_d = ell - k - 2

            for d in range(1, max_d + 1):

                r = k + d
                j = ell - 1 - r

                if j >= 1:
                    rows.append(
                        (k, ell, d)
                    )

    return rows


# ============================================================================
# Training verification
# ============================================================================

def train_test(max_ell):
    rows = admissible_rows(max_ell)

    failures = []

    print("=" * 78)
    print("1. TRAINING CERTIFICATE")
    print("=" * 78)

    for k, ell, d in rows:

        rec = row_record(k, ell, d)

        if not rec["match"]:
            failures.append(rec)

    print(
        f"rows checked = {len(rows)}"
    )
    print(
        f"leading-law failures = {len(failures)}"
    )

    if failures:
        print()
        print("FIRST FAILURES")
        print("-" * 78)

        for rec in failures[:20]:
            print(
                f"({rec['k']},{rec['ell']}) "
                f"d={rec['d']} "
                f"r={rec['r']} "
                f"j={rec['j']}"
            )
            print(
                f"  c_j(N) = "
                f"{sp.factor(rec['coefficient'])}"
            )
            print(
                f"  degree  = {rec['degree']}"
            )
            print(
                f"  lead    = {rec['lead']}"
            )
            print(
                f"  expected= {rec['expected']}"
            )
            print()

    else:
        print("STATUS = PASS")

    return rows, failures


# ============================================================================
# Leave-one-ell-out
# ============================================================================

def leave_one_ell_out(max_ell):
    """
    Since the law is closed-form rather than fitted, LOO is used as a
    consistency sweep: each ell is removed from the reporting set and
    independently recomputed.

    The important test is simply that every held-out ell passes the exact
    leading-law identity.
    """

    print("=" * 78)
    print("2. LEAVE-ONE-ELL-OUT CONSISTENCY")
    print("=" * 78)

    all_rows = admissible_rows(max_ell)

    ells = sorted(
        {ell for _, ell, _ in all_rows}
    )

    total_failures = 0

    for held_out in ells:

        local_failures = 0
        count = 0

        for k, ell, d in all_rows:

            if ell != held_out:
                continue

            count += 1

            rec = row_record(k, ell, d)

            if not rec["match"]:
                local_failures += 1

        total_failures += local_failures

        print(
            f"ell={held_out:2d}: "
            f"rows={count:2d} "
            f"failures={local_failures}"
        )

    print(
        f"LOO failures = {total_failures}"
    )

    return total_failures


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ells):
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    failures = 0
    total = 0

    for ell in ells:

        local_total = 0
        local_fail = 0

        for k in range(1, ell, 2):

            max_d = ell - k - 2

            for d in range(1, max_d + 1):

                rec = row_record(k, ell, d)

                if rec is None:
                    continue

                local_total += 1
                total += 1

                if not rec["match"]:
                    local_fail += 1
                    failures += 1

        print(
            f"ell={ell:2d}: "
            f"rows={local_total:3d} "
            f"failures={local_fail:3d}"
        )

    print()
    print(
        f"forward rows = {total}"
    )
    print(
        f"forward failures = {failures}"
    )

    return failures


# ============================================================================
# Print reference-style examples
# ============================================================================

def examples():
    print("=" * 78)
    print("4. INDEXED COEFFICIENT EXAMPLES")
    print("=" * 78)

    examples_list = [
        (1, 5, 1),
        (3, 7, 1),
        (5, 9, 1),
        (1, 7, 2),
        (3, 9, 2),
        (5, 11, 2),
        (1, 9, 3),
        (3, 11, 4),
        (5, 13, 5),
        (1, 13, 6),
        (3, 17, 8),
        (1, 23, 10),
    ]

    for k, ell, d in examples_list:

        rec = row_record(k, ell, d)

        print(
            f"({k},{ell}) "
            f"d={d} "
            f"r={rec['r']} "
            f"j={rec['j']}"
        )
        print(
            f"  c_j(N) = "
            f"{sp.factor(rec['coefficient'])}"
        )
        print(
            f"  degree  = {rec['degree']}"
        )
        print(
            f"  lead    = {rec['lead']}"
        )
        print(
            f"  expected= {rec['expected']}"
        )
        print(
            f"  match   = {rec['match']}"
        )
        print()


# ============================================================================
# d-summary
# ============================================================================

def d_summary(max_ell):
    print("=" * 78)
    print("5. SUMMARY BY d")
    print("=" * 78)

    rows = admissible_rows(max_ell)

    for d in sorted(
        {d for _, _, d in rows}
    ):

        subset = [
            row for row in rows
            if row[2] == d
        ]

        fail = 0

        for k, ell, _ in subset:
            rec = row_record(
                k,
                ell,
                d,
            )

            if not rec["match"]:
                fail += 1

        print(
            f"d={d:2d}: "
            f"rows={len(subset):3d} "
            f"failures={fail:3d}"
        )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 129")
    print("CORRECT POST-BOUNDARY SINGLE-COEFFICIENT LEADING LAW")
    print("INDEX j = ell - k - d - 1")
    print("FINITE BINOMIAL KERNEL -> ONE NEWTON COEFFICIENT")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 15

    # Use a small forward holdout because exact symbolic expansion is costly.
    FORWARD_ELL = [17]

    # ------------------------------------------------------------------
    # 0. Quotient sanity
    # ------------------------------------------------------------------

    print("=" * 78)
    print("0. QUOTIENT SANITY")
    print("=" * 78)

    sanity_cases = [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]

    quotient_failures = 0

    for k, ell in sanity_cases:

        try:
            quotient_Q(k, ell)
            print(
                f"({k},{ell}) quotient = PASS"
            )
        except Exception as exc:
            quotient_failures += 1
            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = "
        f"{quotient_failures}/{len(sanity_cases)}"
    )

    if quotient_failures:
        raise ArithmeticError(
            "Quotient sanity failed."
        )

    print()

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    _, failures = train_test(
        TRAIN_MAX_ELL
    )

    if failures:
        raise ArithmeticError(
            "Correct-index leading law failed."
        )

    print()

    # ------------------------------------------------------------------
    # Examples
    # ------------------------------------------------------------------

    examples()
    print()

    # ------------------------------------------------------------------
    # LOO
    # ------------------------------------------------------------------

    loo_failures = leave_one_ell_out(
        TRAIN_MAX_ELL
    )

    if loo_failures:
        raise ArithmeticError(
            "LOO consistency failed."
        )

    print()

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    forward_failures = forward_holdout(
        FORWARD_ELL
    )

    if forward_failures:
        raise ArithmeticError(
            "Forward ell holdout failed."
        )

    print()

    # ------------------------------------------------------------------
    # d-summary
    # ------------------------------------------------------------------

    d_summary(
        TRAIN_MAX_ELL
    )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "Experiment 128 established the exact index:"
    )
    print(
        "    j = ell - k - d - 1"
    )
    print()
    print(
        "This experiment uses that index directly and does not"
    )
    print(
        "assume j=d+1."
    )
    print()
    print(
        "The next mathematical question after a PASS is not another"
    )
    print(
        "index search. It is the symbolic derivation of the leading"
    )
    print(
        "coefficient of this exact Newton coefficient."
    )
    print()
    print(
        "The coefficient being tested is:"
    )
    print(
        "    c_{ell-k-d-1}(N)."
    )
    print()
    print(
        "This is the narrowest experiment needed to reconnect"
    )
    print(
        "Experiment 126 with the independently reconstructed"
    )
    print(
        "Newton tensor."
    )
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc),
        )
        sys.exit(1)

