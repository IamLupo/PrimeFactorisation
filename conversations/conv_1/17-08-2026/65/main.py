#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 130
TOP-N EDGE CERTIFICATE FROM THE CORRECT NEWTON INDEX

TARGET:
    j = ell - k - d - 1

CONJECTURED N-DEGREE:
    D = k + floor((d-1)/2)

CONJECTURED LEADING COEFFICIENT:

    d odd:
        L = (-1)^((d-1)/2) * (k + ell)

    d even:
        L = (-1)^(d/2)
            * (k + ell)/2
            * (k - ell + d + 1)

GOAL
----
Use the exact finite binomial kernel and exact Newton recurrence, but:

  * build each (k,ell) quotient only once;
  * convert to (S,N) only once;
  * extract all required post-boundary coefficients for that pair in one
    descending Newton pass;
  * verify both the degree law and the leading-coefficient law;
  * perform leave-one-ell-out and forward-ell checks.

This is the first experiment after 128/129 that tests the COMPLETE
top-degree formula rather than only the leading coefficient.

It deliberately does NOT fit a polynomial in (k,ell,d).

NO FULL C/D TENSOR
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
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
# Exact finite binomial kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact synthetic quotient
# ============================================================================

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

    remainder = sp.expand(R.as_expr())

    if remainder != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(remainder)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric p,q -> S,N
# ============================================================================

def symmetric_to_NS(expr):
    """
    Substitute q=S-p and reduce modulo

        p^2 - S*p + N.

    A symmetric expression must leave no p term.
    """

    expr_sub = sp.expand(
        expr.subs(q, S - p)
    )

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

    p_coeff = sp.expand(
        rem.coeff(p, 1)
    )

    if p_coeff != 0:
        raise ArithmeticError(
            "Symmetric conversion retained p-dependence: "
            + str(sp.factor(p_coeff))
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Newton basis
# ============================================================================

def build_newton_basis(max_j):
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
# Admissible post-boundary rows
# ============================================================================

def target_rows_for_pair(k, ell):
    rows = []

    for d in range(1, ell - k - 1):
        r = k + d
        j = ell - 1 - r

        if j < 1:
            continue

        rows.append(
            (d, r, j)
        )

    return rows


# ============================================================================
# Conjectured degree law
# ============================================================================

def expected_degree(k, d):
    return k + (d - 1) // 2


# ============================================================================
# Conjectured leading coefficient law
# ============================================================================

def expected_lead(k, ell, d):

    if d % 2 == 1:
        return sp.expand(
            (-1) ** ((d - 1) // 2)
            * (k + ell)
        )

    return sp.expand(
        (-1) ** (d // 2)
        * sp.Rational(k + ell, 2)
        * (k - ell + d + 1)
    )


# ============================================================================
# Extract all required coefficients for one (k,ell)
# ============================================================================

def extract_pair_targets(Q_NS, ell, targets):
    """
    One descending Newton pass for a single (k,ell).

    We do not retain the full tensor. We only save coefficients whose
    indices are requested by targets.
    """

    if not targets:
        return {}

    target_by_j = {
        j: d
        for d, r, j in targets
    }

    min_j = min(
        target_by_j
    )

    P = build_newton_basis(
        ell - 1
    )

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    found = {}

    for j in range(
        ell - 1,
        min_j - 1,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(S**j)
        )

        if j in target_by_j:
            found[j] = sp.expand(c)

        if c != 0:
            subtract = sp.Poly(
                sp.expand(c * P[j]),
                S,
                domain=sp.QQ.frac_field(N),
            )

            work = work - subtract

    missing = [
        j for j in target_by_j
        if j not in found
    ]

    if missing:
        raise ArithmeticError(
            f"Failed to extract Newton indices "
            f"{missing} for ell={ell}"
        )

    return found


# ============================================================================
# Analyze one coefficient
# ============================================================================

def analyze_coefficient(
    k,
    ell,
    d,
    r,
    j,
    coefficient,
):
    poly = sp.Poly(
        sp.expand(coefficient),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        degree = None
        lead = sp.Integer(0)
    else:
        degree = int(poly.degree())
        lead = sp.expand(poly.LC())

    expected_D = expected_degree(
        k,
        d
    )

    expected_L = expected_lead(
        k,
        ell,
        d
    )

    degree_ok = (
        degree == expected_D
    )

    lead_ok = (
        sp.expand(
            lead - expected_L
        ) == 0
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "r": r,
        "j": j,
        "coefficient": sp.expand(coefficient),
        "degree": degree,
        "expected_degree": expected_D,
        "lead": lead,
        "expected_lead": expected_L,
        "degree_ok": degree_ok,
        "lead_ok": lead_ok,
        "ok": (
            degree_ok
            and lead_ok
        ),
    }


# ============================================================================
# Build exact dataset once per (k,ell)
# ============================================================================

def build_dataset(max_ell):

    records = []
    pair_count = 0

    print("=" * 78)
    print("BUILDING CACHED EXACT POST-BOUNDARY DATASET")
    print("=" * 78)

    for ell in range(
        3,
        max_ell + 1,
        2
    ):

        for k in range(
            1,
            ell,
            2
        ):

            targets = target_rows_for_pair(
                k,
                ell
            )

            if not targets:
                continue

            pair_count += 1

            Qpq = quotient_Q(
                k,
                ell
            )

            QNS = symmetric_to_NS(
                Qpq
            )

            coefficients = extract_pair_targets(
                QNS,
                ell,
                targets
            )

            for d, r, j in targets:

                rec = analyze_coefficient(
                    k,
                    ell,
                    d,
                    r,
                    j,
                    coefficients[j]
                )

                records.append(
                    rec
                )

    print(
        f"pairs processed = {pair_count}"
    )
    print(
        f"rows            = {len(records)}"
    )
    print()

    return records


# ============================================================================
# Validation
# ============================================================================

def validate_records(records):
    degree_failures = []
    lead_failures = []

    for rec in records:

        if not rec["degree_ok"]:
            degree_failures.append(
                rec
            )

        if not rec["lead_ok"]:
            lead_failures.append(
                rec
            )

    print("=" * 78)
    print("EXACT TOP-N EDGE CERTIFICATE")
    print("=" * 78)

    print(
        f"rows                  = {len(records)}"
    )
    print(
        f"degree-law failures   = {len(degree_failures)}"
    )
    print(
        f"leading-law failures  = {len(lead_failures)}"
    )

    if degree_failures:
        print()
        print("FIRST DEGREE FAILURES")
        print("-" * 78)

        for rec in degree_failures[:15]:
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
                f"  degree={rec['degree']} "
                f"expected={rec['expected_degree']}"
            )
            print()

    if lead_failures:
        print()
        print("FIRST LEADING-COEFFICIENT FAILURES")
        print("-" * 78)

        for rec in lead_failures[:15]:
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
                f"  degree={rec['degree']} "
                f"expected={rec['expected_degree']}"
            )
            print(
                f"  lead={rec['lead']} "
                f"expected={rec['expected_lead']}"
            )
            print()

    return (
        degree_failures,
        lead_failures,
    )


# ============================================================================
# Examples
# ============================================================================

def print_examples(records):

    wanted = {
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
    }

    print("=" * 78)
    print("CERTIFICATE EXAMPLES")
    print("=" * 78)

    for rec in records:

        key = (
            rec["k"],
            rec["ell"],
            rec["d"],
        )

        if key not in wanted:
            continue

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
            f"  degree  = "
            f"{rec['degree']} "
            f"(expected {rec['expected_degree']})"
        )

        print(
            f"  lead    = "
            f"{rec['lead']} "
            f"(expected {rec['expected_lead']})"
        )

        print(
            f"  status  = "
            f"{rec['ok']}"
        )

        print()


# ============================================================================
# Leave-one-ell-out
# ============================================================================

def loo_summary(records):

    print("=" * 78)
    print("LEAVE-ONE-ELL-OUT")
    print("=" * 78)

    ells = sorted(
        {
            rec["ell"]
            for rec in records
        }
    )

    total_failures = 0

    for held_out in ells:

        subset = [
            rec
            for rec in records
            if rec["ell"] == held_out
        ]

        degree_failures = sum(
            not rec["degree_ok"]
            for rec in subset
        )

        lead_failures = sum(
            not rec["lead_ok"]
            for rec in subset
        )

        local = (
            degree_failures
            + lead_failures
        )

        total_failures += local

        print(
            f"ell={held_out:2d}: "
            f"rows={len(subset):3d} "
            f"degree_fail={degree_failures:2d} "
            f"lead_fail={lead_failures:2d}"
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
    print("FORWARD ELL HOLDOUT")
    print("=" * 78)

    all_records = []

    for ell in ells:

        for k in range(
            1,
            ell,
            2
        ):

            targets = target_rows_for_pair(
                k,
                ell
            )

            if not targets:
                continue

            Qpq = quotient_Q(
                k,
                ell
            )

            QNS = symmetric_to_NS(
                Qpq
            )

            coefficients = extract_pair_targets(
                QNS,
                ell,
                targets
            )

            for d, r, j in targets:

                all_records.append(
                    analyze_coefficient(
                        k,
                        ell,
                        d,
                        r,
                        j,
                        coefficients[j]
                    )
                )

    degree_failures = sum(
        not rec["degree_ok"]
        for rec in all_records
    )

    lead_failures = sum(
        not rec["lead_ok"]
        for rec in all_records
    )

    print(
        f"forward rows = {len(all_records)}"
    )
    print(
        f"degree failures = {degree_failures}"
    )
    print(
        f"leading failures = {lead_failures}"
    )

    for ell in ells:

        local = [
            rec
            for rec in all_records
            if rec["ell"] == ell
        ]

        print(
            f"ell={ell:2d}: "
            f"rows={len(local):3d} "
            f"degree_fail={sum(not r['degree_ok'] for r in local):2d} "
            f"lead_fail={sum(not r['lead_ok'] for r in local):2d}"
        )

    return (
        degree_failures,
        lead_failures,
        all_records,
    )


# ============================================================================
# Equivalent closed forms
# ============================================================================

def verify_equivalent_formula(records):

    print("=" * 78)
    print("PARITY-FORMULA EQUIVALENCE CHECK")
    print("=" * 78)

    failures = 0

    for rec in records:

        k = rec["k"]
        ell = rec["ell"]
        d = rec["d"]

        a = sp.Rational(
            k - 1,
            2
        )

        b = sp.Rational(
            ell - 1,
            2
        )

        if d % 2 == 1:

            alt = sp.expand(
                2
                * (-1) ** ((d - 1) // 2)
                * (a + b + 1)
            )

        else:

            alt = sp.expand(
                (-1) ** (d // 2)
                * (a + b + 1)
                * (
                    2*a
                    - 2*b
                    + d
                    + 1
                )
            )

        if sp.expand(
            alt - rec["expected_lead"]
        ) != 0:
            failures += 1

    print(
        f"parity-equivalence failures = "
        f"{failures}/{len(records)}"
    )

    return failures


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 130")
    print("TOP-N EDGE CERTIFICATE FROM THE CORRECT NEWTON INDEX")
    print("j = ell-k-d-1")
    print("D = k + floor((d-1)/2)")
    print("EXACT FINITE BINOMIAL KERNEL")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 17
    FORWARD_ELL = [19]

    # --------------------------------------------------------------
    # 0. Exact sanity
    # --------------------------------------------------------------

    print("=" * 78)
    print("0. BASIC QUOTIENT SANITY")
    print("=" * 78)

    sanity = [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]

    sanity_failures = 0

    for k, ell in sanity:
        try:
            quotient_Q(k, ell)
            print(
                f"({k},{ell}) quotient = PASS"
            )
        except Exception as exc:
            sanity_failures += 1
            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = "
        f"{sanity_failures}/{len(sanity)}"
    )

    if sanity_failures:
        raise ArithmeticError(
            "Basic quotient sanity failed."
        )

    print()

    # --------------------------------------------------------------
    # 1. Training
    # --------------------------------------------------------------

    records = build_dataset(
        TRAIN_MAX_ELL
    )

    degree_failures, lead_failures = validate_records(
        records
    )

    print()

    # --------------------------------------------------------------
    # 2. Examples
    # --------------------------------------------------------------

    print_examples(
        records
    )

    # --------------------------------------------------------------
    # 3. Formula equivalence
    # --------------------------------------------------------------

    print()
    equivalent_failures = verify_equivalent_formula(
        records
    )

    # --------------------------------------------------------------
    # 4. LOO
    # --------------------------------------------------------------

    print()
    loo_failures = loo_summary(
        records
    )

    # --------------------------------------------------------------
    # 5. Forward
    # --------------------------------------------------------------

    print()
    forward_degree_failures, forward_lead_failures, forward_records = (
        forward_holdout(
            FORWARD_ELL
        )
    )

    # --------------------------------------------------------------
    # 6. Final
    # --------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    total_failures = (
        degree_failures
        + lead_failures
        + equivalent_failures
        + loo_failures
        + forward_degree_failures
        + forward_lead_failures
    )

    if total_failures == 0:
        print(
            "STATUS = PASS"
        )
        print()
        print(
            "The exact coefficient"
        )
        print(
            "    c_{ell-k-d-1}(N)"
        )
        print(
            "satisfies the tested degree law"
        )
        print(
            "    deg_N = k + floor((d-1)/2)"
        )
        print(
            "and the tested leading-coefficient law"
        )
        print(
            "on the training grid and forward holdout."
        )
        print()
        print(
            "The next step is now a symbolic extraction"
        )
        print(
            "of the highest-N term from the finite binomial"
        )
        print(
            "kernel itself."
        )
    else:
        print(
            "STATUS = FAIL"
        )
        print(
            f"total diagnostic failures = "
            f"{total_failures}"
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
            str(exc)
        )
        sys.exit(1)

