# ============================================================================
# EXPERIMENT 202
# EXACT pq-KERNEL NORMALIZATION BY PARITY
# ============================================================================
#
# Goal:
#
#   Experiment 201 established experimentally that
#
#       p+q+1 | F(k,l;p,q)
#
#   on the tested odd-k grid exactly when k and l have the same parity.
#
#   Therefore the factorization
#
#       F = (p+q+1) Q(pq,p+q)
#
#   cannot be the universal normalization.
#
# This experiment determines the exact quotient/remainder structure in the
# two parity branches without guessing affine formulas.
#
# It also fixes two implementation problems from Experiment 201:
#
#   * factor detection is done by exact polynomial division;
#   * p <-> q symmetry uses simultaneous replacement.
#
# No previous experiment output is read.
# No unrestricted candidate search is performed.
#
# ============================================================================

import sympy as sp


# ----------------------------------------------------------------------------
# SYMBOLS
# ----------------------------------------------------------------------------

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ----------------------------------------------------------------------------
# EXACT pq KERNEL
# ----------------------------------------------------------------------------

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ----------------------------------------------------------------------------
# EXACT DIVISION BY p+q+1
# ----------------------------------------------------------------------------

def divide_by_linear(k, ell):
    F = sp.Poly(exact_F(k, ell), p, q)
    D = sp.Poly(p + q + 1, p, q)

    Q, R = sp.div(F, D)

    return sp.expand(Q.as_expr()), sp.expand(R.as_expr())


# ----------------------------------------------------------------------------
# TRUE SIMULTANEOUS SWAP p <-> q
# ----------------------------------------------------------------------------

def swap_pq(expr):
    return sp.expand(
        expr.xreplace({
            p: q,
            q: p,
        })
    )


def is_symmetric_pq(expr):
    return sp.expand(expr - swap_pq(expr)) == 0


# ----------------------------------------------------------------------------
# ROOT TEST
# ----------------------------------------------------------------------------

def root_test(k, ell):
    return sp.expand(
        exact_F(k, ell).subs(p, -q - 1)
    )


# ----------------------------------------------------------------------------
# EXACT (1+q)-POWER EXTRACTION
# ----------------------------------------------------------------------------

def q_factor_power(expr):
    """
    Return the exact multiplicity of (q+1) dividing expr.

    This never creates rational functions.
    """
    poly = sp.Poly(sp.expand(expr), q)

    if poly.is_zero:
        return None, sp.Integer(0)

    divisor = sp.Poly(q + 1, q)
    current = poly
    multiplicity = 0

    while True:
        quotient, remainder = sp.div(current, divisor)

        if remainder.is_zero:
            multiplicity += 1
            current = quotient
        else:
            break

    return multiplicity, sp.expand(current.as_expr())


# ----------------------------------------------------------------------------
# POLYNOMIAL CONTENT
# ----------------------------------------------------------------------------

def polynomial_content(expr):
    """
    Primitive integer content of a univariate polynomial in q.
    """
    poly = sp.Poly(sp.expand(expr), q)

    if poly.is_zero:
        return 0

    return int(sp.gcd_list(poly.all_coeffs()))


# ----------------------------------------------------------------------------
# PARITY-DEPENDENT BASIC NORMALIZATION
# ----------------------------------------------------------------------------

def normalization_data(k, ell):
    F = exact_F(k, ell)

    Q, R = divide_by_linear(k, ell)

    if R == 0:
        branch = "DIVISIBLE_BY_(p+q+1)"
    else:
        branch = "NONDIVISIBLE"

    return {
        "F": F,
        "Q": Q,
        "R": R,
        "branch": branch,
    }


# ----------------------------------------------------------------------------
# DIRECT CASE REPORT
# ----------------------------------------------------------------------------

def report_case(k, ell):
    data = normalization_data(k, ell)

    F = data["F"]
    Q = data["Q"]
    R = data["R"]

    print(
        f"k={k:2d} ell={ell:2d} "
        f"same_parity={'YES' if (k-ell) % 2 == 0 else 'NO'}"
    )

    print(
        f"  divisible by p+q+1 = "
        f"{'YES' if R == 0 else 'NO'}"
    )

    if R == 0:

        reconstructed = sp.expand(
            (p + q + 1) * Q
        )

        print(
            f"  quotient reconstruction = "
            f"{'PASS' if sp.expand(F - reconstructed) == 0 else 'FAIL'}"
        )

        print(
            f"  quotient p<->q symmetric = "
            f"{'YES' if is_symmetric_pq(Q) else 'NO'}"
        )

        print(
            f"  quotient degree in p = "
            f"{sp.Poly(Q, p, q).degree(p)}"
        )

        print(
            f"  quotient total degree = "
            f"{sp.Poly(Q, p, q).total_degree()}"
        )

    else:

        print(
            f"  remainder root-test = {root_test(k, ell)}"
        )

        mult, residual = q_factor_power(R)

        if mult is None:
            print("  remainder q+1 multiplicity = ALL")
        else:
            print(
                f"  remainder (q+1)^a multiplicity = {mult}"
            )

            print(
                f"  residual content = "
                f"{polynomial_content(residual)}"
            )

            print(
                f"  residual = {residual}"
            )


# ----------------------------------------------------------------------------
# PARITY CLASS SUMMARY
# ----------------------------------------------------------------------------

def parity_summary():

    print("=" * 78)
    print("1. PARITY-BRANCH SUMMARY")
    print("=" * 78)

    summary = {
        "same": {"total": 0, "divisible": 0, "nondivisible": 0},
        "different": {"total": 0, "divisible": 0, "nondivisible": 0},
    }

    for k in range(1, 14, 2):
        for ell in range(3, 25):

            Q, R = divide_by_linear(k, ell)

            key = "same" if (k - ell) % 2 == 0 else "different"

            summary[key]["total"] += 1

            if R == 0:
                summary[key]["divisible"] += 1
            else:
                summary[key]["nondivisible"] += 1

    for key in ("same", "different"):

        item = summary[key]

        print(
            f"{key:10s}: "
            f"total={item['total']:3d} "
            f"divisible={item['divisible']:3d} "
            f"nondivisible={item['nondivisible']:3d}"
        )


# ----------------------------------------------------------------------------
# REMAINDER SHAPE ANALYSIS
# ----------------------------------------------------------------------------

def remainder_analysis():

    print()
    print("=" * 78)
    print("2. NONDIVISIBLE-BRANCH REMAINDER ANALYSIS")
    print("=" * 78)

    test_cases = [
        (1, 4),
        (1, 6),
        (1, 8),
        (1, 10),
        (3, 8),
        (3, 10),
        (3, 12),
        (5, 10),
        (5, 12),
        (7, 14),
        (9, 20),
        (11, 24),
    ]

    for k, ell in test_cases:

        Q, R = divide_by_linear(k, ell)

        print()
        print(f"k={k}, ell={ell}")

        print(
            f"  R = {R}"
        )

        mult, residual = q_factor_power(R)

        if mult is not None:
            print(
                f"  (q+1)-multiplicity = {mult}"
            )
            print(
                f"  residual = {residual}"
            )


# ----------------------------------------------------------------------------
# COMPARE REMAINDER TO CLOSED EXPRESSIONS
# ----------------------------------------------------------------------------

def remainder_closed_shape():

    print()
    print("=" * 78)
    print("3. EXACT REMAINDER IDENTITY TEST")
    print("=" * 78)

    # From the root calculation:
    #
    #   F(-q-1,q)
    #
    # can be evaluated directly.
    #
    # We test whether the observed remainder is exactly this polynomial,
    # which it must be modulo p+q+1.

    for k, ell in [
        (1, 4),
        (1, 6),
        (3, 8),
        (3, 10),
        (5, 10),
        (7, 14),
    ]:

        _, R = divide_by_linear(k, ell)

        root = root_test(k, ell)

        # R is a polynomial in q.  Verify exact equality.
        ok = sp.expand(R - root) == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"R == F(-q-1,q): "
            f"{'YES' if ok else 'NO'}"
        )


# ----------------------------------------------------------------------------
# QUOTIENT SYMMETRY TEST
# ----------------------------------------------------------------------------

def quotient_symmetry_analysis():

    print()
    print("=" * 78)
    print("4. QUOTIENT SYMMETRY IN THE DIVISIBLE BRANCH")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 5),
        (1, 7),
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 9),
        (5, 11),
        (5, 19),
        (7, 13),
        (7, 15),
        (9, 19),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        Q, R = divide_by_linear(k, ell)

        if R != 0:
            continue

        symmetric = is_symmetric_pq(Q)

        print(
            f"k={k:2d} ell={ell:2d} "
            f"Q symmetric={'YES' if symmetric else 'NO'}"
        )

        if not symmetric:
            swapped = swap_pq(Q)

            difference = sp.expand(Q - swapped)

            print(
                f"  antisymmetry defect = {difference}"
            )


# ----------------------------------------------------------------------------
# ALTERNATIVE NORMALIZATION TEST
# ----------------------------------------------------------------------------

def alternative_normalization():

    print()
    print("=" * 78)
    print("5. DIRECT FACTORIZATION BY (p-q)")
    print("=" * 78)

    cases = [
        (1, 4),
        (1, 6),
        (3, 8),
        (3, 10),
        (5, 10),
        (7, 14),
        (1, 3),
        (3, 7),
        (5, 11),
    ]

    for k, ell in cases:

        Fpoly = sp.Poly(exact_F(k, ell), p, q)

        D = sp.Poly(p - q, p, q)

        Q, R = sp.div(Fpoly, D)

        Q = sp.expand(Q.as_expr())
        R = sp.expand(R.as_expr())

        print(
            f"k={k:2d} ell={ell:2d} "
            f"(p-q) divides = "
            f"{'YES' if R == 0 else 'NO'}"
        )

        if R != 0:
            print(
                f"  remainder = {R}"
            )


# ----------------------------------------------------------------------------
# ANCHOR BRANCH CLASSIFICATION
# ----------------------------------------------------------------------------

def anchor_classification():

    print()
    print("=" * 78)
    print("6. THEOREM ANCHOR BRANCH CLASSIFICATION")
    print("=" * 78)

    anchors = [
        ("A1", 1, 10, 2, 27),
        ("A2", 3, 20, 6, 935),
        ("A3", 7, 40, 15, -1797818),
        ("B1", 3, 7, 3, -3),
        ("B2", 3, 9, 4, 9),
        ("B3", 3, 11, 5, -16),
        ("B4", 5, 11, 5, -5),
        ("B5", 5, 19, 9, -196),
        ("B6", 9, 21, 11, -84),
    ]

    for name, k, ell, s, expected in anchors:

        Q, R = divide_by_linear(k, ell)

        print(
            f"{name:2s}: "
            f"k={k:2d} ell={ell:2d} s={s:2d} "
            f"parity={'same' if (k-ell)%2==0 else 'different'} "
            f"divisible={'YES' if R == 0 else 'NO'} "
            f"expected={expected}"
        )


# ----------------------------------------------------------------------------
# EXACT KERNEL SPECIALIZATION
# ----------------------------------------------------------------------------

def direct_kernel_checks():

    print()
    print("=" * 78)
    print("7. DIRECT SPECIALIZATIONS")
    print("=" * 78)

    # These identities are useful because they do not involve Q at all.

    for k, ell in [
        (1, 4),
        (1, 10),
        (3, 8),
        (5, 10),
    ]:

        F = exact_F(k, ell)

        print()
        print(f"k={k}, ell={ell}")

        print(
            "  F(-q-1,q) ="
        )

        print(
            f"    {root_test(k, ell)}"
        )


# ----------------------------------------------------------------------------
# FINAL DIAGNOSTIC
# ----------------------------------------------------------------------------

def final_diagnostic():

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    same_total = 0
    same_div = 0

    different_total = 0
    different_div = 0

    for k in range(1, 14, 2):
        for ell in range(3, 25):

            _, R = divide_by_linear(k, ell)

            if (k - ell) % 2 == 0:

                same_total += 1

                if R == 0:
                    same_div += 1

            else:

                different_total += 1

                if R == 0:
                    different_div += 1

    print(
        f"same-parity cases       = {same_total}"
    )

    print(
        f"same-parity divisible   = {same_div}"
    )

    print(
        f"different-parity cases  = {different_total}"
    )

    print(
        f"different-parity divide= {different_div}"
    )

    print()

    if same_div == same_total and different_div == 0:

        print(
            "OBSERVED PARITY LAW:"
        )

        print(
            "    p+q+1 divides F exactly on the same-parity branch."
        )

        print()

        print(
            "Therefore the universal object is NOT Q=F/(p+q+1)."
        )

        print(
            "The theorem's coefficient must come from a normalization"
        )

        print(
            "that remains valid on the opposite-parity branch."
        )

        print()

        print(
            "NEXT MATHEMATICAL TARGET:"
        )

        print(
            "derive a single parity-unified symmetric normalization of F"
        )

        print(
            "before attempting any coefficient/index identification."
        )

    else:

        print(
            "The simple parity law is not sufficient."
        )

        print(
            "Further structural analysis of the exact kernel is required."
        )


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 202")
    print("EXACT pq-KERNEL NORMALIZATION BY PARITY")
    print("=" * 78)

    print()
    print(
        "This experiment does not construct Q(N,S) unless the"
    )
    print(
        "necessary divisibility has been independently verified."
    )

    parity_summary()
    remainder_analysis()
    remainder_closed_shape()
    quotient_symmetry_analysis()
    alternative_normalization()
    anchor_classification()
    direct_kernel_checks()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 202")
    print("=" * 78)


if __name__ == "__main__":
    main()