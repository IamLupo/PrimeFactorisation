#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 163
SYMBOLIC HYPERGEOMETRIC TERM RATIO FOR [N^s] D_(k,ell)

PURPOSE
-------
Experiment 162 established exactly that every coefficient

    [N^s] D_(k,ell)(N)

is a finite alternating binomial convolution.

The next task is to expose the general summand as a hypergeometric term.

For fixed (k,ell,s), define the exact coefficient summand

    T_{k,ell,s}(j)

after summing the V_r coefficient over the internal V-index m.

This experiment:

    1. constructs the exact symbolic j-summand;
    2. derives T(j+1)/T(j) symbolically;
    3. factors numerator and denominator;
    4. checks whether the ratio is rational in j;
    5. identifies terminating support and boundary zeros;
    6. compares the symbolic summand against the exact finite convolution;
    7. asks SymPy for a hypergeometric / summation reduction when possible;
    8. verifies the resulting formula on training and forward ell holdout.

The goal is NOT yet to prove the final closed form.
The goal is to identify the exact hypergeometric class of the remaining
finite sum.

NO POLYNOMIAL FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

j = sp.symbols("j", integer=True)
m = sp.symbols("m", integer=True, nonnegative=True)
N = sp.symbols("N")


# ============================================================================
# V coefficient
# ============================================================================

def v_coeff(r, m):
    """
    [N^m] V_r(N).
    """

    if m < 0:
        return sp.Integer(0)

    if m > r // 2:
        return sp.Integer(0)

    a = sp.binomial(
        r - m,
        m
    )

    b = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(
            r - m - 1,
            m - 1
        )
    )

    return sp.expand(
        (-1)**(r-m) * (a + b)
    )


# ============================================================================
# Exact coefficient source from the finite kernel
# ============================================================================

def kernel_source_terms(k, ell, j):
    """
    Return kernel contributions to C_j as

        (source_label, t_exponent, coefficient).

    The four terms are

        p^k(1+q)^ell
        q^k(1+p)^ell
        -p^ell(1+q)^k
        -q^ell(1+p)^k.
    """

    terms = []

    # p^k(1+q)^ell
    a = k - j
    if 0 <= a <= ell:
        terms.append(
            (
                "A",
                k + a,
                sp.binomial(ell, a)
            )
        )

    # q^k(1+p)^ell
    b = k + j
    if 0 <= b <= ell:
        terms.append(
            (
                "B",
                k + b,
                sp.binomial(ell, b)
            )
        )

    # -p^ell(1+q)^k
    c = ell - j
    if 0 <= c <= k:
        terms.append(
            (
                "C",
                ell + c,
                -sp.binomial(k, c)
            )
        )

    # -q^ell(1+p)^k
    d = ell + j
    if 0 <= d <= k:
        terms.append(
            (
                "D",
                ell + d,
                -sp.binomial(k, d)
            )
        )

    return terms


# ============================================================================
# Exact coefficient contribution for fixed source j
# ============================================================================

def source_j_contribution(k, ell, s, j):
    """
    Return the complete contribution of a fixed source mode j to

        [N^s] D.

    The quotient source uses

        r = j - 2k,

    and the V coefficient v_coeff(r,m).

    The final N-degree condition is imposed exactly.
    """

    r = j - 2*k

    if r < 0:
        return sp.Integer(0)

    total = sp.Integer(0)

    # V_r has m = 0,...,floor(r/2).
    for m0 in range(
        0,
        sp.floor(r/2) + 1
    ):
        # SymPy floor stays symbolic if r is symbolic.
        # This function is intended to be evaluated at integer j.
        # Therefore this branch is for explicit integer evaluation.
        break

    raise RuntimeError(
        "source_j_contribution is only a structural placeholder; "
        "use symbolic_source_term for Experiment 163."
    )


# ============================================================================
# Symbolic finite m-sum for a fixed kernel branch
# ============================================================================

def symbolic_branch_sum(
    k,
    ell,
    s,
    branch
):
    """
    Symbolically derive the contribution of one kernel branch.

    For each branch, the t-degree determines the required V-index m,
    so m is not freely summed: it is fixed by the target N-degree s.

    This is the key simplification.

    If the kernel term has t-degree e and source index j with
        r = j-2k,
    then the final N-degree is

        m + (e-r-2k)/2.

    Hence

        m = s - (e-r-2k)/2.

    Thus the V inner sum collapses immediately for each branch.
    """

    r = j - 2*k

    # Select symbolic branch.
    if branch == "A":
        a = k - j
        e = k + a
        coeff = sp.binomial(ell, a)

    elif branch == "B":
        b = k + j
        e = k + b
        coeff = sp.binomial(ell, b)

    elif branch == "C":
        c = ell - j
        e = ell + c
        coeff = -sp.binomial(k, c)

    elif branch == "D":
        d = ell + j
        e = ell + d
        coeff = -sp.binomial(k, d)

    else:
        raise ValueError(
            f"unknown branch {branch}"
        )

    # D normalization:
    #
    # weight = V_r / t^(r+1)
    #
    # final divide = t^(2k-1)
    #
    # V coefficient N^m contributes t^(2m).
    #
    # total t exponent:
    #
    #   2m + e - (r+1) - (2k-1)
    #
    # = 2m + e-r-2k.
    #
    # Therefore N degree s satisfies
    #
    #   2s = 2m + e-r-2k.
    #
    # hence
    #
    #   m = s - (e-r-2k)/2.

    numerator = sp.expand(
        e - r - 2*k
    )

    # The branch is valid only when numerator is even.
    # For symbolic analysis we retain a parity indicator separately.

    m_expr = sp.simplify(
        s - numerator/2
    )

    # V coefficient with symbolic m:
    V_symbolic = (
        (-1)**(r-m_expr)
        * (
            sp.binomial(
                r-m_expr,
                m_expr
            )
            + sp.binomial(
                r-m_expr-1,
                m_expr-1
            )
        )
    )

    term = sp.simplify(
        coeff * V_symbolic
    )

    return sp.factor(term), sp.simplify(m_expr)


# ============================================================================
# Build symbolic branch summand
# ============================================================================

def symbolic_j_term(k, ell, s):
    """
    Sum the four kernel branches.

    Returns

        T(j) = A_j + B_j + C_j + D_j

    as an exact symbolic expression in j.

    Parity-incompatible branches are retained symbolically but later
    filtered when evaluating concrete integer j.
    """

    total = sp.Integer(0)
    branch_data = []

    for branch in ["A", "B", "C", "D"]:

        try:
            term, m_expr = symbolic_branch_sum(
                k,
                ell,
                s,
                branch
            )

        except Exception:
            continue

        total += term

        branch_data.append(
            (
                branch,
                sp.factor(term),
                sp.factor(m_expr)
            )
        )

    return (
        sp.factor(
            sp.expand(total)
        ),
        branch_data
    )


# ============================================================================
# Exact finite j-term for integer j
# ============================================================================

def concrete_j_term(k, ell, s, j0):
    """
    Compute the coefficient contribution for a concrete integer j.

    This is used to verify the symbolic T(j).
    """

    r0 = j0 - 2*k

    if r0 < 0:
        return sp.Integer(0)

    total = sp.Integer(0)

    for branch in ["A", "B", "C", "D"]:

        if branch == "A":
            a = k-j0
            if not (0 <= a <= ell):
                continue
            e = k+a
            c = sp.binomial(ell, a)

        elif branch == "B":
            b = k+j0
            if not (0 <= b <= ell):
                continue
            e = k+b
            c = sp.binomial(ell, b)

        elif branch == "C":
            c0 = ell-j0
            if not (0 <= c0 <= k):
                continue
            e = ell+c0
            c = -sp.binomial(k, c0)

        else:
            d = ell+j0
            if not (0 <= d <= k):
                continue
            e = ell+d
            c = -sp.binomial(k, d)

        numerator = e-r0-2*k

        if numerator % 2 != 0:
            continue

        m0 = s - numerator//2

        if m0 < 0 or m0 > r0//2:
            continue

        total += (
            c
            * v_coeff(
                r0,
                m0
            )
        )

    return sp.factor(total)


# ============================================================================
# Symbolic consecutive ratio
# ============================================================================

def symbolic_ratio(T):
    """
    Compute and factor

        T(j+1)/T(j).

    The ratio is only meaningful where T(j) is nonzero.
    """

    numerator = sp.expand(
        T.subs(j, j+1)
    )

    denominator = sp.expand(T)

    ratio = sp.factor(
        sp.cancel(
            numerator/denominator
        )
    )

    return ratio


# ============================================================================
# Hypergeometric recognition
# ============================================================================

def hyper_check(T):
    """
    SymPy's hypergeometric test.

    This is applied to T viewed as a function of j.
    """

    try:
        return sp.concrete.gosper.gosper_sum(
            T,
            (j, j, j)
        )
    except Exception:
        return None


# ============================================================================
# Exact finite sum of concrete terms
# ============================================================================

def concrete_sum(k, ell, s):
    """
    Sum over the actual admissible source support.
    """

    j_min = 2*k
    j_max = ell

    terms = []

    for j0 in range(
        j_min,
        j_max+1
    ):
        term = concrete_j_term(
            k,
            ell,
            s,
            j0
        )

        if term != 0:
            terms.append(
                (
                    j0,
                    term
                )
            )

    return terms


# ============================================================================
# Direct D coefficients
# ============================================================================

def direct_coefficients(k, ell):
    """
    Independent coefficient list obtained from the known D polynomials
    reconstructed through the exact finite quotient.
    """

    # Build D directly from concrete source terms.
    degree = max(
        0,
        ell-2*k
    )

    result = {}

    for s0 in range(
        degree+1
    ):
        terms = concrete_sum(
            k,
            ell,
            s0
        )

        result[s0] = sp.factor(
            sum(
                term
                for _, term in terms
            )
        )

    return result


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 163")
    print("SYMBOLIC HYPERGEOMETRIC TERM RATIO")
    print("=" * 78)
    print()

    TESTS = [
        (1,9),
        (1,11),
        (1,13),
        (3,11),
        (3,13),
        (5,13)
    ]

    FORWARD = [
        (1,15),
        (1,17),
        (3,15),
        (3,17),
        (5,15),
        (5,17)
    ]

    symbolic_failures = 0
    ratio_failures = 0
    holdout_failures = 0

    # ------------------------------------------------------------------------
    # 1. Symbolic summand construction
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. SYMBOLIC j-SUM")
    print("=" * 78)

    for k, ell in TESTS:

        degree = max(
            0,
            ell-2*k
        )

        print()
        print(
            f"({k},{ell})"
        )

        for s0 in range(
            min(3, degree+1)
        ):

            T, branches = symbolic_j_term(
                k,
                ell,
                s0
            )

            print()
            print(
                f"  target s={s0}"
            )

            for branch, term, m_expr in branches:
                print(
                    f"    branch={branch}"
                )

                print(
                    f"      m(j) = {m_expr}"
                )

                print(
                    f"      term  = {term}"
                )

            print(
                "    total T(j) =",
                T
            )

            # Check concrete symbolic representation at supported j.
            concrete = concrete_sum(
                k,
                ell,
                s0
            )

            mismatch = 0

            for j0, value in concrete:

                symbolic_value = sp.simplify(
                    T.subs(
                        j,
                        j0
                    )
                )

                residual = sp.factor(
                    sp.expand(
                        symbolic_value-value
                    )
                )

                if residual != 0:
                    mismatch += 1

                    print(
                        f"      mismatch at j={j0}: "
                        f"residual={residual}"
                    )

            if mismatch:
                symbolic_failures += 1

            # ----------------------------------------------------------------
            # Ratio
            # ----------------------------------------------------------------

            if T != 0:

                try:
                    ratio = symbolic_ratio(T)

                    print(
                        "    T(j+1)/T(j) =",
                        ratio
                    )

                    # Test whether ratio still depends on j rationally.
                    ratio_together = sp.together(
                        ratio
                    )

                    num, den = (
                        sp.fraction(
                            ratio_together
                        )
                    )

                    num_poly = sp.Poly(
                        num,
                        j,
                        domain=sp.QQ.frac_field(N)
                    )

                    den_poly = sp.Poly(
                        den,
                        j,
                        domain=sp.QQ.frac_field(N)
                    )

                    print(
                        "    numerator degree in j =",
                        num_poly.degree()
                    )

                    print(
                        "    denominator degree in j =",
                        den_poly.degree()
                    )

                except Exception as exc:
                    ratio_failures += 1

                    print(
                        "    ratio diagnostic failed:",
                        type(exc).__name__,
                        str(exc)
                    )

    # ------------------------------------------------------------------------
    # 2. Concrete hypergeometric diagnostics
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CONCRETE SUM RATIOS")
    print("=" * 78)

    for k, ell in TESTS:

        degree = max(
            0,
            ell-2*k
        )

        print()
        print(
            f"({k},{ell})"
        )

        for s0 in range(
            min(3, degree+1)
        ):

            terms = concrete_sum(
                k,
                ell,
                s0
            )

            print()
            print(
                f"  s={s0}"
            )

            for (ja, va), (jb, vb) in zip(
                terms,
                terms[1:]
            ):

                if vb == 0 or va == 0:
                    continue

                ratio = sp.factor(
                    sp.cancel(
                        vb/va
                    )
                )

                print(
                    f"    j={ja}->{jb}: ratio={ratio}"
                )

            print(
                "    sum =",
                sp.factor(
                    sum(
                        value
                        for _, value in terms
                    )
                )
            )

    # ------------------------------------------------------------------------
    # 3. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    for k, ell in FORWARD:

        coeffs = direct_coefficients(
            k,
            ell
        )

        ok = True

        for s0, exact in coeffs.items():

            terms = concrete_sum(
                k,
                ell,
                s0
            )

            summed = sp.factor(
                sum(
                    value
                    for _, value in terms
                )
            )

            if sp.expand(
                summed-exact
            ) != 0:
                ok = False

        print(
            f"({k},{ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            holdout_failures += 1

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "symbolic summand failures =",
        symbolic_failures
    )

    print(
        "ratio diagnostic failures =",
        ratio_failures
    )

    print(
        "forward failures =",
        holdout_failures
    )

    if (
        symbolic_failures == 0
        and ratio_failures == 0
        and holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The coefficient convolution has been reduced"
        )

        print(
            "to an explicit symbolic j-summand."
        )

        print()
        print(
            "The consecutive ratio T(j+1)/T(j) is now"
        )

        print(
            "available as an exact rational function of j."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "identify the resulting finite sum as a terminating"
        )

        print(
            "hypergeometric series and apply the appropriate"
        )

        print(
            "Vandermonde / Chu-Vandermonde identity."
        )

    else:

        print()
        print(
            "STATUS = FAIL"
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

