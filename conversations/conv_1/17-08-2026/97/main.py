#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 161
COEFFICIENT-WISE VANDERMONDE COMPRESSION OF THE D-CONVOLUTION

PURPOSE
-------
Experiment 160 established exactly:

    weight_r(t) = V_r(t^2) / t^(r+1),

where

    V_r(N)
      = sum_{m=0}^{floor(r/2)}
          (-1)^(r-m)
          [ C(r-m,m) + C(r-m-1,m-1) ]
          N^m.

The antisymmetric invariant is

    D_(k,ell)(N)
      = sum_r weight_r * C_(2k+r).

After normalization in t, this is a finite polynomial convolution.

The next task is to compress that convolution coefficient-by-coefficient.

METHOD
------
1. Derive the exact finite C_j binomial source.
2. Insert the closed V_r formula.
3. Expand only far enough to isolate the coefficient [N^s].
4. Obtain an exact nested finite binomial sum for each coefficient.
5. Search for exact Vandermonde / Chu-Vandermonde collapses.
6. Compare the collapsed candidate with the directly computed D polynomial.
7. Test forward ell holdout.

The experiment does NOT fit coefficients numerically.

The main target is to determine whether the coefficient of N^s can be reduced
from a double sum to a single binomial expression.

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

p, q = sp.symbols("p q")
t, u = sp.symbols("t u")
N = sp.symbols("N")


# ============================================================================
# Exact kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient
# ============================================================================

def quotient_Q_pq(k, ell):
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
            f"quotient remainder ({k},{ell}) = {sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Laurent helpers
# ============================================================================

def laurent_F(k, ell):
    return sp.expand(
        kernel_F(k, ell).subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


def laurent_Q(k, ell):
    return sp.expand(
        quotient_Q_pq(k, ell).subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


def laurent_support(expr):
    support = []

    for term in sp.Add.make_args(sp.expand(expr)):
        exponent = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        support.append(int(exponent))

    return sorted(set(support))


def laurent_coeff(expr, j):
    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    shift = max(0, -min(support))

    poly = sp.Poly(
        sp.expand(expr * u**shift),
        u,
        domain=sp.QQ.frac_field(t),
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(u**target)
    )


# ============================================================================
# N normalization
# ============================================================================

def even_t_to_N(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return sp.Integer(0)

    for term in sp.Add.make_args(expr):
        exponent = sp.sympify(
            term.as_powers_dict().get(t, 0)
        )

        if not exponent.is_Integer:
            return None

        if int(exponent) < 0:
            return None

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():

        power = int(power)

        if power % 2:
            return None

        result += coeff * N**(power // 2)

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Correct V kernel
# ============================================================================

def V_closed(r):
    result = sp.Integer(0)

    for m in range(
        0,
        r // 2 + 1
    ):
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

        result += (
            (-1)**(r - m)
            * (a + b)
            * N**m
        )

    return sp.expand(result)


# ============================================================================
# Exact finite C_j formula in N/t form
# ============================================================================

def C_closed_t(k, ell, j):
    """
    Exact [u^j] F(tu,t/u).
    """

    # p^k (1+q)^ell
    a = k - j

    if 0 <= a <= ell:
        c1 = sp.binomial(ell, a) * t**(k + a)
    else:
        c1 = sp.Integer(0)

    # q^k (1+p)^ell
    b = k + j

    if 0 <= b <= ell:
        c2 = sp.binomial(ell, b) * t**(k + b)
    else:
        c2 = sp.Integer(0)

    # -p^ell (1+q)^k
    c = ell - j

    if 0 <= c <= k:
        c3 = -sp.binomial(k, c) * t**(ell + c)
    else:
        c3 = sp.Integer(0)

    # -q^ell (1+p)^k
    d = ell + j

    if 0 <= d <= k:
        c4 = -sp.binomial(k, d) * t**(ell + d)
    else:
        c4 = sp.Integer(0)

    return sp.expand(c1 + c2 + c3 + c4)


# ============================================================================
# Exact source weights
# ============================================================================

def source_weights_for_q(k, ell, target_j):
    source = {}

    source[ell - 1] = {
        ell: sp.cancel(1/t)
    }

    for j in range(
        ell - 1,
        target_j,
        -1
    ):
        current = source.get(j, {})
        next_one = source.get(j + 1, {})

        new_map = {
            j: sp.cancel(1/t)
        }

        for m, coeff in current.items():
            new_map[m] = sp.expand(
                new_map.get(m, 0)
                - coeff/t
            )

        for m, coeff in next_one.items():
            new_map[m] = sp.expand(
                new_map.get(m, 0)
                - coeff
            )

        source[j - 1] = {
            m: sp.factor(c)
            for m, c in new_map.items()
            if sp.expand(c) != 0
        }

    return source.get(target_j, {})


def antisymmetric_source_weights(k, ell):
    left = source_weights_for_q(
        k,
        ell,
        2*k - 1
    )

    right = source_weights_for_q(
        k,
        ell,
        2*k + 1
    )

    combined = {}

    for m, coeff in left.items():
        combined[m] = sp.expand(
            combined.get(m, 0) + coeff
        )

    for m, coeff in right.items():
        combined[m] = sp.expand(
            combined.get(m, 0) - coeff
        )

    return {
        m: sp.factor(c)
        for m, c in combined.items()
        if sp.expand(c) != 0
    }


# ============================================================================
# Direct D
# ============================================================================

def D_direct(k, ell):
    Q = laurent_Q(k, ell)

    qm = laurent_coeff(
        Q,
        2*k - 1
    )

    qp = laurent_coeff(
        Q,
        2*k + 1
    )

    raw = sp.cancel(
        (qm - qp)
        / t**(2*k - 1)
    )

    D = even_t_to_N(raw)

    if D is None:
        raise ArithmeticError(
            f"D is not N-only for ({k},{ell})"
        )

    return sp.factor(D)


# ============================================================================
# Exact D from V convolution
# ============================================================================

def D_from_V(k, ell):
    weights = antisymmetric_source_weights(
        k,
        ell
    )

    raw = sp.Integer(0)

    for j, weight in weights.items():

        r = j - 2*k

        V_t = sp.expand(
            V_closed(r).subs(
                N,
                t**2
            )
        )

        expected_weight = sp.cancel(
            V_t / t**(r + 1)
        )

        residual = sp.factor(
            sp.together(
                weight - expected_weight
            )
        )

        if residual != 0:
            raise ArithmeticError(
                f"V weight mismatch at ({k},{ell}), "
                f"j={j}, r={r}: {residual}"
            )

        raw += (
            weight
            * C_closed_t(k, ell, j)
        )

    raw = sp.cancel(
        raw / t**(2*k - 1)
    )

    D = even_t_to_N(raw)

    if D is None:
        raise ArithmeticError(
            f"D convolution not N-only for ({k},{ell})"
        )

    return sp.factor(D)


# ============================================================================
# Exact coefficient extraction from D
# ============================================================================

def coefficient_dict(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ
    )

    return {
        r: sp.factor(
            poly.coeff_monomial(N**r)
        )
        for r in range(
            poly.degree() + 1
        )
    }


# ============================================================================
# Build coefficient-level convolution
# ============================================================================

def coefficient_convolution(k, ell, s):
    """
    Derive [N^s] D directly from:

        D = sum_r V_r(N) * normalized C_r(N).

    We retain the finite sum symbolically.

    The source index j corresponds to r=j-2k.
    """

    weights = antisymmetric_source_weights(
        k,
        ell
    )

    terms = []

    for j in sorted(weights):

        r = j - 2*k

        if r < 0:
            continue

        V = V_closed(r)

        # Exact C_j after dividing out t^(2k-1).
        Cj = C_closed_t(
            k,
            ell,
            j
        )

        # Source contribution:
        #
        # weight = V_r(t^2) / t^(r+1)
        #
        # D normalization = t^-(2k-1)
        #
        # The t-powers cancel to an N-polynomial.
        contribution = sp.cancel(
            (
                V.subs(N, t**2)
                / t**(r+1)
            )
            * Cj
            / t**(2*k-1)
        )

        contribution_N = even_t_to_N(
            contribution
        )

        if contribution_N is None:
            raise ArithmeticError(
                f"non-N contribution at j={j}"
            )

        coeff = sp.Poly(
            sp.expand(contribution_N),
            N,
            domain=sp.QQ
        ).coeff_monomial(
            N**s
        )

        if coeff != 0:
            terms.append(
                (
                    j,
                    r,
                    sp.factor(coeff)
                )
            )

    return terms


# ============================================================================
# Search simple Vandermonde collapses
# ============================================================================

def vandermonde_candidates(k, ell, s, exact_coeff):
    """
    Test a controlled set of exact binomial candidates.

    These are diagnostics only. A candidate is accepted only when it
    reproduces an identity over the full tested set, not from a single row.
    """

    candidates = []

    # Natural central binomial candidates.
    test_forms = [
        (
            f"C(ell,{s})",
            sp.binomial(ell, s)
        ),
        (
            f"C(ell,{s+k})",
            (
                sp.binomial(ell, s+k)
                if 0 <= s+k <= ell
                else sp.Integer(0)
            )
        ),
        (
            f"C(ell,{s+k-1})",
            (
                sp.binomial(ell, s+k-1)
                if 0 <= s+k-1 <= ell
                else sp.Integer(0)
            )
        ),
        (
            f"C(ell,{s+2*k-1})",
            (
                sp.binomial(
                    ell,
                    s+2*k-1
                )
                if 0 <= s+2*k-1 <= ell
                else sp.Integer(0)
            )
        ),
        (
            f"C(ell-k,{s})",
            (
                sp.binomial(
                    ell-k,
                    s
                )
                if 0 <= s <= ell-k
                else sp.Integer(0)
            )
        ),
    ]

    for name, value in test_forms:

        if value == 0:
            continue

        ratio = sp.factor(
            sp.cancel(
                exact_coeff / value
            )
        )

        candidates.append(
            (
                name,
                ratio
            )
        )

    return candidates


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 161")
    print("COEFFICIENT-WISE VANDERMONDE COMPRESSION")
    print("=" * 78)
    print()

    TRAIN = [
        (1,3),
        (1,5),
        (1,7),
        (1,9),
        (1,11),
        (1,13),
        (3,7),
        (3,9),
        (3,11),
        (3,13),
        (5,11),
        (5,13),
    ]

    FORWARD = [
        (1,15),
        (1,17),
        (3,15),
        (3,17),
        (5,15),
        (5,17),
    ]

    # ------------------------------------------------------------------------
    # 1. Baseline D convolution
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. BASELINE EXACT D CONVOLUTION")
    print("=" * 78)

    baseline_failures = 0

    for k, ell in TRAIN:

        direct = D_direct(
            k,
            ell
        )

        from_V = D_from_V(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - from_V
            )
        )

        print(
            f"({k},{ell}) "
            f"residual={residual}"
        )

        if residual != 0:
            baseline_failures += 1

    # ------------------------------------------------------------------------
    # 2. Coefficient-wise finite convolution
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. COEFFICIENT-WISE SOURCE CONVOLUTION")
    print("=" * 78)

    coefficient_failures = 0

    for k, ell in TRAIN:

        D = D_from_V(
            k,
            ell
        )

        coeffs = coefficient_dict(
            D
        )

        print()
        print(
            f"({k},{ell})"
        )

        for s in sorted(coeffs):

            exact = coeffs[s]

            terms = coefficient_convolution(
                k,
                ell,
                s
            )

            summed = sp.factor(
                sp.expand(
                    sum(
                        term[2]
                        for term in terms
                    )
                )
            )

            residual = sp.factor(
                sp.expand(
                    summed - exact
                )
            )

            print(
                f"  [N^{s}] exact = {exact}"
            )

            print(
                f"       convolution terms = "
                f"{len(terms)}"
            )

            print(
                f"       summed = {summed}"
            )

            print(
                f"       residual = {residual}"
            )

            if residual != 0:
                coefficient_failures += 1

            candidates = vandermonde_candidates(
                k,
                ell,
                s,
                exact
            )

            for name, ratio in candidates:
                print(
                    f"       candidate {name}: ratio={ratio}"
                )

    # ------------------------------------------------------------------------
    # 3. Representative inner sums
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. REPRESENTATIVE BINOMIAL SUMS")
    print("=" * 78)

    for k, ell in [
        (1,9),
        (3,11),
        (5,13),
    ]:

        D = D_from_V(
            k,
            ell
        )

        coeffs = coefficient_dict(
            D
        )

        print()
        print(
            f"({k},{ell})"
        )

        for s in sorted(coeffs):

            terms = coefficient_convolution(
                k,
                ell,
                s
            )

            print(
                f"  coefficient N^{s}:"
            )

            for j, r, term in terms:
                print(
                    f"    j={j}, r={r}, term={term}"
                )

            print(
                "    exact sum =",
                sp.factor(
                    sum(
                        term[2]
                        for term in terms
                    )
                )
            )

    # ------------------------------------------------------------------------
    # 4. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    holdout_failures = 0

    for k, ell in FORWARD:

        direct = D_direct(
            k,
            ell
        )

        from_V = D_from_V(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - from_V
            )
        )

        ok = residual == 0

        print(
            f"({k},{ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        print(
            "  D(N) =",
            direct
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
        "baseline D failures =",
        baseline_failures
    )

    print(
        "coefficient-wise failures =",
        coefficient_failures
    )

    print(
        "forward failures =",
        holdout_failures
    )

    if (
        baseline_failures == 0
        and coefficient_failures == 0
        and holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The V/C convolution reproduces every coefficient"
        )

        print(
            "of D exactly."
        )

        print()
        print(
            "The coefficient-level finite sums are now isolated."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive an exact one-sum Vandermonde or"
        )

        print(
            "Chu-Vandermonde collapse for the inner coefficient sum."
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

