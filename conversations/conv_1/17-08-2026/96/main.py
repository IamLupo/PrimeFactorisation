#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 160
CORRECTED V-KERNEL SOURCE MATCH + EXACT D CONVOLUTION
==============================================================================

PURPOSE
-------
Experiment 159 established the corrected antisymmetric Green kernel

    V_r = W_r - N W_(r-2)

with

    V_0 = 1
    V_1 = -1
    V_(r+1) = -V_r - N V_(r-1),

and the closed form

    V_r(N)
      = sum_{m=0}^{floor(r/2)}
          (-1)^(r-m)
          [ C(r-m,m) + C(r-m-1,m-1) ]
          N^m.

The previous run failed only because the source-weight checker compared
expressions containing t against expressions containing N without imposing

    N = t^2.

For example,

    -(2*t^2-1)/t^3

and

    (1-2*N)/t^3

are exactly identical after N=t^2.

This experiment fixes that normalization.

TARGET
------
1. Verify W recurrence and closed form.
2. Verify V recurrence, W-NW shift, and closed form.
3. Compare actual source weights against

       V_r(t^2) / t^(r+1)

   directly in the t-domain.
4. Reconstruct D from the corrected V kernel and exact finite-binomial
   kernel coefficients.
5. Compare with the independently constructed quotient.
6. Verify D(0)=1.
7. Verify forward ell holdout.
8. Print the exact finite convolution for representative pairs.

IMPORTANT
---------
This experiment does not use polynomial fitting.

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
# Exact finite kernel
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

    remainder = sp.expand(R.as_expr())

    if remainder != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(remainder)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Laurent expressions
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

    for term in sp.Add.make_args(
        sp.expand(expr)
    ):
        exponent = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent in u: {term}"
            )

        support.append(int(exponent))

    return sorted(set(support))


def laurent_coeff(expr, j):
    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    shift = max(
        0,
        -min(support)
    )

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
# Robust t-Laurent decomposition
# ============================================================================

def split_t_laurent(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return 0, sp.Integer(0)

    minimum = None

    for term in sp.Add.make_args(expr):
        exponent = sp.sympify(
            term.as_powers_dict().get(t, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"noninteger exponent of t: {term}"
            )

        exponent = int(exponent)

        if minimum is None or exponent < minimum:
            minimum = exponent

    reduced = sp.expand(
        expr / t**minimum
    )

    return minimum, reduced


def even_t_to_N(expr):
    """
    Convert an ordinary even polynomial in t to N=t^2.
    """

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
        domain=sp.QQ,
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():
        power = int(power)

        if power % 2 != 0:
            return None

        result += coeff * N**(power // 2)

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Universal W kernel
# ============================================================================

def W_closed(r):
    return sp.factor(
        sp.expand(
            sum(
                (-1)**(r-m)
                * sp.binomial(r-m, m)
                * N**m
                for m in range(
                    0,
                    r // 2 + 1
                )
            )
        )
    )


def W_recursive(max_r):
    W = {
        0: sp.Integer(1),
        1: sp.Integer(-1),
    }

    for r in range(1, max_r):
        W[r + 1] = sp.factor(
            sp.expand(
                -W[r]
                - N*W[r-1]
            )
        )

    return W


# ============================================================================
# Correct antisymmetric V kernel
# ============================================================================

def V_recursive(max_r):
    V = {
        0: sp.Integer(1),
        1: sp.Integer(-1),
    }

    if max_r >= 2:
        V[2] = sp.factor(
            1 - 2*N
        )

    for r in range(2, max_r):
        V[r + 1] = sp.factor(
            sp.expand(
                -V[r]
                - N*V[r-1]
            )
        )

    return V


def V_from_W(r, W):
    if r == 0:
        return sp.Integer(1)

    if r == 1:
        return sp.Integer(-1)

    return sp.factor(
        sp.expand(
            W[r]
            - N*W[r-2]
        )
    )


def V_closed(r):
    result = sp.Integer(0)

    for m in range(
        0,
        r // 2 + 1
    ):
        first = sp.binomial(
            r-m,
            m
        )

        second = (
            sp.Integer(0)
            if m == 0
            else sp.binomial(
                r-m-1,
                m-1
            )
        )

        result += (
            (-1)**(r-m)
            * (first + second)
            * N**m
        )

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Source weights for quotient modes
# ============================================================================

def source_weights_for_q(
    k,
    ell,
    target_j
):
    source = {}

    # q_(ell-1) = C_ell / t
    source[ell - 1] = {
        ell: sp.cancel(1/t)
    }

    for j in range(
        ell - 1,
        target_j,
        -1
    ):
        current = source.get(
            j,
            {}
        )

        next_one = source.get(
            j + 1,
            {}
        )

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

    return source.get(
        target_j,
        {}
    )


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
            combined.get(m, 0)
            + coeff
        )

    for m, coeff in right.items():
        combined[m] = sp.expand(
            combined.get(m, 0)
            - coeff
        )

    return {
        m: sp.factor(c)
        for m, c in combined.items()
        if sp.expand(c) != 0
    }


# ============================================================================
# FIXED source-weight comparison
# ============================================================================

def expected_source_weight_t(r):
    """
    Expected source coefficient directly in the t-domain.

        V_r(N=t^2) / t^(r+1)
    """

    V_t = sp.expand(
        V_closed(r).subs(
            N,
            t**2
        )
    )

    return sp.cancel(
        V_t / t**(r+1)
    )


def source_weight_residual(
    k,
    m,
    actual
):
    r = m - 2*k

    if r < 0:
        return None, r

    expected = expected_source_weight_t(
        r
    )

    residual = sp.factor(
        sp.together(
            actual - expected
        )
    )

    return residual, r


# ============================================================================
# Exact finite Laurent coefficient of F
# ============================================================================

def kernel_mode_closed(k, ell, j):
    # p^k (1+q)^ell
    a = k - j

    if 0 <= a <= ell:
        c1 = (
            sp.binomial(ell, a)
            * t**(k+a)
        )
    else:
        c1 = sp.Integer(0)

    # q^k (1+p)^ell
    b = k + j

    if 0 <= b <= ell:
        c2 = (
            sp.binomial(ell, b)
            * t**(k+b)
        )
    else:
        c2 = sp.Integer(0)

    # -p^ell (1+q)^k
    c = ell - j

    if 0 <= c <= k:
        c3 = (
            -sp.binomial(k, c)
            * t**(ell+c)
        )
    else:
        c3 = sp.Integer(0)

    # -q^ell (1+p)^k
    d = ell + j

    if 0 <= d <= k:
        c4 = (
            -sp.binomial(k, d)
            * t**(ell+d)
        )
    else:
        c4 = sp.Integer(0)

    return sp.expand(
        c1 + c2 + c3 + c4
    )


# ============================================================================
# D from corrected V source
# ============================================================================

def D_from_V_source(k, ell):
    weights = antisymmetric_source_weights(
        k,
        ell
    )

    raw = sp.Integer(0)

    for m, actual_weight in weights.items():

        residual, r = source_weight_residual(
            k,
            m,
            actual_weight
        )

        if residual is None:
            raise ArithmeticError(
                f"negative Green index r={r}"
            )

        if residual != 0:
            raise ArithmeticError(
                "source weight mismatch: "
                f"pair=({k},{ell}), "
                f"m={m}, r={r}, "
                f"actual={actual_weight}, "
                f"expected={expected_source_weight_t(r)}, "
                f"residual={residual}"
            )

        raw += (
            actual_weight
            * kernel_mode_closed(
                k,
                ell,
                m
            )
        )

    # D = (q_(2k-1)-q_(2k+1)) / t^(2k-1)
    raw = sp.cancel(
        raw / t**(2*k - 1)
    )

    base, reduced = split_t_laurent(
        raw
    )

    if base != 0:
        raise ArithmeticError(
            f"D retains t-power {base} for ({k},{ell})"
        )

    D = even_t_to_N(
        reduced
    )

    if D is None:
        raise ArithmeticError(
            f"D is not N-only for ({k},{ell})"
        )

    return sp.factor(D)


# ============================================================================
# Direct D from quotient
# ============================================================================

def D_direct(k, ell):
    Q = laurent_Q(
        k,
        ell
    )

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

    base, reduced = split_t_laurent(
        raw
    )

    if base != 0:
        raise ArithmeticError(
            f"direct D retains t-power {base}"
        )

    D = even_t_to_N(
        reduced
    )

    if D is None:
        raise ArithmeticError(
            f"direct D is not N-only"
        )

    return sp.factor(D)


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 160")
    print("CORRECTED V-KERNEL SOURCE MATCH + EXACT D CONVOLUTION")
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

    green_failures = 0
    source_failures = 0
    D_failures = 0
    boundary_failures = 0
    holdout_failures = 0

    # ------------------------------------------------------------------------
    # 1. V certificate
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. V GREEN KERNEL CERTIFICATE")
    print("=" * 78)

    W = W_recursive(20)
    V = V_recursive(20)

    for r in range(21):

        from_W = V_from_W(
            r,
            W
        )

        closed = V_closed(r)

        residual_1 = sp.factor(
            sp.expand(
                V[r] - from_W
            )
        )

        residual_2 = sp.factor(
            sp.expand(
                V[r] - closed
            )
        )

        print(
            f"r={r:2d} "
            f"V={V[r]} "
            f"closed={closed} "
            f"V-(W-NW[-2])={residual_1} "
            f"V-closed={residual_2}"
        )

        if residual_1 != 0 or residual_2 != 0:
            green_failures += 1

    # ------------------------------------------------------------------------
    # 2. Source weights
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SOURCE-WEIGHT MATCHING")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        print()
        print(
            f"({k},{ell})"
        )

        weights = antisymmetric_source_weights(
            k,
            ell
        )

        for m in sorted(weights):

            actual = weights[m]

            residual, r = source_weight_residual(
                k,
                m,
                actual
            )

            if residual is None:
                continue

            expected = expected_source_weight_t(r)

            if residual != 0:
                source_failures += 1

                print(
                    f"  C_{m}: "
                    f"r={r} FAIL"
                )

                print(
                    "      actual   =",
                    actual
                )

                print(
                    "      expected =",
                    expected
                )

                print(
                    "      residual =",
                    residual
                )

            else:

                print(
                    f"  C_{m}: "
                    f"r={r} PASS"
                )

    # ------------------------------------------------------------------------
    # 3. D convolution
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. D FROM CORRECTED V CONVOLUTION")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k + 1:
            continue

        direct = D_direct(
            k,
            ell
        )

        source = D_from_V_source(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source
            )
        )

        boundary = sp.factor(
            direct.subs(
                N,
                0
            )
        )

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  D_direct =",
            direct
        )

        print(
            "  D_Vsource =",
            source
        )

        print(
            "  residual =",
            residual
        )

        print(
            "  D(0) =",
            boundary
        )

        if residual != 0:
            D_failures += 1

        if boundary != 1:
            boundary_failures += 1

        print(
            "  status =",
            "PASS"
            if residual == 0 and boundary == 1
            else "FAIL"
        )

    # ------------------------------------------------------------------------
    # 4. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    for k, ell in FORWARD:

        direct = D_direct(
            k,
            ell
        )

        source = D_from_V_source(
            k,
            ell
        )

        residual = sp.factor(
            sp.expand(
                direct - source
            )
        )

        boundary = sp.factor(
            direct.subs(
                N,
                0
            )
        )

        ok = (
            residual == 0
            and boundary == 1
        )

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
    # 5. Representative convolution
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. REPRESENTATIVE FINITE CONVOLUTION")
    print("=" * 78)

    for k, ell in [
        (1,9),
        (3,11),
        (5,13),
    ]:

        weights = antisymmetric_source_weights(
            k,
            ell
        )

        print()
        print(
            f"({k},{ell})"
        )

        for m in sorted(weights):

            r = m - 2*k

            print(
                f"  V_{r}(N) / t^{r+1} "
                f"* C_{m}"
            )

            print(
                "      V_r =",
                V_closed(r)
            )

            print(
                "      C_m =",
                kernel_mode_closed(
                    k,
                    ell,
                    m
                )
            )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "V kernel failures =",
        green_failures
    )

    print(
        "source-weight failures =",
        source_failures
    )

    print(
        "D convolution failures =",
        D_failures
    )

    print(
        "D(0) failures =",
        boundary_failures
    )

    print(
        "forward failures =",
        holdout_failures
    )

    if (
        green_failures == 0
        and source_failures == 0
        and D_failures == 0
        and boundary_failures == 0
        and holdout_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The t/N normalization is now handled in the same"
        )

        print(
            "variable domain before source comparison."
        )

        print()
        print(
            "The antisymmetric Green kernel exactly generates"
        )

        print(
            "the finite source weights:"
        )

        print(
            "    weight_r = V_r(t^2) / t^(r+1)."
        )

        print()
        print(
            "The same corrected kernel reproduces D_(k,ell)(N)"
        )

        print(
            "exactly on training and forward holdout."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "collapse the exact V/C convolution into a single"
        )

        print(
            "finite binomial or terminating hypergeometric sum."
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

