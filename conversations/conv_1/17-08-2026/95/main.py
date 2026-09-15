#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 159
CORRECTED ANTISYMMETRIC GREEN KERNEL V_r
EXACT D-CONVOLUTION CERTIFICATE

Experiment 158R established:

    W_0 = 1
    W_1 = -1
    W_(r+1) = -W_r - N W_(r-1)

with

    W_r =
      sum_{m=0}^{floor(r/2)}
        (-1)^(r-m) binom(r-m,m) N^m.

The actual antisymmetric source weights are NOT W_r directly.

They satisfy

    V_r = W_r - N W_(r-2),

with

    V_0 = 1
    V_1 = -1
    V_2 = 1 - 2N

and the same second-order operator thereafter:

    V_(r+1) = -V_r - N V_(r-1).

Equivalent closed form:

    V_r =
      sum_{m=0}^{floor(r/2)}
        (-1)^(r-m)
        [
            binom(r-m,m)
          + binom(r-m-1,m-1)
        ]
        N^m,

where the second binomial is interpreted as zero for m=0.

TARGET
------
1. Verify V recurrence.
2. Verify W_r - N W_(r-2) exactly.
3. Verify the closed V finite-binomial formula.
4. Match every antisymmetric source weight from the quotient recurrence.
5. Use V_r to reconstruct D_(k,ell)(N) directly from the finite kernel.
6. Compare against the independently constructed quotient.
7. Test forward ell holdout.
8. Inspect whether the D convolution can be compressed by a terminating
   Vandermonde/Chu-Vandermonde-type summation.

This is NOT polynomial fitting.

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
# Exact kernel / quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q_pq(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q)
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q)
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
                q: t/u
            }
        )
    )


def laurent_Q(k, ell):
    return sp.expand(
        quotient_Q_pq(k, ell).subs(
            {
                p: t*u,
                q: t/u
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
                f"noninteger Laurent exponent in u: {term}"
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
        domain=sp.QQ.frac_field(t)
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(u**target)
    )


# ============================================================================
# Laurent-t normalization
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
                f"noninteger t exponent: {term}"
            )

        exponent = int(exponent)

        if minimum is None or exponent < minimum:
            minimum = exponent

    return minimum, sp.expand(
        expr / t**minimum
    )


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

    return sp.factor(sp.expand(result))


# ============================================================================
# Correct W sequence
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
        1: sp.Integer(-1)
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
# Correct antisymmetric V sequence
# ============================================================================

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


def V_recursive(max_r):
    V = {
        0: sp.Integer(1),
        1: sp.Integer(-1)
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


def V_closed(r):
    result = sp.Integer(0)

    for m in range(
        0,
        r // 2 + 1
    ):
        a = sp.binomial(
            r-m,
            m
        )

        if m == 0:
            b = sp.Integer(0)
        else:
            b = sp.binomial(
                r-m-1,
                m-1
            )

        result += (
            (-1)**(r-m)
            * (a + b)
            * N**m
        )

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Source weights of q_j
# ============================================================================

def source_weights_for_q(
    k,
    ell,
    target_j
):
    source = {}

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

        source[j-1] = {
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
# Match a source weight against V_r
# ============================================================================

def source_weight_match(k, m, coeff):
    r = m - 2*k

    if r < 0:
        return {
            "r": r,
            "valid": False,
            "reason": "negative r"
        }

    # Actual expected source coefficient:
    #
    #     V_r(N) / t^(r+1)
    #
    expected = sp.cancel(
        V_closed(r)
        / t**(r+1)
    )

    residual = sp.factor(
        sp.expand(
            coeff - expected
        )
    )

    return {
        "r": r,
        "expected": expected,
        "residual": residual,
        "valid": residual == 0
    }


# ============================================================================
# Exact finite kernel coefficient
# ============================================================================

def kernel_mode_closed(k, ell, j):
    # p^k(1+q)^ell
    r1 = k-j

    if 0 <= r1 <= ell:
        c1 = (
            sp.binomial(ell, r1)
            * t**(k+r1)
        )
    else:
        c1 = sp.Integer(0)

    # q^k(1+p)^ell
    r2 = k+j

    if 0 <= r2 <= ell:
        c2 = (
            sp.binomial(ell, r2)
            * t**(k+r2)
        )
    else:
        c2 = sp.Integer(0)

    # -p^ell(1+q)^k
    r3 = ell-j

    if 0 <= r3 <= k:
        c3 = (
            -sp.binomial(k, r3)
            * t**(ell+r3)
        )
    else:
        c3 = sp.Integer(0)

    # -q^ell(1+p)^k
    r4 = ell+j

    if 0 <= r4 <= k:
        c4 = (
            -sp.binomial(k, r4)
            * t**(ell+r4)
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

    for m, coeff in weights.items():
        r = m - 2*k

        if r < 0:
            raise ArithmeticError(
                f"unexpected negative Green index {r}"
            )

        # Replace the actual weight by the corrected closed V weight.
        corrected_weight = sp.cancel(
            V_closed(r)
            / t**(r+1)
        )

        # Exact source must match before summing.
        if sp.factor(
            sp.expand(
                coeff - corrected_weight
            )
        ) != 0:
            raise ArithmeticError(
                f"source weight mismatch at m={m}, "
                f"r={r}, pair=({k},{ell})"
            )

        raw += (
            corrected_weight
            * kernel_mode_closed(
                k,
                ell,
                m
            )
        )

    raw = sp.cancel(
        raw
        / t**(2*k - 1)
    )

    base, reduced = split_t_laurent(raw)

    if base != 0:
        raise ArithmeticError(
            f"D retained t power {base} for ({k},{ell})"
        )

    result = even_t_to_N(reduced)

    if result is None:
        raise ArithmeticError(
            f"D retained odd t powers for ({k},{ell})"
        )

    return sp.factor(result)


# ============================================================================
# Direct D
# ============================================================================

def D_direct(k, ell):
    Q = laurent_Q(k, ell)

    qm = laurent_coeff(
        Q,
        2*k-1
    )

    qp = laurent_coeff(
        Q,
        2*k+1
    )

    raw = sp.cancel(
        (qm-qp)
        / t**(2*k-1)
    )

    base, reduced = split_t_laurent(raw)

    if base != 0:
        raise ArithmeticError(
            f"D direct retained t power {base}"
        )

    result = even_t_to_N(reduced)

    if result is None:
        raise ArithmeticError(
            f"D direct retained odd t powers"
        )

    return sp.factor(result)


# ============================================================================
# Convolution coefficient inspection
# ============================================================================

def polynomial_coefficients(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ
    )

    if poly.is_zero:
        return []

    return [
        sp.factor(
            poly.coeff_monomial(N**r)
        )
        for r in range(
            poly.degree()+1
        )
    ]


# ============================================================================
# Search for convolution compression
# ============================================================================

def vandermonde_diagnostic(k, ell, D):
    """
    Compare D coefficients against basic binomial transforms.

    This is only a symbolic diagnostic. It does not define the answer by fit.
    """

    coeffs = polynomial_coefficients(D)

    observations = []

    for r, coeff in enumerate(coeffs):

        candidates = {
            "C(ell,r)": sp.binomial(ell, r),
            "C(ell-1,r)": sp.binomial(ell-1, r),
            "C(ell-2,r)": sp.binomial(ell-2, r),
            "C(ell,k+r)": (
                sp.binomial(ell, k+r)
                if 0 <= k+r <= ell
                else sp.Integer(0)
            ),
            "C(ell,k+r-1)": (
                sp.binomial(ell, k+r-1)
                if 0 <= k+r-1 <= ell
                else sp.Integer(0)
            )
        }

        for label, value in candidates.items():

            if value == 0:
                continue

            ratio = sp.factor(
                sp.cancel(
                    sp.Rational(coeff, value)
                )
            )

            observations.append(
                (
                    r,
                    coeff,
                    label,
                    ratio
                )
            )

    return observations


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 159")
    print("CORRECTED ANTISYMMETRIC GREEN KERNEL")
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

    green_failures = 0
    source_failures = 0
    D_failures = 0
    boundary_failures = 0
    holdout_failures = 0

    # ------------------------------------------------------------------------
    # 1. V recurrence + closed form
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. CORRECTED V GREEN POLYNOMIAL")
    print("=" * 78)

    W = W_recursive(20)
    V = V_recursive(20)

    for r in range(21):

        from_W = V_from_W(
            r,
            W
        )

        recurrence_V = V[r]
        closed_V = V_closed(r)

        residual_1 = sp.factor(
            sp.expand(
                recurrence_V - from_W
            )
        )

        residual_2 = sp.factor(
            sp.expand(
                recurrence_V - closed_V
            )
        )

        print(
            f"r={r:2d} "
            f"V={recurrence_V} "
            f"closed={closed_V} "
            f"V-(W-NW[-2])={residual_1} "
            f"V-closed={residual_2}"
        )

        if residual_1 != 0 or residual_2 != 0:
            green_failures += 1

    # ------------------------------------------------------------------------
    # 2. Source matching
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ACTUAL SOURCE WEIGHTS VS V_r")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k+1:
            continue

        weights = antisymmetric_source_weights(
            k,
            ell
        )

        print()
        print(
            f"({k},{ell})"
        )

        for m in sorted(weights):

            check = source_weight_match(
                k,
                m,
                weights[m]
            )

            r = check["r"]

            if not check["valid"]:

                source_failures += 1

                print(
                    f"  C_{m}: FAIL "
                    f"r={r} "
                    f"actual={weights[m]} "
                    f"expected={check['expected']} "
                    f"residual={check['residual']}"
                )

            else:

                print(
                    f"  C_{m}: "
                    f"r={r} "
                    f"weight={weights[m]} "
                    f"PASS"
                )

    # ------------------------------------------------------------------------
    # 3. Direct D vs V-source D
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. D FROM CORRECTED V SOURCE")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k+1:
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
                direct-source
            )
        )

        boundary = sp.factor(
            direct.subs(N, 0)
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

        if residual == 0 and boundary == 1:
            print(
                "  status = PASS"
            )
        else:
            print(
                "  status = FAIL"
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
                direct-source
            )
        )

        ok = (
            residual == 0
            and sp.factor(
                direct.subs(N, 0)
            ) == 1
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
    # 5. Convolution diagnostics
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. D COEFFICIENT / VANDERMONDE DIAGNOSTIC")
    print("=" * 78)

    for k, ell in TRAIN:

        if ell < 2*k+1:
            continue

        D = D_from_V_source(
            k,
            ell
        )

        print()
        print(
            f"({k},{ell}) D coefficients =",
            polynomial_coefficients(D)
        )

        observations = vandermonde_diagnostic(
            k,
            ell,
            D
        )

        # Print only the first few exact ratios so the output remains
        # manageable.
        for obs in observations[:8]:

            r, coeff, label, ratio = obs

            print(
                f"  r={r}: "
                f"{coeff} / {label} = {ratio}"
            )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "V recurrence/closed-form failures =",
        green_failures
    )

    print(
        "source-weight failures =",
        source_failures
    )

    print(
        "D source residual failures =",
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
            "The corrected antisymmetric Green kernel is exact:"
        )

        print(
            "    V_r = W_r - N W_(r-2)"
        )

        print(
            "and"
        )

        print(
            "    V_r = sum_m (-1)^(r-m)"
        )

        print(
            "          [C(r-m,m)+C(r-m-1,m-1)] N^m."
        )

        print()
        print(
            "The same V kernel reproduces D_(k,ell)(N)"
        )

        print(
            "exactly on training and forward ell holdouts."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "compress the resulting finite convolution for D"
        )

        print(
            "using an exact Vandermonde/Chu-Vandermonde identity."
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

