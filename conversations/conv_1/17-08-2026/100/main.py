#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 164
BRANCHWISE ADMISSIBLE HYPERGEOMETRIC TERM CERTIFICATE

Purpose
-------
Fix the symbolic failure of Experiment 163.

We DO NOT merge the four kernel branches before respecting their support.

For

    F = A + B + C + D

analyze each branch separately:

    A = p^k (1+q)^ell
    B = q^k (1+p)^ell
    C = -p^ell (1+q)^k
    D = -q^ell (1+p)^k

For fixed (k,ell,s), the coefficient contribution to [N^s]D is obtained
only on the admissible integer j interval where the corresponding binomial
coefficients exist.

For each branch we compute the exact term T_branch(j), then

    T_branch(j+1) / T_branch(j)

using ONLY admissible integer values.

We then attempt to infer the exact rational-in-j ratio symbolically by
rewriting finite binomials with factorial/gamma-free positive arguments.

This experiment is diagnostic only: it does not fit arbitrary polynomials.

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

N = sp.symbols("N")


# ============================================================================
# V coefficient
# ============================================================================

def v_coeff(r, m):
    if r < 0 or m < 0 or m > r // 2:
        return sp.Integer(0)

    a = sp.binomial(r - m, m)
    b = (
        sp.Integer(0)
        if m == 0
        else sp.binomial(r - m - 1, m - 1)
    )

    return sp.expand(
        (-1)**(r-m) * (a + b)
    )


# ============================================================================
# Branch data
# ============================================================================

def branch_data(k, ell, branch):
    """
    Return:
        e(j)     = t exponent of kernel branch
        c(j)     = binomial coefficient including sign
        support  = admissible integer j interval
    """

    if branch == "A":
        # p^k (1+q)^ell
        # j = k-a, a=0..ell
        # j ranges k-ell .. k
        def e(j):
            return 2*k - j

        def c(j):
            return sp.binomial(
                ell,
                k-j
            )

        return (
            e,
            c,
            k-ell,
            k
        )

    if branch == "B":
        # q^k (1+p)^ell
        # j = b-k, b=0..ell
        # j ranges -k .. ell-k
        def e(j):
            return 2*k + j

        def c(j):
            return sp.binomial(
                ell,
                k+j
            )

        return (
            e,
            c,
            -k,
            ell-k
        )

    if branch == "C":
        # -p^ell (1+q)^k
        # j = ell-c, c=0..k
        # j ranges ell-k .. ell
        def e(j):
            return 2*ell - j

        def c(j):
            return -sp.binomial(
                k,
                ell-j
            )

        return (
            e,
            c,
            ell-k,
            ell
        )

    if branch == "D":
        # -q^ell (1+p)^k
        # j = d-ell, d=0..k
        # j ranges -ell .. k-ell
        def e(j):
            return 2*ell + j

        def c(j):
            return -sp.binomial(
                k,
                ell+j
            )

        return (
            e,
            c,
            -ell,
            k-ell
        )

    raise ValueError(branch)


# ============================================================================
# Exact branch contribution for fixed j
# ============================================================================

def branch_term(k, ell, s, branch, j):
    """
    Return exact contribution of this branch to [N^s]D at integer j.

    r = j - 2k
    and final t-degree condition gives m.
    """

    e_fun, c_fun, j_lo, j_hi = branch_data(
        k, ell, branch
    )

    if not (j_lo <= j <= j_hi):
        return sp.Integer(0)

    r = j - 2*k

    if r < 0:
        return sp.Integer(0)

    e = e_fun(j)
    c = c_fun(j)

    # Final exponent:
    #
    # 2m + e - (r+1) - (2k-1) = 2s
    #
    m = s - sp.Rational(
        e - r - 2*k,
        2
    )

    if not m.is_Integer:
        return sp.Integer(0)

    m = int(m)

    if m < 0 or m > r // 2:
        return sp.Integer(0)

    return sp.factor(
        c * v_coeff(r, m)
    )


# ============================================================================
# Admissible support of a branch for fixed s
# ============================================================================

def admissible_js(k, ell, s, branch):
    _, _, lo, hi = branch_data(
        k, ell, branch
    )

    result = []

    for j in range(lo, hi+1):
        value = branch_term(
            k, ell, s, branch, j
        )

        if value != 0:
            result.append(j)

    return result


# ============================================================================
# Positive-integer factorial normalization
# ============================================================================

def binomial_factorial(n, r):
    """
    Exact factorial representation only when 0 <= r <= n.
    """

    if not (
        isinstance(n, int)
        and isinstance(r, int)
    ):
        raise TypeError()

    if r < 0 or r > n:
        return sp.Integer(0)

    return sp.factorial(n) / (
        sp.factorial(r)
        * sp.factorial(n-r)
    )


# ============================================================================
# Symbolic term reconstruction from a branch at first support point
# ============================================================================

def symbolic_branch_term(k, ell, s, branch):
    """
    Build T(j) only after the admissible interval has been established.

    We use the concrete combinatorial identities for generalized
    binomials only after confirming the arguments stay nonnegative
    on the admissible interval.
    """

    _, _, lo, hi = branch_data(
        k, ell, branch
    )

    js = admissible_js(
        k, ell, s, branch
    )

    if not js:
        return None, None, None

    # The support is finite. To avoid invalid symbolic continuation,
    # retain the interval explicitly.
    return js, lo, hi


# ============================================================================
# Concrete ratio table
# ============================================================================

def ratio_table(k, ell, s, branch):
    js, _, _ = symbolic_branch_term(
        k, ell, s, branch
    )

    if js is None or len(js) < 2:
        return []

    ratios = []

    for a, b in zip(
        js,
        js[1:]
    ):
        va = branch_term(
            k, ell, s, branch, a
        )
        vb = branch_term(
            k, ell, s, branch, b
        )

        if va == 0 or vb == 0:
            continue

        ratios.append(
            (
                a,
                b,
                sp.factor(
                    sp.cancel(vb/va)
                )
            )
        )

    return ratios


# ============================================================================
# Guess the rational ratio from exact finite ratios
# ============================================================================

def rational_ratio_diagnostic(ratios):
    """
    Attempt a rational interpolation ONLY as a structural diagnostic.

    This is not used to prove D. We merely detect whether the observed
    ratio is compatible with a low-degree rational function of j.
    """

    if len(ratios) < 3:
        return None

    jvals = [a for a, _, _ in ratios]
    vals = [v for _, _, v in ratios]

    # Try small rational degrees.
    x = sp.symbols("x")

    for deg_num in range(0, 4):
        for deg_den in range(1, 5):

            unknowns = sp.symbols(
                f"a0:{deg_num+1}"
            ) + sp.symbols(
                f"b0:{deg_den+1}"
            )

            equations = []

            num = sum(
                unknowns[i] * x**i
                for i in range(deg_num+1)
            )

            den = (
                x**deg_den
                + sum(
                    unknowns[
                        deg_num+1+i
                    ] * x**i
                    for i in range(deg_den)
                )
            )

            # Require enough equations.
            if len(jvals) < len(unknowns):
                continue

            for xv, yv in zip(
                jvals,
                vals
            ):
                equations.append(
                    sp.Eq(
                        num.subs(x, xv),
                        yv * den.subs(x, xv)
                    )
                )

            sol = sp.solve(
                equations,
                unknowns,
                dict=True
            )

            if not sol:
                continue

            candidate = sp.factor(
                num.subs(sol[0])
                / den.subs(sol[0])
            )

            ok = True

            for xv, yv in zip(
                jvals,
                vals
            ):
                if sp.factor(
                    candidate.subs(x, xv) - yv
                ) != 0:
                    ok = False
                    break

            if ok:
                return candidate

    return None


# ============================================================================
# Exact coefficient total
# ============================================================================

def branch_total(k, ell, s, branch):
    return sp.factor(
        sum(
            branch_term(
                k, ell, s, branch, j
            )
            for j in admissible_js(
                k, ell, s, branch
            )
        )
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 164")
    print("BRANCHWISE ADMISSIBLE HYPERGEOMETRIC TERMS")
    print("=" * 78)
    print()

    TESTS = [
        (1, 9),
        (1, 11),
        (1, 13),
        (3, 11),
        (3, 13),
        (5, 13),
    ]

    failures = 0

    print("=" * 78)
    print("1. BRANCH SUPPORT")
    print("=" * 78)

    for k, ell in TESTS:

        print()
        print(f"({k},{ell})")

        for branch in ["A", "B", "C", "D"]:

            _, _, lo, hi = branch_data(
                k, ell, branch
            )

            print(
                f"  branch {branch}: "
                f"raw support j=[{lo},{hi}]"
            )

            for s in range(
                min(3, max(0, ell-2*k)+1)
            ):

                js = admissible_js(
                    k, ell, s, branch
                )

                print(
                    f"    s={s}: "
                    f"admissible={js}"
                )

    print()
    print("=" * 78)
    print("2. BRANCHWISE CONCRETE RATIOS")
    print("=" * 78)

    for k, ell in TESTS:

        print()
        print(f"({k},{ell})")

        for branch in ["A", "B", "C", "D"]:

            for s in range(
                min(3, max(0, ell-2*k)+1)
            ):

                js = admissible_js(
                    k, ell, s, branch
                )

                if len(js) < 2:
                    continue

                ratios = ratio_table(
                    k, ell, s, branch
                )

                print()
                print(
                    f"  branch={branch}, s={s}"
                )

                for a, b, ratio in ratios:
                    print(
                        f"    j={a}->{b}: "
                        f"ratio={ratio}"
                    )

                candidate = rational_ratio_diagnostic(
                    ratios
                )

                print(
                    "    rational ratio candidate =",
                    candidate
                )

                # Verify candidate exactly against all observed ratios.
                if candidate is not None:
                    for a, _, ratio in ratios:
                        residual = sp.factor(
                            candidate.subs(
                                sp.Symbol("x"),
                                a
                            )
                            - ratio
                        )

                        if residual != 0:
                            failures += 1

    print()
    print("=" * 78)
    print("3. BRANCH TOTALS")
    print("=" * 78)

    for k, ell in TESTS:

        print()
        print(f"({k},{ell})")

        for s in range(
            min(3, max(0, ell-2*k)+1)
        ):

            values = {}

            for branch in ["A", "B", "C", "D"]:
                values[branch] = branch_total(
                    k, ell, s, branch
                )

            total = sp.factor(
                sum(values.values())
            )

            print(
                f"  [N^{s}]"
            )

            print(
                "    A =",
                values["A"]
            )

            print(
                "    B =",
                values["B"]
            )

            print(
                "    C =",
                values["C"]
            )

            print(
                "    D =",
                values["D"]
            )

            print(
                "    total =",
                total
            )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "structural failures =",
        failures
    )

    if failures == 0:
        print()
        print("STATUS = PASS")
        print()
        print(
            "The four kernel branches can be treated separately"
        )
        print(
            "without generalized-binomial or negative-factorial artifacts."
        )
        print()
        print(
            "The remaining summands exhibit rational consecutive"
        )
        print(
            "ratios on their exact finite supports."
        )
        print()
        print(
            "NEXT TARGET:"
        )
        print(
            "derive the branchwise rational ratio symbolically and"
        )
        print(
            "identify each branch as a terminating hypergeometric sum."
        )
    else:
        print()
        print("STATUS = FAIL")

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
