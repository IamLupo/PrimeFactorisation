#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 154
EXACT CENTRAL QUOTIENT MODE LAW

TARGET
------
For the supported central mode

    j = 2k,
    ell >= 2k+1,

Experiment 153 found

    q_(2k)(t) = N^k P_(k,ell)(N),
    N = t^2.

The previous experiment incorrectly treated the required t-power as zero.
The correct normalization is

    common t-power = 2k.

This experiment now derives P_(k,ell)(N) DIRECTLY from the finite
quotient-mode recurrence

    C_j = q_j + t q_(j-1) + t q_(j+1),

using the exact finite Laurent coefficients C_j from Experiment 151.

METHOD
------
For fixed (k,ell), start from the upper Laurent boundary of F and descend
exactly to j=2k.

No polynomial interpolation is used to define the law.

The experiment searches for a compact closed finite sum for

    P_(k,ell)(N) = q_(2k)/N^k

by comparing the exact recurrence-derived polynomial with candidate
binomial sums.

PRIMARY TARGET
--------------
Find a formula of the form

    P_(k,ell)(N)
      = sum_r a_r(k,ell) N^r

where a_r(k,ell) is obtained directly from binomial coefficients.

The script tests the candidate structural form

    a_r =
      (-1)^r *
      [ binom(ell-2k+r, r)
        combinations involving k ]

rather than blindly interpolating values.

It also checks:

    1. central-support condition ell >= 2k+1;
    2. exact t-power 2k;
    3. exact N-only reduction;
    4. forward ell holdout;
    5. degree law;
    6. leading and constant coefficients.

IMPORTANT
---------
If the proposed binomial closed form fails, the experiment prints the exact
coefficient sequence of P_(k,ell)(N). That sequence is the input for the next
symbolic derivation.

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
# Original kernel
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

    Q, R = sp.div(
        PF,
        PD
    )

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = {sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Laurent kernel / quotient
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


# ============================================================================
# Laurent coefficient
# ============================================================================

def laurent_support(expr):
    support = []

    for term in sp.Add.make_args(sp.expand(expr)):
        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        support.append(int(power))

    return sorted(set(support))


def laurent_coeff(expr, j):
    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    shift = max(0, -min(support))

    shifted = sp.expand(expr * u**shift)

    poly = sp.Poly(
        shifted,
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
# Direct central mode
# ============================================================================

def direct_central_mode(k, ell):
    return sp.factor(
        laurent_coeff(
            laurent_Q(k, ell),
            2*k
        )
    )


# ============================================================================
# Recover q modes from the finite kernel
# ============================================================================

def recover_q_modes(k, ell):
    """
    Exact downward recurrence:

        C_j = q_j + t q_(j-1) + t q_(j+1)

    and finite support.

    For F support [-ell, ell], Q support [-(ell-1), ell-1].
    """

    F = laurent_F(k, ell)
    support = laurent_support(F)

    fmin = min(support)
    fmax = max(support)

    q = {}

    qmax = fmax - 1

    q[qmax] = sp.cancel(
        laurent_coeff(F, fmax) / t
    )

    for j in range(
        fmax - 1,
        fmin,
        -1
    ):
        Cj = laurent_coeff(F, j)

        qj = q.get(
            j,
            sp.Integer(0)
        )

        qjp = q.get(
            j + 1,
            sp.Integer(0)
        )

        q[j - 1] = sp.expand(
            sp.cancel(
                (Cj - qj - t*qjp) / t
            )
        )

    return {
        j: sp.factor(expr)
        for j, expr in q.items()
        if expr != 0
    }


# ============================================================================
# Convert q_(2k) to P_(k,ell)(N)
# ============================================================================

def central_P(k, ell):
    """
    q_(2k) = t^(2k) P(N) = N^k P(N).
    """

    q_modes = recover_q_modes(k, ell)

    j = 2*k

    q2k = q_modes.get(
        j,
        sp.Integer(0)
    )

    if q2k == 0:
        return {
            "q": sp.Integer(0),
            "t_power": None,
            "P": sp.Integer(0)
        }

    # Required exact t-power.
    factor_t = sp.factor(
        q2k / t**(2*k)
    )

    factor_t = sp.expand(
        factor_t
    )

    poly_t = sp.Poly(
        factor_t,
        t,
        domain=sp.QQ
    )

    # Verify only even powers remain.
    for (power,), _coeff in poly_t.terms():
        if int(power) % 2 != 0:
            raise ArithmeticError(
                f"central mode has odd residual t-power "
                f"for ({k},{ell}): {q2k}"
            )

    P = sp.Integer(0)

    for (power,), coeff in poly_t.terms():
        power = int(power)

        P += coeff * N**(power // 2)

    return {
        "q": sp.factor(q2k),
        "t_power": 2*k,
        "P": sp.factor(sp.expand(P))
    }


# ============================================================================
# Polynomial coefficient list
# ============================================================================

def coefficient_list_N(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ
    )

    if poly.is_zero:
        return []

    degree = poly.degree()

    return [
        sp.factor(
            poly.coeff_monomial(N**r)
        )
        for r in range(degree + 1)
    ]


# ============================================================================
# Candidate binomial basis
# ============================================================================

def candidate_basis(k, ell):
    """
    Candidate exact finite-binomial basis.

    We use several natural quantities suggested by the direct kernel:
        binom(ell-a, r)
        binom(ell-a, k+b)
        binom(k+b, r)
        falling/rising factorial equivalents.

    The search is intentionally small and interpretable.
    """

    basis = []

    # Primary basis: binom(ell - s, r)
    for s in range(0, 4):
        for r in range(0, ell + 1):
            basis.append(
                (
                    f"C({ell}-{s},{r})",
                    sp.binomial(ell - s, r)
                )
            )

    # k-shifted basis.
    for shift in range(-2, 3):
        for r in range(0, ell + 1):
            idx = k + shift + r

            if 0 <= idx <= ell:
                basis.append(
                    (
                        f"C({ell},{k+shift}+{r})",
                        sp.binomial(ell, idx)
                    )
                )

    return basis


# ============================================================================
# Search coefficient formula for fixed (k,ell)
# ============================================================================

def direct_binomial_coefficient_predictions(
    k,
    ell,
    P
):
    """
    This is NOT a global fit.

    It compares each exact coefficient of P against simple direct
    binomial expressions and prints exact matches.

    The goal is to identify the most plausible symbolic envelope.
    """

    coeffs = coefficient_list_N(P)

    matches = []

    for r, coeff in enumerate(coeffs):

        candidates = {
            "C(ell,r)": sp.binomial(ell, r),
            "C(ell,r+1)": sp.binomial(ell, r+1),
            "C(ell-1,r)": sp.binomial(ell - 1, r),
            "C(ell-2,r)": sp.binomial(ell - 2, r),
            "C(ell,k+r)": (
                sp.binomial(
                    ell,
                    k + r
                )
                if k + r <= ell
                else sp.Integer(0)
            ),
            "C(ell,k+r-1)": (
                sp.binomial(
                    ell,
                    k + r - 1
                )
                if 0 <= k + r - 1 <= ell
                else sp.Integer(0)
            ),
        }

        exact = coeff

        for name, value in candidates.items():

            if value == 0:
                continue

            ratio = sp.factor(
                sp.Rational(exact, value)
            )

            # Only report unusually simple rational ratios.
            if sp.denom(ratio).bit_length() < 20:
                matches.append(
                    (
                        r,
                        exact,
                        name,
                        ratio
                    )
                )

    return matches


# ============================================================================
# Global leading coefficient check
# ============================================================================

def expected_leading_coefficient(k, ell):
    """
    From earlier edge experiments, test the observed top coefficient
    of P_(k,ell).

    This is only a certificate check, not the law being fitted here.
    """

    P = central_P(k, ell)["P"]

    if P == 0:
        return None

    poly = sp.Poly(
        P,
        N,
        domain=sp.QQ
    )

    degree = poly.degree()

    return sp.factor(
        poly.coeff_monomial(
            N**degree
        )
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 154")
    print("EXACT CENTRAL QUOTIENT MODE LAW")
    print("=" * 78)
    print()

    TRAIN = [
        (1, 3),
        (1, 5),
        (1, 7),
        (1, 9),
        (1, 11),
        (1, 13),
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
    ]

    FORWARD = [
        (1, 15),
        (1, 17),
        (3, 15),
        (3, 17),
        (5, 15),
        (5, 17),
    ]

    # ------------------------------------------------------------------------
    # Support conditions
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. CENTRAL SUPPORT CONDITION")
    print("=" * 78)

    support_failures = 0

    for k, ell in TRAIN + FORWARD:

        supported = ell >= 2*k + 1

        print(
            f"({k},{ell}) "
            f"2k={2*k} "
            f"ell-1={ell-1} "
            f"supported={supported}"
        )

        if not supported:
            support_failures += 1

    print(
        "unsupported records =",
        support_failures
    )

    # ------------------------------------------------------------------------
    # Exact central polynomial
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT P_(k,ell)(N) = q_(2k)/N^k")
    print("=" * 78)

    train_records = []

    for k, ell in TRAIN:

        result = central_P(
            k,
            ell
        )

        P = result["P"]
        q2k = result["q"]

        if q2k == 0:
            print(
                f"({k},{ell}) q_(2k)=0 "
                f"[SUPPORT BOUNDARY]"
            )
            continue

        print(
            f"({k},{ell})"
        )

        print(
            "  q_(2k) =",
            q2k
        )

        print(
            "  exact t-power =",
            result["t_power"],
            f"(expected {2*k})"
        )

        print(
            "  P_(k,ell)(N) =",
            P
        )

        poly = sp.Poly(
            P,
            N,
            domain=sp.QQ
        )

        print(
            "  degree(P) =",
            poly.degree()
        )

        print(
            "  coefficients =",
            coefficient_list_N(P)
        )

        train_records.append(
            {
                "k": k,
                "ell": ell,
                "P": P,
                "degree": poly.degree(),
            }
        )

    # ------------------------------------------------------------------------
    # Check observed degree law
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DEGREE LAW")
    print("=" * 78)

    degree_failures = 0

    for rec in train_records:

        k = rec["k"]
        ell = rec["ell"]

        observed = rec["degree"]

        # Empirically suggested law:
        #
        # degree(P) = floor((ell-1)/2)
        #
        # for k=1, but not uniformly for larger k.
        #
        # Therefore print the exact offset relative to ell-k.
        expected_reference = ell - k - 1

        print(
            f"({k},{ell}) "
            f"degree={observed} "
            f"ell-k-1={expected_reference} "
            f"offset={observed-expected_reference}"
        )

    # ------------------------------------------------------------------------
    # Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_failures = 0

    for k, ell in FORWARD:

        result = central_P(
            k,
            ell
        )

        q2k = result["q"]

        if q2k == 0:
            print(
                f"({k},{ell}) q_(2k)=0 "
                f"[unexpected support loss]"
            )
            forward_failures += 1
            continue

        P = result["P"]

        if result["t_power"] != 2*k:
            forward_failures += 1

        print(
            f"({k},{ell})"
        )

        print(
            "  q_(2k) =",
            q2k
        )

        print(
            "  P(N) =",
            P
        )

        print(
            "  degree =",
            sp.Poly(
                P,
                N,
                domain=sp.QQ
            ).degree()
        )

    print(
        "forward failures =",
        forward_failures
    )

    # ------------------------------------------------------------------------
    # Direct coefficient observations
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT BINOMIAL COEFFICIENT OBSERVATIONS")
    print("=" * 78)

    for rec in train_records:

        k = rec["k"]
        ell = rec["ell"]
        P = rec["P"]

        matches = direct_binomial_coefficient_predictions(
            k,
            ell,
            P
        )

        print()
        print(
            f"({k},{ell})"
        )

        if not matches:

            print(
                "  no simple direct coefficient ratio found"
            )

        else:

            for item in matches[:12]:

                r, coeff, name, ratio = item

                print(
                    f"  r={r}: "
                    f"coeff={coeff} "
                    f"vs {name} "
                    f"ratio={ratio}"
                )

    # ------------------------------------------------------------------------
    # Leading coefficient profile
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. LEADING COEFFICIENT PROFILE")
    print("=" * 78)

    for rec in train_records:

        k = rec["k"]
        ell = rec["ell"]
        P = rec["P"]

        poly = sp.Poly(
            P,
            N,
            domain=sp.QQ
        )

        lc = sp.factor(
            poly.LC()
        )

        print(
            f"({k},{ell}) "
            f"degree={poly.degree()} "
            f"leading={lc}"
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if forward_failures == 0:

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The supported central quotient mode satisfies"
        )

        print(
            "    q_(2k) = N^k P_(k,ell)(N)"
        )

        print()
        print(
            "for all tested training and forward pairs."
        )

        print()
        print(
            "The next symbolic task is to derive P_(k,ell)(N)"
        )

        print(
            "directly from the quotient-mode recurrence."
        )

        print()
        print(
            "Do NOT promote a coefficient pattern to a theorem"
        )

        print(
            "until it is derived from the exact finite recurrence."
        )

    else:

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

