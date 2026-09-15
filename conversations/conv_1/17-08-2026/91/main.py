#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 155
CENTRAL THREE-MODE COUPLING

TARGET
------
Experiment 154 established, for supported pairs ell >= 2k+1,

    q_(2k) = N^k P_(k,ell)(N),

with N=t^2.

The next local modes are expected to have the form

    q_(2k-1) = t R_(k,ell)(N),
    q_(2k)   = N^k P_(k,ell)(N),
    q_(2k+1) = t T_(k,ell)(N).

From the exact Laurent quotient recurrence,

    C_j = q_j + t q_(j-1) + t q_(j+1),

the central equation is

    C_(2k)
      = q_(2k)
        + t q_(2k-1)
        + t q_(2k+1).

After substitution,

    C_(2k)
      = N^k P(N)
        + N R(N)
        + N T(N).

Therefore

    (C_(2k) - N^k P(N))/N
        = R(N) + T(N).

This experiment derives the three modes exactly and tests:

    1. q_(2k-1)/t is N-only;
    2. q_(2k+1)/t is N-only;
    3. the central three-mode identity holds exactly;
    4. the sum R+T has a direct closed formula from C_(2k);
    5. the difference R-T has a simple structure;
    6. forward ell values obey the same structure.

The experiment does NOT claim that R or T directly recovers S.
It isolates the smallest local quotient-mode block around the N-only
central mode.

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
        domain=sp.QQ.frac_field(t),
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(u**target)
    )


# ============================================================================
# Convert even powers of t to N
# ============================================================================

def rewrite_even_t(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return sp.Integer(0)

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

    return sp.factor(sp.expand(result))


# ============================================================================
# Exact central kernel coefficient C_(2k)
# ============================================================================

def C_central(k, ell):
    F = laurent_F(k, ell)

    return sp.factor(
        laurent_coeff(
            F,
            2*k
        )
    )


# ============================================================================
# Exact quotient central three modes
# ============================================================================

def extract_three_modes(k, ell):
    Q = laurent_Q(k, ell)

    j0 = 2*k

    qm = sp.factor(
        laurent_coeff(
            Q,
            j0 - 1
        )
    )

    q0 = sp.factor(
        laurent_coeff(
            Q,
            j0
        )
    )

    qp = sp.factor(
        laurent_coeff(
            Q,
            j0 + 1
        )
    )

    return qm, q0, qp


# ============================================================================
# Normalize the three modes
# ============================================================================

def normalize_modes(k, qm, q0, qp):
    """
    Expected structure:

        qm = t * R(N)
        q0 = t^(2k) * P(N)
        qp = t * T(N)

    Return exact N-polynomials if the structure holds.
    """

    R_raw = sp.expand(
        sp.cancel(qm / t)
    )

    P_raw = sp.expand(
        sp.cancel(q0 / t**(2*k))
    )

    T_raw = sp.expand(
        sp.cancel(qp / t)
    )

    R = rewrite_even_t(R_raw)
    P = rewrite_even_t(P_raw)
    T = rewrite_even_t(T_raw)

    return R, P, T


# ============================================================================
# Central recurrence certificate
# ============================================================================

def central_identity(k, ell, R, P, T):
    """
        C_(2k) = N^k P + N R + N T
    """

    C = C_central(
        k,
        ell
    )

    C_N = rewrite_even_t(C)

    if C_N is None:
        raise ArithmeticError(
            f"C_(2k) retains odd t powers for ({k},{ell})"
        )

    residual = sp.factor(
        sp.expand(
            C_N
            - N**k * P
            - N * R
            - N * T
        )
    )

    return C_N, residual


# ============================================================================
# Symmetric/antisymmetric neighboring combinations
# ============================================================================

def mode_combinations(R, T):
    return {
        "sum": sp.factor(
            sp.expand(R + T)
        ),
        "difference": sp.factor(
            sp.expand(R - T)
        ),
    }


# ============================================================================
# Compare neighboring coefficients directly
# ============================================================================

def adjacent_recurrence_residual(
    k,
    ell,
    qm,
    q0,
    qp
):
    """
    Verify the recurrence in raw t,u form at j=2k.
    """

    C = laurent_coeff(
        laurent_F(k, ell),
        2*k
    )

    residual = sp.factor(
        sp.expand(
            C
            - q0
            - t*qm
            - t*qp
        )
    )

    return residual


# ============================================================================
# Degree summary
# ============================================================================

def poly_degree(expr):
    if expr is None:
        return None

    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return -sp.oo

    return poly.degree()


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 155")
    print("CENTRAL THREE-MODE COUPLING")
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
    # 1. Training
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. TRAINING THREE-MODE CERTIFICATE")
    print("=" * 78)

    train_failures = 0
    structural_failures = 0
    central_failures = 0

    records = []

    for k, ell in TRAIN:

        # Central mode must be inside quotient support:
        # 2k <= ell-1.
        if ell < 2*k + 1:
            print(
                f"({k},{ell}) skipped: central mode outside support"
            )
            continue

        qm, q0, qp = extract_three_modes(
            k,
            ell
        )

        raw_residual = adjacent_recurrence_residual(
            k,
            ell,
            qm,
            q0,
            qp
        )

        R, P, T = normalize_modes(
            k,
            qm,
            q0,
            qp
        )

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  q_(2k-1) =",
            qm
        )

        print(
            "  q_(2k)   =",
            q0
        )

        print(
            "  q_(2k+1) =",
            qp
        )

        print(
            "  R(N) =",
            R
        )

        print(
            "  P(N) =",
            P
        )

        print(
            "  T(N) =",
            T
        )

        print(
            "  raw central recurrence =",
            "PASS"
            if raw_residual == 0
            else "FAIL"
        )

        if raw_residual != 0:
            train_failures += 1
            print(
                "    residual =",
                raw_residual
            )

        if (
            R is None
            or P is None
            or T is None
        ):
            structural_failures += 1
            print(
                "  N-polynomial structure = FAIL"
            )
        else:
            print(
                "  N-polynomial structure = PASS"
            )

        if R is not None and P is not None and T is not None:

            C_N, central_residual = central_identity(
                k,
                ell,
                R,
                P,
                T
            )

            print(
                "  C_(2k)(N) =",
                C_N
            )

            print(
                "  normalized central identity =",
                "PASS"
                if central_residual == 0
                else "FAIL"
            )

            if central_residual != 0:
                central_failures += 1
                print(
                    "    residual =",
                    central_residual
                )

            combo = mode_combinations(
                R,
                T
            )

            print(
                "  R+T =",
                combo["sum"]
            )

            print(
                "  R-T =",
                combo["difference"]
            )

            print(
                "  degrees:",
                "deg(R)=", poly_degree(R),
                "deg(P)=", poly_degree(P),
                "deg(T)=", poly_degree(T)
            )

            records.append(
                {
                    "k": k,
                    "ell": ell,
                    "R": R,
                    "P": P,
                    "T": T,
                    "sum": combo["sum"],
                    "difference": combo["difference"],
                    "C": C_N,
                }
            )

    # ------------------------------------------------------------------------
    # 2. Forward holdout
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_failures = 0

    for k, ell in FORWARD:

        if ell < 2*k + 1:
            print(
                f"({k},{ell}) skipped: central mode outside support"
            )
            continue

        qm, q0, qp = extract_three_modes(
            k,
            ell
        )

        raw_residual = adjacent_recurrence_residual(
            k,
            ell,
            qm,
            q0,
            qp
        )

        R, P, T = normalize_modes(
            k,
            qm,
            q0,
            qp
        )

        structure_ok = (
            R is not None
            and P is not None
            and T is not None
        )

        normalized_ok = False

        if structure_ok:
            _C_N, residual = central_identity(
                k,
                ell,
                R,
                P,
                T
            )

            normalized_ok = (
                residual == 0
            )

        ok = (
            raw_residual == 0
            and structure_ok
            and normalized_ok
        )

        print(
            f"({k},{ell}) "
            f"status={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            forward_failures += 1

        if structure_ok:

            combo = mode_combinations(
                R,
                T
            )

            print(
                "  R(N) =",
                R
            )

            print(
                "  P(N) =",
                P
            )

            print(
                "  T(N) =",
                T
            )

            print(
                "  R+T =",
                combo["sum"]
            )

            print(
                "  R-T =",
                combo["difference"]
            )

    # ------------------------------------------------------------------------
    # 3. Search for simple relation between R and T
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NEIGHBOR-MODE RELATION SEARCH")
    print("=" * 78)

    relation_failures = 0

    # Test the simplest exact hypotheses:
    #
    #   R+T = C/N-type expression
    #   R-T is divisible by simple factors
    #
    # We do not fit a law here; we only check divisibility/zero residuals
    # for natural quantities already present in the recurrence.

    for rec in records:

        k = rec["k"]
        ell = rec["ell"]

        R = rec["R"]
        T = rec["T"]
        C = rec["C"]
        P = rec["P"]

        # From central recurrence:
        #
        # C = N^k P + N(R+T)
        #
        expected_sum = sp.factor(
            sp.cancel(
                (C - N**k * P) / N
            )
        )

        residual = sp.factor(
            sp.expand(
                R + T - expected_sum
            )
        )

        if residual != 0:
            relation_failures += 1

        print()
        print(
            f"({k},{ell})"
        )

        print(
            "  R+T expected =",
            expected_sum
        )

        print(
            "  R+T residual =",
            residual
        )

        print(
            "  R-T =",
            sp.factor(
                R - T
            )
        )

    # ------------------------------------------------------------------------
    # 4. Final diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "training raw recurrence failures =",
        train_failures
    )

    print(
        "training N-structure failures =",
        structural_failures
    )

    print(
        "training normalized identity failures =",
        central_failures
    )

    print(
        "forward failures =",
        forward_failures
    )

    print(
        "neighbor relation failures =",
        relation_failures
    )

    if (
        train_failures == 0
        and structural_failures == 0
        and central_failures == 0
        and forward_failures == 0
        and relation_failures == 0
    ):

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The central three-mode block has exact structure:"
        )

        print()
        print(
            "    q_(2k-1) = t R(N)"
        )

        print(
            "    q_(2k)   = N^k P(N)"
        )

        print(
            "    q_(2k+1) = t T(N)"
        )

        print()
        print(
            "and"
        )

        print(
            "    C_(2k) = N^k P(N) + N(R(N)+T(N))."
        )

        print()
        print(
            "The next target is to derive a closed formula"
        )

        print(
            "for R-T, because R+T is already forced by C_(2k)."
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

