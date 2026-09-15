#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 150
DIRECT LAURENT-KERNEL / N-ONLY INFORMATION TEST

PURPOSE
-------
Experiment 149 ruled out low-order recurrences in ell whose coefficients
depend only on N.

We now return to the original finite kernel BEFORE introducing S.

Use the exact Laurent substitution

    p = t*u
    q = t*u^(-1)

so that

    N = p*q = t^2
    S = p+q = t*(u+u^(-1)).

For fixed (k,ell),

    F_(k,ell)(t*u,t*u^(-1))

is a finite Laurent polynomial

    sum_j C_j(t) u^j.

This experiment classifies the Laurent coefficients C_j(t).

The central question is:

    Does any nontrivial informative Laurent coefficient depend ONLY
    on t^2=N, with no hidden u-branch information?

We classify every coefficient into:

    A. N-only:
         C_j(t) is a polynomial/rational expression in t^2 only.

    B. t-dependent but branch-free:
         coefficient contains t with odd powers but no u ambiguity.
         Such a coefficient is not directly an N-only invariant.

    C. branch-dependent:
         information is carried by a nonzero Laurent exponent j.

We additionally inspect the coefficient pairing

    C_j(t) u^j + C_(-j)(t) u^(-j),

and verify the symmetry

    C_j = C_(-j)

expected from p <-> q.

The important outcome is NOT whether arbitrary coefficients exist.
The target is an informative quantity that survives as a function of

    N = t^2

without requiring a choice of u.

This is a structural obstruction/certificate experiment.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO POLYNOMIAL FITTING
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
# Direct Laurent substitution
# ============================================================================

def laurent_kernel(k, ell):
    F = kernel_F(k, ell)

    return sp.expand(
        F.subs(
            {
                p: t*u,
                q: t/u,
            }
        )
    )


# ============================================================================
# Laurent support and coefficient extraction
# ============================================================================

def laurent_support(expr):
    expr = sp.expand(expr)

    support = []

    for term in sp.Add.make_args(expr):

        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"non-integer Laurent exponent in term: {term}"
            )

        support.append(
            int(power)
        )

    return sorted(
        set(support)
    )


def laurent_coeff(expr, j):
    """
    Exact coefficient [u^j].
    """

    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    min_power = min(support)
    shift = max(
        0,
        -min_power
    )

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ.frac_field(t),
    )

    wanted = j + shift

    if wanted < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            u**wanted
        )
    )


# ============================================================================
# Classify a coefficient as function of N=t^2
# ============================================================================

def classify_t_polynomial(expr):
    """
    Determine whether expr depends only on t^2.

    Since the Laurent coefficient itself contains no u, this is purely
    a parity question in t.

    Returns:
        N_ONLY
        T_ODD
        ZERO
    """

    expr = sp.expand(expr)

    if expr == 0:
        return "ZERO"

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    for (power,), coeff in poly.terms():

        if power % 2 == 1 and coeff != 0:
            return "T_ODD"

    return "N_ONLY"


# ============================================================================
# Rewrite even t-polynomial as polynomial in N=t^2
# ============================================================================

def rewrite_as_N(expr):
    """
    Replace t^(2r) by N^r exactly.
    """

    N = sp.Symbol("N")

    poly = sp.Poly(
        sp.expand(expr),
        t,
        domain=sp.QQ,
    )

    result = sp.Integer(0)

    for (power,), coeff in poly.terms():

        power = int(power)

        if power % 2 != 0:
            raise ArithmeticError(
                f"odd t-power encountered in N rewrite: {expr}"
            )

        result += coeff * N**(power // 2)

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Branch dependence classification
# ============================================================================

def classify_laurent_term(j, coeff):
    if coeff == 0:
        return "ZERO"

    if j != 0:
        return "BRANCH_DEPENDENT"

    return classify_t_polynomial(
        coeff
    )


# ============================================================================
# Pairing check
# ============================================================================

def check_pair_symmetry(expr):
    support = laurent_support(expr)

    failures = []

    for j in support:

        cj = laurent_coeff(
            expr,
            j
        )

        cminus = laurent_coeff(
            expr,
            -j
        )

        if sp.expand(cj - cminus) != 0:

            failures.append(
                (
                    j,
                    sp.factor(cj),
                    sp.factor(cminus)
                )
            )

    return failures


# ============================================================================
# Test whether the entire Laurent polynomial has only N-only coefficients
# ============================================================================
#
# This is intentionally NOT expected to pass: the whole point of the
# experiment is to identify where hidden branch information enters.
# ============================================================================

def informative_profile(expr):
    support = laurent_support(expr)

    profile = []

    for j in support:

        coeff = laurent_coeff(
            expr,
            j
        )

        cls = classify_laurent_term(
            j,
            coeff
        )

        profile.append(
            {
                "j": j,
                "coeff": sp.factor(coeff),
                "class": cls
            }
        )

    return profile


# ============================================================================
# Sample pairs
# ============================================================================

PAIRS = [
    (1, 3),
    (1, 5),
    (1, 7),
    (1, 9),
    (3, 5),
    (3, 7),
    (3, 9),
    (5, 7),
    (5, 9),
]


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 150")
    print("DIRECT LAURENT-KERNEL / N-ONLY INFORMATION TEST")
    print("=" * 78)
    print()

    total_pair_failures = 0
    total_symmetry_failures = 0
    total_N_only_nonzero = 0
    total_branch_terms = 0
    total_t_odd_terms = 0

    # ------------------------------------------------------------------------
    # 1. Direct Laurent construction
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. DIRECT LAURENT KERNELS")
    print("=" * 78)

    for k, ell in PAIRS:

        print()
        print(
            f"({k},{ell})"
        )

        try:
            R = laurent_kernel(
                k,
                ell
            )
        except Exception as exc:
            print(
                "  construction failure:",
                type(exc).__name__,
                str(exc)
            )
            total_pair_failures += 1
            continue

        print(
            "  support =",
            laurent_support(R)
        )

        print(
            "  R(t,u) ="
        )

        print(
            "   ",
            sp.factor(R)
        )

        # --------------------------------------------------------------------
        # Symmetry
        # --------------------------------------------------------------------

        symmetry_failures = check_pair_symmetry(
            R
        )

        if symmetry_failures:

            total_symmetry_failures += 1

            print(
                "  u <-> u^(-1) symmetry = FAIL"
            )

            for failure in symmetry_failures[:5]:
                print(
                    "    j=",
                    failure[0],
                    "c_j=",
                    failure[1],
                    "c_-j=",
                    failure[2]
                )

        else:

            print(
                "  u <-> u^(-1) symmetry = PASS"
            )

        # --------------------------------------------------------------------
        # Profile
        # --------------------------------------------------------------------

        profile = informative_profile(
            R
        )

        n_only = 0
        branch = 0
        t_odd = 0

        print()
        print(
            "  COEFFICIENT CLASSIFICATION"
        )

        for item in profile:

            j = item["j"]
            coeff = item["coeff"]
            cls = item["class"]

            if cls == "N_ONLY":
                n_only += 1
                total_N_only_nonzero += 1

                print(
                    f"    j={j:>3} "
                    f"[N_ONLY] "
                    f"{coeff} "
                    f"-> "
                    f"{rewrite_as_N(coeff)}"
                )

            elif cls == "T_ODD":
                t_odd += 1
                total_t_odd_terms += 1

                print(
                    f"    j={j:>3} "
                    f"[T_ODD] "
                    f"{coeff}"
                )

            elif cls == "BRANCH_DEPENDENT":
                branch += 1
                total_branch_terms += 1

                print(
                    f"    j={j:>3} "
                    f"[BRANCH] "
                    f"{coeff}"
                )

        print()
        print(
            f"  nonzero N-only coefficients = {n_only}"
        )

        print(
            f"  nonzero t-odd coefficients = {t_odd}"
        )

        print(
            f"  branch-dependent coefficients = {branch}"
        )

        # --------------------------------------------------------------------
        # Check whether any nonzero j=0 coefficient exists.
        # --------------------------------------------------------------------

        c0 = laurent_coeff(
            R,
            0
        )

        if c0 != 0:

            print()
            print(
                "  ZERO-MODE / N-ONLY CANDIDATE"
            )

            print(
                "    C_0(t) =",
                sp.factor(c0)
            )

            if classify_t_polynomial(c0) == "N_ONLY":

                print(
                    "    C_0(N) =",
                    rewrite_as_N(c0)
                )

                print(
                    "    STATUS = N-ONLY CANDIDATE"
                )

            else:

                print(
                    "    STATUS = NOT N-ONLY"
                )

        # --------------------------------------------------------------------
        # Pair invariant
        # --------------------------------------------------------------------
        #
        # For j>0:
        #
        #     c_j u^j + c_j u^-j
        #
        # depends on u+u^-1 unless c_j=0.
        #
        # We report this explicitly.
        # --------------------------------------------------------------------

        positive_support = [
            j
            for j in laurent_support(R)
            if j > 0
        ]

        branch_informative = [
            j for j in positive_support
            if laurent_coeff(R, j) != 0
        ]

        print()
        print(
            "  positive branch modes =",
            branch_informative
        )

    # ------------------------------------------------------------------------
    # 2. Global summary
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. GLOBAL INFORMATION SUMMARY")
    print("=" * 78)

    print(
        "pairs tested =",
        len(PAIRS)
    )

    print(
        "pair construction failures =",
        total_pair_failures
    )

    print(
        "Laurent symmetry failures =",
        total_symmetry_failures
    )

    print(
        "nonzero N-only coefficients =",
        total_N_only_nonzero
    )

    print(
        "nonzero t-odd coefficients =",
        total_t_odd_terms
    )

    print(
        "branch-dependent coefficients =",
        total_branch_terms
    )

    # ------------------------------------------------------------------------
    # 3. Dedicated zero-mode analysis
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ZERO-MODE ANALYSIS")
    print("=" * 78)

    zero_modes = []

    for k, ell in PAIRS:

        R = laurent_kernel(
            k,
            ell
        )

        c0 = laurent_coeff(
            R,
            0
        )

        zero_modes.append(
            (
                k,
                ell,
                sp.factor(c0)
            )
        )

        print(
            f"({k},{ell}) C_0(t) =",
            sp.factor(c0)
        )

        if c0 != 0 and classify_t_polynomial(c0) == "N_ONLY":

            print(
                "   N-only form =",
                rewrite_as_N(c0)
            )

    # ------------------------------------------------------------------------
    # 4. Structural obstruction test
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL OBSTRUCTION TEST")
    print("=" * 78)

    # A Laurent mode j != 0 necessarily carries u^j+u^-j.
    # Therefore it is branch-dependent unless its coefficient vanishes.
    #
    # We determine whether any tested kernel has:
    #
    #   - no positive Laurent modes; OR
    #   - only N-only zero mode.
    #
    pure_N_pairs = 0

    for k, ell, c0 in zero_modes:

        R = laurent_kernel(
            k,
            ell
        )

        positive = [
            j
            for j in laurent_support(R)
            if j > 0
            and laurent_coeff(R, j) != 0
        ]

        if (
            not positive
            and c0 != 0
            and classify_t_polynomial(c0) == "N_ONLY"
        ):

            pure_N_pairs += 1

    print(
        "pure N-only kernels =",
        pure_N_pairs
    )

    # ------------------------------------------------------------------------
    # 5. Final diagnostic
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if total_symmetry_failures != 0:

        print(
            "STATUS = FAIL"
        )

        print(
            "The expected p<->q / u<->u^-1 symmetry failed."
        )

        return

    if pure_N_pairs > 0:

        print(
            "STATUS = N-ONLY SIGNAL FOUND"
        )

        print()
        print(
            "At least one tested kernel has no nonzero branch"
        )

        print(
            "Laurent modes and contains a nonzero N-only zero mode."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "determine whether that zero mode can be generalized"
        )

        print(
            "to a family sufficient for trace recovery."
        )

    elif total_N_only_nonzero > 0:

        print(
            "STATUS = PARTIAL N-ONLY SIGNAL"
        )

        print()
        print(
            "Nonzero N-only Laurent coefficients exist, but"
        )

        print(
            "branch-dependent modes remain."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "test whether the N-only modes form a closed recurrence."
        )

    else:

        print(
            "STATUS = BRANCH DEPENDENCE CONFIRMED"
        )

        print()
        print(
            "No nonzero N-only Laurent coefficient was found"
        )

        print(
            "among the tested kernels."
        )

        print(
            "All informative modes carry nonzero Laurent exponent"
        )

        print(
            "and therefore retain hidden u / p:q branch information."
        )

        print()
        print(
            "This strongly suggests that the KAPPA construction"
        )

        print(
            "does not itself produce an N-only observable."
        )

        print()
        print(
            "NEXT MATHEMATICAL QUESTION:"
        )

        print(
            "look for an external N-computable operation that"
        )

        print(
            "selects or evaluates the hidden Laurent mode."
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

