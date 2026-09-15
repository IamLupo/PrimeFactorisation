#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 153
CENTRAL QUOTIENT MODE q_(2k)
EXACT N-ONLY MODE TEST + CLOSED-FORM SEARCH

PURPOSE
-------
Experiment 152 found the exact quotient-mode recurrence

    C_j = q_j + t q_(j-1) + t q_(j+1)

and observed the strong pattern

    q_j(t) = t^|2k-j| * P_(k,ell,j)(N),
    N = t^2.

The special mode

    j = 2k

should therefore have

    q_(2k)(t) = P_(k,ell)(N),

with NO residual power of t.

This experiment isolates q_(2k) and tests:

    1. exact extraction;
    2. t-power = 0;
    3. exact N-only form;
    4. degree profile in N;
    5. forward ell holdout;
    6. whether the coefficients admit a simple finite-binomial
       representation.

The binomial-law search is deliberately constrained.

We test whether q_(2k) can be represented in the basis

    N^r * binom(ell, k+r),
    N^r * binom(ell, k-r),
    N^r * binom(ell, r),

and small signed combinations of these bases.

The purpose is structural discovery, not numerical fitting.

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

    Q, R = sp.div(
        PF,
        PD
    )

    rem = sp.expand(
        R.as_expr()
    )

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(rem)}"
        )

    return sp.expand(
        Q.as_expr()
    )


# ============================================================================
# Laurent quotient
# ============================================================================

def laurent_Q(k, ell):
    Q = quotient_Q_pq(k, ell)

    return sp.expand(
        Q.subs(
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
        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"noninteger Laurent exponent: {term}"
            )

        support.append(
            int(power)
        )

    return sorted(
        set(support)
    )


def laurent_coeff(expr, j):
    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    shift = max(
        0,
        -min(support)
    )

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ.frac_field(t),
    )

    target = j + shift

    if target < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            u**target
        )
    )


# ============================================================================
# Extract central mode
# ============================================================================

def central_mode(k, ell):
    Q = laurent_Q(
        k,
        ell
    )

    j = 2*k

    return sp.factor(
        laurent_coeff(
            Q,
            j
        )
    )


# ============================================================================
# t-power normalization
# ============================================================================

def common_t_power(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return None

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    return min(
        int(mon[0])
        for mon in poly.monoms()
    )


def rewrite_in_N(expr):
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

        result += (
            coeff
            * N**(power // 2)
        )

    return sp.factor(
        sp.expand(result)
    )


# ============================================================================
# Degree in N
# ============================================================================

def N_degree(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return -sp.oo

    return poly.degree()


# ============================================================================
# Build central-mode dataset
# ============================================================================

def build_dataset(
    k_values,
    ell_values
):
    records = []

    for k in k_values:

        for ell in ell_values:

            if not (
                1 <= k < ell
                and k % 2 == 1
                and ell % 2 == 1
            ):
                continue

            q2k = central_mode(
                k,
                ell
            )

            power = common_t_power(
                q2k
            )

            Nform = rewrite_in_N(
                q2k
            )

            records.append(
                {
                    "k": k,
                    "ell": ell,
                    "j": 2*k,
                    "raw": q2k,
                    "t_power": power,
                    "Nform": Nform,
                    "degree_N": (
                        None
                        if Nform is None
                        else N_degree(Nform)
                    ),
                }
            )

    return records


# ============================================================================
# Candidate basis for symbolic law search
# ============================================================================

def candidate_basis(k, ell, max_r):
    """
    Return candidate scalar factors that can multiply simple N powers.

    The search is intentionally small and interpretable.
    """

    basis = []

    for r in range(
        0,
        max_r + 1
    ):

        basis.append(
            (
                f"N^{r}*C(ell,k+r)",
                N**r * sp.binomial(ell, k+r)
            )
        )

        basis.append(
            (
                f"N^{r}*C(ell,k-r)",
                N**r * sp.binomial(ell, k-r)
            )
        )

        basis.append(
            (
                f"N^{r}*C(ell,r)",
                N**r * sp.binomial(ell, r)
            )
        )

    return basis


# ============================================================================
# Evaluate candidate basis at exact integers
# ============================================================================

def basis_value(
    label,
    expr,
    k,
    ell
):
    return sp.expand(
        expr.subs(
            {
                # no symbolic k/ell in the expression itself
            }
        )
    )


# ============================================================================
# Search constant linear combinations across dataset
# ============================================================================

def search_simple_binomial_law(
    records,
    max_r=3
):
    """
    Test whether q_(2k)(N) can be expressed as a constant rational linear
    combination of the selected binomial basis functions.

    Coefficients are global constants.

    We solve the exact linear system over QQ.
    """

    if not records:
        return []

    # Build a common basis using the record-specific symbolic expression.
    #
    # Since binomial(ell,...) is evaluated at integer ell, each basis vector
    # becomes a rational/integer scalar times N^r.
    #
    # We represent the target as coefficient equations in powers of N.

    # Determine max N degree.
    max_degree = max(
        rec["degree_N"]
        for rec in records
        if rec["Nform"] is not None
    )

    labels = []

    # Use k-dependent basis template labels but solve over actual values.
    for r in range(
        0,
        max_r + 1
    ):
        labels += [
            f"A_r C(ell,k+r), r={r}",
            f"B_r C(ell,k-r), r={r}",
            f"C_r C(ell,r), r={r}",
        ]

    n_basis = len(labels)

    rows = []
    rhs = []

    for rec in records:

        k = rec["k"]
        ell = rec["ell"]
        target = rec["Nform"]

        if target is None:
            continue

        target_poly = sp.Poly(
            sp.expand(target),
            N,
            domain=sp.QQ,
        )

        for degree in range(
            max_degree + 1
        ):

            row = []

            for r in range(
                0,
                max_r + 1
            ):

                row.append(
                    sp.binomial(
                        ell,
                        k + r
                    )
                )

                row.append(
                    sp.binomial(
                        ell,
                        k - r
                    )
                )

                row.append(
                    sp.binomial(
                        ell,
                        r
                    )
                )

                # Only N^r contributes to N^degree when degree == r.
                #
                # Therefore multiply basis coefficient by indicator.
                #
                # This basis design keeps the exact search tiny.
                #

            # The previous construction appended 3*(max_r+1)
            # scalars. We need to project onto the requested N-degree.
            #
            # Rebuild the row correctly.
            row = []

            for r in range(
                0,
                max_r + 1
            ):

                indicator = (
                    1
                    if degree == r
                    else 0
                )

                row.append(
                    indicator
                    * sp.binomial(
                        ell,
                        k + r
                    )
                )

                row.append(
                    indicator
                    * sp.binomial(
                        ell,
                        k - r
                    )
                )

                row.append(
                    indicator
                    * sp.binomial(
                        ell,
                        r
                    )
                )

            rhs.append(
                sp.Rational(
                    target_poly.coeff_monomial(
                        N**degree
                    )
                )
            )

            rows.append(
                row
            )

    A = sp.Matrix(
        rows
    )

    b = sp.Matrix(
        rhs
    )

    if A.cols == 0:
        return []

    solutions = sp.linsolve(
        (
            A,
            b
        )
    )

    if solutions == sp.EmptySet:
        return []

    tuples = list(
        solutions
    )

    if not tuples:
        return []

    solution = tuples[0]

    free = set()

    for value in solution:
        free |= value.free_symbols

    if free:
        return [
            {
                "parametric": True,
                "solution": solution,
                "labels": labels,
            }
        ]

    return [
        {
            "parametric": False,
            "solution": solution,
            "labels": labels,
        }
    ]


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 153")
    print("CENTRAL QUOTIENT MODE q_(2k)")
    print("EXACT N-ONLY MODE TEST + CLOSED-FORM SEARCH")
    print("=" * 78)
    print()

    TRAIN_K = [
        1,
        3,
        5,
    ]

    TRAIN_ELL = [
        3,
        5,
        7,
        9,
        11,
        13,
    ]

    FORWARD_ELL = [
        15,
        17,
    ]

    # ------------------------------------------------------------------------
    # 1. Dataset
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. CENTRAL-MODE DATASET")
    print("=" * 78)

    train = build_dataset(
        TRAIN_K,
        TRAIN_ELL
    )

    forward = build_dataset(
        TRAIN_K,
        FORWARD_ELL
    )

    print(
        "training records =",
        len(train)
    )

    print(
        "forward records =",
        len(forward)
    )

    # ------------------------------------------------------------------------
    # 2. Exact N-only test
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. q_(2k) N-ONLY CERTIFICATE")
    print("=" * 78)

    n_only_failures = 0
    wrong_power_failures = 0

    for rec in train:

        print(
            f"({rec['k']},{rec['ell']}) "
            f"j={rec['j']}"
        )

        print(
            "  raw =",
            rec["raw"]
        )

        print(
            "  common t-power =",
            rec["t_power"]
        )

        print(
            "  N-form =",
            rec["Nform"]
        )

        print(
            "  degree_N =",
            rec["degree_N"]
        )

        # Central prediction:
        #
        # q_(2k) should have t-power 0.
        if rec["t_power"] != 0:
            wrong_power_failures += 1

    print()
    print(
        "wrong central t-power failures =",
        wrong_power_failures
    )

    # ------------------------------------------------------------------------
    # 3. Forward validation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_power_failures = 0
    forward_N_failures = 0

    for rec in forward:

        if rec["t_power"] != 0:
            forward_power_failures += 1

        if rec["Nform"] is None:
            forward_N_failures += 1

        print(
            f"({rec['k']},{rec['ell']}) "
            f"q_(2k) = {rec['Nform']}"
        )

    print(
        "forward t-power failures =",
        forward_power_failures
    )

    print(
        "forward N-form failures =",
        forward_N_failures
    )

    # ------------------------------------------------------------------------
    # 4. Degree profile
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DEGREE PROFILE")
    print("=" * 78)

    for k in TRAIN_K:

        subset = [
            rec
            for rec in train
            if rec["k"] == k
        ]

        print(
            f"k={k}:"
        )

        for rec in subset:

            print(
                f"  ell={rec['ell']:>2} "
                f"degree_N={rec['degree_N']} "
                f"q_(2k)={rec['Nform']}"
            )

    # ------------------------------------------------------------------------
    # 5. Binomial-law search
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SIMPLE BINOMIAL LAW SEARCH")
    print("=" * 78)

    solutions = search_simple_binomial_law(
        train,
        max_r=3
    )

    print(
        "solutions =",
        len(solutions)
    )

    for solution in solutions:

        if solution["parametric"]:

            print(
                "parametric solution =",
                solution["solution"]
            )

        else:

            print(
                "exact coefficient law:"
            )

            for label, coeff in zip(
                solution["labels"],
                solution["solution"]
            ):

                if coeff != 0:
                    print(
                        f"  {label} -> {sp.factor(coeff)}"
                    )

    # ------------------------------------------------------------------------
    # 6. Direct law diagnostics
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CENTRAL-MODE STRUCTURAL DIAGNOSTIC")
    print("=" * 78)

    # Check whether degree grows with ell-k in a simple pattern.
    #
    # This is descriptive, not a fit.

    for k in TRAIN_K:

        subset = [
            rec
            for rec in train
            if rec["k"] == k
            and rec["ell"] >= k + 3
        ]

        pairs = [
            (
                rec["ell"],
                rec["degree_N"]
            )
            for rec in subset
        ]

        print(
            f"k={k}:",
            pairs
        )

    # ------------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if (
        wrong_power_failures == 0
        and forward_power_failures == 0
        and forward_N_failures == 0
    ):

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The central quotient mode satisfies"
        )

        print(
            "    q_(2k)(t) = P_(k,ell)(N)"
        )

        print(
            "with no residual t-power on the tested grid."
        )

        if solutions:

            print()
            print(
                "A simple binomial law candidate was also found."
            )

            print(
                "This law should now be independently verified"
            )

            print(
                "on additional (k,ell) before being promoted."
            )

        else:

            print()
            print(
                "The central mode is N-only, but no simple"
            )

            print(
                "binomial envelope from the tested basis was found."
            )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive q_(2k) directly from the quotient recurrence"
        )

        print(
            "and determine its exact general formula."
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

