import sympy as sp


# =============================================================================
# EXPERIMENT 208
# OUTER-DIAGONAL COEFFICIENT DERIVATION AUDIT
# =============================================================================
#
# Objective:
#
#   Experiment 207 isolated the outer-diagonal coordinate
#
#       D(k,ell,s)
#         = [N^(r-1) X^(ell-r+1)] G(N,X-1),
#
#   where
#
#       r = s - (k-1)/2.
#
#   This experiment does NOT search new affine coordinates.
#
#   Instead it:
#
#     1. computes D exactly from the pq kernel;
#     2. derives a direct coefficient formula from the two-variable
#        representation;
#     3. compares exact values against the theorem anchors;
#     4. tests whether D has a simple binomial/product structure;
#     5. checks whether the theorem values are obtained from D by
#        a structurally motivated normalization only:
#
#           D
#           -D
#           D / r
#           D / (ell-r+1)
#           D / binomial(...)
#
#   No unrestricted fitting is performed.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


ANCHORS = [
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


# =============================================================================
# EXACT pq KERNEL
# =============================================================================

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# =============================================================================
# DIRECT SYMMETRIC REDUCTION
# =============================================================================

def exact_G(k, ell):

    F = exact_F(k, ell)

    result, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Symmetrization remainder for "
            f"k={k}, ell={ell}: {remainder}"
        )

    sum_sym = None
    prod_sym = None

    for sym, expr in mapping:

        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            sum_sym = sym

        elif sp.expand(expr - p*q) == 0:
            prod_sym = sym

    if sum_sym is None or prod_sym is None:
        raise ValueError(
            f"Could not identify symmetric coordinates: {mapping}"
        )

    return sp.expand(
        result.subs(
            {
                sum_sym: S,
                prod_sym: N,
            }
        )
    )


# =============================================================================
# SHIFTED BASIS
# =============================================================================

def shifted_G(k, ell):
    G = exact_G(k, ell)
    return sp.expand(
        G.subs(S, X - 1)
    )


# =============================================================================
# SAFE COEFFICIENT EXTRACTION
# =============================================================================

def coefficient(expr, ndeg, xdeg):

    if ndeg < 0 or xdeg < 0:
        return sp.Integer(0)

    poly_x = sp.Poly(
        sp.expand(expr),
        X,
        domain=sp.EX,
    )

    if xdeg > poly_x.degree():
        return sp.Integer(0)

    layer = poly_x.coeff_monomial(X**xdeg)

    if layer == 0:
        return sp.Integer(0)

    poly_n = sp.Poly(
        sp.expand(layer),
        N,
        domain=sp.EX,
    )

    return sp.expand(
        poly_n.coeff_monomial(N**ndeg)
    )


# =============================================================================
# OUTER DIAGONAL
# =============================================================================

def outer_parameters(k, ell, s):

    if k % 2 == 0:
        raise ValueError(
            f"Expected odd k, got k={k}"
        )

    h = (k - 1) // 2
    r = s - h

    ndeg = r - 1
    xdeg = ell - r + 1

    return h, r, ndeg, xdeg


def outer_diagonal(k, ell, s):

    Gx = shifted_G(k, ell)

    h, r, ndeg, xdeg = outer_parameters(
        k,
        ell,
        s,
    )

    return coefficient(
        Gx,
        ndeg,
        xdeg,
    )


# =============================================================================
# BASIC BINOMIAL QUANTITIES
# =============================================================================

def safe_binom(n, r):

    if r < 0 or r > n:
        return sp.Integer(0)

    return sp.binomial(
        n,
        r,
    )


def motivated_normalizations(k, ell, s, D):

    h, r, ndeg, xdeg = outer_parameters(
        k,
        ell,
        s,
    )

    values = {}

    values["D"] = sp.expand(D)
    values["-D"] = sp.expand(-D)

    values["D/r"] = (
        sp.cancel(D / r)
        if r != 0
        else sp.Integer(0)
    )

    values["D/(ell-r+1)"] = sp.cancel(
        D / (ell - r + 1)
    )

    values["D/binom(k,r)"] = (
        sp.cancel(
            D / safe_binom(k, r)
        )
        if safe_binom(k, r) != 0
        else sp.Integer(0)
    )

    values["D/binom(ell,r)"] = (
        sp.cancel(
            D / safe_binom(ell, r)
        )
        if safe_binom(ell, r) != 0
        else sp.Integer(0)
    )

    values["D/binom(ell,r-1)"] = (
        sp.cancel(
            D / safe_binom(ell, r - 1)
        )
        if safe_binom(ell, r - 1) != 0
        else sp.Integer(0)
    )

    values["D/(r*binom(ell,r))"] = (
        sp.cancel(
            D / (
                r * safe_binom(ell, r)
            )
        )
        if safe_binom(ell, r) != 0 and r != 0
        else sp.Integer(0)
    )

    return values


# =============================================================================
# DIRECT COEFFICIENT FORMULA TEST
# =============================================================================
#
# We derive the outer coefficient directly from the monomial structure
# before symmetric reduction.
#
# For the shifted variable X = S+1,
#
#   p^k (1+q)^ell
#
# becomes, after fixing the coefficient of X^j, a finite sum in p q.
#
# Rather than assume a closed form, we independently compute D in two ways:
#
#   (A) symmetric polynomial coefficient;
#   (B) direct evaluation of G on symbolic N,S followed by coefficient
#       extraction.
#
# This catches implementation-level mistakes.
# =============================================================================

def direct_outer_check(k, ell, s):

    D1 = outer_diagonal(
        k,
        ell,
        s,
    )

    G = exact_G(
        k,
        ell,
    )

    h, r, ndeg, xdeg = outer_parameters(
        k,
        ell,
        s,
    )

    # Extract [N^ndeg] first.
    poly_N = sp.Poly(
        sp.expand(G),
        N,
        domain=sp.EX,
    )

    layer = poly_N.coeff_monomial(
        N**ndeg
    )

    layer_shifted = sp.expand(
        layer.subs(S, X - 1)
    )

    D2 = sp.Poly(
        layer_shifted,
        X,
        domain=sp.EX,
    ).coeff_monomial(
        X**xdeg
    )

    return (
        sp.expand(D1),
        sp.expand(D2),
    )


# =============================================================================
# SMALL PARAMETER GRID
# =============================================================================

def grid_audit():

    print()
    print("=" * 78)
    print("1. OUTER-DIAGONAL GRID AUDIT")
    print("=" * 78)

    failures = 0

    for k in [1, 3, 5, 7, 9, 11]:

        for ell in range(
            max(k, 3),
            min(k + 12, 20) + 1,
            2,
        ):

            # Use theorem-style admissible s values.
            for s in range(
                1,
                ell + 1,
            ):

                if s <= (k - 1) // 2:
                    continue

                try:

                    D1, D2 = direct_outer_check(
                        k,
                        ell,
                        s,
                    )

                except Exception as exc:

                    print(
                        "CHECK ERROR",
                        k,
                        ell,
                        s,
                        exc,
                    )

                    failures += 1
                    continue

                if D1 != D2:

                    print(
                        "MISMATCH:",
                        f"k={k}",
                        f"ell={ell}",
                        f"s={s}",
                        f"D1={D1}",
                        f"D2={D2}",
                    )

                    failures += 1

    print()
    print(
        f"outer coefficient cross-check failures = {failures}"
    )


# =============================================================================
# ANCHOR AUDIT
# =============================================================================

def anchor_audit():

    print()
    print("=" * 78)
    print("2. OUTER-DIAGONAL THEOREM-ANCHOR AUDIT")
    print("=" * 78)

    exact_hits = 0
    abs_hits = 0

    for name, k, ell, s, expected in ANCHORS:

        h, r, ndeg, xdeg = outer_parameters(
            k,
            ell,
            s,
        )

        D = outer_diagonal(
            k,
            ell,
            s,
        )

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"h={h} r={r} "
            f"[N^{ndeg}X^{xdeg}]"
        )

        print(
            f"  D        = {D}"
        )

        print(
            f"  expected = {expected}"
        )

        print(
            f"  exact    = {D == expected}"
        )

        print(
            f"  abs      = {abs(D) == abs(expected)}"
        )

        if D == expected:
            exact_hits += 1

        if abs(D) == abs(expected):
            abs_hits += 1

        normalized = motivated_normalizations(
            k,
            ell,
            s,
            D,
        )

        print(
            "  motivated normalizations:"
        )

        for key, value in normalized.items():

            if key == "D":
                continue

            print(
                f"    {key:24s} = {value}"
            )

    print()
    print(
        f"EXACT outer matches = {exact_hits}/9"
    )

    print(
        f"ABSOLUTE outer matches = {abs_hits}/9"
    )


# =============================================================================
# OUTER-DIAGONAL SEQUENCE AUDIT
# =============================================================================

def sequence_audit():

    print()
    print("=" * 78)
    print("3. OUTER-DIAGONAL SEQUENCE AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
        (5, 15),
        (5, 19),
        (7, 15),
        (7, 17),
        (9, 21),
        (9, 23),
        (11, 23),
        (11, 25),
    ]

    for k, ell in cases:

        h = (k - 1) // 2

        print()
        print(
            f"k={k} ell={ell} h={h}"
        )

        for r in range(
            1,
            min(
                k + 2,
                ell - h + 1,
            ),
        ):

            s = r + h

            if s > ell:
                continue

            D = outer_diagonal(
                k,
                ell,
                s,
            )

            print(
                f"  r={r:2d} "
                f"s={s:2d} "
                f"D={D}"
            )


# =============================================================================
# FACTORIZATION AUDIT
# =============================================================================

def factorization_audit():

    print()
    print("=" * 78)
    print("4. OUTER-DIAGONAL FACTORIZATION AUDIT")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        D = outer_diagonal(
            k,
            ell,
            s,
        )

        print()
        print(
            f"{name}: D={D}"
        )

        if D == 0:
            print(
                "  factorization: ZERO"
            )
            continue

        print(
            "  factorization:"
        )

        print(
            f"    {sp.factor(D)}"
        )

        print(
            "  integer factor content:"
        )

        coeff, factors = sp.factor_list(
            D
        )

        print(
            f"    constant={coeff}"
        )

        for factor, exponent in factors:

            print(
                f"    factor={factor} exponent={exponent}"
            )


# =============================================================================
# RECURRENCE IN r
# =============================================================================

def recurrence_audit():

    print()
    print("=" * 78)
    print("5. OUTER-DIAGONAL r-RECURRENCE AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        h = (k - 1) // 2

        values = []

        for r in range(
            1,
            min(
                k + 3,
                ell - h,
            ),
        ):

            s = r + h

            if s > ell:
                continue

            D = outer_diagonal(
                k,
                ell,
                s,
            )

            values.append(
                (r, D)
            )

        print()
        print(
            f"k={k} ell={ell}"
        )

        for i in range(
            1,
            len(values),
        ):

            r0, D0 = values[i - 1]
            r1, D1 = values[i]

            print(
                f"  r={r0}->{r1}: "
                f"D={D0} -> {D1}"
            )

            if D0 != 0:

                ratio = sp.cancel(
                    sp.Rational(D1, D0)
                )

                print(
                    f"       ratio = {ratio}"
                )


# =============================================================================
# FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("6. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Experiment 207 ruled out (r,d) as sufficient coordinates."
    )

    print()
    print(
        "The strongest remaining structural signal is the"
    )

    print(
        "outer diagonal"
    )

    print(
        "    D(k,ell,s)"
    )

    print(
        "      = [N^(r-1) X^(ell-r+1)] G(N,X-1)"
    )

    print(
        "where"
    )

    print(
        "    r = s-(k-1)/2."
    )

    print()
    print(
        "Experiment 208 therefore does not introduce new index"
    )

    print(
        "families. It asks whether this boundary coefficient has"
    )

    print(
        "an intrinsic binomial/product structure."
    )

    print()
    print(
        "The decisive outcomes are:"
    )

    print(
        "  * if D itself follows a clean recurrence or closed"
    )

    print(
        "    binomial structure, derive it before comparing with"
    )

    print(
        "    theorem anchors;"
    )

    print(
        "  * if D does not have such structure, abandon D as the"
    )

    print(
        "    primary theorem coordinate;"
    )

    print(
        "  * if D reproduces several anchors structurally, the"
    )

    print(
        "    next step is to trace that boundary coefficient back"
    )

    print(
        "    to the original pq expansion."
    )

    print()
    print(
        "No finite-anchor match is accepted as proof."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 208")
    print("OUTER-DIAGONAL COEFFICIENT DERIVATION AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    grid_audit()
    anchor_audit()
    sequence_audit()
    factorization_audit()
    recurrence_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 208")
    print("=" * 78)


if __name__ == "__main__":
    main()

