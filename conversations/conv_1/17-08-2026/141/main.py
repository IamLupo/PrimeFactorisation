import sympy as sp


# =============================================================================
# EXPERIMENT 207
# EFFECTIVE-INDEX / BOUNDARY-DISTANCE AUDIT
# =============================================================================
#
# Goal:
#   Determine whether the theorem index s is shifted by k in the exact
#   symmetric kernel.
#
# No unrestricted affine search.
# No previous output is read.
# No guessed pq -> (k,ell,s) map is used.
#
# We test only index transformations that are directly motivated by:
#
#     h = (k-1)/2
#     r = s-h
#     d = ell-s
#     u = s-h-1
#     v = ell-s-h
#
# and the corresponding boundary distances in the shifted S+1 basis.
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
# SHIFTED S+1 BASIS
# =============================================================================

def shifted_coefficients(k, ell):

    G = exact_G(k, ell)

    shifted = sp.expand(
        G.subs(S, X - 1)
    )

    if shifted == 0:
        return {0: sp.Integer(0)}

    poly = sp.Poly(
        shifted,
        X,
        domain=sp.EX,
    )

    degree = poly.degree()

    if degree is None:
        return {0: sp.Integer(0)}

    degree = int(degree)

    return {
        j: sp.expand(
            poly.coeff_monomial(X**j)
        )
        for j in range(degree + 1)
    }


# =============================================================================
# SAFE N COEFFICIENT
# =============================================================================

def N_coeff(poly_expr, d):

    if d < 0:
        return sp.Integer(0)

    expr = sp.expand(poly_expr)

    if expr == 0:
        return sp.Integer(0)

    poly = sp.Poly(
        expr,
        N,
        domain=sp.EX,
    )

    return sp.expand(
        poly.coeff_monomial(N**d)
    )


def get_coeff(coeffs, j, d):

    if j < 0 or d < 0:
        return sp.Integer(0)

    if j not in coeffs:
        return sp.Integer(0)

    return N_coeff(
        coeffs[j],
        d,
    )


# =============================================================================
# EFFECTIVE INDEX SYSTEM
# =============================================================================

def effective_indices(k, ell, s):

    # k is odd in all theorem anchors.
    h = (k - 1) // 2

    return {
        "s": s,
        "h": h,

        # Direct theorem index
        "r0": s,

        # Shift by half-k
        "r1": s - h,
        "r2": s - h - 1,
        "r3": s + h,
        "r4": s + h + 1,

        # Distance from upper boundary
        "d0": ell - s,
        "d1": ell - s - h,
        "d2": ell - s - h - 1,
        "d3": ell - s + h,
        "d4": ell - s + h + 1,
    }


# =============================================================================
# MOTIVATED COORDINATE FAMILIES
# =============================================================================

def candidate_coordinates(k, ell, s):

    h = (k - 1) // 2

    r = s - h
    d = ell - s

    candidates = []

    # -------------------------------------------------------------------------
    # Family A:
    # effective theorem index r
    # -------------------------------------------------------------------------

    candidates.extend([
        ("A0", d, r),
        ("A1", d + 1, r),
        ("A2", d - 1, r),
        ("A3", d, r - 1),
        ("A4", d, r + 1),
    ])

    # -------------------------------------------------------------------------
    # Family B:
    # upper-boundary distance
    # -------------------------------------------------------------------------

    u = d - h

    candidates.extend([
        ("B0", u, r),
        ("B1", u + 1, r),
        ("B2", u - 1, r),
        ("B3", u, r - 1),
        ("B4", u, r + 1),
    ])

    # -------------------------------------------------------------------------
    # Family C:
    # complementary effective index
    # -------------------------------------------------------------------------

    v = ell - s - h

    candidates.extend([
        ("C0", v, r),
        ("C1", v + 1, r),
        ("C2", v - 1, r),
        ("C3", v, r - 1),
        ("C4", v, r + 1),
    ])

    # -------------------------------------------------------------------------
    # Family D:
    # direct shifted theorem index in N-degree and S+1-degree
    # -------------------------------------------------------------------------

    j = ell - r + 1
    ndeg = r - 1

    candidates.extend([
        ("D0", j, ndeg),
        ("D1", j - 1, ndeg),
        ("D2", j + 1, ndeg),
        ("D3", j, ndeg - 1),
        ("D4", j, ndeg + 1),
    ])

    return candidates


# =============================================================================
# ANCHOR COORDINATE AUDIT
# =============================================================================

def anchor_audit():

    print()
    print("=" * 78)
    print("1. EFFECTIVE INDEX COORDINATE AUDIT")
    print("=" * 78)

    exact_counts = {}
    abs_counts = {}

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(k, ell)

        h = (k - 1) // 2
        r = s - h
        d = ell - s

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"h={h} r=s-h={r} d=ell-s={d} "
            f"expected={expected}"
        )

        found = False

        for label, j, ndeg in candidate_coordinates(
            k,
            ell,
            s,
        ):

            value = get_coeff(
                coeffs,
                j,
                ndeg,
            )

            print(
                f"  {label:4s} "
                f"[N^{ndeg} X^{j}] = {value}"
            )

            exact_counts[label] = (
                exact_counts.get(label, 0)
                + int(value == expected)
            )

            abs_counts[label] = (
                abs_counts.get(label, 0)
                + int(abs(value) == abs(expected))
            )

            if value == expected:
                found = True

        if not found:
            print(
                "  no effective-index coordinate match"
            )

    print()
    print("EXACT MATCH COUNTS")

    for label in sorted(exact_counts):
        print(
            f"  {label:4s}: "
            f"{exact_counts[label]}/9"
        )

    print()
    print("ABSOLUTE-VALUE MATCH COUNTS")

    for label in sorted(abs_counts):
        print(
            f"  {label:4s}: "
            f"{abs_counts[label]}/9"
        )


# =============================================================================
# EFFECTIVE INDEX TABLE
# =============================================================================

def index_table():

    print()
    print("=" * 78)
    print("2. EFFECTIVE INDEX TABLE")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        idx = effective_indices(
            k,
            ell,
            s,
        )

        print()
        print(
            f"{name}: "
            f"k={k:2d} ell={ell:2d} s={s:2d} "
            f"expected={expected}"
        )

        for key, value in idx.items():
            print(
                f"  {key:3s} = {value}"
            )


# =============================================================================
# SAME-EFFECTIVE-INDEX COLLISION AUDIT
# =============================================================================

def collision_audit():

    print()
    print("=" * 78)
    print("3. EFFECTIVE-INDEX COLLISION AUDIT")
    print("=" * 78)

    groups = {}

    for name, k, ell, s, expected in ANCHORS:

        h = (k - 1) // 2

        r = s - h
        d = ell - s

        key = (r, d)

        groups.setdefault(
            key,
            []
        ).append(
            (
                name,
                k,
                ell,
                s,
                expected,
            )
        )

    for key, values in groups.items():

        if len(values) < 2:
            continue

        print()
        print(
            f"(r,d)={key}"
        )

        for item in values:
            print(
                f"  {item}"
            )

        vals = {
            item[-1]
            for item in values
        }

        if len(vals) == 1:
            print(
                "  CONSISTENT: same effective coordinates "
                "give same theorem value"
            )
        else:
            print(
                "  NOT CONSISTENT: theorem value differs"
            )


# =============================================================================
# BOUNDARY-DISTANCE COLLISION AUDIT
# =============================================================================

def boundary_collision_audit():

    print()
    print("=" * 78)
    print("4. BOUNDARY-DISTANCE COLLISION AUDIT")
    print("=" * 78)

    groups = {}

    for name, k, ell, s, expected in ANCHORS:

        h = (k - 1) // 2

        r = s - h
        u = ell - s - h

        key = (r, u)

        groups.setdefault(
            key,
            []
        ).append(
            (
                name,
                k,
                ell,
                s,
                expected,
            )
        )

    for key, values in groups.items():

        if len(values) < 2:
            continue

        print()
        print(
            f"(r,u)={key}"
        )

        for item in values:
            print(
                f"  {item}"
            )

        vals = {
            item[-1]
            for item in values
        }

        if len(vals) == 1:
            print(
                "  CONSISTENT"
            )
        else:
            print(
                "  NOT CONSISTENT"
            )


# =============================================================================
# PARITY COMPARISON
# =============================================================================

def parity_audit():

    print()
    print("=" * 78)
    print("5. EFFECTIVE INDEX PARITY AUDIT")
    print("=" * 78)

    cases = [
        (1, 9, 4),
        (1, 10, 4),
        (3, 9, 4),
        (3, 10, 4),
        (5, 11, 5),
        (5, 12, 5),
        (7, 15, 7),
        (7, 16, 7),
        (9, 21, 10),
        (9, 22, 10),
    ]

    for k, ell, s in cases:

        h = (k - 1) // 2
        r = s - h
        d = ell - s

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        # Most natural effective coordinate.
        j = d
        ndeg = r

        value = get_coeff(
            coeffs,
            j,
            ndeg,
        )

        parity = (
            "same"
            if (k - ell) % 2 == 0
            else "different"
        )

        print(
            f"k={k:2d} ell={ell:2d} s={s:2d} "
            f"h={h:2d} r={r:3d} d={d:3d} "
            f"parity={parity:9s} "
            f"[N^{ndeg}X^{j}]={value}"
        )


# =============================================================================
# SMALL LOCAL TRANSFORMATION AUDIT
# =============================================================================

def local_transform_audit():

    print()
    print("=" * 78)
    print("6. LOCAL EFFECTIVE-INDEX TRANSFORMATION AUDIT")
    print("=" * 78)

    transformations = [
        ("r",     lambda r, d: (r, d)),
        ("r-1",   lambda r, d: (r - 1, d)),
        ("r+1",   lambda r, d: (r + 1, d)),
        ("d",     lambda r, d: (r, d)),
        ("d-1",   lambda r, d: (r, d - 1)),
        ("d+1",   lambda r, d: (r, d + 1)),
        ("r,d-1", lambda r, d: (r, d - 1)),
        ("r-1,d", lambda r, d: (r - 1, d)),
        ("r-1,d-1", lambda r, d: (r - 1, d - 1)),
        ("r+1,d+1", lambda r, d: (r + 1, d + 1)),
    ]

    counts = {
        name: 0
        for name, _ in transformations
    }

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        h = (k - 1) // 2
        r = s - h
        d = ell - s

        print()
        print(
            f"{name}: k={k} ell={ell} "
            f"s={s} expected={expected}"
        )

        for label, transform in transformations:

            ndeg, j = transform(
                r,
                d,
            )

            value = get_coeff(
                coeffs,
                j,
                ndeg,
            )

            print(
                f"  {label:10s} "
                f"[N^{ndeg}X^{j}]={value}"
            )

            if value == expected:
                counts[label] += 1

    print()
    print("TRANSFORMATION MATCH COUNTS")

    for label in counts:
        print(
            f"  {label:10s}: "
            f"{counts[label]}/9"
        )


# =============================================================================
# FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Experiment 206 showed that the direct coordinate"
    )

    print(
        "    [N^(s-1) (S+1)^(ell-s+1)]"
    )

    print(
        "is not stable across the anchors."
    )

    print()
    print(
        "Experiment 207 therefore tests the first natural"
    )

    print(
        "reindexing forced by odd k:"
    )

    print(
        "    h = (k-1)/2"
    )

    print(
        "    r = s-h"
    )

    print(
        "together with the boundary distance"
    )

    print(
        "    d = ell-s."
    )

    print()
    print(
        "The purpose is not to discover an arbitrary formula."
    )

    print(
        "It is to determine whether the theorem index s is"
    )

    print(
        "actually an effective post-boundary index after removing"
    )

    print(
        "the compulsory k-dependent half-width."
    )

    print()
    print(
        "A convincing result requires:"
    )

    print(
        "  1. repeated matches across many anchors;"
    )

    print(
        "  2. consistency when different (k,ell,s) produce"
    )

    print(
        "     the same effective coordinates;"
    )

    print(
        "  3. agreement on both parity branches;"
    )

    print(
        "  4. a derivation from the exact pq kernel."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 207")
    print("EFFECTIVE-INDEX / BOUNDARY-DISTANCE AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    index_table()
    anchor_audit()
    collision_audit()
    boundary_collision_audit()
    parity_audit()
    local_transform_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 207")
    print("=" * 78)


if __name__ == "__main__":
    main()

