import sympy as sp

# ==============================================================================
# EXPERIMENT 240
# EXACT pq KERNEL -> FULL L0/L1/L2/L3 RECONSTRUCTION AUDIT
# ==============================================================================
#
# Standalone main.py
# One file only
# Exact arithmetic over QQ
# The pq-kernel is explicitly defined here
#
# Main goals:
#   1. Verify the proposed exact_F kernel directly.
#   2. Construct G(N,X) with X = p+q+1 and N = pq.
#   3. Extract homogeneous layers L0,L1,L2,L3.
#   4. Verify all previously established interior formulas.
#   5. Verify known L2/L3 boundary corrections.
#   6. Localize the k=13 anomaly.
#
# No fitting of new formulas is performed.
# No L4 analysis is performed.
# ==============================================================================


p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
K, L = sp.symbols("K L")
A = sp.symbols("A")


# ==============================================================================
# EXACT pq KERNEL
# ==============================================================================

def exact_F(k, ell):
    k = int(k)
    ell = int(ell)

    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ==============================================================================
# BASIC SYMMETRIC RECONSTRUCTION
# ==============================================================================

# p+q = X-1
# pq   = N
#
# Any symmetric polynomial in p,q can therefore be rewritten in N,X.
#
# We use the recurrence for power sums:
#   S_0 = 2
#   S_1 = X-1
#   S_m = (X-1) S_{m-1} - N S_{m-2}
#
# A monomial p^i q^j + p^j q^i becomes
#   N^min(i,j) * S_|i-j|
#


def power_sum(m):
    m = int(m)

    if m == 0:
        return sp.Integer(2)

    if m == 1:
        return X - 1

    s0 = sp.Integer(2)
    s1 = X - 1

    for _ in range(2, m + 1):
        s2 = sp.expand((X - 1) * s1 - N * s0)
        s0, s1 = s1, s2

    return sp.expand(s1)


def symmetric_pq_to_NX(expr):
    expr = sp.Poly(sp.expand(expr), p, q)

    terms = {}

    for (ip, iq), coeff in expr.terms():
        if ip < iq:
            continue

        if ip == iq:
            key = (ip, 0)
            terms[key] = terms.get(key, sp.Integer(0)) + coeff
            continue

        # Because the input is symmetric, coefficients of (ip,iq)
        # and (iq,ip) must agree.
        c1 = coeff
        c2 = expr.coeff_monomial(p**iq * q**ip)

        if sp.simplify(c1 - c2) != 0:
            raise ValueError(
                f"Input is not symmetric at exponents ({ip},{iq})."
            )

        d = ip - iq
        key = (iq, d)
        terms[key] = terms.get(key, sp.Integer(0)) + c1

    out = sp.Integer(0)

    for (base, diff), coeff in terms.items():
        if diff == 0:
            out += coeff * N**base
        else:
            out += coeff * N**base * power_sum(diff)

    return sp.expand(out)


# ==============================================================================
# EXACT G(N,X)
# ==============================================================================

def exact_G(k, ell):
    F = exact_F(k, ell)
    G = symmetric_pq_to_NX(F)
    return sp.expand(G)


# ==============================================================================
# HOMOGENEOUS LAYER EXTRACTION
# ==============================================================================

def monomial_degree_NX(term):
    powers = term.as_powers_dict()
    return int(powers.get(N, 0) + powers.get(X, 0))


def homogeneous_component(expr, degree):
    out = sp.Integer(0)

    for term in sp.Poly(sp.expand(expr), N, X).terms():
        (a, b), coeff = term

        if a + b == degree:
            out += coeff * N**a * X**b

    return sp.expand(out)


def coefficient(expr, a, b):
    a = int(a)
    b = int(b)

    poly = sp.Poly(sp.expand(expr), N, X)
    return sp.expand(poly.coeff_monomial(N**a * X**b))


def exact_layers(k, ell):
    G = exact_G(k, ell)

    L0 = homogeneous_component(G, ell)
    L1 = homogeneous_component(G, ell - 1)
    L2 = homogeneous_component(G, ell - 2)
    L3 = homogeneous_component(G, ell - 3)

    return G, L0, L1, L2, L3


# ==============================================================================
# KNOWN INTERIOR PRODUCT LAW
# ==============================================================================

def P_r(k, a, ell, r):
    k = sp.Integer(k)
    a = sp.Integer(a)
    ell = sp.Integer(ell)
    r = int(r)

    out = (
        (-1)**(r + 1)
        * sp.binomial(k + r, a)
        / sp.factorial(r)
    )

    for j in range(1, r):
        out *= ell - a - j

    out *= ((k + r) * ell - k * a) / (k + r)

    return sp.factor(out)


# ==============================================================================
# ESTABLISHED L2 BOUNDARY CORRECTIONS
# ==============================================================================

def delta_L2_k(k):
    k = sp.Integer(k)
    return sp.binomial(k + 2, 3)


def delta_L2_k1(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        -sp.binomial(k + 2, 2) * ell
        + k * (k + 2) * (k + 3) / 2
    )


def delta_L2_k2_candidate_low_k(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        (ell - k - 4)
        * ((2 * k - 1) * ell - k * (k + 3))
        / 2
    )


# ==============================================================================
# ESTABLISHED L3 BOUNDARY CORRECTIONS
# ==============================================================================

def delta_L3_k(k):
    k = sp.Integer(k)
    return sp.binomial(k + 2, 3)


def delta_L3_k1(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        -sp.binomial(k + 2, 2) * ell
        + k * (k + 2) * (k + 3) / 2
    )


def delta_L3_k2_low_k(k, ell):
    k = sp.Integer(k)
    ell = sp.Integer(ell)

    return sp.expand(
        (ell - k - 4)
        * ((2 * k - 1) * ell - k * (k + 3))
        / 2
    )


# ==============================================================================
# LAYER COEFFICIENT HELPERS
# ==============================================================================

def layer_coeff(layer, a, degree):
    return coefficient(layer, a, degree - a)


# ==============================================================================
# 1. KERNEL SANITY AUDIT
# ==============================================================================

def kernel_sanity_audit():
    print("=" * 78)
    print("1. EXACT pq KERNEL SANITY AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 4),
        (3, 7),
        (3, 9),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:
        F = exact_F(k, ell)

        swapped = sp.expand(F.xreplace({p: q, q: p}))
        antisym = sp.expand(F + swapped)

        diagonal = sp.expand(F.subs(q, p))

        ok = (
            sp.simplify(antisym) == 0
            and sp.simplify(diagonal) == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"symmetric=True zero_on_p=q=True "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"kernel sanity failures = {failures}")
    print()


# ==============================================================================
# 2. TOP LAYER AUDIT
# ==============================================================================

def top_layer_audit():
    print("=" * 78)
    print("2. TOP-LAYER AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 10),
        (3, 7),
        (3, 20),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
        (13, 25),
    ]

    failures = 0

    for k, ell in cases:
        _, L0, _, _, _ = exact_layers(k, ell)

        expected = sp.expand(
            -X**(ell - k)
            * ((X + N)**k - N**k)
        )

        residual = sp.expand(L0 - expected)
        ok = sp.simplify(residual) == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"residual={residual} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"top-layer failures = {failures}")
    print()


# ==============================================================================
# 3. INTERIOR L1/L2/L3 PRODUCT-LAW AUDIT
# ==============================================================================

def interior_product_law_audit():
    print("=" * 78)
    print("3. INTERIOR L1/L2/L3 PRODUCT-LAW AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),
        (7, 15),
        (7, 17),
        (9, 21),
        (9, 23),
        (11, 23),
    ]

    failures = 0
    tested = 0

    for k, ell in cases:
        _, _, L1, L2, L3 = exact_layers(k, ell)

        local_failures = 0

        # r=1 corresponds to L1, degree ell-1
        # r=2 corresponds to L2, degree ell-2
        # r=3 corresponds to L3, degree ell-3
        for r, layer in [(1, L1), (2, L2), (3, L3)]:
            for a in range(k):
                actual = layer_coeff(layer, a, ell - r)
                expected = sp.factor(P_r(k, a, ell, r))
                residual = sp.simplify(actual - expected)

                tested += 1

                if residual != 0:
                    local_failures += 1
                    failures += 1

                    print(
                        f"FAIL k={k} ell={ell} r={r} a={a}"
                    )
                    print(f"  actual   = {actual}")
                    print(f"  expected = {expected}")
                    print(f"  residual = {residual}")

        print(
            f"k={k:2d} ell={ell:2d} "
            f"interior failures={local_failures} "
            f"{'PASS' if local_failures == 0 else 'FAIL'}"
        )

    print()
    print(f"tested interior coefficients = {tested}")
    print(f"interior formula failures = {failures}")
    print()


# ==============================================================================
# 4. L2 BOUNDARY LOCALIZATION
# ==============================================================================

def l2_boundary_localization():
    print("=" * 78)
    print("4. L2 BOUNDARY LOCALIZATION")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (5, 13),
        (5, 15),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, L2, _ = exact_layers(k, ell)

        rows = {
            "a=k": k,
            "a=k+1": k + 1,
        }

        print(f"k={k:2d} ell={ell:2d}")

        for label, a in rows.items():
            actual = layer_coeff(L2, a, ell - 2)
            interior = P_r(k, a, ell, 2)

            if a == k:
                correction = delta_L2_k(k)
            else:
                correction = delta_L2_k1(k, ell)

            predicted = sp.expand(interior + correction)
            residual = sp.simplify(actual - predicted)

            print(
                f"  {label:5s} "
                f"actual={actual} "
                f"predicted={predicted} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

            if residual != 0:
                failures += 1

    print()
    print(f"L2 boundary failures = {failures}")
    print()


# ==============================================================================
# 5. L3 KNOWN BOUNDARY AUDIT
# ==============================================================================

def l3_known_boundary_audit():
    print("=" * 78)
    print("5. L3 KNOWN BOUNDARY AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (5, 13),
        (5, 15),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, _, L3 = exact_layers(k, ell)

        print(f"k={k:2d} ell={ell:2d}")

        for label, a, correction in [
            ("a=k", k, delta_L3_k(k)),
            ("a=k+1", k + 1, delta_L3_k1(k, ell)),
        ]:
            actual = layer_coeff(L3, a, ell - 3)
            interior = P_r(k, a, ell, 3)
            predicted = sp.expand(interior + correction)
            residual = sp.simplify(actual - predicted)

            print(
                f"  {label:6s} "
                f"actual={actual} "
                f"predicted={predicted} "
                f"residual={residual} "
                f"{'PASS' if residual == 0 else 'FAIL'}"
            )

            if residual != 0:
                failures += 1

    print()
    print(f"L3 known-boundary failures = {failures}")
    print()


# ==============================================================================
# 6. L3 a=k+2 ANOMALY LOCALIZATION
# ==============================================================================

def l3_k2_localization():
    print("=" * 78)
    print("6. L3 a=k+2 ANOMALY LOCALIZATION")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 13),
        (5, 15),
        (5, 17),
        (7, 17),
        (7, 19),
        (9, 21),
        (9, 23),
        (11, 23),
        (11, 25),
        (13, 25),
        (13, 27),
    ]

    print(
        "Columns: k ell actual interior low-k-candidate anomaly"
    )

    for k, ell in cases:
        _, _, _, _, L3 = exact_layers(k, ell)

        a = k + 2

        actual = layer_coeff(L3, a, ell - 3)
        interior = P_r(k, a, ell, 3)
        candidate = sp.expand(
            interior + delta_L3_k2_low_k(k, ell)
        )
        anomaly = sp.expand(actual - candidate)

        print(
            f"k={k:2d} ell={ell:2d} "
            f"actual={str(actual):>8} "
            f"interior={str(interior):>8} "
            f"candidate={str(candidate):>8} "
            f"anomaly={str(anomaly):>8}"
        )

    print()


# ==============================================================================
# 7. k=13 LOCAL ROW AUDIT
# ==============================================================================

def k13_local_rows():
    print("=" * 78)
    print("7. k=13 LOCAL L2/L3 ROW AUDIT")
    print("=" * 78)

    for ell in [25, 27, 29]:
        _, _, _, L2, L3 = exact_layers(13, ell)

        print(f"k=13 ell={ell}")

        for layer_name, layer, degree in [
            ("L2", L2, ell - 2),
            ("L3", L3, ell - 3),
        ]:
            print(f"  {layer_name}")

            for a in [13, 14, 15]:
                actual = layer_coeff(layer, a, degree)
                print(
                    f"    a={a:2d} "
                    f"[N^{a}X^{degree-a}] = {actual}"
                )

    print()


# ==============================================================================
# 8. FULL LOCAL RECONSTRUCTION FOR k=13
# ==============================================================================

def k13_reconstruction():
    print("=" * 78)
    print("8. k=13 COMPLETE LOCAL RECONSTRUCTION")
    print("=" * 78)

    failures = 0

    for ell in [25, 27, 29]:
        _, _, _, L2, L3 = exact_layers(13, ell)

        print(f"k=13 ell={ell}")

        # L2:
        for label, a, correction in [
            ("L2 a=k", 13, delta_L2_k(13)),
            ("L2 a=k+1", 14, delta_L2_k1(13, ell)),
        ]:
            actual = layer_coeff(L2, a, ell - 2)
            predicted = sp.expand(
                P_r(13, a, ell, 2) + correction
            )
            residual = sp.simplify(actual - predicted)

            print(
                f"  {label:10s} "
                f"actual={actual} "
                f"predicted={predicted} "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

        # L3:
        for label, a, correction in [
            ("L3 a=k", 13, delta_L3_k(13)),
            ("L3 a=k+1", 14, delta_L3_k1(13, ell)),
        ]:
            actual = layer_coeff(L3, a, ell - 3)
            predicted = sp.expand(
                P_r(13, a, ell, 3) + correction
            )
            residual = sp.simplify(actual - predicted)

            print(
                f"  {label:10s} "
                f"actual={actual} "
                f"predicted={predicted} "
                f"residual={residual}"
            )

            if residual != 0:
                failures += 1

    print()
    print(f"k=13 known-formula local failures = {failures}")
    print()


# ==============================================================================
# 9. EXACT KERNEL RECONSTRUCTION SPOT CHECKS
# ==============================================================================

def direct_kernel_spot_checks():
    print("=" * 78)
    print("9. DIRECT KERNEL SPOT CHECKS")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (7, 15),
        (9, 21),
    ]

    failures = 0

    for k, ell in cases:
        F = exact_F(k, ell)
        G = exact_G(k, ell)

        # Verify that converting G back to p,q gives F.
        back = sp.expand(
            G.subs({
                N: p * q,
                X: p + q + 1,
            })
        )

        residual = sp.expand(F - back)
        ok = residual == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction residual={residual} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"kernel reconstruction failures = {failures}")
    print()


# ==============================================================================
# 10. FINAL DIAGNOSTIC
# ==============================================================================

def final_diagnostic():
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment uses the proposed exact pq kernel directly:

    F = p^k(1+q)^ell
      + q^k(1+p)^ell
      - p^ell(1+q)^k
      - q^ell(1+p)^k

The purpose is no longer to infer the kernel.

The decisive hierarchy is:

    exact_F
        ->
    exact_G(N,X)
        ->
    L0,L1,L2,L3
        ->
    interior product law
        ->
    boundary corrections.

The k=13 issue is tested locally rather than through
interpolation.

Interpretation:

  * If kernel reconstruction and all established layers pass,
    the pq kernel is strongly confirmed.

  * If L2 passes but L3 a=k+2 fails,
    the anomaly is genuinely a third-layer boundary effect.

  * If L2 fails at k=13,
    the discrepancy enters before L3.

  * If the a=k and a=k+1 rows pass while a=k+2 fails,
    the boundary transition is localized exactly at a=k+2.

No L4 analysis is performed.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    kernel_sanity_audit()
    top_layer_audit()
    interior_product_law_audit()
    l2_boundary_localization()
    l3_known_boundary_audit()
    l3_k2_localization()
    k13_local_rows()
    k13_reconstruction()
    direct_kernel_spot_checks()
    final_diagnostic()


if __name__ == "__main__":
    main()

