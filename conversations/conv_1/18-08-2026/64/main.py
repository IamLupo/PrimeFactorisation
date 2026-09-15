import random
import sympy as sp

# ============================================================
# EXPERIMENT 67
# RANDOM PRIME PAIR / n=pq / EXACT ALL HOMOGENEOUS LAYERS
#
# This is a standalone replacement.
# No external experiment script is required.
# ============================================================

RNG_SEED = None          # Set an integer for reproducibility, or None.
P_MIN = 100_000
P_MAX = 1_000_000

K = 9
ELL = 16

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")


# ============================================================
# 1. EXACT SUPPLIED KERNEL
# ============================================================

def exact_F_9_16(pv, qv):
    """
    Exact supplied k=9, ell=16 kernel.

    Accepts either symbolic SymPy expressions or integers.
    No floating-point coercion is performed.
    """
    pv = sp.sympify(pv)
    qv = sp.sympify(qv)

    expr = (
        -9*pv**16*qv**8
        -36*pv**16*qv**7
        -84*pv**16*qv**6
        -126*pv**16*qv**5
        -126*pv**16*qv**4
        -84*pv**16*qv**3
        -36*pv**16*qv**2
        -9*pv**16*qv
        -pv**16

        +16*pv**15*qv**9

        +120*pv**14*qv**9
        +560*pv**13*qv**9
        +1820*pv**12*qv**9
        +4368*pv**11*qv**9
        +8008*pv**10*qv**9

        +16*pv**9*qv**15
        +120*pv**9*qv**14
        +560*pv**9*qv**13
        +1820*pv**9*qv**12
        +4368*pv**9*qv**11
        +8008*pv**9*qv**10
        +22880*pv**9*qv**9
        +12870*pv**9*qv**8
        +11440*pv**9*qv**7
        +8008*pv**9*qv**6
        +4368*pv**9*qv**5
        +1820*pv**9*qv**4
        +560*pv**9*qv**3
        +120*pv**9*qv**2
        +16*pv**9*qv
        +pv**9

        -9*pv**8*qv**16
        +12870*pv**8*qv**9

        -36*pv**7*qv**16
        +11440*pv**7*qv**9

        -84*pv**6*qv**16
        +8008*pv**6*qv**9

        -126*pv**5*qv**16
        +4368*pv**5*qv**9

        -126*pv**4*qv**16
        +1820*pv**4*qv**9

        -84*pv**3*qv**16
        +560*pv**3*qv**9

        -36*pv**2*qv**16
        +120*pv**2*qv**9

        -9*pv*qv**16
        +16*pv*qv**9

        -qv**16
        +qv**9
    )

    return sp.expand(expr)


# ============================================================
# 2. RANDOM PRIME GENERATION
# ============================================================

def random_prime(rng, lo, hi):
    while True:
        candidate = rng.randint(lo, hi)
        if sp.isprime(candidate):
            return candidate


# ============================================================
# 3. SYMMETRIC REPRESENTATION
# ============================================================

def symmetric_kernel_from_F(F):
    """
    Rewrite F(p,q) exactly using

        s = p + q
        n = p q

    and then substitute

        s = X - 1.

    Returns G(N,X).
    """

    result = sp.symmetrize(
        sp.expand(F),
        [p, q],
        formal=True,
    )

    if len(result) != 3:
        raise RuntimeError(
            f"Unexpected symmetrize result: {result}"
        )

    symmetric_part, remainder, mapping = result

    remainder = sp.expand(remainder)

    if remainder != 0:
        raise RuntimeError(
            "Kernel is not symmetric in p and q.\n"
            f"Symmetry remainder = {remainder}"
        )

    s = sp.Symbol("s")
    n = sp.Symbol("n")

    substitutions = {}

    for generated, original in mapping:
        original = sp.expand(original)

        if sp.expand(original - (p + q)) == 0:
            substitutions[generated] = s

        elif sp.expand(original - p*q) == 0:
            substitutions[generated] = n

        else:
            raise RuntimeError(
                "Unexpected symmetric-variable mapping:\n"
                f"  {generated} -> {original}"
            )

    G_s_n = sp.expand(
        symmetric_part.subs(substitutions)
    )

    G_NX = sp.expand(
        G_s_n.subs({
            s: X - 1,
            n: N,
        })
    )

    return G_NX


# ============================================================
# 4. EXACT HOMOGENEOUS DECOMPOSITION
# ============================================================

def homogeneous_layers(expr):
    """
    Split a polynomial in N,X according to total degree.

    For example:

        N^2 X^3 -> total degree 5.

    Returns:

        {
            degree: exact homogeneous polynomial
        }

    """
    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
    )

    layers = {}

    for monomial, coefficient in poly.terms():
        n_power, x_power = monomial
        total_degree = n_power + x_power

        term = (
            coefficient
            * N**n_power
            * X**x_power
        )

        layers[total_degree] = sp.expand(
            layers.get(
                total_degree,
                sp.Integer(0),
            ) + term
        )

    return {
        degree: sp.expand(layers[degree])
        for degree in sorted(
            layers,
            reverse=True,
        )
    }


# ============================================================
# 5. EXACT EVALUATION
# ============================================================

def exact_eval(expr, substitutions):
    """
    Exact evaluation with an explicit float guard.
    """
    value = sp.expand(
        expr.subs(substitutions)
    )

    if value.has(sp.Float):
        raise RuntimeError(
            "Floating-point arithmetic detected."
        )

    return value


def assert_no_float(expr, label):
    """
    Reject any floating-point number anywhere in expr.
    """
    if expr.has(sp.Float):
        raise RuntimeError(
            f"Floating-point coefficient detected in {label}"
        )


# ============================================================
# 6. MAIN EXPERIMENT
# ============================================================

def main():

    rng = random.Random(RNG_SEED)

    print("=" * 78)
    print(
        "EXPERIMENT 67 — RANDOM PRIME PAIR / n=pq / "
        "EXACT ALL HOMOGENEOUS LAYERS"
    )
    print("=" * 78)

    # --------------------------------------------------------
    # 0. GENERATED PRIME PAIR
    # --------------------------------------------------------

    p_value = random_prime(
        rng,
        P_MIN,
        P_MAX,
    )

    q_value = random_prime(
        rng,
        P_MIN,
        P_MAX,
    )

    while q_value == p_value:
        q_value = random_prime(
            rng,
            P_MIN,
            P_MAX,
        )

    n_value = p_value * q_value
    S_value = p_value + q_value
    X_value = S_value + 1
    phi_value = (
        p_value - 1
    ) * (
        q_value - 1
    )

    print("\n0. GENERATED PRIME PAIR")
    print(f"  p = {p_value}")
    print(f"  q = {q_value}")
    print(f"  p != q = {p_value != q_value}")
    print(f"  k = {K}")
    print(f"  ell = {ELL}")

    # --------------------------------------------------------
    # 1. SEMIPRIME CONSTRUCTION
    # --------------------------------------------------------

    print("\n1. SEMIPRIME CONSTRUCTION")
    print(f"  n = p*q = {n_value}")
    print(f"  S = p+q = {S_value}")
    print(f"  X = p+q+1 = {X_value}")
    print(
        f"  phi(n) = (p-1)(q-1) = {phi_value}"
    )

    # --------------------------------------------------------
    # 2. EXACT F(p,q)
    # --------------------------------------------------------

    F_symbolic = exact_F_9_16(
        p,
        q,
    )

    assert_no_float(
        F_symbolic,
        "F(p,q)",
    )

    F_value = exact_eval(
        F_symbolic,
        {
            p: p_value,
            q: q_value,
        },
    )

    print("\n2. EXACT KERNEL")
    print(
        "  kernel source = supplied k=9, ell=16 exact kernel"
    )
    print(
        f"  F_{{{K},{ELL}}}(p,q) = {F_symbolic}"
    )
    print(
        f"  F_{{{K},{ELL}}}"
        f"({p_value},{q_value}) = {F_value}"
    )

    # --------------------------------------------------------
    # 3. EXACT G(N,X)
    # --------------------------------------------------------

    G_symbolic = symmetric_kernel_from_F(
        F_symbolic
    )

    assert_no_float(
        G_symbolic,
        "G(N,X)",
    )

    G_value = exact_eval(
        G_symbolic,
        {
            N: n_value,
            X: X_value,
        },
    )

    print("\n3. EXACT SYMMETRIC KERNEL")
    print(
        f"  G_{{{K},{ELL}}}(N,X) = {G_symbolic}"
    )
    print(
        f"  G_{{{K},{ELL}}}"
        f"({n_value},{X_value}) = {G_value}"
    )

    # --------------------------------------------------------
    # 4. ALL HOMOGENEOUS LAYERS
    # --------------------------------------------------------

    layers = homogeneous_layers(
        G_symbolic
    )

    degrees = sorted(
        layers,
        reverse=True,
    )

    max_degree = degrees[0]
    min_degree = degrees[-1]

    print("\n4. ALL EXACT HOMOGENEOUS LAYERS")
    print(
        f"  highest total degree = {max_degree}"
    )
    print(
        f"  lowest total degree = {min_degree}"
    )
    print(
        f"  number of nonzero layers = {len(degrees)}"
    )

    layer_values = {}

    for degree in degrees:

        layer = layers[degree]

        assert_no_float(
            layer,
            f"G^({degree})",
        )

        value = exact_eval(
            layer,
            {
                N: n_value,
                X: X_value,
            },
        )

        layer_values[degree] = value

        print(
            f"\n  degree={degree}: G^({degree})"
        )
        print(
            f"    polynomial = {layer}"
        )
        print(
            f"    value = {value}"
        )

    # --------------------------------------------------------
    # 5. TOP THREE LAYERS
    # --------------------------------------------------------

    top_degrees = degrees[:3]

    top_three_expr = sp.expand(
        sum(
            layers[d]
            for d in top_degrees
        )
    )

    top_three_value = sum(
        layer_values[d]
        for d in top_degrees
    )

    lower_expr = sp.expand(
        G_symbolic
        - top_three_expr
    )

    lower_value = (
        G_value
        - top_three_value
    )

    print("\n5. TOP THREE LAYERS")

    print(
        f"  degrees included = {top_degrees}"
    )

    for degree in top_degrees:
        print(
            f"  G^({degree}) value = "
            f"{layer_values[degree]}"
        )

    print(
        f"  top-3 sum = {top_three_value}"
    )

    print(
        f"  lower-layer residual = {lower_value}"
    )

    print(
        f"  lower-layer residual polynomial = "
        f"{lower_expr}"
    )

    # --------------------------------------------------------
    # 6. EXACT RECONSTRUCTION
    # --------------------------------------------------------

    all_layers_expr = sp.expand(
        sum(
            layers.values()
        )
    )

    all_layers_value = sum(
        layer_values.values()
    )

    lower_layers_expr = sp.expand(
        sum(
            layers[d]
            for d in degrees[3:]
        )
    )

    lower_layers_value = sum(
        layer_values[d]
        for d in degrees[3:]
    )

    print("\n6. RECONSTRUCTION CHECKS")

    print(
        "  F(p,q) == G(pq,p+q+1): "
        f"{F_value == G_value}"
    )

    print(
        "  sum(all homogeneous layers) == G: "
        f"{all_layers_expr == G_symbolic}"
    )

    print(
        "  evaluated sum(all layers) == G value: "
        f"{all_layers_value == G_value}"
    )

    print(
        "  top-3 + lower residual == G: "
        f"{sp.expand(top_three_expr + lower_expr) == G_symbolic}"
    )

    print(
        "  lower residual == sum(all lower layers): "
        f"{sp.expand(lower_expr - lower_layers_expr) == 0}"
    )

    print(
        "  evaluated lower residual == "
        "sum(all lower layers): "
        f"{lower_value == lower_layers_value}"
    )

    # Hard failures.
    if F_value != G_value:
        raise AssertionError(
            "F(p,q) != G(pq,p+q+1)"
        )

    if sp.expand(
        all_layers_expr - G_symbolic
    ) != 0:
        raise AssertionError(
            "Homogeneous layer reconstruction failed"
        )

    if all_layers_value != G_value:
        raise AssertionError(
            "Evaluated homogeneous layer reconstruction failed"
        )

    if sp.expand(
        top_three_expr
        + lower_expr
        - G_symbolic
    ) != 0:
        raise AssertionError(
            "Top-three plus residual reconstruction failed"
        )

    if sp.expand(
        lower_expr
        - lower_layers_expr
    ) != 0:
        raise AssertionError(
            "Lower-layer residual reconstruction failed"
        )

    if lower_value != lower_layers_value:
        raise AssertionError(
            "Evaluated lower-layer reconstruction failed"
        )

    # --------------------------------------------------------
    # 7. DEGREE PROFILE
    # --------------------------------------------------------

    print("\n7. HOMOGENEOUS DEGREE PROFILE")

    for degree in range(
        max_degree,
        min_degree - 1,
        -1,
    ):
        if degree in layers:
            print(
                f"  degree={degree}: PRESENT"
            )
        else:
            print(
                f"  degree={degree}: ZERO"
            )

    # --------------------------------------------------------
    # 8. EXACTNESS CONTROL
    # --------------------------------------------------------

    print("\n8. EXACTNESS CONTROLS")

    contains_float_F = F_symbolic.has(
        sp.Float
    )

    contains_float_G = G_symbolic.has(
        sp.Float
    )

    contains_float_layers = any(
        layer.has(sp.Float)
        for layer in layers.values()
    )

    print(
        f"  F contains Float: "
        f"{contains_float_F}"
    )

    print(
        f"  G contains Float: "
        f"{contains_float_G}"
    )

    print(
        f"  layers contain Float: "
        f"{contains_float_layers}"
    )

    if (
        contains_float_F
        or contains_float_G
        or contains_float_layers
    ):
        raise AssertionError(
            "Floating-point arithmetic entered "
            "the symbolic pipeline"
        )

    # --------------------------------------------------------
    # 9. FINAL NOTE
    # --------------------------------------------------------

    print("\n9. NOTE")
    print(
        "  One fresh semiprime was generated."
    )
    print(
        "  No fitting, regression, ranking, or "
        "statistical analysis was performed."
    )
    print(
        "  All symbolic transformations and "
        "evaluations are exact SymPy arithmetic."
    )
    print(
        "  Floating-point coefficients are explicitly rejected."
    )
    print(
        "  The top three layers are only a partial decomposition."
    )
    print(
        "  The exact residual contains every homogeneous "
        "layer below them."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT COMPLETE — ALL EXACT RECONSTRUCTION CHECKS PASSED")
    print("=" * 78)


if __name__ == "__main__":
    main()