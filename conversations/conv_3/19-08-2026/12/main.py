import sympy as sp


# =============================================================================
# EXPERIMENT 389
# EXACT INSTANCE-LOCAL TRACE / NORM REDUCTION
# =============================================================================

N, K, t = sp.symbols("N K t")


# =============================================================================
# TEST INSTANCES
# =============================================================================

INSTANCES = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
]


# =============================================================================
# LAYER POLYNOMIAL
# =============================================================================

def build_R():

    h16 = (
        -9*t**8
        -36*t**7
        -84*t**6
        -126*t**5
        -126*t**4
        -84*t**3
        -36*t**2
        -9*t
        -1
    )

    h15 = (
        88*t**9
        +396*t**8
        +1164*t**7
        +2226*t**6
        +2898*t**5
        +2604*t**4
        +1596*t**3
        +639*t**2
        +151*t
        +16
    )

    h14 = (
        -276*t**10
        -1380*t**9
        -5460*t**8
        -13560*t**7
        -23058*t**6
        -27510*t**5
        -23100*t**4
        -13410*t**3
        -5135*t**2
        -1169*t
        -120
    )

    numerator = sp.expand(
        h15**2 - h16*h14
    )

    denominator = sp.expand(
        h16**2
    )

    return numerator, denominator


R_NUM, R_DEN = build_R()


# =============================================================================
# EXACT KAPPA VALUE
# =============================================================================

def kappa_value(p, q):

    return (
        1
        - (p**2 - p + 1)
        * (q**2 - q + 1)
    )


# =============================================================================
# EXACT TRUE X / TRUE T
# =============================================================================

def true_values(p, q):

    n = sp.Integer(p * q)

    x_true = sp.Integer(
        p + q + 1
    )

    k_true = sp.Integer(
        kappa_value(p, q)
    )

    t_true = sp.cancel(
        n / x_true
    )

    return n, x_true, k_true, t_true


# =============================================================================
# KAPPA QUADRATIC
# =============================================================================

def kappa_quadratic_numeric(n, k):

    return sp.Poly(
        sp.expand(
            (k + n**2 + 2) * t**2
            - n * (n + 3) * t
            + n**2
        ),
        t,
    )


# =============================================================================
# OTHER KAPPA ROOT
# =============================================================================

def other_kappa_root(
    n,
    k,
    t_true,
):

    a = k + n**2 + 2
    c = n**2

    t_other = sp.cancel(
        c / (a * t_true)
    )

    return t_other


# =============================================================================
# EXACT LAYER VALUE
# =============================================================================

def R_value(t_value):

    return sp.cancel(
        R_NUM.subs(t, t_value)
        / R_DEN.subs(t, t_value)
    )


# =============================================================================
# REDUCE NUMERATOR / DENOMINATOR MODULO THE NUMERIC KAPPA QUADRATIC
#
# IMPORTANT:
# Q is already specialized to concrete integers.
#
# This is tiny compared with symbolic reduction over Q(N,K).
# =============================================================================

def reduce_numeric(expr, Q):

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    rem = poly.rem(Q)

    return sp.expand(rem.as_expr())


# =============================================================================
# LINEAR REPRESENTATION
#
# expr = a*t + b
# =============================================================================

def linear_coefficients(expr):

    poly = sp.Poly(
        expr,
        t,
        domain=sp.QQ,
    )

    return (
        poly.coeff_monomial(t),
        poly.coeff_monomial(1),
    )


# =============================================================================
# INSTANCE-LOCAL TRACE / NORM
#
# After reduction:
#
#   R(t) = (u1*t + u0)/(v1*t + v0)
#
# and Q(t) = a*t^2 + b*t + c.
#
# For the two roots:
#
#   S = t1+t2 = -b/a
#   P = t1*t2 = c/a
#
# We compute:
#
#   Trace(R) = R1 + R2
#   Norm(R)  = R1*R2
#
# using exact integer/rational arithmetic.
# =============================================================================

def trace_norm_instance(
    n,
    k,
):

    Q = kappa_quadratic_numeric(
        n,
        k,
    )

    reduced_num = reduce_numeric(
        R_NUM,
        Q,
    )

    reduced_den = reduce_numeric(
        R_DEN,
        Q,
    )

    u1, u0 = linear_coefficients(
        reduced_num
    )

    v1, v0 = linear_coefficients(
        reduced_den
    )

    a = sp.Integer(
        k + n**2 + 2
    )

    b = sp.Integer(
        -n * (n + 3)
    )

    c = sp.Integer(
        n**2
    )

    S = sp.cancel(
        -b / a
    )

    P = sp.cancel(
        c / a
    )

    denominator_product = sp.cancel(
        v1**2 * P
        + v1*v0*S
        + v0**2
    )

    trace = sp.cancel(
        (
            2*u1*v1*P
            + (u1*v0 + u0*v1)*S
            + 2*u0*v0
        )
        / denominator_product
    )

    norm = sp.cancel(
        (
            u1**2*P
            + u1*u0*S
            + u0**2
        )
        / denominator_product
    )

    return (
        Q,
        reduced_num,
        reduced_den,
        trace,
        norm,
    )


# =============================================================================
# FACTOR SIGNATURE
# =============================================================================

def factor_signature(value):

    value = sp.factor(
        sp.cancel(value)
    )

    num, den = sp.fraction(value)

    num_poly = sp.Poly(
        num,
        sp.Symbol("dummy"),
    ) if False else None

    num_fac = sp.factorint(
        abs(int(sp.numer(value)))
    )

    den_fac = sp.factorint(
        abs(int(sp.denom(value)))
    )

    return (
        len(num_fac),
        len(den_fac),
        max(num_fac.keys()) if num_fac else 1,
        max(den_fac.keys()) if den_fac else 1,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 110)
    print("EXPERIMENT 389 START")
    print("=" * 110)
    print()
    print("EXACT INSTANCE-LOCAL TRACE / NORM REDUCTION")
    print()
    print("No symbolic N,K rational-function reduction is performed.")
    print("Every instance is specialized before reduction.")
    print()

    global_trace_pass = True
    global_norm_pass = True
    global_other_root_pass = True
    global_reduction_pass = True

    trace_sizes = []
    norm_sizes = []

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        n, x_true, k_true, t_true = true_values(
            p,
            q,
        )

        t_other = other_kappa_root(
            n,
            k_true,
            t_true,
        )

        # -------------------------------------------------------------
        # Verify both Kappa roots.
        # -------------------------------------------------------------

        Q = kappa_quadratic_numeric(
            n,
            k_true,
        )

        q_true = sp.expand(
            Q.as_expr().subs(
                t,
                t_true,
            )
        )

        q_other = sp.expand(
            Q.as_expr().subs(
                t,
                t_other,
            )
        )

        true_root_pass = (
            q_true == 0
        )

        other_root_pass = (
            q_other == 0
        )

        # -------------------------------------------------------------
        # Direct branches.
        # -------------------------------------------------------------

        r_true = R_value(
            t_true
        )

        r_other = R_value(
            t_other
        )

        direct_trace = sp.cancel(
            r_true + r_other
        )

        direct_norm = sp.cancel(
            r_true * r_other
        )

        # -------------------------------------------------------------
        # Instance-local reduction.
        # -------------------------------------------------------------

        (
            Q_used,
            reduced_num,
            reduced_den,
            trace,
            norm,
        ) = trace_norm_instance(
            n,
            k_true,
        )

        trace_pass = (
            sp.cancel(trace - direct_trace)
            == 0
        )

        norm_pass = (
            sp.cancel(norm - direct_norm)
            == 0
        )

        # -------------------------------------------------------------
        # Check reduced representation itself.
        #
        # R(t) - reduced fraction must vanish modulo Q.
        # -------------------------------------------------------------

        reduced_difference = sp.cancel(
            (
                R_NUM
                * reduced_den
                - R_DEN
                * reduced_num
            )
        )

        reduction_remainder = sp.rem(
            sp.Poly(
                sp.expand(
                    reduced_difference
                ),
                t,
                domain=sp.QQ,
            ),
            Q,
        )

        reduction_pass = (
            reduction_remainder.as_expr()
            == 0
        )

        # -------------------------------------------------------------
        # Sizes.
        # -------------------------------------------------------------

        trace_num = sp.numer(trace)
        trace_den = sp.denom(trace)

        norm_num = sp.numer(norm)
        norm_den = sp.denom(norm)

        trace_size = (
            len(str(abs(int(trace_num)))),
            len(str(abs(int(trace_den)))),
        )

        norm_size = (
            len(str(abs(int(norm_num)))),
            len(str(abs(int(norm_den)))),
        )

        trace_sizes.append(
            trace_size
        )

        norm_sizes.append(
            norm_size
        )

        # -------------------------------------------------------------
        # Output.
        # -------------------------------------------------------------

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {n}")
        print(f"  X_TRUE  = {x_true}")
        print(f"  K_TRUE  = {k_true}")
        print(f"  t_TRUE  = {t_true}")
        print(f"  t_OTHER = {t_other}")
        print()

        print("KAPPA ROOTS")
        print(
            f"  t_TRUE root  : {true_root_pass}"
        )
        print(
            f"  t_OTHER root : {other_root_pass}"
        )
        print()

        print("INSTANCE-LOCAL REDUCTION")
        print(
            "  reduced numerator degree = "
            f"{sp.degree(reduced_num, t)}"
        )
        print(
            "  reduced denominator degree = "
            f"{sp.degree(reduced_den, t)}"
        )
        print(
            "  reduced numerator terms = "
            f"{len(sp.Poly(reduced_num, t).terms())}"
        )
        print(
            "  reduced denominator terms = "
            f"{len(sp.Poly(reduced_den, t).terms())}"
        )
        print()

        print("REDUCTION IDENTITY")
        print(
            f"  R(t) modulo Q reconstruction : "
            f"{reduction_pass}"
        )
        print()

        print("TRACE")
        print(
            f"  direct R1+R2 == trace : "
            f"{trace_pass}"
        )
        print(
            f"  numerator digits   = {trace_size[0]}"
        )
        print(
            f"  denominator digits = {trace_size[1]}"
        )
        print()

        print("NORM")
        print(
            f"  direct R1*R2 == norm : "
            f"{norm_pass}"
        )
        print(
            f"  numerator digits   = {norm_size[0]}"
        )
        print(
            f"  denominator digits = {norm_size[1]}"
        )
        print()

        print("BRANCH VALUES")
        print(
            f"  R_TRUE  numerator digits = "
            f"{len(str(abs(int(sp.numer(r_true)))))}"
        )
        print(
            f"  R_TRUE  denominator digits = "
            f"{len(str(abs(int(sp.denom(r_true)))))}"
        )
        print(
            f"  R_OTHER numerator digits = "
            f"{len(str(abs(int(sp.numer(r_other)))))}"
        )
        print(
            f"  R_OTHER denominator digits = "
            f"{len(str(abs(int(sp.denom(r_other)))))}"
        )
        print()

        # Keep expressions suppressed.
        print("STRUCTURAL VALUES")
        print(
            f"  trace expression size = "
            f"{len(str(trace))}"
        )
        print(
            f"  norm expression size = "
            f"{len(str(norm))}"
        )
        print()

        global_trace_pass &= trace_pass
        global_norm_pass &= norm_pass
        global_other_root_pass &= other_root_pass
        global_reduction_pass &= reduction_pass

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all Kappa roots pass       = "
        f"{global_other_root_pass}"
    )

    print(
        "  all reductions pass        = "
        f"{global_reduction_pass}"
    )

    print(
        "  all Trace identities pass  = "
        f"{global_trace_pass}"
    )

    print(
        "  all Norm identities pass   = "
        f"{global_norm_pass}"
    )

    print()

    print("TRACE SIZE RANGE")

    print(
        f"  minimum numerator digits = "
        f"{min(x[0] for x in trace_sizes)}"
    )

    print(
        f"  maximum numerator digits = "
        f"{max(x[0] for x in trace_sizes)}"
    )

    print(
        f"  minimum denominator digits = "
        f"{min(x[1] for x in trace_sizes)}"
    )

    print(
        f"  maximum denominator digits = "
        f"{max(x[1] for x in trace_sizes)}"
    )

    print()

    print("NORM SIZE RANGE")

    print(
        f"  minimum numerator digits = "
        f"{min(x[0] for x in norm_sizes)}"
    )

    print(
        f"  maximum numerator digits = "
        f"{max(x[0] for x in norm_sizes)}"
    )

    print(
        f"  minimum denominator digits = "
        f"{min(x[1] for x in norm_sizes)}"
    )

    print(
        f"  maximum denominator digits = "
        f"{max(x[1] for x in norm_sizes)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 389 FINAL STATUS")
    print("=" * 110)
    print()

    print(
        "This experiment deliberately avoids the symbolic "
        "Q(N,K) computation."
    )

    print()
    print(
        "The central objects are:"
    )
    print(
        "  Trace(R) = R(t_TRUE) + R(t_OTHER)"
    )
    print(
        "  Norm(R)  = R(t_TRUE) * R(t_OTHER)"
    )

    print()
    print(
        "If these remain much smaller than the full resultant, "
        "the trace/norm pair is the better object for structural "
        "elimination."
    )

    print()
    print("=" * 110)
    print("EXPERIMENT 389 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
