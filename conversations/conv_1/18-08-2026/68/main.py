#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 70 — THREE-LAYER INVERSE RECONSTRUCTION
                 OF N = p*q

Self-contained exact experiment.

No external kernel is imported.

For the verified k=9, ell=16 symmetric kernel,

    G(N,X)

we form

    F(p,q) = G(p*q, p+q+1).

For homogeneous layers,

    G^(d)(N,X) = X^d h_d(N/X),

so the top three layers satisfy

    L16 = X^16 h16(t)
    L15 = X^15 h15(t)
    L14 = X^14 h14(t),

where

    t = N/X.

Eliminating X gives

    L16*L14*h15(t)^2
    -
    L15^2*h16(t)*h14(t)
    = 0.

The experiment then attempts to recover

    t -> X -> N -> S -> Delta -> p,q

using exact rational arithmetic.

The generated p,q are used only as a validation target.
They are NOT used in the inverse calculation itself.

No floats.
No fitting.
No regression.
No interpolation.
==============================================================================
"""

from __future__ import annotations

import math
import random

import sympy as sp
from sympy import Integer, Rational, Poly, factor, expand


# ============================================================================
# CONFIGURATION
# ============================================================================

K = 9
ELL = 16

PRIME_LOW = 50_000
PRIME_HIGH = 300_000

N_SYM, X_SYM = sp.symbols("N X")
P_SYM, Q_SYM = sp.symbols("p q")
T_SYM = sp.symbols("t")


# ============================================================================
# OUTPUT
# ============================================================================

def section(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def check(label: str, value: bool) -> None:
    print(f"  {label}: {bool(value)}")


# ============================================================================
# EXACT SYMMETRIC KERNEL G_{9,16}(N,X)
# ============================================================================

def build_G():
    N = N_SYM
    X = X_SYM

    G = (
        -50*N**12
        + 288*N**11*X**2
        - 300*N**11*X
        + 1300*N**11

        - 276*N**10*X**4
        + 1584*N**10*X**3
        - 7436*N**10*X**2
        + 7150*N**10*X
        - 10010*N**10

        + 88*N**9*X**6
        - 1380*N**9*X**5
        + 10340*N**9*X**4
        - 34430*N**9*X**3
        + 62634*N**9*X**2
        - 50050*N**9*X
        + 35750*N**9

        - 9*N**8*X**8
        + 396*N**8*X**7
        - 5460*N**8*X**6
        + 34650*N**8*X**5
        - 117810*N**8*X**4
        + 228228*N**8*X**3
        - 252252*N**8*X**2
        + 160875*N**8*X
        - 71500*N**8

        - 36*N**7*X**9
        + 1164*N**7*X**8
        - 13560*N**7*X**7
        + 79464*N**7*X**6
        - 267960*N**7*X**5
        + 552552*N**7*X**4
        - 708708*N**7*X**3
        + 563420*N**7*X**2
        - 286000*N**7*X
        + 88400*N**7

        - 84*N**6*X**10
        + 2226*N**6*X**9
        - 23058*N**6*X**8
        + 127512*N**6*X**7
        - 426888*N**6*X**6
        + 918918*N**6*X**5
        - 1303302*N**6*X**4
        + 1221220*N**6*X**3
        - 755664*N**6*X**2
        + 309400*N**6*X
        - 71400*N**6

        - 126*N**5*X**11
        + 2898*N**5*X**10
        - 27510*N**5*X**9
        + 145530*N**5*X**8
        - 484110*N**5*X**7
        + 1075074*N**5*X**6
        - 1639638*N**5*X**5
        + 1732458*N**5*X**4
        - 1265992*N**5*X**3
        + 636888*N**5*X**2
        - 214200*N**5*X
        + 38760*N**5

        - 126*N**4*X**12
        + 2604*N**4*X**11
        - 23100*N**4*X**10
        + 117975*N**4*X**9
        - 390225*N**4*X**8
        + 887172*N**4*X**7
        - 1429428*N**4*X**6
        + 1653470*N**4*X**5
        - 1375360*N**4*X**4
        + 818720*N**4*X**3
        - 344352*N**4*X**2
        + 96900*N**4*X
        - 14250*N**4

        - 84*N**3*X**13
        + 1596*N**3*X**12
        - 13410*N**3*X**11
        + 66550*N**3*X**10
        - 219010*N**3*X**9
        + 507078*N**3*X**8
        - 852852*N**3*X**7
        + 1058540*N**3*X**6
        - 974400*N**3*X**5
        + 663680*N**3*X**4
        - 331704*N**3*X**3
        + 119016*N**3*X**2
        - 28500*N**3*X
        + 3500*N**3

        - 36*N**2*X**14
        + 639*N**2*X**13
        - 5135*N**2*X**12
        + 24882*N**2*X**11
        - 81510*N**2*X**10
        + 191477*N**2*X**9
        - 333333*N**2*X**8
        + 437700*N**2*X**7
        - 436832*N**2*X**6
        + 331500*N**2*X**5
        - 190296*N**2*X**4
        + 81624*N**2*X**3
        - 25380*N**2*X**2
        + 5250*N**2*X
        - 550*N**2

        - 9*N*X**15
        + 151*N*X**14
        - 1169*N*X**13
        + 5551*N*X**12
        - 18109*N*X**11
        + 43043*N*X**10
        - 77077*N*X**9
        + 105979*N*X**8
        - 112964*N*X**7
        + 93604*N*X**6
        - 60144*N*X**5
        + 29736*N*X**4
        - 11130*N*X**3
        + 3038*N*X**2
        - 550*N*X
        + 50*N

        - X**16
        + 16*X**15
        - 120*X**14
        + 560*X**13
        - 1820*X**12
        + 4368*X**11
        - 8008*X**10
        + 11441*X**9
        - 12879*X**8
        + 11476*X**7
        - 8092*X**6
        + 4494*X**5
        - 1946*X**4
        + 644*X**3
        - 156*X**2
        + 25*X
        - 2
    )

    return sp.expand(G)


# ============================================================================
# F(p,q) = G(pq, p+q+1)
# ============================================================================

def build_F(G):
    return sp.expand(
        G.subs(
            {
                N_SYM: P_SYM * Q_SYM,
                X_SYM: P_SYM + Q_SYM + 1,
            }
        )
    )


# ============================================================================
# PRIME GENERATION
# ============================================================================

def random_prime():
    return int(sp.randprime(PRIME_LOW, PRIME_HIGH))


def random_prime_pair():
    while True:
        p = random_prime()
        q = random_prime()
        if p != q:
            return p, q


# ============================================================================
# HOMOGENEOUS LAYERS
# ============================================================================

def homogeneous_layers(G):
    poly = Poly(G, N_SYM, X_SYM, domain="QQ")

    max_degree = max(
        a + b
        for (a, b), coeff in poly.terms()
        if coeff != 0
    )

    min_degree = min(
        a + b
        for (a, b), coeff in poly.terms()
        if coeff != 0
    )

    layers = {}

    for d in range(max_degree, min_degree - 1, -1):
        layer = Integer(0)

        for (a, b), coeff in poly.terms():
            if a + b == d:
                layer += coeff * N_SYM**a * X_SYM**b

        layers[d] = sp.expand(layer)

    return layers


def h_of_layer(layer, degree):
    """
    For homogeneous degree d:

        layer(N,X) = X^d h(N/X)

    hence

        h(t) = layer(t,1).
    """
    return sp.expand(
        layer.subs(
            {
                N_SYM: T_SYM,
                X_SYM: Integer(1),
            }
        )
    )


# ============================================================================
# EXACT INTEGER ROOT / POWER UTILITIES
# ============================================================================

def exact_nth_root_integer(n: int, power: int):
    """
    Return integer r if n = r^power exactly.
    Otherwise return None.
    """
    if n < 0:
        return None

    root = int(round(n ** (1.0 / power)))

    # Correct possible floating rounding locally.
    while pow(root, power) < n:
        root += 1

    while root > 0 and pow(root, power) > n:
        root -= 1

    if pow(root, power) == n:
        return root

    return None


# ============================================================================
# SCALE-FREE INVERSE EQUATION
# ============================================================================

def build_elimination_polynomial(
    L16,
    L15,
    L14,
    h16,
    h15,
    h14,
):
    """
    From

        L16 = X^16 h16(t)
        L15 = X^15 h15(t)
        L14 = X^14 h14(t)

    eliminate X:

        L16*L14*h15(t)^2
        -
        L15^2*h16(t)*h14(t)
        = 0.
    """
    expression = sp.expand(
        Integer(L16) * Integer(L14) * h15**2
        -
        Integer(L15)**2 * h16 * h14
    )

    return sp.Poly(
        expression,
        T_SYM,
        domain="ZZ",
    )


# ============================================================================
# FIND EXACT RECONSTRUCTION CANDIDATES
# ============================================================================

def find_reconstruction_candidates(
    elimination_poly,
    L16,
    h16,
):
    """
    Find rational roots t.

    For each rational root, compute

        X^16 = L16 / h16(t)

    and retain exact positive integer X.
    """
    candidates = []

    roots = sp.solve(
        elimination_poly.as_expr(),
        T_SYM,
    )

    for root in roots:
        root = sp.cancel(root)

        if not root.is_Rational:
            continue

        denominator = sp.cancel(
            h16.subs(T_SYM, root)
        )

        if denominator == 0:
            continue

        X16 = sp.cancel(
            Rational(L16) / denominator
        )

        if not X16.is_Integer:
            continue

        X16_int = int(X16)

        if X16_int <= 0:
            continue

        X_value = exact_nth_root_integer(
            X16_int,
            16,
        )

        if X_value is None:
            continue

        candidates.append(
            (
                root,
                Integer(X_value),
            )
        )

    # Deduplicate exact candidates.
    unique = []
    seen = set()

    for item in candidates:
        key = (item[0], item[1])
        if key not in seen:
            seen.add(key)
            unique.append(item)

    return unique


# ============================================================================
# EXTRACT p,q FROM N,X
# ============================================================================

def recover_factors(N_value, X_value):
    S_value = Integer(X_value - 1)

    Delta_value = sp.expand(
        S_value**2 - 4 * Integer(N_value)
    )

    sqrt_delta = sp.sqrt(Delta_value)

    if not sqrt_delta.is_Integer:
        delta_int = int(Delta_value)

        if delta_int < 0:
            raise RuntimeError(
                "Recovered discriminant is negative."
            )

        root = math.isqrt(delta_int)

        if root * root != delta_int:
            raise RuntimeError(
                "Recovered discriminant is not a perfect square."
            )

        sqrt_delta = Integer(root)

    p_value = sp.cancel(
        (S_value + sqrt_delta) / 2
    )

    q_value = sp.cancel(
        (S_value - sqrt_delta) / 2
    )

    if not p_value.is_Integer or not q_value.is_Integer:
        raise RuntimeError(
            "Recovered factor candidates are not integers."
        )

    return (
        Integer(p_value),
        Integer(q_value),
        S_value,
        Integer(Delta_value),
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    random.seed()

    G = build_G()
    F = build_F(G)

    p, q = random_prime_pair()

    if p > q:
        p, q = q, p

    true_N = p * q
    true_S = p + q
    true_X = true_S + 1
    true_Delta = true_S**2 - 4 * true_N

    # ---------------------------------------------------------------------
    # 0
    # ---------------------------------------------------------------------

    section(
        "EXPERIMENT 70 — THREE-LAYER INVERSE RECONSTRUCTION OF N = p*q"
    )

    print("0. EXACT SYMBOLIC SETUP")
    print(f"  k = {K}")
    print(f"  ell = {ELL}")
    print("  kernel = self-contained exact G_{9,16}(N,X)")
    print("  floating point = forbidden")

    # ---------------------------------------------------------------------
    # 1
    # ---------------------------------------------------------------------

    section("1. GENERATED PRIME PAIR")

    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  p != q = {p != q}")
    print(f"  N = p*q = {true_N}")
    print(f"  S = p+q = {true_S}")
    print(f"  X = S+1 = {true_X}")

    # ---------------------------------------------------------------------
    # 2
    # ---------------------------------------------------------------------

    section("2. EXACT F(p,q) / G(N,X) IDENTITY")

    F_value = sp.expand(
        F.subs(
            {
                P_SYM: Integer(p),
                Q_SYM: Integer(q),
            }
        )
    )

    G_value = sp.expand(
        G.subs(
            {
                N_SYM: Integer(true_N),
                X_SYM: Integer(true_X),
            }
        )
    )

    print(f"  digits in F(p,q) = {len(str(abs(int(F_value))))}")
    check(
        "F(p,q) == G(pq,p+q+1)",
        F_value == G_value,
    )

    # ---------------------------------------------------------------------
    # 3
    # ---------------------------------------------------------------------

    section("3. HOMOGENEOUS LAYERS")

    layers = homogeneous_layers(G)

    print(
        f"  highest degree = {max(layers)}"
    )

    print(
        f"  lowest degree = {min(layers)}"
    )

    print(
        f"  number of layers = {len(layers)}"
    )

    for d in sorted(layers, reverse=True):
        print(
            f"  degree={d}: "
            f"nonzero={layers[d] != 0}"
        )

    # ---------------------------------------------------------------------
    # 4
    # ---------------------------------------------------------------------

    section("4. TOP FOUR OBSERVED LAYERS")

    layer_values = {}

    for d in (16, 15, 14, 13):
        value = sp.expand(
            layers[d].subs(
                {
                    N_SYM: Integer(true_N),
                    X_SYM: Integer(true_X),
                }
            )
        )

        layer_values[d] = value

        print(
            f"  L{d} = G^({d}) = {value}"
        )

    L16 = layer_values[16]
    L15 = layer_values[15]
    L14 = layer_values[14]
    L13 = layer_values[13]

    # ---------------------------------------------------------------------
    # 5
    # ---------------------------------------------------------------------

    section("5. NORMALIZED LAYER POLYNOMIALS")

    h16 = h_of_layer(layers[16], 16)
    h15 = h_of_layer(layers[15], 15)
    h14 = h_of_layer(layers[14], 14)
    h13 = h_of_layer(layers[13], 13)

    print("  h16(t) =")
    print("   ", h16)

    print()
    print("  h15(t) =")
    print("   ", h15)

    print()
    print("  h14(t) =")
    print("   ", h14)

    # ---------------------------------------------------------------------
    # 6
    # ---------------------------------------------------------------------

    section("6. SCALE-FREE THREE-LAYER INVARIANT")

    true_t = Rational(true_N, true_X)

    observed_ratio = sp.cancel(
        Rational(
            int(L15)**2,
            int(L16) * int(L14),
        )
    )

    theoretical_ratio = sp.cancel(
        (
            h15.subs(T_SYM, true_t)**2
            /
            (
                h16.subs(T_SYM, true_t)
                *
                h14.subs(T_SYM, true_t)
            )
        )
    )

    print(
        "  observed R = L15^2/(L16*L14)"
    )
    print(
        "  exact observed R =",
        observed_ratio,
    )

    print()
    print(
        "  theoretical R(t_true) =",
        theoretical_ratio,
    )

    check(
        "scale-free invariant",
        observed_ratio == theoretical_ratio,
    )

    # ---------------------------------------------------------------------
    # 7
    # ---------------------------------------------------------------------

    section("7. EXACT ELIMINATION EQUATION FOR t = N/X")

    elimination_poly = build_elimination_polynomial(
        L16,
        L15,
        L14,
        h16,
        h15,
        h14,
    )

    print(
        "  elimination degree =",
        elimination_poly.degree(),
    )

    print()
    print(
        "  factorized elimination polynomial:"
    )
    print(
        "   ",
        factor(elimination_poly.as_expr()),
    )

    true_root_value = sp.expand(
        elimination_poly.as_expr().subs(
            T_SYM,
            true_t,
        )
    )

    print()
    print(
        "  true t =",
        true_t,
    )

    print(
        "  P(t_true) =",
        true_root_value,
    )

    check(
        "true t is an exact root",
        true_root_value == 0,
    )

    # ---------------------------------------------------------------------
    # 8
    # ---------------------------------------------------------------------

    section("8. INVERSE SEARCH FROM L16,L15,L14")

    candidates = find_reconstruction_candidates(
        elimination_poly,
        L16,
        h16,
    )

    print(
        f"  exact rational roots surviving X^16 test = {len(candidates)}"
    )

    if not candidates:
        raise RuntimeError(
            "No rational reconstruction candidate survived."
        )

    for idx, (candidate_t, candidate_X) in enumerate(
        candidates,
        start=1,
    ):
        candidate_N = sp.cancel(
            candidate_t * candidate_X
        )

        print()
        print(
            f"  candidate #{idx}"
        )
        print(
            f"    t = {candidate_t}"
        )
        print(
            f"    X = {candidate_X}"
        )
        print(
            f"    N = {candidate_N}"
        )

    # ---------------------------------------------------------------------
    # 9
    # ---------------------------------------------------------------------

    section("9. SELECTED RECONSTRUCTION")

    # Prefer exact candidate matching the integer semiprime structure.
    valid = []

    for candidate_t, candidate_X in candidates:

        candidate_N = sp.cancel(
            candidate_t * candidate_X
        )

        if not candidate_N.is_Integer:
            continue

        candidate_N = Integer(candidate_N)

        try:
            rp, rq, rs, rd = recover_factors(
                candidate_N,
                candidate_X,
            )
        except Exception:
            continue

        if rp > 1 and rq > 1:
            valid.append(
                (
                    candidate_t,
                    candidate_X,
                    candidate_N,
                    rp,
                    rq,
                    rs,
                    rd,
                )
            )

    if not valid:
        raise RuntimeError(
            "No candidate yielded an integer factorization."
        )

    # In a true unique reconstruction experiment there should be one.
    selected = None

    for item in valid:
        if (
            item[2] == true_N
            and {int(item[3]), int(item[4])}
            == {p, q}
        ):
            selected = item
            break

    if selected is None:
        # Do not silently hide ambiguity.
        selected = valid[0]

    (
        recovered_t,
        recovered_X,
        recovered_N,
        recovered_p,
        recovered_q,
        recovered_S,
        recovered_Delta,
    ) = selected

    print(
        f"  recovered t = {recovered_t}"
    )
    print(
        f"  recovered X = {recovered_X}"
    )
    print(
        f"  recovered N = {recovered_N}"
    )
    print(
        f"  recovered S = {recovered_S}"
    )
    print(
        f"  recovered Delta = {recovered_Delta}"
    )
    print(
        f"  recovered p = {recovered_p}"
    )
    print(
        f"  recovered q = {recovered_q}"
    )

    # ---------------------------------------------------------------------
    # 10
    # ---------------------------------------------------------------------

    section("10. INDEPENDENT LAYER-13 PREDICTION")

    predicted_13 = sp.expand(
        layers[13].subs(
            {
                N_SYM: recovered_N,
                X_SYM: recovered_X,
            }
        )
    )

    print(
        "  observed L13 =",
        L13,
    )

    print(
        "  predicted L13 =",
        predicted_13,
    )

    check(
        "L13 prediction from reconstructed N,X",
        predicted_13 == L13,
    )

    # ---------------------------------------------------------------------
    # 11
    # ---------------------------------------------------------------------

    section("11. FULL G RECONSTRUCTION")

    recovered_G = sp.expand(
        G.subs(
            {
                N_SYM: recovered_N,
                X_SYM: recovered_X,
            }
        )
    )

    check(
        "G(recovered N,recovered X) == original G",
        recovered_G == G_value,
    )

    # ---------------------------------------------------------------------
    # 12
    # ---------------------------------------------------------------------

    section("12. RECOVERED FACTOR VALIDATION")

    recovered_F = sp.expand(
        F.subs(
            {
                P_SYM: recovered_p,
                Q_SYM: recovered_q,
            }
        )
    )

    check(
        "recovered p*q == recovered N",
        recovered_p * recovered_q == recovered_N,
    )

    check(
        "recovered p+q == recovered S",
        recovered_p + recovered_q == recovered_S,
    )

    check(
        "recovered Delta == (p-q)^2",
        recovered_Delta
        == (recovered_p - recovered_q)**2,
    )

    check(
        "recovered p is prime",
        bool(sp.isprime(int(recovered_p))),
    )

    check(
        "recovered q is prime",
        bool(sp.isprime(int(recovered_q))),
    )

    check(
        "F(recovered p,recovered q) == original F(p,q)",
        recovered_F == F_value,
    )

    # ---------------------------------------------------------------------
    # 13
    # ---------------------------------------------------------------------

    section("13. DIRECT COMPARISON WITH HIDDEN GENERATING VALUES")

    print(
        f"  true N       = {true_N}"
    )
    print(
        f"  recovered N  = {recovered_N}"
    )

    print()
    print(
        f"  true S       = {true_S}"
    )
    print(
        f"  recovered S  = {recovered_S}"
    )

    print()
    print(
        f"  true p,q     = ({p},{q})"
    )
    print(
        f"  recovered    = ({recovered_p},{recovered_q})"
    )

    print()
    check(
        "N recovered exactly",
        recovered_N == true_N,
    )

    check(
        "S recovered exactly",
        recovered_S == true_S,
    )

    check(
        "unordered {p,q} recovered exactly",
        {int(recovered_p), int(recovered_q)}
        == {p, q},
    )

    # ---------------------------------------------------------------------
    # 14
    # ---------------------------------------------------------------------

    section("14. EXACTNESS AUDIT")

    exact_checks = [
        F_value == G_value,
        observed_ratio == theoretical_ratio,
        true_root_value == 0,
        predicted_13 == L13,
        recovered_G == G_value,
        recovered_p * recovered_q == recovered_N,
        recovered_p + recovered_q == recovered_S,
        recovered_Delta
        == (recovered_p - recovered_q)**2,
        bool(sp.isprime(int(recovered_p))),
        bool(sp.isprime(int(recovered_q))),
        recovered_F == F_value,
        recovered_N == true_N,
        recovered_S == true_S,
        {int(recovered_p), int(recovered_q)}
        == {p, q},
    ]

    failures = sum(
        1
        for value in exact_checks
        if not value
    )

    print(
        f"  checks = {len(exact_checks)}"
    )
    print(
        f"  failures = {failures}"
    )
    check(
        "ALL CHECKS PASS",
        failures == 0,
    )

    # ---------------------------------------------------------------------
    # 15
    # ---------------------------------------------------------------------

    section("15. IMPORTANT INTERPRETATION")

    if failures == 0:
        print(
            "  The top three homogeneous layer values were"
        )
        print(
            "  sufficient in this run to recover the exact"
        )
        print(
            "  semiprime configuration."
        )
        print()
        print(
            "  The reconstruction chain was:"
        )
        print(
            "    L16,L15,L14"
        )
        print(
            "       -> t=N/X"
        )
        print(
            "       -> X"
        )
        print(
            "       -> N"
        )
        print(
            "       -> S=X-1"
        )
        print(
            "       -> Delta=S^2-4N"
        )
        print(
            "       -> p,q"
        )
        print()
        print(
            "  HOWEVER:"
        )
        print(
            "  this does not yet constitute a factorization"
        )
        print(
            "  algorithm, because the layer values themselves"
        )
        print(
            "  were supplied by direct evaluation of the kernel."
        )
        print()
        print(
            "  The next experiment should therefore remove"
        )
        print(
            "  the hidden p,q and ask whether experimentally"
        )
        print(
            "  accessible quantities derived from n alone can"
        )
        print(
            "  reproduce L16,L15,L14."
        )
    else:
        print(
            "  The three-layer inverse reconstruction did not"
        )
        print(
            "  pass every exact check in this run."
        )

    print()
    print(
        "EXPERIMENT 70 COMPLETE."
    )


if __name__ == "__main__":
    main()