from __future__ import annotations

from math import isqrt, gcd


# ==============================================================================
# EXPERIMENT 407
# ==============================================================================
#
# EXACT TRUE/CONJUGATE ROOT CROSS-PRODUCT AUDIT
#
# Established:
#
#   Q_-(Z) = 2Z^2-(N+1-d)Z+2N
#   roots(Q_-) = {p,q}
#
#   Q_+(Z) = 2Z^2-(N+1+d)Z+2N
#   roots(Q_+) = {z1,z2}
#
# with:
#
#   p*q   = N
#   z1*z2 = N
#
#   p+q   = (N+1-d)/2
#   z1+z2 = (N+1+d)/2
#
# Therefore the two branches have identical product but different sums.
#
# Objective:
#
#   1. Compute exact symmetric cross-products:
#
#        (p-z1)(p-z2)
#        (q-z1)(q-z2)
#
#   2. Search for simple closed forms involving N,d,p,q,S.
#
#   3. Compute the four cross-products:
#
#        p*z1, p*z2, q*z1, q*z2
#
#      only through symmetric information; do not numerically approximate
#      irrational conjugate roots.
#
#   4. Construct the cross-product quartic:
#
#        R(W) =
#          (W-p*z1)(W-p*z2)
#          (W-q*z1)(W-q*z2)
#
#      and determine whether its coefficients collapse.
#
#   5. Eliminate p,q in favor of their symmetric data:
#
#        p+q = S
#        pq = N
#
#      and eliminate z1,z2 using:
#
#        z1+z2 = T
#        z1*z2 = N
#
#      where:
#
#        S = (N+1-d)/2
#        T = (N+1+d)/2
#
#   6. Test whether R(W) depends only on N and d.
#
#   7. Substitute:
#
#        d^2 = -4K-3N^2+6N+1
#
#      wherever possible and search for an N,K-only form.
#
#   8. Check whether R(W) has a particularly simple center,
#      e.g. after shifting W by N or N^2.
#
#   9. Test whether the cross-product polynomial is a perfect
#      quadratic in W^2, or has a discriminant/simple factorization.
#
# No resultants.
# No symbolic multivariate factorization.
# Only exact integer arithmetic.
# ==============================================================================


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


def sq_info(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def digits(n: int) -> int:
    return len(str(abs(n)))


def bitlen(n: int) -> int:
    return abs(n).bit_length()


def recover(p: int, q: int):
    N = p * q
    S = p + q
    X = S + 1
    K = -N*N + N*X - X*X + 3*X - 2

    D = -4*K - 3*N*N + 6*N + 1
    ok, d = sq_info(D)

    return N, S, X, K, D, d, ok


# ------------------------------------------------------------------------------
# Polynomial utilities
# Coefficients are stored low -> high.
# ------------------------------------------------------------------------------

def poly_add(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_scale(a, c):
    return [x*c for x in a]


def poly_eval(p, x):
    v = 0
    for c in reversed(p):
        v = v*x + c
    return v


def fmt_poly(p, var="W"):
    terms = []
    for i in range(len(p)-1, -1, -1):
        c = p[i]
        if c == 0:
            continue

        if i == 0:
            term = str(abs(c))
        elif i == 1:
            term = "" if abs(c) == 1 else str(abs(c)) + var
        else:
            term = "" if abs(c) == 1 else str(abs(c)) + var + f"^{i}"

        if not terms:
            terms.append(("-" if c < 0 else "") + term)
        else:
            terms.append((" - " if c < 0 else " + ") + term)

    return "".join(terms) if terms else "0"


def main():

    print("=" * 110)
    print("EXPERIMENT 407 START")
    print("=" * 110)
    print()
    print("EXACT TRUE/CONJUGATE ROOT CROSS-PRODUCT AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        N, S, X, K, D, d, d_ok = recover(p, q)

        local_ok = True

        # True/conjugate root sums.
        #
        # S = p+q = (N+1-d)/2
        # T = z1+z2 = (N+1+d)/2
        #
        # Work with doubled sums to keep everything integral.
        twoS = N + 1 - d
        twoT = N + 1 + d

        if twoS % 2 != 0 or twoT % 2 != 0:
            local_ok = False

        T = twoT // 2

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  K = {K}")
        print()

        print("FIRST RADICAL")
        print(f"  D = {D}")
        print(f"  D = d^2 : {d_ok}")
        print(f"  d = {d}")
        print()

        if not d_ok:
            local_ok = False

        print("BRANCH SUMS")
        print(f"  true branch sum S = {S}")
        print(f"  conjugate branch sum T = {T}")
        print(f"  S+T = {S+T}")
        print(f"  expected N+1 = {N+1}")
        print(f"  S+T identity = {S+T == N+1}")
        print(f"  T-S = {T-S}")
        print(f"  expected d = {d}")
        print(f"  T-S identity = {T-S == d}")
        print()

        if S + T != N + 1 or T - S != d:
            local_ok = False

        # ----------------------------------------------------------------------
        # Symmetric cross-products
        #
        # For any root x of Q_-:
        #
        #   (x-z1)(x-z2)
        #      = x^2 - T*x + N
        #
        # So for x=p,q we get exact integers.
        # ----------------------------------------------------------------------

        cp_p = p*p - T*p + N
        cp_q = q*q - T*q + N

        print("CROSS-DISTANCE PRODUCTS")
        print("  (p-z1)(p-z2) = p^2 - T*p + N")
        print(f"  value = {cp_p}")
        print()
        print("  (q-z1)(q-z2) = q^2 - T*q + N")
        print(f"  value = {cp_q}")
        print()

        # Search for simple factorizations.
        candidates_p = {
            "p*(p-q)": p*(p-q),
            "p*(p-S)": p*(p-S),
            "p*(p-T)": p*(p-T),
            "-d*p": -d*p,
            "-d*q": -d*q,
            "-d*(p)": -d*p,
            "-d*(q)": -d*q,
            "-2*d*p": -2*d*p,
            "-2*d*q": -2*d*q,
        }

        candidates_q = {
            "q*(q-p)": q*(q-p),
            "q*(q-S)": q*(q-S),
            "q*(q-T)": q*(q-T),
            "-d*p": -d*p,
            "-d*q": -d*q,
            "-2*d*p": -2*d*p,
            "-2*d*q": -2*d*q,
        }

        print("SIMPLE CROSS-PRODUCT MATCHES")
        matches = []

        for name, value in candidates_p.items():
            if cp_p == value:
                matches.append(f"(p-z1)(p-z2) = {name}")

        for name, value in candidates_q.items():
            if cp_q == value:
                matches.append(f"(q-z1)(q-z2) = {name}")

        if matches:
            for m in matches:
                print(f"  {m}")
        else:
            print("  no tested simple factorization")
        print()

        # ----------------------------------------------------------------------
        # Cross-product quartic
        #
        # For x=p:
        #
        #   (W-p z1)(W-p z2)
        #
        #   = W^2 - p*T*W + p^2*N
        #
        # Similarly for q.
        #
        # Product gives a quartic in W.
        # ----------------------------------------------------------------------

        Rp = [p*p*N, -p*T, 1]
        Rq = [q*q*N, -q*T, 1]

        R = poly_mul(Rp, Rq)

        print("CROSS-PRODUCT POLYNOMIAL")
        print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
        print(f"  degree = {len(R)-1}")
        print(f"  R(W) = {fmt_poly(R)}")
        print()

        if len(R) != 5 or R[-1] != 1:
            local_ok = False

        # ----------------------------------------------------------------------
        # Derive coefficients symbolically from S,T,N.
        #
        # R(W) =
        #
        # W^4
        # - T(p+q)W^3
        # + N(p^2+q^2)W^2
        # + p q N T (p+q) W
        # + p^2 q^2 N^2
        #
        # Since pq=N and p+q=S:
        #
        # p^2+q^2 = S^2 - 2N
        #
        # W coefficient:
        # p q N T S = N^2 S T
        #
        # constant = N^4.
        #
        # ----------------------------------------------------------------------

        coeff_W4 = 1
        coeff_W3 = -T*S
        coeff_W2 = N*(S*S - 2*N)
        coeff_W1 = N*N*S*T
        coeff_W0 = N**4

        expected_R = [
            coeff_W0,
            coeff_W1,
            coeff_W2,
            coeff_W3,
            coeff_W4,
        ]

        print("SYMMETRIC COEFFICIENT RECONSTRUCTION")
        print(f"  W^4 = {coeff_W4}")
        print(f"  W^3 = {coeff_W3}")
        print(f"  W^2 = {coeff_W2}")
        print(f"  W^1 = {coeff_W1}")
        print(f"  W^0 = {coeff_W0}")
        print()

        print(f"  coefficient reconstruction = {R == expected_R}")
        print()

        if R != expected_R:
            local_ok = False

        # ----------------------------------------------------------------------
        # The key simplification:
        #
        # S*T
        #
        # = ((N+1)^2-d^2)/4
        #
        # and because d^2 is known from N,K, this becomes N,K-only.
        # ----------------------------------------------------------------------

        four_ST = (N+1)*(N+1) - d*d

        print("KEY SYMMETRIC PRODUCT")
        print("  4*S*T = (N+1)^2-d^2")
        print(f"  actual = {4*S*T}")
        print(f"  expected = {four_ST}")
        print(f"  identity = {4*S*T == four_ST}")
        print()

        if 4*S*T != four_ST:
            local_ok = False

        # ----------------------------------------------------------------------
        # N,K-only version of 4ST
        #
        # d^2 = -4K - 3N^2 + 6N + 1
        #
        # => 4ST
        #    = N^2+2N+1 + 4K + 3N^2 - 6N - 1
        #    = 4(N^2-N+K)
        #
        # Therefore:
        #
        #        S*T = N^2-N+K
        #
        # This is a particularly important identity.
        # ----------------------------------------------------------------------

        st_nk = N*N - N + K

        print("N,K-ONLY CROSS-BRANCH PRODUCT")
        print("  S*T = N^2-N+K")
        print(f"  S*T = {S*T}")
        print(f"  N^2-N+K = {st_nk}")
        print(f"  identity = {S*T == st_nk}")
        print()

        if S*T != st_nk:
            local_ok = False

        # ----------------------------------------------------------------------
        # N,K-only cross-polynomial coefficients.
        #
        # W^3 = -S*T
        #      = -(N^2-N+K)
        #
        # W^2 = N(S^2-2N)
        #
        # Need S^2:
        #
        # S = (N+1-d)/2
        #
        # Instead:
        #
        # S^2 = ((N+1)^2 + d^2 - 2d(N+1))/4
        #
        # This still has d linearly, so W^2 may retain a d term.
        #
        # But W^1 depends on ST and therefore is N,K-only.
        # ----------------------------------------------------------------------

        print("N,K-ONLY COEFFICIENT AUDIT")

        nk_W3 = -(N*N - N + K)
        nk_W1 = N*N * st_nk
        nk_W0 = N**4

        print(f"  W^3 coefficient = -(N^2-N+K) = {nk_W3}")
        print(f"  actual W^3 = {R[3]}")
        print(f"  W^3 identity = {R[3] == nk_W3}")
        print()

        print(f"  W^1 coefficient = N^2*(N^2-N+K) = {nk_W1}")
        print(f"  actual W^1 = {R[1]}")
        print(f"  W^1 identity = {R[1] == nk_W1}")
        print()

        print(f"  W^0 coefficient = N^4 = {nk_W0}")
        print(f"  actual W^0 = {R[0]}")
        print(f"  W^0 identity = {R[0] == nk_W0}")
        print()

        if R[3] != nk_W3 or R[1] != nk_W1 or R[0] != nk_W0:
            local_ok = False

        # ----------------------------------------------------------------------
        # Center transformation W -> N/Y or W -> N*Y
        #
        # Since each branch has product N, reciprocal roots may simplify.
        #
        # If x1*x2=N then:
        #
        #   1/x1 + 1/x2 = (x1+x2)/N.
        #
        # Therefore test reciprocal branch sums:
        #
        #   1/p + 1/q = S/N
        #   1/z1 + 1/z2 = T/N
        #
        # and their difference:
        #
        #   1/z1 + 1/z2 - (1/p+1/q) = d/N.
        # ----------------------------------------------------------------------

        print("RECIPROCAL BRANCH RELATION")
        print("  (1/z1 + 1/z2) - (1/p + 1/q) = d/N")
        reciprocal_num = T - S
        print(f"  numerator = {reciprocal_num}")
        print(f"  expected d = {d}")
        print(f"  identity = {reciprocal_num == d}")
        print()

        if reciprocal_num != d:
            local_ok = False

        # ----------------------------------------------------------------------
        # Product of cross-products:
        #
        # (p-z1)(p-z2)(q-z1)(q-z2)
        #
        # This is R(0) = N^4.
        #
        # Check directly in symmetric form.
        # ----------------------------------------------------------------------

        cross_product = cp_p * cp_q

        print("CROSS-DISTANCE PRODUCT")
        print("  (p-z1)(p-z2)(q-z1)(q-z2)")
        print(f"  actual = {cross_product}")
        print(f"  expected N^4 = {N**4}")
        print(f"  identity = {cross_product == N**4}")
        print()

        # Note:
        # cp_p and cp_q are generally negative, so their product is positive.

        if cross_product != N**4:
            local_ok = False

        # ----------------------------------------------------------------------
        # Additional low-complexity candidates for cp_p + cp_q and cp_p*cp_q.
        # ----------------------------------------------------------------------

        print("CROSS-DISTANCE SYMMETRIC INVARIANTS")

        cp_sum = cp_p + cp_q

        # Derived directly:
        # cp_p + cp_q
        #   = p^2+q^2 - T(p+q) + 2N
        #   = S^2-2N - TS + 2N
        #   = S(S-T)
        #   = -S*d
        #
        expected_cp_sum = -S*d

        print(f"  C = cp_p + cp_q = {cp_sum}")
        print(f"  expected -S*d = {expected_cp_sum}")
        print(f"  identity = {cp_sum == expected_cp_sum}")
        print()

        if cp_sum != expected_cp_sum:
            local_ok = False

        print("KEY CROSS-BRANCH IDENTITY")
        print("  (p-z1)(p-z2) + (q-z1)(q-z2) = -S*d")
        print()

        # ----------------------------------------------------------------------
        # Factor the quadratic whose roots are cp_p and cp_q.
        #
        # Their:
        #
        #   sum  = -S*d
        #   product = N^4
        #
        # Therefore:
        #
        #   C^2 + S*d*C + N^4 = 0
        #
        # This may itself expose a new layer.
        # ----------------------------------------------------------------------

        print("CROSS-DISTANCE QUADRATIC")
        print("  C^2 + S*d*C + N^4 = 0")
        print()

        check_cp_p = cp_p*cp_p + S*d*cp_p + N**4
        check_cp_q = cp_q*cp_q + S*d*cp_q + N**4

        print(f"  polynomial(cp_p) = {check_cp_p}")
        print(f"  polynomial(cp_q) = {check_cp_q}")
        print(
            f"  both roots = "
            f"{check_cp_p == 0 and check_cp_q == 0}"
        )
        print()

        if check_cp_p != 0 or check_cp_q != 0:
            local_ok = False

        # ----------------------------------------------------------------------
        # N,K-only version:
        #
        # S*d = ((N+1-d)d)/2
        #
        # This still contains d.
        #
        # Square the coefficient:
        #
        # (S*d)^2 = S^2*d^2
        #
        # which is N,K-only.
        # Audit this explicitly.
        # ----------------------------------------------------------------------

        Sd2 = (S*d)**2
        Sd2_nk = S*S * D

        print("S*d MAGNITUDE SQUARE")
        print("  (S*d)^2 = S^2*d^2")
        print(f"  actual = {Sd2}")
        print(f"  reconstructed = {Sd2_nk}")
        print(f"  identity = {Sd2 == Sd2_nk}")
        print()

        if Sd2 != Sd2_nk:
            local_ok = False

        # ----------------------------------------------------------------------
        # Size
        # ----------------------------------------------------------------------

        print("SIZE DATA")
        print(f"  N digits            = {digits(N)}")
        print(f"  K digits            = {digits(K)}")
        print(f"  d digits            = {digits(d)}")
        print(f"  S digits            = {digits(S)}")
        print(f"  T digits            = {digits(T)}")
        print(f"  cp_p digits         = {digits(cp_p)}")
        print(f"  cp_q digits         = {digits(cp_q)}")
        print(f"  N^4 digits          = {digits(N**4)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

        global_ok = global_ok and local_ok

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("BRANCH DATA")
    print("  p*q = z1*z2 = N")
    print("  (z1+z2)-(p+q) = d")
    print()

    print("NEW CROSS-BRANCH IDENTITY")
    print("  S*T = N^2-N+K")
    print()
    print("  where")
    print("    S = p+q")
    print("    T = z1+z2")
    print()

    print("CROSS-DISTANCE IDENTITY")
    print("  C_p + C_q = -S*d")
    print("  C_p*C_q = N^4")
    print("  therefore")
    print("    C^2 + S*d*C + N^4 = 0")
    print("  for C in {C_p,C_q}.")
    print()

    print("CROSS-PRODUCT POLYNOMIAL")
    print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
    print("  Its coefficients collapse substantially under")
    print("  p+q=S, pq=N, z1+z2=T, z1z2=N.")
    print()

    print("CORE QUESTION")
    print("  Does the cross-branch quadratic expose another")
    print("  low-complexity invariant linking N,K directly?")
    print()

    print("=" * 110)
    print("EXPERIMENT 407 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 407 FINISHED")
    print()


if __name__ == "__main__":
    main()
