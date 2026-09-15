from __future__ import annotations

from math import isqrt


# ==============================================================================
# EXPERIMENT 408
# ==============================================================================
#
# EXACT N,K-ONLY CROSS-BRANCH POLYNOMIAL / DISTANCE QUADRATIC AUDIT
#
# Corrections from EXPERIMENT 407:
#
#   1. For
#
#        R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)
#
#      the W coefficient is NEGATIVE:
#
#        -N^2*S*T
#
#   2. The W^2 coefficient is:
#
#        N(S^2 + T^2 - 2N)
#
#      not merely N(S^2-2N).
#
#   3. For
#
#        C_p = (p-z1)(p-z2)
#        C_q = (q-z1)(q-z2)
#
#      we have:
#
#        C_p = -d*p
#        C_q = -d*q
#
#      hence:
#
#        C_p + C_q = -d*S
#        C_p*C_q   = d^2*N
#
#      and therefore:
#
#        C^2 + d*S*C + d^2*N = 0.
#
#   4. Since
#
#        d^2 = -4K - 3N^2 + 6N + 1,
#
#      the constant term d^2*N is immediately N,K-only.
#
#   5. The full cross-product quartic is expected to reduce to:
#
#      R(W) =
#        W^4
#        -(N^2-N+K)W^3
#        +N(-N^2+2N+1-2K)W^2
#        -N^2(N^2-N+K)W
#        +N^4
#
# Objective:
#
#   - Verify all corrected identities exactly.
#   - Verify the cross-distance quadratic.
#   - Verify the complete quartic as an N,K-only polynomial.
#   - Search for reciprocal / palindromic structure.
#   - Test the transformation W -> N/W.
#   - Test whether R(W) is self-reciprocal after normalization.
#   - Test simple factorizations involving W-N, W+N, W^2-NW+N^2.
#   - Determine whether the cross-branch layer adds genuinely new
#     information or is completely determined by N,K.
#
# No resultants.
# No symbolic multivariate factoring.
# Exact integer arithmetic only.
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


def exact_sqrt(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def digits(n: int) -> int:
    return len(str(abs(n)))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_eval(poly, x):
    v = 0
    for c in reversed(poly):
        v = v * x + c
    return v


def fmt_poly(poly, var="W"):
    terms = []

    for power in range(len(poly) - 1, -1, -1):
        c = poly[power]
        if c == 0:
            continue

        ac = abs(c)

        if power == 0:
            body = str(ac)
        elif power == 1:
            body = var if ac == 1 else f"{ac}{var}"
        else:
            body = var if ac == 1 else f"{ac}{var}"
            body += f"^{power}"

        if not terms:
            terms.append(("-" if c < 0 else "") + body)
        else:
            terms.append((" - " if c < 0 else " + ") + body)

    return "".join(terms) if terms else "0"


def main():

    print("=" * 110)
    print("EXPERIMENT 408 START")
    print("=" * 110)
    print()
    print("EXACT N,K-ONLY CROSS-BRANCH POLYNOMIAL / DISTANCE QUADRATIC AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        N = p * q
        S = p + q
        X = S + 1

        # Same K construction used throughout the experiments.
        K = -N * N + N * X - X * X + 3 * X - 2

        D = -4 * K - 3 * N * N + 6 * N + 1
        d_ok, d = exact_sqrt(D)

        local_ok = d_ok

        T = (N + 1 + d) // 2

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  T = {T}")
        print(f"  K = {K}")
        print()

        print("FIRST RADICAL")
        print(f"  D = {D}")
        print(f"  D = d^2 : {d_ok}")
        print(f"  d = {d}")
        print()

        # ------------------------------------------------------------------
        # True/conjugate branch checks.
        # ------------------------------------------------------------------

        print("BRANCH SUMS")

        branch_sum_check = (S + T == N + 1)
        branch_diff_check = (T - S == d)

        print(f"  S+T = {S+T}")
        print(f"  expected N+1 = {N+1}")
        print(f"  S+T identity = {branch_sum_check}")
        print()

        print(f"  T-S = {T-S}")
        print(f"  expected d = {d}")
        print(f"  T-S identity = {branch_diff_check}")
        print()

        if not branch_sum_check or not branch_diff_check:
            local_ok = False

        # ------------------------------------------------------------------
        # Correct cross-distance quantities.
        # ------------------------------------------------------------------

        cp = p * p - T * p + N
        cq = q * q - T * q + N

        print("CROSS-DISTANCE VALUES")
        print(f"  C_p = (p-z1)(p-z2) = {cp}")
        print(f"  C_q = (q-z1)(q-z2) = {cq}")
        print()

        cp_simple = -d * p
        cq_simple = -d * q

        print("SIMPLE CROSS-DISTANCE IDENTITIES")
        print(f"  C_p = -d*p : {cp == cp_simple}")
        print(f"  C_q = -d*q : {cq == cq_simple}")
        print()

        if cp != cp_simple or cq != cq_simple:
            local_ok = False

        # ------------------------------------------------------------------
        # Cross-distance sum/product.
        # ------------------------------------------------------------------

        Csum = cp + cq
        Cprod = cp * cq

        expected_Csum = -d * S
        expected_Cprod = d * d * N

        print("CROSS-DISTANCE SUM")
        print(f"  C_p+C_q = {Csum}")
        print(f"  -d*S    = {expected_Csum}")
        print(f"  identity = {Csum == expected_Csum}")
        print()

        print("CROSS-DISTANCE PRODUCT")
        print(f"  C_p*C_q = {Cprod}")
        print(f"  d^2*N   = {expected_Cprod}")
        print(f"  identity = {Cprod == expected_Cprod}")
        print()

        if Csum != expected_Csum or Cprod != expected_Cprod:
            local_ok = False

        # ------------------------------------------------------------------
        # Correct cross-distance quadratic.
        # ------------------------------------------------------------------

        q_cp_p = cp * cp + d * S * cp + d * d * N
        q_cp_q = cq * cq + d * S * cq + d * d * N

        print("CORRECTED CROSS-DISTANCE QUADRATIC")
        print("  C^2 + d*S*C + d^2*N = 0")
        print(f"  polynomial(C_p) = {q_cp_p}")
        print(f"  polynomial(C_q) = {q_cp_q}")
        print(
            f"  both roots = "
            f"{q_cp_p == 0 and q_cp_q == 0}"
        )
        print()

        if q_cp_p != 0 or q_cp_q != 0:
            local_ok = False

        # ------------------------------------------------------------------
        # N,K-only coefficient d^2*N.
        # ------------------------------------------------------------------

        d2_nk = -4 * K - 3 * N * N + 6 * N + 1
        nk_Cprod = N * d2_nk

        print("N,K-ONLY CROSS-DISTANCE PRODUCT")
        print("  C_p*C_q = N*(-4K-3N^2+6N+1)")
        print(f"  actual = {Cprod}")
        print(f"  N,K expression = {nk_Cprod}")
        print(f"  identity = {Cprod == nk_Cprod}")
        print()

        if Cprod != nk_Cprod:
            local_ok = False

        # ------------------------------------------------------------------
        # ST = N^2-N+K.
        # ------------------------------------------------------------------

        ST = S * T
        ST_nk = N * N - N + K

        print("BRANCH-SUM PRODUCT")
        print("  S*T = N^2-N+K")
        print(f"  S*T = {ST}")
        print(f"  N^2-N+K = {ST_nk}")
        print(f"  identity = {ST == ST_nk}")
        print()

        if ST != ST_nk:
            local_ok = False

        # ------------------------------------------------------------------
        # Correct R(W).
        #
        # For p:
        #
        #   W^2 - p*T*W + p^2*N
        #
        # For q:
        #
        #   W^2 - q*T*W + q^2*N
        #
        # Product:
        #
        # W^4
        # -ST W^3
        # +N(S^2+T^2-2N)W^2
        # -N^2ST W
        # +N^4
        # ------------------------------------------------------------------

        R1 = [p * p * N, -p * T, 1]
        R2 = [q * q * N, -q * T, 1]
        R = poly_mul(R1, R2)

        expected_W4 = 1
        expected_W3 = -S * T
        expected_W2 = N * (S * S + T * T - 2 * N)
        expected_W1 = -N * N * S * T
        expected_W0 = N ** 4

        expected_R = [
            expected_W0,
            expected_W1,
            expected_W2,
            expected_W3,
            expected_W4,
        ]

        print("CORRECTED CROSS-PRODUCT POLYNOMIAL")
        print("  R(W) = (W-p*z1)(W-p*z2)(W-q*z1)(W-q*z2)")
        print(f"  degree = {len(R)-1}")
        print(f"  R(W) = {fmt_poly(R)}")
        print()

        print("COEFFICIENT RECONSTRUCTION")
        print(f"  W^4 = {R[4]} ; expected {expected_W4}")
        print(f"  W^3 = {R[3]} ; expected {expected_W3}")
        print(f"  W^2 = {R[2]} ; expected {expected_W2}")
        print(f"  W^1 = {R[1]} ; expected {expected_W1}")
        print(f"  W^0 = {R[0]} ; expected {expected_W0}")
        print(f"  exact coefficient reconstruction = {R == expected_R}")
        print()

        if R != expected_R:
            local_ok = False

        # ------------------------------------------------------------------
        # Fully N,K-only form.
        #
        # ST = N^2-N+K
        #
        # S^2+T^2:
        #
        #   (S+T)^2 - 2ST
        #   = (N+1)^2 - 2(N^2-N+K)
        #   = -N^2 + 4N + 1 - 2K
        #
        # Thus:
        #
        #   S^2+T^2-2N
        #   = -N^2 + 2N + 1 - 2K.
        # ------------------------------------------------------------------

        nk_W3 = -(N * N - N + K)
        nk_W2 = N * (-N * N + 2 * N + 1 - 2 * K)
        nk_W1 = -N * N * (N * N - N + K)
        nk_W0 = N ** 4

        nk_R = [
            nk_W0,
            nk_W1,
            nk_W2,
            nk_W3,
            1,
        ]

        print("FULL N,K-ONLY QUARTIC")
        print("  R_NK(W) =")
        print("    W^4")
        print("    -(N^2-N+K)W^3")
        print("    +N(-N^2+2N+1-2K)W^2")
        print("    -N^2(N^2-N+K)W")
        print("    +N^4")
        print()

        print(f"  W^3 N,K identity = {R[3] == nk_W3}")
        print(f"  W^2 N,K identity = {R[2] == nk_W2}")
        print(f"  W^1 N,K identity = {R[1] == nk_W1}")
        print(f"  W^0 N,K identity = {R[0] == nk_W0}")
        print(f"  complete N,K-only identity = {R == nk_R}")
        print()

        if R != nk_R:
            local_ok = False

        # ------------------------------------------------------------------
        # Reciprocal symmetry.
        #
        # Check:
        #
        #   W^4 R(N/W)
        #
        # and compare with a simple transformation of R(W).
        #
        # Also inspect whether:
        #
        #   coefficients satisfy a3/a0 = a1/a4
        #
        # ------------------------------------------------------------------

        reciprocal_scaled = [
            R[4] * N**4,
            R[3] * N**3,
            R[2] * N**2,
            R[1] * N,
            R[0],
        ]

        # Coefficient list corresponds to W^0..W^4 after W -> N/W,
        # multiplied by W^4.

        print("RECIPROCAL TRANSFORMATION")
        print("  Test polynomial W^4*R(N/W)")
        print(
            f"  transformed coefficients (low->high) = "
            f"{reciprocal_scaled}"
        )
        print()

        # Simple comparison:
        # R(N/W) transformation is not expected to equal R(W),
        # but the branch product N may induce a predictable relation.
        #
        # Compare normalized coefficient ratios where integral.
        a0, a1, a2, a3, a4 = R

        ratio_check_1 = a1 * a0 == a3 * a4
        ratio_check_2 = a1 == -N * N * (-a3)
        ratio_check_3 = a0 == N**4 * a4

        print("PALINDROMIC / ANTI-PALINDROMIC AUDIT")
        print(f"  a0 = N^4*a4 : {ratio_check_3}")
        print(f"  a1 = N^2*a3 : {a1 == a3 * N * N}")
        print(
            "  a1 = -N^2*a3 : "
            f"{a1 == -a3 * N * N}"
        )
        print(
            "  a0*a1 = a3*a4 : "
            f"{ratio_check_1}"
        )
        print()

        # The expected relation is:
        #   a1 = N^2*a3
        # because both a1 and a3 are negative.
        if not ratio_check_3 or a1 != a3 * N * N:
            local_ok = False

        # ------------------------------------------------------------------
        # Evaluate special values.
        # ------------------------------------------------------------------

        special = {
            "R(0)": poly_eval(R, 0),
            "R(N)": poly_eval(R, N),
            "R(-N)": poly_eval(R, -N),
            "R(N^2)": poly_eval(R, N*N),
            "R(1)": poly_eval(R, 1),
            "R(-1)": poly_eval(R, -1),
        }

        print("SPECIAL-VALUE AUDIT")
        for name, value in special.items():
            print(f"  {name} = {value}")
        print()

        # ------------------------------------------------------------------
        # Search simple divisibility/factor candidates for R(W).
        # ------------------------------------------------------------------

        print("SIMPLE ROOT / FACTOR CANDIDATES")

        candidates = {
            "W=N": N,
            "W=-N": -N,
            "W=1": 1,
            "W=-1": -1,
            "W=N^2": N*N,
            "W=-N^2": -(N*N),
        }

        any_special_root = False

        for label, w in candidates.items():
            value = poly_eval(R, w)
            hit = (value == 0)
            print(f"  {label:<8} -> root = {hit}")
            any_special_root = any_special_root or hit

        print()

        # ------------------------------------------------------------------
        # Important corrected product:
        #
        # C_p*C_q = d^2*N.
        #
        # Search whether d^2*N has especially compact N,K forms.
        # ------------------------------------------------------------------

        print("CROSS-DISTANCE PRODUCT FORMS")

        compact_d2N = N * (-4*K - 3*N*N + 6*N + 1)

        print(f"  d^2*N = {d*d*N}")
        print(f"  N*(-4K-3N^2+6N+1) = {compact_d2N}")
        print(f"  identity = {d*d*N == compact_d2N}")
        print()

        if d*d*N != compact_d2N:
            local_ok = False

        # ------------------------------------------------------------------
        # Size.
        # ------------------------------------------------------------------

        print("SIZE DATA")
        print(f"  N digits       = {digits(N)}")
        print(f"  K digits       = {digits(K)}")
        print(f"  d digits       = {digits(d)}")
        print(f"  R(W^2) digits  = {digits(R[2])}")
        print(f"  R(W^1) digits  = {digits(R[1])}")
        print(f"  R(W^0) digits  = {digits(R[0])}")
        print(f"  d^2*N digits   = {digits(d*d*N)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

        global_ok = global_ok and local_ok

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("CORRECTED CROSS-DISTANCE STRUCTURE")
    print("  C_p = -d*p")
    print("  C_q = -d*q")
    print("  C_p+C_q = -d*S")
    print("  C_p*C_q = d^2*N")
    print()

    print("CORRECTED CROSS-DISTANCE QUADRATIC")
    print("  C^2 + d*S*C + d^2*N = 0")
    print()

    print("BRANCH-SUM INVARIANT")
    print("  S*T = N^2-N+K")
    print()

    print("FULL N,K-ONLY CROSS-PRODUCT POLYNOMIAL")
    print("  R(W) = W^4")
    print("       -(N^2-N+K)W^3")
    print("       +N(-N^2+2N+1-2K)W^2")
    print("       -N^2(N^2-N+K)W")
    print("       +N^4")
    print()

    print("RECIPROCAL STRUCTURE")
    print("  a0 = N^4*a4")
    print("  a1 = N^2*a3")
    print()

    print("CORE QUESTION")
    print("  Does the complete cross-branch polynomial contain")
    print("  any information beyond the N,K invariants already")
    print("  forced by the two conjugate quadratics?")
    print()

    print("=" * 110)
    print("EXPERIMENT 408 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 408 FINISHED")
    print()


if __name__ == "__main__":
    main()
