from __future__ import annotations

from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 406
# ==============================================================================
#
# EXACT CONJUGATE-ROOT TRANSFORMATION AUDIT
#
# Objective:
#   1. Recover d from N,K.
#   2. Construct the true quadratic:
#
#        Q_-(Z) = 2Z^2-(N+1-d)Z+2N
#
#      with roots p,q.
#
#   3. Construct the conjugate quadratic:
#
#        Q_+(Z) = 2Z^2-(N+1+d)Z+2N
#
#      with formal roots z1,z2.
#
#   4. Determine whether the conjugate roots admit simple
#      algebraic transformations involving:
#
#        p, q, N, S=p+q, g=q-p, d
#
#   5. Test candidate relations including:
#
#        z1*z2 = N
#        z1+z2 = (N+1+d)/2
#
#        z_i * p
#        z_i / p
#        z_i * q
#        z_i / q
#
#        z_i - p
#        z_i - q
#
#        z_i + p
#        z_i + q
#
#        z_i * (p-2)
#        z_i * (q-2)
#
#        z_i / (p-2)
#        z_i / (q-2)
#
#   6. Since Delta_+ is nonsquare in all tested instances,
#      avoid pretending that z1,z2 are rational integers.
#
#   7. Work instead with their symmetric combinations and
#      test whether simple expressions involving sqrt(Delta_+)
#      collapse to known quantities.
#
#   8. Define:
#
#        T = N+1+d
#        H = sqrt(Delta_+)   (symbolic object, not assumed integer)
#
#      so:
#
#        z1 = (T+H)/4
#        z2 = (T-H)/4
#
#   9. Test exact identities for:
#
#        z1*z2
#        (z1-p)(z2-p)
#        (z1-q)(z2-q)
#        (z1-(p-2))(z2-(q-2))
#        (z1-p)(z1-q)
#
#      by eliminating H symbolically using H^2=Delta_+.
#
#  10. Search for especially simple expressions in N,d,p,q.
#
#  11. Determine whether the conjugate branch is merely a
#      Galois/sign-conjugate object or whether it encodes an
#      additional arithmetic transformation.
#
# No resultants.
# No symbolic multivariate factorization.
# Only exact integer arithmetic plus controlled symbolic
# reduction using H^2 = Delta_+.
#
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


def digits(n: int) -> int:
    return len(str(abs(n)))


def bit_length(n: int) -> int:
    return abs(n).bit_length()


def square_info(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def recover_data(p: int, q: int):
    N = p * q
    S = p + q
    X = S + 1

    K = -N**2 + N*X - X**2 + 3*X - 2

    D = -4*K - 3*N**2 + 6*N + 1
    ok, d = square_info(D)

    return {
        "p": p,
        "q": q,
        "N": N,
        "S": S,
        "X": X,
        "K": K,
        "D": D,
        "d": d,
        "d_ok": ok,
    }


def qminus(z: int, N: int, d: int) -> int:
    return 2*z*z - (N+1-d)*z + 2*N


def qplus(z: int, N: int, d: int) -> int:
    return 2*z*z - (N+1+d)*z + 2*N


# ------------------------------------------------------------------------------
# Controlled quadratic-field arithmetic
#
# Represent a + b*H as tuple (a,b), where H^2 = delta.
#
# We only need multiplication and equality after reduction.
# ------------------------------------------------------------------------------


def pair_add(x, y):
    return x[0] + y[0], x[1] + y[1]


def pair_sub(x, y):
    return x[0] - y[0], x[1] - y[1]


def pair_mul(x, y, delta):
    a, b = x
    c, d = y
    return (
        a*c + b*d*delta,
        a*d + b*c,
    )


def pair_scale(x, k: int):
    return x[0] * k, x[1] * k


def pair_zero(x) -> bool:
    return x[0] == 0 and x[1] == 0


def pair_str(x) -> str:
    a, b = x
    if b == 0:
        return str(a)
    if a == 0:
        return f"({b})*H"
    sign = "+" if b > 0 else "-"
    return f"{a} {sign} {abs(b)}*H"


def main():

    print("=" * 110)
    print("EXPERIMENT 406 START")
    print("=" * 110)
    print()
    print("EXACT CONJUGATE-ROOT TRANSFORMATION AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact arithmetic plus H^2-reduction is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        data = recover_data(p, q)

        p = data["p"]
        q = data["q"]
        N = data["N"]
        S = data["S"]
        K = data["K"]
        D = data["D"]
        d = data["d"]

        local_ok = True

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
        print(f"  D = d^2 : {data['d_ok']}")
        print(f"  d = {d}")
        print()

        if not data["d_ok"]:
            local_ok = False

        # ----------------------------------------------------------------------
        # True quadratic
        # ----------------------------------------------------------------------

        qm_p = qminus(p, N, d)
        qm_q = qminus(q, N, d)

        print("TRUE FACTOR QUADRATIC")
        print("  Q_-(Z) = 2Z^2-(N+1-d)Z+2N")
        print(f"  Q_-(p) = {qm_p}")
        print(f"  Q_-(q) = {qm_q}")
        print(f"  exact roots = {qm_p == 0 and qm_q == 0}")
        print()

        if qm_p != 0 or qm_q != 0:
            local_ok = False

        # ----------------------------------------------------------------------
        # Conjugate quadratic
        # ----------------------------------------------------------------------

        delta_plus = (N+1+d)**2 - 16*N
        plus_sq, plus_root = square_info(delta_plus)

        T = N + 1 + d

        print("CONJUGATE QUADRATIC")
        print("  Q_+(Z) = 2Z^2-(N+1+d)Z+2N")
        print(f"  Delta_+ = {delta_plus}")
        print(f"  Delta_+ square = {plus_sq}")
        print()

        if plus_sq:
            print(f"  H = sqrt(Delta_+) = {plus_root}")
        else:
            print("  H = sqrt(Delta_+) is irrational over Q")

        print()
        print("FORMAL CONJUGATE ROOTS")
        print("  z1 = (T+H)/4")
        print("  z2 = (T-H)/4")
        print(f"  T = {T}")
        print()

        # Each root is represented as:
        # z1 = T/4 + H/4
        # z2 = T/4 - H/4
        z1 = (T, 1)   # numerator over 4
        z2 = (T, -1)

        # ----------------------------------------------------------------------
        # Vieta
        # ----------------------------------------------------------------------

        print("CONJUGATE VIETA")
        print("  z1+z2 = T/2")
        print("  z1*z2 = N")

        # Since z1,z2 have denominator 4:
        # numerator product = (T+H)(T-H) = T^2-H^2
        #
        # divide by 16.
        product_num = T*T - delta_plus

        print(f"  T^2 - Delta_+ = {product_num}")
        print(f"  16N = {16*N}")
        print(f"  product identity = {product_num == 16*N}")
        print()

        if product_num != 16*N:
            local_ok = False

        # ----------------------------------------------------------------------
        # Root subtraction products
        #
        # For arbitrary integer c:
        #
        #   (z1-c)(z2-c)
        #     = z1*z2 - c(z1+z2) + c^2
        #
        #     = N - c*T/2 + c^2
        #
        # This removes H entirely.
        # ----------------------------------------------------------------------

        print("SYMMETRIC DISTANCE PRODUCTS")

        test_constants = {
            "p": p,
            "q": q,
            "p-1": p - 1,
            "q-1": q - 1,
            "p-2": p - 2,
            "q-2": q - 2,
            "S": S,
            "d": d,
            "d+1": d + 1,
            "d+3": d + 3,
        }

        for name, c in test_constants.items():

            # Exact rational numerator over 2:
            # product = (2N - c*T + 2*c^2)/2
            num = 2*N - c*T + 2*c*c

            if num % 2 == 0:
                value = num // 2
                print(f"  (z1-{name})(z2-{name}) = {value}")
            else:
                print(
                    f"  (z1-{name})(z2-{name}) = "
                    f"{num}/2"
                )

        print()

        # ----------------------------------------------------------------------
        # Particularly important candidates
        # ----------------------------------------------------------------------

        def distance_product(c):
            return N - c*T/2 + c*c

        # Use denominator-free form.
        print("SPECIAL ROOT RELATIONS")

        # (z1-p)(z2-p)
        num_pp = 2*N - p*T + 2*p*p
        num_qq = 2*N - q*T + 2*q*q

        print(
            "  2*(z1-p)(z2-p) = "
            f"{num_pp}"
        )
        print(
            "  2*(z1-q)(z2-q) = "
            f"{num_qq}"
        )

        # Since p*q=N:
        # polynomial Q+ evaluated at p is:
        # Q+(p) = 2(z1-p)(z2-p)
        #
        # This should exactly match -2dp.
        #
        print()
        print("  Q_+(p) relation:")
        print(f"    Q_+(p) = {qplus(p, N, d)}")
        print(f"    2*(z1-p)(z2-p) = {num_pp}")
        print(
            f"    equality = "
            f"{qplus(p, N, d) == num_pp}"
        )

        print()
        print("  Q_+(q) relation:")
        print(f"    Q_+(q) = {qplus(q, N, d)}")
        print(f"    2*(z1-q)(z2-q) = {num_qq}")
        print(
            f"    equality = "
            f"{qplus(q, N, d) == num_qq}"
        )

        if qplus(p, N, d) != num_pp:
            local_ok = False
        if qplus(q, N, d) != num_qq:
            local_ok = False

        print()

        # ----------------------------------------------------------------------
        # Test whether z1/z2 can equal simple rational quantities.
        #
        # Since Delta_+ is generally nonsquare, exact equality to a rational
        # candidate would require the H coefficient to disappear.
        #
        # We test only candidates where a numerical/simple rational relation
        # would plausibly occur.
        # ----------------------------------------------------------------------

        print("RATIONAL-CANDIDATE ROOT TESTS")

        candidates = {
            "p": p,
            "q": q,
            "p-2": p - 2,
            "q-2": q - 2,
            "S": S,
            "S/2": (S, 2),
            "N/p": (N, p),
            "N/q": (N, q),
            "d": d,
            "d/2": (d, 2),
            "N/d": (N, d),
        }

        for name, value in candidates.items():

            # Rational equality z1 = a/b:
            # (T+H)/4 = a/b
            #
            # This can only hold if:
            #   H = 4a/b - T
            #
            # and therefore H would be rational.
            #
            # So if Delta_+ is nonsquare, all such tests fail automatically.
            #
            if isinstance(value, tuple):
                a, b = value
                # Candidate as reduced rational is not necessary.
                # Compare the required H^2 exactly:
                lhs_num = 4*a - T*b
                rhs = lhs_num * lhs_num
                lhs = delta_plus * b * b
                equal = lhs == rhs
            else:
                a = value
                rhs = (4*a - T) ** 2
                equal = rhs == delta_plus

            print(f"  z1 == {name:<8} : {equal}")

        print()

        # ----------------------------------------------------------------------
        # Test simple transformed quantities involving H.
        #
        # Compute:
        #
        #   4z1 - T = H
        #   T - 4z2 = H
        #
        # Therefore:
        #
        #   4(z1-p) + (p relation)
        #
        # We explicitly check whether H ± simple integer becomes a known
        # multiple of d.
        # ----------------------------------------------------------------------

        print("H-TRANSFORM TESTS")

        known_values = {
            "d": d,
            "2d": 2*d,
            "N": N,
            "S": S,
            "q-p": q-p,
            "p+q": S,
            "d+N": d+N,
            "N-d": N-d,
            "N+1+d": N+1+d,
            "N+1-d": N+1-d,
        }

        # Test Delta_+ against (known + offset)^2.
        offsets = [-4, -2, -1, 0, 1, 2, 4]

        for name, value in known_values.items():
            for off in offsets:
                cand = value + off
                if cand * cand == delta_plus:
                    print(
                        f"  sqrt(Delta_+) = {name} + ({off}) "
                        f"--> EXACT"
                    )

        print("  no exact simple candidate found" if not plus_sq else
              "  Delta_+ itself is square")
        print()

        # ----------------------------------------------------------------------
        # Product/sum transformations
        # ----------------------------------------------------------------------

        print("CONJUGATE-BRANCH SYMMETRIC DATA")

        conj_sum_twice = T
        conj_product = N

        true_sum = S
        true_product = N

        print(f"  true sum              = {true_sum}")
        print(f"  conjugate sum * 2     = {conj_sum_twice}")
        print(f"  true product          = {true_product}")
        print(f"  conjugate product     = {conj_product}")
        print()

        # Difference of the branch sums.
        #
        # ((N+1+d)/2) - S
        #   = (N+1+d - 2S)/2
        #
        # Since d=N-2S+1:
        #   = d
        #
        conj_sum_minus_true = (N + 1 + d - 2*S) // 2

        print("  conjugate sum - true sum")
        print(f"    = {conj_sum_minus_true}")
        print(f"    equals d = {conj_sum_minus_true == d}")
        print()

        if conj_sum_minus_true != d:
            local_ok = False

        # This is an important structural identity.
        print("BRANCH-SUM IDENTITY")
        print("  (z1+z2) - (p+q) = d")
        print(
            f"  actual   = {conj_sum_minus_true}"
        )
        print(
            f"  expected = {d}"
        )
        print(
            f"  identity = {conj_sum_minus_true == d}"
        )
        print()

        # ----------------------------------------------------------------------
        # Quartic factor interpretation
        # ----------------------------------------------------------------------

        print("QUARTIC DECOMPOSITION")

        # True quadratic:
        #   2Z^2 - (N+1-d)Z + 2N
        #
        # Conjugate:
        #   2Z^2 - (N+1+d)Z + 2N
        #
        # Their coefficient-sum structure:
        #
        #   coefficient sums cancel d
        #   coefficient differences isolate 2d
        #
        minus_linear = N + 1 - d
        plus_linear = N + 1 + d

        print(f"  minus linear coefficient magnitude = {minus_linear}")
        print(f"  plus  linear coefficient magnitude = {plus_linear}")
        print(f"  sum  = {minus_linear + plus_linear}")
        print(f"  expected 2(N+1) = {2*(N+1)}")
        print(f"  difference = {plus_linear - minus_linear}")
        print(f"  expected 2d = {2*d}")
        print()

        if minus_linear + plus_linear != 2*(N+1):
            local_ok = False

        if plus_linear - minus_linear != 2*d:
            local_ok = False

        # ----------------------------------------------------------------------
        # Arithmetic signatures
        # ----------------------------------------------------------------------

        print("ARITHMETIC SIGNATURE")
        print(f"  gcd(Delta_+, N)   = {gcd(delta_plus, N)}")
        print(f"  gcd(Delta_+, d)   = {gcd(delta_plus, d)}")
        print(f"  gcd(Delta_+, D)   = {gcd(delta_plus, D)}")
        print(f"  gcd(Delta_+, S)   = {gcd(delta_plus, S)}")
        print(f"  gcd(Delta_+, q-p) = {gcd(delta_plus, q-p)}")
        print()

        # ----------------------------------------------------------------------
        # Size
        # ----------------------------------------------------------------------

        print("SIZE DATA")
        print(f"  N digits            = {digits(N)}")
        print(f"  K digits            = {digits(K)}")
        print(f"  d digits            = {digits(d)}")
        print(f"  Delta_+ digits      = {digits(delta_plus)}")
        print(f"  Delta_+ bit-length  = {bit_length(delta_plus)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

        global_ok = global_ok and local_ok

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("TRUE BRANCH")
    print("  Q_-(Z) has roots p,q")
    print()

    print("CONJUGATE BRANCH")
    print("  Q_+(Z) is obtained from Q_-(Z) by d -> -d")
    print("  z1*z2 = N")
    print("  z1+z2 = (N+1+d)/2")
    print()

    print("BRANCH-SUM RELATION")
    print("  (z1+z2) - (p+q) = d")
    print()

    print("QUARTIC INTERPRETATION")
    print("  P(Z) = Q_-(Z) Q_+(Z)")
    print("  the two branches have the same product N")
    print("  their root sums differ by exactly d")
    print()

    print("CORE QUESTION")
    print(
        "  Does the conjugate branch admit a simple transformation "
        "to p,q or to another familiar factor quantity?"
    )
    print()

    print("=" * 110)
    print("EXPERIMENT 406 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 406 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
