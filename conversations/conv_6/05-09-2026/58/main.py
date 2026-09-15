# ================================================================
# START EXPERIMENT 141
# Exact double-radix cancellation / mixed-difference structure
# ================================================================

import math
import random
import time

try:
    import numpy as np
except ImportError:
    np = None


R1 = [17, 43, 59]
R2 = [19, 37, 61]

BIT_SIZES = [30, 36, 42, 48, 54]


# ------------------------------------------------
# Miller-Rabin
# ------------------------------------------------

def is_probable_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37
    ]

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic for the ranges we use here.
    witnesses = [2, 3, 5, 7, 11, 13, 17]

    for a in witnesses:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


def random_prime(bits):
    while True:
        x = random.getrandbits(bits)
        x |= (1 << (bits - 1))
        x |= 1

        if is_probable_prime(x):
            return x


def random_semiprime(bits):
    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:
        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p != q:
            return p, q


# ------------------------------------------------
# Utility
# ------------------------------------------------

def fmt_matrix(M):
    rows = []
    for row in M:
        rows.append(
            "["
            + ", ".join(str(x) for x in row)
            + "]"
        )
    return "\n".join(rows)


def matrix_norm(M):
    if np is not None:
        return float(np.linalg.norm(np.array(M, dtype=float)))

    s = 0.0
    for row in M:
        for x in row:
            s += float(x) * float(x)

    return math.sqrt(s)


def singular_values(M):
    if np is None:
        return None

    A = np.array(M, dtype=float)
    return np.linalg.svd(A, compute_uv=False)


def rank3_minor_values(M):
    """
    For a 3x3 matrix return all 2x2 minors.
    """
    vals = []

    for i in range(3):
        for j in range(i + 1, 3):
            for k in range(3):
                for l in range(k + 1, 3):

                    v = (
                        M[i][k] * M[j][l]
                        - M[i][l] * M[j][k]
                    )

                    vals.append(v)

    return vals


def max_abs(M):
    return max(abs(x) for row in M for x in row)


# ================================================================
# EXPERIMENT
# ================================================================

print("=" * 72)
print("START EXPERIMENT 141")
print("Exact double-radix cancellation / mixed-difference structure")
print("=" * 72)
print()
print("R1 =", R1)
print("R2 =", R2)
print()


for bits in BIT_SIZES:

    print("=" * 72)
    print(f"GENERATING {bits} -BIT SEMIPRIME")
    print("=" * 72)

    p, q = random_semiprime(bits)
    n = p * q

    print()
    print("-" * 72)
    print(f"n bits = {n.bit_length()}")
    print(f"n      = {n}")
    print(f"true p = {p}")
    print(f"true q = {q}")
    print()

    # ------------------------------------------------------------
    # Factor decomposition in every radix.
    # ------------------------------------------------------------

    k = []
    a = []

    for r in R1:
        ki = p // r
        ai = p % r

        k.append(ki)
        a.append(ai)

    ell = []
    b = []

    for s in R2:
        lj = q // s
        bj = q % s

        ell.append(lj)
        b.append(bj)

    print("TRUE k   =", k)
    print("TRUE a   =", a)
    print("TRUE ell =", ell)
    print("TRUE b   =", b)
    print()

    # ------------------------------------------------------------
    # Quotient and carry matrices.
    # ------------------------------------------------------------

    Q = []
    E = []

    for i, r in enumerate(R1):

        qrow = []
        erow = []

        for j, s in enumerate(R2):

            Rij = r * s

            qij = n // Rij
            eij = qij - k[i] * ell[j]

            qrow.append(qij)
            erow.append(eij)

        Q.append(qrow)
        E.append(erow)

    print("Q MATRIX")
    print(fmt_matrix(Q))
    print()

    print("E MATRIX")
    print(fmt_matrix(E))
    print()

    # ------------------------------------------------------------
    # First-order radix cancellation:
    #
    # X[i][j] = r_i Q_ij - r_0 Q_0j
    #
    # This removes the dominant K term in the row direction.
    # ------------------------------------------------------------

    X = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            value = (
                r * Q[i][j]
                - R1[0] * Q[0][j]
            )

            row.append(value)

        X.append(row)

    print("FIRST-ORDER ROW-CANCELLED MATRIX X")
    print(fmt_matrix(X))
    print()

    # ------------------------------------------------------------
    # First-order column cancellation.
    #
    # Y[i][j] = s_j Q_ij - s_0 Q_i0
    # ------------------------------------------------------------

    Y = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            value = (
                s * Q[i][j]
                - R2[0] * Q[i][0]
            )

            row.append(value)

        Y.append(row)

    print("FIRST-ORDER COLUMN-CANCELLED MATRIX Y")
    print(fmt_matrix(Y))
    print()

    # ------------------------------------------------------------
    # Exact double cancellation.
    #
    # D[i][j] =
    #     r_i s_j Q_ij
    #   - r_i s_0 Q_i0
    #   - r_0 s_j Q_0j
    #   + r_0 s_0 Q_00
    #
    # ------------------------------------------------------------

    D = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            value = (
                r * s * Q[i][j]
                - r * R2[0] * Q[i][0]
                - R1[0] * s * Q[0][j]
                + R1[0] * R2[0] * Q[0][0]
            )

            row.append(value)

        D.append(row)

    print("DOUBLE-RADIX ANNIHILATOR D")
    print(fmt_matrix(D))
    print()

    # ------------------------------------------------------------
    # Exact residue-difference outer product.
    #
    # A_i = a_0 - a_i
    # B_j = b_0 - b_j
    #
    # P[i][j] = A_i B_j
    #
    # ------------------------------------------------------------

    A = [
        a[0] - ai
        for ai in a
    ]

    B = [
        b[0] - bj
        for bj in b
    ]

    P = []

    for Ai in A:

        row = []

        for Bj in B:
            row.append(Ai * Bj)

        P.append(row)

    print("RESIDUE-DIFFERENCE OUTER PRODUCT P = A B^T")
    print("A =", A)
    print("B =", B)
    print(fmt_matrix(P))
    print()

    # ------------------------------------------------------------
    # Carry-only double-cancellation term.
    #
    # G = D - P
    #
    # Exact identity:
    #
    # D = P + G
    #
    # ------------------------------------------------------------

    G = []

    for i in range(3):

        row = []

        for j in range(3):

            gij = D[i][j] - P[i][j]
            row.append(gij)

        G.append(row)

    print("CARRY DOUBLE-CANCELLATION MATRIX G = D - P")
    print(fmt_matrix(G))
    print()

    # ------------------------------------------------------------
    # Exact identity check.
    # ------------------------------------------------------------

    identity_ok = True

    for i in range(3):
        for j in range(3):

            lhs = D[i][j]
            rhs = P[i][j] + G[i][j]

            if lhs != rhs:
                identity_ok = False

    print(
        "EXACT D = P + G IDENTITY:",
        "PASS" if identity_ok else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Direct remainder interpretation.
    #
    # Let
    #
    # t_ij = n mod (r_i s_j)
    #
    # Since
    #
    # r_i s_j Q_ij = n - t_ij,
    #
    # the entire D matrix can also be written as a mixed
    # difference of known remainders.
    # ------------------------------------------------------------

    T = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            t = n % (r * s)
            row.append(t)

        T.append(row)

    D_from_T = []

    for i in range(3):

        row = []

        for j in range(3):

            value = (
                T[i][0]
                + T[0][j]
                - T[i][j]
                - T[0][0]
            )

            row.append(value)

        D_from_T.append(row)

    print("REMAINDER MATRIX T = n mod (r_i s_j)")
    print(fmt_matrix(T))
    print()

    print("D reconstructed only from remainders")
    print(fmt_matrix(D_from_T))
    print()

    remainder_identity_ok = (D == D_from_T)

    print(
        "D = MIXED DIFFERENCE OF REMAINDERS:",
        "PASS" if remainder_identity_ok else "FAIL"
    )
    print()

    # ------------------------------------------------------------
    # Compare magnitudes.
    # ------------------------------------------------------------

    print("MAGNITUDE COMPARISON")

    print(f"max |Q| = {max_abs(Q)}")
    print(f"max |E| = {max_abs(E)}")
    print(f"max |X| = {max_abs(X)}")
    print(f"max |Y| = {max_abs(Y)}")
    print(f"max |D| = {max_abs(D)}")
    print(f"max |P| = {max_abs(P)}")
    print(f"max |G| = {max_abs(G)}")
    print()

    # ------------------------------------------------------------
    # Singular values.
    #
    # If G unexpectedly becomes low rank, that would be highly
    # interesting.
    # ------------------------------------------------------------

    if np is not None:

        sv_Q = singular_values(Q)
        sv_E = singular_values(E)
        sv_X = singular_values(X)
        sv_Y = singular_values(Y)
        sv_D = singular_values(D)
        sv_P = singular_values(P)
        sv_G = singular_values(G)

        print("SINGULAR VALUES")
        print("Q:", sv_Q)
        print("E:", sv_E)
        print("X:", sv_X)
        print("Y:", sv_Y)
        print("D:", sv_D)
        print("P:", sv_P)
        print("G:", sv_G)
        print()

        print("RELATIVE SECOND SINGULAR VALUE")
        print(f"Q : {sv_Q[1] / sv_Q[0]:.12e}")
        print(f"E : {sv_E[1] / sv_E[0]:.12e}")
        print(f"X : {sv_X[1] / sv_X[0]:.12e}")
        print(f"Y : {sv_Y[1] / sv_Y[0]:.12e}")
        print(f"D : {sv_D[1] / sv_D[0]:.12e}")
        print(f"P : {sv_P[1] / sv_P[0]:.12e}")
        print(f"G : {sv_G[1] / sv_G[0]:.12e}")
        print()

    # ------------------------------------------------------------
    # Exact 2x2 minors.
    #
    # P must have rank 1, so all its 2x2 minors vanish.
    # Test whether D or G unexpectedly inherit rank structure.
    # ------------------------------------------------------------

    minors_Q = rank3_minor_values(Q)
    minors_E = rank3_minor_values(E)
    minors_D = rank3_minor_values(D)
    minors_P = rank3_minor_values(P)
    minors_G = rank3_minor_values(G)

    print("MAX ABSOLUTE 2x2 MINOR")

    print(
        "Q:",
        max(abs(x) for x in minors_Q)
    )

    print(
        "E:",
        max(abs(x) for x in minors_E)
    )

    print(
        "D:",
        max(abs(x) for x in minors_D)
    )

    print(
        "P:",
        max(abs(x) for x in minors_P)
    )

    print(
        "G:",
        max(abs(x) for x in minors_G)
    )

    print()

    # ------------------------------------------------------------
    # A deliberately small brute-force test over residue
    # difference vectors.
    #
    # The question here is NOT to factor n.
    #
    # It asks:
    #
    #   Does the observable D sharply constrain A and B?
    #
    # We enumerate all possible A-vectors and see whether the
    # resulting P=A B^T can be close to D for some admissible B.
    #
    # This is intentionally only a diagnostic.
    # ------------------------------------------------------------

    A_candidates = []

    for aa0 in range(-(R1[0] - 1), R1[0]):

        for aa1 in range(-(R1[1] - 1), R1[1]):

            for aa2 in range(-(R1[2] - 1), R1[2]):

                # A_i = a0 - ai, with a0 and ai both valid residues.
                #
                # Necessary condition:
                # there exists a0 such that ai = a0 - Ai is valid
                # for every radix.
                #
                valid = False

                for a0_test in range(R1[0]):

                    ai0 = a0_test - aa0
                    ai1 = a0_test - aa1
                    ai2 = a0_test - aa2

                    if (
                        0 <= aa0 <= R1[0] - 1
                        and
                        - (R1[1] - 1) <= aa1 <= R1[1] - 1
                        and
                        - (R1[2] - 1) <= aa2 <= R1[2] - 1
                    ):
                        # Only use exact residue feasibility when
                        # all implied residues are legal.
                        #
                        # a0 is allowed to differ from ai by A_i.
                        #
                        if (
                            0 <= ai0 < R1[0]
                            and
                            0 <= ai1 < R1[1]
                            and
                            0 <= ai2 < R1[2]
                        ):
                            valid = True
                            break

                if valid:
                    A_candidates.append(
                        (aa0, aa1, aa2)
                    )

    # Same construction for B.
    B_candidates = []

    for bb0 in range(-(R2[0] - 1), R2[0]):

        for bb1 in range(-(R2[1] - 1), R2[1]):

            for bb2 in range(-(R2[2] - 1), R2[2]):

                valid = False

                for b0_test in range(R2[0]):

                    bj0 = b0_test - bb0
                    bj1 = b0_test - bb1
                    bj2 = b0_test - bb2

                    if (
                        0 <= bj0 < R2[0]
                        and
                        0 <= bj1 < R2[1]
                        and
                        0 <= bj2 < R2[2]
                    ):
                        valid = True
                        break

                if valid:
                    B_candidates.append(
                        (bb0, bb1, bb2)
                    )

    true_A = tuple(A)
    true_B = tuple(B)

    print()
    print("RESIDUE-DIFFERENCE STATE COUNTS")
    print("A candidates =", len(A_candidates))
    print("B candidates =", len(B_candidates))
    print("true A =", true_A)
    print("true B =", true_B)
    print()

    print(
        "TRUE A PRESENT =",
        true_A in A_candidates
    )

    print(
        "TRUE B PRESENT =",
        true_B in B_candidates
    )

    print()

    # ------------------------------------------------------------
    # Test whether the double observable D itself gives a very
    # small residual after subtracting the true outer product.
    # ------------------------------------------------------------

    ratio = (
        matrix_norm(G) / matrix_norm(P)
        if matrix_norm(P) != 0
        else float("inf")
    )

    print("TRUE RESIDUE PRODUCT VS CARRY TERM")
    print(f"||P|| = {matrix_norm(P):.12e}")
    print(f"||G|| = {matrix_norm(G):.12e}")
    print(f"||G|| / ||P|| = {ratio:.12e}")
    print()

    # ------------------------------------------------------------
    # Runtime marker.
    # ------------------------------------------------------------

    print("-" * 72)
    print(f"COMPLETED {bits}-BIT CASE")
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 141")
print("=" * 72)
