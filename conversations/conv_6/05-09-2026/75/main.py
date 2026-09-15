# ================================================================
# START EXPERIMENT 157B
# CRT-compressed global C consistency
#
# Fixes Experiment 157:
#
# The previous experiment scanned every p mod M:
#
#     M = product(r_i) * product(s_j)
#
# This is unnecessarily expensive.
#
# Here we enumerate ONLY the possible residue vectors
#
#     a_i = p mod r_i
#
# for the selected rows.
#
# CRT reconstructs p mod R, where
#
#     R = product(r_i).
#
# Then:
#
#     q = n * p^{-1} mod S
#
# where
#
#     S = product(s_j).
#
# This reduces a block from O(R*S) states to O(R) states.
#
# We then test the complete C pattern on the selected block.
#
# ================================================================

import math
import random
import time


R1 = [
    17,
    43,
    59,
    71,
    83,
]

R2 = [
    19,
    37,
    61,
    73,
    89,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]

SAMPLES_PER_SIZE = 3


# ================================================================
# PRIMALITY
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in [
        2, 3, 5, 7, 11, 13, 17,
    ]:

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

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p != q:
            return p, q


# ================================================================
# CRT
# ================================================================

def crt_pairwise(
    residues,
    moduli,
):

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli
    ):

        # x + M*t == a (mod m)
        #
        # M*t == a-x (mod m)

        rhs = (
            a - x
        ) % m

        inv = pow(
            M % m,
            -1,
            m
        )

        t = (
            rhs * inv
        ) % m

        x += M * t
        M *= m

    return (
        x % M,
        M
    )


# ================================================================
# TRUE C MATRIX
# ================================================================

def build_true_C(
    p,
    q,
):

    k = []
    a = []

    for r in R1:

        k.append(
            p // r
        )

        a.append(
            p % r
        )

    ell = []
    b = []

    for s in R2:

        ell.append(
            q // s
        )

        b.append(
            q % s
        )

    C = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            beta = (
                (k[i] * b[j])
                % s
            )

            alpha = (
                (ell[j] * a[i])
                % r
            )

            numerator = (
                r * beta
                + s * alpha
                + a[i] * b[j]
            )

            cij = (
                numerator
                // (r * s)
            )

            row.append(cij)

        C.append(row)

    return C, a, b


# ================================================================
# PRODUCT
# ================================================================

def prod(values):

    x = 1

    for v in values:
        x *= v

    return x


# ================================================================
# ALL RESIDUE VECTORS
# ================================================================

def residue_vectors(
    moduli
):

    vector = [0] * len(moduli)

    def rec(
        pos
    ):

        if pos == len(moduli):

            yield tuple(vector)
            return

        m = moduli[pos]

        for x in range(m):

            vector[pos] = x

            yield from rec(
                pos + 1
            )

    yield from rec(0)


# ================================================================
# COMPUTE C CELL FROM GLOBAL p,q RESIDUES
# ================================================================

def cell_C_from_global(
    pmod_R,
    qmod_S,
    r,
    s,
):

    # ------------------------------------------------------------
    # Recover p modulo r*s.
    #
    # Since r divides R and s divides S, combine:
    #
    #     p = pmod_R (mod R)
    #
    # but pmod_R alone isn't p mod r*s.
    #
    # We instead reconstruct k mod s using the exact identity
    #
    #     p = r*k+a.
    #
    # To obtain k mod s we need p mod r*s.
    #
    # Therefore use the fact that q is forced by n modulo S,
    # and derive the missing local quotient residue through the
    # factorization congruence.
    #
    # For p mod r = a and q mod s = b, we have:
    #
    #     pq = n (mod r*s)
    #
    # Since a and b are known, derive the quotient residues.
    # ------------------------------------------------------------

    a = pmod_R % r
    b = qmod_S % s

    # We cannot recover k mod s from pmod_R alone.
    # The selected block intentionally gives only row residues.
    #
    # So this helper is replaced by the exact local consistency
    # calculation below.

    return a, b


# ================================================================
# TEST GLOBAL RESIDUE STATE
#
# Given:
#
#     p mod R
#
# compute:
#
#     q mod S = n * inverse(p) mod S
#
# Then C cells are evaluated by constructing the relevant local
# quotient residues using CRT over the row and column modulus.
#
# ================================================================

def evaluate_block_state(
    n,
    pmod,
    R,
    S,
    rows,
    cols,
    C,
):

    if math.gcd(
        pmod,
        R
    ) != 1:

        return False

    # ------------------------------------------------------------
    # Since p is only represented modulo R, we need p modulo r*s
    # for each cell.
    #
    # We can recover it because p is known modulo R, and the
    # missing modulo-s information is not determined.
    #
    # Therefore enumerate only the k_i mod s_j values induced by
    # the factorization congruence.
    #
    # For a selected block, we can instead enumerate the possible
    # b_j values and use pq=n (mod r_i*s_j).
    #
    # This keeps the state bounded by product(s_j).
    # ------------------------------------------------------------

    qmod = (
        n
        * pow(
            pmod,
            -1,
            S
        )
    ) % S

    for i in rows:

        r = R1[i]

        a = pmod % r

        for j in cols:

            s = R2[j]

            b = qmod % s

            nmod = n % (
                r * s
            )

            target = (
                nmod
                - a * b
            ) % (
                r * s
            )

            # ----------------------------------------------------
            # Need x = k_i mod s and y = ell_j mod r such that
            #
            #     r*b*x + s*a*y == target (mod r*s).
            #
            # Since:
            #
            #     r*b*x == target (mod s)
            #
            # and:
            #
            #     s*a*y == target (mod r)
            #
            # solve independently.
            # ----------------------------------------------------

            A1 = (
                r * b
            ) % s

            B1 = (
                target % s
            )

            g1 = math.gcd(
                A1,
                s
            )

            if B1 % g1 != 0:
                return False

            m1 = s // g1

            if m1 == 1:

                xs = range(g1)

            else:

                inv1 = pow(
                    A1 // g1,
                    -1,
                    m1
                )

                x0 = (
                    (B1 // g1)
                    * inv1
                ) % m1

                xs = [
                    (
                        x0
                        + t * m1
                    ) % s
                    for t in range(g1)
                ]

            A2 = (
                s * a
            ) % r

            B2 = (
                target % r
            )

            g2 = math.gcd(
                A2,
                r
            )

            if B2 % g2 != 0:
                return False

            m2 = r // g2

            if m2 == 1:

                ys = range(g2)

            else:

                inv2 = pow(
                    A2 // g2,
                    -1,
                    m2
                )

                y0 = (
                    (B2 // g2)
                    * inv2
                ) % m2

                ys = [
                    (
                        y0
                        + t * m2
                    ) % r
                    for t in range(g2)
                ]

            good = False

            for x in xs:

                beta = (
                    x * b
                ) % s

                for y in ys:

                    alpha = (
                        y * a
                    ) % r

                    cij = (
                        (
                            r * beta
                            + s * alpha
                            + a * b
                        )
                        // (r * s)
                    )

                    if cij == C[i][j]:

                        good = True
                        break

                if good:
                    break

            if not good:
                return False

    return True


# ================================================================
# BLOCK SEARCH
#
# Enumerate p modulo R for the selected rows.
#
# This is dramatically smaller than scanning p modulo R*S.
# ================================================================

def search_block(
    n,
    p,
    q,
    C,
    rows,
    cols,
):

    row_moduli = [
        R1[i]
        for i in rows
    ]

    col_moduli = [
        R2[j]
        for j in cols
    ]

    R = prod(
        row_moduli
    )

    S = prod(
        col_moduli
    )

    survivors = []

    tested = 0

    invertible = 0

    # ------------------------------------------------------------
    # True state.
    # ------------------------------------------------------------

    true_p = p % R

    true_q = q % S

    # ------------------------------------------------------------
    # Enumerate p residue vectors equivalently by pmod modulo R.
    #
    # There are R possibilities.
    # ------------------------------------------------------------

    for pmod in range(
        1,
        R
    ):

        tested += 1

        if math.gcd(
            pmod,
            R
        ) != 1:

            continue

        invertible += 1

        # --------------------------------------------------------
        # q modulo S forced by pq=n mod S.
        #
        # Important:
        # pmod is invertible modulo S only if we know p modulo S.
        # We do NOT know p modulo S from pmod.
        #
        # Therefore this direct reduction is still incomplete.
        #
        # We return here deliberately so Experiment 157B can
        # expose the missing degree of freedom instead of silently
        # making an invalid inference.
        # --------------------------------------------------------

        # No valid global q reconstruction from p mod R alone.

        break

    # Signal that this formulation is invalid.
    return (
        R,
        S,
        tested,
        invertible,
        [],
        true_p,
        true_q,
    )


# ================================================================
# MAIN
# ================================================================

print("=" * 72)
print("START EXPERIMENT 157B")
print("CRT-compressed global C consistency")
print("=" * 72)
print()

print(
    "IMPORTANT: This version deliberately refuses to make the"
)
print(
    "invalid inference q = n*p^{-1} mod S from p mod R alone."
)
print()

print("=" * 72)
print("MATHEMATICAL CHECK")
print("=" * 72)
print()

print(
    "A p-residue modulo R does NOT determine p modulo S."
)

print(
    "Therefore q modulo S cannot be recovered from p modulo R."
)

print(
    "The full p modulo R*S state used by Experiment 157 was"
)

print(
    "actually necessary for that particular global parameterization."
)

print()

print(
    "The correct compression must instead exploit the shared"
)

print(
    "k_i / ell_j quotient variables directly."
)

print()

print("=" * 72)
print("FINISHED EXPERIMENT 157B")
print("=" * 72)
