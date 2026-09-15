#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 637
==========================================================================================

SYMMETRIC / DISCRIMINANT gcd DEPTH THEOREM

Goal
----

Experiment 636 identified the following candidate symmetric construction.

Let

    S = p + q
    D = p - q
    n = p q
    Delta = D^2 = S^2 - 4n

FRAME A:

    U = S
    R = Delta - 36
      = S^2 - 4n - 36

FRAME B:

    U = S - 2
    R = Delta - 16
      = S^2 - 4n - 16

Candidate law:

    depth = v2(gcd(U,R))

The purpose of this experiment is NOT to search arbitrary formulas.

Instead it verifies, exactly and separately:

    1. U/R identities
    2. v2(gcd(U,R)) == depth
    3. relation to d = gcd(A,B)
    4. relation to gcd(X,Y)
    5. normalized coprime residual proof
    6. equal/unequal valuation branches
    7. degenerate A=0 and B=0 branches
    8. symbolic algebraic identities
    9. possible odd gcd contamination
   10. exact final theorem

No level-by-level traversal is required for the main theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from collections import Counter, defaultdict


# ==============================================================================
# CONSTANTS
# ==============================================================================

INF = 10**9

PRIME_LIMIT = 6000

# Use enough primes to reproduce the existing ~329k semiprime domain.
# Adjust upward if necessary, but keep this deterministic.
EXAMPLE_LIMIT = 20


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str

    A: int
    B: int

    X: int
    Y: int

    depth: int


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def v2(x: int) -> int:
    """
    2-adic valuation.

    v2(0) is represented by INF.
    """
    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    """
    Positive odd part of x.
    """
    x = abs(x)

    if x == 0:
        return 0

    return x >> v2(x)


def normalize_pair(A: int, B: int):
    """
    Return:

        d = gcd(A,B)
        a = A/d
        b = B/d

    for nonzero pair.

    For degenerate pair, return None.
    """
    if A == 0 and B == 0:
        return None

    d = gcd(abs(A), abs(B))

    a = A // d
    b = B // d

    return d, a, b


def gcd0(a: int, b: int) -> int:
    """
    Symmetric gcd convention.
    """
    return gcd(abs(a), abs(b))


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        sieve[0] = 0

    if limit >= 1:
        sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# FRAME
# ==============================================================================

def frame_from_n(n: int) -> str:
    """
    Existing convention:

        n == 1 mod 4 -> Frame A
        n == 3 mod 4 -> Frame B
    """
    r = n & 3

    if r == 1:
        return "A"

    if r == 3:
        return "B"

    raise ValueError(
        f"Unexpected odd semiprime residue n mod 4={r}"
    )


# ==============================================================================
# RAW FACTOR RESIDUALS
# ==============================================================================

def residuals(frame: str, p: int, q: int):
    if frame == "A":
        return p - 3, q + 3

    if frame == "B":
        return p + 1, q - 3

    raise ValueError(frame)


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(frame: str, A: int, B: int):
    if frame == "A":
        num_x = B - A
        num_y = A + B

    elif frame == "B":
        num_x = 3 * A - B
        num_y = 3 * A + B

    else:
        raise ValueError(frame)

    if num_x % 2 != 0 or num_y % 2 != 0:
        raise ArithmeticError(
            f"Non-integral global X,Y: frame={frame} A={A} B={B}"
        )

    return num_x // 2, num_y // 2


# ==============================================================================
# TRUE GLOBAL DEPTH
# ==============================================================================

def global_depth(X: int, Y: int) -> int:
    """
    Established global law:

        depth = 1 + v2(gcd(X,Y))

    Wait:
        Previous experiments established exactly:

            depth = 1 + v2(gcd(X,Y))

    However examples show n=21:

            gcd(X,Y)=4
            v2=2
            depth=3.

    Therefore this function uses that established law.
    """
    g = gcd0(X, Y)

    if g == 0:
        raise ArithmeticError("gcd(X,Y)=0")

    return 1 + v2(g)


# ==============================================================================
# SYMMETRIC / DISCRIMINANT CONSTRUCTION
# ==============================================================================

def symmetric_data(n: int, p: int, q: int, frame: str):
    S = p + q
    D = p - q

    Delta = D * D

    if frame == "A":
        U = S
        R = Delta - 36
        sigma = 0
        kappa = 36

    elif frame == "B":
        U = S - 2
        R = Delta - 16
        sigma = 2
        kappa = 16

    else:
        raise ValueError(frame)

    # Equivalent discriminant-only formula.
    R2 = S * S - 4 * n - kappa

    return {
        "S": S,
        "D": D,
        "Delta": Delta,
        "U": U,
        "R": R,
        "R2": R2,
        "sigma": sigma,
        "kappa": kappa,
    }


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(primes: list[int]) -> list[State]:
    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:

            n = p * q

            # Only odd semiprimes.
            if (n & 1) == 0:
                continue

            frame = frame_from_n(n)

            A, B = residuals(
                frame,
                p,
                q,
            )

            X, Y = global_xy(
                frame,
                A,
                B,
            )

            depth = global_depth(
                X,
                Y,
            )

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    X=X,
                    Y=Y,
                    depth=depth,
                )
            )

    return states


# ==============================================================================
# TEST FRAME
# ==============================================================================

def test_frame_domain(states):
    print("=" * 90)
    print("TEST 0: FRAME DOMAIN")
    print("=" * 90)

    failures = 0

    checked = 0

    for s in states:
        checked += 1

        if s.frame == "A":
            expected = 1
        elif s.frame == "B":
            expected = 3
        else:
            failures += 1
            continue

        if s.n % 4 != expected:
            failures += 1

    print(
        f"checked={checked} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 1:
# EXACT U/R IDENTITIES
# ==============================================================================

def test_exact_ur_identities(states):
    print("=" * 90)
    print("TEST 1: EXACT U/R IDENTITIES")
    print("=" * 90)

    failures = 0

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        S = info["S"]
        D = info["D"]

        U = info["U"]
        R = info["R"]
        R2 = info["R2"]

        # Discriminant identity.
        if D * D != S * S - 4 * s.n:
            failures += 1

        # R definition.
        if s.frame == "A":
            if U != S:
                failures += 1

            if R != D * D - 36:
                failures += 1

        elif s.frame == "B":
            if U != S - 2:
                failures += 1

            if R != D * D - 16:
                failures += 1

        # Quadratic reconstruction.
        if R != R2:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2:
# CANDIDATE DEPTH LAW
# ==============================================================================

def test_candidate_depth(states):
    print("=" * 90)
    print("TEST 2: DEPTH = v2(gcd(U,R))")
    print("=" * 90)

    failures = 0
    shown = 0

    distribution = Counter()

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        g = gcd0(U, R)

        if g == 0:
            failures += 1
            continue

        predicted = v2(g)

        distribution[predicted] += 1

        if predicted != s.depth:

            failures += 1

            if shown < 20:
                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"U={U}",
                    f"R={R}",
                    f"gcd={g}",
                    f"v2g={predicted}",
                    f"depth={s.depth}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("depth/gcd valuation distribution:")

    for k in sorted(distribution):
        print(
            f"    {k}: {distribution[k]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 3:
# EXACT RELATION TO gcd(X,Y)
# ==============================================================================

def test_gcd_xy_relation(states):
    print("=" * 90)
    print("TEST 3: SYMMETRIC gcd vs gcd(X,Y)")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        g_sym = gcd0(U, R)
        g_xy = gcd0(s.X, s.Y)

        if v2(g_sym) != v2(g_xy):
            failures += 1

            if shown < 20:
                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"g_sym={g_sym}",
                    f"g_xy={g_xy}",
                    f"v2sym={v2(g_sym)}",
                    f"v2xy={v2(g_xy)}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 4:
# SYMMETRIC gcd / RESIDUAL gcd
# ==============================================================================

def test_residual_gcd_relation(states):
    print("=" * 90)
    print("TEST 4: gcd(U,R) vs d=gcd(A,B)")
    print("=" * 90)

    failures = 0

    classes = Counter()
    examples = {}

    for s in states:

        if s.A == 0 and s.B == 0:
            continue

        d = gcd(abs(s.A), abs(s.B))

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        g = gcd0(
            info["U"],
            info["R"],
        )

        alpha = v2(s.A)
        beta = v2(s.B)

        relation = (
            "equal-v2"
            if alpha == beta
            else "unequal-v2"
        )

        ratio_key = None

        if d != 0 and g % d == 0:
            ratio_key = g // d

        classes[
            (s.frame, relation, ratio_key)
        ] += 1

        examples.setdefault(
            (s.frame, relation, ratio_key),
            s,
        )

        # We do NOT expect the exact integer gcd to equal d.
        #
        # We only expect:
        #
        #     v2(g) = depth
        #
        # and:
        #
        #     v2(g)
        #       =
        #     v2(d) + equality.
        #
        if v2(g) != v2(d) + (alpha == beta):
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    symmetric-gcd relationship classes:")

    for key in sorted(classes, key=str):
        print(
            f"        {key}: {classes[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 5:
# NORMALIZED COPRIME RESIDUAL DERIVATION
# ==============================================================================

def test_normalized_coprime_pair(states):
    print("=" * 90)
    print("TEST 5: NORMALIZED COPRIME RESIDUAL PAIR")
    print("=" * 90)

    failures = 0

    parity_counts = Counter()

    for s in states:

        if s.A == 0 or s.B == 0:
            continue

        d, a, b = normalize_pair(
            s.A,
            s.B,
        )

        if gcd(abs(a), abs(b)) != 1:
            failures += 1
            continue

        parity = (
            a & 1,
            b & 1,
        )

        parity_counts[
            (s.frame, parity)
        ] += 1

        # (0,0) impossible for coprime normalized pair.
        if parity == (0, 0):
            failures += 1

    print(
        f"checked={sum(parity_counts.values())} failures={failures}"
    )

    print()
    print("    normalized parity classes:")

    for key in sorted(
        parity_counts,
        key=str,
    ):
        print(
            f"        {key}: {parity_counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 6:
# DERIVE U/R IN TERMS OF d,a,b
# ==============================================================================

def test_symbolic_normalized_form(states):
    print("=" * 90)
    print("TEST 6: NORMALIZED U/R ALGEBRA")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        if s.A == 0 or s.B == 0:
            continue

        d, a, b = normalize_pair(
            s.A,
            s.B,
        )

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        if s.frame == "A":

            # A=d*a
            # B=d*b
            #
            # S=A+B
            #
            # therefore
            #
            # U=d(a+b)
            #
            # D=B-A=d(b-a)
            #
            # R=D^2-36
            #  = d^2(b-a)^2 - 36

            expected_U = d * (a + b)

            expected_R = (
                d * d * (b - a) ** 2
                - 36
            )

        else:

            # B-frame:
            #
            # U=S-2
            #
            # S = p+q
            #
            # A=p+1
            # B=q-3
            #
            # p=A-1
            # q=B+3
            #
            # S=A+B+2
            #
            # so U=A+B=d(a+b).
            #
            # D=p-q
            #  = A-B-4
            #
            # therefore:
            #
            # R=(A-B-4)^2-16
            #
            #   = (A-B)(A-B-8)
            #
            #   = d(a-b)*(d(a-b)-8)

            expected_U = d * (a + b)

            expected_R = (
                (d * (a - b) - 4) ** 2
                - 16
            )

        if U != expected_U:
            failures += 1

            if shown < 10:
                print(
                    "    U mismatch:",
                    s,
                    expected_U,
                )

                shown += 1

        if R != expected_R:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 7:
# EXACT PARITY -> v2(gcd(U,R))
# ==============================================================================

def test_parity_formula(states):
    print("=" * 90)
    print("TEST 7: PARITY -> SYMMETRIC gcd VALUATION")
    print("=" * 90)

    failures = 0

    counts = Counter()

    for s in states:

        if s.A == 0 or s.B == 0:
            continue

        d, a, b = normalize_pair(
            s.A,
            s.B,
        )

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        g = gcd0(
            info["U"],
            info["R"],
        )

        expected = (
            v2(d)
            + (a % 2 == 1 and b % 2 == 1)
        )

        actual = v2(g)

        counts[
            (s.frame, a & 1, b & 1)
        ] += 1

        if actual != expected:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    exact parity law:")

    for key in sorted(
        counts,
        key=str,
    ):
        print(
            f"        {key}: {counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 8:
# DEGENERATE A=0
# ==============================================================================

def test_A_zero(states):
    print("=" * 90)
    print("TEST 8: FRAME-A A=0 BRANCH")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        if s.frame != "A" or s.A != 0:
            continue

        checked += 1

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        g = gcd0(
            info["U"],
            info["R"],
        )

        # Here:
        #
        # A=0
        # B=q+3
        #
        # X=B/2
        # Y=B/2
        #
        # gcd(X,Y)=|B|/2
        #
        # depth=1+v2(B/2)=v2(B).
        #
        # U=B
        #
        # R=B^2-36
        #
        # gcd(U,R) has the same v2 as B.

        if v2(g) != v2(s.B):
            failures += 1

    print(
        f"checked={checked} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 9:
# DEGENERATE B=0
# ==============================================================================

def test_B_zero(states):
    print("=" * 90)
    print("TEST 9: FRAME-B B=0 BRANCH")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        if s.frame != "B" or s.B != 0:
            continue

        checked += 1

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        g = gcd0(
            info["U"],
            info["R"],
        )

        # Expected:
        #
        # p=q=3
        # n=9
        #
        # U=p+q-2=4
        #
        # R=(p-q)^2-16=-16
        #
        # gcd(U,R)=4
        #
        # depth=2.

        if v2(g) != s.depth:
            failures += 1

    print(
        f"checked={checked} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 10:
# ODD PART OF SYMMETRIC gcd
# ==============================================================================

def test_odd_part(states):
    print("=" * 90)
    print("TEST 10: ODD PART OF gcd(U,R)")
    print("=" * 90)

    """
    This test deliberately does NOT assume that the entire integer
    gcd(U,R) equals gcd(X,Y).

    We only verify that the odd part is irrelevant to depth:

        depth = v2(gcd(U,R)).

    We report the odd-part distribution to see whether a cleaner
    exact integer gcd identity is hiding here.
    """

    odd_parts = Counter()
    failures = 0

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        g = gcd0(
            info["U"],
            info["R"],
        )

        odd = odd_part(g)

        odd_parts[
            (s.frame, odd)
        ] += 1

        if v2(g) != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    most common odd parts:")

    common = odd_parts.most_common(30)

    for (frame, odd), count in common:
        print(
            f"        frame={frame} odd={odd} count={count}"
        )

    print()

    return failures


# ==============================================================================
# TEST 11:
# DISCRIMINANT-ONLY FORM
# ==============================================================================

def test_discriminant_only(states):
    print("=" * 90)
    print("TEST 11: DISCRIMINANT-ONLY CONSTRUCTION")
    print("=" * 90)

    failures = 0

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        # R is reconstructed from S and n only.
        #
        # No explicit D is needed.
        R_from_symmetric = (
            info["S"] * info["S"]
            - 4 * s.n
            - info["kappa"]
        )

        if R != R_from_symmetric:
            failures += 1
            continue

        g = gcd0(
            U,
            R_from_symmetric,
        )

        if v2(g) != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 12:
# EXACT FINAL THEOREM
# ==============================================================================

def test_final_theorem(states):
    print("=" * 90)
    print("TEST 12: FINAL SYMMETRIC DEPTH THEOREM")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        symmetric_depth = v2(
            gcd0(U, R)
        )

        global_depth_value = global_depth(
            s.X,
            s.Y,
        )

        if symmetric_depth != global_depth_value:
            failures += 1

            if shown < 20:
                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"U={U}",
                    f"R={R}",
                    f"v2g={symmetric_depth}",
                    f"global={global_depth_value}",
                    f"stored={s.depth}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(states):
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = [
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        77,
        93,
        141,
        183,
        213,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = by_n.get(n)

        if s is None:
            continue

        info = symmetric_data(
            s.n,
            s.p,
            s.q,
            s.frame,
        )

        U = info["U"]
        R = info["R"]

        g = gcd0(
            U,
            R,
        )

        d = gcd0(
            s.A,
            s.B,
        )

        print()
        print(
            f"n={s.n} p={s.p} q={s.q} frame={s.frame}"
        )

        print(
            f"    A={s.A} B={s.B}"
        )

        print(
            f"    X={s.X} Y={s.Y}"
        )

        print(
            f"    d=gcd(A,B)={d}"
        )

        print(
            f"    S={info['S']} D={info['D']}"
        )

        print(
            f"    Delta={info['Delta']}"
        )

        print(
            f"    U={U}"
        )

        print(
            f"    R={R}"
        )

        print(
            f"    gcd(U,R)={g}"
        )

        print(
            f"    v2(gcd(U,R))={v2(g)}"
        )

        print(
            f"    v2(gcd(X,Y))={v2(gcd0(s.X,s.Y))}"
        )

        print(
            f"    depth={s.depth}"
        )


# ==============================================================================
# FINAL SYMBOLIC SUMMARY
# ==============================================================================

def print_summary():
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
The tested symmetric construction is:

    S = p + q
    Delta = (p-q)^2
          = S^2 - 4n.

FRAME A:

    U = S
    R = Delta - 36

      = S^2 - 4n - 36.

FRAME B:

    U = S - 2
    R = Delta - 16

      = S^2 - 4n - 16.


Candidate theorem:

    depth = v2(gcd(U,R)).

Equivalently:

FRAME A:

    depth =
        v2(
            gcd(
                p+q,
                (p-q)^2 - 36
            )
        )

    = v2(
        gcd(
            p+q,
            (p+q)^2 - 4n - 36
        )
      ).

FRAME B:

    depth =
        v2(
            gcd(
                p+q-2,
                (p-q)^2 - 16
            )
        )

    = v2(
        gcd(
            p+q-2,
            (p+q)^2 - 4n - 16
        )
      ).


Residual interpretation:

    d = gcd(A,B)

and for nonzero A,B:

    depth =
        v2(d)
        + [v2(A)=v2(B)].

The symmetric gcd should therefore package both pieces:

    v2(d)
        +
    equality bit.

This experiment does not search arbitrary n+c forms.

It directly tests whether that information is exactly
encoded by one symmetric gcd involving:

    sum
    discriminant
    n
    frame.


IMPORTANT:

The result of this experiment is NOT assumed in advance.

If TEST 2 passes for all states, then the candidate becomes
an exact experimental theorem over the complete generated
domain.

If TEST 2 or TEST 12 fails, the first counterexamples printed
will show exactly which branch breaks.


Potential next step after a PASS:

    reduce p+q itself.

The remaining problem would then be:

    Can S=p+q be obtained from the available invariant
    structure without first knowing p and q?

That is a much more focused question than the previous
blind n-only valuation searches.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 637 START")
    print("=" * 90)
    print()
    print("SYMMETRIC / DISCRIMINANT gcd DEPTH THEOREM")
    print()

    # --------------------------------------------------------------------------
    # PRIME SIEVE
    # --------------------------------------------------------------------------

    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    # --------------------------------------------------------------------------
    # SEMIPRIME GENERATION
    # --------------------------------------------------------------------------

    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    # --------------------------------------------------------------------------
    # TESTS
    # --------------------------------------------------------------------------

    total_failures = 0

    total_failures += test_frame_domain(
        states
    )

    total_failures += test_exact_ur_identities(
        states
    )

    total_failures += test_candidate_depth(
        states
    )

    total_failures += test_gcd_xy_relation(
        states
    )

    total_failures += test_residual_gcd_relation(
        states
    )

    total_failures += test_normalized_coprime_pair(
        states
    )

    total_failures += test_symbolic_normalized_form(
        states
    )

    total_failures += test_parity_formula(
        states
    )

    total_failures += test_A_zero(
        states
    )

    total_failures += test_B_zero(
        states
    )

    total_failures += test_odd_part(
        states
    )

    total_failures += test_discriminant_only(
        states
    )

    total_failures += test_final_theorem(
        states
    )

    # --------------------------------------------------------------------------
    # EXAMPLES
    # --------------------------------------------------------------------------

    print_examples(
        states
    )

    # --------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------

    print_summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 637 FINISHED")
    print("=" * 90)

    print()
    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:
        print(
            "STATUS=ALL TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
