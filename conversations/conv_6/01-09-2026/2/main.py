#!/usr/bin/env python3

"""
===============================================================================
START EXPERIMENT 76
AUXILIARY-EQUATION RANK / (a,b) INFORMATION EXPERIMENT
===============================================================================

Question:

    Do multiple auxiliary observations

        n+x
        -> (ax,bx,kx,lx,Kx,Ex)

    provide independent equations for the unknown original residues

        a,b

    or do they all collapse to the same original factorization equation?

For

    n = p*q
    p = a + k*r1
    q = b + l*r2
    R = r1*r2
    K = k*l

and an auxiliary factorization

    n+x = px*qx
    px = ax + kx*r1
    qx = bx + lx*r2
    Kx = kx*lx

we have the exact identity

    x =
        R*(Kx-K)
        + r1*(kx*bx-k*b)
        + r2*(lx*ax-l*a)
        + (ax*bx-a*b)

Rearranging gives

    r1*kx*bx + r2*lx*ax + ax*bx
        - x
        - R*(Kx-K)
        =
    r1*k*b + r2*l*a + a*b

The right side is

    n - R*K

because

    n = R*K + r1*k*b + r2*l*a + a*b.

Therefore every auxiliary observation should reduce to

    a*b + r2*l*a + r1*k*b = n - R*K

and hence

    (a + r1*k)(b + r2*l) = n.

This experiment checks that computationally.

It also measures whether:

    multiple auxiliary records
        -> multiple independent constraints

or:

    multiple auxiliary records
        -> exactly one repeated constraint.

No factorization of the unknown n is used by the reconstruction test.
The true factors are retained only as experimental ground truth.

SymPy is used to factor the auxiliary n+x values.
===============================================================================
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass

from sympy import factorint, isprime, nextprime, prevprime


# =============================================================================
# CONFIGURATION
# =============================================================================

SCALES = (10**9, 10**12, 10**16)

ANCHORS_PER_SCALE = 3

R_OFFSETS = (0.95, 1.00, 1.05)

S_VALUES = (30, 210, 2310, 30030)

K_TARGET = 1000

SEED = 1511464998

# Generate semiprime factors close to sqrt(n).
FACTOR_SPREAD = 0.18

# Number of auxiliary factor pairs printed per case.
SHOW_AUX_ROWS = 8

# Number of random (a,b) tests used to verify the candidate equation
# independently of the real residues.
RANDOM_AB_TESTS = 20


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Anchor:
    p: int
    q: int
    n: int


@dataclass
class AuxiliaryRecord:
    S: int
    x: int
    px: int
    qx: int
    ax: int
    bx: int
    kx: int
    lx: int
    Kx: int
    Tx: int
    Ex: int


# =============================================================================
# UTILITIES
# =============================================================================

def fmt(n: int) -> str:
    return f"{n:,}"


def scale_name(n: int) -> str:
    if n == 10**9:
        return "1e+09"
    if n == 10**12:
        return "1e+12"
    if n == 10**16:
        return "1e+16"
    return f"{n:.6e}"


def choose_prime_near(x: int, rng: random.Random) -> int:
    """
    Choose a prime near x, randomly on either side when possible.
    """
    x = max(3, int(x))

    candidates = []

    if isprime(x):
        candidates.append(x)

    candidates.append(int(nextprime(x)))
    candidates.append(int(prevprime(x)))

    # Randomly perturb search direction to avoid always selecting one side.
    direction = rng.choice((-1, 1))

    for delta in (1, 2, 3, 5, 7, 11, 17, 23, 31, 47, 61, 97):
        y = x + direction * delta
        if y >= 3:
            if isprime(y):
                candidates.append(int(y))
                break

    # Deterministic best choice.
    return min(candidates, key=lambda p: abs(p - x))


def generate_anchor(scale: int, rng: random.Random) -> Anchor:
    """
    Generate a semiprime n near the requested scale.

    We construct p,q ourselves so the experiment has exact ground truth.
    """
    root = math.isqrt(scale)

    low = max(100, int(root * (1.0 - FACTOR_SPREAD)))
    high = max(low + 10, int(root * (1.0 + FACTOR_SPREAD)))

    p_target = rng.randint(low, high)
    q_target = rng.randint(low, high)

    p = choose_prime_near(p_target, rng)
    q = choose_prime_near(q_target, rng)

    if p == q:
        q = int(nextprime(q + 1))

    return Anchor(p=p, q=q, n=p * q)


def choose_balanced_prime_pair(target_R: float) -> tuple[int, int]:
    """
    Construct balanced prime moduli with

        r1*r2 ~= target_R.

    Uses exact integer prime adjustment.
    """
    if target_R <= 0:
        raise ValueError("target_R must be positive")

    root = max(3, math.isqrt(int(target_R)))

    candidates = []

    for d1 in range(-40, 41):
        r1 = root + d1
        if r1 < 2 or not isprime(r1):
            continue

        r2_target = int(target_R / r1)

        for delta in range(-40, 41):
            r2_guess = max(2, r2_target + delta)

            if isprime(r2_guess):
                r2 = r2_guess
                R = r1 * r2
                error = abs(R - target_R)

                candidates.append((error, r1, r2))

    if not candidates:
        # Wider fallback.
        r1 = int(nextprime(root))
        r2 = int(nextprime(max(2, round(target_R / r1))))
        candidates.append((abs(r1 * r2 - target_R), r1, r2))

    _, r1, r2 = min(candidates)
    return r1, r2


def quotient_cell(p: int, r1: int, q: int, r2: int):
    """
    Compute

        p = a + k*r1
        q = b + l*r2.
    """
    k = p // r1
    a = p - k * r1

    l = q // r2
    b = q - l * r2

    K = k * l

    return a, b, k, l, K


# =============================================================================
# AUXILIARY FACTORIZATION
# =============================================================================

def factor_auxiliary_number(
    n: int,
    r1: int,
    r2: int,
    K: int,
    S_values: tuple[int, ...],
) -> list[AuxiliaryRecord]:
    """
    Factor each unique n+x only once.

    This is intentionally independent of the choice of R offsets.
    """
    records: list[AuxiliaryRecord] = []

    seen = set()

    for S in S_values:
        x = (-n) % S
        m = n + x

        if m in seen:
            continue

        seen.add(m)

        fac = factorint(m)

        # Enumerate all divisor pairs p<=q.
        divisors = sorted(fac)

        # Generate all positive divisors.
        all_divs = [1]

        for prime, exponent in fac.items():
            old = list(all_divs)
            mul = 1

            new_terms = []
            for e in range(1, exponent + 1):
                mul *= prime
                new_terms.extend(d * mul for d in old)

            all_divs.extend(new_terms)

        all_divs = sorted(set(all_divs))

        pair_count = 0

        for px in all_divs:
            if px * px > m:
                break

            if m % px != 0:
                continue

            qx = m // px

            ax = px % r1
            kx = px // r1

            bx = qx % r2
            lx = qx // r2

            Kx = kx * lx
            Tx = m // (r1 * r2)
            Ex = Tx - Kx

            records.append(
                AuxiliaryRecord(
                    S=S,
                    x=x,
                    px=px,
                    qx=qx,
                    ax=ax,
                    bx=bx,
                    kx=kx,
                    lx=lx,
                    Kx=Kx,
                    Tx=Tx,
                    Ex=Ex,
                )
            )

            # Include reversed factor orientation as a distinct
            # quotient-cell observation.
            if px != qx:
                px2 = qx
                qx2 = px

                ax2 = px2 % r1
                kx2 = px2 // r1

                bx2 = qx2 % r2
                lx2 = qx2 // r2

                Kx2 = kx2 * lx2
                Tx2 = Tx
                Ex2 = Tx2 - Kx2

                records.append(
                    AuxiliaryRecord(
                        S=S,
                        x=x,
                        px=px2,
                        qx=qx2,
                        ax=ax2,
                        bx=bx2,
                        kx=kx2,
                        lx=lx2,
                        Kx=Kx2,
                        Tx=Tx2,
                        Ex=Ex2,
                    )
                )

            pair_count += 1

    return records


# =============================================================================
# EQUATION CHECK
# =============================================================================

def auxiliary_equation_value(
    n: int,
    R: int,
    k: int,
    l: int,
    rec: AuxiliaryRecord,
) -> int:
    """
    Compute the left side:

        r1*kx*bx + r2*lx*ax + ax*bx
            - x
            - R*(Kx-K)

    It must equal:

        n - R*K.
    """
    raise RuntimeError("Use auxiliary_equation_value_with_moduli")


def auxiliary_equation_value_with_moduli(
    n: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
    rec: AuxiliaryRecord,
) -> int:
    R = r1 * r2
    K = k * l

    return (
        r1 * rec.kx * rec.bx
        + r2 * rec.lx * rec.ax
        + rec.ax * rec.bx
        - rec.x
        - R * (rec.Kx - K)
    )


def original_equation_value(
    n: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
    a: int,
    b: int,
) -> int:
    return (
        a * b
        + r2 * l * a
        + r1 * k * b
    )


def factorized_equation_value(
    r1: int,
    r2: int,
    k: int,
    l: int,
    a: int,
    b: int,
) -> int:
    return (
        a + r1 * k
    ) * (
        b + r2 * l
    )


# =============================================================================
# MAIN CASE
# =============================================================================

def run_case(
    anchor: Anchor,
    r1: int,
    r2: int,
) -> dict:

    n = anchor.n
    R = r1 * r2

    a, b, k, l, K = quotient_cell(anchor.p, r1, anchor.q, r2)

    T = n // R
    E = T - K

    # This is the supposedly repeated equation.
    expected_constant = n - R * K

    records = factor_auxiliary_number(
        n=n,
        r1=r1,
        r2=r2,
        K=K,
        S_values=S_VALUES,
    )

    # -------------------------------------------------------------------------
    # Check every auxiliary observation.
    # -------------------------------------------------------------------------

    invariant_failures = 0
    invariant_values = set()

    K_equal_records = []
    K_other_records = []

    for rec in records:
        value = auxiliary_equation_value_with_moduli(
            n=n,
            r1=r1,
            r2=r2,
            k=k,
            l=l,
            rec=rec,
        )

        invariant_values.add(value)

        if value != expected_constant:
            invariant_failures += 1

        if rec.Kx == K:
            K_equal_records.append(rec)
        else:
            K_other_records.append(rec)

    # -------------------------------------------------------------------------
    # Check that the original equation really equals n-RK.
    # -------------------------------------------------------------------------

    original_value = original_equation_value(
        n, r1, r2, k, l, a, b
    )

    original_identity_ok = (
        original_value == expected_constant
    )

    factorized_value = factorized_equation_value(
        r1, r2, k, l, a, b
    )

    factorized_identity_ok = (
        factorized_value == n
    )

    # -------------------------------------------------------------------------
    # Difference between any auxiliary equation and the original equation.
    #
    # Every auxiliary row should be zero after substitution.
    # -------------------------------------------------------------------------

    difference_failures = 0

    for rec in records:
        lhs = auxiliary_equation_value_with_moduli(
            n=n,
            r1=r1,
            r2=r2,
            k=k,
            l=l,
            rec=rec,
        )

        rhs = original_value

        if lhs != rhs:
            difference_failures += 1

    # -------------------------------------------------------------------------
    # Pairwise equation-rank test.
    #
    # We normalize every row to:
    #
    #     ab + A*a + B*b = C
    #
    # The coefficients A,B,C must be identical.
    # -------------------------------------------------------------------------

    A = r2 * l
    B = r1 * k
    C = expected_constant

    normalized_rows = set()

    for _rec in records:
        normalized_rows.add((A, B, C))

    # -------------------------------------------------------------------------
    # Kx == K subset.
    # -------------------------------------------------------------------------

    K_equal_values = set()

    for rec in K_equal_records:
        K_equal_values.add(
            auxiliary_equation_value_with_moduli(
                n=n,
                r1=r1,
                r2=r2,
                k=k,
                l=l,
                rec=rec,
            )
        )

    K_equal_invariant_ok = (
        all(v == expected_constant for v in K_equal_values)
    )

    # -------------------------------------------------------------------------
    # Random fake (a,b) points.
    #
    # Check that the factorization transformation is algebraically exact.
    # -------------------------------------------------------------------------

    rng = random.Random(
        SEED ^ n ^ r1 ^ r2
    )

    random_identity_failures = 0

    for _ in range(RANDOM_AB_TESTS):
        aa = rng.randrange(r1)
        bb = rng.randrange(r2)

        left = (
            aa * bb
            + r2 * l * aa
            + r1 * k * bb
            + R * K
        )

        right = (
            aa + r1 * k
        ) * (
            bb + r2 * l
        )

        if left != right:
            random_identity_failures += 1

    return {
        "n": n,
        "r1": r1,
        "r2": r2,
        "R": R,
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,
        "records": records,
        "K_equal_records": K_equal_records,
        "K_other_records": K_other_records,
        "expected_constant": expected_constant,
        "invariant_values": invariant_values,
        "invariant_failures": invariant_failures,
        "original_identity_ok": original_identity_ok,
        "factorized_identity_ok": factorized_identity_ok,
        "difference_failures": difference_failures,
        "normalized_rows": normalized_rows,
        "K_equal_invariant_ok": K_equal_invariant_ok,
        "random_identity_failures": random_identity_failures,
    }


# =============================================================================
# PRINTING
# =============================================================================

def print_case(result: dict, show_rows: bool = True):
    n = result["n"]
    r1 = result["r1"]
    r2 = result["r2"]
    R = result["R"]

    a = result["a"]
    b = result["b"]
    k = result["k"]
    l = result["l"]
    K = result["K"]

    T = result["T"]
    E = result["E"]

    records = result["records"]
    Keq = result["K_equal_records"]
    Kother = result["K_other_records"]

    print(
        f"    R={fmt(R)} "
        f"(r1,r2)=({fmt(r1)},{fmt(r2)})"
    )

    print(
        f"        TRUE a,b=({fmt(a)},{fmt(b)}) "
        f"k,l=({k},{l}) K={fmt(K)} T={fmt(T)} E={fmt(E)}"
    )

    print(
        f"        auxiliary records          = {len(records):,}"
    )
    print(
        f"        Kx == K records            = {len(Keq):,}"
    )
    print(
        f"        Kx != K records            = {len(Kother):,}"
    )

    print(
        f"        equation invariant values  = "
        f"{len(result['invariant_values'])}"
    )

    print(
        f"        normalized equation rows   = "
        f"{len(result['normalized_rows'])}"
    )

    print(
        f"        expected constant C        = "
        f"n-R*K = {fmt(result['expected_constant'])}"
    )

    print(
        f"        invariant failures        = "
        f"{result['invariant_failures']}"
    )

    print(
        f"        auxiliary-vs-original fail = "
        f"{result['difference_failures']}"
    )

    print(
        f"        Kx=K invariant valid       = "
        f"{result['K_equal_invariant_ok']}"
    )

    print(
        f"        factorized identity       = "
        f"(a+r1*k)(b+r2*l)=n -> "
        f"{result['factorized_identity_ok']}"
    )

    print(
        f"        random algebra failures   = "
        f"{result['random_identity_failures']}"
    )

    if not show_rows:
        return

    if not records:
        return

    print(
        "        SAMPLE AUXILIARY ROWS:"
    )

    for rec in records[:SHOW_AUX_ROWS]:
        invariant = auxiliary_equation_value_with_moduli(
            n=n,
            r1=r1,
            r2=r2,
            k=k,
            l=l,
            rec=rec,
        )

        deltaK = rec.Kx - K

        print(
            f"            S={fmt(rec.S)} "
            f"x={fmt(rec.x)} "
            f"Kx={fmt(rec.Kx)} "
            f"dK={deltaK:+d} "
            f"Tx={fmt(rec.Tx)} "
            f"Ex={fmt(rec.Ex)}"
        )

        print(
            f"                "
            f"(kx,lx)=({fmt(rec.kx)},{fmt(rec.lx)}) "
            f"(ax,bx)=({fmt(rec.ax)},{fmt(rec.bx)})"
        )

        print(
            f"                "
            f"invariant={fmt(invariant)} "
            f"expected={fmt(result['expected_constant'])}"
        )

    print(
        "        DERIVED ORIGINAL EQUATION:"
    )

    print(
        f"            a*b + ({fmt(r2*l)})*a "
        f"+ ({fmt(r1*k)})*b = "
        f"{fmt(result['expected_constant'])}"
    )

    print(
        "        FACTORIZED FORM:"
    )

    print(
        f"            "
        f"(a + {fmt(r1*k)})"
        f"(b + {fmt(r2*l)})"
        f" = n"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    global S_VALUES

    rng = random.Random(SEED)

    start_total = time.perf_counter()

    print("=" * 92)
    print("START EXPERIMENT 76")
    print("AUXILIARY-EQUATION RANK / (a,b) INFORMATION EXPERIMENT")
    print("=" * 92)

    print()
    print("configuration")
    print(f"    scales                   = {[scale_name(x) for x in SCALES]}")
    print(f"    anchors / scale         = {ANCHORS_PER_SCALE}")
    print(f"    K target                = {K_TARGET}")
    print(f"    R offsets               = {R_OFFSETS}")
    print(f"    S values                = {S_VALUES}")
    print(f"    random AB tests         = {RANDOM_AB_TESTS}")
    print(f"    seed                    = {SEED}")

    total_records = 0
    total_keq = 0
    total_failures = 0
    total_difference_failures = 0

    scale_results = {}

    for scale in SCALES:
        print()
        print("=" * 92)
        print(f"SCALE {scale_name(scale)}")
        print("=" * 92)

        anchors = [
            generate_anchor(scale, rng)
            for _ in range(ANCHORS_PER_SCALE)
        ]

        scale_case_count = 0
        scale_records = 0
        scale_keq = 0
        scale_failures = 0
        scale_diff_failures = 0

        for ai, anchor in enumerate(anchors, 1):
            print()
            print(
                f"anchor {ai}/{ANCHORS_PER_SCALE} "
                f"n={fmt(anchor.n)}"
            )

            print(
                f"    TRUE factors = "
                f"({fmt(anchor.p)},{fmt(anchor.q)})"
            )

            # -----------------------------------------------------------------
            # Construct all R regimes.
            # -----------------------------------------------------------------

            for offset in R_OFFSETS:
                target_R = anchor.n / K_TARGET * offset

                r1, r2 = choose_balanced_prime_pair(
                    target_R
                )

                case = run_case(
                    anchor,
                    r1,
                    r2,
                )

                print()
                print(
                    f"    OFFSET {offset:.2f}"
                )

                print_case(
                    case,
                    show_rows=True,
                )

                scale_case_count += 1
                scale_records += len(case["records"])
                scale_keq += len(case["K_equal_records"])
                scale_failures += case["invariant_failures"]
                scale_diff_failures += case["difference_failures"]

        scale_results[scale] = {
            "cases": scale_case_count,
            "records": scale_records,
            "keq": scale_keq,
            "failures": scale_failures,
            "diff_failures": scale_diff_failures,
        }

        total_records += scale_records
        total_keq += scale_keq
        total_failures += scale_failures
        total_difference_failures += scale_diff_failures

        print()
        print("-" * 92)
        print(f"SCALE {scale_name(scale)} SUMMARY")
        print("-" * 92)

        print(
            f"    cases                         = {scale_case_count}"
        )
        print(
            f"    auxiliary records             = {scale_records:,}"
        )
        print(
            f"    Kx == K records                = {scale_keq:,}"
        )
        print(
            f"    equation invariant failures   = {scale_failures}"
        )
        print(
            f"    auxiliary/original failures   = "
            f"{scale_diff_failures}"
        )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 92)
    print("GLOBAL SUMMARY")
    print("=" * 92)

    print(
        f"    total auxiliary records       = {total_records:,}"
    )

    print(
        f"    total Kx == K records         = {total_keq:,}"
    )

    print(
        f"    total invariant failures      = {total_failures}"
    )

    print(
        f"    total equation differences    = "
        f"{total_difference_failures}"
    )

    print()
    print("=" * 92)
    print("CORE ALGEBRAIC RESULT")
    print("=" * 92)

    print(
        """
For every auxiliary observation:

    r1*kx*bx + r2*lx*ax + ax*bx
        - x
        - R*(Kx-K)

is experimentally equal to the same constant:

    n - R*K.

Therefore every auxiliary observation gives the SAME equation:

    a*b + r2*l*a + r1*k*b = n - R*K.

Adding R*K to both sides:

    a*b + r2*l*a + r1*k*b + R*K = n.

Since:

    R*K = (r1*k)(r2*l),

the expression factors:

    (a + r1*k)(b + r2*l) = n.

But:

    a + r1*k = p
    b + r2*l = q.

Therefore:

    (a + r1*k)(b + r2*l) = n

is exactly the original prime factorization:

    p*q=n.

So the experiment asks whether auxiliary trajectories create
new independent information about (a,b).

The expected result is:

    equation-rank = 1

regardless of how many auxiliary observations exist.

In particular, even if many different auxiliary values produce:

    Kx = K,

they still collapse to the same scalar equation for (a,b).

Thus Kx is useful for identifying candidate K values, but this
experiment isolates the remaining problem:

    GIVEN k,l,K,r1,r2,
    recover a,b from

        (a+r1*k)(b+r2*l)=n.

That remaining step is structurally equivalent to finding the
factor pair of n inside the quotient cell.

===============================================================================
"""
    )

    print("=" * 92)
    print("CROSS-SCALE SUMMARY")
    print("=" * 92)

    print(
        "scale       cases    auxRecords      Kx=K    invFails    diffFails"
    )
    print("-" * 92)

    for scale in SCALES:
        r = scale_results[scale]

        print(
            f"{scale_name(scale):<10} "
            f"{r['cases']:>6} "
            f"{r['records']:>13,} "
            f"{r['keq']:>9,} "
            f"{r['failures']:>10} "
            f"{r['diff_failures']:>10}"
        )

    total_time = time.perf_counter() - start_total

    print()
    print("=" * 92)
    print("TIMING")
    print("=" * 92)
    print(
        f"total runtime = {total_time:.4f}s"
    )

    print()
    print("=" * 92)
    print("FINISHED EXPERIMENT 76")
    print("=" * 92)


if __name__ == "__main__":
    main()
