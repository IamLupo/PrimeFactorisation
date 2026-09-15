#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 597
# ==============================================================================
# ACTIVE A/B RESIDUE RECURSION + FUNCTION COEFFICIENT RECURSION
#
# Goals:
#
#   1. Discover the active A/B residue at each 2-adic level.
#   2. Determine the residue recurrence.
#   3. Determine the A/B residue separation.
#   4. Show the exact coefficient recursion of the A and B frames.
#   5. Verify same-n coordinate halving.
#   6. Verify that the factor functions remain unchanged after
#      coordinate renormalization.
#
# No files.
# No external data.
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_EXAMPLES = True
EXAMPLE_COUNT = 8


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):

        if not sieve[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        sieve[start:limit + 1:p] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_semiprimes(limit, sieve):

    primes = [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > limit:
            break

        max_q = limit // p

        for q in primes[i:]:

            if q > max_q:
                break

            result.append(
                (p * q, p, q)
            )

    return result


# ==============================================================================
# FRAME A
# ==============================================================================

def frame_A(p: int, q: int, z: int):

    k = 1 << (z - 1)
    d = 2 * k

    x_num = q - p + 6
    y_num = p + q

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    x = x_num // d
    y = y_num // d

    return x, y


# ==============================================================================
# FRAME B
# ==============================================================================

def frame_B(p: int, q: int, z: int):

    k = 1 << (z - 1)
    d = 2 * k

    x_num = 3 * p - q + 6
    y_num = 3 * p + q

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    x = x_num // d
    y = y_num // d

    return x, y


# ==============================================================================
# ACTIVE RESIDUES
# ==============================================================================

def discover_active(semiprimes):

    active = {
        z: {
            "A": defaultdict(int),
            "B": defaultdict(int),
        }
        for z in range(
            MIN_Z,
            MAX_Z + 1,
        )
    }

    for n, p, q in semiprimes:

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            modulus = 1 << z
            residue = n % modulus

            if frame_A(p, q, z) is not None:
                active[z]["A"][residue] += 1

            if frame_B(p, q, z) is not None:
                active[z]["B"][residue] += 1

    return active


# ==============================================================================
# UNIQUE ACTIVE RESIDUE
# ==============================================================================

def unique_residue(active, z, frame):

    values = active[z][frame]

    if len(values) != 1:
        return None

    return next(iter(values))


# ==============================================================================
# ACTIVE RESIDUE REPORT
# ==============================================================================

def report_active(active):

    print()
    print("=" * 90)
    print("ACTIVE RESIDUES")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        modulus = 1 << z

        print()
        print(
            f"z={z} mod={modulus}"
        )

        for frame in ("A", "B"):

            entries = active[z][frame]

            if not entries:

                print(
                    f"    {frame}: NONE"
                )

                continue

            print(
                f"    {frame}: "
                f"{sorted(entries.items())}"
            )


# ==============================================================================
# RESIDUE RECURSION
# ==============================================================================

def report_residue_recursion(active):

    print()
    print("=" * 90)
    print("ACTIVE RESIDUE RECURSION")
    print("=" * 90)

    for frame in ("A", "B"):

        print()
        print(
            f"FRAME {frame}"
        )

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            residue = unique_residue(
                active,
                z,
                frame,
            )

            print(
                f"    z={z}: {residue}"
            )

        print()
        print(
            "    transitions:"
        )

        for z in range(
            MIN_Z,
            MAX_Z,
        ):

            r1 = unique_residue(
                active,
                z,
                frame,
            )

            r2 = unique_residue(
                active,
                z + 1,
                frame,
            )

            if (
                r1 is None
                or
                r2 is None
            ):
                continue

            delta = r2 - r1
            correction = r2 - 2 * r1

            print(
                f"        z={z}: "
                f"{r1} -> {r2} "
                f"delta={delta} "
                f"(r_next-2r={correction})"
            )


# ==============================================================================
# RECURRENCE TEST
# ==============================================================================

def test_recurrence(active, frame):

    print()
    print("=" * 90)
    print(
        f"RECURRENCE TEST FRAME {frame}"
    )
    print("=" * 90)

    pairs = []

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        r1 = unique_residue(
            active,
            z,
            frame,
        )

        r2 = unique_residue(
            active,
            z + 1,
            frame,
        )

        if (
            r1 is not None
            and
            r2 is not None
        ):
            pairs.append(
                (z, r1, r2)
            )

    if not pairs:

        print(
            "    insufficient data"
        )

        return

    candidates = {
        "2r-1": lambda r: 2 * r - 1,
        "2r+1": lambda r: 2 * r + 1,
        "2r-3": lambda r: 2 * r - 3,
        "2r+3": lambda r: 2 * r + 3,
        "r+6": lambda r: r + 6,
        "r-6": lambda r: r - 6,
    }

    for name, func in candidates.items():

        exact = all(
            func(r1) == r2
            for _, r1, r2 in pairs
        )

        print(
            f"    {name:<8}: "
            f"{'EXACT' if exact else 'no'}"
        )

    delta_values = sorted(
        {
            r2 - r1
            for _, r1, r2 in pairs
        }
    )

    doubled_values = sorted(
        {
            r2 - 2 * r1
            for _, r1, r2 in pairs
        }
    )

    print()
    print(
        f"    additive deltas = {delta_values}"
    )
    print(
        f"    doubled offsets = {doubled_values}"
    )


# ==============================================================================
# A/B RELATION
# ==============================================================================

def report_A_B_relation(active):

    print()
    print("=" * 90)
    print("A/B RESIDUE RELATION")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        A = unique_residue(
            active,
            z,
            "A",
        )

        B = unique_residue(
            active,
            z,
            "B",
        )

        print()
        print(
            f"z={z}"
        )

        print(
            f"    A = {A}"
        )

        print(
            f"    B = {B}"
        )

        if (
            A is not None
            and
            B is not None
        ):

            print(
                f"    B-A = {B-A}"
            )

            print(
                f"    A+B = {A+B}"
            )

            print(
                f"    B-2A = {B-2*A}"
            )

            print(
                f"    A-2B = {A-2*B}"
            )


# ==============================================================================
# CHILD BIT
# ==============================================================================

def report_child_bits(active):

    print()
    print("=" * 90)
    print("CHILD BIT OF ACTIVE RESIDUE")
    print("=" * 90)

    for frame in ("A", "B"):

        print()
        print(
            f"FRAME {frame}"
        )

        for z in range(
            MIN_Z,
            MAX_Z,
        ):

            r1 = unique_residue(
                active,
                z,
                frame,
            )

            r2 = unique_residue(
                active,
                z + 1,
                frame,
            )

            if (
                r1 is None
                or
                r2 is None
            ):
                continue

            step = 1 << z
            quotient = (r2 - r1)

            if quotient % step == 0:
                bit = quotient // step
            else:
                bit = "non-binary"

            print(
                f"    z={z}: "
                f"{r1} -> {r2} "
                f"child_bit={bit}"
            )


# ==============================================================================
# A-FRAME COEFFICIENTS
# ==============================================================================

def frame_A_coefficients(z):

    k = 1 << (z - 1)

    # p = -k*x + k*y + 3
    # q =  k*x + k*y - 3

    return (
        (-k, k, 3),
        (k, k, -3),
    )


# ==============================================================================
# B-FRAME COEFFICIENTS
# ==============================================================================

def frame_B_coefficients(z):

    k = Fraction(
        1 << (z - 1),
        1,
    )

    # 3p = -k*x + k*y + 3
    # q  =  k*x + k*y - 3
    #
    # Therefore:
    #
    # p = (-k/3)x + (k/3)y + 1
    # q =  kx + ky - 3

    return (
        (
            -k / 3,
            k / 3,
            Fraction(1),
        ),
        (
            k,
            k,
            Fraction(-3),
        ),
    )


# ==============================================================================
# FUNCTION COEFFICIENT RECURSION
# ==============================================================================

def report_function_recursion():

    print()
    print("=" * 90)
    print("FUNCTION COEFFICIENT RECURSION")
    print("=" * 90)

    # --------------------------------------------------------------------------
    # FRAME A
    # --------------------------------------------------------------------------

    print()
    print("FRAME A")

    previous = None

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        p_coeff, q_coeff = (
            frame_A_coefficients(z)
        )

        print()
        print(
            f"    z={z}"
        )

        print(
            f"        p="
            f"{p_coeff}"
        )

        print(
            f"        q="
            f"{q_coeff}"
        )

        if previous is not None:

            pp, pq = previous

            print(
                "        transition:"
            )

            print(
                f"            p_x: "
                f"{pp[0]} -> {p_coeff[0]}"
            )

            print(
                f"            p_y: "
                f"{pp[1]} -> {p_coeff[1]}"
            )

            print(
                f"            p_c: "
                f"{pp[2]} -> {p_coeff[2]}"
            )

            print(
                f"            q_x: "
                f"{pq[0]} -> {q_coeff[0]}"
            )

            print(
                f"            q_y: "
                f"{pq[1]} -> {q_coeff[1]}"
            )

            print(
                f"            q_c: "
                f"{pq[2]} -> {q_coeff[2]}"
            )

        previous = (
            p_coeff,
            q_coeff,
        )

    # --------------------------------------------------------------------------
    # FRAME B
    # --------------------------------------------------------------------------

    print()
    print("FRAME B")

    previous = None

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        p_coeff, q_coeff = (
            frame_B_coefficients(z)
        )

        print()
        print(
            f"    z={z}"
        )

        print(
            f"        p="
            f"{p_coeff}"
        )

        print(
            f"        q="
            f"{q_coeff}"
        )

        if previous is not None:

            pp, pq = previous

            print(
                "        transition:"
            )

            print(
                f"            p_x: "
                f"{pp[0]} -> {p_coeff[0]}"
            )

            print(
                f"            p_y: "
                f"{pp[1]} -> {p_coeff[1]}"
            )

            print(
                f"            p_c: "
                f"{pp[2]} -> {p_coeff[2]}"
            )

            print(
                f"            q_x: "
                f"{pq[0]} -> {q_coeff[0]}"
            )

            print(
                f"            q_y: "
                f"{pq[1]} -> {q_coeff[1]}"
            )

            print(
                f"            q_c: "
                f"{pq[2]} -> {q_coeff[2]}"
            )

        previous = (
            p_coeff,
            q_coeff,
        )


# ==============================================================================
# SYMBOLIC CONJUGACY
# ==============================================================================

def report_symbolic_conjugacy():

    print()
    print("=" * 90)
    print("SYMBOLIC LEVEL CONJUGACY")
    print("=" * 90)

    print(
        """
FRAME A

    p_z(x,y)
      = -k x + k y + 3

    q_z(x,y)
      =  k x + k y - 3

with:

    k = 2^(z-1).

At the next level:

    k' = 2k

and:

    x' = x/2
    y' = y/2.

Therefore:

    p_(z+1)(x/2,y/2)
      = -2k(x/2) + 2k(y/2) + 3
      = -kx + ky + 3
      = p_z(x,y).

Likewise:

    q_(z+1)(x/2,y/2)
      = q_z(x,y).

The same argument applies to the B frame.
"""
    )

    print()
    print("FRAME A:")
    print(
        "    p_(z+1)(x/2,y/2) = p_z(x,y)"
    )
    print(
        "    q_(z+1)(x/2,y/2) = q_z(x,y)"
    )

    print()
    print("FRAME B:")
    print(
        "    p_(z+1)(x/2,y/2) = p_z(x,y)"
    )
    print(
        "    q_(z+1)(x/2,y/2) = q_z(x,y)"
    )


# ==============================================================================
# SAME-N HALVING TEST
# ==============================================================================

def report_halving(semiprimes):

    print()
    print("=" * 90)
    print("SAME-N COORDINATE HALVING")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        A_count = 0
        A_fail = 0

        B_count = 0
        B_fail = 0

        for n, p, q in semiprimes:

            parent_A = frame_A(
                p,
                q,
                z,
            )

            child_A = frame_A(
                p,
                q,
                z + 1,
            )

            if (
                parent_A is not None
                and
                child_A is not None
            ):

                A_count += 1

                if (
                    child_A[0] * 2
                    != parent_A[0]
                    or
                    child_A[1] * 2
                    != parent_A[1]
                ):

                    A_fail += 1

            parent_B = frame_B(
                p,
                q,
                z,
            )

            child_B = frame_B(
                p,
                q,
                z + 1,
            )

            if (
                parent_B is not None
                and
                child_B is not None
            ):

                B_count += 1

                if (
                    child_B[0] * 2
                    != parent_B[0]
                    or
                    child_B[1] * 2
                    != parent_B[1]
                ):

                    B_fail += 1

        print()
        print(
            f"z={z}"
        )

        print(
            f"    A: "
            f"transitions={A_count} "
            f"failures={A_fail}"
        )

        print(
            f"    B: "
            f"transitions={B_count} "
            f"failures={B_fail}"
        )


# ==============================================================================
# SAMPLE TRANSITIONS
# ==============================================================================

def report_examples(semiprimes):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print("EXAMPLE SAME-N RENORMALIZATIONS")
    print("=" * 90)

    shown = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        count = 0

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        for n, p, q in semiprimes:

            parent = None
            child = None
            frame = None

            a0 = frame_A(
                p,
                q,
                z,
            )

            a1 = frame_A(
                p,
                q,
                z + 1,
            )

            if (
                a0 is not None
                and
                a1 is not None
            ):

                parent = a0
                child = a1
                frame = "A"

            else:

                b0 = frame_B(
                    p,
                    q,
                    z,
                )

                b1 = frame_B(
                    p,
                    q,
                    z + 1,
                )

                if (
                    b0 is not None
                    and
                    b1 is not None
                ):

                    parent = b0
                    child = b1
                    frame = "B"

            if parent is None:
                continue

            print(
                f"    {frame}: "
                f"n={n} "
                f"{parent} -> {child}"
            )

            count += 1

            if count >= EXAMPLE_COUNT:
                break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 597 START")
    print("=" * 90)

    print()
    print(
        "FAST ACTIVE A/B RESIDUE RECURSION"
    )

    # --------------------------------------------------------------------------
    # Sieve
    # --------------------------------------------------------------------------

    print()
    print("[1] PRIME SIEVE")

    sieve = prime_sieve(
        MAX_N
    )

    # --------------------------------------------------------------------------
    # Semiprimes
    # --------------------------------------------------------------------------

    print()
    print("[2] SEMIPRIME GENERATION")

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes="
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Active residues
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] ACTIVE RESIDUES"
    )

    active = discover_active(
        semiprimes
    )

    report_active(
        active
    )

    # --------------------------------------------------------------------------
    # Residue recursion
    # --------------------------------------------------------------------------

    report_residue_recursion(
        active
    )

    # --------------------------------------------------------------------------
    # Recurrence tests
    # --------------------------------------------------------------------------

    test_recurrence(
        active,
        "A",
    )

    test_recurrence(
        active,
        "B",
    )

    # --------------------------------------------------------------------------
    # A/B relation
    # --------------------------------------------------------------------------

    report_A_B_relation(
        active
    )

    # --------------------------------------------------------------------------
    # Child bits
    # --------------------------------------------------------------------------

    report_child_bits(
        active
    )

    # --------------------------------------------------------------------------
    # Function coefficients
    # --------------------------------------------------------------------------

    report_function_recursion()

    # --------------------------------------------------------------------------
    # Symbolic conjugacy
    # --------------------------------------------------------------------------

    report_symbolic_conjugacy()

    # --------------------------------------------------------------------------
    # Same-n halving
    # --------------------------------------------------------------------------

    report_halving(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Examples
    # --------------------------------------------------------------------------

    report_examples(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Final audit
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print("FINAL AUDIT")
    print("=" * 90)

    print(
        """
The active-residue structure is now separated into three layers.

1. RESIDUE RECURSION

       r_z
        ->
       r_(z+1)

2. COORDINATE RENORMALIZATION

       x_(z+1) = x_z / 2
       y_(z+1) = y_z / 2

3. FUNCTION CONJUGACY

       F_(z+1)(x/2,y/2)
           =
       F_z(x,y)

For FRAME A:

       p = -kx + ky + 3
       q =  kx + ky - 3

For FRAME B:

       p = -(k/3)x + (k/3)y + 1
       q =  kx + ky - 3

where:

       k = 2^(z-1).

The important question is therefore no longer whether the
coefficients grow with z -- they necessarily double because
k doubles.

The real question is whether the ACTIVE RESIDUE RECURRENCE
selects exactly which frame remains valid at the next level.

That is the next structural layer to attack.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 597 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()