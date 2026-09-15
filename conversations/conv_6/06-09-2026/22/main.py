#!/usr/bin/env python3

"""
START EXPERIMENT 178

C-VALUE RECOVERY WITHOUT FACTORS

Goal
----

Experiment 177B showed that a small number of observed C-cell values can
make direct p-CRT factor recovery very efficient.

The missing question is:

    Can useful C-cell information be obtained from n alone?

This experiment computes all possible C-values for every (r,s) cell using
ONLY n and the radices.

For a cell (r,s):

    a = p mod r
    b = q mod s
    z = a*b

and

    beta  = (n-z) * r^{-1} mod s
    alpha = (n-z) * s^{-1} mod r

so

    C =
        floor((r*beta + s*alpha + z)/(r*s))

depends only on n,r,s,z.

No p or q are required for the C-value enumeration.

The hidden p,q are used ONLY for benchmark verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

SHOW_CELLS = 25

MAX_C_VALUES = 16

MAX_STATES = 2_000_000

MAX_EXACT_TESTS = 2_000_000

MAX_C_COMBINATIONS = 10_000

MAX_GROUPS = 100


# ============================================================
# NUMBER THEORY
# ============================================================

def is_prime(n):
    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    ]

    for p in small:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:
        if n % d == 0:
            return False

        d += step
        step = 6 - step

    return True


def egcd(a, b):
    """
    Extended Euclidean algorithm.

    Returns g,x,y satisfying:

        a*x + b*y = g
    """

    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(
        b,
        a % b,
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def inv_mod(a, m):
    g, x, _ = egcd(a, m)

    if g != 1:
        raise ValueError(
            f"{a} is not invertible mod {m}; gcd={g}"
        )

    return x % m


def crt_pair(a1, m1, a2, m2):
    """
    Combine:

        x = a1 mod m1
        x = a2 mod m2

    assuming gcd(m1,m2)=1.
    """

    t = (
        (a2 - a1)
        * inv_mod(m1 % m2, m2)
    ) % m2

    x = a1 + m1 * t

    m = m1 * m2

    return x % m, m


# ============================================================
# TEST SEMIPRIME GENERATION
# ============================================================

def make_semiprime(bits):

    low = 1 << (bits // 2 - 1)

    high = 1 << (bits // 2 + 1)

    for p in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(p):
            continue

        target = (
            (1 << bits)
            // p
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            q = target + delta

            if q <= p:
                continue

            if not is_prime(q):
                continue

            n = p * q

            if n.bit_length() == bits:
                return p, q, n

    raise RuntimeError(
        f"Could not generate {bits}-bit semiprime"
    )


# ============================================================
# C VALUE FROM n ONLY
# ============================================================

def c_value_from_z(
    n,
    r,
    s,
    z,
    inv_r_mod_s,
    inv_s_mod_r,
):
    """
    Compute the C value from:

        n,r,s,z

    where z = a*b.
    """

    beta = (
        (n - z)
        * inv_r_mod_s
    ) % s

    alpha = (
        (n - z)
        * inv_s_mod_r
    ) % r

    numerator = (
        r * beta
        + s * alpha
        + z
    )

    denominator = r * s

    return numerator // denominator


# ============================================================
# POSSIBLE C VALUES
# ============================================================

def possible_c_values(
    n,
    r,
    s,
):

    inv_r_mod_s = inv_mod(
        r % s,
        s,
    )

    inv_s_mod_r = inv_mod(
        s % r,
        r,
    )

    values = set()

    for a in range(1, r):

        for b in range(1, s):

            z = a * b

            c = c_value_from_z(
                n,
                r,
                s,
                z,
                inv_r_mod_s,
                inv_s_mod_r,
            )

            values.add(c)

    return sorted(values)


# ============================================================
# TRUE C
# BENCHMARK ONLY
# ============================================================

def true_c3(
    p,
    q,
    r,
    s,
):

    a = p % r
    b = q % s

    k = p // r
    ell = q // s

    beta = (
        k * b
    ) % s

    alpha = (
        ell * a
    ) % r

    numerator = (
        r * beta
        + s * alpha
        + a * b
    )

    denominator = r * s

    return numerator // denominator


# ============================================================
# BUILD RELATION FOR SPECIFIC C
# ============================================================

def build_relation_for_c(
    n,
    r,
    s,
    c_target,
):

    inv_r_mod_s = inv_mod(
        r % s,
        s,
    )

    inv_s_mod_r = inv_mod(
        s % r,
        r,
    )

    relation = []

    for a in range(1, r):

        for b in range(1, s):

            z = a * b

            c = c_value_from_z(
                n,
                r,
                s,
                z,
                inv_r_mod_s,
                inv_s_mod_r,
            )

            if c == c_target:
                relation.append(
                    (a, b)
                )

    return relation


# ============================================================
# CELL RELATION -> P CRT CLASSES
# ============================================================

def cell_to_p_classes(
    n,
    r,
    s,
    relation,
):

    inv_b = {}

    classes = []

    for a, b in relation:

        if b not in inv_b:

            inv_b[b] = inv_mod(
                b,
                s,
            )

        p_mod_s = (
            n
            * inv_b[b]
        ) % s

        p_mod_rs, modulus = crt_pair(
            a,
            r,
            p_mod_s,
            s,
        )

        classes.append(
            (
                p_mod_rs,
                modulus,
            )
        )

    return classes


# ============================================================
# DIRECT P SEARCH
# ============================================================

def search_cells(
    n,
    cell_specs,
):

    sqrt_n = math.isqrt(n)

    states = {
        (0, 1)
    }

    generated = 0

    exact_tests = 0

    start = time.perf_counter()

    for spec in cell_specs:

        r = spec["r"]
        s = spec["s"]
        c = spec["c"]

        relation = build_relation_for_c(
            n,
            r,
            s,
            c,
        )

        classes = cell_to_p_classes(
            n,
            r,
            s,
            relation,
        )

        next_states = {}

        for P, M in states:

            for p2, m2 in classes:

                P3, M3 = crt_pair(
                    P,
                    M,
                    p2,
                    m2,
                )

                smallest = P3

                if smallest == 0:
                    smallest = M3

                if smallest > sqrt_n:
                    continue

                next_states[
                    (P3, M3)
                ] = True

                generated += 1

                if len(next_states) > MAX_STATES:

                    return {
                        "found": None,
                        "generated": generated,
                        "states": len(next_states),
                        "exact_tests": exact_tests,
                        "time": (
                            time.perf_counter()
                            - start
                        ),
                        "aborted": "MAX_STATES",
                    }

        states = set(next_states)

        if not states:

            return {
                "found": None,
                "generated": generated,
                "states": 0,
                "exact_tests": exact_tests,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "aborted": None,
            }

        remaining = set()

        for P, M in states:

            if M <= sqrt_n:

                remaining.add(
                    (P, M)
                )

                continue

            candidate = P

            if candidate == 0:
                candidate = M

            if not (
                1 < candidate <= sqrt_n
            ):
                continue

            exact_tests += 1

            if exact_tests > MAX_EXACT_TESTS:

                return {
                    "found": None,
                    "generated": generated,
                    "states": len(states),
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": "MAX_EXACT_TESTS",
                }

            if n % candidate == 0:

                q = n // candidate

                return {
                    "found": (
                        candidate,
                        q,
                    ),
                    "generated": generated,
                    "states": len(states),
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": None,
                }

        states = remaining

    # --------------------------------------------------------
    # Final enumeration
    # --------------------------------------------------------

    for P, M in states:

        if P == 0:
            P = M

        if P > sqrt_n:
            continue

        count = (
            (sqrt_n - P)
            // M
        )

        for t in range(
            count + 1
        ):

            candidate = (
                P
                + t * M
            )

            if candidate < 2:
                continue

            exact_tests += 1

            if exact_tests > MAX_EXACT_TESTS:

                return {
                    "found": None,
                    "generated": generated,
                    "states": len(states),
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": "MAX_EXACT_TESTS",
                }

            if n % candidate == 0:

                q = n // candidate

                return {
                    "found": (
                        candidate,
                        q,
                    ),
                    "generated": generated,
                    "states": len(states),
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": None,
                }

    return {
        "found": None,
        "generated": generated,
        "states": len(states),
        "exact_tests": exact_tests,
        "time": (
            time.perf_counter()
            - start
        ),
        "aborted": None,
    }


# ============================================================
# BUILD CELL INFORMATION
# ============================================================

def build_cells(n):

    cells = []

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            values = possible_c_values(
                n,
                r,
                s,
            )

            cells.append({
                "i": i,
                "j": j,
                "r": r,
                "s": s,
                "values": values,
                "entropy": len(values),
            })

    cells.sort(
        key=lambda c: (
            c["entropy"],
            c["r"] * c["s"],
        )
    )

    return cells


# ============================================================
# BUILD CANDIDATE CELL GROUPS
# ============================================================

def build_groups(cells):

    useful = [
        c
        for c in cells
        if c["entropy"] <= MAX_C_VALUES
    ]

    groups = []

    # --------------------------------------------------------
    # Single cells
    # --------------------------------------------------------

    for c in useful:

        groups.append([
            c
        ])

        if len(groups) >= MAX_GROUPS:
            return groups

    # --------------------------------------------------------
    # Independent pairs
    # --------------------------------------------------------

    limit = min(
        len(useful),
        SHOW_CELLS,
    )

    for i in range(limit):

        c1 = useful[i]

        for j in range(
            i + 1,
            limit,
        ):

            c2 = useful[j]

            if c1["i"] == c2["i"]:
                continue

            if c1["j"] == c2["j"]:
                continue

            groups.append([
                c1,
                c2,
            ])

            if len(groups) >= MAX_GROUPS:
                return groups

    return groups


# ============================================================
# RUN INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(
        f"START INSTANCE {bits}-BIT"
    )
    print("=" * 72)

    p, q, n = make_semiprime(
        bits
    )

    sqrt_n = math.isqrt(n)

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {sqrt_n}"
    )

    # Hidden benchmark values.
    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    # --------------------------------------------------------
    # C-value enumeration from n alone
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cells = build_cells(n)

    c_time = (
        time.perf_counter()
        - t0
    )

    print()
    print(
        f"C-value enumeration = "
        f"{c_time:.6f}s"
    )

    # --------------------------------------------------------
    # Verify true C values are contained
    # --------------------------------------------------------

    true_missing = []

    for c in cells:

        actual_c = true_c3(
            p,
            q,
            c["r"],
            c["s"],
        )

        c["true_c"] = actual_c

        if actual_c not in c["values"]:

            true_missing.append(
                (
                    c["r"],
                    c["s"],
                    actual_c,
                )
            )

    if true_missing:

        print()
        print(
            "ERROR: true C value was not "
            "present in candidate set."
        )

        for item in true_missing:

            print(
                f"    cell=({item[0]},{item[1]}) "
                f"true_c={item[2]}"
            )

        print()
        print(
            f"FINISHED INSTANCE {bits}-BIT"
        )

        return

    # --------------------------------------------------------
    # Print entropy table
    # --------------------------------------------------------

    print()
    print(
        "Cells ordered by number "
        "of possible C-values:"
    )

    for c in cells[:SHOW_CELLS]:

        print(
            f"  ({c['r']:2d},{c['s']:2d}) "
            f"C-count={c['entropy']:2d} "
            f"C-values={c['values']} "
            f"true={c['true_c']}"
        )

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    unique_count = sum(
        1
        for c in cells
        if c["entropy"] == 1
    )

    two_to_four_count = sum(
        1
        for c in cells
        if 2 <= c["entropy"] <= 4
    )

    five_to_eight_count = sum(
        1
        for c in cells
        if 5 <= c["entropy"] <= 8
    )

    max_value_count = sum(
        1
        for c in cells
        if c["entropy"] <= MAX_C_VALUES
    )

    print()
    print(
        f"UNIQUE cells = "
        f"{unique_count}/25"
    )

    print(
        f"2..4-value cells = "
        f"{two_to_four_count}/25"
    )

    print(
        f"5..8-value cells = "
        f"{five_to_eight_count}/25"
    )

    print(
        f"<= {MAX_C_VALUES}-value cells = "
        f"{max_value_count}/25"
    )

    # --------------------------------------------------------
    # Build candidate groups
    # --------------------------------------------------------

    groups = build_groups(
        cells
    )

    print()
    print(
        f"Candidate cell groups = "
        f"{len(groups)}"
    )

    # --------------------------------------------------------
    # Factor-blind search
    # --------------------------------------------------------

    for group_no, group in enumerate(
        groups,
        1,
    ):

        c_combination_count = math.prod(
            c["entropy"]
            for c in group
        )

        print()
        print(
            f"GROUP {group_no}: "
            f"cells="
            + ", ".join(
                f"({c['r']},{c['s']})"
                for c in group
            )
        )

        print(
            f"    C combinations = "
            f"{c_combination_count}"
        )

        if (
            c_combination_count
            > MAX_C_COMBINATIONS
        ):

            print(
                "    SKIP: too many "
                "C assignments"
            )

            continue

        combinations = [
            []
        ]

        for c in group:

            expanded = []

            for partial in combinations:

                for value in c["values"]:

                    expanded.append(
                        partial
                        + [value]
                    )

            combinations = expanded

        # ----------------------------------------------------
        # Test every candidate C assignment.
        # ----------------------------------------------------

        for combo in combinations:

            specs = []

            for c, cv in zip(
                group,
                combo,
            ):

                specs.append({
                    "r": c["r"],
                    "s": c["s"],
                    "c": cv,
                })

            result = search_cells(
                n,
                specs,
            )

            if result["found"]:

                fp, fq = result["found"]

                if fp > fq:
                    fp, fq = fq, fp

                print()
                print(
                    "*** FACTOR FOUND ***"
                )

                print(
                    f"    C assignment = "
                    f"{combo}"
                )

                print(
                    f"    p = {fp}"
                )

                print(
                    f"    q = {fq}"
                )

                print(
                    f"    correct = "
                    f"{fp * fq == n}"
                )

                print(
                    f"    generated = "
                    f"{result['generated']}"
                )

                print(
                    f"    exact tests = "
                    f"{result['exact_tests']}"
                )

                print(
                    f"    search time = "
                    f"{result['time']:.6f}s"
                )

                print()
                print(
                    f"FINISHED INSTANCE "
                    f"{bits}-BIT"
                )

                return

            if result["aborted"]:
                break

    print()
    print(
        "No factor recovered from "
        "factor-blind C candidates."
    )

    print()
    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "START EXPERIMENT 178"
    )

    print()

    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    print(
        f"BITS = {BITS}"
    )

    print(
        f"MAX_C_VALUES = "
        f"{MAX_C_VALUES}"
    )

    print(
        f"MAX_STATES = "
        f"{MAX_STATES}"
    )

    print(
        f"MAX_EXACT_TESTS = "
        f"{MAX_EXACT_TESTS}"
    )

    print()

    for bits in BITS:
        run_instance(bits)

    print()
    print(
        "FINISHED EXPERIMENT 178"
    )


if __name__ == "__main__":
    main()