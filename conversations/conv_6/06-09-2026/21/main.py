#!/usr/bin/env python3

"""
START EXPERIMENT 177B

Direct P-CRT attack using the CORRECT c3 cell targets.

Important distinction:

The benchmark generator knows the hidden p,q and therefore knows the
observed c3 value of every (r,s) cell.

The SEARCH itself does not know p or q.

For each cell:

    p = r*k + a
    q = s*l + b

    beta  = (k*b) mod s
    alpha = (l*a) mod r

    c3 = floor((r*beta + s*alpha + a*b)/(r*s))

For the true residues (a,b), c3 is the observed cell value.

From pq=n:

    beta  ≡ (n - ab) * r^{-1} (mod s)
    alpha ≡ (n - ab) * s^{-1} (mod r)

Therefore once c3 is known, the valid (a,b) pairs can be found
using ONLY n,r,s,c3.

For every valid pair:

    p ≡ a             (mod r)
    p ≡ n*b^{-1}      (mod s)

giving a direct CRT class for p.

The experiment then searches combinations of independent cells.

"""

import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54, 60, 66, 72, 78, 84]

START_CELL_COUNT = 20

MAX_STATES = 2_000_000
MAX_EXACT_TESTS = 2_000_000

# We only need a few cells before their combined modulus
# should become enormous.
CELL_DEPTHS = [1, 2, 3, 4, 5]


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
    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(b, a % b)

    return g, y1, x1 - (a // b) * y1


def inv_mod(a, m):
    g, x, _ = egcd(a, m)

    if g != 1:
        raise ValueError(
            f"{a} has no inverse mod {m}"
        )

    return x % m


def crt_pair(a1, m1, a2, m2):
    """
    Combine:

        x = a1 mod m1
        x = a2 mod m2

    with gcd(m1,m2)=1.
    """

    t = (
        (a2 - a1)
        * inv_mod(m1 % m2, m2)
    ) % m2

    x = a1 + m1 * t
    m = m1 * m2

    return x % m, m


# ============================================================
# TRUE CELL CARRY -- ONLY USED TO BUILD THE BENCHMARK
# ============================================================

def true_c3(p, q, r, s):
    """
    Compute the actual c3 for this cell from the hidden factors.

    This is benchmark-side information only.
    """

    a = p % r
    b = q % s

    k = p // r
    ell = q // s

    beta = (k * b) % s
    alpha = (ell * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    return c3


# ============================================================
# CELL RELATION
# ============================================================

def build_cell_relation(n, r, s, c3_target):
    """
    Construct all (a,b) satisfying:

        C(n,r,s,a,b) = c3_target

    without using p or q.
    """

    inv_r_s = inv_mod(r % s, s)
    inv_s_r = inv_mod(s % r, r)

    relation = []

    for a in range(1, r):
        for b in range(1, s):

            z = a * b

            beta = (
                (n - z)
                * inv_r_s
            ) % s

            alpha = (
                (n - z)
                * inv_s_r
            ) % r

            c3 = (
                r * beta
                + s * alpha
                + z
            ) // (r * s)

            if c3 == c3_target:
                relation.append((a, b))

    return relation


# ============================================================
# BUILD BENCHMARK CELLS
# ============================================================

def build_cells(n, p, q):

    cells = []

    for i, r in enumerate(R1):
        for j, s in enumerate(R2):

            target = true_c3(
                p,
                q,
                r,
                s,
            )

            relation = build_cell_relation(
                n,
                r,
                s,
                target,
            )

            true_a = p % r
            true_b = q % s

            true_present = (
                (true_a, true_b)
                in relation
            )

            cells.append({
                "i": i,
                "j": j,
                "r": r,
                "s": s,
                "c3": target,
                "relation": relation,
                "size": len(relation),
                "true_pair": (
                    true_a,
                    true_b,
                ),
                "true_present": true_present,
            })

    return cells


# ============================================================
# CONVERT CELL PAIRS -> P CRT CLASSES
# ============================================================

def cell_p_classes(n, cell):

    r = cell["r"]
    s = cell["s"]

    classes = []

    for a, b in cell["relation"]:

        # From pq=n:

        # p*b = n mod s
        #
        # p = n*b^-1 mod s

        p_mod_s = (
            n
            * inv_mod(b, s)
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
                a,
                b,
            )
        )

    return classes


# ============================================================
# DIRECT P SEARCH
# ============================================================

def search_cells(n, cells):

    sqrt_n = math.isqrt(n)

    states = {
        (0, 1)
    }

    generated = 0
    exact_tests = 0

    start = time.perf_counter()

    for depth, cell in enumerate(cells, 1):

        classes = cell_p_classes(
            n,
            cell,
        )

        new_states = {}

        for P, M in states:

            for qP, qM, _, _ in classes:

                P2, M2 = crt_pair(
                    P,
                    M,
                    qP,
                    qM,
                )

                # No positive solution <= sqrt(n)
                # if smallest positive representative
                # is already too large.

                smallest = P2

                if smallest == 0:
                    smallest = M2

                if smallest > sqrt_n:
                    continue

                new_states[
                    (P2, M2)
                ] = True

                generated += 1

                if len(new_states) > MAX_STATES:

                    return {
                        "found": None,
                        "states": len(new_states),
                        "generated": generated,
                        "exact_tests": exact_tests,
                        "time": (
                            time.perf_counter()
                            - start
                        ),
                        "aborted": "MAX_STATES",
                    }

        states = set(new_states)

        print(
            f"        depth={depth} "
            f"cell=({cell['r']},{cell['s']}) "
            f"c3={cell['c3']} "
            f"relation={cell['size']} "
            f"states={len(states)}"
        )

        # ----------------------------------------------------
        # Once M > sqrt(n), at most one p <= sqrt(n)
        # exists in a state.
        # ----------------------------------------------------

        remaining = set()

        for P, M in states:

            if M <= sqrt_n:
                remaining.add((P, M))
                continue

            p_candidate = P

            if p_candidate == 0:
                p_candidate = M

            if not (
                1 < p_candidate <= sqrt_n
            ):
                continue

            exact_tests += 1

            if exact_tests > MAX_EXACT_TESTS:

                return {
                    "found": None,
                    "states": len(states),
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": "MAX_EXACT_TESTS",
                }

            if n % p_candidate == 0:

                q_candidate = (
                    n // p_candidate
                )

                return {
                    "found": (
                        p_candidate,
                        q_candidate,
                    ),
                    "states": len(states),
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": None,
                }

        states = remaining

        if not states:
            return {
                "found": None,
                "states": 0,
                "generated": generated,
                "exact_tests": exact_tests,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "aborted": None,
            }

    # --------------------------------------------------------
    # Final enumeration.
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
                P + t * M
            )

            if candidate < 2:
                continue

            exact_tests += 1

            if exact_tests > MAX_EXACT_TESTS:

                return {
                    "found": None,
                    "states": len(states),
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": "MAX_EXACT_TESTS",
                }

            if n % candidate == 0:

                q_candidate = (
                    n // candidate
                )

                return {
                    "found": (
                        candidate,
                        q_candidate,
                    ),
                    "states": len(states),
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted": None,
                }

    return {
        "found": None,
        "states": len(states),
        "generated": generated,
        "exact_tests": exact_tests,
        "time": (
            time.perf_counter()
            - start
        ),
        "aborted": None,
    }


# ============================================================
# CHOOSE INDEPENDENT R1/R2 CELLS
# ============================================================

def make_cell_sets(cells):

    cells = sorted(
        cells,
        key=lambda c: (
            c["size"],
            c["r"] * c["s"],
        ),
    )

    result = []

    def rec(chosen, start):

        if len(chosen) >= 1:

            result.append(
                list(chosen)
            )

        if len(chosen) == 5:
            return

        for idx in range(
            start,
            len(cells),
        ):

            c = cells[idx]

            if any(
                c["i"] == x["i"]
                or c["j"] == x["j"]
                for x in chosen
            ):
                continue

            rec(
                chosen + [c],
                idx + 1,
            )

            if len(result) >= 100:
                return

    rec([], 0)

    result.sort(
        key=lambda group: (
            math.prod(
                c["size"]
                for c in group
            ),
            -math.prod(
                c["r"] * c["s"]
                for c in group
            ),
        )
    )

    return result[:100]


# ============================================================
# TEST INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

    # --------------------------------------------------------
    # Generate balanced semiprime.
    # --------------------------------------------------------

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

    p = None
    q = None

    for x in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(x):
            continue

        target = (
            (1 << bits)
            // x
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            y = target + delta

            if y <= x:
                continue

            if is_prime(y):

                p = x
                q = y
                break

        if p is not None:
            break

    n = p * q

    sqrt_n = math.isqrt(n)

    print(f"bits(n) = {n.bit_length()}")
    print(f"sqrt(n) = {sqrt_n}")

    # Hidden only for benchmark verification.
    print(f"hidden p = {p}")
    print(f"hidden q = {q}")

    # --------------------------------------------------------
    # Build the actual C-cell observations.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cells = build_cells(
        n,
        p,
        q,
    )

    build_time = (
        time.perf_counter()
        - t0
    )

    print(
        f"C-cell build time = "
        f"{build_time:.6f}s"
    )

    # --------------------------------------------------------
    # Verify that every true pair survived.
    # --------------------------------------------------------

    failed = [
        c
        for c in cells
        if not c["true_present"]
    ]

    if failed:

        print(
            "ERROR: true pair missing "
            f"from {len(failed)} cells"
        )

        for c in failed:
            print(
                c["r"],
                c["s"],
                c["c3"],
                c["true_pair"],
            )

        print(
            f"FINISHED INSTANCE {bits}-BIT"
        )

        return

    # --------------------------------------------------------
    # Rarest cells.
    # --------------------------------------------------------

    rare = sorted(
        cells,
        key=lambda c: c["size"],
    )

    print()
    print("Rarest cells:")

    for c in rare[
        :START_CELL_COUNT
    ]:

        print(
            f"  ({c['r']:2d},{c['s']:2d}) "
            f"c3={c['c3']} "
            f"relation={c['size']:4d} "
            f"true={c['true_pair']}"
        )

    # --------------------------------------------------------
    # Independent cell sets.
    # --------------------------------------------------------

    sets = make_cell_sets(
        cells
    )

    print()
    print(
        f"Candidate cell groups = "
        f"{len(sets)}"
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    for set_no, group in enumerate(
        sets,
        1,
    ):

        if len(group) > 5:
            continue

        print()
        print(
            f"GROUP {set_no}: "
            f"depth={len(group)}"
        )

        modulus = 1

        for c in group:

            modulus *= (
                c["r"]
                * c["s"]
            )

            print(
                f"    ({c['r']},{c['s']}) "
                f"c3={c['c3']} "
                f"relation={c['size']}"
            )

        print(
            f"    combined modulus = "
            f"{modulus}"
        )

        result = search_cells(
            n,
            group,
        )

        print(
            f"    generated   = "
            f"{result['generated']}"
        )

        print(
            f"    states      = "
            f"{result['states']}"
        )

        print(
            f"    exact tests = "
            f"{result['exact_tests']}"
        )

        print(
            f"    time        = "
            f"{result['time']:.6f}s"
        )

        if result["aborted"]:
            print(
                f"    aborted     = "
                f"{result['aborted']}"
            )

        if result["found"]:

            fp, fq = result["found"]

            if fp > fq:
                fp, fq = fq, fp

            print()
            print("*** FACTOR FOUND ***")
            print(f"p = {fp}")
            print(f"q = {fq}")
            print(
                f"correct = "
                f"{fp * fq == n}"
            )

            print()
            print(
                f"FINISHED INSTANCE {bits}-BIT"
            )

            return

    print()
    print(
        "No factor recovered."
    )

    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print("START EXPERIMENT 177B")
    print()

    print(f"R1 = {R1}")
    print(f"R2 = {R2}")
    print(f"BITS = {BITS}")
    print(f"MAX_STATES = {MAX_STATES}")
    print(f"MAX_EXACT_TESTS = {MAX_EXACT_TESTS}")

    for bits in BITS:
        run_instance(bits)

    print()
    print("FINISHED EXPERIMENT 177B")


if __name__ == "__main__":
    main()
