#!/usr/bin/env python3

"""
START EXPERIMENT 179B

GLOBAL C-MATRIX CONSISTENCY
WITH GENERALIZED CRT

Experiment 179 established the idea of solving globally over the
(a,b) residue structure, but its final CRT stage incorrectly assumed
that all cell moduli were pairwise coprime.

That is false.

For example:

    17*19 = 323
    17*37 = 629

and

    gcd(323,629) = 17.

Therefore ordinary CRT cannot be used.

This experiment fixes that by using generalized CRT:

    x ≡ a (mod m)
    x ≡ b (mod n)

is solvable iff

    a ≡ b (mod gcd(m,n)).

When compatible, the merged modulus is lcm(m,n).

The CSP itself remains factor-blind.

The hidden p,q are used only for verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

MAX_SOLUTIONS = 100_000

MAX_NODE_COUNT = 2_000_000

MAX_EXACT_TESTS = 2_000_000

# Stop after this many seconds in the global CSP phase.
MAX_GLOBAL_TIME = 30.0


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

    g, x, _ = egcd(
        a,
        m,
    )

    if g != 1:
        raise ValueError(
            f"{a} not invertible mod {m}; gcd={g}"
        )

    return x % m


# ============================================================
# GENERALIZED CRT
# ============================================================

def crt_pair_general(
    a1,
    m1,
    a2,
    m2,
):
    """
    Generalized CRT.

    Solve:

        x ≡ a1 (mod m1)
        x ≡ a2 (mod m2)

    Returns:

        (x, lcm(m1,m2))

    or:

        None

    if the congruences are incompatible.
    """

    g = math.gcd(
        m1,
        m2,
    )

    # Necessary compatibility condition.
    if (a2 - a1) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    # m1r and m2r are now coprime.
    #
    # a1 + m1*t = a2 (mod m2)
    #
    # m1r*t = (a2-a1)/g (mod m2r)

    rhs = (
        (a2 - a1)
        // g
    ) % m2r

    if m2r == 1:
        t = 0
    else:
        inv = inv_mod(
            m1r % m2r,
            m2r,
        )

        t = (
            rhs * inv
        ) % m2r

    x = (
        a1
        + m1 * t
    )

    lcm = (
        m1
        * m2r
    )

    return (
        x % lcm,
        lcm,
    )


# ============================================================
# TEST SEMIPRIME
# ============================================================

def make_semiprime(bits):

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

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
        f"Unable to generate {bits}-bit semiprime"
    )


# ============================================================
# C-VALUE
# ============================================================

def c_value(
    n,
    r,
    s,
    z,
    inv_r_s,
    inv_s_r,
):

    beta = (
        (n - z)
        * inv_r_s
    ) % s

    alpha = (
        (n - z)
        * inv_s_r
    ) % r

    return (
        r * beta
        + s * alpha
        + z
    ) // (
        r * s
    )


# ============================================================
# CELL CANDIDATES
# ============================================================

def build_cell_candidates(
    n,
    r,
    s,
):

    inv_r_s = inv_mod(
        r % s,
        s,
    )

    inv_s_r = inv_mod(
        s % r,
        r,
    )

    candidates = []

    for a in range(
        1,
        r,
    ):

        for b in range(
            1,
            s,
        ):

            z = a * b

            c = c_value(
                n,
                r,
                s,
                z,
                inv_r_s,
                inv_s_r,
            )

            candidates.append({
                "a": a,
                "b": b,
                "z": z,
                "c": c,
            })

    return candidates


# ============================================================
# ALL CELLS
# ============================================================

def build_all_cells(n):

    cells = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            candidates = (
                build_cell_candidates(
                    n,
                    r,
                    s,
                )
            )

            cells[(i, j)] = {
                "i": i,
                "j": j,
                "r": r,
                "s": s,
                "candidates": candidates,
            }

    return cells


# ============================================================
# TRUE C
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

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (
        r * s
    )


# ============================================================
# INITIAL A/B DOMAINS
# ============================================================

def project_domains(
    cells,
):

    A = {}

    B = {}

    for i, r in enumerate(R1):

        domain = set(
            range(
                1,
                r,
            )
        )

        for j in range(
            len(R2)
        ):

            cell = cells[
                (i, j)
            ]

            possible = {
                x["a"]
                for x in cell["candidates"]
            }

            domain &= possible

        A[i] = domain

    for j, s in enumerate(R2):

        domain = set(
            range(
                1,
                s,
            )
        )

        for i in range(
            len(R1)
        ):

            cell = cells[
                (i, j)
            ]

            possible = {
                x["b"]
                for x in cell["candidates"]
            }

            domain &= possible

        B[j] = domain

    return A, B


# ============================================================
# CELL FILTER
# ============================================================

def valid_candidates(
    cell,
    A,
    B,
):

    i = cell["i"]
    j = cell["j"]

    return [
        x
        for x in cell["candidates"]
        if x["a"] in A[i]
        and x["b"] in B[j]
    ]


# ============================================================
# ARC CONSISTENCY
# ============================================================

def arc_consistency(
    cells,
    A,
    B,
):

    changed = True

    while changed:

        changed = False

        for (
            i,
            j,
        ), cell in cells.items():

            valid = valid_candidates(
                cell,
                A,
                B,
            )

            if not valid:
                return False

            new_a = {
                x["a"]
                for x in valid
            }

            new_b = {
                x["b"]
                for x in valid
            }

            old_a = A[i]
            old_b = B[j]

            reduced_a = (
                old_a
                & new_a
            )

            reduced_b = (
                old_b
                & new_b
            )

            if reduced_a != old_a:

                A[i] = reduced_a
                changed = True

            if reduced_b != old_b:

                B[j] = reduced_b
                changed = True

            if not A[i] or not B[j]:
                return False

    return True


# ============================================================
# VERIFY ASSIGNMENT
# ============================================================

def assignment_to_residues(
    assignment,
):

    A = {}
    B = {}

    for (
        i,
        j,
    ), pair in assignment.items():

        a, b = pair

        if i in A and A[i] != a:
            return None

        if j in B and B[j] != b:
            return None

        A[i] = a
        B[j] = b

    return A, B


# ============================================================
# GENERALIZED CRT FACTOR SEARCH
# ============================================================

def try_factor_from_assignment(
    n,
    assignment,
):

    sqrt_n = math.isqrt(n)

    # --------------------------------------------------------
    # Convert every cell into a p congruence:
    #
    #     p ≡ a mod r
    #
    # and:
    #
    #     p ≡ n*b^-1 mod s
    #
    # giving:
    #
    #     p ≡ P mod lcm(r,s)=r*s
    #
    # for that individual cell.
    # --------------------------------------------------------

    constraints = []

    for (
        i,
        j,
    ), pair in assignment.items():

        r = R1[i]
        s = R2[j]

        a, b = pair

        p_mod_s = (
            n
            * inv_mod(
                b,
                s,
            )
        ) % s

        # r and s are distinct primes,
        # so THIS individual CRT is ordinary CRT.
        individual = crt_pair_general(
            a,
            r,
            p_mod_s,
            s,
        )

        if individual is None:

            return None, {
                "generated": 0,
                "exact_tests": 0,
                "time": 0.0,
                "aborted": "CELL_CONTRADICTION",
            }

        p_class, p_modulus = individual

        constraints.append(
            (
                p_class,
                p_modulus,
                i,
                j,
            )
        )

    # --------------------------------------------------------
    # Apply strongest constraints first.
    #
    # A small modulus is weaker, but constraints with unusual
    # overlap can cause early contradictions. Sorting by modulus
    # descending tends to expose those sooner.
    # --------------------------------------------------------

    constraints.sort(
        key=lambda x:
        -x[1]
    )

    states = {
        (0, 1)
    }

    generated = 0

    exact_tests = 0

    start = time.perf_counter()

    for (
        p_class,
        p_modulus,
        _,
        _,
    ) in constraints:

        next_states = {}

        for P, M in states:

            merged = crt_pair_general(
                P,
                M,
                p_class,
                p_modulus,
            )

            # Shared radices may cause a gcd > 1.
            # Incompatible states are simply discarded.
            if merged is None:
                continue

            P2, M2 = merged

            smallest = P2

            if smallest == 0:
                smallest = M2

            if smallest > sqrt_n:
                continue

            next_states[
                (P2, M2)
            ] = True

            generated += 1

            if generated > MAX_NODE_COUNT:

                return None, {
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted":
                        "MAX_NODE_COUNT",
                }

        states = set(
            next_states
        )

        if not states:

            return None, {
                "generated": generated,
                "exact_tests": exact_tests,
                "time": (
                    time.perf_counter()
                    - start
                ),
                "aborted": None,
            }

        # ----------------------------------------------------
        # Exact test when modulus > sqrt(n).
        # ----------------------------------------------------

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

            if (
                exact_tests
                > MAX_EXACT_TESTS
            ):

                return None, {
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted":
                        "MAX_EXACT_TESTS",
                }

            if n % candidate == 0:

                q = n // candidate

                return (
                    (
                        candidate,
                        q,
                    ),
                    {
                        "generated": generated,
                        "exact_tests": exact_tests,
                        "time": (
                            time.perf_counter()
                            - start
                        ),
                        "aborted": None,
                    },
                )

        states = remaining

    # --------------------------------------------------------
    # Final enumeration.
    # --------------------------------------------------------

    for P, M in states:

        if P == 0:
            P = M

        if P > sqrt_n:
            continue

        count = (
            sqrt_n - P
        ) // M

        for t in range(
            count + 1
        ):

            candidate = (
                P + t * M
            )

            if candidate < 2:
                continue

            exact_tests += 1

            if (
                exact_tests
                > MAX_EXACT_TESTS
            ):

                return None, {
                    "generated": generated,
                    "exact_tests": exact_tests,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                    "aborted":
                        "MAX_EXACT_TESTS",
                }

            if n % candidate == 0:

                q = n // candidate

                return (
                    (
                        candidate,
                        q,
                    ),
                    {
                        "generated": generated,
                        "exact_tests": exact_tests,
                        "time": (
                            time.perf_counter()
                            - start
                        ),
                        "aborted": None,
                    },
                )

    return None, {
        "generated": generated,
        "exact_tests": exact_tests,
        "time": (
            time.perf_counter()
            - start
        ),
        "aborted": None,
    }


# ============================================================
# GLOBAL CSP
# ============================================================

def solve_global_csp(
    n,
    cells,
    A,
    B,
):

    nodes = 0

    solutions = 0

    started = time.perf_counter()

    def recurse(
        A_local,
        B_local,
        assignment,
    ):

        nonlocal nodes
        nonlocal solutions

        if (
            time.perf_counter()
            - started
            > MAX_GLOBAL_TIME
        ):
            return

        nodes += 1

        if nodes > MAX_NODE_COUNT:
            return

        if solutions >= MAX_SOLUTIONS:
            return

        # ----------------------------------------------------
        # Propagate domains.
        # ----------------------------------------------------

        if not arc_consistency(
            cells,
            A_local,
            B_local,
        ):
            return

        # ----------------------------------------------------
        # Select unassigned cell with the fewest candidates.
        # ----------------------------------------------------

        best_key = None

        best_candidates = None

        for key, cell in cells.items():

            if key in assignment:
                continue

            vc = valid_candidates(
                cell,
                A_local,
                B_local,
            )

            if not vc:
                return

            if (
                best_candidates is None
                or len(vc)
                < len(best_candidates)
            ):

                best_key = key
                best_candidates = vc

        # ----------------------------------------------------
        # Full assignment.
        # ----------------------------------------------------

        if best_key is None:

            solutions += 1

            yield dict(
                assignment
            )

            return

        i, j = best_key

        # ----------------------------------------------------
        # Try candidate pair.
        # ----------------------------------------------------

        for candidate in best_candidates:

            if (
                time.perf_counter()
                - started
                > MAX_GLOBAL_TIME
            ):
                return

            A2 = {
                k: set(v)
                for k, v in A_local.items()
            }

            B2 = {
                k: set(v)
                for k, v in B_local.items()
            }

            assignment2 = dict(
                assignment
            )

            assignment2[
                best_key
            ] = (
                candidate["a"],
                candidate["b"],
            )

            # Restrict the shared row and column.
            A2[i] = {
                candidate["a"]
            }

            B2[j] = {
                candidate["b"]
            }

            for result in recurse(
                A2,
                B2,
                assignment2,
            ):

                yield result

    for solution in recurse(
        {
            k: set(v)
            for k, v in A.items()
        },
        {
            k: set(v)
            for k, v in B.items()
        },
        {},
    ):

        yield solution

        if (
            solutions
            >= MAX_SOLUTIONS
        ):
            return


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

    # Hidden only for verification.
    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    # --------------------------------------------------------
    # Build cell candidate relations.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    cells = build_all_cells(
        n
    )

    build_time = (
        time.perf_counter()
        - t0
    )

    print()
    print(
        f"cell construction = "
        f"{build_time:.6f}s"
    )

    # --------------------------------------------------------
    # Initial domains.
    # --------------------------------------------------------

    A, B = project_domains(
        cells
    )

    # --------------------------------------------------------
    # AC-3.
    # --------------------------------------------------------

    t0 = time.perf_counter()

    consistent = arc_consistency(
        cells,
        A,
        B,
    )

    ac_time = (
        time.perf_counter()
        - t0
    )

    print()
    print(
        f"initial AC-3 = "
        f"{ac_time:.6f}s"
    )

    if not consistent:

        print(
            "GLOBAL CSP is inconsistent."
        )

        print()
        print(
            f"FINISHED INSTANCE "
            f"{bits}-BIT"
        )

        return

    # --------------------------------------------------------
    # Print domains.
    # --------------------------------------------------------

    print()
    print(
        "A domain sizes:"
    )

    for i, r in enumerate(R1):

        print(
            f"  A[{i}] mod {r}: "
            f"{len(A[i])}"
        )

    print()
    print(
        "B domain sizes:"
    )

    for j, s in enumerate(R2):

        print(
            f"  B[{j}] mod {s}: "
            f"{len(B[j])}"
        )

    # --------------------------------------------------------
    # True residue survival.
    # --------------------------------------------------------

    print()
    print(
        "True residue survival:"
    )

    for i, r in enumerate(R1):

        value = p % r

        print(
            f"  A[{i}] = {value}: "
            f"{value in A[i]}"
        )

    for j, s in enumerate(R2):

        value = q % s

        print(
            f"  B[{j}] = {value}: "
            f"{value in B[j]}"
        )

    # --------------------------------------------------------
    # Global CSP search.
    # --------------------------------------------------------

    print()
    print(
        "Starting GLOBAL CSP search..."
    )

    t0 = time.perf_counter()

    solution_count = 0

    factor_found = None

    for assignment in solve_global_csp(
        n,
        cells,
        A,
        B,
    ):

        solution_count += 1

        print()
        print(
            f"GLOBAL SOLUTION "
            f"{solution_count}"
        )

        residue_result = (
            assignment_to_residues(
                assignment
            )
        )

        if residue_result is None:

            print(
                "    inconsistent assignment"
            )

            continue

        recovered_A, recovered_B = (
            residue_result
        )

        print(
            f"    assigned cells = "
            f"{len(assignment)}/25"
        )

        print(
            f"    A = "
            f"{[recovered_A[i] for i in range(5)]}"
        )

        print(
            f"    B = "
            f"{[recovered_B[j] for j in range(5)]}"
        )

        # ----------------------------------------------------
        # Compare residue vector with hidden benchmark.
        # ----------------------------------------------------

        true_a = [
            p % r
            for r in R1
        ]

        true_b = [
            q % s
            for s in R2
        ]

        print(
            f"    true A = {true_a}"
        )

        print(
            f"    true B = {true_b}"
        )

        # ----------------------------------------------------
        # Generalized CRT factor recovery.
        # ----------------------------------------------------

        factor, stats = (
            try_factor_from_assignment(
                n,
                assignment,
            )
        )

        print(
            f"    CRT generated = "
            f"{stats['generated']}"
        )

        print(
            f"    exact tests = "
            f"{stats['exact_tests']}"
        )

        print(
            f"    CRT time = "
            f"{stats['time']:.6f}s"
        )

        if stats["aborted"]:

            print(
                f"    CRT status = "
                f"{stats['aborted']}"
            )

        if factor is not None:

            fp, fq = factor

            if fp > fq:
                fp, fq = fq, fp

            factor_found = (
                fp,
                fq,
            )

            print()
            print(
                "*** FACTOR FOUND ***"
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

            break

    total_time = (
        time.perf_counter()
        - t0
    )

    print()
    print(
        f"GLOBAL solutions tested = "
        f"{solution_count}"
    )

    print(
        f"GLOBAL search time = "
        f"{total_time:.6f}s"
    )

    if factor_found is None:

        print()
        print(
            "No factor recovered from "
            "global C consistency."
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
        "START EXPERIMENT 179B"
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
        f"MAX_SOLUTIONS = "
        f"{MAX_SOLUTIONS}"
    )

    print(
        f"MAX_NODE_COUNT = "
        f"{MAX_NODE_COUNT}"
    )

    print(
        f"MAX_EXACT_TESTS = "
        f"{MAX_EXACT_TESTS}"
    )

    print(
        f"MAX_GLOBAL_TIME = "
        f"{MAX_GLOBAL_TIME}"
    )

    print()

    # Quick self-test of generalized CRT.
    print(
        "Generalized CRT self-test:"
    )

    test = crt_pair_general(
        1,
        17,
        1,
        17 * 37,
    )

    print(
        f"    compatible test = {test}"
    )

    test2 = crt_pair_general(
        1,
        17,
        2,
        17 * 37,
    )

    print(
        f"    incompatible test = {test2}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 179B"
    )


if __name__ == "__main__":
    main()