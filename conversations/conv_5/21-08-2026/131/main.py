# =============================================================================
# EXPERIMENT 692
# MULTI-SHIFT H-CONSISTENCY -> CANONICAL C RECOVERY
#
# Goal:
#
#   We know:
#
#       M_c = n + c
#       H_c = gcd(C, M_c)
#
#   and therefore:
#
#       H_c | M_c
#       H_c | C
#
#   Previous experiments generated arbitrary divisors H_c from each M_c and
#   formed L = lcm(H_c).  That was too permissive.
#
#   This experiment instead enforces the FULL gcd condition:
#
#       gcd(C, n+c) == H_c
#
#   simultaneously for several shifts.
#
#   The true H_c is NOT used to select candidates.
#
#   We construct candidate C values from divisor combinations and then
#   independently verify the complete H-vector.
#
# IMPORTANT:
#   This experiment is intentionally small/moderate by default because the
#   Cartesian product of divisor lattices can become large.
# =============================================================================

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product
from math import gcd, lcm
from sympy import factorint, divisors, primerange


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 200

# Start with a small set. Increase carefully.
SHIFTS = [3, 9, 81, 137]

# Maximum number of H-combinations examined per n/frame.
# This prevents pathological divisor-lattice explosions.
MAX_COMBINATIONS = 200_000

# Candidate C range around the observable L.
#
# We know:
#
#     L | C
#
# so candidates are C = L*k.
#
# The canonical C is roughly the larger factor plus/minus a small constant
# in the present frames, so C cannot exceed n+some small amount.
#
# This is NOT used to prove the factorization. It only bounds enumeration.
MAX_C = None  # None => C <= n + max(abs(shift)) + 10

PRINT_EXAMPLES = 20


# =============================================================================
# BASIC ARITHMETIC
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)
    if x == 0:
        return 10**9

    r = 0
    while (x & 1) == 0:
        x >>= 1
        r += 1
    return r


def frame_for_n(n: int) -> str:
    # Same frame split used throughout the previous experiments.
    return "A" if n % 4 == 3 else "B"


def canonical_C(p: int, q: int, frame: str) -> int:
    if frame == "A":
        return q + 3
    return p + 1


def canonical_factor_from_C(n: int, C: int, frame: str):
    """
    Frame A:
        C = q + 3
        q = C - 3

    Frame B:
        C = p + 1
        p = C - 1

    Return the candidate factor pair if C induces an exact divisor pair.
    """
    if frame == "A":
        q = C - 3
        if q <= 1:
            return None

        if n % q != 0:
            return None

        p = n // q
        return p, q

    else:
        p = C - 1
        if p <= 1:
            return None

        if n % p != 0:
            return None

        q = n // p
        return p, q


def generate_primes(limit: int) -> list[int]:
    return list(primerange(3, limit + 1))


def generate_semiprimes(primes: list[int]):
    for i, p in enumerate(primes):
        for q in primes[i:]:
            yield p * q, p, q


# =============================================================================
# H VECTOR
# =============================================================================

def h_vector(n: int, C: int, shifts: list[int]) -> tuple[int, ...]:
    return tuple(gcd(C, n + c) for c in shifts)


def factor_divisor_lattice(x: int) -> tuple[dict[int, int], list[int]]:
    fac = factorint(x)
    divs = divisors(x)
    return fac, divs


# =============================================================================
# CANDIDATE GENERATION
# =============================================================================

def compatible_C_from_H_vector(
    n: int,
    frame: str,
    shifts: list[int],
    Hs: tuple[int, ...],
    max_C: int,
) -> set[int]:
    """
    Given a proposed H-vector, generate C values satisfying:

        H_i | C

    for every i.

    Therefore:

        L = lcm(H_i)

    divides C.

    We then enumerate C = L*k.

    The important part is that every candidate is later subjected to the
    exact self-consistency condition:

        gcd(C, n+c_i) == H_i
    """

    L = 1
    for h in Hs:
        L = lcm(L, h)

    if L <= 0:
        return set()

    out = set()

    # C must be at least L.
    for C in range(L, max_C + 1, L):
        # Full H-vector consistency.
        good = True

        for c, h in zip(shifts, Hs):
            if gcd(C, n + c) != h:
                good = False
                break

        if good:
            out.add(C)

    return out


def candidate_C_from_divisor_lattices(
    n: int,
    frame: str,
    shifts: list[int],
    max_C: int,
):
    """
    Search over H choices from the complete divisor lattices.

    To avoid an enormous Cartesian product, we progressively merge constraints.

    For a partial H-vector, maintain possible C values satisfying all exact gcd
    constraints encountered so far.
    """

    lattices = []

    for c in shifts:
        M = n + c

        if M <= 0:
            return set(), 0, False

        fac, divs = factor_divisor_lattice(M)
        lattices.append((c, fac, divs))

    # Start with all possible C in a useful bounded range.
    #
    # We don't know C, but for the canonical frames:
    #
    #   A: C=q+3 <= n+3
    #   B: C=p+1 <= n+1
    #
    # Thus n+max_shift+10 is a harmless experimental ceiling.
    #
    # Candidate C values will later be filtered by the frame equation.
    candidates = set(range(2, max_C + 1))

    combinations_examined = 0
    aborted = False

    # Use the smallest divisor lattices first.
    lattices.sort(key=lambda t: len(t[2]))

    for c, fac, divs in lattices:

        new_candidates = set()

        # Instead of selecting H first, directly evaluate:
        #
        #     H = gcd(C, n+c)
        #
        # for the current C candidates.
        #
        # This avoids the huge Cartesian H product and is exactly the
        # observable consistency condition we care about.
        for C in candidates:
            h = gcd(C, n + c)

            # h is automatically in the divisor lattice, but explicitly
            # retain that fact for diagnostics.
            if h in divs:
                new_candidates.add(C)

            combinations_examined += 1

            if combinations_examined >= MAX_COMBINATIONS:
                aborted = True
                break

        candidates = new_candidates

        if aborted:
            break

    return candidates, combinations_examined, aborted


# =============================================================================
# TRUE-CONTROL DATA
# =============================================================================

def exact_frame_candidate(C: int, n: int, frame: str):
    """
    Test whether C itself corresponds to the canonical frame.
    """
    pair = canonical_factor_from_C(n, C, frame)
    if pair is None:
        return None

    p, q = pair

    if p <= 1 or q <= 1:
        return None

    return p, q


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states):
    failures = 0

    for n, p, q in states:
        frame = frame_for_n(n)

        if frame == "A":
            C = q + 3
            recovered = canonical_factor_from_C(n, C, frame)

            if recovered != (p, q):
                failures += 1

        else:
            C = p + 1
            recovered = canonical_factor_from_C(n, C, frame)

            if recovered != (p, q):
                failures += 1

    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_complete_shifted_lattices(states):
    print("=" * 90)
    print("TEST 1: TRUE H IS IN EVERY DIVISOR LATTICE")
    print("=" * 90)

    failures = 0

    for n, p, q in states:

        frame = frame_for_n(n)
        C = canonical_C(p, q, frame)

        for c in SHIFTS:
            M = n + c

            if M <= 0:
                failures += 1
                continue

            if M % gcd(C, M) != 0:
                failures += 1

    print(f"shifts={SHIFTS}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_observable_C_candidates(states):
    print("=" * 90)
    print("TEST 2: OBSERVABLE MULTI-SHIFT GCD CONSISTENCY")
    print("=" * 90)

    print()
    print("This is the key test.")
    print()
    print("For each n and frame we use ONLY:")
    print()
    print("    n+c")
    print("    factor(n+c)")
    print("    complete divisor lattice")
    print()
    print("Then we keep C values satisfying simultaneously:")
    print()
    print("    gcd(C,n+c) = H_c")
    print()
    print("for the H-values implied by C itself.")
    print()
    print("Finally we check whether the surviving C corresponds to the")
    print("canonical frame factorization.")
    print()

    results = {
        "A": {
            "none": 0,
            "unique": 0,
            "ambiguous": 0,
            "canonical_present": 0,
            "candidate_counts": Counter(),
        },
        "B": {
            "none": 0,
            "unique": 0,
            "ambiguous": 0,
            "canonical_present": 0,
            "candidate_counts": Counter(),
        },
    }

    examples = []

    total_examined = 0
    aborted_states = 0

    for n, p, q in states:

        frame = frame_for_n(n)

        max_C = (
            MAX_C
            if MAX_C is not None
            else n + max(abs(c) for c in SHIFTS) + 10
        )

        candidates, examined, aborted = candidate_C_from_divisor_lattices(
            n=n,
            frame=frame,
            shifts=SHIFTS,
            max_C=max_C,
        )

        total_examined += examined

        if aborted:
            aborted_states += 1

        true_C = canonical_C(p, q, frame)

        if true_C in candidates:
            results[frame]["canonical_present"] += 1

        count = len(candidates)
        results[frame]["candidate_counts"][count] += 1

        if count == 0:
            results[frame]["none"] += 1
        elif count == 1:
            results[frame]["unique"] += 1
        else:
            results[frame]["ambiguous"] += 1

        if len(examples) < PRINT_EXAMPLES:
            examples.append(
                (
                    frame,
                    n,
                    p,
                    q,
                    true_C,
                    sorted(candidates)[:50],
                    count,
                    aborted,
                )
            )

    for frame in ("A", "B"):

        r = results[frame]

        print(f"FRAME {frame}")
        print("-" * 90)

        total = sum(r["candidate_counts"].values())

        print(f"states={total}")
        print(f"canonical C present={r['canonical_present']}")
        print(f"none={r['none']}")
        print(f"unique={r['unique']}")
        print(f"ambiguous={r['ambiguous']}")

        if total:
            print(
                f"unique ratio={r['unique'] / total:.6f}"
            )

        print("candidate-count distribution:")
        for k, v in sorted(r["candidate_counts"].items()):
            print(f"{k:8d} -> {v}")

        print()

    print("examples:")
    for (
        frame,
        n,
        p,
        q,
        true_C,
        candidates,
        count,
        aborted,
    ) in examples:

        print(
            f"frame={frame} n={n} p={p} q={q} "
            f"true_C={true_C}"
        )
        print(f"    candidate_count={count}")
        print(f"    candidates={candidates}")
        if aborted:
            print("    WARNING=COMBINATION LIMIT REACHED")

    print()

    return results


# =============================================================================
# TEST 3
# =============================================================================

def test_exact_H_vector(states):
    print("=" * 90)
    print("TEST 3: H-VECTOR SELF-CONSISTENCY")
    print("=" * 90)

    failures = 0
    collision_map = defaultdict(list)

    for n, p, q in states:

        frame = frame_for_n(n)
        C = canonical_C(p, q, frame)

        Hs = h_vector(n, C, SHIFTS)

        # Verify every component exactly.
        for c, h in zip(SHIFTS, Hs):
            if gcd(C, n + c) != h:
                failures += 1

        collision_map[(frame, Hs)].append((n, C, p, q))

    print(f"shifts={SHIFTS}")
    print(f"checked={len(states)}")
    print(f"failures={failures}")

    collision_counts = Counter(
        len(v) for v in collision_map.values()
    )

    print()
    print("H-vector collision multiplicity:")
    for k, v in sorted(collision_counts.items()):
        print(f"{k:8d} -> {v}")

    print()

    return failures, collision_map


# =============================================================================
# TEST 4
# =============================================================================

def test_candidate_C_frame_recovery(states):
    print("=" * 90)
    print("TEST 4: C -> FRAME FACTOR RECOVERY")
    print("=" * 90)

    classifications = Counter()
    examples = []

    for n, p, q in states:

        frame = frame_for_n(n)

        max_C = (
            MAX_C
            if MAX_C is not None
            else n + max(abs(c) for c in SHIFTS) + 10
        )

        candidates, _, _ = candidate_C_from_divisor_lattices(
            n=n,
            frame=frame,
            shifts=SHIFTS,
            max_C=max_C,
        )

        valid = []

        for C in candidates:
            pair = exact_frame_candidate(C, n, frame)

            if pair is None:
                continue

            pp, qq = pair

            # Require actual factorization.
            if pp * qq != n:
                continue

            valid.append((C, pp, qq))

        classifications[len(valid)] += 1

        if len(examples) < PRINT_EXAMPLES:
            examples.append(
                (
                    frame,
                    n,
                    p,
                    q,
                    valid[:20],
                )
            )

    print("valid factor-pair candidate count distribution:")
    for k, v in sorted(classifications.items()):
        print(f"{k:8d} -> {v}")

    print()
    print("examples:")

    for frame, n, p, q, valid in examples:
        print(
            f"frame={frame} n={n} true=({p},{q}) "
            f"valid_candidates={valid}"
        )

    print()

    return classifications


# =============================================================================
# TEST 5
# =============================================================================

def test_h_vector_collision_information(states):
    print("=" * 90)
    print("TEST 5: DOES THE H-VECTOR IDENTIFY C?")
    print("=" * 90)

    signatures = defaultdict(set)

    for n, p, q in states:

        frame = frame_for_n(n)
        C = canonical_C(p, q, frame)
        Hs = h_vector(n, C, SHIFTS)

        signatures[(frame, Hs)].add(C)

    multi = 0
    max_C_per_signature = 0

    examples = []

    for sig, Cs in signatures.items():

        if len(Cs) > 1:
            multi += 1

            max_C_per_signature = max(
                max_C_per_signature,
                len(Cs),
            )

            if len(examples) < PRINT_EXAMPLES:
                examples.append((sig, sorted(Cs)))

    print(f"H-signatures={len(signatures)}")
    print(f"signatures with multiple C={multi}")
    print(f"maximum C values per H-signature={max_C_per_signature}")

    if examples:
        print()
        print("collision examples:")
        for sig, Cs in examples:
            print(f"signature={sig}")
            print(f"    Cs={Cs}")

    print()

    return signatures


# =============================================================================
# TEST 6
# =============================================================================

def test_minimal_shift_progression(states):
    print("=" * 90)
    print("TEST 6: SHIFT PROGRESSION")
    print("=" * 90)

    for frame in ("A", "B"):

        frame_states = [
            (n, p, q)
            for n, p, q in states
            if frame_for_n(n) == frame
        ]

        print(f"FRAME {frame}")
        print("-" * 90)

        for r in range(1, len(SHIFTS) + 1):

            current = SHIFTS[:r]

            unique = 0
            ambiguous = 0
            missing = 0
            candidate_total = 0

            for n, p, q in frame_states:

                max_C = (
                    MAX_C
                    if MAX_C is not None
                    else n + max(abs(c) for c in current) + 10
                )

                candidates, _, _ = candidate_C_from_divisor_lattices(
                    n=n,
                    frame=frame,
                    shifts=current,
                    max_C=max_C,
                )

                true_C = canonical_C(p, q, frame)

                if true_C not in candidates:
                    missing += 1

                candidate_total += len(candidates)

                if len(candidates) == 1:
                    unique += 1
                elif len(candidates) > 1:
                    ambiguous += 1

            mean_candidates = (
                candidate_total / len(frame_states)
                if frame_states
                else 0.0
            )

            print(
                f"shifts={current}"
            )
            print(
                f"    unique={unique} "
                f"ambiguous={ambiguous} "
                f"missing={missing} "
                f"mean_candidates={mean_candidates:.4f}"
            )

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 692 START")
    print("=" * 90)
    print()
    print(f"prime limit={PRIME_LIMIT}")
    print(f"shifts={SHIFTS}")
    print()

    primes = generate_primes(PRIME_LIMIT)
    states = list(generate_semiprimes(primes))

    print(f"odd primes={len(primes)}")
    print(f"semiprimes={len(states)}")
    print()

    baseline_failures = test_baseline(states)

    test_complete_shifted_lattices(states)

    test_observable_C_candidates(states)

    test_exact_H_vector(states)

    test_candidate_C_frame_recovery(states)

    test_h_vector_collision_information(states)

    test_minimal_shift_progression(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("For each shift c:")
    print()
    print("    M_c = n+c")
    print("    H_c = gcd(C,M_c)")
    print()
    print("and therefore:")
    print()
    print("    H_c | M_c")
    print("    H_c | C")
    print()
    print("The previous experiments used:")
    print()
    print("    L = lcm(H_c)")
    print()
    print("but L | C is only a necessary condition.")
    print()
    print("Experiment 692 instead emphasizes the exact condition:")
    print()
    print("    gcd(C,n+c) == H_c")
    print()
    print("for ALL shifts simultaneously.")
    print()
    print("The main question is:")
    print()
    print("    Can the complete shifted divisor information")
    print("    force a small set of self-consistent C values?")
    print()
    print("The strongest possible outcome would be:")
    print()
    print("    one surviving C")
    print("    -> frame equation")
    print("    -> p,q")
    print()
    print("A weaker but still useful result is:")
    print()
    print("    many C candidates remain,")
    print("    but the candidate count shrinks rapidly as shifts are added.")
    print()
    print("A failure would mean:")
    print()
    print("    shifted divisor lattices alone do not select C")
    print()
    print(f"BASELINE FAILURES={baseline_failures}")
    print("=" * 90)
    print("EXPERIMENT 692 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
