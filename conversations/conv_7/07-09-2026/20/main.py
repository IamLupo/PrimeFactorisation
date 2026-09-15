from math import isqrt


WINDOW = 100

#
# Maximum displacement hypothesized around a nearby factor.
#
DISPLACEMENT_LIMIT = 10


TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
    (47, 83),
    (53, 89),
    (59, 97),
    (61, 101),
    (67, 103),
    (71, 107),
]


# ============================================================
# NUMBER THEORY
# ============================================================

def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def prime_factor_pairs(n):
    """
    Return all prime factor pairs p <= q of n.
    """
    result = []

    d = 2

    while d * d <= n:

        if n % d == 0:

            q = n // d

            if is_prime(d) and is_prime(q):
                result.append((d, q))

        d += 1

    return result


def gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)


# ============================================================
# CANDIDATE GENERATION
# ============================================================

def candidates_from_neighbor(n, px, qx):
    """
    Use a nearby factor px of n+x.

    Hypothesis:

        p = px + d

    for small integer d.

    For every candidate p:

        p | n

    is required.

    IMPORTANT:
    This function knows only n, px, qx and the displacement
    bound. It does not use the target factorization.
    """

    candidates = []

    for d in range(
        -DISPLACEMENT_LIMIT,
        DISPLACEMENT_LIMIT + 1
    ):

        p_candidate = px + d

        if p_candidate <= 1:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = n // p_candidate

        candidates.append(
            (
                p_candidate,
                q_candidate,
                d,
            )
        )

    return candidates


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 17")
print()

print(
    "n\tp\tq\t"
    "neighbor_x\t"
    "neighbor_px\t"
    "neighbor_qx\t"
    "true_dp\t"
    "true_dq\t"
    "candidate_count\t"
    "target_found\t"
    "candidate_p_values"
)


summary = []


for p, q in TEST_CASES:

    n = p * q

    total_neighbor_states = 0
    useful_neighbors = 0
    target_found_count = 0

    min_candidate_count = None
    best_neighbor = None

    for x in range(
        -WINDOW,
        WINDOW + 1
    ):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = prime_factor_pairs(nx)

        for px, qx in pairs:

            total_neighbor_states += 1

            #
            # Ground truth only for measuring
            # the displacement.
            #
            true_dp = px - p
            true_dq = qx - q

            #
            # Generate candidates using only
            # nearby factor px.
            #
            candidates = candidates_from_neighbor(
                n,
                px,
                qx
            )

            candidate_count = len(candidates)

            target_found = any(
                cp == p
                for cp, cq, d in candidates
            )

            if target_found:

                target_found_count += 1

            if candidate_count > 0:

                useful_neighbors += 1

                if (
                    min_candidate_count is None
                    or candidate_count < min_candidate_count
                ):
                    min_candidate_count = candidate_count
                    best_neighbor = (
                        x,
                        px,
                        qx,
                        true_dp,
                        true_dq,
                        candidates,
                    )

            candidate_values = [
                cp
                for cp, cq, d in candidates
            ]

            print(
                f"{n}\t"
                f"{p}\t"
                f"{q}\t"
                f"{x}\t"
                f"{px}\t"
                f"{qx}\t"
                f"{true_dp}\t"
                f"{true_dq}\t"
                f"{candidate_count}\t"
                f"{target_found}\t"
                f"{candidate_values}"
            )

    #
    # ========================================================
    # SUMMARY FOR THIS n
    # ========================================================
    #

    print()
    print(
        f"SUMMARY n={n}"
    )

    if best_neighbor is None:

        print(
            "  No nearby factor produced a candidate."
        )

        summary.append(
            (
                n,
                p,
                q,
                total_neighbor_states,
                useful_neighbors,
                target_found_count,
                None,
                None,
            )
        )

        continue

    bx, bpx, bqx, bdp, bdq, bcandidates = (
        best_neighbor
    )

    print(
        f"  Best neighbor: x={bx}, "
        f"px={bpx}, qx={bqx}"
    )

    print(
        f"  True displacement: "
        f"dp={bdp}, dq={bdq}"
    )

    print(
        f"  Candidate count: "
        f"{len(bcandidates)}"
    )

    print(
        f"  Candidates: "
        f"{bcandidates}"
    )

    print(
        f"  Target found by some neighbor: "
        f"{target_found_count > 0}"
    )

    summary.append(
        (
            n,
            p,
            q,
            total_neighbor_states,
            useful_neighbors,
            target_found_count,
            len(bcandidates),
            bx,
        )
    )


print()
print("GLOBAL SUMMARY")
print("-" * 120)

print(
    "n\tp\tq\t"
    "neighbor_states\t"
    "useful_neighbors\t"
    "target_hits\t"
    "best_candidate_count\t"
    "best_x"
)

for row in summary:

    print(
        f"{row[0]}\t"
        f"{row[1]}\t"
        f"{row[2]}\t"
        f"{row[3]}\t"
        f"{row[4]}\t"
        f"{row[5]}\t"
        f"{row[6]}\t"
        f"{row[7]}"
    )


print()
print("INTERSECTION TEST")
print("-" * 120)

print(
    "n\t"
    "intersection_size\t"
    "intersection\t"
    "contains_true_p"
)

for p, q in TEST_CASES:

    n = p * q

    intersection = None

    observations = 0

    for x in range(
        -WINDOW,
        WINDOW + 1
    ):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = prime_factor_pairs(nx)

        for px, qx in pairs:

            observations += 1

            candidates = candidates_from_neighbor(
                n,
                px,
                qx
            )

            candidate_set = {
                cp
                for cp, cq, d in candidates
            }

            if intersection is None:
                intersection = candidate_set
            else:
                intersection &= candidate_set

    if intersection is None:
        intersection = set()

    print(
        f"{n}\t"
        f"{len(intersection)}\t"
        f"{sorted(intersection)}\t"
        f"{p in intersection}"
    )


print()
print("FINISHED EXPERIMENT 17")
