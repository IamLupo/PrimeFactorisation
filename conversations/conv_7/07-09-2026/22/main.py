from math import inf


WINDOW = 100
DISPLACEMENT_LIMIT = 50

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
# BASIC PRIME / FACTOR ROUTINES
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
    Return unordered prime*prime factorizations of n.
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


# ============================================================
# ALGEBRAIC DISPLACEMENT SOLVER
# ============================================================

def displacement_solutions(px, qx, x, limit):
    """
    Solve the exact target equation

        -x = px*e + qx*d + d*e

    with

        |d| <= limit
        |e| <= limit

    Rearrangement:

        e = (-x - qx*d) / (px+d)

    This uses ONLY:
        px, qx, x

    and does NOT test divisibility of the original n.
    """

    solutions = []

    for d in range(-limit, limit + 1):

        p_candidate = px + d

        if p_candidate == 0:
            continue

        numerator = -x - qx * d

        if numerator % p_candidate != 0:
            continue

        e = numerator // p_candidate

        if abs(e) > limit:
            continue

        q_candidate = qx + e

        solutions.append(
            (
                d,
                e,
                p_candidate,
                q_candidate,
            )
        )

    return solutions


# ============================================================
# RANK / UNIQUENESS
# ============================================================

def displacement_norm(d, e):
    """
    Primary notion of "small displacement".
    """
    return abs(d) + abs(e)


def displacement_maxnorm(d, e):
    return max(abs(d), abs(e))


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 19")
print()

print(
    "n\tp\tq\t"
    "x\tpx\tqx\t"
    "true_d\ttrue_e\t"
    "solution_count\t"
    "true_found\t"
    "true_rank_L1\t"
    "true_rank_Linf\t"
    "unique_solution\t"
    "smallest_L1\t"
    "smallest_Linf\t"
    "smallest_pair"
)

summary = []


for p, q in TEST_CASES:

    n = p * q

    neighbor_states = 0
    true_found_count = 0
    unique_count = 0

    best_global_L1 = inf
    best_global_Linf = inf

    true_ranks = []

    examples = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        factors = prime_factor_pairs(nx)

        for px, qx in factors:

            neighbor_states += 1

            #
            # Ground truth.
            #
            true_d = p - px
            true_e = q - qx

            solutions = displacement_solutions(
                px,
                qx,
                x,
                DISPLACEMENT_LIMIT
            )

            #
            # Sort by L1 displacement.
            #
            sorted_L1 = sorted(
                solutions,
                key=lambda z: (
                    displacement_norm(z[0], z[1]),
                    displacement_maxnorm(z[0], z[1]),
                    abs(z[0]),
                    abs(z[1]),
                )
            )

            #
            # Sort by Linf displacement.
            #
            sorted_Linf = sorted(
                solutions,
                key=lambda z: (
                    displacement_maxnorm(z[0], z[1]),
                    displacement_norm(z[0], z[1]),
                    abs(z[0]),
                    abs(z[1]),
                )
            )

            true_index_L1 = None
            true_index_Linf = None

            for i, (d, e, pc, qc) in enumerate(sorted_L1):

                if d == true_d and e == true_e:
                    true_index_L1 = i + 1
                    break

            for i, (d, e, pc, qc) in enumerate(sorted_Linf):

                if d == true_d and e == true_e:
                    true_index_Linf = i + 1
                    break

            true_found = (
                true_index_L1 is not None
            )

            unique_solution = (
                len(solutions) == 1
                and true_found
            )

            if true_found:
                true_found_count += 1

                true_ranks.append(
                    (
                        true_index_L1,
                        true_index_Linf,
                    )
                )

            if unique_solution:
                unique_count += 1

            #
            # Record global best algebraic solutions.
            #
            if sorted_L1:

                best = sorted_L1[0]

                best_L1 = displacement_norm(
                    best[0],
                    best[1]
                )

                best_Linf = displacement_maxnorm(
                    best[0],
                    best[1]
                )

                best_global_L1 = min(
                    best_global_L1,
                    best_L1
                )

                best_global_Linf = min(
                    best_global_Linf,
                    best_Linf
                )

            #
            # Print only informative cases:
            #
            #   - true solution is present
            #   - and either unique or very high ranking
            #
            if true_found and (
                unique_solution
                or true_index_L1 <= 3
            ):
                smallest = sorted_L1[0]

                examples.append(
                    (
                        x,
                        px,
                        qx,
                        true_d,
                        true_e,
                        len(solutions),
                        true_index_L1,
                        true_index_Linf,
                        unique_solution,
                        smallest,
                    )
                )

    #
    # Aggregate ranks.
    #
    if true_ranks:

        avg_rank_L1 = (
            sum(r[0] for r in true_ranks)
            / len(true_ranks)
        )

        avg_rank_Linf = (
            sum(r[1] for r in true_ranks)
            / len(true_ranks)
        )

        best_rank_L1 = min(
            r[0] for r in true_ranks
        )

        best_rank_Linf = min(
            r[1] for r in true_ranks
        )

    else:
        avg_rank_L1 = inf
        avg_rank_Linf = inf
        best_rank_L1 = -1
        best_rank_Linf = -1

    found_pct = (
        100.0 * true_found_count / neighbor_states
        if neighbor_states
        else 0.0
    )

    unique_pct = (
        100.0 * unique_count / neighbor_states
        if neighbor_states
        else 0.0
    )

    summary.append(
        (
            n,
            neighbor_states,
            true_found_count,
            unique_count,
            found_pct,
            unique_pct,
            avg_rank_L1,
            avg_rank_Linf,
            best_rank_L1,
            best_rank_Linf,
        )
    )

    print()
    print(
        f"SUMMARY n={n}"
    )

    print(
        f"  nearby states = {neighbor_states}"
    )

    print(
        f"  true displacement found = "
        f"{true_found_count}"
    )

    print(
        f"  true displacement coverage = "
        f"{found_pct:.4f}%"
    )

    print(
        f"  unique algebraic solution = "
        f"{unique_count}"
    )

    print(
        f"  unique-solution rate = "
        f"{unique_pct:.4f}%"
    )

    if true_ranks:

        print(
            f"  average true L1 rank = "
            f"{avg_rank_L1:.4f}"
        )

        print(
            f"  average true Linf rank = "
            f"{avg_rank_Linf:.4f}"
        )

        print(
            f"  best true L1 rank = "
            f"{best_rank_L1}"
        )

        print(
            f"  best true Linf rank = "
            f"{best_rank_Linf}"
        )

    if best_global_L1 != inf:

        print(
            f"  smallest observed L1 = "
            f"{best_global_L1}"
        )

        print(
            f"  smallest observed Linf = "
            f"{best_global_Linf}"
        )

    #
    # Show at most 10 informative examples.
    #
    if examples:

        print()
        print(
            "  INFORMATIVE EXAMPLES"
        )
        print(
            "  x\tpx\tqx\t"
            "true_d\ttrue_e\t"
            "count\trank_L1\trank_Linf\t"
            "unique\tsmallest"
        )

        for row in examples[:10]:

            print(
                "  "
                f"{row[0]}\t"
                f"{row[1]}\t"
                f"{row[2]}\t"
                f"{row[3]}\t"
                f"{row[4]}\t"
                f"{row[5]}\t"
                f"{row[6]}\t"
                f"{row[7]}\t"
                f"{row[8]}\t"
                f"{row[9]}"
            )


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("FINAL SUMMARY")
print("-" * 120)

print(
    "n\t"
    "neighbors\t"
    "true_found\t"
    "unique\t"
    "true_found_pct\t"
    "unique_pct\t"
    "avg_L1_rank\t"
    "avg_Linf_rank\t"
    "best_L1_rank\t"
    "best_Linf_rank"
)

for (
    n,
    neighbors,
    true_found,
    unique,
    found_pct,
    unique_pct,
    avg_L1,
    avg_Linf,
    best_L1,
    best_Linf,
) in summary:

    avg_L1_text = (
        f"{avg_L1:.4f}"
        if avg_L1 != inf
        else "NA"
    )

    avg_Linf_text = (
        f"{avg_Linf:.4f}"
        if avg_Linf != inf
        else "NA"
    )

    print(
        f"{n}\t"
        f"{neighbors}\t"
        f"{true_found}\t"
        f"{unique}\t"
        f"{found_pct:.4f}\t"
        f"{unique_pct:.4f}\t"
        f"{avg_L1_text}\t"
        f"{avg_Linf_text}\t"
        f"{best_L1}\t"
        f"{best_Linf}"
    )


print()
print("FINISHED EXPERIMENT 19")
