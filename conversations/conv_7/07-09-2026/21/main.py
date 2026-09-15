from math import isqrt


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
# BASIC NUMBER THEORY
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
# SOLVE DISPLACEMENT EQUATION
# ============================================================

def predicted_displacements(px, qx, x):
    """
    Solve

        x = px*e + qx*d + d*e

    for small integer d,e.

    Rearrange:

        x = qx*d + e*(px+d)

    Therefore:

        e = (x - qx*d)/(px+d)

    for every d where px+d != 0.

    IMPORTANT:
    This prediction step does NOT use n.
    """

    candidates = []

    for d in range(
        -DISPLACEMENT_LIMIT,
        DISPLACEMENT_LIMIT + 1
    ):

        p_candidate = px + d

        if p_candidate == 0:
            continue

        numerator = (
            x
            - qx * d
        )

        if numerator % p_candidate != 0:
            continue

        e = numerator // p_candidate

        candidates.append(
            (
                d,
                e,
                px + d,
                qx + e,
            )
        )

    return candidates


# ============================================================
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 18")
print()

print(
    "n\tp\tq\t"
    "x\tpx\tqx\t"
    "true_d\ttrue_e\t"
    "prediction_count\t"
    "true_pair_predicted\t"
    "predicted_pairs"
)

summary = []


for p, q in TEST_CASES:

    n = p * q

    total_neighbors = 0
    prediction_events = 0
    true_predicted = 0

    for x in range(
        -WINDOW,
        WINDOW + 1
    ):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        for px, qx in prime_factor_pairs(nx):

            total_neighbors += 1

            #
            # Ground truth ONLY for evaluation.
            #
            true_d = p - px
            true_e = q - qx

            predictions = predicted_displacements(
                px,
                qx,
                x
            )

            prediction_count = len(
                predictions
            )

            true_pair_predicted = any(
                d == true_d
                and e == true_e
                for d, e, pc, qc
                in predictions
            )

            if prediction_count > 0:
                prediction_events += 1

            if true_pair_predicted:
                true_predicted += 1

            print(
                f"{n}\t"
                f"{p}\t"
                f"{q}\t"
                f"{x}\t"
                f"{px}\t"
                f"{qx}\t"
                f"{true_d}\t"
                f"{true_e}\t"
                f"{prediction_count}\t"
                f"{true_pair_predicted}\t"
                f"{predictions}"
            )

    print()
    print(
        f"SUMMARY n={n}"
    )
    print(
        f"  nearby states = {total_neighbors}"
    )
    print(
        f"  states with predictions = "
        f"{prediction_events}"
    )
    print(
        f"  states containing true displacement = "
        f"{true_predicted}"
    )

    if total_neighbors:
        print(
            f"  prediction coverage = "
            f"{100.0 * prediction_events / total_neighbors:.4f}%"
        )

        print(
            f"  true displacement coverage = "
            f"{100.0 * true_predicted / total_neighbors:.4f}%"
        )

    summary.append(
        (
            n,
            total_neighbors,
            prediction_events,
            true_predicted,
        )
    )


print()
print("SUMMARY TABLE")
print("-" * 100)

print(
    "n\t"
    "neighbor_states\t"
    "prediction_states\t"
    "true_displacement_predicted\t"
    "prediction_coverage_pct\t"
    "true_coverage_pct"
)

for n, total, predicted, true_predicted in summary:

    prediction_pct = (
        100.0 * predicted / total
        if total else 0.0
    )

    true_pct = (
        100.0 * true_predicted / total
        if total else 0.0
    )

    print(
        f"{n}\t"
        f"{total}\t"
        f"{predicted}\t"
        f"{true_predicted}\t"
        f"{prediction_pct:.4f}\t"
        f"{true_pct:.4f}"
    )


print()
print("FINISHED EXPERIMENT 18")
