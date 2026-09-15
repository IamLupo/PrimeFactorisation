from math import isqrt

WINDOW = 100

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
    Return all factor pairs (p, q) where both are prime.
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


def factor_sum(p, q):
    return p + q


def discriminant_gap(n, s):
    """
    For a factorization n = pq with s = p+q:

        s^2 - 4n = (p-q)^2
    """
    return s * s - 4 * n


def trajectory_for_target(n):
    """
    Collect nearby prime-factorized neighbors.

    Returns:
        (x, px, qx, Sx, Dx)
    where

        n+x = px*qx
        Sx   = px+qx
        Dx   = Sx^2 - 4(n+x) = (px-qx)^2
    """
    states = []

    for x in range(-WINDOW, WINDOW + 1):
        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        factors = prime_factor_pairs(nx)

        for px, qx in factors:
            sx = factor_sum(px, qx)
            dx = discriminant_gap(nx, sx)

            states.append(
                {
                    "x": x,
                    "px": px,
                    "qx": qx,
                    "sx": sx,
                    "dx": dx,
                }
            )

    return states


def local_prediction_errors(states, target_s):
    """
    Measure how well the simplest local transformations
    predict the unknown target sum S.

    No divisibility tests on n are used here.

    Candidate predictors:

        1. Sx
        2. Sx - x
        3. Sx + x
        4. Sx - round(x / Sx)
        5. Sx + round(x / Sx)
        6. Sx - round(Dx / (2*Sx))
        7. Sx + round(Dx / (2*Sx))

    These are deliberately simple deterministic transforms.
    """

    rows = []

    for state in states:
        x = state["x"]
        sx = state["sx"]
        dx = state["dx"]

        predictions = {}

        predictions["Sx"] = sx
        predictions["Sx_minus_x"] = sx - x
        predictions["Sx_plus_x"] = sx + x

        if sx != 0:
            q = round(x / sx)

            predictions["Sx_minus_x_over_Sx"] = sx - q
            predictions["Sx_plus_x_over_Sx"] = sx + q

            q2 = round(dx / (2 * sx))

            predictions["Sx_minus_D_over_2S"] = sx - q2
            predictions["Sx_plus_D_over_2S"] = sx + q2

        errors = {
            name: abs(pred - target_s)
            for name, pred in predictions.items()
        }

        rows.append(
            {
                **state,
                "predictions": predictions,
                "errors": errors,
            }
        )

    return rows


print("START EXPERIMENT 20")
print()

print(
    "n\tp\tq\tS\t"
    "neighbors\t"
    "best_Sx_error\t"
    "best_Sx_predictor\t"
    "exact_prediction_count"
)

all_summary = []

for p, q in TEST_CASES:
    n = p * q
    target_s = p + q

    states = trajectory_for_target(n)

    scored = local_prediction_errors(
        states,
        target_s
    )

    predictor_names = [
        "Sx",
        "Sx_minus_x",
        "Sx_plus_x",
        "Sx_minus_x_over_Sx",
        "Sx_plus_x_over_Sx",
        "Sx_minus_D_over_2S",
        "Sx_plus_D_over_2S",
    ]

    predictor_errors = {
        name: []
        for name in predictor_names
    }

    exact_counts = {
        name: 0
        for name in predictor_names
    }

    for row in scored:
        for name in predictor_names:
            error = row["errors"].get(name)

            if error is None:
                continue

            predictor_errors[name].append(error)

            if error == 0:
                exact_counts[name] += 1

    best_predictor = None
    best_error = None

    for name in predictor_names:
        errors = predictor_errors[name]

        if not errors:
            continue

        current_best = min(errors)

        if best_error is None or current_best < best_error:
            best_error = current_best
            best_predictor = name

    total_exact = sum(exact_counts.values())

    print(
        f"{n}\t{p}\t{q}\t{target_s}\t"
        f"{len(states)}\t"
        f"{best_error}\t"
        f"{best_predictor}\t"
        f"{total_exact}"
    )

    print()
    print(f"DETAIL n={n}")
    print(f"  target S = {target_s}")
    print(f"  neighbor states = {len(states)}")
    print()

    print("  PREDICTOR PERFORMANCE")

    for name in predictor_names:
        errors = predictor_errors[name]

        if not errors:
            continue

        exact = exact_counts[name]

        average_error = sum(errors) / len(errors)
        median_sorted = sorted(errors)
        median_error = median_sorted[len(median_sorted) // 2]
        minimum_error = min(errors)

        print(
            f"  {name}: "
            f"min={minimum_error}, "
            f"median={median_error}, "
            f"avg={average_error:.4f}, "
            f"exact={exact}"
        )

    print()
    print("  CLOSEST PREDICTIONS")

    closest = sorted(
        scored,
        key=lambda row: min(row["errors"].values())
    )[:10]

    print(
        "  x\tpx\tqx\tSx\t"
        "best_predictor\t"
        "prediction\t"
        "target_S\t"
        "error"
    )

    for row in closest:
        best_name = min(
            row["errors"],
            key=row["errors"].get
        )

        best_pred = row["predictions"][best_name]
        best_err = row["errors"][best_name]

        print(
            f"  {row['x']}\t"
            f"{row['px']}\t"
            f"{row['qx']}\t"
            f"{row['sx']}\t"
            f"{best_name}\t"
            f"{best_pred}\t"
            f"{target_s}\t"
            f"{best_err}"
        )

    all_summary.append(
        {
            "n": n,
            "p": p,
            "q": q,
            "target_s": target_s,
            "neighbors": len(states),
            "best_error": best_error,
            "best_predictor": best_predictor,
            "total_exact": total_exact,
        }
    )

print()
print("FINAL SUMMARY")
print("-" * 110)

print(
    "n\t"
    "target_S\t"
    "neighbors\t"
    "best_error\t"
    "best_predictor\t"
    "total_exact_predictions"
)

for row in all_summary:
    print(
        f"{row['n']}\t"
        f"{row['target_s']}\t"
        f"{row['neighbors']}\t"
        f"{row['best_error']}\t"
        f"{row['best_predictor']}\t"
        f"{row['total_exact']}"
    )

print()
print("FINISHED EXPERIMENT 20")
