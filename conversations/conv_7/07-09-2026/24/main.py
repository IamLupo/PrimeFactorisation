from math import inf


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
    result = []

    d = 2

    while d * d <= n:
        if n % d == 0:
            q = n // d

            if is_prime(d) and is_prime(q):
                result.append((d, q))

        d += 1

    return result


def collect_neighbors(n):
    """
    Collect prime-factorized neighbors:

        n + x = px*qx

    and their factor sums:

        Sx = px + qx
    """

    states = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        for px, qx in prime_factor_pairs(nx):

            sx = px + qx

            states.append(
                {
                    "x": x,
                    "px": px,
                    "qx": qx,
                    "sx": sx,
                }
            )

    return states


def choose_closest(states, side=None):
    """
    Choose exactly one neighbor.

    side=None:
        closest |x| overall

    side="negative":
        closest negative x

    side="positive":
        closest positive x
    """

    candidates = states

    if side == "negative":
        candidates = [
            s for s in states
            if s["x"] < 0
        ]

    elif side == "positive":
        candidates = [
            s for s in states
            if s["x"] > 0
        ]

    if not candidates:
        return None

    return min(
        candidates,
        key=lambda s: abs(s["x"])
    )


def interpolate_to_zero(left, right):
    """
    Linear interpolation of S(x) to x=0.

    left.x < 0
    right.x > 0

    S_hat(0) =
        (S_left * x_right - S_right * x_left)
        / (x_right - x_left)
    """

    xl = left["x"]
    xr = right["x"]

    sl = left["sx"]
    sr = right["sx"]

    return (
        sl * xr - sr * xl
    ) / (
        xr - xl
    )


def nearest_integer(x):
    """
    Symmetric nearest integer rounding.
    """
    return int(x + 0.5)


print("START EXPERIMENT 21")
print()

print(
    "n\tp\tq\tS\t"
    "nearest_x\tnearest_Sx\tnearest_error\t"
    "neg_x\tneg_Sx\tneg_error\t"
    "pos_x\tpos_Sx\tpos_error\t"
    "interp_prediction\tinterp_error\t"
    "two_sided"
)

summary = []


for p, q in TEST_CASES:

    n = p * q
    target_s = p + q

    states = collect_neighbors(n)

    nearest = choose_closest(states)
    negative = choose_closest(states, "negative")
    positive = choose_closest(states, "positive")

    nearest_error = (
        abs(nearest["sx"] - target_s)
        if nearest else None
    )

    negative_error = (
        abs(negative["sx"] - target_s)
        if negative else None
    )

    positive_error = (
        abs(positive["sx"] - target_s)
        if positive else None
    )

    interpolation = None
    interpolation_error = None

    if negative is not None and positive is not None:

        interpolation = interpolate_to_zero(
            negative,
            positive
        )

        interpolation_error = abs(
            interpolation - target_s
        )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{target_s}\t"
        f"{nearest['x'] if nearest else 'NA'}\t"
        f"{nearest['sx'] if nearest else 'NA'}\t"
        f"{nearest_error if nearest_error is not None else 'NA'}\t"
        f"{negative['x'] if negative else 'NA'}\t"
        f"{negative['sx'] if negative else 'NA'}\t"
        f"{negative_error if negative_error is not None else 'NA'}\t"
        f"{positive['x'] if positive else 'NA'}\t"
        f"{positive['sx'] if positive else 'NA'}\t"
        f"{positive_error if positive_error is not None else 'NA'}\t"
        f"{interpolation if interpolation is not None else 'NA'}\t"
        f"{interpolation_error if interpolation_error is not None else 'NA'}\t"
        f"{'YES' if interpolation is not None else 'NO'}"
    )

    print()
    print(f"DETAIL n={n}")
    print(f"  target S = {target_s}")
    print(f"  factorized neighbors = {len(states)}")
    print()

    if nearest:
        print(
            "  CLOSEST NEIGHBOR"
        )
        print(
            f"    x  = {nearest['x']}"
        )
        print(
            f"    px = {nearest['px']}"
        )
        print(
            f"    qx = {nearest['qx']}"
        )
        print(
            f"    Sx = {nearest['sx']}"
        )
        print(
            f"    error = {nearest_error}"
        )

    print()

    if negative:
        print(
            "  CLOSEST NEGATIVE NEIGHBOR"
        )
        print(
            f"    x  = {negative['x']}"
        )
        print(
            f"    px = {negative['px']}"
        )
        print(
            f"    qx = {negative['qx']}"
        )
        print(
            f"    Sx = {negative['sx']}"
        )
        print(
            f"    error = {negative_error}"
        )

    print()

    if positive:
        print(
            "  CLOSEST POSITIVE NEIGHBOR"
        )
        print(
            f"    x  = {positive['x']}"
        )
        print(
            f"    px = {positive['px']}"
        )
        print(
            f"    qx = {positive['qx']}"
        )
        print(
            f"    Sx = {positive['sx']}"
        )
        print(
            f"    error = {positive_error}"
        )

    print()

    if interpolation is not None:
        print(
            "  TWO-SIDED INTERPOLATION"
        )
        print(
            f"    x_minus = {negative['x']}"
        )
        print(
            f"    S_minus = {negative['sx']}"
        )
        print(
            f"    x_plus  = {positive['x']}"
        )
        print(
            f"    S_plus  = {positive['sx']}"
        )
        print(
            f"    predicted S(0) = {interpolation:.10f}"
        )
        print(
            f"    target S       = {target_s}"
        )
        print(
            f"    error           = {interpolation_error:.10f}"
        )

    print()

    summary.append(
        {
            "n": n,
            "target_s": target_s,
            "nearest_error": nearest_error,
            "negative_error": negative_error,
            "positive_error": positive_error,
            "interpolation_error": interpolation_error,
        }
    )


print()
print("FINAL SUMMARY")
print("-" * 120)

print(
    "n\t"
    "S\t"
    "nearest_error\t"
    "negative_error\t"
    "positive_error\t"
    "interpolation_error"
)

for row in summary:

    interp = row["interpolation_error"]

    print(
        f"{row['n']}\t"
        f"{row['target_s']}\t"
        f"{row['nearest_error']}\t"
        f"{row['negative_error']}\t"
        f"{row['positive_error']}\t"
        f"{interp if interp is not None else 'NA'}"
    )


# Aggregate statistics

nearest_values = [
    r["nearest_error"]
    for r in summary
    if r["nearest_error"] is not None
]

negative_values = [
    r["negative_error"]
    for r in summary
    if r["negative_error"] is not None
]

positive_values = [
    r["positive_error"]
    for r in summary
    if r["positive_error"] is not None
]

interpolation_values = [
    r["interpolation_error"]
    for r in summary
    if r["interpolation_error"] is not None
]


print()
print("AGGREGATE RESULTS")

if nearest_values:
    print(
        f"  nearest-neighbor mean error = "
        f"{sum(nearest_values) / len(nearest_values):.6f}"
    )

    print(
        f"  nearest-neighbor exact = "
        f"{sum(e == 0 for e in nearest_values)} / "
        f"{len(nearest_values)}"
    )

if negative_values:
    print(
        f"  negative-side mean error = "
        f"{sum(negative_values) / len(negative_values):.6f}"
    )

    print(
        f"  negative-side exact = "
        f"{sum(e == 0 for e in negative_values)} / "
        f"{len(negative_values)}"
    )

if positive_values:
    print(
        f"  positive-side mean error = "
        f"{sum(positive_values) / len(positive_values):.6f}"
    )

    print(
        f"  positive-side exact = "
        f"{sum(e == 0 for e in positive_values)} / "
        f"{len(positive_values)}"
    )

if interpolation_values:
    print(
        f"  interpolation mean error = "
        f"{sum(interpolation_values) / len(interpolation_values):.6f}"
    )

    print(
        f"  interpolation exact = "
        f"{sum(e == 0 for e in interpolation_values)} / "
        f"{len(interpolation_values)}"
    )

    print(
        f"  interpolation <= 1 = "
        f"{sum(e <= 1 for e in interpolation_values)} / "
        f"{len(interpolation_values)}"
    )


print()
print("FINISHED EXPERIMENT 21")
