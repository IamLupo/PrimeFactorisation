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
    """
    Return all prime factor pairs (p,q), p <= q,
    for n.
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
# EXPERIMENT
# ============================================================

print("START EXPERIMENT 16")
print()

print(
    "n\tp\tq\t"
    "x\tpx\tqx\t"
    "dp=dpx\t"
    "dq=dqx\t"
    "linear_term\t"
    "cross_term\t"
    "reconstructed_x\t"
    "exact"
)


summary = []


for p, q in TEST_CASES:

    n = p * q

    rows = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = prime_factor_pairs(nx)

        for px, qx in pairs:

            #
            # Displacement of each factor.
            #
            dp = px - p
            dq = qx - q

            #
            # Exact identity:
            #
            # x =
            #     p*dq
            #   + q*dp
            #   + dp*dq
            #
            linear_term = (
                p * dq
                + q * dp
            )

            cross_term = (
                dp * dq
            )

            reconstructed_x = (
                linear_term
                + cross_term
            )

            exact = (
                reconstructed_x == x
            )

            rows.append(
                (
                    x,
                    px,
                    qx,
                    dp,
                    dq,
                    linear_term,
                    cross_term,
                    reconstructed_x,
                    exact,
                )
            )

            print(
                f"{n}\t"
                f"{p}\t"
                f"{q}\t"
                f"{x}\t"
                f"{px}\t"
                f"{qx}\t"
                f"{dp}\t"
                f"{dq}\t"
                f"{linear_term}\t"
                f"{cross_term}\t"
                f"{reconstructed_x}\t"
                f"{exact}"
            )

    #
    # ========================================================
    # SUMMARY STATISTICS
    # ========================================================
    #

    if not rows:
        continue

    exact_count = sum(
        row[-1]
        for row in rows
    )

    #
    # Cases where the factor displacement is small.
    #
    small_displacement = [
        row for row in rows
        if abs(row[3]) <= 5
        and abs(row[4]) <= 5
    ]

    #
    # Cases where both factors move in the same direction.
    #
    same_direction = [
        row for row in rows
        if row[3] != 0
        and row[4] != 0
        and (
            (row[3] > 0 and row[4] > 0)
            or
            (row[3] < 0 and row[4] < 0)
        )
    ]

    #
    # Cases where the displacement is opposite.
    #
    opposite_direction = [
        row for row in rows
        if row[3] != 0
        and row[4] != 0
        and (
            (row[3] > 0 and row[4] < 0)
            or
            (row[3] < 0 and row[4] > 0)
        )
    ]

    #
    # Track gcd(x,n).
    #
    gcd_rows = []

    for row in rows:

        x = row[0]

        a = abs(x)
        b = n

        while b:
            a, b = b, a % b

        gcd_rows.append(
            (x, a)
        )

    nontrivial_gcd = [
        item for item in gcd_rows
        if item[1] > 1
    ]

    summary.append(
        (
            n,
            p,
            q,
            len(rows),
            exact_count,
            len(small_displacement),
            len(same_direction),
            len(opposite_direction),
            len(nontrivial_gcd),
        )
    )


print()
print("SUMMARY")
print("-" * 120)

print(
    "n\tp\tq\t"
    "neighbor_states\t"
    "identity_exact\t"
    "small_displacement\t"
    "same_direction\t"
    "opposite_direction\t"
    "gcd(x,n)>1"
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
        f"{row[7]}\t"
        f"{row[8]}"
    )


print()
print("SMALL-DISPLACEMENT DETAILS")
print("-" * 120)

for p, q in TEST_CASES:

    n = p * q

    found = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        for px, qx in prime_factor_pairs(nx):

            dp = px - p
            dq = qx - q

            if abs(dp) <= 5 and abs(dq) <= 5:

                found.append(
                    (
                        x,
                        px,
                        qx,
                        dp,
                        dq,
                        dp * dq,
                    )
                )

    print()
    print(
        f"n={n}, target=(p={p},q={q})"
    )

    if not found:
        print("  None")
        continue

    for row in found:
        print(
            "  x=%4d px=%4d qx=%4d "
            "dp=%4d dq=%4d dp*dq=%5d"
            % row
        )


print()
print("GCD CHECK")
print("-" * 120)

for p, q in TEST_CASES:

    n = p * q

    useful = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        a = abs(x)
        b = n

        while b:
            a, b = b, a % b

        if a > 1:
            useful.append(
                (x, a)
            )

    print(
        f"n={n}: "
        f"nontrivial gcd offsets="
        f"{useful}"
    )


print()
print("FINISHED EXPERIMENT 16")
