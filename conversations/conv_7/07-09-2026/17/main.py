from math import isqrt


R1 = 11
R2 = 13

WINDOW = 50

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


def factor_pairs(n):
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


# ============================================================
# MODULAR COORDINATES
# ============================================================

def coordinate(x, r):
    return x % r, x // r


def state(p, q):
    """
    Full modular state of p,q.
    """

    a, k = coordinate(p, R1)
    b, l = coordinate(q, R2)

    return {
        "p": p,
        "q": q,

        "a": a,
        "b": b,

        "k": k,
        "l": l,

        "K": k * l,

        "S": p + q,
    }


# ============================================================
# TARGET ANALYSIS
# ============================================================

print("START EXPERIMENT 14")
print()

print(
    "n\tp\tq\t"
    "a\tk\tb\tl\tK\t"
    "neighbors\t"
    "same_p_cell\t"
    "same_q_cell\t"
    "same_both_cells\t"
    "exact_p\t"
    "exact_q\t"
    "exact_both\t"
    "same_p_cell_pct\t"
    "same_q_cell_pct\t"
    "same_both_pct"
)

for p, q in TEST_CASES:

    n = p * q

    target = state(p, q)

    target_p_cell = (
        target["a"],
        target["k"],
    )

    target_q_cell = (
        target["b"],
        target["l"],
    )

    neighbors = 0

    same_p_cell = 0
    same_q_cell = 0
    same_both_cells = 0

    exact_p = 0
    exact_q = 0
    exact_both = 0

    #
    # Also retain the actual trajectories so that
    # we can print them later.
    #
    trajectory = []

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = factor_pairs(nx)

        for px, qx in pairs:

            neighbors += 1

            nearby = state(px, qx)

            p_cell = (
                nearby["a"],
                nearby["k"],
            )

            q_cell = (
                nearby["b"],
                nearby["l"],
            )

            if p_cell == target_p_cell:
                same_p_cell += 1

            if q_cell == target_q_cell:
                same_q_cell += 1

            if (
                p_cell == target_p_cell
                and q_cell == target_q_cell
            ):
                same_both_cells += 1

            if px == p:
                exact_p += 1

            if qx == q:
                exact_q += 1

            if (
                (px == p and qx == q)
                or
                (px == q and qx == p)
            ):
                exact_both += 1

            trajectory.append(
                (
                    x,
                    px,
                    qx,
                    nearby["a"],
                    nearby["k"],
                    nearby["b"],
                    nearby["l"],
                    nearby["K"],
                )
            )

    p_pct = (
        100.0 * same_p_cell / neighbors
        if neighbors else 0.0
    )

    q_pct = (
        100.0 * same_q_cell / neighbors
        if neighbors else 0.0
    )

    both_pct = (
        100.0 * same_both_cells / neighbors
        if neighbors else 0.0
    )

    print(
        f"{n}\t"
        f"{p}\t"
        f"{q}\t"
        f"{target['a']}\t"
        f"{target['k']}\t"
        f"{target['b']}\t"
        f"{target['l']}\t"
        f"{target['K']}\t"
        f"{neighbors}\t"
        f"{same_p_cell}\t"
        f"{same_q_cell}\t"
        f"{same_both_cells}\t"
        f"{exact_p}\t"
        f"{exact_q}\t"
        f"{exact_both}\t"
        f"{p_pct:.4f}\t"
        f"{q_pct:.4f}\t"
        f"{both_pct:.4f}"
    )


print()
print("COORDINATE TRAJECTORIES")
print("-" * 120)

for p, q in TEST_CASES:

    n = p * q

    target = state(p, q)

    print()
    print(
        f"n={n} "
        f"target=("
        f"p={p}, q={q}; "
        f"a={target['a']}, k={target['k']}; "
        f"b={target['b']}, l={target['l']}; "
        f"K={target['K']}"
        f")"
    )

    print(
        "x\tpx\tqx\ta_x\tk_x\tb_x\tl_x\tK_x"
    )

    for x in range(-WINDOW, WINDOW + 1):

        if x == 0:
            continue

        nx = n + x

        if nx <= 1:
            continue

        pairs = factor_pairs(nx)

        for px, qx in pairs:

            s = state(px, qx)

            print(
                f"{x}\t"
                f"{px}\t"
                f"{qx}\t"
                f"{s['a']}\t"
                f"{s['k']}\t"
                f"{s['b']}\t"
                f"{s['l']}\t"
                f"{s['K']}"
            )


print()
print("FINISHED EXPERIMENT 14")