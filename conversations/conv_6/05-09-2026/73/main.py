# ================================================================
# START EXPERIMENT 156
# Information-ranked C-cell attack
#
# Experiment 155 showed:
#
#   * every C-cell gives some information;
#   * arc consistency usually performs almost no extra propagation;
#   * some C=2 cells are dramatically more restrictive than ordinary
#     C=0/C=1 cells.
#
# This experiment asks:
#
#   Can we exploit ONLY the most informative cells?
#
# For every cell (i,j), define
#
#     rho_ij = |R_ij| / (r_i*s_j)
#
# and information
#
#     I_ij = -log2(rho_ij).
#
# We rank cells by I_ij and progressively add them.
#
# We separately test:
#
#   1. strongest individual cells;
#   2. strongest cells per row;
#   3. strongest cells per column;
#   4. strongest globally ranked cells.
#
# For each selected set, we project the binary relations onto
# A_i and/or B_j and measure how much the shared-residue domains
# shrink.
#
# IMPORTANT:
#
# This remains an ORACLE experiment.
# The true C pattern is supplied from the generated factors.
#
# It is NOT a factorization algorithm yet.
#
# No sqrt(n) search.
# No p/q scan.
# No k/ell scan.
#
# ================================================================

import math
import random
import time


# ================================================================
# CONFIGURATION
# ================================================================

R1 = [
    17,
    43,
    59,
    71,
    83,
]

R2 = [
    19,
    37,
    61,
    73,
    89,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]

SAMPLES_PER_SIZE = 3


# Number of globally strongest cells to progressively add.
GLOBAL_STEPS = [
    1,
    2,
    3,
    4,
    5,
    6,
    8,
    10,
    12,
    15,
    20,
    25,
]


# ================================================================
# PRIMALITY
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in [
        2, 3, 5, 7, 11, 13, 17,
    ]:

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):

            x = (
                x * x
            ) % n

            if x == n - 1:

                composite = False
                break

        if composite:
            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p != q:
            return p, q


# ================================================================
# FACTOR DATA
# ================================================================

def factor_data(p, q):

    k = []
    a = []

    for r in R1:

        k.append(
            p // r
        )

        a.append(
            p % r
        )

    ell = []
    b = []

    for s in R2:

        ell.append(
            q // s
        )

        b.append(
            q % s
        )

    return k, a, ell, b


# ================================================================
# TRUE C MATRIX
# ================================================================

def build_true_C(
    p,
    q,
):

    k, a, ell, b = factor_data(
        p,
        q
    )

    C = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            beta = (
                (k[i] * b[j])
                % s
            )

            alpha = (
                (ell[j] * a[i])
                % r
            )

            numerator = (
                r * beta
                + s * alpha
                + a[i] * b[j]
            )

            cij = (
                numerator
                // (r * s)
            )

            row.append(cij)

        C.append(row)

    return C, a, b


# ================================================================
# LINEAR CONGRUENCE
#
# Solve:
#
#     A*x == B (mod M)
#
# Returns all x in [0,M-1].
# ================================================================

def solve_linear_congruence(
    A,
    B,
    M,
):

    A %= M
    B %= M

    g = math.gcd(
        A,
        M
    )

    if B % g != 0:
        return []

    A1 = A // g
    B1 = B // g
    M1 = M // g

    if M1 == 1:
        return list(range(M))

    inv = pow(
        A1,
        -1,
        M1
    )

    x0 = (
        B1 * inv
    ) % M1

    return [
        (
            x0
            + t * M1
        ) % M
        for t in range(g)
    ]


# ================================================================
# BUILD ONE CELL RELATION
#
# A relation contains pairs (a,b) that are compatible with:
#
#   n mod (r*s)
#
# and the observed C value.
# ================================================================

def build_cell_relation(
    n,
    r,
    s,
    observed_c,
):

    relation = set()

    n_mod = (
        n % (r * s)
    )

    for a in range(r):

        for b in range(s):

            target = (
                n_mod
                - a * b
            ) % (r * s)

            # ----------------------------------------------------
            # x = k mod s
            #
            # r*b*x == target (mod s)
            # ----------------------------------------------------

            xs = solve_linear_congruence(
                (r * b) % s,
                target % s,
                s,
            )

            if not xs:
                continue

            # ----------------------------------------------------
            # y = ell mod r
            #
            # s*a*y == target (mod r)
            # ----------------------------------------------------

            ys = solve_linear_congruence(
                (s * a) % r,
                target % r,
                r,
            )

            if not ys:
                continue

            compatible = False

            for x in xs:

                beta = (
                    x * b
                ) % s

                for y in ys:

                    alpha = (
                        y * a
                    ) % r

                    numerator = (
                        r * beta
                        + s * alpha
                        + a * b
                    )

                    cij = (
                        numerator
                        // (r * s)
                    )

                    if cij == observed_c:

                        compatible = True
                        break

                if compatible:
                    break

            if compatible:

                relation.add(
                    (a, b)
                )

    return relation


# ================================================================
# INFORMATION METRIC
# ================================================================

def information_bits(
    relation_size,
    total_size,
):

    if relation_size <= 0:
        return float("inf")

    return math.log2(
        total_size
        / relation_size
    )


# ================================================================
# RELATION PROJECTIONS
# ================================================================

def relation_a_projection(
    relation
):

    return {
        a
        for a, b in relation
    }


def relation_b_projection(
    relation
):

    return {
        b
        for a, b in relation
    }


# ================================================================
# DOMAIN PRODUCT
# ================================================================

def product_size(
    domains
):

    x = 1

    for d in domains:

        x *= len(d)

    return x


# ================================================================
# SELECTED-CELL DOMAIN TEST
#
# A cell constrains both variables:
#
#   A_i
#   B_j
#
# We repeatedly intersect each endpoint domain with the projection
# of every selected relation touching that variable.
#
# This is deliberately simpler than full arc consistency.
#
# ================================================================

def projected_domains(
    selected_cells,
    relations,
):

    domains_a = [
        set(range(r))
        for r in R1
    ]

    domains_b = [
        set(range(s))
        for s in R2
    ]

    changed = True
    iterations = 0

    while changed:

        changed = False
        iterations += 1

        for i, j in selected_cells:

            relation = relations[i][j]

            pa = relation_a_projection(
                relation
            )

            pb = relation_b_projection(
                relation
            )

            new_a = (
                domains_a[i]
                & pa
            )

            new_b = (
                domains_b[j]
                & pb
            )

            if new_a != domains_a[i]:

                domains_a[i] = new_a
                changed = True

            if new_b != domains_b[j]:

                domains_b[j] = new_b
                changed = True

            # ----------------------------------------------------
            # Now filter the relation itself according to the
            # current endpoint domains.
            #
            # This creates genuine propagation between A_i and B_j.
            # ----------------------------------------------------

            filtered_relation = {
                (a, b)
                for a, b in relation
                if (
                    a in domains_a[i]
                    and b in domains_b[j]
                )
            }

            if filtered_relation != relation:

                relations[i][j] = (
                    filtered_relation
                )

                new_a = (
                    domains_a[i]
                    & relation_a_projection(
                        filtered_relation
                    )
                )

                new_b = (
                    domains_b[j]
                    & relation_b_projection(
                        filtered_relation
                    )
                )

                if new_a != domains_a[i]:

                    domains_a[i] = new_a
                    changed = True

                if new_b != domains_b[j]:

                    domains_b[j] = new_b
                    changed = True

    return (
        domains_a,
        domains_b,
        iterations,
    )


# ================================================================
# NOTE:
#
# projected_domains mutates its relation dictionary.
#
# We therefore clone only the selected relation sets before each
# independent test.
# ================================================================

def clone_relations(
    relations
):

    return [
        [
            set(relations[i][j])
            for j in range(len(R2))
        ]
        for i in range(len(R1))
    ]


# ================================================================
# CELL RANKING
# ================================================================

def rank_cells(
    relations
):

    cells = []

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            size = len(
                relations[i][j]
            )

            raw = r * s

            info = information_bits(
                size,
                raw
            )

            cells.append({
                "i": i,
                "j": j,
                "r": r,
                "s": s,
                "c": None,
                "size": size,
                "raw": raw,
                "fraction": (
                    size / raw
                ),
                "info": info,
            })

    cells.sort(
        key=lambda x: (
            -x["info"],
            x["i"],
            x["j"],
        )
    )

    return cells


# ================================================================
# MAIN
# ================================================================

print("=" * 72)
print("START EXPERIMENT 156")
print("Information-ranked C-cell attack")
print("=" * 72)
print()

print("R1 =", R1)
print("R2 =", R2)
print(
    "SAMPLES PER BIT SIZE =",
    SAMPLES_PER_SIZE
)
print()

aggregate = {
    "samples": 0,
    "true_survived": 0,
    "strongest_1_true": 0,
    "strongest_5_true": 0,
    "strongest_10_true": 0,
}


for bits in BIT_SIZES:

    print("=" * 72)
    print(
        f"BIT SIZE = {bits}"
    )
    print("=" * 72)
    print()

    for sample in range(
        1,
        SAMPLES_PER_SIZE + 1
    ):

        t0 = time.perf_counter()

        p, q = random_semiprime(
            bits
        )

        n = p * q

        C, true_a, true_b = build_true_C(
            p,
            q
        )

        print("-" * 72)

        print(
            f"SAMPLE {sample}/{SAMPLES_PER_SIZE}"
        )

        print(
            "n bits =",
            n.bit_length()
        )

        print(
            "n      =",
            n
        )

        print(
            "p      =",
            p
        )

        print(
            "q      =",
            q
        )

        print()

        print(
            "TRUE a =",
            true_a
        )

        print(
            "TRUE b =",
            true_b
        )

        print()

        print(
            "C MATRIX"
        )

        for row in C:

            print(row)

        print()

        # --------------------------------------------------------
        # Build relations.
        # --------------------------------------------------------

        relations = [
            [None for _ in R2]
            for _ in R1
        ]

        cells = []

        for i, r in enumerate(R1):

            for j, s in enumerate(R2):

                relation = build_cell_relation(
                    n,
                    r,
                    s,
                    C[i][j],
                )

                relations[i][j] = relation

                raw = r * s

                info = information_bits(
                    len(relation),
                    raw
                )

                cells.append({
                    "i": i,
                    "j": j,
                    "r": r,
                    "s": s,
                    "c": C[i][j],
                    "size": len(relation),
                    "raw": raw,
                    "fraction": len(relation) / raw,
                    "info": info,
                })

        cells.sort(
            key=lambda x: (
                -x["info"],
                x["i"],
                x["j"],
            )
        )

        aggregate[
            "samples"
        ] += 1

        # --------------------------------------------------------
        # Print ranked cells.
        # --------------------------------------------------------

        print(
            "GLOBAL CELL INFORMATION RANK"
        )

        print(
            "rank  cell      C   states/raw"
            "            fraction"
            "          information"
        )

        print(
            "-" * 72
        )

        for rank, cell in enumerate(
            cells,
            start=1
        ):

            print(
                f"{rank:4d}  "
                f"({cell['i']},{cell['j']})"
                f"       "
                f"{cell['c']}   "
                f"{cell['size']:5d}/"
                f"{cell['raw']:5d}      "
                f"{cell['fraction']:.6f}      "
                f"{cell['info']:.6f} bits"
            )

        print()

        # --------------------------------------------------------
        # Check true pair survives every cell.
        # --------------------------------------------------------

        true_all = True

        for cell in cells:

            i = cell["i"]
            j = cell["j"]

            if (
                true_a[i],
                true_b[j]
            ) not in relations[i][j]:

                true_all = False

        if true_all:

            aggregate[
                "true_survived"
            ] += 1

        print(
            "TRUE PAIR IN ALL CELL RELATIONS:",
            true_all
        )

        print()

        # --------------------------------------------------------
        # GLOBAL PROGRESSIVE TEST
        # --------------------------------------------------------

        print(
            "GLOBAL TOP-K PROPAGATION"
        )

        print(
            "K    selected"
            "    A-product"
            "          B-product"
            "          A-frac"
            "          B-frac"
            "          true-A"
            " true-B"
        )

        print(
            "-" * 72
        )

        for K in GLOBAL_STEPS:

            if K > len(cells):
                continue

            selected = [
                (
                    cells[t]["i"],
                    cells[t]["j"]
                )
                for t in range(K)
            ]

            test_relations = clone_relations(
                relations
            )

            domains_a, domains_b, iterations = (
                projected_domains(
                    selected,
                    test_relations
                )
            )

            raw_a = math.prod(
                R1
            )

            raw_b = math.prod(
                R2
            )

            final_a = product_size(
                domains_a
            )

            final_b = product_size(
                domains_b
            )

            true_a_ok = all(
                true_a[i]
                in domains_a[i]
                for i in range(len(R1))
            )

            true_b_ok = all(
                true_b[j]
                in domains_b[j]
                for j in range(len(R2))
            )

            print(
                f"{K:2d}   "
                f"{K:8d}   "
                f"{final_a:12d}   "
                f"{final_b:12d}   "
                f"{final_a/raw_a:.6e}   "
                f"{final_b/raw_b:.6e}   "
                f"{str(true_a_ok):>5}   "
                f"{str(true_b_ok):>5}"
            )

            if K == 1 and \
                    true_a_ok and \
                    true_b_ok:

                aggregate[
                    "strongest_1_true"
                ] += 1

            if K == 5 and \
                    true_a_ok and \
                    true_b_ok:

                aggregate[
                    "strongest_5_true"
                ] += 1

            if K == 10 and \
                    true_a_ok and \
                    true_b_ok:

                aggregate[
                    "strongest_10_true"
                ] += 1

        print()

        # --------------------------------------------------------
        # PER-ROW STRONGEST CELLS
        #
        # Select the single most informative cell touching each row.
        # --------------------------------------------------------

        row_best = {}

        for cell in cells:

            i = cell["i"]

            if i not in row_best:

                row_best[i] = cell

        row_selected = [
            (
                row_best[i]["i"],
                row_best[i]["j"]
            )
            for i in range(len(R1))
        ]

        test_relations = clone_relations(
            relations
        )

        row_domains_a, row_domains_b, _ = (
            projected_domains(
                row_selected,
                test_relations
            )
        )

        print(
            "BEST CELL PER ROW"
        )

        print(
            "selected cells =",
            row_selected
        )

        print(
            "A domains:"
        )

        for i in range(len(R1)):

            print(
                f"  a[{i}] = "
                f"{len(row_domains_a[i])}/"
                f"{R1[i]} "
                f"true="
                f"{true_a[i] in row_domains_a[i]}"
            )

        print()

        print(
            "B domains:"
        )

        for j in range(len(R2)):

            print(
                f"  b[{j}] = "
                f"{len(row_domains_b[j])}/"
                f"{R2[j]} "
                f"true="
                f"{true_b[j] in row_domains_b[j]}"
            )

        print()

        # --------------------------------------------------------
        # PER-COLUMN STRONGEST CELLS
        # --------------------------------------------------------

        col_best = {}

        for cell in cells:

            j = cell["j"]

            if j not in col_best:

                col_best[j] = cell

        col_selected = [
            (
                col_best[j]["i"],
                col_best[j]["j"]
            )
            for j in range(len(R2))
        ]

        test_relations = clone_relations(
            relations
        )

        col_domains_a, col_domains_b, _ = (
            projected_domains(
                col_selected,
                test_relations
            )
        )

        print(
            "BEST CELL PER COLUMN"
        )

        print(
            "selected cells =",
            col_selected
        )

        print(
            "A domains:"
        )

        for i in range(len(R1)):

            print(
                f"  a[{i}] = "
                f"{len(col_domains_a[i])}/"
                f"{R1[i]} "
                f"true="
                f"{true_a[i] in col_domains_a[i]}"
            )

        print()

        print(
            "B domains:"
        )

        for j in range(len(R2)):

            print(
                f"  b[{j}] = "
                f"{len(col_domains_b[j])}/"
                f"{R2[j]} "
                f"true="
                f"{true_b[j] in col_domains_b[j]}"
            )

        print()

        # --------------------------------------------------------
        # TWO-DIRECTION UNION:
        #
        # strongest cell for every row PLUS strongest cell for
        # every column.
        #
        # Duplicate cells are removed.
        # --------------------------------------------------------

        union_selected = list(
            dict.fromkeys(
                row_selected
                + col_selected
            )
        )

        test_relations = clone_relations(
            relations
        )

        union_a, union_b, union_iterations = (
            projected_domains(
                union_selected,
                test_relations
            )
        )

        print(
            "BEST ROW + BEST COLUMN UNION"
        )

        print(
            "selected cells =",
            union_selected
        )

        print(
            "iterations =",
            union_iterations
        )

        print()

        print(
            "A domains:"
        )

        for i in range(len(R1)):

            print(
                f"  a[{i}] = "
                f"{len(union_a[i])}/"
                f"{R1[i]} "
                f"true="
                f"{true_a[i] in union_a[i]}"
            )

        print()

        print(
            "B domains:"
        )

        for j in range(len(R2)):

            print(
                f"  b[{j}] = "
                f"{len(union_b[j])}/"
                f"{R2[j]} "
                f"true="
                f"{true_b[j] in union_b[j]}"
            )

        print()

        # --------------------------------------------------------
        # CHECK THE STRONGEST CELL.
        # --------------------------------------------------------

        strongest = cells[0]

        print(
            "STRONGEST SINGLE CELL"
        )

        print(
            f"cell = ({strongest['i']},"
            f"{strongest['j']})"
        )

        print(
            "C =",
            strongest["c"]
        )

        print(
            "states =",
            strongest["size"]
        )

        print(
            "raw =",
            strongest["raw"]
        )

        print(
            "fraction =",
            f"{strongest['fraction']:.12e}"
        )

        print(
            "information =",
            f"{strongest['info']:.12e}",
            "bits"
        )

        print()

        # --------------------------------------------------------
        # C=0 / C=1 / C=2 information summary.
        # --------------------------------------------------------

        buckets = {
            0: [],
            1: [],
            2: [],
        }

        for cell in cells:

            buckets[
                cell["c"]
            ].append(
                cell["info"]
            )

        print(
            "INFORMATION BY C VALUE"
        )

        for c in [0, 1, 2]:

            values = buckets[c]

            if not values:

                print(
                    f"C={c}: no cells"
                )

                continue

            print(
                f"C={c}: "
                f"count={len(values)}, "
                f"min={min(values):.6f}, "
                f"mean="
                f"{sum(values)/len(values):.6f}, "
                f"max={max(values):.6f}"
            )

        print()

        # --------------------------------------------------------
        # Runtime
        # --------------------------------------------------------

        elapsed = (
            time.perf_counter()
            - t0
        )

        print(
            "SAMPLE RUNTIME =",
            f"{elapsed:.6f} s"
        )

        print()


# ================================================================
# AGGREGATE
# ================================================================

print("=" * 72)
print("EXPERIMENT 156 AGGREGATE")
print("=" * 72)
print()

print(
    "samples =",
    aggregate["samples"]
)

print(
    "true pair survived all relations =",
    aggregate["true_survived"],
    "/",
    aggregate["samples"]
)

print()

print(
    "top-1 test true survival =",
    aggregate["strongest_1_true"],
    "/",
    aggregate["samples"]
)

print(
    "top-5 test true survival =",
    aggregate["strongest_5_true"],
    "/",
    aggregate["samples"]
)

print(
    "top-10 test true survival =",
    aggregate["strongest_10_true"],
    "/",
    aggregate["samples"]
)

print()

print(
    "The key quantities are:"
)

print(
    "  information(cell)"
)

print(
    "  final A/B domain fractions"
)

print(
    "  whether the strongest cells produce singleton residues"
)

print(
    "  whether adding cells beats ordinary arc consistency"
)

print()

print("=" * 72)
print("FINISHED EXPERIMENT 156")
print("=" * 72)
