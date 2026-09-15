# ================================================================
# START EXPERIMENT 158
# Global radix-residue CSP
#
# Global variables:
#
#     X_m = p mod m
#
# for every radix m.
#
# For each cell (r_i, s_j), the exact carry C_ij becomes a binary
# constraint on:
#
#     X_{r_i}, X_{s_j}
#
# because the pair determines p mod (r_i*s_j), while
#
#     q mod m = n * X_m^{-1} mod m.
#
# This experiment measures:
#
#   * binary relation sizes
#   * AC-3 propagation
#   * remaining global Cartesian state count
#   * exact CSP solution count
#   * whether the true residue vector is recovered
#
# There is no sqrt(n) search and no p/q scan.
#
# ================================================================

import math
import random
import time


# ================================================================
# CONFIG
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

ALL_RADICES = R1 + R2

RADIX_INDEX = {
    m: i
    for i, m in enumerate(ALL_RADICES)
}

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]

SAMPLES_PER_SIZE = 3

MAX_SOLUTIONS = 1_000_000


# ================================================================
# PRIMALITY
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for p in small_primes:

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
        2, 3, 5, 7, 11, 13, 17
    ]:

        if a >= n:
            continue

        x = pow(
            a,
            d,
            n
        )

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
        x |= (
            1 << (bits - 1)
        )

        if is_probable_prime(x):
            return x


def random_semiprime(bits):

    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:

        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p != q:
            return p, q


# ================================================================
# CRT FOR COPRIME MODULI
# ================================================================

def crt2(
    a,
    m,
    b,
    n,
):

    # Solve:
    #
    #     x = a (mod m)
    #     x = b (mod n)
    #
    # with gcd(m,n)=1.

    if math.gcd(m, n) != 1:

        raise ValueError(
            f"CRT moduli not coprime: {m}, {n}"
        )

    inv = pow(
        m % n,
        -1,
        n
    )

    t = (
        (b - a)
        * inv
    ) % n

    return (
        a + m * t
    ) % (
        m * n
    )


# ================================================================
# TRUE C MATRIX
#
# THIS IS THE FUNCTION NAME USED THROUGHOUT THE SCRIPT.
# ================================================================

def true_C_matrix(
    p,
    q,
):

    k = [
        p // r
        for r in R1
    ]

    a = [
        p % r
        for r in R1
    ]

    ell = [
        q // s
        for s in R2
    ]

    b = [
        q % s
        for s in R2
    ]

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
                // (
                    r * s
                )
            )

            row.append(cij)

        C.append(row)

    return C


# ================================================================
# TRUE GLOBAL P-RESIDUE VECTOR
# ================================================================

def true_residue_vector(p):

    return tuple(
        p % m
        for m in ALL_RADICES
    )


# ================================================================
# q mod m
#
# Since p*q = n mod m and p mod m is nonzero:
#
#     q = n * p^{-1} mod m
# ================================================================

def q_residue(
    n,
    x,
    m,
):

    if math.gcd(
        x,
        m
    ) != 1:

        return None

    return (
        n
        * pow(
            x,
            -1,
            m
        )
    ) % m


# ================================================================
# CELL C FROM GLOBAL P RESIDUES
#
# Given:
#
#     xr = p mod r
#     xs = p mod s
#
# recover:
#
#     p mod rs
#     k mod s
#
# and:
#
#     q mod r
#     q mod s
#     ell mod r
#
# then calculate the exact C.
# ================================================================

def cell_c_from_global_residues(
    n,
    r,
    s,
    xr,
    xs,
):

    # ------------------------------------------------------------
    # p side
    # ------------------------------------------------------------

    a = xr

    p_rs = crt2(
        xr,
        r,
        xs,
        s
    )

    k_mod_s = (
        (
            p_rs
            - a
        )
        // r
    ) % s

    # ------------------------------------------------------------
    # q side
    # ------------------------------------------------------------

    yr = q_residue(
        n,
        xr,
        r
    )

    ys = q_residue(
        n,
        xs,
        s
    )

    if yr is None or ys is None:
        return None

    b = ys

    q_rs = crt2(
        yr,
        r,
        ys,
        s
    )

    ell_mod_r = (
        (
            q_rs
            - b
        )
        // s
    ) % r

    # ------------------------------------------------------------
    # final carry
    # ------------------------------------------------------------

    beta = (
        k_mod_s
        * b
    ) % s

    alpha = (
        ell_mod_r
        * a
    ) % r

    numerator = (
        r * beta
        + s * alpha
        + a * b
    )

    cij = (
        numerator
        // (
            r * s
        )
    )

    return cij


# ================================================================
# BUILD ONE BINARY CONSTRAINT
#
# relation = {(xr,xs) : C(xr,xs) == observed_c}
# ================================================================

def build_constraint(
    n,
    r,
    s,
    observed_c,
):

    relation = set()

    for xr in range(
        1,
        r
    ):

        for xs in range(
            1,
            s
        ):

            cij = cell_c_from_global_residues(
                n,
                r,
                s,
                xr,
                xs
            )

            if cij == observed_c:

                relation.add(
                    (
                        xr,
                        xs
                    )
                )

    return relation


# ================================================================
# GET RELATION IN xi -> xj ORIENTATION
# ================================================================

def get_relation(
    constraints,
    xi,
    xj,
):

    if (
        xi,
        xj
    ) in constraints:

        return constraints[
            (xi, xj)
        ]

    if (
        xj,
        xi
    ) not in constraints:

        raise KeyError(
            f"Missing constraint for "
            f"{xi},{xj}"
        )

    return {
        (
            b,
            a
        )
        for a, b in constraints[
            (xj, xi)
        ]
    }


# ================================================================
# DOMAIN PRODUCT
# ================================================================

def domain_product(
    domains
):

    result = 1

    for d in domains:
        result *= len(d)

    return result


# ================================================================
# BUILD NEIGHBORS
# ================================================================

def build_neighbors(
    constraints
):

    neighbors = {
        i: set()
        for i in range(
            len(ALL_RADICES)
        )
    }

    for i, j in constraints:

        neighbors[i].add(j)
        neighbors[j].add(i)

    return neighbors


# ================================================================
# AC-3 REVISION
# ================================================================

def revise(
    domains,
    constraints,
    xi,
    xj,
):

    relation = get_relation(
        constraints,
        xi,
        xj
    )

    old_domain = domains[xi]

    new_domain = set()

    for value_i in old_domain:

        supported = False

        for a, b in relation:

            if (
                a == value_i
                and
                b in domains[xj]
            ):

                supported = True
                break

        if supported:

            new_domain.add(
                value_i
            )

    if new_domain == old_domain:

        return False

    domains[xi] = new_domain

    return True


# ================================================================
# AC-3
# ================================================================

def ac3(
    domains,
    constraints
):

    neighbors = build_neighbors(
        constraints
    )

    queue = []

    for i, j in constraints:

        queue.append(
            (
                i,
                j
            )
        )

        queue.append(
            (
                j,
                i
            )
        )

    revisions = 0

    while queue:

        xi, xj = queue.pop(0)

        revisions += 1

        changed = revise(
            domains,
            constraints,
            xi,
            xj
        )

        if not changed:
            continue

        if not domains[xi]:

            return revisions

        for xk in neighbors[xi]:

            if xk == xj:
                continue

            queue.append(
                (
                    xk,
                    xi
                )
            )

    return revisions


# ================================================================
# CHOOSE DFS VARIABLE
# ================================================================

def choose_variable(
    domains,
    assigned,
    neighbors
):

    candidates = []

    for i in range(
        len(ALL_RADICES)
    ):

        if i in assigned:
            continue

        candidates.append(
            (
                len(domains[i]),
                -len(neighbors[i]),
                i
            )
        )

    return min(
        candidates
    )[2]


# ================================================================
# TEST VALUE AGAINST ASSIGNED NEIGHBORS
# ================================================================

def value_consistent(
    var,
    value,
    assigned,
    constraints,
    neighbors
):

    for nb in neighbors[var]:

        if nb not in assigned:
            continue

        if var < nb:

            relation = constraints[
                (
                    var,
                    nb
                )
            ]

            pair = (
                value,
                assigned[nb]
            )

        else:

            relation = constraints[
                (
                    nb,
                    var
                )
            ]

            pair = (
                assigned[nb],
                value
            )

        if pair not in relation:

            return False

    return True


# ================================================================
# DFS ENUMERATION
# ================================================================

def enumerate_solutions(
    domains,
    constraints,
    max_solutions
):

    working_domains = [
        set(d)
        for d in domains
    ]

    neighbors = build_neighbors(
        constraints
    )

    solutions = []

    nodes = 0
    limit_hit = False

    def dfs(
        assigned
    ):

        nonlocal nodes
        nonlocal limit_hit

        if len(solutions) >= max_solutions:

            limit_hit = True
            return

        nodes += 1

        # --------------------------------------------------------
        # Complete assignment
        # --------------------------------------------------------

        if len(assigned) == len(
            ALL_RADICES
        ):

            solution = tuple(
                assigned[i]
                for i in range(
                    len(ALL_RADICES)
                )
            )

            solutions.append(
                solution
            )

            return

        # --------------------------------------------------------
        # Pick variable
        # --------------------------------------------------------

        var = choose_variable(
            working_domains,
            assigned,
            neighbors
        )

        # --------------------------------------------------------
        # Branch
        # --------------------------------------------------------

        for value in sorted(
            working_domains[var]
        ):

            if not value_consistent(
                var,
                value,
                assigned,
                constraints,
                neighbors
            ):
                continue

            assigned[var] = value

            backups = {}

            valid = True

            # ----------------------------------------------------
            # Forward-check neighbors.
            # ----------------------------------------------------

            for nb in neighbors[var]:

                if nb in assigned:
                    continue

                if nb not in backups:

                    backups[nb] = set(
                        working_domains[nb]
                    )

                if var < nb:

                    relation = constraints[
                        (
                            var,
                            nb
                        )
                    ]

                    allowed = {
                        b
                        for a, b in relation
                        if (
                            a == value
                            and
                            b in working_domains[nb]
                        )
                    }

                else:

                    relation = constraints[
                        (
                            nb,
                            var
                        )
                    ]

                    allowed = {
                        a
                        for a, b in relation
                        if (
                            b == value
                            and
                            a in working_domains[nb]
                        )
                    }

                working_domains[nb] &= allowed

                if not working_domains[nb]:

                    valid = False
                    break

            if valid:

                dfs(
                    assigned
                )

            # ----------------------------------------------------
            # Restore
            # ----------------------------------------------------

            for nb, backup in backups.items():

                working_domains[nb] = backup

            del assigned[var]

            if limit_hit:
                return

    dfs({})

    return (
        solutions,
        nodes,
        limit_hit
    )


# ================================================================
# PRINT DOMAINS
# ================================================================

def print_domains(
    title,
    domains
):

    print(title)

    for i, m in enumerate(
        ALL_RADICES
    ):

        values = sorted(
            domains[i]
        )

        if len(values) <= 20:

            shown = values

        else:

            shown = (
                values[:10]
                + ["..."]
                + values[-5:]
            )

        print(
            f"  m={m:2d}: "
            f"{len(values):3d}/"
            f"{m - 1:3d} "
            f"{shown}"
        )

    print()


# ================================================================
# RECONSTRUCT p MOD FULL PRODUCT
# ================================================================

def reconstruct_full_residue(
    solution
):

    x = solution[0]
    M = ALL_RADICES[0]

    for i in range(
        1,
        len(ALL_RADICES)
    ):

        x = crt2(
            x,
            M,
            solution[i],
            ALL_RADICES[i]
        )

        M *= ALL_RADICES[i]

    return (
        x,
        M
    )


# ================================================================
# MAIN EXPERIMENT
# ================================================================

print("=" * 72)
print("START EXPERIMENT 158")
print("Global radix-residue CSP")
print("=" * 72)
print()

print(
    "R1 =",
    R1
)

print(
    "R2 =",
    R2
)

print(
    "ALL RADICES =",
    ALL_RADICES
)

print()

aggregate_samples = 0
aggregate_true_survival = 0
aggregate_unique = 0
aggregate_limited = 0


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

        C = true_C_matrix(
            p,
            q
        )

        true_x = true_residue_vector(
            p
        )

        aggregate_samples += 1

        print("-" * 72)

        print(
            f"SAMPLE "
            f"{sample}/{SAMPLES_PER_SIZE}"
        )

        print(
            "n bits =",
            n.bit_length()
        )

        print(
            "n =",
            n
        )

        print(
            "true p =",
            p
        )

        print(
            "true q =",
            q
        )

        print()

        print(
            "TRUE p RESIDUES"
        )

        for i, m in enumerate(
            ALL_RADICES
        ):

            print(
                f"  p mod {m:2d} =",
                true_x[i]
            )

        print()

        print(
            "C MATRIX"
        )

        for row in C:
            print(row)

        print()

        # --------------------------------------------------------
        # Build constraints
        # --------------------------------------------------------

        constraints = {}

        relation_rows = []

        print(
            "BUILDING CONSTRAINTS"
        )

        for i, r in enumerate(R1):

            for j, s in enumerate(R2):

                xi = RADIX_INDEX[r]
                xj = RADIX_INDEX[s]

                relation = build_constraint(
                    n,
                    r,
                    s,
                    C[i][j]
                )

                constraints[
                    (xi, xj)
                ] = relation

                raw = (
                    (r - 1)
                    * (s - 1)
                )

                count = len(
                    relation
                )

                if count == 0:

                    info = float("inf")

                else:

                    info = math.log2(
                        raw / count
                    )

                true_pair = (
                    true_x[xi],
                    true_x[xj]
                )

                true_ok = (
                    true_pair
                    in relation
                )

                relation_rows.append(
                    {
                        "r": r,
                        "s": s,
                        "c": C[i][j],
                        "count": count,
                        "raw": raw,
                        "info": info,
                        "true": true_ok,
                    }
                )

        # --------------------------------------------------------
        # Correct sorting implementation.
        # --------------------------------------------------------

        relation_rows.sort(
            key=lambda item: (
                -item["info"],
                item["r"],
                item["s"],
            )
        )

        print()

        print(
            "CONSTRAINT INFORMATION RANK"
        )

        print(
            "rank cell       C     states/raw"
            "        fraction       info"
        )

        print(
            "-" * 72
        )

        for rank, item in enumerate(
            relation_rows,
            start=1
        ):

            count = item["count"]
            raw = item["raw"]
            info = item["info"]

            if info == float("inf"):

                info_text = "INF"

            else:

                info_text = (
                    f"{info:.6f}"
                )

            print(
                f"{rank:4d} "
                f"({item['r']:2d},"
                f"{item['s']:2d}) "
                f" C={item['c']} "
                f"{count:5d}/"
                f"{raw:5d} "
                f"{count/raw:.6f} "
                f"{info_text}"
            )

        print()

        true_all = all(
            item["true"]
            for item in relation_rows
        )

        print(
            "TRUE PAIR SURVIVES ALL "
            "CONSTRAINTS =",
            true_all
        )

        print()

        # --------------------------------------------------------
        # Initial domains.
        #
        # Since every radix is prime and p is a generated prime
        # larger than all radices, p mod m != 0.
        # --------------------------------------------------------

        domains = [
            set(
                range(
                    1,
                    m
                )
            )
            for m in ALL_RADICES
        ]

        raw_states = domain_product(
            domains
        )

        print(
            "RAW GLOBAL CSP STATE COUNT =",
            raw_states
        )

        print()

        # --------------------------------------------------------
        # AC-3
        # --------------------------------------------------------

        ac_revisions = ac3(
            domains,
            constraints
        )

        print(
            "AC-3 REVISIONS =",
            ac_revisions
        )

        print()

        print_domains(
            "DOMAINS AFTER AC-3",
            domains
        )

        ac_states = domain_product(
            domains
        )

        print(
            "AC-3 GLOBAL STATE COUNT =",
            ac_states
        )

        print(
            "AC-3 REDUCTION FRACTION =",
            f"{ac_states/raw_states:.12e}"
        )

        print()

        true_ac = all(
            true_x[i]
            in domains[i]
            for i in range(
                len(ALL_RADICES)
            )
        )

        print(
            "TRUE STATE SURVIVES AC-3 =",
            true_ac
        )

        if true_ac:

            aggregate_true_survival += 1

        print()

        singleton_count = sum(
            len(d) == 1
            for d in domains
        )

        empty_count = sum(
            len(d) == 0
            for d in domains
        )

        print(
            "SINGLETON VARIABLES =",
            singleton_count,
            "/",
            len(ALL_RADICES)
        )

        print(
            "EMPTY VARIABLES =",
            empty_count
        )

        print()

        # --------------------------------------------------------
        # Exact CSP enumeration
        # --------------------------------------------------------

        if empty_count == 0:

            solutions, nodes, limit_hit = (
                enumerate_solutions(
                    domains,
                    constraints,
                    MAX_SOLUTIONS
                )
            )

            print(
                "GLOBAL DFS"
            )

            print(
                "  nodes visited =",
                nodes
            )

            print(
                "  solutions found =",
                len(solutions)
            )

            print(
                "  solution limit hit =",
                limit_hit
            )

            print()

            if limit_hit:

                aggregate_limited += 1

            true_tuple = tuple(
                true_x
            )

            true_found = (
                true_tuple
                in solutions
            )

            print(
                "  TRUE GLOBAL RESIDUE "
                "VECTOR FOUND =",
                true_found
            )

            print()

            if len(solutions) == 1:

                aggregate_unique += 1

                unique = solutions[0]

                p_mod_M, M = (
                    reconstruct_full_residue(
                        unique
                    )
                )

                print(
                    "UNIQUE GLOBAL RESIDUE"
                )

                print(
                    "  p mod M =",
                    p_mod_M
                )

                print(
                    "  true p mod M =",
                    p % M
                )

                print(
                    "  M =",
                    M
                )

                print(
                    "  MATCH =",
                    p_mod_M == p % M
                )

                print()

            elif len(solutions) <= 20:

                print(
                    "GLOBAL SOLUTIONS"
                )

                for sol in solutions:

                    marker = ""

                    if sol == true_tuple:

                        marker = " <-- TRUE"

                    print(
                        "  ",
                        sol,
                        marker
                    )

                print()

        else:

            print(
                "GLOBAL DFS SKIPPED: "
                "EMPTY DOMAIN"
            )

            print()

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
print("EXPERIMENT 158 AGGREGATE")
print("=" * 72)
print()

print(
    "samples =",
    aggregate_samples
)

print(
    "true state survived AC-3 =",
    aggregate_true_survival,
    "/",
    aggregate_samples
)

print(
    "unique CSP solutions =",
    aggregate_unique,
    "/",
    aggregate_samples
)

print(
    "solution-limit cases =",
    aggregate_limited,
    "/",
    aggregate_samples
)

print()

print(
    "The important result is the number of globally consistent"
)

print(
    "radix-residue vectors remaining after AC-3 / DFS."
)

print()

print("=" * 72)
print("FINISHED EXPERIMENT 158")
print("=" * 72)