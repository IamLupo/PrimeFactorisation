#!/usr/bin/env python3

import math
import random
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple


# ========================================================================
# START EXPERIMENT 160
# Adaptive relation-join / variable elimination
# ========================================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48]
SAMPLES_PER_SIZE = 1

# Hard safety boundary.
#
# We deliberately refuse to materialize a relation larger than this.
MAX_TUPLES = 500_000

# Number of different starting cells to test.
SEED_COUNT = 3


# ========================================================================
# Prime generation
# ========================================================================

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    ]

    for p in small:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic Miller-Rabin bases for uint64.
    bases = [
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    ]

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(lo: int, hi: int) -> int:
    while True:
        x = random.randrange(lo, hi)

        if x % 2 == 0:
            x += 1

        if x >= hi:
            continue

        if is_probable_prime(x):
            return x


def generate_balanced_semiprime(
    bits: int,
) -> Tuple[int, int, int]:

    p_bits = bits // 2
    q_bits = bits - p_bits

    p_lo = 1 << (p_bits - 1)
    p_hi = 1 << p_bits

    q_lo = 1 << (q_bits - 1)
    q_hi = 1 << q_bits

    while True:
        p = random_prime(p_lo, p_hi)
        q = random_prime(q_lo, q_hi)

        if p == q:
            continue

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ========================================================================
# Exact residue / carry machinery
# ========================================================================

def crt_pair(
    x_r: int,
    r: int,
    x_s: int,
    s: int,
) -> int:
    """
    Reconstruct p mod (r*s) from:
        p mod r = x_r
        p mod s = x_s
    """

    t = (
        (x_s - x_r)
        * pow(r, -1, s)
    ) % s

    return x_r + r * t


def q_residue(
    n: int,
    p_residue: int,
    modulus: int,
) -> int:

    return (
        n
        * pow(p_residue, -1, modulus)
    ) % modulus


def c_from_residues(
    n: int,
    r: int,
    s: int,
    x_r: int,
    x_s: int,
) -> int:
    """
    Compute the exact c3 carry from only

        p mod r
        p mod s

    because q mod r and q mod s are determined by n=pq.
    """

    # ------------------------------------------------------------
    # Recover p mod r*s.
    # ------------------------------------------------------------

    p_mod_rs = crt_pair(
        x_r,
        r,
        x_s,
        s,
    )

    a = x_r

    k_mod_s = (
        (p_mod_rs - a)
        // r
    ) % s

    # ------------------------------------------------------------
    # Recover q residues.
    # ------------------------------------------------------------

    y_r = q_residue(
        n,
        x_r,
        r,
    )

    y_s = q_residue(
        n,
        x_s,
        s,
    )

    q_mod_rs = crt_pair(
        y_r,
        r,
        y_s,
        s,
    )

    b = y_s

    l_mod_r = (
        (q_mod_rs - b)
        // s
    ) % r

    # ------------------------------------------------------------
    # Exact C carry.
    # ------------------------------------------------------------

    alpha = (
        l_mod_r * a
    ) % r

    beta = (
        k_mod_s * b
    ) % s

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


def exact_c_from_pq(
    p: int,
    q: int,
    r: int,
    s: int,
) -> int:

    k, a = divmod(p, r)
    l, b = divmod(q, s)

    beta = (k * b) % s
    alpha = (l * a) % r

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


def build_c_matrix(
    p: int,
    q: int,
) -> List[List[int]]:

    return [
        [
            exact_c_from_pq(
                p,
                q,
                r,
                s,
            )
            for s in R2
        ]
        for r in R1
    ]


# ========================================================================
# Relation representation
# ========================================================================

@dataclass
class BinaryRelation:
    var1: int
    var2: int

    tuples: Set[Tuple[int, int]]

    # Number of allowed states.
    raw_size: int = 0

    @property
    def size(self) -> int:
        return len(self.tuples)

    @property
    def fraction(self) -> float:
        if self.raw_size == 0:
            return 0.0

        return self.size / self.raw_size

    @property
    def information(self) -> float:
        if self.size == 0:
            return float("inf")

        return -math.log2(
            self.fraction
        )


@dataclass
class Factor:
    """
    A higher-order relation over several residue variables.
    """

    variables: Tuple[int, ...]
    tuples: Set[Tuple[int, ...]]

    @property
    def size(self) -> int:
        return len(self.tuples)


# ========================================================================
# Build all binary relations
# ========================================================================

def build_relations(
    n: int,
    C: List[List[int]],
) -> Dict[Tuple[int, int], BinaryRelation]:

    relations = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            observed_c = C[i][j]

            allowed = set()

            for x_r in range(1, r):

                for x_s in range(1, s):

                    c = c_from_residues(
                        n,
                        r,
                        s,
                        x_r,
                        x_s,
                    )

                    if c == observed_c:
                        allowed.add(
                            (x_r, x_s)
                        )

            relations[(r, s)] = BinaryRelation(
                var1=r,
                var2=s,
                tuples=allowed,
                raw_size=(r - 1) * (s - 1),
            )

    return relations


# ========================================================================
# Relation statistics
# ========================================================================

def relation_key(
    relation: BinaryRelation,
) -> Tuple[float, int, int]:

    return (
        relation.information,
        relation.var1,
        relation.var2,
    )


def rank_relations(
    relations: Dict[Tuple[int, int], BinaryRelation],
) -> List[BinaryRelation]:

    return sorted(
        relations.values(),
        key=relation_key,
        reverse=True,
    )


# ========================================================================
# Variable domains
# ========================================================================

def domain_size(var: int) -> int:
    return var - 1


def raw_factor_size(
    variables: Tuple[int, ...],
) -> int:

    value = 1

    for var in variables:
        value *= domain_size(var)

    return value


def factor_information(
    factor: Factor,
) -> float:

    raw = raw_factor_size(
        factor.variables
    )

    if factor.size == 0:
        return float("inf")

    return math.log2(
        raw / factor.size
    )


# ========================================================================
# True-state checks
# ========================================================================

def true_assignment(
    p: int,
    variables: Tuple[int, ...],
) -> Tuple[int, ...]:

    return tuple(
        p % v
        for v in variables
    )


def true_survives(
    factor: Factor,
    p: int,
) -> bool:

    assignment = true_assignment(
        p,
        factor.variables,
    )

    return assignment in factor.tuples


# ========================================================================
# Join helpers
# ========================================================================

def common_variables(
    factor: Factor,
    relation: BinaryRelation,
) -> List[int]:

    current = set(
        factor.variables
    )

    return [
        v
        for v in (
            relation.var1,
            relation.var2,
        )
        if v in current
    ]


def candidate_join_size_one_common(
    factor: Factor,
    relation: BinaryRelation,
    common_var: int,
) -> int:
    """
    Exact natural-join cardinality estimate.

    Since relation is binary and exactly one variable is shared:

        |F JOIN R|
        =
        sum_x deg_F(x) * deg_R(x)
    """

    pos = factor.variables.index(
        common_var
    )

    current_degree = defaultdict(int)

    for tup in factor.tuples:
        current_degree[tup[pos]] += 1

    if relation.var1 == common_var:
        common_pos = 0
        new_pos = 1
    else:
        common_pos = 1
        new_pos = 0

    relation_degree = defaultdict(int)

    for tup in relation.tuples:
        relation_degree[
            tup[common_pos]
        ] += 1

    total = 0

    for value, degree in current_degree.items():
        total += (
            degree
            * relation_degree.get(
                value,
                0,
            )
        )

    return total


def expected_join_size_one_common(
    factor: Factor,
    relation: BinaryRelation,
    common_var: int,
) -> float:
    """
    Independence prediction.

    Only one new variable is introduced.
    """

    d_common = domain_size(
        common_var
    )

    new_var = (
        relation.var2
        if relation.var1 == common_var
        else relation.var1
    )

    d_new = domain_size(
        new_var
    )

    candidate_fraction = (
        relation.size
        / (d_common * d_new)
    )

    return (
        factor.size
        * candidate_fraction
    )


def join_one_common(
    factor: Factor,
    relation: BinaryRelation,
    common_var: int,
) -> Factor:

    pos = factor.variables.index(
        common_var
    )

    # The new variable is added to the end.
    new_var = (
        relation.var2
        if relation.var1 == common_var
        else relation.var1
    )

    # ------------------------------------------------------------
    # Index binary relation by common value.
    # ------------------------------------------------------------

    if relation.var1 == common_var:
        common_pos = 0
        new_pos = 1
    else:
        common_pos = 1
        new_pos = 0

    rel_index = defaultdict(list)

    for tup in relation.tuples:
        rel_index[
            tup[common_pos]
        ].append(
            tup[new_pos]
        )

    # ------------------------------------------------------------
    # Natural join.
    # ------------------------------------------------------------

    new_tuples = set()

    for tup in factor.tuples:

        common_value = tup[pos]

        for new_value in rel_index.get(
            common_value,
            (),
        ):

            new_tuples.add(
                tup + (new_value,)
            )

    return Factor(
        variables=
            factor.variables
            + (new_var,),
        tuples=new_tuples,
    )


def filter_factor(
    factor: Factor,
    relation: BinaryRelation,
) -> Factor:
    """
    Both variables of the binary relation already exist
    in the higher-order factor.

    The relation therefore acts as a pure filter.
    """

    p1 = factor.variables.index(
        relation.var1
    )

    p2 = factor.variables.index(
        relation.var2
    )

    rel = relation.tuples

    new_tuples = {
        tup
        for tup in factor.tuples
        if (
            tup[p1],
            tup[p2],
        ) in rel
    }

    return Factor(
        variables=factor.variables,
        tuples=new_tuples,
    )


def filter_expected_size(
    factor: Factor,
    relation: BinaryRelation,
) -> float:

    return (
        factor.size
        * relation.fraction
    )


# ========================================================================
# Adaptive path search
# ========================================================================

def choose_next_relation(
    factor: Factor,
    relations: Dict[Tuple[int, int], BinaryRelation],
    used: Set[Tuple[int, int]],
):
    """
    Adaptive selection:

    1. Prefer a relation whose BOTH variables already exist.
       Such a relation is a pure filter and may collapse the factor.

    2. Otherwise consider relations with exactly one shared variable.

    Among one-shared relations, prefer the one with the smallest
    exact predicted join size, with information as a tie breaker.
    """

    current_vars = set(
        factor.variables
    )

    # ------------------------------------------------------------
    # Pure filters.
    # ------------------------------------------------------------

    filters = []

    for key, relation in relations.items():

        if key in used:
            continue

        shared = common_variables(
            factor,
            relation,
        )

        if len(shared) == 2:

            filters.append(
                relation
            )

    if filters:

        filters.sort(
            key=lambda rel: (
                rel.fraction,
                -rel.information,
            )
        )

        return filters[0], "filter"

    # ------------------------------------------------------------
    # One-variable joins.
    # ------------------------------------------------------------

    candidates = []

    for key, relation in relations.items():

        if key in used:
            continue

        shared = common_variables(
            factor,
            relation,
        )

        if len(shared) != 1:
            continue

        common_var = shared[0]

        predicted = candidate_join_size_one_common(
            factor,
            relation,
            common_var,
        )

        candidates.append(
            (
                predicted,
                -relation.information,
                relation,
                common_var,
            )
        )

    if not candidates:
        return None, None

    candidates.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    _, _, relation, common_var = candidates[0]

    return (
        relation,
        common_var,
    )


# ========================================================================
# Single adaptive run
# ========================================================================

def run_path(
    n: int,
    p: int,
    q: int,
    relations: Dict[Tuple[int, int], BinaryRelation],
    seed: BinaryRelation,
    seed_number: int,
):

    print()
    print(
        "======================================================================="
    )
    print(
        f"SEED {seed_number}: "
        f"({seed.var1},{seed.var2}) "
        f"C relation"
    )
    print(
        "======================================================================="
    )

    print(
        f"seed states = {seed.size}/{seed.raw_size}"
    )

    print(
        f"seed information = "
        f"{seed.information:.6f} bits"
    )

    factor = Factor(
        variables=(
            seed.var1,
            seed.var2,
        ),
        tuples=set(seed.tuples),
    )

    used = {
        (
            seed.var1,
            seed.var2,
        )
    }

    print(
        f"true state survives = "
        f"{true_survives(factor, p)}"
    )

    step = 0

    while True:

        choice, mode = choose_next_relation(
            factor,
            relations,
            used,
        )

        if choice is None:
            print()
            print(
                "NO MORE CONNECTED RELATIONS"
            )
            break

        relation = choice

        key = (
            relation.var1,
            relation.var2,
        )

        step += 1

        print()
        print(
            f"STEP {step}"
        )

        print(
            f"  candidate = "
            f"({relation.var1},{relation.var2})"
        )

        print(
            f"  candidate relation size = "
            f"{relation.size}/{relation.raw_size}"
        )

        print(
            f"  candidate information = "
            f"{relation.information:.6f} bits"
        )

        before_size = factor.size

        # ------------------------------------------------------------
        # Filter: both variables already exist.
        # ------------------------------------------------------------

        if mode == "filter":

            expected = filter_expected_size(
                factor,
                relation,
            )

            print(
                f"  operation = FILTER"
            )

            print(
                f"  current states = "
                f"{before_size}"
            )

            print(
                f"  independent expected = "
                f"{expected:.3f}"
            )

            factor = filter_factor(
                factor,
                relation,
            )

            actual = factor.size

            print(
                f"  actual states = "
                f"{actual}"
            )

            if expected > 0:
                print(
                    f"  correlation ratio = "
                    f"{actual / expected:.6f}"
                )

        # ------------------------------------------------------------
        # Join: exactly one variable shared.
        # ------------------------------------------------------------

        else:

            common_var = mode

            predicted = candidate_join_size_one_common(
                factor,
                relation,
                common_var,
            )

            expected = expected_join_size_one_common(
                factor,
                relation,
                common_var,
            )

            print(
                f"  operation = JOIN"
            )

            print(
                f"  shared variable = "
                f"{common_var}"
            )

            print(
                f"  exact predicted size = "
                f"{predicted}"
            )

            print(
                f"  independence expected = "
                f"{expected:.3f}"
            )

            if predicted > MAX_TUPLES:

                print()
                print(
                    "  *** JOIN CAP REACHED ***"
                )

                print(
                    f"  predicted {predicted:,} tuples"
                )

                print(
                    f"  maximum {MAX_TUPLES:,} tuples"
                )

                print()
                print(
                    "  STOPPING THIS PATH"
                )

                break

            factor = join_one_common(
                factor,
                relation,
                common_var,
            )

            actual = factor.size

            print(
                f"  actual states = "
                f"{actual}"
            )

            print(
                f"  variables now = "
                f"{factor.variables}"
            )

            if expected > 0:
                print(
                    f"  correlation ratio = "
                    f"{actual / expected:.6f}"
                )

        # ------------------------------------------------------------
        # Global statistics.
        # ------------------------------------------------------------

        used.add(key)

        raw = raw_factor_size(
            factor.variables
        )

        info = factor_information(
            factor
        )

        print(
            f"  raw state space = "
            f"{raw}"
        )

        print(
            f"  information = "
            f"{info:.6f} bits"
        )

        print(
            f"  surviving fraction = "
            f"{factor.size / raw:.12e}"
        )

        survives = true_survives(
            factor,
            p,
        )

        print(
            f"  TRUE STATE SURVIVES = "
            f"{survives}"
        )

        if factor.size == 0:

            print()
            print(
                "  *** EMPTY RELATION ***"
            )

            break

        if factor.size > MAX_TUPLES:

            print()
            print(
                "  *** MATERIALIZED RELATION ABOVE CAP ***"
            )

            break

    # --------------------------------------------------------------------
    # Final result.
    # --------------------------------------------------------------------

    print()
    print(
        "PATH SUMMARY"
    )

    print(
        f"  final variables = "
        f"{factor.variables}"
    )

    print(
        f"  final states = "
        f"{factor.size}"
    )

    print(
        f"  final raw states = "
        f"{raw_factor_size(factor.variables)}"
    )

    print(
        f"  final information = "
        f"{factor_information(factor):.6f} bits"
    )

    print(
        f"  TRUE STATE SURVIVES = "
        f"{true_survives(factor, p)}"
    )

    return factor


# ========================================================================
# Main sample
# ========================================================================

def run_sample(
    bits: int,
    sample_index: int,
):

    print()
    print(
        "------------------------------------------------------------------------"
    )

    print(
        f"SAMPLE {sample_index}/{SAMPLES_PER_SIZE}"
    )

    n, p, q = generate_balanced_semiprime(
        bits
    )

    print(
        f"n bits = {bits}"
    )

    print(
        f"n = {n}"
    )

    print(
        f"true p = {p}"
    )

    print(
        f"true q = {q}"
    )

    print()
    print(
        "TRUE p RESIDUES"
    )

    for r in R1:
        print(
            f"  p mod {r} = {p % r}"
        )

    for s in R2:
        print(
            f"  p mod {s} = {p % s}"
        )

    C = build_c_matrix(
        p,
        q,
    )

    print()
    print(
        "C MATRIX"
    )

    for row in C:
        print(row)

    print()
    print(
        "BUILDING EXACT BINARY RELATIONS"
    )

    relations = build_relations(
        n,
        C,
    )

    print("DONE")

    # ------------------------------------------------------------
    # Rank cells by information.
    # ------------------------------------------------------------

    ranked = rank_relations(
        relations
    )

    print()
    print(
        "ALL CELL INFORMATION"
    )

    print(
        "rank cell       C     states/raw       "
        "fraction       information"
    )

    print(
        "------------------------------------------------------------------------"
    )

    for rank, rel in enumerate(
        ranked,
        start=1,
    ):

        i = R1.index(rel.var1)
        j = R2.index(rel.var2)

        print(
            f"{rank:4d} "
            f"({rel.var1},{rel.var2}) "
            f"C={C[i][j]} "
            f"{rel.size:6d}/{rel.raw_size:<6d} "
            f"{rel.fraction:.8f} "
            f"{rel.information:10.6f}"
        )

    # ------------------------------------------------------------
    # Run adaptive paths from strongest cells.
    # ------------------------------------------------------------

    for seed_number, seed in enumerate(
        ranked[:SEED_COUNT],
        start=1,
    ):

        run_path(
            n,
            p,
            q,
            relations,
            seed,
            seed_number,
        )


# ========================================================================
# Driver
# ========================================================================

def main():

    random.seed(160)

    print("=" * 72)
    print("START EXPERIMENT 160")
    print("Adaptive relation-join / variable elimination")
    print("=" * 72)

    print()
    print(f"R1 = {R1}")
    print(f"R2 = {R2}")

    print()
    print(
        f"MAX_TUPLES = {MAX_TUPLES:,}"
    )

    print(
        f"SEED_COUNT = {SEED_COUNT}"
    )

    for bits in BIT_SIZES:

        print()
        print("=" * 72)
        print(
            f"BIT SIZE = {bits}"
        )
        print("=" * 72)

        for sample_index in range(
            1,
            SAMPLES_PER_SIZE + 1,
        ):

            run_sample(
                bits,
                sample_index,
            )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 160")
    print("=" * 72)


if __name__ == "__main__":
    main()
