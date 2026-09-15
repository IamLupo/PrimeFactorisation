#!/usr/bin/env python3

import math
import random
from collections import Counter
from typing import Dict, List, Set, Tuple


# ========================================================================
# START EXPERIMENT 161
# Carry relation geometry and degree structure
# ========================================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 2

TOP_CELLS = 10
TOP_DETAILED = 3


# ========================================================================
# Prime generation
# ========================================================================

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    ]

    for p in small_primes:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic Miller-Rabin for 64-bit integers.
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
        a %= n

        if a == 0:
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

        p = random_prime(
            p_lo,
            p_hi,
        )

        q = random_prime(
            q_lo,
            q_hi,
        )

        if p == q:
            continue

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ========================================================================
# Exact carry
# ========================================================================

def exact_c(
    p: int,
    q: int,
    r: int,
    s: int,
) -> int:

    k, a = divmod(p, r)
    ell, b = divmod(q, s)

    beta = (k * b) % s
    alpha = (ell * a) % r

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
            exact_c(
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
# CRT / global residue reconstruction
# ========================================================================

def crt_pair(
    x_r: int,
    r: int,
    x_s: int,
    s: int,
) -> int:
    """
    Reconstruct z mod r*s from

        z mod r = x_r
        z mod s = x_s

    with r,s distinct primes.
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


def c_from_global_residues(
    n: int,
    r: int,
    s: int,
    x: int,
    y: int,
) -> int:
    """
    Compute the exact c3 carry using only

        x = p mod r
        y = p mod s.

    Since n = p*q, q mod r and q mod s are determined
    by the inverse of p modulo those radices.
    """

    # ------------------------------------------------------------
    # Recover p mod r*s.
    # ------------------------------------------------------------

    p_mod_rs = crt_pair(
        x,
        r,
        y,
        s,
    )

    a = x

    k_mod_s = (
        (p_mod_rs - a)
        // r
    ) % s

    # ------------------------------------------------------------
    # Recover q mod r and q mod s.
    # ------------------------------------------------------------

    q_r = q_residue(
        n,
        x,
        r,
    )

    q_s = q_residue(
        n,
        y,
        s,
    )

    q_mod_rs = crt_pair(
        q_r,
        r,
        q_s,
        s,
    )

    b = q_s

    ell_mod_r = (
        (q_mod_rs - b)
        // s
    ) % r

    # ------------------------------------------------------------
    # Exact c3.
    # ------------------------------------------------------------

    alpha = (
        ell_mod_r * a
    ) % r

    beta = (
        k_mod_s * b
    ) % s

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


# ========================================================================
# Relation construction
# ========================================================================

def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> Set[Tuple[int, int]]:

    relation: Set[Tuple[int, int]] = set()

    for x in range(1, r):

        for y in range(1, s):

            c = c_from_global_residues(
                n,
                r,
                s,
                x,
                y,
            )

            if c == target_c:
                relation.add(
                    (x, y)
                )

    return relation


# ========================================================================
# Degree analysis
# ========================================================================

def entropy_from_counts(
    counts: List[int],
) -> float:

    total = sum(counts)

    if total == 0:
        return 0.0

    h = 0.0

    for count in counts:

        if count <= 0:
            continue

        p = count / total

        h -= p * math.log2(p)

    return h


def degree_analysis(
    relation: Set[Tuple[int, int]],
    r: int,
    s: int,
) -> Dict:

    dx = Counter()
    dy = Counter()

    neighbors_x: Dict[int, Set[int]] = {}
    neighbors_y: Dict[int, Set[int]] = {}

    for x, y in relation:

        dx[x] += 1
        dy[y] += 1

        neighbors_x.setdefault(
            x,
            set(),
        ).add(y)

        neighbors_y.setdefault(
            y,
            set(),
        ).add(x)

    x_degrees = [
        dx[x]
        for x in range(1, r)
    ]

    y_degrees = [
        dy[y]
        for y in range(1, s)
    ]

    active_x = [
        d
        for d in x_degrees
        if d > 0
    ]

    active_y = [
        d
        for d in y_degrees
        if d > 0
    ]

    # ------------------------------------------------------------
    # Conditional entropy:
    #
    # H(Y|X) = average log2(degree_X)
    # H(X|Y) = average log2(degree_Y)
    # ------------------------------------------------------------

    total = len(relation)

    h_y_given_x = 0.0

    for x in range(1, r):

        d = dx[x]

        if d:
            h_y_given_x += (
                (d / total)
                * math.log2(d)
            )

    h_x_given_y = 0.0

    for y in range(1, s):

        d = dy[y]

        if d:
            h_x_given_y += (
                (d / total)
                * math.log2(d)
            )

    # ------------------------------------------------------------
    # Basic degree stats.
    # ------------------------------------------------------------

    def stats(values: List[int]) -> Dict:

        nonzero = [
            d
            for d in values
            if d > 0
        ]

        if not nonzero:
            return {
                "min": 0,
                "max": 0,
                "mean": 0.0,
                "median": 0.0,
                "zero": len(values),
            }

        ordered = sorted(nonzero)
        n = len(ordered)

        if n % 2:
            median = ordered[n // 2]
        else:
            median = (
                ordered[n // 2 - 1]
                + ordered[n // 2]
            ) / 2.0

        return {
            "min": min(nonzero),
            "max": max(nonzero),
            "mean": sum(nonzero) / n,
            "median": median,
            "zero": len(values) - n,
        }

    sx = stats(x_degrees)
    sy = stats(y_degrees)

    # ------------------------------------------------------------
    # Fraction of active rows with small degree.
    # ------------------------------------------------------------

    def branch_fraction(
        values: List[int],
        threshold: int,
    ) -> float:

        if not values:
            return 0.0

        return (
            sum(
                d <= threshold
                for d in values
            )
            / len(values)
        )

    branch_scores_x = {
        t: branch_fraction(
            active_x,
            t,
        )
        for t in [1, 2, 3, 4]
    }

    branch_scores_y = {
        t: branch_fraction(
            active_y,
            t,
        )
        for t in [1, 2, 3, 4]
    }

    return {
        "dx": dx,
        "dy": dy,
        "neighbors_x": neighbors_x,
        "neighbors_y": neighbors_y,
        "x_stats": sx,
        "y_stats": sy,
        "h_y_given_x": h_y_given_x,
        "h_x_given_y": h_x_given_y,
        "x_function_score": branch_scores_x[1],
        "y_function_score": branch_scores_y[1],
        "branch_scores_x": branch_scores_x,
        "branch_scores_y": branch_scores_y,
    }


# ========================================================================
# Geometry
# ========================================================================

def contiguous_runs(
    values: List[int],
) -> int:

    if not values:
        return 0

    values = sorted(values)

    runs = 1

    for previous, current in zip(
        values,
        values[1:],
    ):

        if current != previous + 1:
            runs += 1

    return runs


def best_arithmetic_step(
    values: List[int],
) -> Tuple[int, float]:

    if len(values) < 2:
        return 0, 1.0

    values = sorted(values)

    differences = Counter(
        b - a
        for a, b in zip(
            values,
            values[1:],
        )
    )

    step, count = differences.most_common(1)[0]

    score = (
        count
        / (len(values) - 1)
    )

    return step, score


def geometry_analysis(
    neighbors: Dict[int, Set[int]],
) -> Dict:

    run_counts = []
    step_scores = []

    for values in neighbors.values():

        vals = sorted(values)

        if not vals:
            continue

        run_counts.append(
            contiguous_runs(vals)
        )

        _, score = best_arithmetic_step(
            vals
        )

        step_scores.append(score)

    return {
        "mean_runs": (
            sum(run_counts)
            / len(run_counts)
            if run_counts
            else 0.0
        ),
        "max_runs": (
            max(run_counts)
            if run_counts
            else 0
        ),
        "mean_best_step_score": (
            sum(step_scores)
            / len(step_scores)
            if step_scores
            else 0.0
        ),
    }


# ========================================================================
# Marginal degree entropy
# ========================================================================

def marginal_degree_entropy(
    relation: Set[Tuple[int, int]],
    r: int,
    s: int,
) -> Dict:

    dx = Counter(
        x
        for x, _ in relation
    )

    dy = Counter(
        y
        for _, y in relation
    )

    x_counts = [
        dx[x]
        for x in range(1, r)
    ]

    y_counts = [
        dy[y]
        for y in range(1, s)
    ]

    return {
        "x": entropy_from_counts(
            x_counts
        ),
        "y": entropy_from_counts(
            y_counts
        ),
    }


# ========================================================================
# Complete cell analysis
# ========================================================================

def analyse_cells(
    n: int,
    p: int,
    q: int,
):

    C = build_c_matrix(
        p,
        q,
    )

    records = []

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            target_c = C[i][j]

            relation = build_relation(
                n,
                r,
                s,
                target_c,
            )

            raw = (
                (r - 1)
                * (s - 1)
            )

            fraction = (
                len(relation)
                / raw
            )

            information = (
                -math.log2(fraction)
            )

            true_pair = (
                p % r,
                p % s,
            )

            survives = (
                true_pair in relation
            )

            degrees = degree_analysis(
                relation,
                r,
                s,
            )

            geometry_x = geometry_analysis(
                degrees["neighbors_x"]
            )

            geometry_y = geometry_analysis(
                degrees["neighbors_y"]
            )

            marginal_entropy = (
                marginal_degree_entropy(
                    relation,
                    r,
                    s,
                )
            )

            records.append(
                {
                    "i": i,
                    "j": j,
                    "r": r,
                    "s": s,
                    "c": target_c,
                    "relation": relation,
                    "raw": raw,
                    "size": len(relation),
                    "fraction": fraction,
                    "information": information,
                    "true_survives": survives,
                    "degrees": degrees,
                    "geometry_x": geometry_x,
                    "geometry_y": geometry_y,
                    "marginal_entropy": marginal_entropy,
                }
            )

    return C, records


# ========================================================================
# Printing cell table
# ========================================================================

def print_cell_table(
    records,
):

    ranked = sorted(
        records,
        key=lambda rec: rec["information"],
        reverse=True,
    )

    print()
    print(
        "CELL RELATIONS RANKED BY INFORMATION"
    )

    print(
        "rank cell       C   states/raw       "
        "info      H(Y|X)   H(X|Y)   maxX maxY"
    )

    print(
        "------------------------------------------------------------------------"
    )

    for rank, rec in enumerate(
        ranked[:TOP_CELLS],
        start=1,
    ):

        degrees = rec["degrees"]

        print(
            f"{rank:4d} "
            f"({rec['r']:2d},{rec['s']:2d}) "
            f"C={rec['c']} "
            f"{rec['size']:6d}/{rec['raw']:<6d} "
            f"{rec['information']:8.4f} "
            f"{degrees['h_y_given_x']:9.4f} "
            f"{degrees['h_x_given_y']:9.4f} "
            f"{degrees['x_stats']['max']:4d} "
            f"{degrees['y_stats']['max']:4d}"
        )

    return ranked


# ========================================================================
# Detailed relation printer
# ========================================================================

def print_detailed_relation(
    rec,
):

    r = rec["r"]
    s = rec["s"]
    c = rec["c"]

    degrees = rec["degrees"]

    print()
    print(
        "------------------------------------------------------------------------"
    )

    print(
        f"DETAILED RELATION ({r},{s}) C={c}"
    )

    print(
        "------------------------------------------------------------------------"
    )

    print(
        f"states = "
        f"{rec['size']}/{rec['raw']}"
    )

    print(
        f"fraction = "
        f"{rec['fraction']:.12e}"
    )

    print(
        f"information = "
        f"{rec['information']:.6f} bits"
    )

    print(
        f"H(Y|X) = "
        f"{degrees['h_y_given_x']:.6f} bits"
    )

    print(
        f"H(X|Y) = "
        f"{degrees['h_x_given_y']:.6f} bits"
    )

    print()
    print(
        "X DEGREE STATISTICS"
    )

    print(
        f"  min active degree = "
        f"{degrees['x_stats']['min']}"
    )

    print(
        f"  max degree         = "
        f"{degrees['x_stats']['max']}"
    )

    print(
        f"  mean degree        = "
        f"{degrees['x_stats']['mean']:.6f}"
    )

    print(
        f"  median degree      = "
        f"{degrees['x_stats']['median']:.6f}"
    )

    print(
        f"  zero-degree X      = "
        f"{degrees['x_stats']['zero']}"
    )

    print()
    print(
        "Y DEGREE STATISTICS"
    )

    print(
        f"  min active degree = "
        f"{degrees['y_stats']['min']}"
    )

    print(
        f"  max degree         = "
        f"{degrees['y_stats']['max']}"
    )

    print(
        f"  mean degree        = "
        f"{degrees['y_stats']['mean']:.6f}"
    )

    print(
        f"  median degree      = "
        f"{degrees['y_stats']['median']:.6f}"
    )

    print(
        f"  zero-degree Y      = "
        f"{degrees['y_stats']['zero']}"
    )

    print()
    print(
        "BRANCH SCORES"
    )

    for threshold in [1, 2, 3, 4]:

        print(
            f"  X degree <= {threshold}: "
            f"{degrees['branch_scores_x'][threshold] * 100:.2f}%"
        )

        print(
            f"  Y degree <= {threshold}: "
            f"{degrees['branch_scores_y'][threshold] * 100:.2f}%"
        )

    print()
    print(
        "GEOMETRY"
    )

    print(
        f"  X mean contiguous runs = "
        f"{rec['geometry_x']['mean_runs']:.6f}"
    )

    print(
        f"  X max contiguous runs = "
        f"{rec['geometry_x']['max_runs']}"
    )

    print(
        f"  X mean best arithmetic-step score = "
        f"{rec['geometry_x']['mean_best_step_score']:.6f}"
    )

    print(
        f"  Y mean contiguous runs = "
        f"{rec['geometry_y']['mean_runs']:.6f}"
    )

    print(
        f"  Y max contiguous runs = "
        f"{rec['geometry_y']['max_runs']}"
    )

    print(
        f"  Y mean best arithmetic-step score = "
        f"{rec['geometry_y']['mean_best_step_score']:.6f}"
    )

    print()
    print(
        "X -> Y ALLOWED SETS"
    )

    # IMPORTANT:
    #
    # neighbors_x.items() returns:
    #
    #     x, values
    #
    # not a five-element tuple.
    #
    # This is the bug fixed from the previous script.

    for x, values in sorted(
        degrees["neighbors_x"].items()
    ):

        vals = sorted(values)

        runs = contiguous_runs(
            vals
        )

        step, step_score = (
            best_arithmetic_step(vals)
        )

        if len(vals) <= 16:

            display = str(vals)

        else:

            display = (
                str(vals[:7])
                + " ... "
                + str(vals[-5:])
            )

        print(
            f"  x={x:3d} "
            f"degree={len(vals):3d} "
            f"runs={runs:2d} "
            f"step={step:3d} "
            f"step-score={step_score:.3f} "
            f"y={display}"
        )

    print()
    print(
        "Y -> X ALLOWED SETS"
    )

    for y, values in sorted(
        degrees["neighbors_y"].items()
    ):

        vals = sorted(values)

        runs = contiguous_runs(
            vals
        )

        step, step_score = (
            best_arithmetic_step(vals)
        )

        if len(vals) <= 16:

            display = str(vals)

        else:

            display = (
                str(vals[:7])
                + " ... "
                + str(vals[-5:])
            )

        print(
            f"  y={y:3d} "
            f"degree={len(vals):3d} "
            f"runs={runs:2d} "
            f"step={step:3d} "
            f"step-score={step_score:.3f} "
            f"x={display}"
        )

    true_pair = (
        None,
        None,
    )

    print()


# ========================================================================
# Global summary for one sample
# ========================================================================

def summarize_sample(
    records,
):

    print()
    print(
        "GLOBAL CELL SUMMARY"
    )

    mean_info = (
        sum(
            rec["information"]
            for rec in records
        )
        / len(records)
    )

    max_info = max(
        rec["information"]
        for rec in records
    )

    mean_h_yx = (
        sum(
            rec["degrees"]["h_y_given_x"]
            for rec in records
        )
        / len(records)
    )

    mean_h_xy = (
        sum(
            rec["degrees"]["h_x_given_y"]
            for rec in records
        )
        / len(records)
    )

    mean_runs_x = (
        sum(
            rec["geometry_x"]["mean_runs"]
            for rec in records
        )
        / len(records)
    )

    mean_runs_y = (
        sum(
            rec["geometry_y"]["mean_runs"]
            for rec in records
        )
        / len(records)
    )

    print(
        f"mean cell information = "
        f"{mean_info:.6f} bits"
    )

    print(
        f"max cell information = "
        f"{max_info:.6f} bits"
    )

    print(
        f"mean H(Y|X) = "
        f"{mean_h_yx:.6f} bits"
    )

    print(
        f"mean H(X|Y) = "
        f"{mean_h_xy:.6f} bits"
    )

    print(
        f"mean X contiguous runs = "
        f"{mean_runs_x:.6f}"
    )

    print(
        f"mean Y contiguous runs = "
        f"{mean_runs_y:.6f}"
    )

    failed = sum(
        not rec["true_survives"]
        for rec in records
    )

    print(
        f"TRUE RELATION SURVIVAL FAILURES = "
        f"{failed}"
    )


# ========================================================================
# Cross-size summary
# ========================================================================

def print_cross_size_summary(
    summaries,
):

    print()
    print(
        "======================================================================="
    )

    print(
        "CARRY-RELATION GEOMETRY ACROSS BIT SIZES"
    )

    print(
        "======================================================================="
    )

    print(
        "bits   meanInfo   maxInfo   "
        "meanH(Y|X)   meanH(X|Y)   "
        "meanRunsX   meanRunsY"
    )

    print(
        "------------------------------------------------------------------------"
    )

    for bits, records_list in summaries:

        records = [
            rec
            for sample_records in records_list
            for rec in sample_records
        ]

        mean_info = (
            sum(
                rec["information"]
                for rec in records
            )
            / len(records)
        )

        max_info = max(
            rec["information"]
            for rec in records
        )

        mean_h_yx = (
            sum(
                rec["degrees"]["h_y_given_x"]
                for rec in records
            )
            / len(records)
        )

        mean_h_xy = (
            sum(
                rec["degrees"]["h_x_given_y"]
                for rec in records
            )
            / len(records)
        )

        mean_runs_x = (
            sum(
                rec["geometry_x"]["mean_runs"]
                for rec in records
            )
            / len(records)
        )

        mean_runs_y = (
            sum(
                rec["geometry_y"]["mean_runs"]
                for rec in records
            )
            / len(records)
        )

        print(
            f"{bits:4d} "
            f"{mean_info:10.5f} "
            f"{max_info:9.5f} "
            f"{mean_h_yx:12.5f} "
            f"{mean_h_xy:12.5f} "
            f"{mean_runs_x:11.5f} "
            f"{mean_runs_y:11.5f}"
        )


# ========================================================================
# Main sample runner
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

    for m in R1 + R2:

        print(
            f"  p mod {m} = {p % m}"
        )

    C, records = analyse_cells(
        n,
        p,
        q,
    )

    print()
    print(
        "C MATRIX"
    )

    for row in C:
        print(row)

    ranked = print_cell_table(
        records
    )

    summarize_sample(
        records
    )

    print()
    print(
        "DETAILED STRONGEST RELATIONS"
    )

    for rec in ranked[:TOP_DETAILED]:

        print_detailed_relation(
            rec
        )

    return records


# ========================================================================
# Main
# ========================================================================

def main():

    random.seed(161)

    print("=" * 72)
    print("START EXPERIMENT 161")
    print("Carry relation geometry and degree structure")
    print("=" * 72)

    print()
    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    summaries = []

    for bits in BIT_SIZES:

        print()
        print("=" * 72)
        print(
            f"BIT SIZE = {bits}"
        )
        print("=" * 72)

        records_list = []

        for sample_index in range(
            1,
            SAMPLES_PER_SIZE + 1,
        ):

            records = run_sample(
                bits,
                sample_index,
            )

            records_list.append(
                records
            )

        summaries.append(
            (
                bits,
                records_list,
            )
        )

    print_cross_size_summary(
        summaries
    )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 161")
    print("=" * 72)


if __name__ == "__main__":
    main()