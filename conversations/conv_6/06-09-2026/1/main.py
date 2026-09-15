#!/usr/bin/env python3

import math
import random
import itertools
from typing import Dict, List, Tuple, Set


# ========================================================================
# START EXPERIMENT 159
# Higher-order radix/carry constraint composition
# ========================================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 2

BLOCK_SHAPES = [
    (2, 2),
    (2, 3),
    (3, 2),
    (3, 3),
]

TOP_K = 10


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

    # Deterministic enough for our small experimental range.
    bases = [2, 3, 5, 7, 11, 13, 17]

    for a in bases:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


def random_prime(lo: int, hi: int) -> int:
    if lo % 2 == 0:
        lo += 1

    while True:
        x = random.randrange(lo, hi)

        if x % 2 == 0:
            x += 1

        if x >= hi:
            continue

        if is_probable_prime(x):
            return x


def generate_balanced_semiprime(bits: int) -> Tuple[int, int, int]:
    """
    Generate p*q with approximately balanced prime factors.
    """

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
# Exact carry decomposition
# ========================================================================

def cell_components(
    p: int,
    q: int,
    r: int,
    s: int,
) -> Tuple[int, int, int, int, int, int, int]:
    """
    Exact decomposition for one radix pair.

    p = r*k + a
    q = s*l + b

    c1 = floor(k*b/s)
    c2 = floor(l*a/r)

    beta  = (k*b) mod s
    alpha = (l*a) mod r

    c3 =
        floor((r*beta + s*alpha + a*b)/(r*s))

    E = c1 + c2 + c3
    """

    k, a = divmod(p, r)
    l, b = divmod(q, s)

    c1 = (k * b) // s
    c2 = (l * a) // r

    beta = (k * b) % s
    alpha = (l * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    return (
        k,
        l,
        a,
        b,
        alpha,
        beta,
        c3,
    )


def exact_c(
    p: int,
    q: int,
    r: int,
    s: int,
) -> int:

    (
        _k,
        _l,
        _a,
        _b,
        _alpha,
        _beta,
        c3,
    ) = cell_components(p, q, r, s)

    return c3


def build_c_matrix(
    p: int,
    q: int,
) -> List[List[int]]:

    return [
        [
            exact_c(p, q, r, s)
            for s in R2
        ]
        for r in R1
    ]


# ========================================================================
# Relation construction
# ========================================================================

# Relation:
#
#   relation[i][j][x]
#
# is a bit-mask containing all y in 1..s-1 such that
#
#   C_ij(x,y) == observed_C_ij
#
# where
#
#   x = p mod r_i
#   y = p mod s_j
#
# q mod r/s is recovered through
#
#   q = n * p^{-1} mod modulus
#
# so the cell is indeed determined by the two p residues.


def q_residue(
    n: int,
    x: int,
    modulus: int,
) -> int:

    return (n * pow(x, -1, modulus)) % modulus


def crt_pair(
    x_r: int,
    r: int,
    x_s: int,
    s: int,
) -> int:
    """
    Reconstruct p mod r*s from:
        p mod r = x_r
        p mod s = x_s

    r,s are distinct primes.
    """

    # p = x_r + r*t
    #
    # x_r + r*t = x_s (mod s)
    #
    # r*t = x_s - x_r (mod s)

    t = (
        (x_s - x_r)
        * pow(r, -1, s)
    ) % s

    return x_r + r * t


def cell_c_from_residues(
    n: int,
    r: int,
    s: int,
    x_r: int,
    x_s: int,
) -> int:

    p_mod_rs = crt_pair(
        x_r,
        r,
        x_s,
        s,
    )

    # Since p ≡ p_mod_rs mod r*s,
    #
    # k = (p-a)/r
    #
    # modulo s is determined by p mod r*s.

    a = x_r

    k_mod_s = (
        (p_mod_rs - a)
        // r
    ) % s

    # q residues are determined from n/p mod modulus.

    y_r = q_residue(n, x_r, r)
    y_s = q_residue(n, x_s, s)

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

    alpha = (l_mod_r * a) % r
    beta = (k_mod_s * b) % s

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    return c3


def build_relations(
    n: int,
    C: List[List[int]],
) -> List[List[Dict[int, int]]]:

    relations = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            observed = C[i][j]

            mapping: Dict[int, int] = {}

            for x in range(1, r):

                mask = 0

                for y in range(1, s):

                    c = cell_c_from_residues(
                        n,
                        r,
                        s,
                        x,
                        y,
                    )

                    if c == observed:
                        mask |= 1 << (y - 1)

                mapping[x] = mask

            row.append(mapping)

        relations.append(row)

    return relations


# ========================================================================
# Utilities
# ========================================================================

def bit_values(mask: int, modulus: int) -> List[int]:
    values = []

    for y in range(1, modulus):

        if mask & (1 << (y - 1)):
            values.append(y)

    return values


def bit_count(mask: int) -> int:
    return mask.bit_count()


def full_mask(modulus: int) -> int:
    return (1 << (modulus - 1)) - 1


def block_raw_state_count(
    rows: Tuple[int, ...],
    cols: Tuple[int, ...],
) -> int:

    value = 1

    for i in rows:
        value *= R1[i] - 1

    for j in cols:
        value *= R2[j] - 1

    return value


def block_joint_information(
    satisfying: int,
    raw: int,
) -> float:

    if satisfying <= 0:
        return float("inf")

    return -math.log2(satisfying / raw)


# ========================================================================
# Cell-level information
# ========================================================================

def cell_information(
    relations,
) -> List[Tuple[float, int, int, int, int]]:

    results = []

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            allowed = 0

            for x in range(1, r):

                allowed += bit_count(
                    relations[i][j][x]
                )

            raw = (r - 1) * (s - 1)

            fraction = allowed / raw

            info = (
                float("inf")
                if allowed == 0
                else -math.log2(fraction)
            )

            results.append(
                (
                    info,
                    i,
                    j,
                    allowed,
                    raw,
                )
            )

    results.sort(reverse=True)

    return results


# ========================================================================
# Higher-order block evaluator
# ========================================================================

def evaluate_block(
    relations,
    rows: Tuple[int, ...],
    cols: Tuple[int, ...],
    true_p: int,
    true_q: int,
) -> Dict:

    raw = block_raw_state_count(
        rows,
        cols,
    )

    satisfying = 0
    surviving_row_tuples = 0

    # True p residues.
    true_row_values = tuple(
        true_p % R1[i]
        for i in rows
    )

    true_col_values = tuple(
        true_p % R2[j]
        for j in cols
    )

    true_survives = True

    # --------------------------------------------------------------------
    # Enumerate ONLY the selected R1 variables.
    #
    # For every selected row tuple:
    #
    #   each R2 variable is constrained independently
    #
    # by intersecting the masks from the selected rows.
    #
    # This is exact block composition and avoids Cartesian enumeration
    # of the complete 2m/large global CSP.
    # --------------------------------------------------------------------

    row_domains = [
        range(1, R1[i])
        for i in rows
    ]

    for row_values in itertools.product(*row_domains):

        masks = {}

        impossible = False

        for col_pos, j in enumerate(cols):

            s = R2[j]

            mask = full_mask(s)

            for row_pos, i in enumerate(rows):

                x = row_values[row_pos]

                mask &= relations[i][j][x]

                if mask == 0:
                    impossible = True
                    break

            if impossible:
                break

            masks[j] = mask

        if impossible:
            continue

        surviving_row_tuples += 1

        ways = 1

        for j in cols:
            ways *= bit_count(masks[j])

        satisfying += ways

        # Check true configuration.
        if row_values == true_row_values:

            for j in cols:

                true_y = true_p % R2[j]

                if not (
                    masks[j]
                    & (1 << (true_y - 1))
                ):
                    true_survives = False

    fraction = satisfying / raw

    info = (
        float("inf")
        if satisfying == 0
        else -math.log2(fraction)
    )

    return {
        "rows": rows,
        "cols": cols,
        "raw": raw,
        "satisfying": satisfying,
        "fraction": fraction,
        "information": info,
        "surviving_row_tuples": surviving_row_tuples,
        "true_survives": true_survives,
    }


# ========================================================================
# Expected additive cell information
# ========================================================================

def additive_cell_information(
    relations,
    C,
    rows,
    cols,
) -> float:

    total = 0.0

    for i in rows:

        for j in cols:

            r = R1[i]
            s = R2[j]

            allowed = 0

            for x in range(1, r):
                allowed += bit_count(
                    relations[i][j][x]
                )

            raw = (r - 1) * (s - 1)

            if allowed == 0:
                return float("inf")

            total += -math.log2(
                allowed / raw
            )

    return total


# ========================================================================
# Run all block shapes
# ========================================================================

def evaluate_all_blocks(
    relations,
    C,
    p,
    q,
) -> Dict[Tuple[int, int], List[Dict]]:

    results = {}

    for shape in BLOCK_SHAPES:

        nr, nc = shape

        shape_results = []

        for rows in itertools.combinations(
            range(len(R1)),
            nr,
        ):

            for cols in itertools.combinations(
                range(len(R2)),
                nc,
            ):

                result = evaluate_block(
                    relations,
                    rows,
                    cols,
                    p,
                    q,
                )

                result["additive_cell_information"] = (
                    additive_cell_information(
                        relations,
                        C,
                        rows,
                        cols,
                    )
                )

                result["synergy"] = (
                    result["information"]
                    - result["additive_cell_information"]
                )

                shape_results.append(result)

        shape_results.sort(
            key=lambda x: x["information"],
            reverse=True,
        )

        results[shape] = shape_results

    return results


# ========================================================================
# Pretty printing
# ========================================================================

def format_indices(
    rows,
    cols,
) -> str:

    r_names = [
        str(R1[i])
        for i in rows
    ]

    c_names = [
        str(R2[j])
        for j in cols
    ]

    return (
        "("
        + ",".join(r_names)
        + ") x ("
        + ",".join(c_names)
        + ")"
    )


def print_block_summary(
    results,
    shape,
):

    print()
    print(
        "------------------------------------------------------------------------"
    )

    print(
        f"TOP {TOP_K} BLOCKS FOR SHAPE {shape[0]}x{shape[1]}"
    )

    print(
        "------------------------------------------------------------------------"
    )

    print(
        "rank block                         info      "
        "sum-cell   synergy     fraction        states/raw"
    )

    print(
        "------------------------------------------------------------------------"
    )

    for rank, result in enumerate(
        results[shape][:TOP_K],
        start=1,
    ):

        block_name = format_indices(
            result["rows"],
            result["cols"],
        )

        print(
            f"{rank:4d} "
            f"{block_name:<30} "
            f"{result['information']:9.4f} "
            f"{result['additive_cell_information']:9.4f} "
            f"{result['synergy']:9.4f} "
            f"{result['fraction']:.6e} "
            f"{result['satisfying']}/{result['raw']}"
        )

        if not result["true_survives"]:
            print(
                "      WARNING: TRUE STATE DOES NOT SURVIVE"
            )


def print_shape_extrema(
    results,
    shape,
):

    data = results[shape]

    best = data[0]

    least = data[-1]

    print()
    print(
        f"SHAPE {shape[0]}x{shape[1]} SUMMARY"
    )

    print(
        f"  best information        = "
        f"{best['information']:.6f} bits"
    )

    print(
        f"  best block              = "
        f"{format_indices(best['rows'], best['cols'])}"
    )

    print(
        f"  best satisfying states  = "
        f"{best['satisfying']}"
    )

    print(
        f"  best raw states         = "
        f"{best['raw']}"
    )

    print(
        f"  best fraction           = "
        f"{best['fraction']:.6e}"
    )

    print(
        f"  best additive cell info = "
        f"{best['additive_cell_information']:.6f} bits"
    )

    print(
        f"  best synergy            = "
        f"{best['synergy']:.6f} bits"
    )

    print(
        f"  weakest information     = "
        f"{least['information']:.6f} bits"
    )


# ========================================================================
# Main experiment
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

    print(
        f"n bits = {bits}"
    )

    n, p, q = generate_balanced_semiprime(bits)

    print(f"n = {n}")
    print(f"true p = {p}")
    print(f"true q = {q}")

    print()
    print("TRUE p RESIDUES")

    for m in R1 + R2:
        print(
            f"  p mod {m} = {p % m}"
        )

    C = build_c_matrix(p, q)

    print()
    print("C MATRIX")

    for row in C:
        print(row)

    print()
    print("BUILDING HIGHER-ORDER RELATIONS")

    relations = build_relations(
        n,
        C,
    )

    print("DONE")

    # --------------------------------------------------------------------
    # Cell information baseline
    # --------------------------------------------------------------------

    cell_info = cell_information(
        relations
    )

    print()
    print(
        "TOP CELL-LEVEL INFORMATION"
    )

    print(
        "rank cell       C     info"
    )

    print(
        "------------------------------------------------------------------------"
    )

    for rank, (
        info,
        i,
        j,
        allowed,
        raw,
    ) in enumerate(
        cell_info[:TOP_K],
        start=1,
    ):

        print(
            f"{rank:4d} "
            f"({R1[i]},{R2[j]}) "
            f"C={C[i][j]} "
            f"{info:9.6f}"
        )

    # --------------------------------------------------------------------
    # Higher-order blocks
    # --------------------------------------------------------------------

    print()
    print(
        "EVALUATING HIGHER-ORDER BLOCK COMPOSITIONS"
    )

    results = evaluate_all_blocks(
        relations,
        C,
        p,
        q,
    )

    # --------------------------------------------------------------------
    # Print summaries
    # --------------------------------------------------------------------

    for shape in BLOCK_SHAPES:

        print_shape_extrema(
            results,
            shape,
        )

        print_block_summary(
            results,
            shape,
        )

    # --------------------------------------------------------------------
    # Critical comparison
    # --------------------------------------------------------------------

    print()
    print(
        "======================================================================="
    )

    print(
        "COMPOSITION TEST"
    )

    print(
        "======================================================================="
    )

    for shape in BLOCK_SHAPES:

        best = results[shape][0]

        ratio = (
            best["information"]
            / max(
                1e-12,
                best["additive_cell_information"],
            )
        )

        print()
        print(
            f"{shape[0]}x{shape[1]} block:"
        )

        print(
            f"  joint information       = "
            f"{best['information']:.6f} bits"
        )

        print(
            f"  additive cell estimate  = "
            f"{best['additive_cell_information']:.6f} bits"
        )

        print(
            f"  joint/additive ratio    = "
            f"{ratio:.6f}"
        )

        print(
            f"  synergy                 = "
            f"{best['synergy']:.6f} bits"
        )

        print(
            f"  true state survives     = "
            f"{best['true_survives']}"
        )

        print(
            f"  block                   = "
            f"{format_indices(best['rows'], best['cols'])}"
        )

    # --------------------------------------------------------------------
    # Best 3x3 block
    # --------------------------------------------------------------------

    best33 = results[(3, 3)][0]

    print()
    print(
        "======================================================================="
    )

    print(
        "BEST 3x3 BLOCK"
    )

    print(
        "======================================================================="
    )

    print(
        f"block = "
        f"{format_indices(best33['rows'], best33['cols'])}"
    )

    print(
        f"raw states = "
        f"{best33['raw']}"
    )

    print(
        f"satisfying states = "
        f"{best33['satisfying']}"
    )

    print(
        f"fraction = "
        f"{best33['fraction']:.12e}"
    )

    print(
        f"information = "
        f"{best33['information']:.6f} bits"
    )

    print(
        f"additive cell information = "
        f"{best33['additive_cell_information']:.6f} bits"
    )

    print(
        f"synergy = "
        f"{best33['synergy']:.6f} bits"
    )

    print(
        f"surviving row tuples = "
        f"{best33['surviving_row_tuples']}"
    )

    print(
        f"true state survives = "
        f"{best33['true_survives']}"
    )


# ========================================================================
# Driver
# ========================================================================

def main():

    random.seed(159)

    print("=" * 72)
    print("START EXPERIMENT 159")
    print("Higher-order radix/carry constraint composition")
    print("=" * 72)

    print()
    print(f"R1 = {R1}")
    print(f"R2 = {R2}")

    for bits in BIT_SIZES:

        print()
        print("=" * 72)
        print(f"BIT SIZE = {bits}")
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
    print("FINISHED EXPERIMENT 159")
    print("=" * 72)


if __name__ == "__main__":
    main()
