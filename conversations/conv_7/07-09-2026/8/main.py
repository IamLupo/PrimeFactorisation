from collections import defaultdict
from sympy import factorint


# ============================================================
# SETTINGS
# ============================================================

TEST_CASES = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
    (37, 71),
    (41, 73),
    (43, 79),
]

R1 = 11
R2 = 13

# Hard output limits
MAX_OUTPUT_LINES = 2000
MAX_RAW_STATES_PER_CASE = 100


# ============================================================
# OUTPUT CONTROL
# ============================================================

output_lines = 0


def emit(text=""):
    """
    Print while enforcing the global output limit.
    """
    global output_lines

    if output_lines >= MAX_OUTPUT_LINES:
        return False

    print(text)
    output_lines += 1
    return True


# ============================================================
# MACMAHON M2
# ============================================================

def macmahon_M2_terms(n):
    """
    Generate all M2 states

        m1*s1 + m2*s2 = n

    with

        0 < s1 < s2
        m1,m2 > 0
    """

    for s1 in range(1, n):

        max_m1 = n // s1

        for s2 in range(s1 + 1, n + 1):

            for m1 in range(1, max_m1 + 1):

                remainder = n - m1 * s1

                if remainder <= 0:
                    continue

                if remainder % s2 != 0:
                    continue

                m2 = remainder // s2

                if m2 <= 0:
                    continue

                yield {
                    "s1": s1,
                    "s2": s2,
                    "m1": m1,
                    "m2": m2,
                    "weight": m1 * m2,
                }


# ============================================================
# QUOTIENT / RESIDUE
# ============================================================

def qr_coordinates(x, r):
    """
    x = a + k*r
    """

    return x % r, x // r


# ============================================================
# FACTOR COORDINATES
# ============================================================

def factor_coordinates(p, q, r1, r2):

    a, k = qr_coordinates(p, r1)
    b, l = qr_coordinates(q, r2)

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
# START
# ============================================================

emit("START EXPERIMENT 5")
emit()

for p, q in TEST_CASES:

    n = p * q

    coords = factor_coordinates(p, q, R1, R2)

    # --------------------------------------------------------
    # Generate states ONCE
    # --------------------------------------------------------

    terms = list(macmahon_M2_terms(n))

    M2 = sum(t["weight"] for t in terms)

    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    emit("=" * 80)
    emit(f"n = {n}")
    emit(f"factorization = {factorint(n)}")

    emit(
        f"factor coordinates: "
        f"p={p} -> a={coords['a']}, k={coords['k']}; "
        f"q={q} -> b={coords['b']}, l={coords['l']}"
    )

    emit(
        f"S={coords['S']}  "
        f"K={coords['K']}"
    )

    emit(f"M2(n) = {M2}")
    emit(f"number of M2 states = {len(terms)}")

    # ========================================================
    # RAW SAMPLE ONLY
    # ========================================================

    emit()
    emit("RAW M2 SAMPLE")
    emit(
        "s1 s2 m1 m2 weight | "
        "s1%r1 s2%r1 | "
        "m1%r1 m2%r1 | "
        "s1%r2 s2%r2 | "
        "m1%r2 m2%r2"
    )

    raw_count = 0

    for t in terms:

        if raw_count >= MAX_RAW_STATES_PER_CASE:
            break

        s1 = t["s1"]
        s2 = t["s2"]
        m1 = t["m1"]
        m2 = t["m2"]
        w = t["weight"]

        emit(
            f"{s1} {s2} {m1} {m2} {w} | "
            f"{s1 % R1} {s2 % R1} | "
            f"{m1 % R1} {m2 % R1} | "
            f"{s1 % R2} {s2 % R2} | "
            f"{m1 % R2} {m2 % R2}"
        )

        raw_count += 1

    if len(terms) > MAX_RAW_STATES_PER_CASE:
        emit(
            f"... "
            f"{len(terms) - MAX_RAW_STATES_PER_CASE} "
            f"additional states omitted ..."
        )

    # ========================================================
    # AGGREGATE A
    # PART-SIZE RESIDUES
    # ========================================================

    part_residue = defaultdict(lambda: {
        "states": 0,
        "weight": 0,
    })

    for t in terms:

        key = (
            t["s1"] % R1,
            t["s2"] % R1,
            t["s1"] % R2,
            t["s2"] % R2,
        )

        part_residue[key]["states"] += 1
        part_residue[key]["weight"] += t["weight"]

    emit()
    emit("AGGREGATE A: PART-SIZE RESIDUES")
    emit(
        "s1_r1 s2_r1 s1_r2 s2_r2 states weight"
    )

    for key in sorted(part_residue):

        if output_lines >= MAX_OUTPUT_LINES:
            break

        s1r1, s2r1, s1r2, s2r2 = key
        d = part_residue[key]

        emit(
            f"{s1r1} {s2r1} "
            f"{s1r2} {s2r2} "
            f"{d['states']} "
            f"{d['weight']}"
        )

    # ========================================================
    # AGGREGATE B
    # MULTIPLICITY RESIDUES
    # ========================================================

    multiplicity_residue = defaultdict(lambda: {
        "states": 0,
        "weight": 0,
    })

    for t in terms:

        key = (
            t["m1"] % R1,
            t["m2"] % R1,
            t["m1"] % R2,
            t["m2"] % R2,
        )

        multiplicity_residue[key]["states"] += 1
        multiplicity_residue[key]["weight"] += t["weight"]

    emit()
    emit("AGGREGATE B: MULTIPLICITY RESIDUES")
    emit(
        "m1_r1 m2_r1 m1_r2 m2_r2 states weight"
    )

    for key in sorted(multiplicity_residue):

        if output_lines >= MAX_OUTPUT_LINES:
            break

        m1r1, m2r1, m1r2, m2r2 = key
        d = multiplicity_residue[key]

        emit(
            f"{m1r1} {m2r1} "
            f"{m1r2} {m2r2} "
            f"{d['states']} "
            f"{d['weight']}"
        )

    # ========================================================
    # AGGREGATE C
    # PART-SIZE QUOTIENT CELLS
    # ========================================================

    part_quotient = defaultdict(lambda: {
        "states": 0,
        "weight": 0,
    })

    for t in terms:

        s1 = t["s1"]
        s2 = t["s2"]

        _, k1 = qr_coordinates(s1, R1)
        _, k2 = qr_coordinates(s2, R2)

        key = (k1, k2)

        part_quotient[key]["states"] += 1
        part_quotient[key]["weight"] += t["weight"]

    emit()
    emit("AGGREGATE C: PART-SIZE QUOTIENT CELLS")
    emit("k1 k2 states weight")

    for key in sorted(part_quotient):

        if output_lines >= MAX_OUTPUT_LINES:
            break

        k1, k2 = key
        d = part_quotient[key]

        emit(
            f"{k1} {k2} "
            f"{d['states']} "
            f"{d['weight']}"
        )

    # ========================================================
    # AGGREGATE D
    # MULTIPLICITY QUOTIENT CELLS
    # ========================================================

    multiplicity_quotient = defaultdict(lambda: {
        "states": 0,
        "weight": 0,
    })

    for t in terms:

        m1 = t["m1"]
        m2 = t["m2"]

        _, k1 = qr_coordinates(m1, R1)
        _, k2 = qr_coordinates(m2, R2)

        key = (k1, k2)

        multiplicity_quotient[key]["states"] += 1
        multiplicity_quotient[key]["weight"] += t["weight"]

    emit()
    emit("AGGREGATE D: MULTIPLICITY QUOTIENT CELLS")
    emit("k1 k2 states weight")

    for key in sorted(multiplicity_quotient):

        if output_lines >= MAX_OUTPUT_LINES:
            break

        k1, k2 = key
        d = multiplicity_quotient[key]

        emit(
            f"{k1} {k2} "
            f"{d['states']} "
            f"{d['weight']}"
        )

    # ========================================================
    # SPECIAL CHECK
    # ========================================================

    matching_states = []

    for t in terms:

        s1 = t["s1"]
        s2 = t["s2"]
        m1 = t["m1"]
        m2 = t["m2"]

        if (
            p in (s1, s2, m1, m2)
            or
            q in (s1, s2, m1, m2)
        ):
            matching_states.append(t)

    emit()
    emit(
        "SPECIAL CHECK: STATES CONTAINING p OR q"
    )
    emit(
        f"matching states = {len(matching_states)}"
    )

    # Only print the first few matching states.
    for t in matching_states[:20]:

        emit(
            f"{t['s1']} "
            f"{t['s2']} "
            f"{t['m1']} "
            f"{t['m2']} "
            f"{t['weight']}"
        )

    # ========================================================
    # TARGET
    # ========================================================

    emit()
    emit("FACTOR COORDINATE TARGET")
    emit(f"p = {p}")
    emit(f"q = {q}")
    emit(f"a = {coords['a']}")
    emit(f"b = {coords['b']}")
    emit(f"k = {coords['k']}")
    emit(f"l = {coords['l']}")
    emit(f"K = {coords['K']}")
    emit(f"S = {coords['S']}")

    if output_lines >= MAX_OUTPUT_LINES:
        break

    emit()


# ============================================================
# END
# ============================================================

emit(
    f"Output lines used: "
    f"{output_lines}/{MAX_OUTPUT_LINES}"
)

emit("FINISHED EXPERIMENT 5")