from math import comb

print("=" * 78)
print("EXPERIMENT 198")
print("EXACT pq-LEVEL SUMMAND EXTRACTION AUDIT")
print("=" * 78)
print()
print("This experiment is designed to locate the FIRST place where")
print("p/q-dependence enters the exact construction.")
print()
print("It does NOT search arbitrary closed forms.")
print("It does NOT enumerate unrestricted affine candidates.")
print()


# ============================================================================
# CONFIGURATION
# ============================================================================

# Set this to True ONLY after inserting the real exact pq-level summand
# construction in exact_pq_summands().
CONNECTED_TO_EXACT_PQ = False


# ============================================================================
# KNOWN REDUCED-STATE ANCHORS FROM EXPERIMENTS 191-194
# ============================================================================
#
# These are preserved as reduced-state anchors only.
# No p,q pair is fabricated for them.

REDUCED_ANCHORS = [
    ("A1", 1, 10, 2, 27),
    ("A2", 3, 20, 6, 935),
    ("A3", 7, 40, 15, -1797818),
    ("B1", 3, 7, 3, -3),
    ("B2", 3, 9, 4, 9),
    ("B3", 3, 11, 5, -16),
    ("B4", 5, 11, 5, -5),
    ("B5", 5, 19, 9, -196),
    ("B6", 9, 21, 11, -84),
]


# ============================================================================
# ACTUAL pq CASES ALREADY USED IN EXPERIMENT 195
# ============================================================================
#
# These are genuine p,q pairs from the earlier pq reduction audit.

PQ_CASES = [
    (3, 5),
    (3, 7),
    (3, 11),
    (3, 13),
    (3, 17),
    (3, 19),
    (3, 23),
    (3, 29),
    (3, 31),

    (5, 7),
    (5, 11),
    (5, 13),
    (5, 17),
    (5, 19),
    (5, 23),
    (5, 29),
    (5, 31),

    (7, 11),
    (7, 13),
    (7, 17),
    (7, 19),
    (7, 23),
    (7, 29),
    (7, 31),

    (11, 13),
    (11, 17),
    (11, 19),
    (11, 23),
    (11, 29),
    (11, 31),

    (13, 17),
    (13, 19),
    (13, 23),
    (13, 29),
    (13, 31),

    (17, 19),
    (17, 23),
    (17, 29),
    (17, 31),

    (19, 23),
    (19, 29),
    (19, 31),

    (23, 29),
    (23, 31),

    (29, 31),
]


# ============================================================================
# KNOWN EXPERIMENT-195 TEMPORARY REDUCTION MAP
# ============================================================================
#
# This is included ONLY so the script can reproduce the old diagnostic.
# It is NOT claimed to be the theoretical pq -> (k,ell,s) map.
#
# Observed from Experiment 195:
#
# p=3  -> k=1
# p=5  -> k=3
# p=7  -> k=5
# p=11 -> k=9
# etc.
#
# For the q values used there, the old reduction produced:
# ell = 2*floor(q/2) +  ?  etc.
#
# Rather than reconstruct that heuristic here, we preserve the actual
# reduction values explicitly for the pq cases from Experiment 195.

REDUCED_MAP = {
    (3, 5): (1, 8, 2),
    (3, 7): (1, 10, 3),
    (3, 11): (1, 14, 5),
    (3, 13): (1, 16, 6),
    (3, 17): (1, 20, 8),
    (3, 19): (1, 22, 9),
    (3, 23): (1, 26, 11),
    (3, 29): (1, 32, 14),
    (3, 31): (1, 34, 15),

    (5, 7): (3, 12, 3),
    (5, 11): (3, 16, 5),
    (5, 13): (3, 18, 6),
    (5, 17): (3, 22, 8),
    (5, 19): (3, 24, 9),
    (5, 23): (3, 28, 11),
    (5, 29): (3, 34, 14),
    (5, 31): (3, 36, 15),

    (7, 11): (5, 18, 5),
    (7, 13): (5, 20, 6),
    (7, 17): (5, 24, 8),
    (7, 19): (5, 26, 9),
    (7, 23): (5, 30, 11),
    (7, 29): (5, 36, 14),
    (7, 31): (5, 38, 15),

    (11, 13): (9, 24, 6),
    (11, 17): (9, 28, 8),
    (11, 19): (9, 30, 9),
    (11, 23): (9, 34, 11),
    (11, 29): (9, 40, 14),
    (11, 31): (9, 42, 15),

    (13, 17): (11, 30, 8),
    (13, 19): (11, 32, 9),
    (13, 23): (11, 36, 11),
    (13, 29): (11, 42, 14),
    (13, 31): (11, 44, 15),

    (17, 19): (15, 36, 9),
    (17, 23): (15, 40, 11),
    (17, 29): (15, 46, 14),
    (17, 31): (15, 48, 15),

    (19, 23): (17, 42, 11),
    (19, 29): (17, 48, 14),
    (19, 31): (17, 50, 15),

    (23, 29): (21, 52, 14),
    (23, 31): (21, 54, 15),

    (29, 31): (27, 60, 15),
}


# ============================================================================
# BASIC HELPERS
# ============================================================================

def C(n, r):
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def k_from_p(p):
    return (p - 1) // 2


def reduced_state_from_map(p, q):
    return REDUCED_MAP.get((p, q))


# ============================================================================
# REFERENCE CONVOLUTION
# ============================================================================

def reference_convolution(k, ell, s):
    total = 0

    for j in range(s + 1):
        total += (
            (-1) ** j
            * C(ell, j)
            * C(ell - j, s - j)
        )

    return total


# ============================================================================
# EXACT pq-LEVEL SUMMAND FUNCTION
# ============================================================================
#
# THIS IS THE ONLY PART THAT MUST EVENTUALLY BE CONNECTED TO THE REAL
# MATHEMATICAL pq CONSTRUCTION.
#
# The important point is that we do NOT pretend the old reduced convolution
# is the pq expression.
#
# When connected, return:
#
# [
#     {
#         "index": j,
#         "value": term_value,
#         "factors": {
#             "sign": ...,
#             "factor_1": ...,
#             "factor_2": ...,
#             ...
#         }
#     },
#     ...
# ]
#
# Do not collapse the sum here.

def exact_pq_summands(p, q):
    if not CONNECTED_TO_EXACT_PQ:
        return None

    # ========================================================================
    # REPLACE THIS BODY WITH THE CURRENT EXACT n=pq SUMMAND CONSTRUCTION
    # ========================================================================
    #
    # Example:
    #
    # terms = []
    #
    # for j in range(j_min, j_max + 1):
    #
    #     factor_a = ...
    #     factor_b = ...
    #     sign = ...
    #
    #     value = sign * factor_a * factor_b
    #
    #     terms.append({
    #         "index": j,
    #         "value": value,
    #         "factors": {
    #             "sign": sign,
    #             "factor_a": factor_a,
    #             "factor_b": factor_b,
    #         }
    #     })
    #
    # return terms
    #
    # ========================================================================

    raise RuntimeError(
        "CONNECTED_TO_EXACT_PQ=True but exact_pq_summands() "
        "has not been connected to the real pq construction."
    )


# ============================================================================
# 1. REDUCED ANCHOR BASELINE
# ============================================================================

print("=" * 78)
print("1. REDUCED ANCHOR BASELINE")
print("=" * 78)

for label, k, ell, s, expected in REDUCED_ANCHORS:
    current = reference_convolution(k, ell, s)

    status = "PASS" if current == expected else "FAIL"

    print(
        f"{label}: "
        f"k={k:2d} ell={ell:3d} s={s:2d} "
        f"expected={expected:12d} "
        f"reference={current:12d} "
        f"{status}"
    )


# ============================================================================
# 2. pq -> REDUCED-STATE AUDIT
# ============================================================================

print()
print("=" * 78)
print("2. pq -> REDUCED-STATE AUDIT")
print("=" * 78)

for p, q in PQ_CASES:

    state = reduced_state_from_map(p, q)

    if state is None:
        print(
            f"p={p:2d} q={q:2d}: REDUCED MAP MISSING"
        )
        continue

    k, ell, s = state
    value = reference_convolution(k, ell, s)

    print(
        f"p={p:2d} q={q:2d} "
        f"n={p*q:4d} "
        f"-> k={k:2d} ell={ell:3d} s={s:2d} "
        f"reference={value:10d}"
    )


# ============================================================================
# 3. EXACT PQ CONNECTION STATUS
# ============================================================================

print()
print("=" * 78)
print("3. EXACT pq CONNECTION STATUS")
print("=" * 78)

if not CONNECTED_TO_EXACT_PQ:
    print()
    print("STATUS: EXACT pq SUMMANDS NOT YET CONNECTED")
    print()
    print("This run is intentionally stopping before making any")
    print("mathematical claim about the pq-level summands.")
    print()
    print("The next required insertion is the real exact n=pq")
    print("summand construction inside exact_pq_summands().")
else:
    print()
    print("STATUS: EXACT pq SUMMANDS CONNECTED")


# ============================================================================
# 4. EXACT TERM EXTRACTION, IF CONNECTED
# ============================================================================

if CONNECTED_TO_EXACT_PQ:

    print()
    print("=" * 78)
    print("4. EXACT pq-LEVEL TERM EXTRACTION")
    print("=" * 78)

    for p, q in PQ_CASES:

        terms = exact_pq_summands(p, q)

        if terms is None:
            continue

        total = sum(t["value"] for t in terms)

        print()
        print(
            f"p={p:2d} q={q:2d} n={p*q:4d}"
        )
        print(
            f"number of terms = {len(terms)}"
        )
        print(
            f"exact pq sum    = {total}"
        )

        for t in terms:
            print(
                f"  j={t['index']:3d} "
                f"value={t['value']:12d} "
                f"factors={t.get('factors', {})}"
            )


# ============================================================================
# 5. TERM-LEVEL k DEPENDENCE
# ============================================================================

if CONNECTED_TO_EXACT_PQ:

    print()
    print("=" * 78)
    print("5. TERM-LEVEL k-DEPENDENCE")
    print("=" * 78)

    # Group by q so p changes while q is held fixed.
    by_q = {}

    for p, q in PQ_CASES:
        by_q.setdefault(q, []).append((p, q))

    for q, cases in sorted(by_q.items()):

        if len(cases) < 2:
            continue

        print()
        print(f"q={q}")

        signatures = []

        for p, _ in cases:

            terms = exact_pq_summands(p, q)

            signature = tuple(
                (
                    t["index"],
                    t["value"],
                    tuple(
                        sorted(
                            t.get("factors", {}).items()
                        )
                    )
                )
                for t in terms
            )

            signatures.append(signature)

            print(
                f"  p={p:2d} "
                f"k={k_from_p(p):2d} "
                f"terms={len(terms):3d}"
            )

        identical = len(set(signatures)) == 1

        print(
            f"  exact summand structure identical = "
            f"{identical}"
        )


# ============================================================================
# 6. FIRST k-DEPENDENT COMPONENT
# ============================================================================

if CONNECTED_TO_EXACT_PQ:

    print()
    print("=" * 78)
    print("6. FIRST k-DEPENDENT COMPONENT")
    print("=" * 78)

    print()
    print(
        "For each pq summand, inspect factor labels and identify"
    )
    print(
        "which factor changes when p changes."
    )

    for q, cases in sorted(by_q.items()):

        if len(cases) < 2:
            continue

        p0, _ = cases[0]
        p1, _ = cases[1]

        t0 = exact_pq_summands(p0, q)
        t1 = exact_pq_summands(p1, q)

        print()
        print(
            f"q={q}: compare p={p0} and p={p1}"
        )

        max_len = max(len(t0), len(t1))

        for i in range(max_len):

            if i >= len(t0) or i >= len(t1):
                print(
                    f"  TERM COUNT / RANGE CHANGES at position {i}"
                )
                continue

            a = t0[i]
            b = t1[i]

            changed = []

            if a["index"] != b["index"]:
                changed.append("index")

            if a["value"] != b["value"]:
                changed.append("value")

            fa = a.get("factors", {})
            fb = b.get("factors", {})

            all_keys = sorted(set(fa) | set(fb))

            for key in all_keys:
                if fa.get(key) != fb.get(key):
                    changed.append(key)

            print(
                f"  term {i:3d}: "
                f"changed={changed}"
            )


# ============================================================================
# 7. COLLAPSE-LOSS COMPARISON
# ============================================================================

if CONNECTED_TO_EXACT_PQ:

    print()
    print("=" * 78)
    print("7. COLLAPSE-LOSS COMPARISON")
    print("=" * 78)

    for p, q in PQ_CASES:

        state = reduced_state_from_map(p, q)

        if state is None:
            continue

        k, ell, s = state

        pq_terms = exact_pq_summands(p, q)

        pq_total = sum(
            t["value"]
            for t in pq_terms
        )

        reduced_total = reference_convolution(
            k,
            ell,
            s
        )

        print(
            f"p={p:2d} q={q:2d} "
            f"pq_exact={pq_total:12d} "
            f"reduced={reduced_total:12d} "
            f"difference={pq_total - reduced_total:12d}"
        )


# ============================================================================
# 8. FINAL DIAGNOSTIC
# ============================================================================

print()
print("=" * 78)
print("FINAL DIAGNOSTIC")
print("=" * 78)
print()

print("The logical target is:")
print()
print("    exact n=pq construction")
print("              |")
print("              v")
print("      exact summand T_j(p,q)")
print("              |")
print("              v")
print("      admissible j-range")
print("              |")
print("              v")
print("      reduction to (k,ell,s)")
print("              |")
print("              v")
print("         convolution")
print()

if not CONNECTED_TO_EXACT_PQ:
    print("CURRENT STATUS:")
    print()
    print("  The reduced convolution is still the only connected object.")
    print("  No pq-level mathematical conclusion is being claimed.")
    print()
    print("NEXT STEP:")
    print()
    print("  Insert the actual exact n=pq summand formula into")
    print("  exact_pq_summands().")
    print()
    print("Then set:")
    print("  CONNECTED_TO_EXACT_PQ = True")
    print()
else:
    print("CURRENT STATUS:")
    print()
    print("  Exact pq-level summands are connected.")
    print("  Inspect the first factor/index/range change across p.")
    print()

print("=" * 78)
print("END EXPERIMENT 198")
print("=" * 78)