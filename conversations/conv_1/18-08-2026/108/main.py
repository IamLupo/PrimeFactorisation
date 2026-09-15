"""
==============================================================================
EXPERIMENT 111 — EXACT MINIMAL CHANNEL SUPPORT LAW / BRANCH-THRESHOLD AUDIT
==============================================================================

Observed support:

    A: c=0
    B: c=1

Candidate unified law:

    s(p,c) = max(p-a, floor((p+b*c+d)/2))

and the corresponding ceil form.

Goals:

    1. Find the minimal exact integer parameters.
    2. Check uniqueness/minimality.
    3. Compare floor and ceil representations.
    4. Verify the compressed law against every observed point.
    5. Determine where the two branches exchange dominance.
    6. Check whether max(0, ...) is actually unnecessary.

Exact integer arithmetic only.
No floating point.
No extrapolation.
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass


# ============================================================================
# OBSERVED DATA
# ============================================================================

DATA = {
    "A-even": {
        "c": 0,
        "points": [(0, 0), (2, 1), (4, 2), (6, 4), (8, 6)],
    },
    "A-odd": {
        "c": 0,
        "points": [(1, 1), (3, 2), (5, 3), (7, 5)],
    },
    "B-even": {
        "c": 1,
        "points": [(0, 0), (2, 1), (4, 2), (6, 4)],
    },
    "B-odd": {
        "c": 1,
        "points": [(1, 0), (3, 1), (5, 3), (7, 5)],
    },
}


# ============================================================================
# INTEGER HELPERS
# ============================================================================

def floor_div(a: int, b: int) -> int:
    return a // b


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


# ============================================================================
# CANDIDATE
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    mode: str
    a: int
    b: int
    d: int

    def value(self, p: int, c: int) -> int:
        linear = p - self.a
        numerator = p + self.b * c + self.d

        if self.mode == "floor":
            half = floor_div(numerator, 2)
        elif self.mode == "ceil":
            half = ceil_div(numerator, 2)
        else:
            raise ValueError(self.mode)

        return max(linear, half)

    def value_with_zero(self, p: int, c: int) -> int:
        return max(0, self.value(p, c))

    def complexity(self) -> tuple[int, int, int, str]:
        """
        Prefer:
          1. smallest total coefficient magnitude;
          2. smallest number of nonzero corrections;
          3. shortest description;
          4. lexicographic stability.
        """
        magnitude = abs(self.a) + abs(self.b) + abs(self.d)
        nonzero = sum(
            int(x != 0)
            for x in (self.a, self.b, self.d)
        )
        return (
            magnitude,
            nonzero,
            len(self.text()),
            self.text(),
        )

    def text(self) -> str:
        if self.a == 0:
            left = "p"
        elif self.a > 0:
            left = f"p-{self.a}"
        else:
            left = f"p+{-self.a}"

        terms = ["p"]

        if self.b != 0:
            if self.b == 1:
                terms.append("+c")
            elif self.b == -1:
                terms.append("-c")
            elif self.b > 0:
                terms.append(f"+{self.b}c")
            else:
                terms.append(f"{self.b}c")

        if self.d != 0:
            if self.d > 0:
                terms.append(f"+{self.d}")
            else:
                terms.append(str(self.d))

        numerator = "".join(terms)

        fn = "floor" if self.mode == "floor" else "ceil"

        return (
            f"max({left}, {fn}(({numerator})/2))"
        )


# ============================================================================
# EXACT VALIDATION
# ============================================================================

def exact(candidate: Candidate, include_zero: bool = False) -> bool:
    for block in DATA.values():
        c = block["c"]

        for p, observed in block["points"]:
            predicted = (
                candidate.value_with_zero(p, c)
                if include_zero
                else candidate.value(p, c)
            )

            if predicted != observed:
                return False

    return True


# ============================================================================
# SEARCH
# ============================================================================

def search(limit: int = 6) -> list[Candidate]:
    found: list[Candidate] = []

    for mode in ("floor", "ceil"):
        for a in range(-limit, limit + 1):
            for b in range(-limit, limit + 1):
                for d in range(-limit, limit + 1):
                    candidate = Candidate(
                        mode=mode,
                        a=a,
                        b=b,
                        d=d,
                    )

                    if exact(candidate):
                        found.append(candidate)

    found.sort(key=lambda x: x.complexity())
    return found


# ============================================================================
# DIRECT KNOWN LAW
# ============================================================================

def known_law(p: int, c: int) -> int:
    return max(
        p - 2,
        floor_div(p + 1 - c, 2),
    )


def known_law_exact() -> bool:
    for block in DATA.values():
        c = block["c"]

        for p, observed in block["points"]:
            if known_law(p, c) != observed:
                return False

    return True


# ============================================================================
# EQUIVALENT REPRESENTATIONS
# ============================================================================

def equivalent_floor_ceil(p: int, c: int) -> bool:
    left = ceil_div(p - c, 2)
    right = floor_div(p + 1 - c, 2)

    return left == right


def all_equivalent_forms_exact() -> bool:
    for p in range(0, 50):
        for c in (0, 1):
            if not equivalent_floor_ceil(p, c):
                return False

    return True


# ============================================================================
# ZERO TERM REDUNDANCY
# ============================================================================

def zero_redundant() -> bool:
    """
    Check whether the known law is already nonnegative on the entire
    nonnegative integer p-domain for c in {0,1}.
    """
    for p in range(0, 100):
        for c in (0, 1):
            if known_law(p, c) < 0:
                return False

    return True


# ============================================================================
# BRANCH COMPARISON
# ============================================================================

def branch_values(p: int, c: int) -> tuple[int, int]:
    left = p - 2
    right = floor_div(p + 1 - c, 2)
    return left, right


def dominant_branch(p: int, c: int) -> str:
    left, right = branch_values(p, c)

    if left > right:
        return "linear"
    if right > left:
        return "half"
    return "tie"


def print_branch_table() -> None:
    print()
    print("  p  c   p-2   floor((p+1-c)/2)   dominant")
    print("  " + "-" * 48)

    for c in (0, 1):
        for p in range(0, 11):
            left, right = branch_values(p, c)
            dom = dominant_branch(p, c)

            print(
                f"  {p:2d} {c:2d} "
                f"{left:5d} "
                f"{right:18d}   "
                f"{dom}"
            )


# ============================================================================
# SYMBOLIC SWITCH ANALYSIS
# ============================================================================

def threshold_audit() -> None:
    print()
    print("  Threshold condition:")
    print()
    print("      p - 2 >= floor((p+1-c)/2)")
    print()
    print("  Exact observed transition:")
    print()

    for c in (0, 1):
        print(f"  channel c={c}")

        last_half = None

        for p in range(0, 12):
            dom = dominant_branch(p, c)

            if dom == "half":
                last_half = p

            print(
                f"    p={p}: {dom}"
            )

        print(
            f"    last half-dominant p = {last_half}"
        )


# ============================================================================
# VALIDATION REPORT
# ============================================================================

def validate_candidate(candidate: Candidate) -> bool:
    all_ok = True

    for name, block in DATA.items():
        c = block["c"]

        print(f"  {name}")

        for p, observed in block["points"]:
            predicted = candidate.value(p, c)
            ok = predicted == observed
            all_ok = all_ok and ok

            print(
                f"    p={p}: "
                f"observed={observed} "
                f"predicted={predicted} "
                f"exact={ok}"
            )

    return all_ok


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print(
        "EXPERIMENT 111 — EXACT MINIMAL CHANNEL SUPPORT LAW / "
        "BRANCH-THRESHOLD AUDIT"
    )
    print("=" * 78)

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("1. OBSERVED SUPPORT")
    print("=" * 78)

    for name, block in DATA.items():
        print(
            f"  {name}: c={block['c']}"
        )

        for p, s in block["points"]:
            print(
                f"    p={p}: s={s}"
            )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("2. EXACT SEARCH")
    print("=" * 78)

    candidates = search(limit=6)

    print(
        f"  exact candidate count = {len(candidates)}"
    )

    for candidate in candidates[:30]:
        print(
            f"    {candidate.text()}"
            f"    score={candidate.complexity()}"
        )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("3. MINIMAL EXACT CANDIDATE")
    print("=" * 78)

    if candidates:
        best = candidates[0]

        print(
            f"  best = {best.text()}"
        )

        print(
            f"  mode={best.mode}"
        )
        print(
            f"  a={best.a}"
        )
        print(
            f"  b={best.b}"
        )
        print(
            f"  d={best.d}"
        )
        print(
            f"  complexity={best.complexity()}"
        )
    else:
        best = None
        print(
            "  no exact candidate found"
        )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("4. KNOWN COMPRESSED LAW")
    print("=" * 78)

    print(
        "  s(p,c) = max(p-2, floor((p+1-c)/2))"
    )

    print(
        f"  exact={known_law_exact()}"
    )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("5. EQUIVALENT CEILING FORM")
    print("=" * 78)

    print(
        "  ceil((p-c)/2) == floor((p+1-c)/2)"
    )

    print(
        f"  exact_for_p>=0,c={{{{0,1}}}} = "
        f"{all_equivalent_forms_exact()}"
    )

    print()
    print(
        "  Therefore the equivalent form is:"
    )
    print(
        "      s(p,c) = max(p-2, ceil((p-c)/2))"
    )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("6. ZERO-TERM REDUNDANCY")
    print("=" * 78)

    print(
        "  Is max(0, ...) necessary?"
    )
    print(
        f"  answer = {not zero_redundant()}"
    )

    print(
        "  zero term redundant on p>=0,c in {0,1} = "
        f"{zero_redundant()}"
    )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("7. BRANCH TABLE")
    print("=" * 78)

    print_branch_table()

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("8. BRANCH-THRESHOLD AUDIT")
    print("=" * 78)

    threshold_audit()

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("9. BEST-CANDIDATE VALIDATION")
    print("=" * 78)

    best_ok = False

    if best is not None:
        best_ok = validate_candidate(best)

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("10. SECTOR EQUIVALENCE")
    print("=" * 78)

    sector_forms = {
        "A-even": lambda p: max(
            p - 2,
            floor_div(p, 2),
        ),
        "A-odd": lambda p: max(
            p - 2,
            ceil_div(p, 2),
        ),
        "B-even": lambda p: max(
            p - 2,
            floor_div(p, 2),
        ),
        "B-odd": lambda p: max(
            p - 2,
            floor_div(p - 1, 2),
        ),
    }

    for name, block in DATA.items():
        c = block["c"]

        eq = True

        for p, observed in block["points"]:
            common = known_law(p, c)
            sector = sector_forms[name](p)

            if common != sector:
                eq = False

        print(
            f"  {name}: compressed_equals_sector={eq}"
        )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 110 found the minimal-looking channel interaction

      c = 0  for A
      c = 1  for B

  and an exact candidate

      s(p,c)
        = max(
            p - 2,
            floor((p + 1 - c)/2)
          ).

  Because

      ceil((p-c)/2)
        = floor((p+1-c)/2)

  for integer p,c, the equivalent form is

      s(p,c)
        = max(
            p - 2,
            ceil((p-c)/2)
          ).

  This reproduces all four observed sectors:

      c=0:
          max(p-2, ceil(p/2))
          -> A-even/A-odd according to parity.

      c=1:
          max(p-2, ceil((p-1)/2))
          -> B-even/B-odd according to parity.

  The key question now is no longer support interpolation.

  It is whether this compressed law has a natural interpretation as
  the maximum of two competing endpoint constraints:

      linear branch:
          p - 2

      half-index branch:
          ceil((p-c)/2).

  The branch-threshold audit identifies precisely where the finite
  support switches from one mechanism to the other.

  Everything is exact over Z.
  No floating point.
  No recurrence search.
  No extrapolation.
        """
    .strip())

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    direct_ok = known_law_exact()
    equivalent_ok = all_equivalent_forms_exact()

    failures = 0

    if not direct_ok:
        failures += 1

    if not equivalent_ok:
        failures += 1

    if best is not None and not best_ok:
        failures += 1

    print(
        f"  minimal_search_found={bool(candidates)}"
    )
    print(
        f"  known_law_exact={direct_ok}"
    )
    print(
        f"  floor_ceil_equivalence={equivalent_ok}"
    )
    print(
        f"  best_candidate_exact={best_ok}"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 111 COMPLETE")


if __name__ == "__main__":
    main()

