"""
==============================================================================
EXPERIMENT 110 — EXACT CHANNEL/PARITY SUPPORT-LAW MINIMIZATION
==============================================================================

PURPOSE

Experiment 109 showed that m=8/7 plus parity does not produce a common
affine-delta law of the tested form.

Here we introduce the smallest additional discrete variable:

    c = 0  for A
    c = 1  for B

and parity

    e = 0  even
    e = 1  odd.

We search exact laws of the form

    s(p,c,e) =
        max(
            0,
            p - a(c,e),
            floor((p + d(c,e))/2)
        )

and the corresponding ceil version.

The affine functions are

    a(c,e) = a0 + a1*c + a2*e + a3*c*e
    d(c,e) = d0 + d1*c + d2*e + d3*c*e

This is deliberately small.

No floating point.
No recurrence search.
No extrapolation.
==============================================================================

OBSERVED DATA

A-even:
    p = 0,2,4,6,8
    s = 0,1,2,4,6

A-odd:
    p = 1,3,5,7
    s = 1,2,3,5

B-even:
    p = 0,2,4,6
    s = 0,1,2,4

B-odd:
    p = 1,3,5,7
    s = 0,1,3,5

==============================================================================
"""

from __future__ import annotations

from itertools import product
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Exact observed data
# ---------------------------------------------------------------------------

DATA = [
    # name, channel c, parity e, points
    (
        "A-even",
        0,
        0,
        [(0, 0), (2, 1), (4, 2), (6, 4), (8, 6)],
    ),
    (
        "A-odd",
        0,
        1,
        [(1, 1), (3, 2), (5, 3), (7, 5)],
    ),
    (
        "B-even",
        1,
        0,
        [(0, 0), (2, 1), (4, 2), (6, 4)],
    ),
    (
        "B-odd",
        1,
        1,
        [(1, 0), (3, 1), (5, 3), (7, 5)],
    ),
]


# ---------------------------------------------------------------------------
# Exact integer floor / ceil division
# ---------------------------------------------------------------------------

def floor_div(a: int, b: int) -> int:
    return a // b


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def max3(a: int, b: int, c: int) -> int:
    return max(a, b, c)


# ---------------------------------------------------------------------------
# Affine function in channel/parity
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AffineCE:
    c0: int
    c1: int
    c2: int
    c3: int

    def value(self, c: int, e: int) -> int:
        return (
            self.c0
            + self.c1 * c
            + self.c2 * e
            + self.c3 * c * e
        )

    def complexity(self) -> int:
        return (
            abs(self.c0)
            + abs(self.c1)
            + abs(self.c2)
            + abs(self.c3)
        )

    def text(self, symbol: str) -> str:
        parts = []

        if self.c0 != 0:
            parts.append(str(self.c0))

        if self.c1 != 0:
            parts.append(
                f"{self.c1:+d}c"
            )

        if self.c2 != 0:
            parts.append(
                f"{self.c2:+d}e"
            )

        if self.c3 != 0:
            parts.append(
                f"{self.c3:+d}ce"
            )

        if not parts:
            return "0"

        out = parts[0]

        for part in parts[1:]:
            if part.startswith("+"):
                out += " + " + part[1:]
            elif part.startswith("-"):
                out += " - " + part[1:]
            else:
                out += " + " + part

        if symbol:
            return out

        return out


# ---------------------------------------------------------------------------
# Candidate law
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Candidate:
    mode: str
    a: AffineCE
    d: AffineCE

    def evaluate(self, p: int, c: int, e: int) -> int:
        a = self.a.value(c, e)
        d = self.d.value(c, e)

        left = p - a

        if self.mode == "floor":
            half = floor_div(p + d, 2)
        elif self.mode == "ceil":
            half = ceil_div(p + d, 2)
        else:
            raise ValueError(self.mode)

        return max3(0, left, half)

    def complexity(self) -> tuple[int, int, int]:
        """
        Smaller is preferred.

        First:
            total affine complexity

        Second:
            number of interaction terms ce actually used

        Third:
            total textual complexity
        """
        interaction_count = (
            int(self.a.c3 != 0) +
            int(self.d.c3 != 0)
        )

        total = self.a.complexity() + self.d.complexity()

        text_len = len(self.text())

        return (
            total,
            interaction_count,
            text_len,
        )

    def text(self) -> str:
        mode_text = (
            "floor((p + D(c,e))/2)"
            if self.mode == "floor"
            else "ceil((p + D(c,e))/2)"
        )

        return (
            "max(0, "
            f"p - ({self.a.text('') }), "
            f"{mode_text.replace('D(c,e)', self.d.text(''))}"
            ")"
        )


# ---------------------------------------------------------------------------
# Search space
# ---------------------------------------------------------------------------

COEFF_RANGE = range(-4, 5)


def all_affines() -> list[AffineCE]:
    out = []

    for coeffs in product(
        COEFF_RANGE,
        repeat=4,
    ):
        a = AffineCE(*coeffs)

        # Keep the search finite but not huge.
        if a.complexity() <= 6:
            out.append(a)

    return out


AFFINES = all_affines()


# ---------------------------------------------------------------------------
# Exact validation
# ---------------------------------------------------------------------------

def is_exact(candidate: Candidate) -> bool:
    for _, c, e, points in DATA:
        for p, observed in points:
            predicted = candidate.evaluate(p, c, e)

            if predicted != observed:
                return False

    return True


def error_table(candidate: Candidate):
    result = {}

    for name, c, e, points in DATA:
        rows = []

        for p, observed in points:
            predicted = candidate.evaluate(p, c, e)

            rows.append(
                (
                    p,
                    observed,
                    predicted,
                    predicted == observed,
                )
            )

        result[name] = rows

    return result


# ---------------------------------------------------------------------------
# Search exact candidates
# ---------------------------------------------------------------------------

def search_candidates() -> list[Candidate]:
    found = []

    for mode in ("floor", "ceil"):
        for a in AFFINES:
            for d in AFFINES:
                candidate = Candidate(
                    mode=mode,
                    a=a,
                    d=d,
                )

                if is_exact(candidate):
                    found.append(candidate)

    found.sort(key=lambda x: x.complexity())

    return found


# ---------------------------------------------------------------------------
# More restricted families
# ---------------------------------------------------------------------------

def search_no_interaction() -> list[Candidate]:
    """
    Require c*e terms to vanish.

        a = a0 + a1*c + a2*e
        d = d0 + d1*c + d2*e
    """
    found = []

    reduced = []

    for a in AFFINES:
        if a.c3 == 0:
            reduced.append(a)

    for d in AFFINES:
        if d.c3 == 0:
            for mode in ("floor", "ceil"):
                candidate = Candidate(mode, a, d)

                if is_exact(candidate):
                    found.append(candidate)

    found.sort(key=lambda x: x.complexity())

    return found


def search_channel_base_plus_parity_shift() -> list[Candidate]:
    """
    Particularly simple form:

        a(c,e) = A0 + A1*c + A2*e
        d(c,e) = D0 + D1*c + D2*e

    with no interaction term.
    """
    return search_no_interaction()


# ---------------------------------------------------------------------------
# Direct laws suggested by the data
# ---------------------------------------------------------------------------

def direct_sector_formula(
    p: int,
    c: int,
    e: int,
) -> int:
    """
    Exact piecewise law inferred directly from the four sectors:

        A-even: max(p-2, floor(p/2))
        A-odd : max(p-2, ceil(p/2))
        B-even: max(p-2, floor(p/2))
        B-odd : max(0, p-2)
    """

    if c == 0 and e == 0:
        return max(
            0,
            p - 2,
            floor_div(p, 2),
        )

    if c == 0 and e == 1:
        return max(
            0,
            p - 2,
            ceil_div(p, 2),
        )

    if c == 1 and e == 0:
        return max(
            0,
            p - 2,
            floor_div(p, 2),
        )

    if c == 1 and e == 1:
        return max(
            0,
            p - 2,
        )

    raise ValueError("invalid c/e")


def direct_formula_exact() -> bool:
    for _, c, e, points in DATA:
        for p, observed in points:
            if direct_sector_formula(p, c, e) != observed:
                return False

    return True


# ---------------------------------------------------------------------------
# Search a slightly richer compressed family
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RichCandidate:
    mode: str
    a: int
    b: int
    delta: int

    def evaluate(
        self,
        p: int,
        c: int,
        e: int,
    ) -> int:
        """
        Search

            max(
                0,
                p - a - b*c,
                floor((p + delta + c*e)/2)
            )

        and its ceil version.
        """
        left = p - self.a - self.b * c
        numerator = p + self.delta + c * e

        if self.mode == "floor":
            half = floor_div(numerator, 2)
        else:
            half = ceil_div(numerator, 2)

        return max3(0, left, half)

    def text(self) -> str:
        fn = "floor" if self.mode == "floor" else "ceil"

        return (
            f"max(0, p-{self.a}-{self.b}c, "
            f"{fn}((p+{self.delta}+ce)/2))"
        )


def search_rich() -> list[RichCandidate]:
    found = []

    for mode in ("floor", "ceil"):
        for a in range(-4, 5):
            for b in range(-4, 5):
                for delta in range(-4, 5):
                    candidate = RichCandidate(
                        mode,
                        a,
                        b,
                        delta,
                    )

                    ok = True

                    for _, c, e, points in DATA:
                        for p, observed in points:
                            if (
                                candidate.evaluate(
                                    p, c, e
                                )
                                != observed
                            ):
                                ok = False
                                break

                        if not ok:
                            break

                    if ok:
                        found.append(candidate)

    found.sort(
        key=lambda x: (
            abs(x.a) + abs(x.b) + abs(x.delta),
            x.text(),
        )
    )

    return found


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

def heading(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:

    print("=" * 78)
    print("EXPERIMENT 110 — EXACT CHANNEL/PARITY SUPPORT-LAW MINIMIZATION")
    print("=" * 78)

    heading("1. OBSERVED SUPPORT")

    for name, c, e, points in DATA:
        print(
            f"  {name}: c={c} parity={e}"
        )

        for p, s in points:
            print(
                f"    p={p}: s={s}"
            )

    heading("2. DIRECT SECTOR LAW")

    print(
        "  A-even : max(0, p-2, floor(p/2))"
    )
    print(
        "  A-odd  : max(0, p-2, ceil(p/2))"
    )
    print(
        "  B-even : max(0, p-2, floor(p/2))"
    )
    print(
        "  B-odd  : max(0, p-2)"
    )

    print(
        f"  exact={direct_formula_exact()}"
    )

    heading("3. SMALL NO-INTERACTION AFFINE SEARCH")

    no_interaction = search_no_interaction()

    print(
        f"  exact candidates = {len(no_interaction)}"
    )

    for candidate in no_interaction[:20]:
        print(
            f"    {candidate.text()}"
            f"    score={candidate.complexity()}"
        )

    heading("4. GENERAL CHANNEL/PARITY AFFINE SEARCH")

    candidates = search_candidates()

    print(
        f"  exact candidates = {len(candidates)}"
    )

    for candidate in candidates[:30]:
        print(
            f"    {candidate.text()}"
            f"    score={candidate.complexity()}"
        )

    heading("5. RICH COMPRESSED SEARCH")

    rich = search_rich()

    print(
        f"  exact candidates = {len(rich)}"
    )

    for candidate in rich[:30]:
        print(
            f"    {candidate.text()}"
        )

    heading("6. BEST GENERAL CANDIDATE VALIDATION")

    if candidates:
        best = candidates[0]

        print(
            f"  selected = {best.text()}"
        )

        for name, c, e, points in DATA:
            print(f"  {name}")

            for p, observed in points:
                predicted = best.evaluate(
                    p, c, e
                )

                print(
                    f"    p={p}: "
                    f"observed={observed} "
                    f"predicted={predicted} "
                    f"exact={predicted == observed}"
                )
    else:
        print(
            "  No exact candidate in this family."
        )

    heading("7. RICH-CANDIDATE VALIDATION")

    if rich:
        best_rich = rich[0]

        print(
            f"  selected = {best_rich.text()}"
        )

        all_ok = True

        for name, c, e, points in DATA:
            print(f"  {name}")

            for p, observed in points:
                predicted = best_rich.evaluate(
                    p, c, e
                )

                ok = predicted == observed
                all_ok = all_ok and ok

                print(
                    f"    p={p}: "
                    f"observed={observed} "
                    f"predicted={predicted} "
                    f"exact={ok}"
                )

        print(
            f"  rich_candidate_exact={all_ok}"
        )
    else:
        print(
            "  No exact rich candidate."
        )

    heading("8. SUPPORT-LAW DIFFERENCE TABLE")

    for name, c, e, points in DATA:
        print(f"  {name}")

        for p, observed in points:
            predicted = direct_sector_formula(
                p, c, e
            )

            print(
                f"    p={p}: "
                f"s={observed} "
                f"direct={predicted} "
                f"difference={predicted - observed}"
            )

    heading("9. STRUCTURAL INTERPRETATION")

    print(
        """
  Experiment 109 showed that m=8/7 plus parity cannot encode the
  support law in the tested affine-delta family.

  The exact sector data instead gives:

      A-even:
          max(0, p-2, floor(p/2))

      A-odd:
          max(0, p-2, ceil(p/2))

      B-even:
          max(0, p-2, floor(p/2))

      B-odd:
          max(0, p-2).

  Experiment 110 asks whether this apparent four-sector behavior can
  still be compressed by introducing only the binary channel variable

      c = 0  (A)
      c = 1  (B)

  and binary parity variable

      e = 0  (even)
      e = 1  (odd).

  In particular, we search for a single exact formula in c and e,
  rather than four independently fitted formulas.

  A successful low-complexity candidate would show that the apparent
  channel split is merely a compact discrete interaction.

  A failure would establish that the B-odd sector is genuinely a
  separate boundary regime at the tested level.

  Everything is exact over the integers.
  No floating point.
  No extrapolation.
  No recurrence search.
        """.strip()
    )

    heading("10. FINAL EXACTNESS")

    direct_ok = direct_formula_exact()
    has_general = bool(candidates)
    has_rich = bool(rich)

    print(
        f"  direct_sector_law_exact = {direct_ok}"
    )
    print(
        f"  general_affine_exact_candidate = {has_general}"
    )
    print(
        f"  rich_exact_candidate = {has_rich}"
    )
    print(
        f"  general_candidate_count = {len(candidates)}"
    )
    print(
        f"  rich_candidate_count = {len(rich)}"
    )

    failures = 0 if direct_ok else 1

    print(
        f"  failures = {failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS = {failures == 0}"
    )

    print()
    print("EXPERIMENT 110 COMPLETE")


if __name__ == "__main__":
    main()

