"""
==============================================================================
EXPERIMENT 109 — EXACT MINIMAL UNIFIED SUPPORT-LAW SEARCH
==============================================================================

Goal
----
Starting from the exact terminal multiplicities observed in Experiment 108,
search for the simplest unified discrete law

    s(p, m, parity)

using only small affine expressions together with floor/ceil/max.

This is a finite structural audit only.

No floating point.
No recurrence search.
No extrapolation.
==============================================================================

Data
----
A channel:
    m = 8
    K = 6

B channel:
    m = 7
    K = 5

parity:
    0 = even
    1 = odd

Observed terminal multiplicities:

A-even:
    p=0 -> 0
    p=2 -> 1
    p=4 -> 2
    p=6 -> 4
    p=8 -> 6

A-odd:
    p=1 -> 1
    p=3 -> 2
    p=5 -> 3
    p=7 -> 5

B-even:
    p=0 -> 0
    p=2 -> 1
    p=4 -> 2
    p=6 -> 4

B-odd:
    p=1 -> 0
    p=3 -> 1
    p=5 -> 3
    p=7 -> 5
==============================================================================

Candidate families
------------------

We test expressions of the form

    max(p - a, floor((p + delta) / 2))
    max(p - a, ceil ((p + delta) / 2))

where

    a, delta

are small affine functions of

    m, parity

and, separately, small direct candidates involving m.

The objective is NOT to select an arbitrary exact formula when many exist.
Instead we score each exact law by:

    1. number of free parameters;
    2. expression complexity;
    3. whether it uses parity explicitly;
    4. whether the same coefficients work for both channels.

==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Callable, Iterable


# ---------------------------------------------------------------------------
# Exact observed data
# ---------------------------------------------------------------------------

DATA = {
    "A-even": {
        "m": 8,
        "K": 6,
        "parity": 0,
        "points": [(0, 0), (2, 1), (4, 2), (6, 4), (8, 6)],
    },
    "A-odd": {
        "m": 8,
        "K": 6,
        "parity": 1,
        "points": [(1, 1), (3, 2), (5, 3), (7, 5)],
    },
    "B-even": {
        "m": 7,
        "K": 5,
        "parity": 0,
        "points": [(0, 0), (2, 1), (4, 2), (6, 4)],
    },
    "B-odd": {
        "m": 7,
        "K": 5,
        "parity": 1,
        "points": [(1, 0), (3, 1), (5, 3), (7, 5)],
    },
}


# ---------------------------------------------------------------------------
# Exact floor / ceil division for integers
# ---------------------------------------------------------------------------

def floor_div(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError("b must be positive")
    return a // b


def ceil_div(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError("b must be positive")
    return -((-a) // b)


def max2(a: int, b: int) -> int:
    return a if a >= b else b


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def exact_on_all_points(
    law: Callable[[int, int, int], int],
) -> bool:
    for sector in DATA.values():
        m = sector["m"]
        parity = sector["parity"]
        for p, observed in sector["points"]:
            predicted = law(p, m, parity)
            if predicted != observed:
                return False
    return True


def sector_error_table(
    law: Callable[[int, int, int], int],
) -> dict[str, list[tuple[int, int, int]]]:
    out: dict[str, list[tuple[int, int, int]]] = {}

    for name, sector in DATA.items():
        m = sector["m"]
        parity = sector["parity"]
        rows = []

        for p, observed in sector["points"]:
            predicted = law(p, m, parity)
            rows.append((p, observed, predicted))

        out[name] = rows

    return out


def complexity_score(text: str) -> tuple[int, int]:
    """
    Small deterministic textual complexity score.

    First component:
        number of operators / structural tokens.

    Second component:
        string length.

    This is only a tie-breaker among exact formulas.
    """
    token_count = sum(text.count(tok) for tok in (
        "max", "floor", "ceil", "+", "-", "*"
    ))
    return token_count, len(text)


# ---------------------------------------------------------------------------
# Affine parameterization
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Affine:
    """
    a(m, parity) = c0 + cm*m + cp*parity
    """
    c0: int
    cm: int
    cp: int

    def __call__(self, m: int, parity: int) -> int:
        return self.c0 + self.cm * m + self.cp * parity

    def text(self, variable_m: str = "m", variable_p: str = "e") -> str:
        parts: list[str] = []

        if self.c0:
            parts.append(str(self.c0))

        if self.cm:
            term = f"{abs(self.cm)}{variable_m}"
            if self.cm < 0:
                term = "-" + term
            parts.append(term)

        if self.cp:
            term = f"{abs(self.cp)}{variable_p}"
            if self.cp < 0:
                term = "-" + term
            parts.append(term)

        if not parts:
            return "0"

        result = parts[0]

        for part in parts[1:]:
            if part.startswith("-"):
                result += " - " + part[1:]
            else:
                result += " + " + part

        return result


# ---------------------------------------------------------------------------
# Candidate construction
# ---------------------------------------------------------------------------

AFFINE_RANGE = range(-4, 5)


def generate_affines() -> list[Affine]:
    out = []

    for c0, cm, cp in product(
        AFFINE_RANGE,
        AFFINE_RANGE,
        AFFINE_RANGE,
    ):
        # Keep the search compact by rejecting unnecessarily large forms.
        complexity = abs(c0) + abs(cm) + abs(cp)

        if complexity <= 5:
            out.append(Affine(c0, cm, cp))

    return out


AFFINES = generate_affines()


@dataclass(frozen=True)
class Candidate:
    kind: str
    a: Affine
    delta: Affine

    def evaluate(self, p: int, m: int, parity: int) -> int:
        a_val = self.a(m, parity)
        d_val = self.delta(m, parity)

        left = p - a_val

        if self.kind == "floor":
            right = floor_div(p + d_val, 2)
        elif self.kind == "ceil":
            right = ceil_div(p + d_val, 2)
        else:
            raise ValueError(f"unknown kind: {self.kind}")

        return max2(left, right)

    def text(self) -> str:
        a_text = self.a.text()
        d_text = self.delta.text()

        right = f"floor((p + ({d_text}))/2)"
        if self.kind == "ceil":
            right = f"ceil((p + ({d_text}))/2)"

        return f"max(p - ({a_text}), {right})"

    def score(self) -> tuple[int, int, int, int]:
        """
        Lower is better.

        Prefer:
            - fewer explicit parity dependencies;
            - smaller affine coefficient complexity;
            - shorter resulting formula.
        """
        parity_use = (
            abs(self.a.cp) +
            abs(self.delta.cp)
        )

        affine_size = (
            abs(self.a.c0) +
            abs(self.a.cm) +
            abs(self.a.cp) +
            abs(self.delta.c0) +
            abs(self.delta.cm) +
            abs(self.delta.cp)
        )

        text_len = len(self.text())

        return (
            parity_use,
            affine_size,
            text_len,
            0 if self.kind == "floor" else 1,
        )


def candidate_is_exact(c: Candidate) -> bool:
    return exact_on_all_points(c.evaluate)


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def search_candidates() -> list[Candidate]:
    exact: list[Candidate] = []

    for kind in ("floor", "ceil"):
        for a in AFFINES:
            for delta in AFFINES:
                c = Candidate(kind, a, delta)

                if candidate_is_exact(c):
                    exact.append(c)

    exact.sort(key=lambda c: c.score())
    return exact


# ---------------------------------------------------------------------------
# More restricted structural families
# ---------------------------------------------------------------------------

def search_delta_only() -> list[tuple[str, int]]:
    """
    Search

        max(p-2, floor((p+d)/2))
        max(p-2, ceil ((p+d)/2))

    with constant integer d.
    """
    out = []

    for kind in ("floor", "ceil"):
        for d in range(-8, 9):
            if exact_on_all_points(
                lambda p, m, parity, kind=kind, d=d:
                max2(
                    p - 2,
                    floor_div(p + d, 2)
                    if kind == "floor"
                    else ceil_div(p + d, 2),
                )
            ):
                out.append((kind, d))

    return out


def search_m_only() -> list[tuple[str, int, int, int]]:
    """
    Search

        max(p - (m + a),
            floor((p + b*m + c)/2))

    and the ceil analogue.
    """
    out = []

    for kind in ("floor", "ceil"):
        for a in range(-8, 9):
            for b in range(-3, 4):
                for c in range(-8, 9):

                    def law(
                        p: int,
                        m: int,
                        parity: int,
                        kind=kind,
                        a=a,
                        b=b,
                        c=c,
                    ) -> int:
                        left = p - (m + a)
                        numerator = p + b * m + c

                        if kind == "floor":
                            right = floor_div(numerator, 2)
                        else:
                            right = ceil_div(numerator, 2)

                        return max2(left, right)

                    if exact_on_all_points(law):
                        out.append((kind, a, b, c))

    return out


# ---------------------------------------------------------------------------
# Direct validation of the preferred four-sector assignment
# ---------------------------------------------------------------------------

def preferred_law(p: int, m: int, parity: int) -> int:
    """
    Exact law corresponding to the selected Experiment-108 assignment:

        A-even: delta = 0
        A-odd : delta = 1
        B-even: delta = 0
        B-odd : delta = -1

    The sector can be encoded through the channel offset m:

        A has m=8
        B has m=7

    and the parity contribution.

    We test a few compact equivalent parameterizations separately.
    """
    # Recover the channel:
    if m == 8:
        delta = 0 + parity
    elif m == 7:
        delta = -1 + parity
    else:
        raise ValueError("Unexpected m")

    return max2(p - 2, floor_div(p + delta, 2))


def validate_preferred() -> bool:
    return exact_on_all_points(preferred_law)


# ---------------------------------------------------------------------------
# Attempt a simple m/parity compression of delta
# ---------------------------------------------------------------------------

def search_delta_affine_in_m_parity() -> list[Affine]:
    """
    Search delta(m,e) = c0 + cm*m + cp*e.
    """
    exact = []

    for aff in AFFINES:
        def law(
            p: int,
            m: int,
            parity: int,
            aff=aff,
        ) -> int:
            delta = aff(m, parity)
            return max2(
                p - 2,
                floor_div(p + delta, 2),
            )

        if exact_on_all_points(law):
            exact.append(aff)

    return exact


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def print_heading(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def main() -> None:
    print("=" * 78)
    print("EXPERIMENT 109 — EXACT MINIMAL UNIFIED SUPPORT-LAW SEARCH")
    print("=" * 78)

    print_heading("1. EXACT OBSERVED SUPPORT DATA")

    for name, sector in DATA.items():
        print(
            f"  {name}: m={sector['m']} "
            f"K={sector['K']} parity={sector['parity']}"
        )

        for p, s in sector["points"]:
            print(f"    p={p}: s={s}")

    print_heading("2. DIRECT DELTA-ONLY SEARCH")

    delta_only = search_delta_only()

    if delta_only:
        for kind, d in delta_only:
            print(f"  {kind}, delta={d}: EXACT")
    else:
        print("  NONE")

    print_heading("3. AFFINE DELTA(m,parity) SEARCH")

    affine_delta = search_delta_affine_in_m_parity()

    if affine_delta:
        for aff in affine_delta[:50]:
            print(
                f"  delta(m,e) = {aff.text()} : EXACT"
            )

        if len(affine_delta) > 50:
            print(
                f"  ... {len(affine_delta) - 50} more exact affine laws"
            )
    else:
        print("  NONE")

    print_heading("4. MINIMAL AFFINE DELTA SEARCH WITH p-2 BASE")

    if affine_delta:
        best = sorted(
            affine_delta,
            key=lambda a: (
                abs(a.c0) + abs(a.cm) + abs(a.cp),
                abs(a.cm),
                abs(a.cp),
                abs(a.c0),
            ),
        )[0]

        print(f"  best delta(m,e) = {best.text()}")

        def best_law(p: int, m: int, parity: int) -> int:
            d = best(m, parity)
            return max2(p - 2, floor_div(p + d, 2))

        print(
            f"  globally exact = {exact_on_all_points(best_law)}"
        )

    print_heading("5. GENERAL SMALL AFFINE a(m,e), delta(m,e) SEARCH")

    exact_candidates = search_candidates()

    print(f"  exact candidate count = {len(exact_candidates)}")

    if exact_candidates:
        print("  best exact candidates:")

        for c in exact_candidates[:20]:
            print(
                f"    {c.text()}"
                f"    score={c.score()}"
            )

    print_heading("6. BEST CANDIDATE VALIDATION")

    if exact_candidates:
        best = exact_candidates[0]

        print(f"  selected = {best.text()}")

        all_exact = True

        for name, sector in DATA.items():
            m = sector["m"]
            parity = sector["parity"]

            print(f"  {name}")

            for p, observed in sector["points"]:
                predicted = best.evaluate(p, m, parity)
                ok = (predicted == observed)
                all_exact = all_exact and ok

                print(
                    f"    p={p}: "
                    f"observed={observed} "
                    f"predicted={predicted} "
                    f"exact={ok}"
                )

        print(f"  best_candidate_exact={all_exact}")
    else:
        all_exact = False
        print("  no exact candidate found")

    print_heading("7. PREFERRED EXPERIMENT-108 LAW")

    preferred_ok = validate_preferred()

    print(
        "  s(p)=max(p-2, floor((p+delta(m,e))/2))"
    )
    print(
        "  delta = 0+e for m=8"
    )
    print(
        "  delta = -1+e for m=7"
    )
    print(f"  exact={preferred_ok}")

    for name, sector in DATA.items():
        print(f"  {name}")

        m = sector["m"]
        parity = sector["parity"]

        if m == 8:
            delta = parity
        else:
            delta = -1 + parity

        print(f"    delta={delta}")

        for p, observed in sector["points"]:
            predicted = preferred_law(p, m, parity)

            print(
                f"    p={p}: "
                f"observed={observed} "
                f"predicted={predicted} "
                f"exact={predicted == observed}"
            )

    print_heading("8. m-ONLY STRUCTURAL SEARCH")

    m_only = search_m_only()

    if m_only:
        for item in m_only[:30]:
            print(
                "  kind=%s  a=%s  b=%s  c=%s"
                % item
            )

        if len(m_only) > 30:
            print(
                f"  ... {len(m_only) - 30} more exact m-only laws"
            )
    else:
        print("  NONE")

    print_heading("9. SECTOR ERROR SANITY")

    for name, rows in sector_error_table(
        exact_candidates[0].evaluate
        if exact_candidates
        else preferred_law
    ).items():

        bad = [
            (p, obs, pred)
            for p, obs, pred in rows
            if obs != pred
        ]

        print(
            f"  {name}: "
            f"bad_points={bad}"
        )

    print_heading("10. STRUCTURAL INTERPRETATION")

    print(
        """
  The observed terminal multiplicity is an integer-valued finite support
  quantity.  Rather than interpolating it as an arbitrary polynomial,
  this experiment searches a deliberately small family of discrete laws.

  The principal family is

      s(p,m,e)
        = max(
            p - a(m,e),
            floor((p + delta(m,e))/2)
          )

  together with the corresponding ceil form.

  The important question is whether the same affine functions of

      m = endpoint index
      e = parity

  explain BOTH channels.

  An exact common affine law is stronger than four separately fitted
  channel formulas because the channel distinction is compressed into
  the endpoint parameter and parity.

  This remains a finite identity audit on the observed support data.

  No statement is made about unseen values of m or K.
        """.strip()
    )

    print_heading("11. FINAL EXACTNESS")

    preferred = validate_preferred()
    candidate_exists = bool(exact_candidates)

    print(f"  preferred_law_exact = {preferred}")
    print(f"  exact_affine_candidate_exists = {candidate_exists}")
    print(f"  total_exact_candidates = {len(exact_candidates)}")

    failures = 0 if preferred and candidate_exists else 1

    print(f"  failures = {failures}")

    if failures == 0:
        print("  ALL BASIC CHECKS PASS = True")
    else:
        print("  ALL BASIC CHECKS PASS = False")

    print()
    print("EXPERIMENT 109 COMPLETE")


if __name__ == "__main__":
    main()

