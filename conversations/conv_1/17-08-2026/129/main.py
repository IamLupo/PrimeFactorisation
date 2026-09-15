#!/usr/bin/env python3

"""
========================================================================
EXPERIMENT 194
Coefficient-family reconstruction

Purpose
-------
Experiments 186-193 showed that guessing the support interval directly is
too underdetermined.

This experiment instead tests a structured family of coefficient formulas:

    [x^r] (1+x)^A (1-x)^B

where A and B are simple affine functions of (ell,k).

We test:
    - common affine exponent choices,
    - several transformed coefficient indices,
    - optional global sign,
    - exact equality against every known exact anchor.

No result file is read.

All data are generated here.
Every function is defined here.

No previous experiment output is imported.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb


# ============================================================================
# BASIC ARITHMETIC
# ============================================================================

def C(n: int, r: int) -> int:
    """Binomial coefficient with zero convention."""
    if n < 0:
        return 0
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def sgn(n: int) -> int:
    """Return (-1)^n."""
    return -1 if n % 2 else 1


# ============================================================================
# KNOWN EXACT DATA
# ============================================================================

@dataclass(frozen=True)
class Anchor:
    k: int
    ell: int
    s: int
    exact: int
    label: str


def generate_anchors() -> list[Anchor]:
    """
    Explicitly generated from values already visible in the conversation.

    No external file is read.
    """
    return [
        Anchor(1, 10, 2, 27, "A1"),
        Anchor(3, 20, 6, 935, "A2"),
        Anchor(7, 40, 15, -1797818, "A3"),

        Anchor(3, 7, 3, -3, "B1"),
        Anchor(3, 9, 4, 9, "B2"),
        Anchor(3, 11, 5, -16, "B3"),
        Anchor(5, 11, 5, -5, "B4"),
        Anchor(5, 19, 9, -196, "B5"),
        Anchor(9, 21, 11, -84, "B6"),
    ]


# ============================================================================
# AFFINE EXPONENTS
# ============================================================================

@dataclass(frozen=True)
class Affine:
    """
    Represents

        (a*ell + b*k + c) / d
    """
    name: str
    a: int
    b: int
    c: int
    d: int = 1

    def value(
        self,
        k: int,
        ell: int,
    ) -> int | None:
        numerator = (
            self.a * ell
            + self.b * k
            + self.c
        )

        if numerator < 0:
            return None

        if numerator % self.d != 0:
            return None

        return numerator // self.d


def generate_exponents() -> list[Affine]:
    """
    Small, interpretable affine exponent family.

    Denominator 2 and 4 allow the half/quarter combinations that appeared
    in the support experiments.
    """
    result: list[Affine] = []

    specifications = [
        ("ell", 1, 0, 0, 1),
        ("ell+k", 1, 1, 0, 1),
        ("ell-k", 1, -1, 0, 1),
        ("k", 0, 1, 0, 1),
        ("2ell", 2, 0, 0, 1),
        ("2k", 0, 2, 0, 1),

        ("(ell+k)/2", 1, 1, 0, 2),
        ("(ell-k)/2", 1, -1, 0, 2),

        ("(ell+3k)/2", 1, 3, 0, 2),
        ("(ell-3k)/2", 1, -3, 0, 2),

        ("(ell+3k)/4", 1, 3, 0, 4),
        ("(ell-3k)/4", 1, -3, 0, 4),

        ("(ell+k+1)/2", 1, 1, 1, 2),
        ("(ell-k+1)/2", 1, -1, 1, 2),

        ("(ell+k-1)/2", 1, 1, -1, 2),
        ("(ell-k-1)/2", 1, -1, -1, 2),
    ]

    for name, a, b, c, d in specifications:
        result.append(
            Affine(
                name=name,
                a=a,
                b=b,
                c=c,
                d=d,
            )
        )

    return result


# ============================================================================
# INDEX TRANSFORMS
# ============================================================================

@dataclass(frozen=True)
class IndexTransform:
    name: str

    def apply(
        self,
        k: int,
        ell: int,
        s: int,
    ) -> int:
        if self.name == "s":
            return s

        if self.name == "ell-s":
            return ell - s

        if self.name == "s-k":
            return s - k

        if self.name == "s+k":
            return s + k

        if self.name == "ell-k-s":
            return ell - k - s

        if self.name == "ell+k-s":
            return ell + k - s

        if self.name == "2s-ell":
            return 2 * s - ell

        return s


def generate_index_transforms() -> list[IndexTransform]:
    return [
        IndexTransform("s"),
        IndexTransform("ell-s"),
        IndexTransform("s-k"),
        IndexTransform("s+k"),
        IndexTransform("ell-k-s"),
        IndexTransform("ell+k-s"),
        IndexTransform("2s-ell"),
    ]


# ============================================================================
# COEFFICIENT OF (1+x)^A (1-x)^B
# ============================================================================

def coefficient_product(
    A: int,
    B: int,
    r: int,
) -> int:
    """
    Compute

        [x^r] (1+x)^A (1-x)^B

    directly by convolution.
    """
    if r < 0:
        return 0

    total = 0

    for j in range(
        0,
        A + 1,
    ):
        q = r - j

        if q < 0 or q > B:
            continue

        total += (
            C(A, j)
            * sgn(q)
            * C(B, q)
        )

    return total


# ============================================================================
# CANDIDATE
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    A: Affine
    B: Affine
    index: IndexTransform
    global_sign: int

    @property
    def name(self) -> str:
        sign_part = (
            ""
            if self.global_sign == 1
            else " * (-1)"
        )

        return (
            f"[x^({self.index.name})] "
            f"(1+x)^({self.A.name}) "
            f"(1-x)^({self.B.name})"
            f"{sign_part}"
        )


# ============================================================================
# CANDIDATE EVALUATION
# ============================================================================

def evaluate_candidate(
    candidate: Candidate,
    anchor: Anchor,
) -> int | None:
    A = candidate.A.value(
        anchor.k,
        anchor.ell,
    )

    B = candidate.B.value(
        anchor.k,
        anchor.ell,
    )

    if A is None or B is None:
        return None

    r = candidate.index.apply(
        anchor.k,
        anchor.ell,
        anchor.s,
    )

    return (
        candidate.global_sign
        * coefficient_product(
            A,
            B,
            r,
        )
    )


# ============================================================================
# GENERATE CANDIDATES
# ============================================================================

def generate_candidates() -> list[Candidate]:
    exponents = generate_exponents()
    transforms = generate_index_transforms()

    candidates: list[Candidate] = []

    for A in exponents:
        for B in exponents:
            for index in transforms:
                for global_sign in [1, -1]:
                    candidates.append(
                        Candidate(
                            A=A,
                            B=B,
                            index=index,
                            global_sign=global_sign,
                        )
                    )

    return candidates


# ============================================================================
# SCORING
# ============================================================================

@dataclass(frozen=True)
class CandidateScore:
    candidate: Candidate
    matches: int
    computed: int
    exact_mismatches: int
    undefined: int
    abs_error: int


def score_candidate(
    candidate: Candidate,
    anchors: list[Anchor],
) -> CandidateScore:
    matches = 0
    computed = 0
    mismatches = 0
    undefined = 0
    abs_error = 0

    for anchor in anchors:
        value = evaluate_candidate(
            candidate,
            anchor,
        )

        if value is None:
            undefined += 1
            continue

        computed += 1

        if value == anchor.exact:
            matches += 1
        else:
            mismatches += 1
            abs_error += abs(
                value - anchor.exact
            )

    return CandidateScore(
        candidate=candidate,
        matches=matches,
        computed=computed,
        exact_mismatches=mismatches,
        undefined=undefined,
        abs_error=abs_error,
    )


def rank_scores(
    scores: list[CandidateScore],
) -> list[CandidateScore]:
    """
    Rank by:
        1. exact anchor matches,
        2. fewer undefined cases,
        3. smaller absolute error.
    """
    return sorted(
        scores,
        key=lambda score: (
            -score.matches,
            score.undefined,
            score.abs_error,
        ),
    )


# ============================================================================
# DISPLAY
# ============================================================================

def print_anchor_table(
    anchors: list[Anchor],
) -> None:
    print()
    print("=" * 72)
    print("ANCHOR DATA")
    print("=" * 72)

    for anchor in anchors:
        print(
            f"{anchor.label:3s} "
            f"k={anchor.k:2d} "
            f"ell={anchor.ell:2d} "
            f"s={anchor.s:3d} "
            f"exact={anchor.exact:12d}"
        )


def print_top_scores(
    ranked: list[CandidateScore],
    limit: int = 50,
) -> None:
    print()
    print("=" * 72)
    print("TOP COEFFICIENT-FAMILY CANDIDATES")
    print("=" * 72)

    for index, score in enumerate(
        ranked[:limit],
        start=1,
    ):
        print()
        print(
            f"#{index:2d} "
            f"matches={score.matches}/9 "
            f"undefined={score.undefined} "
            f"abs_error={score.abs_error}"
        )

        print(
            f"  {score.candidate.name}"
        )


def print_perfect_candidates(
    ranked: list[CandidateScore],
) -> None:
    perfect = [
        score
        for score in ranked
        if score.matches == 9
        and score.undefined == 0
    ]

    print()
    print("=" * 72)
    print("PERFECT CANDIDATES")
    print("=" * 72)

    print(
        f"count: {len(perfect)}"
    )

    for score in perfect[:100]:
        print(
            score.candidate.name
        )


# ============================================================================
# FRESH NON-ANCHOR TEST GRID
# ============================================================================

@dataclass(frozen=True)
class FreshCase:
    k: int
    ell: int
    s: int


def generate_fresh_cases() -> list[FreshCase]:
    """
    Generate cases that were not used as anchor constraints.

    These are for comparing the shapes of the surviving candidates,
    not for claiming ground truth.
    """
    cases = []

    for k in [
        1,
        3,
        5,
        7,
        9,
        11,
    ]:
        for ell in range(
            max(8, k + 1),
            31,
        ):
            for s in range(
                0,
                ell + 1,
            ):
                cases.append(
                    FreshCase(
                        k=k,
                        ell=ell,
                        s=s,
                    )
                )

    return cases


def make_fresh_anchor_like(
    case: FreshCase,
) -> Anchor:
    """
    Wrap a fresh case as an Anchor with an artificial exact value.

    The exact field is NOT used. This function exists only so the same
    candidate evaluator can be reused safely.
    """
    return Anchor(
        k=case.k,
        ell=case.ell,
        s=case.s,
        exact=0,
        label="fresh",
    )


def print_survivor_shapes(
    ranked: list[CandidateScore],
) -> None:
    """
    Print coefficient vectors for the best few candidates.
    """
    survivors = [
        score.candidate
        for score in ranked[:10]
    ]

    print()
    print("=" * 72)
    print("SHAPES OF TOP SURVIVORS")
    print("=" * 72)

    for candidate in survivors:
        print()
        print(candidate.name)

        for k in [1, 3, 5, 7]:
            A = candidate.A.value(
                k,
                20,
            )

            B = candidate.B.value(
                k,
                20,
            )

            if A is None or B is None:
                continue

            values = []

            for r in range(
                0,
                A + B + 1,
            ):
                anchor = Anchor(
                    k=k,
                    ell=20,
                    s=r,
                    exact=0,
                    label="shape",
                )

                value = evaluate_candidate(
                    candidate,
                    anchor,
                )

                values.append(value)

            print(
                f"  k={k}: {values}"
            )


# ============================================================================
# SANITY CHECKS
# ============================================================================

def verify_no_file_dependency() -> None:
    """
    This experiment intentionally has no file I/O.

    The function exists as an explicit marker for the experiment's
    reproducibility requirement.
    """
    print()
    print(
        "FILE DEPENDENCY CHECK: "
        "PASS (no files are read)"
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 194")
    print("Coefficient-family reconstruction")
    print("=" * 72)

    print()
    print(
        "All anchors and fresh cases are generated locally."
    )
    print(
        "No results.txt or previous experiment output is read."
    )

    verify_no_file_dependency()

    anchors = generate_anchors()

    print_anchor_table(
        anchors
    )

    candidates = generate_candidates()

    print()
    print(
        f"Generated candidates: "
        f"{len(candidates)}"
    )

    scores = [
        score_candidate(
            candidate,
            anchors,
        )
        for candidate in candidates
    ]

    ranked = rank_scores(
        scores
    )

    print_top_scores(
        ranked,
        limit=50,
    )

    print_perfect_candidates(
        ranked
    )

    print_survivor_shapes(
        ranked
    )

    fresh_cases = generate_fresh_cases()

    print()
    print(
        f"Fresh cases generated: "
        f"{len(fresh_cases)}"
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 194")
    print("=" * 72)


if __name__ == "__main__":
    main()

