#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 364R — EXACT CROSS-TRANSFORM Z-CONSISTENCY / COMMON-ANNIHILATOR
                      COMPATIBILITY AUDIT
==============================================================================

Purpose
-------
Experiment 363R discovered seven exact primitive integer annihilators after
canonical finite-difference transformations. Each one forces a rational
prediction for the SAME missing source value

    Z = Q_3(5) = Q(2,3).

Those forced values were:

    H10 / transport_3
    H10 / lower_transport_3
    H01 / horizontal_width3
    H01 / diamond_5
    H11 / vertical_width3
    H02 / transport_3
    H02 / lower_transport_3

364R asks the next logically stronger question:

    Are these discoveries mutually compatible?

Audits
------

1. Exact reconstruction of all discovered Z-values.

2. Exact equality matrix.

3. Pairwise differences.

4. Rational compatibility classes.

5. Integer compatibility.

6. Small-denominator compatibility.

7. Cross-transform intersection:
       can two or more independent transformed laws hold simultaneously?

8. Common-law intersection:
       does every discovered law force the same Z?

9. Exact primitive relation cross-check.

10. Optional symbolic residual verification for each relation.

IMPORTANT
---------
No predicted Z is inserted into the source table.

The seven values remain model consequences only.

No missing source value is treated as observed.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.

Output is deliberately compact.
"""


from __future__ import annotations

import math
import sys
from collections import defaultdict

import sympy as sp


# ============================================================================
# SYMBOL
# ============================================================================

Z = sp.Symbol("Z")


# ============================================================================
# DISCOVERIES FROM 363R
# ============================================================================

DISCOVERIES = [
    {
        "name": "H10 / transport_3",
        "relation": [
            439971449708397458,
            172645072515828862,
            -1320984232593184369,
        ],
        "forced_Z":
            sp.Rational(
                -3191330642348143851622391184,
                1320984232593184369,
            ),
    },

    {
        "name": "H10 / lower_transport_3",
        "relation": [
            39782170988615124,
            -86322536257914431,
            -240269201882469729,
        ],
        "forced_Z":
            sp.Rational(
                -392992110572009658827397946,
                240269201882469729,
            ),
    },

    {
        "name": "H01 / horizontal_width3",
        "relation": [
            2392847188642012,
            -43785379160637941,
            1750377005640251053,
        ],
        "forced_Z":
            sp.Rational(
                -50926698968037761526399730,
                1750377005640251053,
            ),
    },

    {
        "name": "H01 / diamond_5",
        "relation": [
            519086597060649301029083460844675708,
            -572993045081469360489454673312828768,
            328024347471513762409086608714698424,
            -281452566868358673109594874417536446,
            5342832968590334980318521932388509869,
        ],
        "forced_Z":
            sp.Rational(
                732126390015674197979598697500986324215505626,
                5342832968590334980318521932388509869,
            ),
    },

    {
        "name": "H11 / vertical_width3",
        "relation": [
            52652305849810441236,
            13646355081654935184,
            12359557018829527669,
        ],
        "forced_Z":
            sp.Rational(
                -42367633024449172779805183330,
                12359557018829527669,
            ),
    },

    {
        "name": "H02 / transport_3",
        "relation": [
            21743049162319153956,
            9380422859347821304,
            79729788761536106653,
        ],
        "forced_Z":
            sp.Rational(
                -1636649764753763519072419014,
                79729788761536106653,
            ),
    },

    {
        "name": "H02 / lower_transport_3",
        "relation": [
            72877453359532340,
            -9380422859347821304,
            -4210540774600727807,
        ],
        "forced_Z":
            sp.Rational(
                -42490000284788818260262826,
                4210540774600727807,
            ),
    },
]


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def is_integer(value):
    value = sp.Rational(value)
    return value.q == 1


def denominator(value):
    return int(
        sp.denom(
            sp.Rational(value)
        )
    )


def numerator(value):
    return int(
        sp.numer(
            sp.Rational(value)
        )
    )


def factor_integer(value):
    value = abs(int(value))

    if value == 0:
        return {}

    return sp.factorint(value)


def primitive_integer_relation(vector):
    values = [
        sp.Rational(v)
        for v in vector
    ]

    lcm = 1

    for value in values:
        lcm = sp.ilcm(
            lcm,
            int(value.q),
        )

    integers = [
        int(value * lcm)
        for value in values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        integers = [
            value // g
            for value in integers
        ]

    for value in integers:
        if value != 0 and value < 0:
            integers = [
                -v
                for v in integers
            ]
            break

    return integers


def relation_symbolic_residual(
    relation,
    symbolic_position,
):
    """
    Given a primitive relation c_0,...,c_{m-1},
    treat the designated entry as Z.
    """

    residual = 0

    for index, coefficient in enumerate(
        relation
    ):

        if index == symbolic_position:
            residual += (
                sp.Integer(coefficient)
                * Z
            )
        else:
            residual += sp.Integer(
                coefficient
            )

    return clean(
        residual
    )


# ============================================================================
# 1. DISCOVERY NORMALIZATION
# ============================================================================

def discovery_audit():

    print()
    print("=" * 78)
    print(
        "1. EXACT DISCOVERY NORMALIZATION"
    )
    print("=" * 78)

    for item in DISCOVERIES:

        value = clean(
            item["forced_Z"]
        )

        relation = primitive_integer_relation(
            item["relation"]
        )

        print()
        print(
            "  {}".format(
                item["name"]
            )
        )

        print(
            "    forced_Z={}".format(
                value
            )
        )

        print(
            "    numerator={}".format(
                numerator(value)
            )
        )

        print(
            "    denominator={}".format(
                denominator(value)
            )
        )

        print(
            "    integer={}".format(
                is_integer(value)
            )
        )

        print(
            "    relation_primitive={}".format(
                relation
            )
        )


# ============================================================================
# 2. EXACT EQUALITY MATRIX
# ============================================================================

def equality_matrix():

    print()
    print("=" * 78)
    print(
        "2. EXACT FORCED-Z EQUALITY MATRIX"
    )
    print("=" * 78)

    n = len(DISCOVERIES)

    matrix = []

    for i in range(n):

        row = []

        for j in range(n):

            zi = sp.Rational(
                DISCOVERIES[i]["forced_Z"]
            )

            zj = sp.Rational(
                DISCOVERIES[j]["forced_Z"]
            )

            row.append(
                zi == zj
            )

        matrix.append(row)

    print(
        "  labels={}".format(
            [
                item["name"]
                for item in DISCOVERIES
            ]
        )
    )

    for i, row in enumerate(matrix):

        print(
            "  {}: {}".format(
                i,
                row,
            )
        )

    return matrix


# ============================================================================
# 3. PAIRWISE DIFFERENCES
# ============================================================================

def pairwise_difference_audit():

    print()
    print("=" * 78)
    print(
        "3. PAIRWISE FORCED-Z DIFFERENCE AUDIT"
    )
    print("=" * 78)

    nonzero_count = 0
    equal_pairs = []

    for i in range(
        len(DISCOVERIES)
    ):

        for j in range(
            i + 1,
            len(DISCOVERIES),
        ):

            zi = sp.Rational(
                DISCOVERIES[i]["forced_Z"]
            )

            zj = sp.Rational(
                DISCOVERIES[j]["forced_Z"]
            )

            difference = clean(
                zi - zj
            )

            if difference == 0:

                equal_pairs.append(
                    (
                        DISCOVERIES[i]["name"],
                        DISCOVERIES[j]["name"],
                    )
                )

            else:

                nonzero_count += 1

    print(
        "  equal_pairs={}".format(
            equal_pairs
        )
    )

    print(
        "  nonzero_pairwise_differences={}".format(
            nonzero_count
        )
    )

    return equal_pairs, nonzero_count


# ============================================================================
# 4. COMPATIBILITY CLASSES
# ============================================================================

def compatibility_classes():

    print()
    print("=" * 78)
    print(
        "4. EXACT RATIONAL COMPATIBILITY CLASSES"
    )
    print("=" * 78)

    classes = defaultdict(list)

    for item in DISCOVERIES:

        value = clean(
            item["forced_Z"]
        )

        classes[value].append(
            item["name"]
        )

    for index, (
        value,
        names,
    ) in enumerate(
        classes.items()
    ):

        print()
        print(
            "  class_{}:".format(
                index
            )
        )

        print(
            "    Z={}".format(
                value
            )
        )

        print(
            "    members={}".format(
                names
            )
        )

    print()
    print(
        "  distinct_Z_values={}".format(
            len(classes)
        )
    )

    return classes


# ============================================================================
# 5. INTEGER FEASIBILITY
# ============================================================================

def integer_feasibility():

    print()
    print("=" * 78)
    print(
        "5. INTEGER-COMPATIBILITY AUDIT"
    )
    print("=" * 78)

    integer_discoveries = []

    for item in DISCOVERIES:

        value = sp.Rational(
            item["forced_Z"]
        )

        if is_integer(value):
            integer_discoveries.append(
                item["name"]
            )

    print(
        "  integer_compatible_discoveries={}".format(
            integer_discoveries
        )
    )

    print(
        "  integer_compatible_count={}".format(
            len(integer_discoveries)
        )
    )

    all_integer = (
        len(integer_discoveries)
        == len(DISCOVERIES)
    )

    print(
        "  all_discoveries_integer_compatible={}".format(
            all_integer
        )
    )

    return integer_discoveries


# ============================================================================
# 6. DENOMINATOR PROFILE
# ============================================================================

def denominator_profile():

    print()
    print("=" * 78)
    print(
        "6. FORCED-Z DENOMINATOR PROFILE"
    )
    print("=" * 78)

    denominators = []

    for item in DISCOVERIES:

        value = sp.Rational(
            item["forced_Z"]
        )

        d = denominator(value)

        denominators.append(d)

        print()
        print(
            "  {}".format(
                item["name"]
            )
        )

        print(
            "    denominator={}".format(
                d
            )
        )

        print(
            "    factorization={}".format(
                factor_integer(d)
            )
        )

    return denominators


# ============================================================================
# 7. SHARED DENOMINATOR / COMMON VALUE AUDIT
# ============================================================================

def common_rational_audit():

    print()
    print("=" * 78)
    print(
        "7. COMMON-RATIONAL / COMMON-DENOMINATOR AUDIT"
    )
    print("=" * 78)

    values = [
        sp.Rational(
            item["forced_Z"]
        )
        for item in DISCOVERIES
    ]

    common_denominator = 1

    for value in values:

        common_denominator = sp.ilcm(
            common_denominator,
            int(value.q),
        )

    scaled_numerators = [
        int(
            value
            * common_denominator
        )
        for value in values
    ]

    numerator_gcd = 0

    for value in scaled_numerators:
        numerator_gcd = math.gcd(
            numerator_gcd,
            abs(value),
        )

    print(
        "  least_common_denominator={}".format(
            common_denominator
        )
    )

    print(
        "  scaled_numerators={}".format(
            scaled_numerators
        )
    )

    print(
        "  gcd_of_scaled_numerators={}".format(
            numerator_gcd
        )
    )

    # A common rational value exists iff all values are equal.
    common_value = (
        values[0]
        if all(
            value == values[0]
            for value in values[1:]
        )
        else None
    )

    print(
        "  common_forced_value={}".format(
            common_value
        )
    )

    return (
        common_denominator,
        common_value,
    )


# ============================================================================
# 8. TRANSFORM-LEVEL CONSISTENCY
# ============================================================================

def transform_consistency():

    print()
    print("=" * 78)
    print(
        "8. TRANSFORM-LEVEL CONSISTENCY AUDIT"
    )
    print("=" * 78)

    groups = defaultdict(list)

    for item in DISCOVERIES:

        transform = item["name"].split(
            " / "
        )[0]

        groups[transform].append(
            item
        )

    for transform, items in groups.items():

        values = [
            sp.Rational(
                item["forced_Z"]
            )
            for item in items
        ]

        identical = all(
            value == values[0]
            for value in values[1:]
        )

        print()
        print(
            "  transform={}:".format(
                transform
            )
        )

        print(
            "    discovery_count={}".format(
                len(items)
            )
        )

        print(
            "    forced_Z_values={}".format(
                values
            )
        )

        print(
            "    internally_consistent={}".format(
                identical
            )
        )

    return groups


# ============================================================================
# 9. PAIRWISE COMMON-LAW AUDIT
# ============================================================================

def pairwise_common_law_audit():

    print()
    print("=" * 78)
    print(
        "9. PAIRWISE COMMON-LAW INTERSECTION AUDIT"
    )
    print("=" * 78)

    compatible_pairs = []

    for i in range(
        len(DISCOVERIES)
    ):

        for j in range(
            i + 1,
            len(DISCOVERIES),
        ):

            zi = sp.Rational(
                DISCOVERIES[i]["forced_Z"]
            )

            zj = sp.Rational(
                DISCOVERIES[j]["forced_Z"]
            )

            if zi == zj:

                compatible_pairs.append(
                    (
                        DISCOVERIES[i]["name"],
                        DISCOVERIES[j]["name"],
                        zi,
                    )
                )

    if compatible_pairs:

        for pair in compatible_pairs:

            print(
                "  compatible_pair={}".format(
                    pair
                )
            )

    else:

        print(
            "  compatible_pairs=[]"
        )

    print(
        "  compatible_pair_count={}".format(
            len(compatible_pairs)
        )
    )

    return compatible_pairs


# ============================================================================
# 10. SYMBOLIC RELATION CROSS-CHECK
# ============================================================================

def symbolic_relation_audit():

    print()
    print("=" * 78)
    print(
        "10. SYMBOLIC RELATION CROSS-CHECK"
    )
    print("=" * 78)

    failures = 0

    for item in DISCOVERIES:

        relation = item["relation"]

        residual = relation_symbolic_residual(
            relation,
            len(relation) - 1,
        )

        forced = clean(
            item["forced_Z"]
        )

        evaluated = clean(
            residual.subs(
                Z,
                forced,
            )
        )

        print()
        print(
            "  {}".format(
                item["name"]
            )
        )

        print(
            "    residual={}".format(
                residual
            )
        )

        print(
            "    residual_at_forced_Z={}".format(
                evaluated
            )
        )

        if evaluated != 0:
            failures += 1

    print()
    print(
        "  cross_check_failures={}".format(
            failures
        )
    )

    return failures


# ============================================================================
# 11. STRONG VERDICT
# ============================================================================

def strong_verdict(
    classes,
    integer_discoveries,
):

    print()
    print("=" * 78)
    print(
        "11. STRONG CROSS-DISCOVERY VERDICT"
    )
    print("=" * 78)

    discovery_count = len(
        DISCOVERIES
    )

    distinct_count = len(
        classes
    )

    if discovery_count == 0:

        verdict = (
            "NO_DISCOVERIES"
        )

    elif distinct_count == 1:

        if len(integer_discoveries) == discovery_count:
            verdict = (
                "COMMON_INTEGER_Z_SURVIVES"
            )
        else:
            verdict = (
                "COMMON_RATIONAL_Z_BUT_NOT_INTEGER"
            )

    else:

        verdict = (
            "CROSS-DISCOVERY_INCONSISTENCY"
        )

    print(
        "  discovery_count={}".format(
            discovery_count
        )
    )

    print(
        "  distinct_forced_Z_count={}".format(
            distinct_count
        )
    )

    print(
        "  verdict={}".format(
            verdict
        )
    )

    if verdict == (
        "CROSS-DISCOVERY_INCONSISTENCY"
    ):

        print(
            "  interpretation="
            "the transformed integer laws cannot all describe the same "
            "missing source cell simultaneously"
        )

    elif verdict == (
        "COMMON_INTEGER_Z_SURVIVES"
    ):

        print(
            "  interpretation="
            "multiple independent transformed laws agree on one integer Z"
        )

    elif verdict == (
        "COMMON_RATIONAL_Z_BUT_NOT_INTEGER"
    ):

        print(
            "  interpretation="
            "multiple transformed laws agree, but only at a non-integer Z"
        )

    return verdict


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 364R — EXACT CROSS-TRANSFORM Z-CONSISTENCY / "
        "COMMON-ANNIHILATOR COMPATIBILITY AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "  discovery_count={}".format(
            len(DISCOVERIES)
        )
    )

    discovery_audit()

    equality = equality_matrix()

    equal_pairs, nonzero_count = (
        pairwise_difference_audit()
    )

    classes = compatibility_classes()

    integer_discoveries = (
        integer_feasibility()
    )

    denominator_profile()

    common_rational_audit()

    transform_groups = (
        transform_consistency()
    )

    compatible_pairs = (
        pairwise_common_law_audit()
    )

    cross_check_failures = (
        symbolic_relation_audit()
    )

    verdict = strong_verdict(
        classes,
        integer_discoveries,
    )

    # =========================================================================
    # FINAL EXACTNESS
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "12. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  discoveries_reused_exactly={}".format(
            True
        )
    )

    print(
        "  equality_audit_completed=True"
    )

    print(
        "  pairwise_difference_audit_completed=True"
    )

    print(
        "  compatibility_class_audit_completed=True"
    )

    print(
        "  integer_compatibility_audit_completed=True"
    )

    print(
        "  denominator_profile_completed=True"
    )

    print(
        "  transform_internal_consistency_completed=True"
    )

    print(
        "  pairwise_common_law_audit_completed=True"
    )

    print(
        "  symbolic_relation_cross_check_failures={}".format(
            cross_check_failures
        )
    )

    print(
        "  compatible_pair_count={}".format(
            len(compatible_pairs)
        )
    )

    print(
        "  distinct_forced_Z_count={}".format(
            len(classes)
        )
    )

    print(
        "  strong_verdict={}".format(
            verdict
        )
    )

    print(
        "  missing_value_inserted=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures={}".format(
            cross_check_failures
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            cross_check_failures == 0
        )
    )

    print()
    print(
        "EXPERIMENT 364R COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc
            )
        )

        raise
