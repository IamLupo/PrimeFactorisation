#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;

struct CaseData {
    u64 p = 0;
    u64 q = 0;

    std::vector<u64> q_digits;
    std::vector<u64> radix;
    std::vector<u64> weights;
    std::vector<u64> powers;

    u64 tuple_count = 0;
    u64 interval_count = 0;
};

bool safe_add(u64 a, u64 b, u64& out) {
    if (b > UINT64_MAX - a) {
        return false;
    }

    out = a + b;
    return true;
}

bool safe_mul(u64 a, u64 b, u64& out) {
    if (a != 0 && b > UINT64_MAX / a) {
        return false;
    }

    out = a * b;
    return true;
}

bool build_case(
    u64 p,
    u64 q,
    CaseData& c
) {
    c = CaseData{};

    c.p = p;
    c.q = q;

    if (p < 2 || q == 0) {
        return false;
    }

    u64 temp = q;

    while (true) {
        const u64 d = temp % p;

        u64 radix_value;

        if (!safe_add(d, 1, radix_value)) {
            return false;
        }

        c.q_digits.push_back(d);
        c.radix.push_back(radix_value);

        temp /= p;

        if (temp == 0) {
            break;
        }
    }

    c.weights.resize(c.radix.size());
    c.powers.resize(c.radix.size());

    c.weights[0] = 1;
    c.powers[0] = 1;

    for (size_t i = 1; i < c.radix.size(); ++i) {
        if (!safe_mul(
                c.weights[i - 1],
                c.radix[i - 1],
                c.weights[i]
            )) {
            return false;
        }

        if (!safe_mul(
                c.powers[i - 1],
                p,
                c.powers[i]
            )) {
            return false;
        }
    }

    c.tuple_count = 1;

    for (u64 r : c.radix) {
        if (!safe_mul(
                c.tuple_count,
                r,
                c.tuple_count
            )) {
            return false;
        }
    }

    c.interval_count =
        c.tuple_count - 1;

    return true;
}

bool unrank(
    u64 k,
    const CaseData& c,
    u64& j
) {
    if (k >= c.tuple_count) {
        return false;
    }

    j = 0;

    for (size_t i = 0;
         i < c.radix.size();
         ++i) {

        const u64 d =
            (k / c.weights[i]) %
            c.radix[i];

        u64 contribution;

        if (!safe_mul(
                d,
                c.powers[i],
                contribution
            )) {
            return false;
        }

        if (!safe_add(
                j,
                contribution,
                j
            )) {
            return false;
        }
    }

    return true;
}

void print_case_state(
    const CaseData& c,
    const char* label
) {
    std::cout << label << '\n';

    std::cout
        << "p=" << c.p
        << " q=" << c.q
        << '\n';

    std::cout
        << "q_digits:";

    for (u64 x : c.q_digits) {
        std::cout << ' ' << x;
    }

    std::cout << '\n';

    std::cout
        << "radix:";

    for (u64 x : c.radix) {
        std::cout << ' ' << x;
    }

    std::cout << '\n';

    std::cout
        << "weights:";

    for (u64 x : c.weights) {
        std::cout << ' ' << x;
    }

    std::cout << '\n';

    std::cout
        << "powers:";

    for (u64 x : c.powers) {
        std::cout << ' ' << x;
    }

    std::cout << '\n';

    std::cout
        << "tuple_count="
        << c.tuple_count
        << '\n';

    std::cout
        << "interval_count="
        << c.interval_count
        << "\n\n";
}

bool compare_state(
    const CaseData& a,
    const CaseData& b
) {
    return a.p == b.p &&
           a.q == b.q &&
           a.q_digits == b.q_digits &&
           a.radix == b.radix &&
           a.weights == b.weights &&
           a.powers == b.powers &&
           a.tuple_count == b.tuple_count &&
           a.interval_count == b.interval_count;
}

bool test_repeated_case() {
    CaseData a;
    CaseData b;

    if (!build_case(5, 100, a)) {
        return false;
    }

    if (!build_case(5, 100, b)) {
        return false;
    }

    print_case_state(
        a,
        "FIRST BUILD"
    );

    print_case_state(
        b,
        "SECOND BUILD"
    );

    const bool state_equal =
        compare_state(a, b);

    std::cout
        << "repeated_state_equal="
        << state_equal
        << '\n';

    const u64 expected[] = {
        0,
        25,
        50,
        75
    };

    bool unrank_pass = true;

    for (u64 k = 0;
         k < 4;
         ++k) {

        u64 ja;
        u64 jb;

        const bool oka =
            unrank(k, a, ja);

        const bool okb =
            unrank(k, b, jb);

        const bool pass =
            oka &&
            okb &&
            ja == expected[k] &&
            jb == expected[k];

        std::cout
            << "k=" << k
            << " ja=" << ja
            << " jb=" << jb
            << " expected=" << expected[k]
            << " pass=" << pass
            << '\n';

        if (!pass) {
            unrank_pass = false;
        }
    }

    std::cout
        << "unrank_pass="
        << unrank_pass
        << '\n';

    return state_equal &&
           unrank_pass;
}

bool test_multiple_builds() {
    bool pass = true;

    for (int n = 0; n < 1000; ++n) {
        CaseData c;

        if (!build_case(
                5,
                100,
                c
            )) {
            pass = false;
            break;
        }

        const u64 expected[] = {
            0,
            25,
            50,
            75
        };

        for (u64 k = 0;
             k < 4;
             ++k) {

            u64 j;

            if (!unrank(
                    k,
                    c,
                    j
                )) {
                pass = false;
                break;
            }

            if (j != expected[k]) {
                std::cout
                    << "FAIL iteration="
                    << n
                    << " k="
                    << k
                    << " j="
                    << j
                    << " expected="
                    << expected[k]
                    << '\n';

                pass = false;
                break;
            }
        }

        if (!pass) {
            break;
        }
    }

    std::cout
        << "1000_rebuilds_pass="
        << pass
        << '\n';

    return pass;
}

int main() {
    std::cout
        << "START EXPERIMENT 222R\n\n";

    const bool repeated_pass =
        test_repeated_case();

    std::cout << '\n';

    const bool rebuild_pass =
        test_multiple_builds();

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << (
            repeated_pass &&
            rebuild_pass
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 222R\n";

    return (
        repeated_pass &&
        rebuild_pass
    ) ? 0 : 1;
}
