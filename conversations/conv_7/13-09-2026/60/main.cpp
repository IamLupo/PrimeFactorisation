#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Interval {
    u64 left;
    u64 right;
};

struct CaseData {
    u64 p;
    u64 a;
    u64 b;
    u64 z;
    u64 e;
    u64 s0;
    u64 q;
    u64 m;
};

bool safe_mul(u64 a, u64 b, u64& out) {
    if (a != 0 && b > UINT64_MAX / a) {
        return false;
    }
    out = a * b;
    return true;
}

bool safe_add(u64 a, u64 b, u64& out) {
    if (b > UINT64_MAX - a) {
        return false;
    }
    out = a + b;
    return true;
}

bool pow_u64(u64 base, u64 exp, u64& out) {
    out = 1;
    for (u64 i = 0; i < exp; ++i) {
        if (!safe_mul(out, base, out)) {
            return false;
        }
    }
    return true;
}

u64 digit(u64 x, u64 p, u64 i) {
    while (i > 0) {
        x /= p;
        --i;
    }
    return x % p;
}

bool leq_p(u64 j, u64 q, u64 p) {
    while (j != 0 || q != 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
    }

    return true;
}

bool build_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    CaseData& c
) {
    if (p < 2) {
        return false;
    }

    if (b > p - 2) {
        return false;
    }

    u64 b1;
    if (!safe_add(b, 1, b1)) {
        return false;
    }

    u64 pa;
    if (!pow_u64(p, a, pa)) {
        return false;
    }

    u64 s0;
    if (!safe_mul(b1, pa, s0)) {
        return false;
    }

    u64 e;
    if (!safe_add(a, z, e)) {
        return false;
    }

    if (!safe_add(e, 1, e)) {
        return false;
    }

    u64 pe;
    if (!pow_u64(p, e, pe)) {
        return false;
    }

    u64 qpe;
    if (!safe_mul(q, pe, qpe)) {
        return false;
    }

    u64 mp1;
    if (!safe_add(s0, qpe, mp1)) {
        return false;
    }

    if (mp1 == 0) {
        return false;
    }

    c.p = p;
    c.a = a;
    c.b = b;
    c.z = z;
    c.e = e;
    c.s0 = s0;
    c.q = q;
    c.m = mp1 - 1;

    return true;
}

std::vector<Interval> generate_intervals(const CaseData& c) {
    const u64 p = c.p;
    const u64 q = c.q;
    const u64 e = c.e;
    const u64 s0 = c.s0;

    u64 pe;
    pow_u64(p, e, pe);

    std::vector<Interval> intervals;

    for (u64 j = 0; j < q; ++j) {
        if (!leq_p(j, q, p)) {
            continue;
        }

        u64 first;
        if (!safe_mul(j, pe, first)) {
            continue;
        }

        if (!safe_add(s0, first, first)) {
            continue;
        }

        u64 r = 0;
        u64 jj = j;
        u64 qq = q;

        while (jj % p == qq % p) {
            ++r;
            jj /= p;
            qq /= p;

            if (jj == 0 && qq == 0) {
                break;
            }
        }

        u64 pr;
        if (!pow_u64(p, r, pr)) {
            continue;
        }

        u64 jr = j;
        for (u64 i = 0; i < r; ++i) {
            jr /= p;
        }

        u64 prefix = jr;
        if (!safe_add(prefix, 1, prefix)) {
            continue;
        }

        u64 power_er;
        if (!safe_mul(pe, pr, power_er)) {
            continue;
        }

        u64 end_plus_1;
        if (!safe_mul(prefix, power_er, end_plus_1)) {
            continue;
        }

        if (end_plus_1 == 0) {
            continue;
        }

        u64 right = end_plus_1 - 1;

        intervals.push_back({first, right});
    }

    std::sort(
        intervals.begin(),
        intervals.end(),
        [](const Interval& x, const Interval& y) {
            if (x.left != y.left) {
                return x.left < y.left;
            }
            return x.right < y.right;
        }
    );

    return intervals;
}

bool check_case(const CaseData& c, bool print_details) {
    const auto intervals = generate_intervals(c);

    const u64 I = static_cast<u64>(intervals.size());

    bool pass = true;

    if (I == 0) {
        return false;
    }

    // 1. First HIT starts exactly at s0.
    const bool first_start_ok = intervals.front().left == c.s0;

    // 2. Last HIT ends exactly at m-s0.
    bool last_end_ok = false;

    if (c.m >= c.s0) {
        const u64 expected_last_right = c.m - c.s0;
        last_end_ok = intervals.back().right == expected_last_right;
    }

    // 3. Every internal gap has length s0.
    bool internal_ok = true;

    u64 internal_gap_checks = 0;
    u64 internal_gap_failures = 0;

    for (size_t i = 1; i < intervals.size(); ++i) {
        ++internal_gap_checks;

        const u64 prev_right = intervals[i - 1].right;
        const u64 next_left = intervals[i].left;

        if (next_left <= prev_right) {
            internal_ok = false;
            ++internal_gap_failures;
            continue;
        }

        const u64 gap = next_left - prev_right - 1;

        if (gap != c.s0) {
            internal_ok = false;
            ++internal_gap_failures;
        }
    }

    // 4. Explicit alternating coverage.
    bool partition_ok = true;

    // First MISS block: [0, s0-1].
    if (c.s0 > 0) {
        if (intervals.front().left != c.s0) {
            partition_ok = false;
        }
    }

    // Each internal boundary must be:
    // next.left = previous.right + s0 + 1.
    for (size_t i = 1; i < intervals.size(); ++i) {
        u64 expected_left;

        if (!safe_add(intervals[i - 1].right, c.s0, expected_left)) {
            partition_ok = false;
            break;
        }

        if (!safe_add(expected_left, 1, expected_left)) {
            partition_ok = false;
            break;
        }

        if (intervals[i].left != expected_left) {
            partition_ok = false;
            break;
        }
    }

    // Final MISS block must have length s0.
    bool final_block_ok = false;

    if (intervals.back().right <= c.m) {
        const u64 tail = c.m - intervals.back().right;
        final_block_ok = tail == c.s0;
    }

    // 5. Complete coverage count.
    u64 total_length = 0;
    for (const auto& iv : intervals) {
        const u64 len = iv.right - iv.left + 1;
        total_length += len;
    }

    u64 miss_blocks = 0;
    if (I <= UINT64_MAX - 1) {
        miss_blocks = I + 1;
    } else {
        pass = false;
    }

    u64 expected_miss = 0;
    bool expected_miss_ok = safe_mul(c.s0, miss_blocks, expected_miss);

    u64 domain_size = c.m + 1;

    u64 expected_hit = 0;
    bool expected_hit_ok =
        expected_miss_ok && expected_miss <= domain_size;

    if (expected_hit_ok) {
        expected_hit = domain_size - expected_miss;
    }

    const bool hit_count_ok =
        expected_hit_ok && total_length == expected_hit;

    // 6. Exact alternating partition:
    // MISS, HIT, MISS, HIT, ..., MISS.
    bool exact_partition = first_start_ok &&
                           last_end_ok &&
                           internal_ok &&
                           partition_ok &&
                           final_block_ok &&
                           hit_count_ok;

    pass = pass && exact_partition;

    if (print_details) {
        std::cout << "CASE\n";
        std::cout << "p=" << c.p
                  << " a=" << c.a
                  << " b=" << c.b
                  << " z=" << c.z
                  << " e=" << c.e
                  << " s0=" << c.s0
                  << " q=" << c.q << '\n';

        std::cout << "m=" << c.m << '\n';
        std::cout << "intervals=" << I << '\n';

        std::cout << "first_start=" << intervals.front().left << '\n';
        std::cout << "expected_first_start=" << c.s0 << '\n';
        std::cout << "first_start_pass=" << first_start_ok << '\n';

        std::cout << "last_right=" << intervals.back().right << '\n';
        std::cout << "expected_last_right=" << (c.m - c.s0) << '\n';
        std::cout << "last_end_pass=" << last_end_ok << '\n';

        std::cout << "internal_gap_checks=" << internal_gap_checks << '\n';
        std::cout << "internal_gap_failures=" << internal_gap_failures << '\n';
        std::cout << "internal_gap_pass=" << internal_ok << '\n';

        std::cout << "partition_pass=" << partition_ok << '\n';
        std::cout << "final_block_pass=" << final_block_ok << '\n';

        std::cout << "hit_length=" << total_length << '\n';
        std::cout << "expected_hit_length=" << expected_hit << '\n';
        std::cout << "hit_count_pass=" << hit_count_ok << '\n';

        std::cout << "final_pass=" << exact_partition << '\n';

        const size_t sample_count = std::min<size_t>(10, intervals.size());

        for (size_t i = 0; i < sample_count; ++i) {
            std::cout << "sample[" << i << "]="
                      << "[" << intervals[i].left
                      << "," << intervals[i].right << "]\n";
        }

        std::cout << '\n';
    }

    return pass;
}

u64 random_bounded(std::mt19937_64& rng, u64 lo, u64 hi) {
    std::uniform_int_distribution<u64> dist(lo, hi);
    return dist(rng);
}

bool run_random_suite(
    u64 p,
    int cases,
    u64 max_a,
    u64 max_b_or_z,
    u64 max_q_digits,
    std::mt19937_64& rng
) {
    int failures = 0;
    int skipped = 0;

    for (int i = 0; i < cases; ++i) {
        const u64 a = random_bounded(rng, 0, max_a);

        const u64 b = random_bounded(rng, 0, p - 2);

        const u64 z = random_bounded(rng, 0, max_b_or_z);

        u64 p_pow_q;
        if (!pow_u64(p, max_q_digits, p_pow_q)) {
            ++skipped;
            continue;
        }

        if (p_pow_q <= 2) {
            ++skipped;
            continue;
        }

        const u64 q = random_bounded(rng, 1, p_pow_q - 1);

        CaseData c;

        if (!build_case(p, a, b, z, q, c)) {
            ++skipped;
            continue;
        }

        // Avoid cases where q generates an astronomically large
        // interval set; this experiment is about exact partitioning.
        if (q > 2000000) {
            ++skipped;
            continue;
        }

        if (!check_case(c, false)) {
            ++failures;

            if (failures <= 3) {
                std::cout << "RANDOM FAILURE\n";
                std::cout << "p=" << c.p
                          << " a=" << c.a
                          << " b=" << c.b
                          << " z=" << c.z
                          << " e=" << c.e
                          << " s0=" << c.s0
                          << " q=" << c.q
                          << " m=" << c.m << "\n\n";
            }
        }
    }

    std::cout << "RANDOM SUITE\n";
    std::cout << "p=" << p
              << " cases=" << cases
              << " failures=" << failures
              << " skipped=" << skipped << '\n';
    std::cout << "suite_pass=" << (failures == 0) << "\n\n";

    return failures == 0;
}

int main() {
    std::cout << "START EXPERIMENT 216\n\n";

    bool all_pass = true;

    {
        CaseData c;
        build_case(2, 0, 0, 2, 63, c);
        all_pass = check_case(c, true) && all_pass;
    }

    {
        CaseData c;
        build_case(3, 0, 0, 3, 80, c);
        all_pass = check_case(c, true) && all_pass;
    }

    {
        CaseData c;
        build_case(5, 0, 0, 2, 100, c);
        all_pass = check_case(c, true) && all_pass;
    }

    {
        CaseData c;
        build_case(3, 3, 2, 1, 80, c);
        all_pass = check_case(c, true) && all_pass;
    }

    // Large dense binary case.
    {
        CaseData c;
        build_case(2, 0, 0, 39, 1048575, c);
        all_pass = check_case(c, false) && all_pass;
    }

    // Large dense base-5 case; chosen to stay inside uint64_t.
    {
        CaseData c;
        build_case(5, 0, 0, 14, 390624, c);
        all_pass = check_case(c, false) && all_pass;
    }

    // Sparse case.
    {
        CaseData c;
        build_case(2, 0, 0, 19, 1073741825ULL, c);
        all_pass = check_case(c, true) && all_pass;
    }

    std::mt19937_64 rng(216216216ULL);

    all_pass =
        run_random_suite(2, 300, 10, 10, 18, rng) &&
        all_pass;

    all_pass =
        run_random_suite(3, 300, 8, 8, 11, rng) &&
        all_pass;

    all_pass =
        run_random_suite(5, 300, 6, 6, 9, rng) &&
        all_pass;

    std::cout << "THEOREM TARGET\n";
    std::cout << "MISS_0 = [0,s0-1]\n";
    std::cout << "MISS_j = [E_{j-1}+1, x_j-1]\n";
    std::cout << "x_j = E_{j-1}+s0+1\n";
    std::cout << "last_right = m-s0\n";
    std::cout << "every MISS block has length s0\n";
    std::cout << "complete partition = MISS-HIT-MISS-HIT-...-HIT-MISS\n\n";

    std::cout << "OVERALL_PASS=" << all_pass << '\n';
    std::cout << "FINISHED EXPERIMENT 216\n";

    return all_pass ? 0 : 1;
}
