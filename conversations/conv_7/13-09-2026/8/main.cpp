#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 lo;
    u64 hi;
};

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

bool lucas_hit(u64 p, u64 m, u64 t) {
    while (m > 0 || t > 0) {
        const u64 md = m % p;
        const u64 td = t % p;

        if (td > md) {
            return true;
        }

        m /= p;
        t /= p;
    }

    return false;
}

std::vector<Interval> build_intervals(u64 p, u64 m) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = lucas_hit(p, m, t);

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            intervals.push_back({start, t - 1});
            inside = false;
        }
    }

    if (inside) {
        intervals.push_back({start, m});
    }

    return intervals;
}

u64 start_gcd(const std::vector<Interval>& intervals) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        g = std::gcd(
            g,
            intervals[i].lo - intervals[i - 1].lo
        );
    }

    return g;
}

int power_exponent(u64 value, u64 p) {
    int e = 0;

    if (value == 0) {
        return 0;
    }

    while (value % p == 0) {
        value /= p;
        ++e;
    }

    if (value != 1) {
        return -1;
    }

    return e;
}

std::vector<int> digits_lsf(u64 n, u64 p) {
    std::vector<int> digits;

    while (n > 0) {
        digits.push_back(
            static_cast<int>(n % p)
        );

        n /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_string(
    const std::vector<int>& digits
) {
    std::string s;

    for (auto it = digits.rbegin();
         it != digits.rend();
         ++it) {

        s += std::to_string(*it);
    }

    return s;
}

std::string digits_lsf_string(
    const std::vector<int>& digits
) {
    std::string s = "[";

    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i != 0) {
            s += ",";
        }

        s += std::to_string(digits[i]);
    }

    s += "]";

    return s;
}

/*
    Generate every digit word of exactly 'length' digits
    in base p, including leading zeroes.

    The most-significant digit is forced to be non-zero.
*/
void generate_words(
    u64 p,
    int position,
    int length,
    std::vector<int>& digits,
    std::vector<std::vector<int>>& words
) {
    if (position == length) {
        if (digits[0] != 0) {
            words.push_back(digits);
        }

        return;
    }

    for (int d = 0; d < static_cast<int>(p); ++d) {
        digits[position] = d;

        generate_words(
            p,
            position + 1,
            length,
            digits,
            words
        );
    }
}

u64 digits_to_value(
    const std::vector<int>& digits,
    u64 p
) {
    u64 value = 0;

    for (int digit : digits) {
        value = value * p +
                static_cast<u64>(digit);
    }

    return value;
}

void classify_binary_words() {
    std::cout << "\nPHASE 1: COMPLETE BINARY WORDS\n";

    const u64 p = 2;

    for (int length = 2; length <= 9; ++length) {
        std::cout << "\nBINARY LENGTH " << length << "\n";

        std::vector<std::vector<int>> words;

        std::vector<int> digits(length);

        generate_words(
            p,
            0,
            length,
            digits,
            words
        );

        std::vector<std::vector<std::string>> groups(10);

        for (const auto& word_msb : words) {
            const u64 m = digits_to_value(word_msb, p);

            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            const u64 g = start_gcd(intervals);

            const int e = power_exponent(g, p);

            if (e >= 0 &&
                e < static_cast<int>(groups.size())) {

                groups[e].push_back(
                    digits_string(
                        digits_lsf(m, p)
                    )
                );
            }
        }

        for (int e = 1;
             e < static_cast<int>(groups.size());
             ++e) {

            if (groups[e].empty()) {
                continue;
            }

            std::cout
                << "e=" << e
                << " count=" << groups[e].size()
                << "\n";

            std::cout << "patterns=";

            const std::size_t limit =
                std::min<std::size_t>(
                    groups[e].size(),
                    80
                );

            for (std::size_t i = 0;
                 i < limit;
                 ++i) {

                if (i != 0) {
                    std::cout << ",";
                }

                std::cout << groups[e][i];
            }

            if (groups[e].size() > limit) {
                std::cout << ",...";
            }

            std::cout << "\n";
        }
    }
}

void classify_ternary_words() {
    std::cout << "\nPHASE 2: COMPLETE TERNARY WORDS\n";

    const u64 p = 3;

    for (int length = 2; length <= 6; ++length) {
        std::cout << "\nTERNARY LENGTH " << length << "\n";

        std::vector<std::vector<int>> words;

        std::vector<int> digits(length);

        generate_words(
            p,
            0,
            length,
            digits,
            words
        );

        std::vector<std::vector<std::string>> groups(10);

        for (const auto& word_msb : words) {
            const u64 m =
                digits_to_value(word_msb, p);

            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            const u64 g = start_gcd(intervals);

            const int e =
                power_exponent(g, p);

            if (e >= 0 &&
                e < static_cast<int>(groups.size())) {

                groups[e].push_back(
                    digits_string(
                        digits_lsf(m, p)
                    )
                );
            }
        }

        for (int e = 1;
             e < static_cast<int>(groups.size());
             ++e) {

            if (groups[e].empty()) {
                continue;
            }

            std::cout
                << "e=" << e
                << " count=" << groups[e].size()
                << "\n";

            std::cout << "patterns=";

            const std::size_t limit =
                std::min<std::size_t>(
                    groups[e].size(),
                    100
                );

            for (std::size_t i = 0;
                 i < limit;
                 ++i) {

                if (i != 0) {
                    std::cout << ",";
                }

                std::cout << groups[e][i];
            }

            if (groups[e].size() > limit) {
                std::cout << ",...";
            }

            std::cout << "\n";
        }
    }
}

void test_candidate_digit_statistics() {
    std::cout
        << "\nPHASE 3: DIGIT-STATISTIC CORRELATION\n";

    struct Stats {
        int leading_pminus1 = 0;
        int trailing_zeroes = 0;
        int first_zero_to_nonzero = -1;
        int first_nonzero_to_lower = -1;
        int first_change = -1;
        int last_non_pminus1 = -1;
    };

    auto compute_stats =
        [](const std::vector<int>& d, int p) {

            Stats s;

            /*
                d is least-significant first.
            */

            int n =
                static_cast<int>(d.size());

            while (
                s.leading_pminus1 < n &&
                d[s.leading_pminus1] == p - 1
            ) {
                ++s.leading_pminus1;
            }

            while (
                s.trailing_zeroes < n &&
                d[s.trailing_zeroes] == 0
            ) {
                ++s.trailing_zeroes;
            }

            for (int i = 1; i < n; ++i) {
                if (d[i - 1] == 0 &&
                    d[i] > 0) {

                    if (s.first_zero_to_nonzero == -1) {
                        s.first_zero_to_nonzero = i;
                    }
                }

                if (d[i] < d[i - 1]) {
                    if (s.first_nonzero_to_lower == -1) {
                        s.first_nonzero_to_lower = i;
                    }
                }

                if (d[i] != d[i - 1]) {
                    if (s.first_change == -1) {
                        s.first_change = i;
                    }
                }
            }

            for (int i = n - 1; i >= 0; --i) {
                if (d[i] != p - 1) {
                    s.last_non_pminus1 = i;
                    break;
                }
            }

            return s;
        };

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        std::cout << "\nSTATISTICS p=" << p << "\n";

        u64 total = 0;
        u64 exact_zero_to_nonzero = 0;
        u64 exact_trailing_zeroes = 0;
        u64 exact_last_non_p1 = 0;
        u64 exact_first_change = 0;

        for (u64 m = 1; m <= 3000; ++m) {
            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 g =
                start_gcd(intervals);

            const int e =
                power_exponent(g, p);

            const auto digits =
                digits_lsf(m, p);

            const Stats s =
                compute_stats(digits, static_cast<int>(p));

            if (s.first_zero_to_nonzero == e) {
                ++exact_zero_to_nonzero;
            }

            if (s.trailing_zeroes == e) {
                ++exact_trailing_zeroes;
            }

            if (s.last_non_pminus1 + 1 == e) {
                ++exact_last_non_p1;
            }

            if (s.first_change == e) {
                ++exact_first_change;
            }
        }

        std::cout
            << "multi_interval_cases="
            << total << "\n";

        std::cout
            << "first_zero_to_nonzero_exact="
            << exact_zero_to_nonzero
            << "/" << total << "\n";

        std::cout
            << "trailing_zeroes_exact="
            << exact_trailing_zeroes
            << "/" << total << "\n";

        std::cout
            << "last_non_pminus1_exact="
            << exact_last_non_p1
            << "/" << total << "\n";

        std::cout
            << "first_change_exact="
            << exact_first_change
            << "/" << total << "\n";
    }
}

void run_direct_examples() {
    std::cout
        << "\nPHASE 4: DIRECT EXAMPLES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 6},
        {2, 12},
        {2, 13},
        {2, 20},
        {2, 24},
        {2, 25},
        {2, 27},
        {2, 28},
        {2, 40},
        {2, 48},
        {2, 96},
        {3, 18},
        {3, 23},
        {3, 54},
        {5, 50},
        {5, 59},
        {5, 75},
        {7, 98},
        {11, 242}
    };

    for (const auto& c : cases) {
        const auto intervals =
            build_intervals(c.p, c.m);

        const u64 g =
            start_gcd(intervals);

        const int e =
            power_exponent(g, c.p);

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_lsf_string(
                digits_lsf(c.m, c.p)
            )
            << " intervals="
            << intervals.size()
            << " G=" << g
            << " e=" << e
            << "\n";

        std::cout << "starts=";

        const std::size_t limit =
            std::min<std::size_t>(
                intervals.size(),
                12
            );

        for (std::size_t i = 0;
             i < limit;
             ++i) {

            if (i != 0) {
                std::cout << ",";
            }

            std::cout << intervals[i].lo;
        }

        if (intervals.size() > limit) {
            std::cout << ",...";
        }

        std::cout << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 164;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    classify_binary_words();
    classify_ternary_words();
    test_candidate_digit_statistics();
    run_direct_examples();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
