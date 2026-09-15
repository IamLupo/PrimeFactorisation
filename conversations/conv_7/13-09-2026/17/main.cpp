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

struct Structure {
    int a;
    int z;
    int e;

    u64 b;
    u64 s0;
    u64 modulus;
    u64 q;
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

bool lucas_hit(
    u64 p,
    u64 m,
    u64 t
) {
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

Structure analyze_structure(
    u64 p,
    u64 m
) {
    Structure s{};

    u64 x = m;

    while (x > 0 &&
           x % p == p - 1) {

        ++s.a;
        x /= p;
    }

    s.b = x % p;

    u64 pa = 1;

    for (int i = 0;
         i < s.a;
         ++i) {
        pa *= p;
    }

    x /= p;

    while (x > 0 &&
           x % p == 0) {

        ++s.z;
        x /= p;
    }

    s.e =
        s.a +
        1 +
        s.z;

    s.s0 =
        (s.b + 1) * pa;

    s.modulus = 1;

    for (int i = 0;
         i < s.e;
         ++i) {
        s.modulus *= p;
    }

    s.q =
        m / s.modulus;

    return s;
}

bool digitwise_leq(
    u64 p,
    u64 a,
    u64 b
) {
    while (a > 0 || b > 0) {
        const u64 ad = a % p;
        const u64 bd = b % p;

        if (ad > bd) {
            return false;
        }

        a /= p;
        b /= p;
    }

    return true;
}

/*
    All predicted interval starts.
*/
std::vector<u64> predicted_starts(
    u64 p,
    const Structure& s
) {
    std::vector<u64> starts;

    for (u64 j = 0;
         j < s.q;
         ++j) {

        if (!digitwise_leq(
                p,
                j,
                s.q
            )) {
            continue;
        }

        starts.push_back(
            s.s0 +
            j * s.modulus
        );
    }

    return starts;
}

/*
    Actual interval construction.

    Used only for small exhaustive validation.
*/
std::vector<Interval> actual_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1;
         t <= m;
         ++t) {

        const bool hit =
            lucas_hit(
                p,
                m,
                t
            );

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            result.push_back({
                start,
                t - 1
            });

            inside = false;
        }
    }

    if (inside) {
        result.push_back({
            start,
            m
        });
    }

    return result;
}

/*
    For a predicted start s, find its actual endpoint.

    This is deliberately only used for small cases.
*/
u64 actual_endpoint(
    u64 p,
    u64 m,
    u64 start
) {
    u64 t = start;

    while (t <= m &&
           lucas_hit(p, m, t)) {
        ++t;
    }

    return t - 1;
}

/*
    The low e digits of an interval start are fixed.
    Therefore the first failure after a start should occur
    when the increment changes the relevant digit pattern.

    We experimentally derive the endpoint by looking at
    the first non-HIT after the start.
*/
u64 predicted_endpoint_candidate(
    u64 p,
    u64 m,
    u64 start
) {
    /*
        Current hypothesis:

        interval length is determined by the largest
        low-order digit block that remains inside the
        HIT condition.

        We calculate the endpoint directly from the
        predecessor digit condition.

        This function is intentionally kept separate so
        a failed formula is easy to diagnose.
    */

    u64 t = start;

    while (t < m) {
        if (!lucas_hit(
                p,
                m,
                t + 1
            )) {
            break;
        }

        ++t;
    }

    return t;
}

/*
    Generate a compact digit representation.
*/
std::vector<int> digits(
    u64 x,
    u64 p
) {
    std::vector<int> d;

    while (x > 0) {
        d.push_back(
            static_cast<int>(x % p)
        );

        x /= p;
    }

    if (d.empty()) {
        d.push_back(0);
    }

    return d;
}

std::string digits_string(
    const std::vector<int>& d
) {
    std::string s = "[";

    for (std::size_t i = 0;
         i < d.size();
         ++i) {

        if (i != 0) {
            s += ",";
        }

        s += std::to_string(d[i]);
    }

    s += "]";

    return s;
}

std::string intervals_string(
    const std::vector<Interval>& intervals,
    std::size_t limit = 20
) {
    std::string s;

    const std::size_t n =
        std::min(
            intervals.size(),
            limit
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        if (i != 0) {
            s += " ";
        }

        s += "[";
        s += std::to_string(intervals[i].lo);
        s += ",";
        s += std::to_string(intervals[i].hi);
        s += "]";
    }

    if (intervals.size() > limit) {
        s += " ...";
    }

    return s;
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: COMPLETE INTERVAL START VALIDATION\n";

    u64 total = 0;
    u64 start_pass = 0;
    u64 count_pass = 0;
    u64 endpoint_pass = 0;
    u64 complete_pass = 0;

    u64 failures = 0;

    for (u64 p = 2;
         p <= 31;
         ++p) {

        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1;
             m <= 1500;
             ++m) {

            const Structure s =
                analyze_structure(
                    p,
                    m
                );

            const auto actual =
                actual_intervals(
                    p,
                    m
                );

            if (actual.size() < 2) {
                continue;
            }

            ++total;

            const auto predicted =
                predicted_starts(
                    p,
                    s
                );

            bool starts_ok = true;

            if (predicted.size() !=
                actual.size()) {

                starts_ok = false;
            } else {

                for (std::size_t i = 0;
                     i < actual.size();
                     ++i) {

                    if (actual[i].lo !=
                        predicted[i]) {

                        starts_ok = false;
                        break;
                    }
                }
            }

            if (starts_ok) {
                ++start_pass;
            }

            if (predicted.size() ==
                actual.size()) {

                ++count_pass;
            }

            /*
                Endpoint validation.

                Currently this simply establishes the
                actual endpoint function for every predicted
                start. This is the data from which the next
                endpoint formula can be inferred.
            */
            bool all_endpoint_hits = true;

            for (u64 start : predicted) {
                const u64 endpoint =
                    actual_endpoint(
                        p,
                        m,
                        start
                    );

                if (endpoint < start) {
                    all_endpoint_hits = false;
                    break;
                }
            }

            if (all_endpoint_hits) {
                ++endpoint_pass;
            }

            if (starts_ok &&
                all_endpoint_hits) {

                ++complete_pass;
            } else {
                ++failures;

                if (failures <= 20) {
                    std::cout
                        << "\nFAILURE\n"
                        << "p=" << p
                        << " m=" << m
                        << " digits="
                        << digits_string(
                            digits(m, p)
                        )
                        << "\n"
                        << "actual="
                        << intervals_string(actual)
                        << "\n"
                        << "predicted_starts=";

                    for (u64 x : predicted) {
                        std::cout
                            << x << " ";
                    }

                    std::cout << "\n";
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "multi_interval_cases="
        << total
        << "\n"
        << "start_set_pass="
        << start_pass
        << "/" << total
        << "\n"
        << "count_pass="
        << count_pass
        << "/" << total
        << "\n"
        << "endpoint_validity_pass="
        << endpoint_pass
        << "/" << total
        << "\n"
        << "complete_pass="
        << complete_pass
        << "/" << total
        << "\n"
        << "failures="
        << failures
        << "\n";
}

void run_endpoint_pattern_phase() {
    std::cout
        << "\nPHASE 2: ENDPOINT PATTERNS\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 10},
        {2, 12},
        {2, 20},
        {2, 24},
        {2, 28},
        {2, 40},
        {2, 48},

        {3, 22},
        {3, 41},
        {3, 59},
        {3, 68},
        {3, 125},
        {3, 377},

        {5, 59},
        {5, 94},
        {5, 194},
        {5, 219}
    };

    for (const auto& c : cases) {
        const Structure s =
            analyze_structure(
                c.p,
                c.m
            );

        const auto intervals =
            actual_intervals(
                c.p,
                c.m
            );

        if (intervals.size() < 2) {
            continue;
        }

        std::cout
            << "\np=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_string(
                digits(c.m, c.p)
            )
            << "\n";

        std::cout
            << "s0=" << s.s0
            << " modulus=" << s.modulus
            << " q=" << s.q
            << "\n";

        std::cout
            << "intervals="
            << intervals_string(
                intervals
            )
            << "\n";

        std::cout
            << "lengths=";

        for (std::size_t i = 0;
             i < intervals.size() &&
             i < 20;
             ++i) {

            if (i != 0) {
                std::cout << ",";
            }

            std::cout
                << intervals[i].hi -
                intervals[i].lo +
                1;
        }

        std::cout
            << "\n";
    }
}

void run_large_start_phase() {
    std::cout
        << "\nPHASE 3: LARGE START FORMULA\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37, 43,
        53, 67, 101,
        127, 211, 431,
        1009, 2003
    };

    constexpr int CASES = 50000;

    u64 state =
        0x1722026ULL;

    u64 tested = 0;
    u64 passed = 0;
    u64 failed = 0;

    for (int i = 0;
         i < CASES;
         ++i) {

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 p =
            primes[
                state %
                primes.size()
            ];

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 m =
            1 +
            state %
            1000000000000000000ULL;

        const Structure s =
            analyze_structure(
                p,
                m
            );

        const auto starts =
            predicted_starts(
                p,
                s
            );

        if (starts.empty()) {
            continue;
        }

        ++tested;

        bool ok = true;

        for (u64 start : starts) {
            if (!lucas_hit(
                    p,
                    m,
                    start
                )) {

                ok = false;
                break;
            }

            if (start > 0 &&
                lucas_hit(
                    p,
                    m,
                    start - 1
                )) {

                ok = false;
                break;
            }
        }

        if (ok) {
            ++passed;
        } else {
            ++failed;

            if (failed <= 10) {
                std::cout
                    << "FAILURE "
                    << "p=" << p
                    << " m=" << m
                    << "\n";
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "tested="
        << tested
        << "\n"
        << "start_formula_pass="
        << passed
        << "/" << tested
        << "\n"
        << "failures="
        << failed
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 172;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_endpoint_pattern_phase();
    run_large_start_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
