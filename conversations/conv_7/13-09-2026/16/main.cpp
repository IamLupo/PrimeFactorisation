#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;

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

/*
    j <=_p q

    means every base-p digit of j is <=
    the corresponding digit of q.
*/
bool digitwise_leq(
    u64 p,
    u64 j,
    u64 q
) {
    while (j > 0 || q > 0) {
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

/*
    Collect actual HIT interval starts.
*/
std::vector<u64> collect_starts(
    u64 p,
    u64 m
) {
    std::vector<u64> starts;

    bool previous_hit = false;

    for (u64 t = 1;
         t <= m;
         ++t) {

        const bool hit =
            lucas_hit(
                p,
                m,
                t
            );

        if (hit && !previous_hit) {
            starts.push_back(t);
        }

        previous_hit = hit;
    }

    return starts;
}

std::vector<u64> actual_indices(
    u64 p,
    u64 m,
    const Structure& s
) {
    const auto starts =
        collect_starts(
            p,
            m
        );

    std::vector<u64> result;

    for (u64 start : starts) {
        result.push_back(
            (start - s.s0) /
            s.modulus
        );
    }

    return result;
}

/*
    Predicted index set:

        j < q

        and

        j <=_p q
*/
std::vector<u64> predicted_indices(
    u64 p,
    const Structure& s
) {
    std::vector<u64> result;

    for (u64 j = 0;
         j < s.q;
         ++j) {

        if (digitwise_leq(
                p,
                j,
                s.q
            )) {

            result.push_back(j);
        }
    }

    return result;
}

std::string vector_string(
    const std::vector<u64>& values,
    std::size_t limit = 40
) {
    std::string s = "[";

    const std::size_t n =
        std::min(
            values.size(),
            limit
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        if (i != 0) {
            s += ",";
        }

        s +=
            std::to_string(values[i]);
    }

    if (values.size() > limit) {
        s += ",...";
    }

    s += "]";

    return s;
}

bool vectors_equal(
    const std::vector<u64>& a,
    const std::vector<u64>& b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0;
         i < a.size();
         ++i) {

        if (a[i] != b[i]) {
            return false;
        }
    }

    return true;
}

void print_failure(
    u64 p,
    u64 m,
    const Structure& s,
    const std::vector<u64>& actual,
    const std::vector<u64>& predicted
) {
    std::cout
        << "\nFAILURE\n";

    std::cout
        << "p=" << p
        << " m=" << m
        << "\n";

    std::cout
        << "a=" << s.a
        << " b=" << s.b
        << " z=" << s.z
        << " e=" << s.e
        << "\n";

    std::cout
        << "s0=" << s.s0
        << " modulus=" << s.modulus
        << " q=" << s.q
        << "\n";

    std::cout
        << "actual="
        << vector_string(actual)
        << "\n";

    std::cout
        << "predicted="
        << vector_string(predicted)
        << "\n";
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: COMPLETE INDEX-SET TEST\n";

    u64 total = 0;
    u64 pass = 0;
    u64 fail = 0;

    u64 decomposition_pass = 0;
    u64 decomposition_fail = 0;

    u64 count_pass = 0;

    for (u64 p = 2;
         p <= 37;
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
                actual_indices(
                    p,
                    m,
                    s
                );

            if (actual.size() < 2) {
                continue;
            }

            ++total;

            /*
                First verify the exact decomposition

                    m = q*p^e + s0 - 1
            */
            const bool decomposition =
                m ==
                s.q * s.modulus +
                s.s0 -
                1;

            if (decomposition) {
                ++decomposition_pass;
            } else {
                ++decomposition_fail;
            }

            /*
                Compare the complete index set.
            */
            const auto predicted =
                predicted_indices(
                    p,
                    s
                );

            if (vectors_equal(
                    actual,
                    predicted
                )) {

                ++pass;
            } else {
                ++fail;

                if (fail <= 20) {
                    print_failure(
                        p,
                        m,
                        s,
                        actual,
                        predicted
                    );
                }
            }

            /*
                Count prediction:

                Number of digitwise subnumbers of q is

                    product_i (q_i + 1)

                but j=q itself is excluded, hence:

                    predicted_count
                    =
                    product_i(q_i+1) - 1

                because j<q.
            */
            u64 count = 1;

            u64 x = s.q;

            while (x > 0) {
                const u64 digit =
                    x % p;

                x /= p;

                count *=
                    digit + 1;
            }

            if (count > 0) {
                --count;
            }

            if (count ==
                predicted.size()) {

                ++count_pass;
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n";

    std::cout
        << "multi_start_cases="
        << total
        << "\n";

    std::cout
        << "complete_index_set_pass="
        << pass
        << "/" << total
        << "\n";

    std::cout
        << "complete_index_set_fail="
        << fail
        << "\n";

    std::cout
        << "decomposition_pass="
        << decomposition_pass
        << "/" << total
        << "\n";

    std::cout
        << "decomposition_fail="
        << decomposition_fail
        << "\n";

    std::cout
        << "digit_count_formula_pass="
        << count_pass
        << "/" << total
        << "\n";
}

void run_examples() {
    std::cout
        << "\nPHASE 2: INDEX-SET EXAMPLES\n";

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

        const auto actual =
            actual_indices(
                c.p,
                c.m,
                s
            );

        if (actual.size() < 2) {
            continue;
        }

        const auto predicted =
            predicted_indices(
                c.p,
                s
            );

        std::cout
            << "\np=" << c.p
            << " m=" << c.m
            << " e=" << s.e
            << " s0=" << s.s0
            << " p^e=" << s.modulus
            << " q=" << s.q
            << "\n";

        std::cout
            << "actual="
            << vector_string(actual)
            << "\n";

        std::cout
            << "predicted="
            << vector_string(predicted)
            << "\n";

        std::cout
            << "match="
            << (
                vectors_equal(
                    actual,
                    predicted
                )
                    ? "YES"
                    : "NO"
            )
            << "\n";
    }
}

/*
    Large validation does not enumerate j.

    It randomly samples j values from [0,q) and verifies that
    the digitwise condition exactly matches whether

        s0 + j*p^e

    is a HIT interval start.

    The actual start condition is:

        current point = HIT
        predecessor = MISS
*/
void run_large_random_phase() {
    std::cout
        << "\nPHASE 3: LARGE RANDOM INDEX TEST\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009,
        2003, 4001
    };

    constexpr int CASES = 50000;
    constexpr int J_PER_CASE = 10;

    u64 state =
        0x1712026ULL;

    u64 tested = 0;
    u64 pass = 0;
    u64 fail = 0;

    for (int c = 0;
         c < CASES;
         ++c) {

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

        if (s.q < 2) {
            continue;
        }

        for (int k = 0;
             k < J_PER_CASE;
             ++k) {

            state =
                state *
                6364136223846793005ULL +
                1442695040888963407ULL;

            const u64 j =
                state % s.q;

            const u64 start =
                s.s0 +
                j * s.modulus;

            const bool actual_start =
                lucas_hit(
                    p,
                    m,
                    start
                ) &&
                !lucas_hit(
                    p,
                    m,
                    start - 1
                );

            const bool predicted =
                digitwise_leq(
                    p,
                    j,
                    s.q
                );

            ++tested;

            if (actual_start ==
                predicted) {

                ++pass;
            } else {
                ++fail;

                if (fail <= 20) {
                    std::cout
                        << "FAILURE "
                        << "p=" << p
                        << " m=" << m
                        << " q=" << s.q
                        << " j=" << j
                        << " actual="
                        << actual_start
                        << " predicted="
                        << predicted
                        << "\n";
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n";

    std::cout
        << "tested_index_positions="
        << tested
        << "\n";

    std::cout
        << "predicate_match="
        << pass
        << "/" << tested
        << "\n";

    std::cout
        << "failures="
        << fail
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 171;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_examples();
    run_large_random_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
