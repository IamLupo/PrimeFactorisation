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
    int exponent;
    u64 first_start;
    u64 modulus;
};

bool is_prime(u64 n) {
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) return false;
    }

    return true;
}

bool lucas_hit(u64 p, u64 m, u64 t) {
    while (m > 0 || t > 0) {
        const u64 md = m % p;
        const u64 td = t % p;

        if (td > md) return true;

        m /= p;
        t /= p;
    }

    return false;
}

Structure analyze_structure(u64 p, u64 m) {
    Structure s{};

    u64 x = m;

    while (x > 0 && x % p == p - 1) {
        ++s.a;
        x /= p;
    }

    if (x == 0) {
        s.exponent = std::max(1, s.a);
        s.first_start = 1;
        s.modulus = 1;

        for (int i = 0; i < s.exponent; ++i) {
            s.modulus *= p;
        }

        return s;
    }

    const u64 b = x % p;

    u64 pa = 1;

    for (int i = 0; i < s.a; ++i) {
        pa *= p;
    }

    x /= p;

    while (x > 0 && x % p == 0) {
        ++s.z;
        x /= p;
    }

    s.exponent = s.a + 1 + s.z;
    s.first_start = (b + 1) * pa;

    s.modulus = 1;

    for (int i = 0; i < s.exponent; ++i) {
        s.modulus *= p;
    }

    return s;
}

std::vector<u64> collect_starts(
    u64 p,
    u64 m
) {
    std::vector<u64> starts;

    bool previous_hit = false;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = lucas_hit(p, m, t);

        if (hit && !previous_hit) {
            starts.push_back(t);
        }

        previous_hit = hit;
    }

    return starts;
}

std::vector<u64> index_sequence(
    const std::vector<u64>& starts,
    u64 first,
    u64 modulus
) {
    std::vector<u64> indices;

    for (u64 s : starts) {
        if (s < first ||
            (s - first) % modulus != 0) {
            continue;
        }

        indices.push_back(
            (s - first) / modulus
        );
    }

    return indices;
}

std::string vector_string(
    const std::vector<u64>& values,
    std::size_t limit = 30
) {
    std::string s = "[";

    const std::size_t n =
        std::min(values.size(), limit);

    for (std::size_t i = 0; i < n; ++i) {
        if (i != 0) s += ",";
        s += std::to_string(values[i]);
    }

    if (values.size() > limit) {
        s += ",...";
    }

    s += "]";

    return s;
}

void run_exhaustive() {
    std::cout
        << "\nPHASE 1: INDEX-LATTICE VALIDATION\n";

    u64 total = 0;
    u64 residue_pass = 0;
    u64 first_pass = 0;
    u64 gcd_index_pass = 0;
    u64 monotone_pass = 0;

    for (u64 p = 2; p <= 37; ++p) {
        if (!is_prime(p)) continue;

        for (u64 m = 1; m <= 1500; ++m) {
            const Structure s =
                analyze_structure(p, m);

            const auto starts =
                collect_starts(p, m);

            if (starts.size() < 2) continue;

            ++total;

            bool residues_ok = true;

            for (u64 start : starts) {
                if (start < s.first_start ||
                    (start - s.first_start) %
                        s.modulus != 0) {
                    residues_ok = false;
                    break;
                }
            }

            if (residues_ok) {
                ++residue_pass;
            }

            if (starts.front() ==
                s.first_start) {
                ++first_pass;
            }

            const auto indices =
                index_sequence(
                    starts,
                    s.first_start,
                    s.modulus
                );

            u64 g = 0;

            for (std::size_t i = 1;
                 i < indices.size();
                 ++i) {
                g = std::gcd(
                    g,
                    indices[i] -
                    indices[i - 1]
                );
            }

            if (g == 1) {
                ++gcd_index_pass;
            }

            bool strictly_increasing = true;

            for (std::size_t i = 1;
                 i < indices.size();
                 ++i) {

                if (indices[i] <=
                    indices[i - 1]) {

                    strictly_increasing = false;
                    break;
                }
            }

            if (strictly_increasing) {
                ++monotone_pass;
            }

            if (!residues_ok ||
                starts.front() !=
                    s.first_start ||
                g != 1 ||
                !strictly_increasing) {

                if (total <= 20) {
                    std::cout
                        << "FAILURE "
                        << "p=" << p
                        << " m=" << m
                        << " starts="
                        << vector_string(starts)
                        << " indices="
                        << vector_string(indices)
                        << "\n";
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "multi_start_cases="
        << total << "\n"
        << "residue_pass="
        << residue_pass
        << "/" << total << "\n"
        << "first_start_pass="
        << first_pass
        << "/" << total << "\n"
        << "index_gcd_1="
        << gcd_index_pass
        << "/" << total << "\n"
        << "strictly_increasing_indices="
        << monotone_pass
        << "/" << total << "\n";
}

void run_examples() {
    std::cout
        << "\nPHASE 2: INDEX EXAMPLES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 12},
        {2, 20},
        {2, 24},
        {2, 28},
        {2, 40},
        {2, 48},
        {3, 41},
        {3, 59},
        {3, 125},
        {3, 377},
        {5, 59},
        {5, 94}
    };

    for (const auto& c : cases) {
        const Structure s =
            analyze_structure(
                c.p,
                c.m
            );

        const auto starts =
            collect_starts(
                c.p,
                c.m
            );

        const auto indices =
            index_sequence(
                starts,
                s.first_start,
                s.modulus
            );

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " e=" << s.exponent
            << " s0=" << s.first_start
            << " modulus=" << s.modulus
            << "\n";

        std::cout
            << "starts="
            << vector_string(starts)
            << "\n";

        std::cout
            << "indices="
            << vector_string(indices)
            << "\n";
    }
}

void run_large_direct() {
    std::cout
        << "\nPHASE 3: LARGE DIRECT BOUNDARY TEST\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37, 43,
        53, 67, 101, 127,
        211, 431, 1009
    };

    u64 state = 0x1692026ULL;

    constexpr int CASES = 10000;

    u64 tested = 0;
    u64 failures = 0;

    for (int i = 0; i < CASES; ++i) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 p =
            primes[state % primes.size()];

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 m =
            1 + state % 1000000000000ULL;

        const Structure s =
            analyze_structure(p, m);

        /*
            We can test the first predicted start
            and one complete period later without
            enumerating all starts.
        */
        const bool first =
            lucas_hit(
                p,
                m,
                s.first_start
            );

        const bool before =
            s.first_start > 0 &&
            !lucas_hit(
                p,
                m,
                s.first_start - 1
            );

        if (!first || !before) {
            ++failures;

            if (failures <= 10) {
                std::cout
                    << "FAIL_FIRST "
                    << "p=" << p
                    << " m=" << m
                    << "\n";
            }

            continue;
        }

        ++tested;
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "tested="
        << tested
        << "\n"
        << "failures="
        << failures
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 169;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive();
    run_examples();
    run_large_direct();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
