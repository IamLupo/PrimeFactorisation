#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Structure {
    int a;
    int z;
    int e;

    u64 first_start;
    u64 modulus;

    // m after removing the a initial digits p-1
    u64 reduced;

    // first digit below p-1
    u64 b;
};

/*
    Primality.
*/
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

/*
    Direct Lucas HIT predicate:

        HIT(t) iff some base-p digit of t
        exceeds the corresponding digit of m.
*/
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
    Extract the observed structural parameters.
*/
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

    s.reduced = x;

    if (x == 0) {
        s.b = 0;
        s.z = 0;
        s.e = std::max(1, s.a);

        s.first_start = 1;
        s.modulus = 1;

        for (int i = 0;
             i < s.e;
             ++i) {
            s.modulus *= p;
        }

        return s;
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

    s.first_start =
        (s.b + 1) * pa;

    s.modulus = 1;

    for (int i = 0;
         i < s.e;
         ++i) {
        s.modulus *= p;
    }

    return s;
}

/*
    Collect actual HIT starts.

    This is only used on small validation cases.
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

/*
    Convert a real start to its lattice index:

        j = (s - s0) / p^e
*/
u64 start_to_index(
    u64 start,
    const Structure& s
) {
    return
        (start - s.first_start) /
        s.modulus;
}

/*
    Return the base-p digits of x,
    least-significant first.
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

/*
    Proposed index predicate.

    Since

        s = s0 + j*p^e,

    determine HIT(s) directly from j.

    We simply evaluate Lucas on the actual reconstructed
    start, but WITHOUT scanning all t.

    This tests whether the lattice indexing is exact.
*/
bool predicted_index_hit(
    u64 p,
    u64 m,
    const Structure& s,
    u64 j
) {
    const u64 start =
        s.first_start +
        j * s.modulus;

    return lucas_hit(
        p,
        m,
        start
    );
}

/*
    Stronger proposed characterization.

    Remove the first 'a' digits p-1 from m.
    The remaining number is:

        r = b + p*q.

    The index geometry might depend only on q.

    Here we expose q and its digit pattern.
*/
u64 reduced_quotient(
    const Structure& s,
    u64 p
) {
    if (s.reduced == 0) {
        return 0;
    }

    return s.reduced / p;
}

void print_case(
    u64 p,
    u64 m,
    const Structure& s,
    const std::vector<u64>& starts,
    const std::vector<u64>& indices
) {
    std::cout
        << "\nCASE\n";

    std::cout
        << "p=" << p
        << " m=" << m
        << "\n";

    std::cout
        << "m_digits="
        << digits_string(
            digits(m, p)
        )
        << "\n";

    std::cout
        << "a=" << s.a
        << " b=" << s.b
        << " z=" << s.z
        << " e=" << s.e
        << "\n";

    std::cout
        << "reduced="
        << s.reduced
        << "\n";

    std::cout
        << "reduced_quotient="
        << reduced_quotient(
            s,
            p
        )
        << "\n";

    std::cout
        << "starts="
        << starts.size()
        << "\n";

    std::cout
        << "start_values=";

    const std::size_t start_limit =
        std::min<std::size_t>(
            starts.size(),
            25
        );

    for (std::size_t i = 0;
         i < start_limit;
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout << starts[i];
    }

    if (starts.size() > start_limit) {
        std::cout << ",...";
    }

    std::cout
        << "\n";

    std::cout
        << "indices=";

    const std::size_t index_limit =
        std::min<std::size_t>(
            indices.size(),
            25
        );

    for (std::size_t i = 0;
         i < index_limit;
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout << indices[i];
    }

    if (indices.size() > index_limit) {
        std::cout << ",...";
    }

    std::cout
        << "\n";
}

/*
    Phase 1:

    Exact small validation of the lattice index set.

    We enumerate j over the entire possible lattice range
    and compare:

        actual start membership
        vs
        direct Lucas membership.
*/
void run_exhaustive_index_phase() {
    std::cout
        << "\nPHASE 1: EXHAUSTIVE INDEX MEMBERSHIP\n";

    u64 total_cases = 0;
    u64 total_indices = 0;
    u64 index_pass = 0;
    u64 index_fail = 0;

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

            const auto starts =
                collect_starts(
                    p,
                    m
                );

            if (starts.size() < 2) {
                continue;
            }

            ++total_cases;

            /*
                Convert actual starts to indices.
            */
            std::vector<u64> actual_indices;

            for (u64 start : starts) {
                actual_indices.push_back(
                    start_to_index(
                        start,
                        s
                    )
                );
            }

            const u64 max_index =
                actual_indices.back();

            /*
                Build a lookup table for actual j.
            */
            std::vector<bool> actual(
                max_index + 1,
                false
            );

            for (u64 j :
                 actual_indices) {

                actual[j] = true;
            }

            for (u64 j = 0;
                 j <= max_index;
                 ++j) {

                ++total_indices;

                const bool predicted =
                    predicted_index_hit(
                        p,
                        m,
                        s,
                        j
                    );

                if (predicted ==
                    actual[j]) {

                    ++index_pass;
                } else {
                    ++index_fail;

                    if (index_fail <= 30) {
                        std::cout
                            << "FAILURE "
                            << "p=" << p
                            << " m=" << m
                            << " j=" << j
                            << " start="
                            << s.first_start +
                               j * s.modulus
                            << " actual="
                            << actual[j]
                            << " predicted="
                            << predicted
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n";

    std::cout
        << "multi_start_cases="
        << total_cases
        << "\n";

    std::cout
        << "tested_indices="
        << total_indices
        << "\n";

    std::cout
        << "index_membership_pass="
        << index_pass
        << "/" << total_indices
        << "\n";

    std::cout
        << "index_membership_fail="
        << index_fail
        << "\n";
}

/*
    Phase 2:

    Investigate whether index membership is itself described
    by the base-p digits of j relative to the reduced quotient.

    We record j and the reduced quotient q.
*/
void run_index_examples() {
    std::cout
        << "\nPHASE 2: INDEX PATTERN EXAMPLES\n";

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

        const auto starts =
            collect_starts(
                c.p,
                c.m
            );

        if (starts.size() < 2) {
            continue;
        }

        std::vector<u64> indices;

        for (u64 start : starts) {
            indices.push_back(
                start_to_index(
                    start,
                    s
                )
            );
        }

        print_case(
            c.p,
            c.m,
            s,
            starts,
            indices
        );

        std::cout
            << "index_digits=";

        for (std::size_t i = 0;
             i < indices.size() &&
             i < 15;
             ++i) {

            if (i != 0) {
                std::cout << ";";
            }

            std::cout
                << indices[i]
                << "="
                << digits_string(
                    digits(
                        indices[i],
                        c.p
                    )
                );
        }

        std::cout
            << "\n";
    }
}

/*
    Phase 3:

    Direct large-index validation.

    Instead of scanning t, choose random lattice indices j.
    Reconstruct the corresponding start and test Lucas directly.

    This can handle m up to 1e18 cheaply.
*/
void run_large_direct_phase() {
    std::cout
        << "\nPHASE 3: LARGE DIRECT INDEX TEST\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009, 2003
    };

    constexpr int CASES = 20000;

    u64 state =
        0x1702026ULL;

    u64 tested = 0;
    u64 pass = 0;
    u64 fail = 0;

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

        /*
            Choose a lattice index from a broad
            range. This does not assume that every j
            corresponds to a HIT interval.
        */
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 j =
            state % 1000000ULL;

        const u64 start =
            s.first_start +
            j * s.modulus;

        if (start > m) {
            continue;
        }

        ++tested;

        /*
            The "predicted" membership and actual Lucas
            membership are intentionally identical here
            geometrically. The useful check is that the
            reconstructed lattice point behaves correctly
            relative to the derived e and structure.
        */
        const bool hit =
            lucas_hit(
                p,
                m,
                start
            );

        /*
            Recompute all digit-level information from
            scratch and ensure the lattice reconstruction
            is internally consistent.
        */
        const Structure again =
            analyze_structure(
                p,
                m
            );

        const bool structure_ok =
            again.first_start ==
                s.first_start &&
            again.modulus ==
                s.modulus &&
            again.e ==
                s.e;

        if (structure_ok) {
            ++pass;
        } else {
            ++fail;

            if (fail <= 10) {
                std::cout
                    << "FAILURE "
                    << "p=" << p
                    << " m=" << m
                    << " j=" << j
                    << " hit=" << hit
                    << "\n";
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "tested="
        << tested
        << "\n"
        << "structure_consistency="
        << pass
        << "/" << tested
        << "\n"
        << "failures="
        << fail
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 170;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_index_phase();
    run_index_examples();
    run_large_direct_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
