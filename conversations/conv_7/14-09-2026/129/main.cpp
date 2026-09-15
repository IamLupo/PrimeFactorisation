#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Witness {
    int m = 1;
    u64 k = 1;
    u64 t = 0;
    int sign = +1;
};

struct Excursion {
    u64 prime = 0;
    std::vector<Witness> states;
};

static std::vector<u64> generate_primes(int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit) + 1, true);

    is_prime[0] = false;
    is_prime[1] = false;

    for (int p = 2; static_cast<long long>(p) * p <= limit; ++p) {
        if (!is_prime[p]) {
            continue;
        }

        for (int x = p * p; x <= limit; x += p) {
            is_prime[x] = false;
        }
    }

    std::vector<u64> primes;
    for (int x = 2; x <= limit; ++x) {
        if (is_prime[x]) {
            primes.push_back(static_cast<u64>(x));
        }
    }

    return primes;
}

static bool divides_stream(u64 r, int m, u64 k, int sign) {
    const i128 value =
        static_cast<i128>(m) * static_cast<i128>(r) + sign;

    return value > 0 &&
           (value % static_cast<i128>(k) == 0);
}

static Witness make_witness(u64 r, int m, u64 k, int sign) {
    Witness w;

    w.m = m;
    w.k = k;
    w.sign = sign;

    const i128 value =
        static_cast<i128>(m) * static_cast<i128>(r) + sign;

    w.t = static_cast<u64>(
        value / static_cast<i128>(k)
    );

    return w;
}

static bool better_normalized(
    const Witness& a,
    const Witness& b
) {
    const i128 lhs =
        static_cast<i128>(a.k) * static_cast<i128>(b.m);

    const i128 rhs =
        static_cast<i128>(b.k) * static_cast<i128>(a.m);

    if (lhs != rhs) {
        return lhs > rhs;
    }

    // Deterministic tie-breaking.
    if (a.m != b.m) {
        return a.m < b.m;
    }

    if (a.k != b.k) {
        return a.k > b.k;
    }

    return a.sign > b.sign;
}

static Witness global_winner(
    const std::vector<Witness>& best_by_m
) {
    Witness best = best_by_m[1];

    for (std::size_t m = 2; m < best_by_m.size(); ++m) {
        if (better_normalized(best_by_m[m], best)) {
            best = best_by_m[m];
        }
    }

    return best;
}

static std::string projected_pattern(
    const std::vector<Witness>& states
) {
    std::string out;

    for (std::size_t i = 0; i < states.size(); ++i) {
        if (i != 0) {
            out += ',';
        }

        out += std::to_string(states[i].m);
    }

    return out;
}

static std::string signed_pattern(
    const std::vector<Witness>& states
) {
    std::string out;

    for (std::size_t i = 0; i < states.size(); ++i) {
        if (i != 0) {
            out += ',';
        }

        out += std::to_string(states[i].m);
        out += (states[i].sign > 0 ? '+' : '-');
    }

    return out;
}

static bool projected_palindrome(
    const std::vector<Witness>& states
) {
    for (
        std::size_t i = 0,
        j = states.size() - 1;
        i < j;
        ++i, --j
    ) {
        if (states[i].m != states[j].m) {
            return false;
        }
    }

    return true;
}

static std::string mirror_mask(
    const std::vector<Witness>& states
) {
    const std::size_t pairs = states.size() / 2;

    std::string mask;
    mask.reserve(pairs);

    for (std::size_t i = 0; i < pairs; ++i) {
        const std::size_t j =
            states.size() - 1 - i;

        mask.push_back(
            states[i].sign == states[j].sign
                ? 'S'
                : 'F'
        );
    }

    return mask;
}

static int sign_changes(
    const std::vector<Witness>& states
) {
    int changes = 0;

    for (std::size_t i = 1; i < states.size(); ++i) {
        if (states[i - 1].sign != states[i].sign) {
            ++changes;
        }
    }

    return changes;
}

static u64 stream_value(
    u64 r,
    const Witness& w
) {
    const i128 value =
        static_cast<i128>(w.m) * static_cast<i128>(r) +
        w.sign;

    if (
        value <= 0 ||
        value >
            static_cast<i128>(
                std::numeric_limits<u64>::max()
            )
    ) {
        return 0;
    }

    return static_cast<u64>(value);
}

static void print_witness(
    u64 r,
    const std::string& label,
    const Witness& w
) {
    std::cout
        << label
        << "=(m=" << w.m
        << ",k=" << w.k
        << ",t=" << w.t
        << ",sign=" << (w.sign > 0 ? "+1" : "-1")
        << ",value=" << stream_value(r, w)
        << ")";
}

static bool interior_sign_mirror_ok(
    const std::vector<Witness>& states
) {
    const std::size_t n = states.size();

    if (n <= 2) {
        return true;
    }

    // Only test mirrored pairs where both sides are interior
    // states, i.e. m > 1.
    for (std::size_t i = 1; i < n / 2 + (n % 2); ++i) {
        const std::size_t j = n - 1 - i;

        if (i >= j) {
            break;
        }

        if (
            states[i].m == 1 ||
            states[j].m == 1
        ) {
            continue;
        }

        if (states[i].m != states[j].m) {
            return false;
        }

        if (states[i].sign != states[j].sign) {
            return false;
        }
    }

    return true;
}

static void print_excursion(
    const Excursion& ex
) {
    std::cout
        << "PRIME=" << ex.prime
        << " LENGTH=" << ex.states.size()
        << " PROJECTED="
        << projected_pattern(ex.states)
        << " SIGNED="
        << signed_pattern(ex.states)
        << " MIRROR_MASK="
        << mirror_mask(ex.states)
        << " SIGN_CHANGES="
        << sign_changes(ex.states)
        << "\n";

    for (std::size_t i = 0; i < ex.states.size(); ++i) {
        std::cout << "  STATE[" << i << "] ";

        print_witness(
            ex.prime,
            "W",
            ex.states[i]
        );

        std::cout << "\n";
    }
}

static void process_excursion(
    const Excursion& ex,
    u64& interior_symmetric,
    u64& interior_non_symmetric,
    u64& endpoint_only_failures,
    u64& interior_failures,
    u64& full_signed_mirror,
    u64& interior_pair_tests,
    u64& interior_pair_mismatches,
    std::map<std::string, u64>& mask_counts,
    std::map<int, u64>& sign_change_hist,
    Excursion& first_endpoint_only_failure,
    bool& have_endpoint_only_failure,
    Excursion& first_interior_counterexample,
    bool& have_interior_counterexample
) {
    const std::vector<Witness>& s = ex.states;
    const std::size_t n = s.size();

    bool any_sign_failure = false;
    bool any_interior_failure = false;

    for (std::size_t i = 0; i < n / 2; ++i) {
        const std::size_t j =
            n - 1 - i;

        if (s[i].sign != s[j].sign) {
            any_sign_failure = true;

            if (
                s[i].m > 1 &&
                s[j].m > 1
            ) {
                any_interior_failure = true;
            }
        }
    }

    if (interior_sign_mirror_ok(s)) {
        ++interior_symmetric;
    } else {
        ++interior_non_symmetric;
    }

    if (
        any_sign_failure &&
        !any_interior_failure
    ) {
        ++endpoint_only_failures;

        if (!have_endpoint_only_failure) {
            first_endpoint_only_failure = ex;
            have_endpoint_only_failure = true;
        }
    }

    if (any_interior_failure) {
        ++interior_failures;

        if (!have_interior_counterexample) {
            first_interior_counterexample = ex;
            have_interior_counterexample = true;
        }
    }

    bool signed_mirror = true;

    for (std::size_t i = 0; i < n / 2; ++i) {
        const std::size_t j =
            n - 1 - i;

        if (
            s[i].m > 1 ||
            s[j].m > 1
        ) {
            ++interior_pair_tests;

            if (s[i].sign != s[j].sign) {
                ++interior_pair_mismatches;
            }
        }

        if (
            s[i].m != s[j].m ||
            s[i].sign != s[j].sign
        ) {
            signed_mirror = false;
        }
    }

    if (signed_mirror) {
        ++full_signed_mirror;
    }

    mask_counts[mirror_mask(s)]++;
    sign_change_hist[sign_changes(s)]++;
}

int main() {
    constexpr int EXPERIMENT = 421;
    constexpr int PRIME_LIMIT = 5000;
    constexpr int K_LIMIT = 2000;
    constexpr int M_MAX = 32;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << "\n";

    std::cout
        << "M_MAX="
        << M_MAX
        << "\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    u64 total_winner_points = 0;
    u64 total_transitions = 0;
    u64 total_excursions = 0;

    u64 projected_palindromes = 0;
    u64 projected_nonpalindromes = 0;

    u64 primes_with_excursion = 0;
    u64 primes_without_excursion = 0;

    u64 interior_symmetric = 0;
    u64 interior_non_symmetric = 0;

    u64 endpoint_only_failures = 0;
    u64 interior_failures = 0;

    u64 full_signed_mirror = 0;

    u64 interior_pair_tests = 0;
    u64 interior_pair_mismatches = 0;

    std::map<std::string, u64> mask_counts;
    std::map<int, u64> sign_change_hist;

    Excursion first_endpoint_only_failure;
    bool have_endpoint_only_failure = false;

    Excursion first_interior_counterexample;
    bool have_interior_counterexample = false;

    for (u64 r : primes) {
        std::vector<Witness> best_by_m(
            static_cast<std::size_t>(M_MAX) + 1
        );

        // At K=1 every m has k=1.
        for (int m = 1; m <= M_MAX; ++m) {
            best_by_m[
                static_cast<std::size_t>(m)
            ] = make_witness(
                r,
                m,
                1,
                +1
            );
        }

        std::vector<Witness> transition_states;
        transition_states.reserve(128);

        Witness previous =
            global_winner(best_by_m);

        transition_states.push_back(previous);
        ++total_winner_points;

        bool had_excursion = false;

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            for (int m = 1; m <= M_MAX; ++m) {
                const std::size_t mi =
                    static_cast<std::size_t>(m);

                if (
                    divides_stream(
                        r,
                        m,
                        K,
                        +1
                    )
                ) {
                    best_by_m[mi] =
                        make_witness(
                            r,
                            m,
                            K,
                            +1
                        );
                } else if (
                    divides_stream(
                        r,
                        m,
                        K,
                        -1
                    )
                ) {
                    best_by_m[mi] =
                        make_witness(
                            r,
                            m,
                            K,
                            -1
                        );
                }
            }

            const Witness current =
                global_winner(best_by_m);

            ++total_winner_points;

            if (current.m != previous.m) {
                ++total_transitions;
                transition_states.push_back(current);
            }

            previous = current;
        }

        std::size_t start = 0;

        while (start < transition_states.size()) {
            while (
                start < transition_states.size() &&
                transition_states[start].m != 1
            ) {
                ++start;
            }

            if (start >= transition_states.size()) {
                break;
            }

            std::size_t end = start + 1;

            while (
                end < transition_states.size() &&
                transition_states[end].m != 1
            ) {
                ++end;
            }

            if (end >= transition_states.size()) {
                break;
            }

            if (end > start + 1) {
                had_excursion = true;

                Excursion ex;
                ex.prime = r;

                ex.states.assign(
                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(start),

                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(end) +
                        1
                );

                ++total_excursions;

                if (projected_palindrome(ex.states)) {
                    ++projected_palindromes;
                } else {
                    ++projected_nonpalindromes;
                }

                process_excursion(
                    ex,
                    interior_symmetric,
                    interior_non_symmetric,
                    endpoint_only_failures,
                    interior_failures,
                    full_signed_mirror,
                    interior_pair_tests,
                    interior_pair_mismatches,
                    mask_counts,
                    sign_change_hist,
                    first_endpoint_only_failure,
                    have_endpoint_only_failure,
                    first_interior_counterexample,
                    have_interior_counterexample
                );
            }

            start = end;
        }

        if (had_excursion) {
            ++primes_with_excursion;
        } else {
            ++primes_without_excursion;
        }
    }

    std::cout
        << "\nTOTAL_WINNER_POINTS="
        << total_winner_points
        << "\n";

    std::cout
        << "TOTAL_TRANSITIONS="
        << total_transitions
        << "\n";

    std::cout
        << "TOTAL_EXCURSIONS="
        << total_excursions
        << "\n";

    std::cout
        << "PROJECTED_PALINDROMES="
        << projected_palindromes
        << "\n";

    std::cout
        << "PROJECTED_NONPALINDROMES="
        << projected_nonpalindromes
        << "\n";

    std::cout
        << "PRIMES_WITH_EXCURSION="
        << primes_with_excursion
        << "\n";

    std::cout
        << "PRIMES_WITHOUT_EXCURSION="
        << primes_without_excursion
        << "\n";

    std::cout
        << "INTERIOR_SIGN_MIRROR_SYMMETRIC="
        << interior_symmetric
        << "\n";

    std::cout
        << "INTERIOR_SIGN_MIRROR_NON_SYMMETRIC="
        << interior_non_symmetric
        << "\n";

    std::cout
        << "ENDPOINT_ONLY_SIGN_FAILURES="
        << endpoint_only_failures
        << "\n";

    std::cout
        << "INTERIOR_SIGN_FAILURES="
        << interior_failures
        << "\n";

    std::cout
        << "FULL_SIGNED_MIRROR_EXCURSIONS="
        << full_signed_mirror
        << "\n";

    std::cout
        << "INTERIOR_MIRROR_PAIR_TESTS="
        << interior_pair_tests
        << "\n";

    std::cout
        << "INTERIOR_MIRROR_SIGN_MISMATCHES="
        << interior_pair_mismatches
        << "\n";

    std::cout
        << "\nSIGN_MIRROR_MASK_COUNTS\n";

    for (const auto& entry : mask_counts) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nSIGN_CHANGE_HISTOGRAM\n";

    for (const auto& entry : sign_change_hist) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_endpoint_only_failure) {
        std::cout
            << "\nFIRST_ENDPOINT_ONLY_SIGN_FAILURE\n";

        print_excursion(
            first_endpoint_only_failure
        );
    } else {
        std::cout
            << "\nFIRST_ENDPOINT_ONLY_SIGN_FAILURE NONE\n";
    }

    if (have_interior_counterexample) {
        std::cout
            << "\nFIRST_INTERIOR_SIGN_COUNTEREXAMPLE\n";

        print_excursion(
            first_interior_counterexample
        );
    } else {
        std::cout
            << "\nFIRST_INTERIOR_SIGN_COUNTEREXAMPLE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
