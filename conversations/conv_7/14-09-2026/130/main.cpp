#include <algorithm>
#include <cstdint>
#include <iostream>
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
    std::vector<bool> is_prime(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    is_prime[0] = false;
    is_prime[1] = false;

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!is_prime[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            is_prime[x] = false;
        }
    }

    std::vector<u64> primes;
    primes.reserve(static_cast<std::size_t>(limit));

    for (int x = 2; x <= limit; ++x) {
        if (is_prime[x]) {
            primes.push_back(static_cast<u64>(x));
        }
    }

    return primes;
}

static bool divides_stream(
    u64 r,
    int m,
    u64 k,
    int sign
) {
    const i128 value =
        static_cast<i128>(m) *
        static_cast<i128>(r) +
        static_cast<i128>(sign);

    if (value <= 0) {
        return false;
    }

    return value %
        static_cast<i128>(k) == 0;
}

static Witness make_witness(
    u64 r,
    int m,
    u64 k,
    int sign
) {
    Witness w;

    w.m = m;
    w.k = k;
    w.sign = sign;

    const i128 value =
        static_cast<i128>(m) *
        static_cast<i128>(r) +
        static_cast<i128>(sign);

    w.t = static_cast<u64>(
        value /
        static_cast<i128>(k)
    );

    return w;
}

static bool better_normalized(
    const Witness& a,
    const Witness& b
) {
    const i128 lhs =
        static_cast<i128>(a.k) *
        static_cast<i128>(b.m);

    const i128 rhs =
        static_cast<i128>(b.k) *
        static_cast<i128>(a.m);

    if (lhs != rhs) {
        return lhs > rhs;
    }

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

    for (
        std::size_t m = 2;
        m < best_by_m.size();
        ++m
    ) {
        if (
            better_normalized(
                best_by_m[m],
                best
            )
        ) {
            best = best_by_m[m];
        }
    }

    return best;
}

static bool projected_palindrome(
    const std::vector<Witness>& states
) {
    if (states.empty()) {
        return true;
    }

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

static bool exact_factor_reversal(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.sign == b.sign &&
        a.k == b.t &&
        a.t == b.k;
}

static bool unordered_factor_match(
    const Witness& a,
    const Witness& b
) {
    if (a.m != b.m) {
        return false;
    }

    return
        (a.k == b.k && a.t == b.t) ||
        (a.k == b.t && a.t == b.k);
}

static bool same_arithmetic_value(
    u64 r,
    const Witness& a,
    const Witness& b
) {
    const i128 va =
        static_cast<i128>(a.m) *
        static_cast<i128>(r) +
        static_cast<i128>(a.sign);

    const i128 vb =
        static_cast<i128>(b.m) *
        static_cast<i128>(r) +
        static_cast<i128>(b.sign);

    return va == vb;
}

static u64 arithmetic_difference(
    u64 r,
    const Witness& a,
    const Witness& b
) {
    const i128 va =
        static_cast<i128>(a.m) *
        static_cast<i128>(r) +
        static_cast<i128>(a.sign);

    const i128 vb =
        static_cast<i128>(b.m) *
        static_cast<i128>(r) +
        static_cast<i128>(b.sign);

    const i128 d =
        va >= vb ? va - vb : vb - va;

    return static_cast<u64>(d);
}

static u64 abs_diff(
    u64 a,
    u64 b
) {
    return a >= b ? a - b : b - a;
}

static void print_witness(
    u64 r,
    const std::string& label,
    const Witness& w
) {
    const i128 value =
        static_cast<i128>(w.m) *
        static_cast<i128>(r) +
        static_cast<i128>(w.sign);

    std::cout
        << label
        << "=(m=" << w.m
        << ",k=" << w.k
        << ",t=" << w.t
        << ",sign="
        << (w.sign > 0 ? "+1" : "-1")
        << ",value=";

    // m <= 32 and r <= 5000, so this always fits in u64.
    std::cout << static_cast<u64>(value);

    std::cout << ")";
}

static void print_excursion(
    const Excursion& ex
) {
    std::cout
        << "PRIME="
        << ex.prime
        << " LENGTH="
        << ex.states.size()
        << "\n";

    for (
        std::size_t i = 0;
        i < ex.states.size();
        ++i
    ) {
        std::cout
            << "  STATE["
            << i
            << "] ";

        print_witness(
            ex.prime,
            "W",
            ex.states[i]
        );

        std::cout << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 422;
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

    u64 total_excursions = 0;
    u64 total_mirrored_pairs = 0;

    u64 projected_palindromes = 0;

    u64 exact_factor_reversals = 0;
    u64 unordered_factor_matches = 0;

    u64 same_stream_pairs = 0;
    u64 different_stream_pairs = 0;

    u64 exact_reversal_same_stream = 0;
    u64 exact_reversal_different_stream = 0;

    u64 only_k_matches = 0;
    u64 only_t_matches = 0;
    u64 neither_factor_matches = 0;

    u64 arithmetic_equal_pairs = 0;
    u64 arithmetic_difference_two = 0;
    u64 arithmetic_other_difference = 0;

    u64 factor_reversal_failures = 0;
    u64 same_stream_reversal_failures = 0;

    std::map<u64, u64> factor_error_histogram;

    Excursion first_factor_failure;
    bool have_first_factor_failure = false;

    Excursion first_same_stream_failure;
    bool have_first_same_stream_failure = false;

    for (u64 r : primes) {
        std::vector<Witness> best_by_m(
            static_cast<std::size_t>(M_MAX) + 1
        );

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

            if (current.m != previous.m) {
                transition_states.push_back(
                    current
                );
            }

            previous = current;
        }

        std::size_t start = 0;

        while (
            start < transition_states.size()
        ) {
            while (
                start < transition_states.size() &&
                transition_states[start].m != 1
            ) {
                ++start;
            }

            if (
                start >=
                transition_states.size()
            ) {
                break;
            }

            std::size_t end = start + 1;

            while (
                end < transition_states.size() &&
                transition_states[end].m != 1
            ) {
                ++end;
            }

            if (
                end >=
                transition_states.size()
            ) {
                break;
            }

            if (end > start + 1) {
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

                if (
                    projected_palindrome(
                        ex.states
                    )
                ) {
                    ++projected_palindromes;
                }

                const std::size_t n =
                    ex.states.size();

                for (
                    std::size_t i = 0;
                    i < n / 2;
                    ++i
                ) {
                    const std::size_t j =
                        n - 1 - i;

                    const Witness& a =
                        ex.states[i];

                    const Witness& b =
                        ex.states[j];

                    ++total_mirrored_pairs;

                    const bool same_stream =
                        a.m == b.m &&
                        a.sign == b.sign;

                    if (same_stream) {
                        ++same_stream_pairs;
                    } else {
                        ++different_stream_pairs;
                    }

                    if (
                        same_arithmetic_value(
                            r,
                            a,
                            b
                        )
                    ) {
                        ++arithmetic_equal_pairs;
                    } else {
                        const u64 d =
                            arithmetic_difference(
                                r,
                                a,
                                b
                            );

                        if (d == 2) {
                            ++arithmetic_difference_two;
                        } else {
                            ++arithmetic_other_difference;
                        }
                    }

                    if (
                        exact_factor_reversal(
                            a,
                            b
                        )
                    ) {
                        ++exact_factor_reversals;

                        if (same_stream) {
                            ++exact_reversal_same_stream;
                        } else {
                            ++exact_reversal_different_stream;
                        }
                    } else {
                        ++factor_reversal_failures;

                        const u64 e1 =
                            abs_diff(a.k, b.t);

                        const u64 e2 =
                            abs_diff(a.t, b.k);

                        ++factor_error_histogram[
                            e1 + e2
                        ];

                        if (
                            !have_first_factor_failure
                        ) {
                            first_factor_failure = ex;
                            have_first_factor_failure = true;
                        }

                        if (same_stream) {
                            ++same_stream_reversal_failures;

                            if (
                                !have_first_same_stream_failure
                            ) {
                                first_same_stream_failure = ex;
                                have_first_same_stream_failure = true;
                            }
                        }
                    }

                    if (
                        unordered_factor_match(
                            a,
                            b
                        )
                    ) {
                        ++unordered_factor_matches;
                    }

                    const bool k_match =
                        a.k == b.t;

                    const bool t_match =
                        a.t == b.k;

                    if (k_match && !t_match) {
                        ++only_k_matches;
                    }

                    if (t_match && !k_match) {
                        ++only_t_matches;
                    }

                    if (!k_match && !t_match) {
                        ++neither_factor_matches;
                    }
                }
            }

            start = end;
        }
    }

    std::cout
        << "\nTOTAL_EXCURSIONS="
        << total_excursions
        << "\n";

    std::cout
        << "PROJECTED_PALINDROMES="
        << projected_palindromes
        << "\n";

    std::cout
        << "TOTAL_MIRRORED_PAIRS="
        << total_mirrored_pairs
        << "\n";

    std::cout
        << "SAME_ARITHMETIC_STREAM_PAIRS="
        << same_stream_pairs
        << "\n";

    std::cout
        << "DIFFERENT_ARITHMETIC_STREAM_PAIRS="
        << different_stream_pairs
        << "\n";

    std::cout
        << "EXACT_FACTOR_REVERSALS="
        << exact_factor_reversals
        << "\n";

    std::cout
        << "UNORDERED_FACTOR_MATCHES="
        << unordered_factor_matches
        << "\n";

    std::cout
        << "EXACT_REVERSALS_SAME_STREAM="
        << exact_reversal_same_stream
        << "\n";

    std::cout
        << "EXACT_REVERSALS_DIFFERENT_STREAM="
        << exact_reversal_different_stream
        << "\n";

    std::cout
        << "ONLY_K_MATCHES="
        << only_k_matches
        << "\n";

    std::cout
        << "ONLY_T_MATCHES="
        << only_t_matches
        << "\n";

    std::cout
        << "NEITHER_FACTOR_MATCHES="
        << neither_factor_matches
        << "\n";

    std::cout
        << "ARITHMETIC_EQUAL_PAIRS="
        << arithmetic_equal_pairs
        << "\n";

    std::cout
        << "ARITHMETIC_DIFFERENCE_TWO="
        << arithmetic_difference_two
        << "\n";

    std::cout
        << "ARITHMETIC_OTHER_DIFFERENCE="
        << arithmetic_other_difference
        << "\n";

    std::cout
        << "FACTOR_REVERSAL_FAILURES="
        << factor_reversal_failures
        << "\n";

    std::cout
        << "SAME_STREAM_REVERSAL_FAILURES="
        << same_stream_reversal_failures
        << "\n";

    std::cout
        << "\nFACTOR_REVERSAL_ERROR_HISTOGRAM\n";

    for (
        const auto& entry :
        factor_error_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_factor_failure) {
        std::cout
            << "\nFIRST_FACTOR_REVERSAL_FAILURE\n";

        print_excursion(
            first_factor_failure
        );
    } else {
        std::cout
            << "\nFIRST_FACTOR_REVERSAL_FAILURE NONE\n";
    }

    if (have_first_same_stream_failure) {
        std::cout
            << "\nFIRST_SAME_STREAM_REVERSAL_FAILURE\n";

        print_excursion(
            first_same_stream_failure
        );
    } else {
        std::cout
            << "\nFIRST_SAME_STREAM_REVERSAL_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}