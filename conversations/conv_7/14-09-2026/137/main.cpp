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

    for (int x = 2; x <= limit; ++x) {
        if (is_prime[x]) {
            primes.push_back(
                static_cast<u64>(x)
            );
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

    return
        value % static_cast<i128>(k) == 0;
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

static bool exact_same_state(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.k == b.k &&
        a.t == b.t &&
        a.sign == b.sign;
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

static bool exact_palindrome(
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
        if (
            !exact_same_state(
                states[i],
                states[j]
            )
        ) {
            return false;
        }
    }

    return true;
}

static bool mirrored_state_path(
    const std::vector<Witness>& states
) {
    const std::size_t n = states.size();

    if (n < 2) {
        return true;
    }

    for (
        std::size_t i = 0;
        i < n / 2;
        ++i
    ) {
        const std::size_t j =
            n - 1 - i;

        if (
            !exact_same_state(
                states[i],
                states[j]
            )
        ) {
            return false;
        }
    }

    return true;
}

static bool projected_path_is_out_and_back(
    const std::vector<Witness>& states
) {
    if (states.size() < 3) {
        return true;
    }

    const std::size_t n =
        states.size();

    const std::size_t center =
        (n - 1) / 2;

    for (
        std::size_t i = 0;
        i <= center;
        ++i
    ) {
        const std::size_t j =
            n - 1 - i;

        if (
            states[i].m !=
            states[j].m
        ) {
            return false;
        }
    }

    return true;
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
        << ",value="
        << static_cast<u64>(value)
        << ")";
}

static void print_path(
    u64 r,
    const std::vector<Witness>& path
) {
    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        std::cout
            << "  K[" << i << "] ";

        print_witness(
            r,
            "W",
            path[i]
        );

        std::cout << "\n";
    }
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

static std::string projected_string(
    const std::vector<Witness>& path
) {
    std::string out;

    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        if (!out.empty()) {
            out += ",";
        }

        out += std::to_string(path[i].m);
    }

    return out;
}

static std::string signed_string(
    const std::vector<Witness>& path
) {
    std::string out;

    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        if (!out.empty()) {
            out += ",";
        }

        out += std::to_string(path[i].m);

        if (path[i].sign > 0) {
            out += "+";
        } else {
            out += "-";
        }
    }

    return out;
}

int main() {
    constexpr int EXPERIMENT = 429;
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

    u64 total_hidden_exact_changes = 0;
    u64 total_hidden_m_changes = 0;

    u64 exact_path_palindromes = 0;
    u64 projected_hidden_path_palindromes = 0;

    u64 projected_out_and_back = 0;

    u64 exact_path_nonpalindromes = 0;
    u64 projected_hidden_path_nonpalindromes = 0;

    u64 hidden_exact_change_only = 0;
    u64 hidden_m_change_only = 0;

    u64 full_hidden_paths = 0;

    u64 max_exact_changes = 0;
    u64 max_m_changes = 0;
    u64 max_path_length = 0;

    std::map<u64, u64>
        exact_change_histogram;

    std::map<u64, u64>
        m_change_histogram;

    std::map<u64, u64>
        path_length_histogram;

    Excursion first_exact_nonpalindrome;
    bool have_first_exact_nonpalindrome = false;

    Excursion first_projected_nonpalindrome;
    bool have_first_projected_nonpalindrome = false;

    Excursion first_exact_only_change;
    bool have_first_exact_only_change = false;

    Excursion first_m_change;
    bool have_first_m_change = false;

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

        /*
         * Full winner state at every K.
         */
        std::vector<Witness> winner_by_k(
            static_cast<std::size_t>(K_LIMIT) + 1
        );

        winner_by_k[1] =
            global_winner(best_by_m);

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
                    const Witness candidate =
                        make_witness(
                            r,
                            m,
                            K,
                            +1
                        );

                    if (
                        candidate.k >
                        best_by_m[mi].k
                    ) {
                        best_by_m[mi] =
                            candidate;
                    }
                }

                if (
                    divides_stream(
                        r,
                        m,
                        K,
                        -1
                    )
                ) {
                    const Witness candidate =
                        make_witness(
                            r,
                            m,
                            K,
                            -1
                        );

                    if (
                        candidate.k >
                        best_by_m[mi].k
                    ) {
                        best_by_m[mi] =
                            candidate;
                    }
                }
            }

            winner_by_k[
                static_cast<std::size_t>(K)
            ] =
                global_winner(best_by_m);
        }

        /*
         * Project the complete path onto m and
         * find the 1 -> ... -> 1 excursions.
         */
        std::vector<Witness> transition_states;
        transition_states.reserve(256);

        transition_states.push_back(
            winner_by_k[1]
        );

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            const Witness& prev =
                winner_by_k[
                    static_cast<std::size_t>(K - 1)
                ];

            const Witness& curr =
                winner_by_k[
                    static_cast<std::size_t>(K)
                ];

            if (curr.m != prev.m) {
                transition_states.push_back(
                    curr
                );
            }
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
                start >= transition_states.size()
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
                end >= transition_states.size()
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
                    !projected_palindrome(
                        ex.states
                    )
                ) {
                    start = end;
                    continue;
                }

                ++total_mirrored_pairs;

                /*
                 * Endpoints are transition states.
                 * Recover their actual K coordinates.
                 */
                const u64 left_k =
                    ex.states.front().k;

                const u64 right_k =
                    ex.states.back().k;

                if (
                    right_k <= left_k
                ) {
                    start = end;
                    continue;
                }

                /*
                 * Full winner path:
                 *
                 * K = left_k ... right_k
                 */
                std::vector<Witness> full_path;

                full_path.reserve(
                    static_cast<std::size_t>(
                        right_k - left_k + 1
                    )
                );

                for (
                    u64 K = left_k;
                    K <= right_k;
                    ++K
                ) {
                    full_path.push_back(
                        winner_by_k[
                            static_cast<std::size_t>(K)
                        ]
                    );
                }

                const u64 exact_changes =
                    full_path.size() <= 1
                        ? 0
                        : static_cast<u64>(
                            full_path.size() - 1
                        );

                u64 m_changes = 0;

                for (
                    std::size_t i = 1;
                    i < full_path.size();
                    ++i
                ) {
                    if (
                        full_path[i].m !=
                        full_path[i - 1].m
                    ) {
                        ++m_changes;
                    }
                }

                /*
                 * Convert exact winner path into its
                 * state-change path only.
                 */
                std::vector<Witness> state_change_path;
                state_change_path.reserve(
                    full_path.size()
                );

                if (!full_path.empty()) {
                    state_change_path.push_back(
                        full_path.front()
                    );

                    for (
                        std::size_t i = 1;
                        i < full_path.size();
                        ++i
                    ) {
                        if (
                            !exact_same_state(
                                full_path[i],
                                full_path[i - 1]
                            )
                        ) {
                            state_change_path.push_back(
                                full_path[i]
                            );
                        }
                    }
                }

                /*
                 * We compare the entire K-by-K winner
                 * path with its reverse.
                 */
                const bool exact_pal =
                    exact_palindrome(
                        full_path
                    );

                const bool projected_pal =
                    projected_palindrome(
                        full_path
                    );

                const bool out_and_back =
                    projected_path_is_out_and_back(
                        full_path
                    );

                if (exact_pal) {
                    ++exact_path_palindromes;
                } else {
                    ++exact_path_nonpalindromes;

                    if (
                        !have_first_exact_nonpalindrome
                    ) {
                        first_exact_nonpalindrome =
                            ex;

                        have_first_exact_nonpalindrome =
                            true;
                    }
                }

                if (projected_pal) {
                    ++projected_hidden_path_palindromes;
                } else {
                    ++projected_hidden_path_nonpalindromes;

                    if (
                        !have_first_projected_nonpalindrome
                    ) {
                        first_projected_nonpalindrome =
                            ex;

                        have_first_projected_nonpalindrome =
                            true;
                    }
                }

                if (out_and_back) {
                    ++projected_out_and_back;
                }

                ++full_hidden_paths;

                total_hidden_exact_changes +=
                    exact_changes;

                total_hidden_m_changes +=
                    m_changes;

                if (exact_changes > max_exact_changes) {
                    max_exact_changes =
                        exact_changes;
                }

                if (m_changes > max_m_changes) {
                    max_m_changes =
                        m_changes;
                }

                if (
                    full_path.size() >
                    max_path_length
                ) {
                    max_path_length =
                        full_path.size();
                }

                exact_change_histogram[
                    exact_changes
                ]++;

                m_change_histogram[
                    m_changes
                ]++;

                path_length_histogram[
                    static_cast<u64>(
                        full_path.size()
                    )
                ]++;

                if (
                    exact_changes >
                    m_changes
                ) {
                    ++hidden_exact_change_only;

                    if (
                        !have_first_exact_only_change
                    ) {
                        first_exact_only_change =
                            ex;

                        have_first_exact_only_change =
                            true;
                    }
                }

                if (m_changes > 0) {
                    ++hidden_m_change_only;

                    if (!have_first_m_change) {
                        first_m_change = ex;
                        have_first_m_change = true;
                    }
                }

                /*
                 * Print nothing here; only first failures
                 * are retained.
                 */
                (void)state_change_path;
            }

            start = end;
        }
    }

    std::cout
        << "\nTOTAL_EXCURSIONS="
        << total_excursions
        << "\n";

    std::cout
        << "FULL_HIDDEN_PATHS="
        << full_hidden_paths
        << "\n";

    std::cout
        << "TOTAL_MIRRORED_PAIRS="
        << total_mirrored_pairs
        << "\n";

    std::cout
        << "\nHIDDEN_PATH_COUNTS\n";

    std::cout
        << "TOTAL_EXACT_STATE_CHANGES="
        << total_hidden_exact_changes
        << "\n";

    std::cout
        << "TOTAL_M_CHANGES="
        << total_hidden_m_changes
        << "\n";

    std::cout
        << "MAX_EXACT_STATE_CHANGES="
        << max_exact_changes
        << "\n";

    std::cout
        << "MAX_M_CHANGES="
        << max_m_changes
        << "\n";

    std::cout
        << "MAX_FULL_PATH_LENGTH="
        << max_path_length
        << "\n";

    std::cout
        << "\nPATH_SYMMETRY\n";

    std::cout
        << "EXACT_PATH_PALINDROMES="
        << exact_path_palindromes
        << "\n";

    std::cout
        << "EXACT_PATH_NONPALINDROMES="
        << exact_path_nonpalindromes
        << "\n";

    std::cout
        << "PROJECTED_HIDDEN_PATH_PALINDROMES="
        << projected_hidden_path_palindromes
        << "\n";

    std::cout
        << "PROJECTED_HIDDEN_PATH_NONPALINDROMES="
        << projected_hidden_path_nonpalindromes
        << "\n";

    std::cout
        << "PROJECTED_OUT_AND_BACK="
        << projected_out_and_back
        << "\n";

    std::cout
        << "\nCHANGE_RELATION\n";

    std::cout
        << "EXACT_CHANGES_EXCEED_M_CHANGES="
        << hidden_exact_change_only
        << "\n";

    std::cout
        << "PATHS_WITH_M_CHANGES="
        << hidden_m_change_only
        << "\n";

    std::cout
        << "\nEXACT_CHANGE_HISTOGRAM\n";

    for (
        const auto& entry :
        exact_change_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM_CHANGE_HISTOGRAM\n";

    for (
        const auto& entry :
        m_change_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFULL_PATH_LENGTH_HISTOGRAM\n";

    for (
        const auto& entry :
        path_length_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_exact_nonpalindrome) {
        std::cout
            << "\nFIRST_EXACT_PATH_NONPALINDROME\n";

        print_excursion(
            first_exact_nonpalindrome
        );
    } else {
        std::cout
            << "\nFIRST_EXACT_PATH_NONPALINDROME NONE\n";
    }

    if (have_first_projected_nonpalindrome) {
        std::cout
            << "\nFIRST_PROJECTED_HIDDEN_PATH_NONPALINDROME\n";

        print_excursion(
            first_projected_nonpalindrome
        );
    } else {
        std::cout
            << "\nFIRST_PROJECTED_HIDDEN_PATH_NONPALINDROME NONE\n";
    }

    if (have_first_exact_only_change) {
        std::cout
            << "\nFIRST_EXACT_ONLY_CHANGE\n";

        print_excursion(
            first_exact_only_change
        );
    } else {
        std::cout
            << "\nFIRST_EXACT_ONLY_CHANGE NONE\n";
    }

    if (have_first_m_change) {
        std::cout
            << "\nFIRST_M_CHANGE\n";

        print_excursion(
            first_m_change
        );
    } else {
        std::cout
            << "\nFIRST_M_CHANGE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
