#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

using namespace std;

static const int PRIME_LIMIT = 100000;
static const int J_LIMIT = 97;

struct Factorization
{
    vector<int> primes;
    vector<int> multiplicities;
};

struct SequenceStats
{
    long long failures = 0;
    long long total_omega = 0;

    long long pair_sequences = 0;
    long long pure_sequences = 0;
    long long max_length_sequences = 0;

    long long index_one = 0;
    long long index_greater_one = 0;

    long long distinct_exponent_one = 0;
    long long distinct_exponent_two = 0;
    long long distinct_exponent_three_or_more = 0;

    long long total_index = 0;
};

vector<int> build_primes(int limit);
vector<int> build_spf(int limit);

Factorization factorize(int n, const vector<int>& spf);

vector<int> factorize_distinct(int n);

int mod_pow(int a, int e, int mod);
bool is_prime(int n);
bool is_primitive_root(int g, int p);
int primitive_root_prime(int p);
vector<int> build_discrete_log_table(int g, int p);

vector<int> all_divisors_from_factorization(
    const Factorization& f
);

bool has_proper_one_residue_divisor(
    int n,
    int j,
    const vector<int>& divisors
);

bool has_proper_zero_sum_subset(
    const vector<int>& exponents,
    int modulus
);

int sequence_index(
    const vector<int>& exponents,
    int modulus
);

int distinct_value_count(
    const vector<int>& values
);

bool all_equal(
    const vector<int>& values
);

void print_factorization(
    const Factorization& f
);

void print_sequence(
    const vector<int>& values
);

void print_histogram(
    const unordered_map<int, long long>& histogram
);

void main_experiment();

int main()
{
    cout << "START EXPERIMENT 498\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 498\n";
    return 0;
}

vector<int> build_primes(int limit)
{
    vector<bool> composite(limit + 1, false);
    vector<int> primes;

    for (int i = 2; i <= limit; ++i)
    {
        if (composite[i])
        {
            continue;
        }

        primes.push_back(i);

        if (1LL * i * i <= limit)
        {
            for (int j = i * i; j <= limit; j += i)
            {
                composite[j] = true;
            }
        }
    }

    return primes;
}

vector<int> build_spf(int limit)
{
    vector<int> spf(limit + 1, 0);

    for (int i = 2; i <= limit; ++i)
    {
        if (spf[i] != 0)
        {
            continue;
        }

        spf[i] = i;

        if (1LL * i * i <= limit)
        {
            for (int j = i * i; j <= limit; j += i)
            {
                if (spf[j] == 0)
                {
                    spf[j] = i;
                }
            }
        }
    }

    return spf;
}

Factorization factorize(int n, const vector<int>& spf)
{
    Factorization result;

    while (n > 1)
    {
        const int p = spf[n];

        int exponent = 0;

        while (n % p == 0)
        {
            n /= p;
            ++exponent;
        }

        result.primes.push_back(p);
        result.multiplicities.push_back(exponent);
    }

    return result;
}

vector<int> factorize_distinct(int n)
{
    vector<int> factors;

    for (int p = 2; 1LL * p * p <= n; ++p)
    {
        if (n % p != 0)
        {
            continue;
        }

        factors.push_back(p);

        while (n % p == 0)
        {
            n /= p;
        }
    }

    if (n > 1)
    {
        factors.push_back(n);
    }

    return factors;
}

int mod_pow(int a, int e, int mod)
{
    long long result = 1;
    long long base = a % mod;

    while (e > 0)
    {
        if (e & 1)
        {
            result = (result * base) % mod;
        }

        base = (base * base) % mod;
        e >>= 1;
    }

    return static_cast<int>(result);
}

bool is_prime(int n)
{
    if (n < 2)
    {
        return false;
    }

    if (n % 2 == 0)
    {
        return n == 2;
    }

    for (int d = 3; 1LL * d * d <= n; d += 2)
    {
        if (n % d == 0)
        {
            return false;
        }
    }

    return true;
}

bool is_primitive_root(int g, int p)
{
    if (p == 2)
    {
        return g == 1;
    }

    const int phi = p - 1;
    const vector<int> factors =
        factorize_distinct(phi);

    for (int q : factors)
    {
        if (mod_pow(g, phi / q, p) == 1)
        {
            return false;
        }
    }

    return true;
}

int primitive_root_prime(int p)
{
    if (p == 2)
    {
        return 1;
    }

    for (int g = 2; g < p; ++g)
    {
        if (is_primitive_root(g, p))
        {
            return g;
        }
    }

    return -1;
}

vector<int> build_discrete_log_table(
    int g,
    int p
)
{
    vector<int> log_table(p, -1);

    if (p == 2)
    {
        log_table[1] = 0;
        return log_table;
    }

    long long value = 1;

    for (int e = 0; e < p - 1; ++e)
    {
        log_table[static_cast<int>(value)] = e;
        value = (value * g) % p;
    }

    return log_table;
}

vector<int> all_divisors_from_factorization(
    const Factorization& f
)
{
    vector<int> divisors = {1};

    for (size_t i = 0; i < f.primes.size(); ++i)
    {
        const int p = f.primes[i];
        const int exponent = f.multiplicities[i];

        vector<int> next;

        int power = 1;

        for (int e = 0; e <= exponent; ++e)
        {
            for (int d : divisors)
            {
                next.push_back(d * power);
            }

            power *= p;
        }

        divisors.swap(next);
    }

    sort(divisors.begin(), divisors.end());

    return divisors;
}

bool has_proper_one_residue_divisor(
    int n,
    int j,
    const vector<int>& divisors
)
{
    for (int d : divisors)
    {
        if (d <= 1 || d >= n)
        {
            continue;
        }

        if (d % j == 1)
        {
            return true;
        }
    }

    return false;
}

bool has_proper_zero_sum_subset(
    const vector<int>& exponents,
    int modulus
)
{
    const int omega =
        static_cast<int>(exponents.size());

    if (omega <= 1)
    {
        return false;
    }

    vector<vector<bool>> possible(
        omega + 1,
        vector<bool>(modulus, false)
    );

    possible[0][0] = true;

    for (int exponent : exponents)
    {
        vector<vector<bool>> next = possible;

        for (int count = 0; count < omega; ++count)
        {
            for (int sum = 0; sum < modulus; ++sum)
            {
                if (!possible[count][sum])
                {
                    continue;
                }

                const int new_sum =
                    (sum + exponent) % modulus;

                next[count + 1][new_sum] = true;
            }
        }

        possible.swap(next);
    }

    for (int count = 1; count < omega; ++count)
    {
        if (possible[count][0])
        {
            return true;
        }
    }

    return false;
}

int sequence_index(
    const vector<int>& exponents,
    int modulus
)
{
    /*
     * For a sequence in C_modulus, multiply all
     * terms by every unit u modulo modulus.
     *
     * For the resulting representatives in
     * {1,...,modulus-1}, calculate:
     *
     *      sum / modulus.
     *
     * The minimum is the cyclic zero-sum index.
     */

    int best_index = 1000000000;

    for (int u = 1; u < modulus; ++u)
    {
        if (gcd(u, modulus) != 1)
        {
            continue;
        }

        int sum = 0;

        for (int exponent : exponents)
        {
            const int x =
                (u * exponent) % modulus;

            /*
             * x cannot be zero because exponent is a
             * nonzero element of C_modulus and u is a unit.
             */
            if (x == 0)
            {
                continue;
            }

            sum += x;
        }

        if (sum % modulus != 0)
        {
            /*
             * Should never happen for a zero-sum sequence.
             */
            continue;
        }

        best_index =
            min(best_index, sum / modulus);
    }

    return best_index;
}

int distinct_value_count(
    const vector<int>& values
)
{
    if (values.empty())
    {
        return 0;
    }

    vector<int> copy = values;

    sort(copy.begin(), copy.end());

    copy.erase(
        unique(copy.begin(), copy.end()),
        copy.end()
    );

    return static_cast<int>(copy.size());
}

bool all_equal(
    const vector<int>& values
)
{
    if (values.empty())
    {
        return true;
    }

    return all_of(
        values.begin(),
        values.end(),
        [&](int x)
        {
            return x == values.front();
        }
    );
}

void print_factorization(
    const Factorization& f
)
{
    for (size_t i = 0; i < f.primes.size(); ++i)
    {
        if (i != 0)
        {
            cout << "*";
        }

        cout << f.primes[i];

        if (f.multiplicities[i] > 1)
        {
            cout << "^" << f.multiplicities[i];
        }
    }
}

void print_sequence(
    const vector<int>& values
)
{
    cout << "[";

    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0)
        {
            cout << ",";
        }

        cout << values[i];
    }

    cout << "]";
}

void print_histogram(
    const unordered_map<int, long long>& histogram
)
{
    vector<pair<int, long long>> values(
        histogram.begin(),
        histogram.end()
    );

    sort(
        values.begin(),
        values.end(),
        [](const auto& a, const auto& b)
        {
            return a.first < b.first;
        }
    );

    for (const auto& [key, value] : values)
    {
        cout << key << ":" << value << " ";
    }

    cout << "\n";
}

void main_experiment()
{
    const auto start =
        chrono::high_resolution_clock::now();

    const vector<int> primes =
        build_primes(PRIME_LIMIT);

    const int max_n =
        J_LIMIT * PRIME_LIMIT + 1;

    const vector<int> spf =
        build_spf(max_n);

    cout << "PRIME_COUNT="
         << primes.size() << "\n";

    cout << "J_LIMIT="
         << J_LIMIT << "\n";

    cout << "MAX_JP_PLUS_1="
         << max_n << "\n";

    cout << "SPF_READY=1\n";

    long long total_failures = 0;
    long long total_omega = 0;

    long long total_pair_sequences = 0;
    long long total_pure_sequences = 0;
    long long total_max_length = 0;

    long long total_index_one = 0;
    long long total_index_greater_one = 0;
    long long total_index = 0;

    long long total_distinct_one = 0;
    long long total_distinct_two = 0;
    long long total_distinct_three_plus = 0;

    unordered_map<int, long long> global_omega_hist;
    unordered_map<int, long long> global_index_hist;
    unordered_map<int, long long> global_distinct_hist;

    for (int j = 2; j <= J_LIMIT; ++j)
    {
        if (!is_prime(j))
        {
            continue;
        }

        const int modulus = j - 1;

        const int primitive_root =
            primitive_root_prime(j);

        if (primitive_root < 0)
        {
            cout << "PRIMITIVE_ROOT_FAILURE J="
                 << j
                 << "\n";

            return;
        }

        const vector<int> log_table =
            build_discrete_log_table(
                primitive_root,
                j
            );

        SequenceStats stats;

        unordered_map<int, long long> omega_hist;
        unordered_map<int, long long> index_hist;
        unordered_map<int, long long> distinct_hist;

        int printed_index_gt_one = 0;
        int printed_pure = 0;
        int printed_nonpure_max = 0;

        for (int p : primes)
        {
            const long long n64 =
                1LL * j * p + 1;

            if (n64 > max_n)
            {
                break;
            }

            const int n =
                static_cast<int>(n64);

            if (spf[n] == n)
            {
                continue;
            }

            const Factorization f =
                factorize(n, spf);

            const vector<int> divisors =
                all_divisors_from_factorization(f);

            const bool qualifying =
                has_proper_one_residue_divisor(
                    n,
                    j,
                    divisors
                );

            if (qualifying)
            {
                continue;
            }

            vector<int> residues;
            vector<int> exponents;

            for (size_t i = 0; i < f.primes.size(); ++i)
            {
                const int residue =
                    f.primes[i] % j;

                const int exponent =
                    log_table[residue];

                if (exponent < 0)
                {
                    cout << "DISCRETE_LOG_FAILURE"
                         << " J=" << j
                         << " P=" << p
                         << " FACTOR="
                         << f.primes[i]
                         << "\n";

                    return;
                }

                for (int e = 0;
                     e < f.multiplicities[i];
                     ++e)
                {
                    residues.push_back(residue);
                    exponents.push_back(exponent);
                }
            }

            /*
             * Reconfirm that this really is a minimal
             * zero-sum sequence.
             */
            int total_exponent = 0;

            for (int x : exponents)
            {
                total_exponent += x;
                total_exponent %= modulus;
            }

            if (total_exponent != 0)
            {
                cout << "ZERO_SUM_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << "\n";

                return;
            }

            if (has_proper_zero_sum_subset(
                    exponents,
                    modulus))
            {
                cout << "MINIMALITY_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << "\n";

                return;
            }

            ++stats.failures;
            ++total_failures;

            const int omega =
                static_cast<int>(exponents.size());

            ++stats.total_omega;
            ++total_omega;

            ++omega_hist[omega];
            ++global_omega_hist[omega];

            const int distinct =
                distinct_value_count(exponents);

            ++distinct_hist[distinct];
            ++global_distinct_hist[distinct];

            const bool pair =
                omega == 2;

            const bool pure =
                all_equal(exponents);

            const bool max_length =
                omega == modulus;

            if (pair)
            {
                ++stats.pair_sequences;
                ++total_pair_sequences;
            }

            if (pure)
            {
                ++stats.pure_sequences;
                ++total_pure_sequences;
            }

            if (max_length)
            {
                ++stats.max_length_sequences;
                ++total_max_length;

                if (!pure &&
                    printed_nonpure_max < 5)
                {
                    cout << "NONPURE_MAX_LENGTH"
                         << " J=" << j
                         << " P=" << p
                         << " N=" << n
                         << " MODULUS=" << modulus
                         << " EXPONENTS=";

                    print_sequence(exponents);

                    cout << "\n";

                    ++printed_nonpure_max;
                }
            }

            if (distinct == 1)
            {
                ++stats.distinct_exponent_one;
                ++total_distinct_one;
            }
            else if (distinct == 2)
            {
                ++stats.distinct_exponent_two;
                ++total_distinct_two;
            }
            else
            {
                ++stats.distinct_exponent_three_or_more;
                ++total_distinct_three_plus;
            }

            const int index =
                sequence_index(
                    exponents,
                    modulus
                );

            if (index == 1)
            {
                ++stats.index_one;
                ++total_index_one;
            }
            else
            {
                ++stats.index_greater_one;
                ++total_index_greater_one;

                if (printed_index_gt_one < 10)
                {
                    cout << "INDEX_GT_ONE"
                         << " J=" << j
                         << " P=" << p
                         << " N=" << n
                         << " MODULUS=" << modulus
                         << " INDEX=" << index
                         << " OMEGA=" << omega
                         << " EXPONENTS=";

                    print_sequence(exponents);

                    cout << "\n";

                    ++printed_index_gt_one;
                }
            }

            stats.total_index += index;
            total_index += index;
            ++index_hist[index];
            ++global_index_hist[index];

            if (pure && printed_pure < 3)
            {
                cout << "PURE_SEQUENCE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << " MODULUS=" << modulus
                     << " OMEGA=" << omega
                     << " EXPONENTS=";

                print_sequence(exponents);

                cout << " INDEX=" << index
                     << "\n";

                ++printed_pure;
            }
        }

        const double avg_omega =
            stats.failures == 0
                ? 0.0
                : static_cast<double>(
                      stats.total_omega
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double avg_index =
            stats.failures == 0
                ? 0.0
                : static_cast<double>(
                      stats.total_index
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double pair_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.pair_sequences
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double pure_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.pure_sequences
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double index_one_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.index_one
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double max_length_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.max_length_sequences
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        cout << "J=" << j
             << " MODULUS=" << modulus
             << " FAILURES=" << stats.failures
             << " AVG_OMEGA=" << avg_omega
             << " PAIRS=" << stats.pair_sequences
             << " PAIR_PERCENT=" << pair_percent
             << " PURE=" << stats.pure_sequences
             << " PURE_PERCENT=" << pure_percent
             << " INDEX_ONE=" << stats.index_one
             << " INDEX_ONE_PERCENT="
             << index_one_percent
             << " INDEX_GT_ONE="
             << stats.index_greater_one
             << " MAX_LENGTH="
             << stats.max_length_sequences
             << " MAX_LENGTH_PERCENT="
             << max_length_percent
             << " AVG_INDEX="
             << avg_index
             << "\n";

        cout << "J=" << j
             << " DISTINCT_EXPONENT_HIST=";

        print_histogram(distinct_hist);

        cout << "J=" << j
             << " INDEX_HIST=";

        print_histogram(index_hist);

        cout << "J=" << j
             << " OMEGA_HIST=";

        print_histogram(omega_hist);
    }

    cout << "\nTOTAL_FAILURES="
         << total_failures
         << "\n";

    cout << "TOTAL_PAIR_SEQUENCES="
         << total_pair_sequences
         << "\n";

    cout << "TOTAL_PURE_SEQUENCES="
         << total_pure_sequences
         << "\n";

    cout << "TOTAL_MAX_LENGTH="
         << total_max_length
         << "\n";

    cout << "TOTAL_INDEX_ONE="
         << total_index_one
         << "\n";

    cout << "TOTAL_INDEX_GT_ONE="
         << total_index_greater_one
         << "\n";

    cout << "TOTAL_DISTINCT_EXPONENT_ONE="
         << total_distinct_one
         << "\n";

    cout << "TOTAL_DISTINCT_EXPONENT_TWO="
         << total_distinct_two
         << "\n";

    cout << "TOTAL_DISTINCT_EXPONENT_THREE_PLUS="
         << total_distinct_three_plus
         << "\n";

    const double global_avg_omega =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_omega) /
              static_cast<double>(total_failures);

    const double global_avg_index =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_index) /
              static_cast<double>(total_failures);

    cout << "GLOBAL_AVG_OMEGA="
         << global_avg_omega
         << "\n";

    cout << "GLOBAL_AVG_INDEX="
         << global_avg_index
         << "\n";

    cout << "GLOBAL_OMEGA_HIST=";
    print_histogram(global_omega_hist);

    cout << "GLOBAL_INDEX_HIST=";
    print_histogram(global_index_hist);

    cout << "GLOBAL_DISTINCT_EXPONENT_HIST=";
    print_histogram(global_distinct_hist);

    cout << "MINIMAL_ZERO_SUM_VERIFIED=1\n";

    const auto end =
        chrono::high_resolution_clock::now();

    const double elapsed_ms =
        chrono::duration<double, milli>(
            end - start
        ).count();

    cout << "ELAPSED_TIME_MS="
         << elapsed_ms
         << "\n";
}