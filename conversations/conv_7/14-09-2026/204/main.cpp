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

struct OmegaStats
{
    long long count = 0;
    long long index_one = 0;
    long long index_two = 0;
    long long index_three = 0;
    long long index_four_plus = 0;

    long long total_omega = 0;
    long long max_index = 0;

    long long distinct_one = 0;
    long long distinct_two = 0;
    long long distinct_three_plus = 0;
};

vector<int> build_primes(int limit);
vector<int> build_spf(int limit);

Factorization factorize(int n, const vector<int>& spf);

bool is_prime(int n);
vector<int> factorize_distinct(int n);

int mod_pow(int a, int e, int mod);
bool is_primitive_root(int g, int p);
int primitive_root_prime(int p);

vector<int> build_discrete_log_table(
    int g,
    int p
);

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

void print_factorization(
    const Factorization& f
);

void print_sequence(
    const vector<int>& values
);

void print_omega_table(
    const vector<OmegaStats>& stats
);

void main_experiment();

int main()
{
    cout << "START EXPERIMENT 499\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 499\n";
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

Factorization factorize(
    int n,
    const vector<int>& spf
)
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

int mod_pow(
    int a,
    int e,
    int mod
)
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

bool is_primitive_root(
    int g,
    int p
)
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
    vector<int> logs(p, -1);

    if (p == 2)
    {
        logs[1] = 0;
        return logs;
    }

    long long value = 1;

    for (int e = 0; e < p - 1; ++e)
    {
        logs[static_cast<int>(value)] = e;
        value = (value * g) % p;
    }

    return logs;
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
    int best = 1000000000;

    for (int u = 1; u < modulus; ++u)
    {
        if (gcd(u, modulus) != 1)
        {
            continue;
        }

        int sum = 0;

        for (int exponent : exponents)
        {
            int x =
                (u * exponent) % modulus;

            if (x == 0)
            {
                return -1;
            }

            sum += x;
        }

        if (sum % modulus != 0)
        {
            continue;
        }

        best = min(best, sum / modulus);
    }

    return best;
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

void print_omega_table(
    const vector<OmegaStats>& stats
)
{
    cout << "\nOMEGA_CLASSIFICATION\n";

    for (size_t omega = 0;
         omega < stats.size();
         ++omega)
    {
        const OmegaStats& s = stats[omega];

        if (s.count == 0)
        {
            continue;
        }

        const double index_one_percent =
            100.0 *
            static_cast<double>(s.index_one) /
            static_cast<double>(s.count);

        const double index_two_percent =
            100.0 *
            static_cast<double>(s.index_two) /
            static_cast<double>(s.count);

        const double index_three_percent =
            100.0 *
            static_cast<double>(s.index_three) /
            static_cast<double>(s.count);

        const double index_four_plus_percent =
            100.0 *
            static_cast<double>(s.index_four_plus) /
            static_cast<double>(s.count);

        const double avg_index =
            static_cast<double>(s.index_one +
                                2 * s.index_two +
                                3 * s.index_three +
                                4 * s.index_four_plus) /
            static_cast<double>(s.count);

        cout << "OMEGA=" << omega
             << " COUNT=" << s.count
             << " INDEX1=" << s.index_one
             << " INDEX1_PERCENT="
             << index_one_percent
             << " INDEX2=" << s.index_two
             << " INDEX2_PERCENT="
             << index_two_percent
             << " INDEX3=" << s.index_three
             << " INDEX3_PERCENT="
             << index_three_percent
             << " INDEX4PLUS="
             << s.index_four_plus
             << " INDEX4PLUS_PERCENT="
             << index_four_plus_percent
             << " MAX_INDEX="
             << s.max_index
             << " AVG_INDEX="
             << avg_index
             << "\n";
    }
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

    /*
     * We need room for omega larger than 21 just in case.
     */
    vector<OmegaStats> omega_stats(128);

    long long total_failures = 0;
    long long total_omega = 0;
    long long total_index = 0;

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

        long long failures = 0;

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

            if (has_proper_one_residue_divisor(
                    n,
                    j,
                    divisors))
            {
                continue;
            }

            vector<int> exponents;

            for (size_t i = 0;
                 i < f.primes.size();
                 ++i)
            {
                const int residue =
                    f.primes[i] % j;

                const int exponent =
                    log_table[residue];

                if (exponent <= 0)
                {
                    cout << "INVALID_EXPONENT"
                         << " J=" << j
                         << " P=" << p
                         << " FACTOR="
                         << f.primes[i]
                         << " EXPONENT="
                         << exponent
                         << "\n";

                    return;
                }

                for (int e = 0;
                     e < f.multiplicities[i];
                     ++e)
                {
                    exponents.push_back(exponent);
                }
            }

            int sum = 0;

            for (int exponent : exponents)
            {
                sum += exponent;
                sum %= modulus;
            }

            if (sum != 0)
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

            const int omega =
                static_cast<int>(exponents.size());

            const int index =
                sequence_index(
                    exponents,
                    modulus
                );

            if (index <= 0)
            {
                cout << "INDEX_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << "\n";

                return;
            }

            const int distinct =
                distinct_value_count(exponents);

            OmegaStats& s =
                omega_stats[omega];

            ++s.count;
            ++s.total_omega;

            s.max_index =
                max<long long>(
                    s.max_index,
                    index
                );

            if (index == 1)
            {
                ++s.index_one;
            }
            else if (index == 2)
            {
                ++s.index_two;
            }
            else if (index == 3)
            {
                ++s.index_three;
            }
            else
            {
                ++s.index_four_plus;
            }

            if (distinct == 1)
            {
                ++s.distinct_one;
            }
            else if (distinct == 2)
            {
                ++s.distinct_two;
            }
            else
            {
                ++s.distinct_three_plus;
            }

            ++failures;
            ++total_failures;

            total_omega += omega;
            total_index += index;

            /*
             * Print the first few omega=3 cases and the first
             * few index>1 cases. These are the most useful
             * examples for identifying an exact rule.
             */
            if (omega == 3 &&
                s.count <= 3)
            {
                cout << "OMEGA3_EXAMPLE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << " INDEX=" << index
                     << " EXPONENTS=";

                print_sequence(exponents);

                cout << "\n";
            }

            if (index > 1 &&
                index <= 2 &&
                s.index_two <= 3)
            {
                cout << "INDEX2_EXAMPLE"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << " OMEGA=" << omega
                     << " DISTINCT="
                     << distinct
                     << " INDEX=" << index
                     << " EXPONENTS=";

                print_sequence(exponents);

                cout << "\n";
            }
        }

        const double percent =
            failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(failures) /
                  static_cast<double>(primes.size());

        cout << "J=" << j
             << " MODULUS=" << modulus
             << " FAILURES=" << failures
             << " FAILURE_PERCENT_OF_PRIMES="
             << percent
             << "\n";
    }

    print_omega_table(omega_stats);

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

    cout << "\nTOTAL_FAILURES="
         << total_failures
         << "\n";

    cout << "TOTAL_OMEGA="
         << total_omega
         << "\n";

    cout << "TOTAL_INDEX="
         << total_index
         << "\n";

    cout << "GLOBAL_AVG_OMEGA="
         << global_avg_omega
         << "\n";

    cout << "GLOBAL_AVG_INDEX="
         << global_avg_index
         << "\n";

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
