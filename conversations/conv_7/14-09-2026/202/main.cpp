#include <algorithm>
#include <chrono>
#include <cmath>
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

struct FailureStats
{
    long long failures = 0;
    long long zero_sum_ok = 0;
    long long minimal_zero_sum_ok = 0;
    long long length_bound_ok = 0;
    long long max_length_hits = 0;
    long long total_omega = 0;
};

vector<int> build_primes(int limit);
vector<int> build_spf(int limit);

Factorization factorize(int n, const vector<int>& spf);

bool is_prime(int n);

vector<int> factorize_small(int n);

int mod_pow(int a, int e, int mod);

bool is_primitive_root(int g, int prime);

int primitive_root_prime(int prime);

vector<int> build_discrete_log_table(int primitive_root, int prime);

vector<int> all_divisors_from_factorization(const Factorization& f);

bool has_proper_one_residue_divisor(
    int n,
    int j,
    const vector<int>& divisors
);

bool has_proper_zero_sum_subset(
    const vector<int>& exponents,
    int modulus
);

void print_histogram(
    const unordered_map<int, long long>& histogram
);

void print_example(
    int j,
    int p,
    int n,
    const Factorization& f,
    const vector<int>& residues,
    const vector<int>& exponents
);

void main_experiment();

int main()
{
    cout << "START EXPERIMENT 496\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 496\n";
    return 0;
}

vector<int> build_primes(int limit)
{
    vector<bool> composite(limit + 1, false);
    vector<int> primes;

    for (int i = 2; i <= limit; ++i)
    {
        if (!composite[i])
        {
            primes.push_back(i);

            if (1LL * i * i <= limit)
            {
                for (int j = i * i; j <= limit; j += i)
                {
                    composite[j] = true;
                }
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
        if (spf[i] == 0)
        {
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

vector<int> factorize_small(int n)
{
    vector<int> factors;

    for (int p = 2; 1LL * p * p <= n; ++p)
    {
        if (n % p == 0)
        {
            factors.push_back(p);

            while (n % p == 0)
            {
                n /= p;
            }
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

bool is_primitive_root(int g, int prime)
{
    const int phi = prime - 1;
    const vector<int> factors = factorize_small(phi);

    for (int q : factors)
    {
        if (mod_pow(g, phi / q, prime) == 1)
        {
            return false;
        }
    }

    return true;
}

int primitive_root_prime(int prime)
{
    for (int g = 2; g < prime; ++g)
    {
        if (is_primitive_root(g, prime))
        {
            return g;
        }
    }

    return -1;
}

vector<int> build_discrete_log_table(int primitive_root, int prime)
{
    vector<int> log_table(prime, -1);

    long long value = 1;

    for (int exponent = 0; exponent < prime - 1; ++exponent)
    {
        log_table[static_cast<int>(value)] = exponent;
        value = (value * primitive_root) % prime;
    }

    return log_table;
}

vector<int> all_divisors_from_factorization(const Factorization& f)
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
    const int omega = static_cast<int>(exponents.size());

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

void print_example(
    int j,
    int p,
    int n,
    const Factorization& f,
    const vector<int>& residues,
    const vector<int>& exponents
)
{
    cout << "EXAMPLE"
         << " J=" << j
         << " P=" << p
         << " N=" << n
         << " FACTORS=";

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

    cout << " RESIDUES=[";

    for (size_t i = 0; i < residues.size(); ++i)
    {
        if (i != 0)
        {
            cout << ",";
        }

        cout << residues[i];
    }

    cout << "] EXPONENTS=[";

    for (size_t i = 0; i < exponents.size(); ++i)
    {
        if (i != 0)
        {
            cout << ",";
        }

        cout << exponents[i];
    }

    cout << "]\n";
}

void main_experiment()
{
    const auto start = chrono::high_resolution_clock::now();

    const vector<int> primes = build_primes(PRIME_LIMIT);
    const int max_n = J_LIMIT * PRIME_LIMIT + 1;
    const vector<int> spf = build_spf(max_n);

    cout << "PRIME_COUNT=" << primes.size() << "\n";
    cout << "J_LIMIT=" << J_LIMIT << "\n";
    cout << "MAX_JP_PLUS_1=" << max_n << "\n";
    cout << "SPF_READY=1\n";

    long long total_failures = 0;
    long long total_zero_sum_ok = 0;
    long long total_minimal_ok = 0;
    long long total_length_bound_ok = 0;
    long long total_max_length_hits = 0;
    long long total_omega = 0;

    unordered_map<int, long long> global_omega_histogram;

    for (int j = 2; j <= J_LIMIT; ++j)
    {
        if (!is_prime(j))
        {
            continue;
        }

        const int modulus = j - 1;
        const int primitive_root = primitive_root_prime(j);

        if (primitive_root < 0)
        {
            cout << "PRIMITIVE_ROOT_FAILURE J=" << j << "\n";
            return;
        }

        const vector<int> log_table =
            build_discrete_log_table(
                primitive_root,
                j
            );

        FailureStats stats;

        unordered_map<int, long long> omega_histogram;

        int printed_examples = 0;

        for (int p : primes)
        {
            const long long n64 = 1LL * j * p + 1;

            if (n64 > max_n)
            {
                break;
            }

            const int n = static_cast<int>(n64);

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

            ++stats.failures;
            ++total_failures;

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
                         << " FACTOR=" << f.primes[i]
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

            const int omega =
                static_cast<int>(exponents.size());

            ++stats.total_omega;
            ++total_omega;
            ++omega_histogram[omega];
            ++global_omega_histogram[omega];

            int total_exponent = 0;

            for (int exponent : exponents)
            {
                total_exponent =
                    (total_exponent + exponent) % modulus;
            }

            const bool zero_sum_ok =
                (total_exponent == 0);

            if (zero_sum_ok)
            {
                ++stats.zero_sum_ok;
                ++total_zero_sum_ok;
            }

            const bool proper_zero_sum =
                has_proper_zero_sum_subset(
                    exponents,
                    modulus
                );

            const bool minimal_zero_sum_ok =
                zero_sum_ok && !proper_zero_sum;

            if (minimal_zero_sum_ok)
            {
                ++stats.minimal_zero_sum_ok;
                ++total_minimal_ok;
            }

            const bool length_bound_ok =
                omega <= modulus;

            if (length_bound_ok)
            {
                ++stats.length_bound_ok;
                ++total_length_bound_ok;
            }

            if (omega == modulus)
            {
                ++stats.max_length_hits;
                ++total_max_length_hits;
            }

            if (printed_examples < 5 &&
                (omega == modulus ||
                 omega >= modulus - 1))
            {
                print_example(
                    j,
                    p,
                    n,
                    f,
                    residues,
                    exponents
                );

                ++printed_examples;
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

        const double zero_sum_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.zero_sum_ok
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        const double minimal_percent =
            stats.failures == 0
                ? 0.0
                : 100.0 *
                  static_cast<double>(
                      stats.minimal_zero_sum_ok
                  ) /
                  static_cast<double>(
                      stats.failures
                  );

        cout << "J=" << j
             << " PRIMITIVE_ROOT=" << primitive_root
             << " MODULUS=" << modulus
             << " FAILURES=" << stats.failures
             << " ZERO_SUM_OK=" << stats.zero_sum_ok
             << " ZERO_SUM_PERCENT="
             << zero_sum_percent
             << " MINIMAL_ZERO_SUM_OK="
             << stats.minimal_zero_sum_ok
             << " MINIMAL_PERCENT="
             << minimal_percent
             << " LENGTH_BOUND_OK="
             << stats.length_bound_ok
             << " MAX_LENGTH_HITS="
             << stats.max_length_hits
             << " AVG_OMEGA="
             << avg_omega
             << "\n";

        cout << "J=" << j << " OMEGA_HIST=";
        print_histogram(omega_histogram);
    }

    cout << "\nGLOBAL_FAILURE_OMEGA_HIST=";
    print_histogram(global_omega_histogram);

    cout << "TOTAL_FAILURES="
         << total_failures << "\n";

    cout << "TOTAL_ZERO_SUM_OK="
         << total_zero_sum_ok << "\n";

    cout << "TOTAL_MINIMAL_ZERO_SUM_OK="
         << total_minimal_ok << "\n";

    cout << "TOTAL_LENGTH_BOUND_OK="
         << total_length_bound_ok << "\n";

    cout << "TOTAL_MAX_LENGTH_HITS="
         << total_max_length_hits << "\n";

    const double global_avg_omega =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_omega) /
              static_cast<double>(total_failures);

    cout << "GLOBAL_AVG_OMEGA="
         << global_avg_omega << "\n";

    cout << "ZERO_SUM_SANITY="
         << (total_zero_sum_ok == total_failures ? 1 : 0)
         << "\n";

    cout << "MINIMAL_ZERO_SUM_SANITY="
         << (total_minimal_ok == total_failures ? 1 : 0)
         << "\n";

    cout << "DAVELENGTH_SANITY="
         << (total_length_bound_ok == total_failures ? 1 : 0)
         << "\n";

    const auto end = chrono::high_resolution_clock::now();

    const double elapsed_ms =
        chrono::duration<double, milli>(
            end - start
        ).count();

    cout << "ELAPSED_TIME_MS="
         << elapsed_ms << "\n";
}
