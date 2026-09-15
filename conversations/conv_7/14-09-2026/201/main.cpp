#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

static const int PRIME_LIMIT = 100000;
static const int J_LIMIT = 100;

struct Factorization
{
    vector<int> primes;
    vector<int> multiplicities;
};

struct FailureRecord
{
    int p;
    int j;
    int n;
    Factorization factorization;
    vector<int> residues;
    int omega;
    int distinct_prime_count;
    int subgroup_size;
    int max_element_order;
    bool omega_leq_subgroup;
    bool omega_eq_subgroup;
};

vector<int> build_primes(int limit);
vector<int> build_spf(int limit);
Factorization factorize(int n, const vector<int>& spf);
int mod_mul(int a, int b, int mod);
int mod_pow(int a, int e, int mod);
bool is_unit(int a, int mod);
vector<int> units_mod(int mod);
vector<int> generated_subgroup(const vector<int>& generators, int mod);
int element_order(int a, int mod);
int max_generator_order(const vector<int>& generators, int mod);
bool has_qualifying_factor(int n, int j, const vector<int>& spf);
bool proper_one_residue_divisor_exists(int n, int j, const vector<int>& divisors);
vector<int> all_divisors_from_factorization(const Factorization& f);
void print_factorization(const Factorization& f);
void print_residues(const vector<int>& residues);
void print_histogram(const unordered_map<int, long long>& histogram);
void print_failure_examples(const vector<FailureRecord>& records, int j);
void main_experiment();

int main()
{
    cout << "START EXPERIMENT 495\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 495\n";
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

int mod_mul(int a, int b, int mod)
{
    return static_cast<int>((1LL * a * b) % mod);
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

bool is_unit(int a, int mod)
{
    return gcd(a, mod) == 1;
}

vector<int> units_mod(int mod)
{
    vector<int> units;

    for (int x = 1; x < mod; ++x)
    {
        if (is_unit(x, mod))
        {
            units.push_back(x);
        }
    }

    return units;
}

vector<int> generated_subgroup(const vector<int>& generators, int mod)
{
    vector<int> subgroup;
    vector<bool> present(mod, false);

    if (mod == 1)
    {
        return {0};
    }

    present[1 % mod] = true;
    subgroup.push_back(1 % mod);

    bool changed = true;

    while (changed)
    {
        changed = false;

        const int current_size = static_cast<int>(subgroup.size());

        for (int i = 0; i < current_size; ++i)
        {
            for (int g : generators)
            {
                const int x1 = mod_mul(subgroup[i], g, mod);
                const int x2 = mod_mul(subgroup[i], mod_pow(g, -1, mod), mod);

                if (!present[x1])
                {
                    present[x1] = true;
                    subgroup.push_back(x1);
                    changed = true;
                }

                if (!present[x2])
                {
                    present[x2] = true;
                    subgroup.push_back(x2);
                    changed = true;
                }
            }
        }
    }

    sort(subgroup.begin(), subgroup.end());
    return subgroup;
}

int element_order(int a, int mod)
{
    if (!is_unit(a, mod))
    {
        return 0;
    }

    int value = 1;

    for (int order = 1; order <= mod; ++order)
    {
        value = mod_mul(value, a, mod);

        if (value == 1)
        {
            return order;
        }
    }

    return 0;
}

int max_generator_order(const vector<int>& generators, int mod)
{
    int maximum = 0;

    for (int g : generators)
    {
        maximum = max(maximum, element_order(g, mod));
    }

    return maximum;
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

bool proper_one_residue_divisor_exists(int n, int j, const vector<int>& divisors)
{
    for (int d : divisors)
    {
        if (d <= 1 || d >= n)
        {
            continue;
        }

        if (d % j == 1 % j)
        {
            return true;
        }
    }

    return false;
}

bool has_qualifying_factor(int n, int j, const vector<int>& spf)
{
    const Factorization f = factorize(n, spf);
    const vector<int> divisors = all_divisors_from_factorization(f);

    return proper_one_residue_divisor_exists(n, j, divisors);
}

void print_factorization(const Factorization& f)
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

void print_residues(const vector<int>& residues)
{
    cout << "[";

    for (size_t i = 0; i < residues.size(); ++i)
    {
        if (i != 0)
        {
            cout << ",";
        }

        cout << residues[i];
    }

    cout << "]";
}

void print_histogram(const unordered_map<int, long long>& histogram)
{
    vector<pair<int, long long>> entries(histogram.begin(), histogram.end());

    sort(
        entries.begin(),
        entries.end(),
        [](const auto& a, const auto& b)
        {
            return a.first < b.first;
        }
    );

    for (const auto& [key, value] : entries)
    {
        cout << key << ":" << value << " ";
    }

    cout << "\n";
}

void print_failure_examples(const vector<FailureRecord>& records, int j)
{
    cout << "\nFAILURE_EXAMPLES_J=" << j << "\n";

    int printed = 0;

    for (const FailureRecord& r : records)
    {
        if (r.j != j)
        {
            continue;
        }

        cout << "P=" << r.p
             << " N=" << r.n
             << " FACTORS=";

        print_factorization(r.factorization);

        cout << " RESIDUES=";
        print_residues(r.residues);

        cout << " OMEGA=" << r.omega
             << " DISTINCT=" << r.distinct_prime_count
             << " SUBGROUP_SIZE=" << r.subgroup_size
             << " MAX_ELEMENT_ORDER=" << r.max_element_order
             << " OMEGA_LEQ_SUBGROUP=" << (r.omega_leq_subgroup ? 1 : 0)
             << " OMEGA_EQ_SUBGROUP=" << (r.omega_eq_subgroup ? 1 : 0)
             << "\n";

        ++printed;

        if (printed >= 10)
        {
            break;
        }
    }

    if (printed == 0)
    {
        cout << "NONE\n";
    }
}

void main_experiment()
{
    const auto start = chrono::high_resolution_clock::now();

    const vector<int> primes = build_primes(PRIME_LIMIT);
    const vector<int> spf = build_spf(J_LIMIT * PRIME_LIMIT + 1);

    cout << "PRIME_COUNT=" << primes.size() << "\n";
    cout << "J_LIMIT=" << J_LIMIT << "\n";
    cout << "MAX_JP_PLUS_1=" << J_LIMIT * PRIME_LIMIT + 1 << "\n";
    cout << "SPF_READY=1\n";

    long long total_failures = 0;
    long long total_failure_omega = 0;
    long long total_failure_distinct = 0;
    long long total_failure_subgroup = 0;
    long long total_failure_max_order = 0;

    long long omega_leq_subgroup = 0;
    long long omega_eq_subgroup = 0;

    long long omega_lt_subgroup = 0;
    long long subgroup_lt_omega = 0;

    unordered_map<int, long long> global_omega_histogram;
    unordered_map<int, long long> global_subgroup_histogram;
    unordered_map<int, long long> global_max_order_histogram;

    unordered_map<int, long long> failure_count_by_j;
    unordered_map<int, long long> equality_count_by_j;
    unordered_map<int, long long> strict_count_by_j;

    vector<FailureRecord> failure_records;

    int total_composite = 0;
    int total_qualifying = 0;

    for (int j = 2; j <= J_LIMIT; ++j)
    {
        long long composite_count = 0;
        long long qualifying_count = 0;
        long long failures = 0;

        long long failure_omega_sum = 0;
        long long failure_distinct_sum = 0;
        long long failure_subgroup_sum = 0;
        long long failure_max_order_sum = 0;

        long long failure_omega_leq = 0;
        long long failure_omega_eq = 0;

        unordered_map<int, long long> omega_histogram;
        unordered_map<int, long long> subgroup_histogram;
        unordered_map<int, long long> max_order_histogram;

        for (int p : primes)
        {
            const long long n64 = 1LL * j * p + 1;

            if (n64 > J_LIMIT * PRIME_LIMIT + 1)
            {
                continue;
            }

            const int n = static_cast<int>(n64);

            if (spf[n] == n)
            {
                continue;
            }

            ++composite_count;
            ++total_composite;

            const Factorization f = factorize(n, spf);
            const vector<int> divisors = all_divisors_from_factorization(f);

            const bool qualifying =
                proper_one_residue_divisor_exists(n, j, divisors);

            if (qualifying)
            {
                ++qualifying_count;
                ++total_qualifying;
                continue;
            }

            ++failures;
            ++total_failures;

            const int omega =
                accumulate(
                    f.multiplicities.begin(),
                    f.multiplicities.end(),
                    0
                );

            const int distinct_prime_count =
                static_cast<int>(f.primes.size());

            vector<int> generators;

            for (int prime_factor : f.primes)
            {
                const int residue = prime_factor % j;

                if (is_unit(residue, j))
                {
                    generators.push_back(residue);
                }
                else
                {
                    cerr << "ERROR_NONUNIT_RESIDUE"
                         << " P=" << p
                         << " J=" << j
                         << " N=" << n
                         << " PRIME_FACTOR=" << prime_factor
                         << "\n";

                    return;
                }
            }

            vector<int> residues;

            for (size_t i = 0; i < f.primes.size(); ++i)
            {
                for (int e = 0; e < f.multiplicities[i]; ++e)
                {
                    residues.push_back(f.primes[i] % j);
                }
            }

            const vector<int> subgroup =
                generated_subgroup(generators, j);

            const int subgroup_size =
                static_cast<int>(subgroup.size());

            const int max_order =
                max_generator_order(generators, j);

            const bool leq =
                omega <= subgroup_size;

            const bool equal =
                omega == subgroup_size;

            if (leq)
            {
                ++omega_leq_subgroup;
                ++failure_omega_leq;
            }
            else
            {
                ++subgroup_lt_omega;
            }

            if (equal)
            {
                ++omega_eq_subgroup;
                ++failure_omega_eq;
            }

            if (omega < subgroup_size)
            {
                ++omega_lt_subgroup;
            }

            failure_omega_sum += omega;
            failure_distinct_sum += distinct_prime_count;
            failure_subgroup_sum += subgroup_size;
            failure_max_order_sum += max_order;

            ++failure_count_by_j[j];
            if (equal)
            {
                ++equality_count_by_j[j];
            }
            if (omega < subgroup_size)
            {
                ++strict_count_by_j[j];
            }

            ++omega_histogram[omega];
            ++subgroup_histogram[subgroup_size];
            ++max_order_histogram[max_order];

            ++global_omega_histogram[omega];
            ++global_subgroup_histogram[subgroup_size];
            ++global_max_order_histogram[max_order];

            FailureRecord record;
            record.p = p;
            record.j = j;
            record.n = n;
            record.factorization = f;
            record.residues = residues;
            record.omega = omega;
            record.distinct_prime_count = distinct_prime_count;
            record.subgroup_size = subgroup_size;
            record.max_element_order = max_order;
            record.omega_leq_subgroup = leq;
            record.omega_eq_subgroup = equal;

            failure_records.push_back(record);
        }

        const double qualify_percent =
            composite_count == 0
                ? 0.0
                : 100.0 * static_cast<double>(qualifying_count)
                    / static_cast<double>(composite_count);

        const double avg_failure_omega =
            failures == 0
                ? 0.0
                : static_cast<double>(failure_omega_sum)
                    / static_cast<double>(failures);

        const double avg_failure_distinct =
            failures == 0
                ? 0.0
                : static_cast<double>(failure_distinct_sum)
                    / static_cast<double>(failures);

        const double avg_failure_subgroup =
            failures == 0
                ? 0.0
                : static_cast<double>(failure_subgroup_sum)
                    / static_cast<double>(failures);

        const double avg_failure_max_order =
            failures == 0
                ? 0.0
                : static_cast<double>(failure_max_order_sum)
                    / static_cast<double>(failures);

        cout << "J=" << j
             << " COMPOSITE=" << composite_count
             << " QUALIFYING=" << qualifying_count
             << " FAILURES=" << failures
             << " QUALIFY_PERCENT=" << qualify_percent
             << " AVG_FAILURE_OMEGA=" << avg_failure_omega
             << " AVG_FAILURE_DISTINCT=" << avg_failure_distinct
             << " AVG_FAILURE_SUBGROUP=" << avg_failure_subgroup
             << " AVG_FAILURE_MAX_ORDER=" << avg_failure_max_order
             << " OMEGA_LEQ_SUBGROUP="
             << failure_omega_leq
             << " OMEGA_EQ_SUBGROUP="
             << failure_omega_eq
             << "\n";

        cout << "J=" << j << " OMEGA_HIST=";
        print_histogram(omega_histogram);

        cout << "J=" << j << " SUBGROUP_HIST=";
        print_histogram(subgroup_histogram);

        cout << "J=" << j << " MAX_ORDER_HIST=";
        print_histogram(max_order_histogram);

        print_failure_examples(failure_records, j);
    }

    cout << "\nGLOBAL_FAILURE_OMEGA_HIST=";
    print_histogram(global_omega_histogram);

    cout << "GLOBAL_FAILURE_SUBGROUP_HIST=";
    print_histogram(global_subgroup_histogram);

    cout << "GLOBAL_FAILURE_MAX_ORDER_HIST=";
    print_histogram(global_max_order_histogram);

    cout << "\nTOTAL_COMPOSITE=" << total_composite << "\n";
    cout << "TOTAL_QUALIFYING=" << total_qualifying << "\n";
    cout << "TOTAL_FAILURES=" << total_failures << "\n";

    cout << "OMEGA_LEQ_SUBGROUP=" << omega_leq_subgroup << "\n";
    cout << "OMEGA_EQ_SUBGROUP=" << omega_eq_subgroup << "\n";
    cout << "OMEGA_LT_SUBGROUP=" << omega_lt_subgroup << "\n";
    cout << "SUBGROUP_LT_OMEGA=" << subgroup_lt_omega << "\n";

    const double avg_omega =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_failure_omega)
                / static_cast<double>(total_failures);

    const double avg_distinct =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_failure_distinct)
                / static_cast<double>(total_failures);

    const double avg_subgroup =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_failure_subgroup)
                / static_cast<double>(total_failures);

    const double avg_max_order =
        total_failures == 0
            ? 0.0
            : static_cast<double>(total_failure_max_order)
                / static_cast<double>(total_failures);

    cout << "GLOBAL_AVG_FAILURE_OMEGA=" << avg_omega << "\n";
    cout << "GLOBAL_AVG_FAILURE_DISTINCT=" << avg_distinct << "\n";
    cout << "GLOBAL_AVG_FAILURE_SUBGROUP=" << avg_subgroup << "\n";
    cout << "GLOBAL_AVG_FAILURE_MAX_ORDER=" << avg_max_order << "\n";

    cout << "\nJ_FAILURE_COUNTS\n";
    for (int j = 2; j <= J_LIMIT; ++j)
    {
        const auto it = failure_count_by_j.find(j);

        if (it == failure_count_by_j.end())
        {
            continue;
        }

        const long long failures = it->second;
        const long long equalities = equality_count_by_j[j];
        const long long strict = strict_count_by_j[j];

        cout << "J=" << j
             << " FAILURES=" << failures
             << " OMEGA_EQ_SUBGROUP=" << equalities
             << " OMEGA_LT_SUBGROUP=" << strict
             << "\n";
    }

    cout << "\nSANITY_CHECKS\n";

    bool sanity_pass = true;

    for (const FailureRecord& r : failure_records)
    {
        if (!r.omega_leq_subgroup)
        {
            cout << "BOUND_FAILURE"
                 << " P=" << r.p
                 << " J=" << r.j
                 << " N=" << r.n
                 << " OMEGA=" << r.omega
                 << " SUBGROUP_SIZE=" << r.subgroup_size
                 << "\n";

            sanity_pass = false;
            break;
        }
    }

    cout << "OMEGA_LEQ_SUBGROUP_SANITY="
         << (sanity_pass ? 1 : 0)
         << "\n";

    const auto end = chrono::high_resolution_clock::now();
    const double elapsed_ms =
        chrono::duration<double, milli>(end - start).count();

    cout << "ELAPSED_TIME_MS=" << elapsed_ms << "\n";
}
