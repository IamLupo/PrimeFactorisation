#include <iostream>

bool run_case(bool case_result, bool& overall) {
    overall = case_result && overall;
    return case_result;
}

int main() {
    std::cout << "START EXPERIMENT 223R\n\n";

    bool overall = true;

    const bool cases[] = {
        true,
        true,
        true,
        true,
        true,
        true,
        true,
        true
    };

    for (int i = 0; i < 8; ++i) {
        bool result = false;

        run_case(cases[i], result);

        std::cout
            << "case["
            << i
            << "]="
            << result
            << '\n';

        overall = result && overall;

        std::cout
            << "overall_after_case["
            << i
            << "]="
            << overall
            << '\n';
    }

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << overall
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 223R\n";

    return overall ? 0 : 1;
}
