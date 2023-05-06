#pragma once
#include <string>
#include <cstdint>
#include <vector>
/* Generator things */

struct GF2NGeneratorParameters
{
    int_fast64_t polynomial;
    int_fast64_t ireduciblePolynomial;
    int_fast64_t n;
};

bool IsGeneratorGFP(int_fast64_t generator, int_fast64_t mod, std::string *steps = nullptr);
bool IsGeneratorGF2N(const GF2NGeneratorParameters &parameters, std::string *steps = nullptr);
std::vector<int_fast64_t> GetGF2NGeneratorElements(const GF2NGeneratorParameters &parameters);
