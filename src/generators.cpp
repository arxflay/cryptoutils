#include "generators.h"
#include "algos.h"
#include <set>
#include <fmt/core.h>
#include <cmath>

bool IsGeneratorGFP(int_fast64_t generator, int_fast64_t mod, std::string *steps)
{
    if (mod <= generator)
        throw std::runtime_error("mod <= generator");

    if (steps)
    {
        *steps = fmt::format("{}^0 = {}\n", generator, 1);
        *steps += fmt::format("{}^1 = {}\n", generator, generator);
    }

    std::set<int_fast64_t> numbers;
    numbers.insert(1);
    numbers.insert(generator);

    int_fast64_t number = generator;

    for (int_fast64_t i = 2; i < (mod - 1); i++)
    {
        number = (number * generator) % mod;
        if (steps)
            *steps += fmt::format("{}^{} = {}\n", generator, i, number);

        if (numbers.insert(number).second == false)
            return false;
    }

    return true;
}

bool IsGeneratorGF2N(const GF2NGeneratorParameters &parameters, std::string *steps)
{
    if (parameters.n <=1)
        return true;
    else if (parameters.ireduciblePolynomial <= parameters.polynomial)
        throw std::runtime_error("ireduciblePolynomial <= polynomial");

    int_fast64_t number = parameters.polynomial;
    int_fast64_t mod = (int_fast64_t)std::pow(2, parameters.n);

    std::string polynomialStr = ConvertNumberToBinary(parameters.polynomial, parameters.n);
    
    if (steps)
    {
        *steps = fmt::format("{}^0 = {}\n", polynomialStr, ConvertNumberToBinary(1, parameters.n));
        *steps += fmt::format("{}^1 = {}\n", polynomialStr, polynomialStr);
    }

    std::set<int_fast64_t> numbers;
    numbers.insert(1);
    numbers.insert(parameters.polynomial);

    //start from 2 because first two numbers are already known
    //mod - 1 because first number ^0 is already known, otherwise sequence will repeat
    for (int_fast64_t i = 2; i < mod - 1; i++)
    {
        number = MultiplyBinary(number, parameters.polynomial);
        number = ReduceGF2N(parameters, number);

        if (steps)
            *steps += fmt::format("{}^{} = {}\n", polynomialStr, i, ConvertNumberToBinary(number, parameters.n));
        
        if (numbers.insert(number).second == false)
            return false;
    }

    return true;
}

std::vector<int_fast64_t> GetGF2NGeneratorElements(const GF2NGeneratorParameters &parameters)
{
    int_fast64_t number = parameters.polynomial;
    int_fast64_t mod = (int_fast64_t)std::pow(2, parameters.n);

    std::vector<int_fast64_t> numbers;
    numbers.push_back(1);
    numbers.push_back(parameters.polynomial);

    //start from 2 because first two numbers are already known
    //mod - 1 because first number ^0 is already known, otherwise sequence will repeat
    for (int_fast64_t i = 2; i < mod - 1; i++)
    {
        number = MultiplyBinary(number, parameters.polynomial);
        while(number >= (1 << parameters.n))
        {
            int_fast64_t msb = GetMSB(number);
            number ^= parameters.ireduciblePolynomial << (msb - (parameters.n + 1));
        }        
        numbers.push_back(number);
    }

    return numbers;
}

