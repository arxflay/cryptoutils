#include "algos.h"

#include "generators.h"
#include <iostream>
#include <array>
#include <cmath>
#include <utility>
#include <algorithm>
#include <iostream>
#include <exception>
#include <fmt/core.h>

std::tuple<int_fast64_t, int_fast64_t> ExtendedGCDRec(int_fast64_t a, int_fast64_t b, int_fast64_t x1, int_fast64_t x2, int_fast64_t y1, int_fast64_t y2)
{
    if (b == 0)
        return std::tie(a, x2);
    
    int_fast64_t q = a / b;

    return ExtendedGCDRec(b, a % b, x2 - q * x1, x1, y2 - q * y1, y1);
}

std::tuple<int_fast64_t, int_fast64_t> ExtendedGCD(int_fast64_t a, int_fast64_t b)
{
    if (a < b)
        std::swap(a, b);

    return ExtendedGCDRec(a, b, 1, 0, 0, 1);
}

int_fast64_t InverseMod(int_fast64_t num, int_fast64_t mod)
{
    int_fast64_t inverse = std::get<1>(ExtendedGCD(num, mod));
    return inverse >= 0 ? inverse : mod + inverse;
}

int_fast64_t GCD(int_fast64_t num, int_fast64_t mod)
{
    return std::get<0>(ExtendedGCD(num, mod));
}

int_fast64_t ModExp(int_fast64_t num, int_fast64_t power, int_fast64_t mod)
{
    int_fast64_t prev = num % mod;
    int_fast64_t result = (power & 1) ? prev : 1;

    for (int_fast64_t i = 1; i < 64; i++)
    {
        prev = (prev * prev) % mod;
        if (power & (1LL << i))
            result = (prev * result) % mod;
    }

    return result;
}

int_fast64_t PositiveMod(int_fast64_t num, int_fast64_t mod)
{
    return (num >= 0) ? (num % mod) : ((num % mod) + mod);
}

/*returns MSB from 1*/
int_fast64_t GetMSB(int_fast64_t num)
{
    int_fast64_t firstByte = 0xFF;
    int_fast64_t index = 0;

    for (int_fast64_t byte = 7; byte >= 0; byte--)
    {
        int_fast64_t byteOfInterest = firstByte << (8 * byte);
        if (!(num & byteOfInterest))
            continue;
        
        num = (num & byteOfInterest) >> (8 * byte);
        for (int_fast64_t pos = 7; pos >= 0; pos--)
        {
            if (num & (1 << pos))
            {
                index = 8 * byte + pos + 1;
                break;
            }
        }
    }

    return index;
}

int_fast64_t ReduceGF2N(const GF2NGeneratorParameters &parameters, int_fast64_t polynomial)
{
    while(polynomial >= (1 << parameters.n))
    {
        int_fast64_t msb = GetMSB(polynomial);
        polynomial ^= parameters.ireduciblePolynomial << (msb - (parameters.n + 1));
    }

    return polynomial;
}

int_fast64_t MultiplyBinary(int_fast64_t number, int_fast64_t multiplicant)
{
    int_fast64_t multiplicantMSB = GetMSB(multiplicant);
    int_fast64_t out = 0;
    for (int_fast64_t i = 0; i < multiplicantMSB; i++)
        if (multiplicant & (1 << i))
            out ^= (number << i);
    
    return out;
}

std::string ConvertNumberToBinary(int_fast64_t number, int_fast64_t minLen)
{
    if (minLen == 0)  
        minLen = GetMSB(number);

    std::string out;

    for (int_fast64_t i = 0; i < minLen; i++)
         out = ((!!(number & (1 << i))) ? "1" : "0") + out;

    return out;
}

int_fast64_t ConvertBinaryToNumber(std::string_view binary)
{
    int_fast64_t result = 0;
    int strLen = static_cast<int>(binary.length());
    if (strLen == 0)
        return 0;

    for (int i = strLen - 1; i >= 0; i--)
        if (binary[i] == '1')
            result |= (1 << ((strLen - 1) - i));

    return result;
}

std::tuple<int_fast64_t, int_fast64_t> DoFermantFactorization(int_fast64_t number, std::string *steps)
{
    int_fast64_t t = (int_fast64_t)std::ceil(std::sqrt(number));
    
    if (steps)
        *steps = fmt::format("t0 = sqrt({}) = {}\n", number, t);
    
    if (t * t == number)
        return std::tie(t, t);

    int_fast64_t s = 0;
    int counter = 0;
    while (true)
    {
        int_fast64_t z = (t * t) - number;
        s = (int_fast64_t)std::floor(std::sqrt(z));
        
        if  (steps)
            *steps = fmt::format("z = {}^2 - {} = {}\n", t, number, z);
        
        if (s * s == z)
            break;

        t = t + 1;

        if (steps)
            *steps = fmt::format("\nt{} = {} + 1 = {}\n", ++counter, t - 1, t);
    }
    
    return std::make_tuple<int_fast64_t, int_fast64_t>(t + s, t - s);
}

std::tuple<int_fast64_t, int_fast64_t> DoRhoFactorization(int_fast64_t number, int_fast64_t rndNumber, std::string *steps)
{
    int_fast64_t x1 = ((rndNumber * rndNumber) - 1) % number;
    int_fast64_t x2 = ((x1 * x1) - 1) % number;
    int_fast64_t result = 0;

    int counter = 3;

    while(true)
    {
        if (steps)
            *steps += fmt::format("i = {}\nx{} = {}\nx{} = {}\n", (counter / 2) - 1, counter - 2, x1, counter - 1, x2);

        int_fast64_t xDelta = std::abs(x2 - x1);
        result = GCD(xDelta, number);

        if (steps)
        {
            *steps += fmt::format("|x{} - x{}| = |{} - {}| = |{}|\n", counter - 1, counter - 2, x2, x1, xDelta);
            *steps += fmt::format("gcd({}, {}) = {}\n\n", xDelta, number, result);
        }

        if (result > 1)
            break;

        x1 = ((x1 * x1) - 1) % number;
        int_fast64_t x3 = ((x2 * x2) - 1) % number;
        x2 = ((x3 * x3) - 1) % number;
        
        counter += 2;
    }


    return std::make_tuple(result, number / result);
}

LehmanPeraltResult LehmanPeralt(const std::vector<int_fast64_t> &numbers, int_fast64_t examinedNumber, std::string *steps)
{
    if (numbers.size() == 0)
        return {};

    LehmanPeraltResult out;

    int_fast64_t power = (examinedNumber - 1) / 2;
    
    for (int_fast64_t i = 0; i < static_cast<int_fast64_t>(numbers.size()); i++)
    {
        int_fast64_t result = ModExp(numbers[i], power, examinedNumber);

        if (result == 1)
            out.flags |= LehmanPeraltFlags::ONE_OCCURED;
        else if (result == examinedNumber - 1)
        {
            out.flags |= LehmanPeraltFlags::MINUS_ONE_OCCURED;
            result = -1;
        }
        else
            out.flags |= LehmanPeraltFlags::NOT_PRIME;

        if (steps)
            *steps += fmt::format("d{} = {} ^ {} mod {} = {}\n", i + 1, numbers[i], power, examinedNumber, result);
    }

    if (!(out.flags & LehmanPeraltFlags::NOT_PRIME))
        out.chance = 1.0 - std::pow(2.0, -(double)numbers.size());

    return out;
}


