#pragma once
#include <tuple>
#include <string>
#include <vector>

/*Common algorithms */

std::tuple<int_fast64_t, int_fast64_t> ExtendedGCD(int_fast64_t a, int_fast64_t b);

int_fast64_t InverseMod(int_fast64_t num, int_fast64_t mod);
int_fast64_t GCD(int_fast64_t num, int_fast64_t mod);

int_fast64_t ModExp(int_fast64_t num, int_fast64_t power, int_fast64_t mod);
int_fast64_t Mod(int_fast64_t num, int_fast64_t mod);
int_fast64_t PositiveMod(int_fast64_t num, int_fast64_t mod);
int_fast64_t GetMSB(int_fast64_t num);

struct GF2NGeneratorParameters;
int_fast64_t ReduceGF2N(const GF2NGeneratorParameters &parameters, int_fast64_t polynomial);

int_fast64_t MultiplyBinary(int_fast64_t number, int_fast64_t multiplicant);
std::string ConvertNumberToBinary(int_fast64_t number, int_fast64_t minLen);
int_fast64_t ConvertBinaryToNumber(std::string_view binary);

/*Factorization methods */

std::tuple<int_fast64_t, int_fast64_t> DoFermantFactorization(int_fast64_t number, std::string *steps = nullptr);
std::tuple<int_fast64_t, int_fast64_t> DoRhoFactorization(int_fast64_t number, int_fast64_t rndNumber, std::string *steps = nullptr);

/*Lehman peralt */

enum LehmanPeraltFlags
{
    ONE_OCCURED = 0x01,
    MINUS_ONE_OCCURED = 0x02,
    NOT_PRIME = 0x04
};

struct LehmanPeraltResult
{
    double chance;
    uint8_t flags;
};

LehmanPeraltResult LehmanPeralt(const std::vector<int_fast64_t> &numbers, int_fast64_t examinedNumber, std::string *steps = nullptr);



