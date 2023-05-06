#pragma once
#include <vector>
#include <cstdint>
#include <string>


struct GF2NGeneratorParameters;

/* Basic information for all EC */

struct ECPoint
{
    int_fast64_t x;
    int_fast64_t y;
    bool operator==(const ECPoint &p) const noexcept { return x == p.x && y == p.y; }
    bool operator!=(const ECPoint &p) const noexcept { return x != p.x || y != p.y; }
};

struct ECCurve
{
    int_fast64_t a;
    int_fast64_t b;
    int_fast64_t p;
};

/*EC y^2 mod p = (x^3 + ax + b) mod p*/

bool ECAlignsOn(const ECCurve &curve, const ECPoint &p);
ECPoint ECDouble(const ECCurve &curve, const ECPoint &p, std::string *steps = nullptr);
ECPoint ECSum(const ECCurve &curve, const ECPoint &p, const ECPoint &q, std::string *steps = nullptr);

/*EC GF2N y^2 + xy = x^3 + ax^2 + b*/

bool ECAlignsOnGF2N(const ECCurve &curve, const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps = nullptr);
ECPoint ECDoubleGF2N(const ECCurve &curve, const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps = nullptr);
ECPoint ECSumGF2N(const ECCurve &curve, const ECPoint &p, const ECPoint &q, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps = nullptr);
ECPoint ECMultiplyGF2N(const ECCurve &curve, const ECPoint &p, int_fast64_t scalar, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps = nullptr);

/*Formatting functions for GF2N Curve*/

void ECPointToStrGF2N(const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, std::string &xStr, std::string &yStr);
void ECCurveToStrGF2N(const ECCurve &curve, const std::vector<int_fast64_t> &generatorPoints, std::string &curveA, std::string &curveB);
int_fast64_t ECPointGetPowerGF2N(const std::vector<int_fast64_t> &generatorPoints, int_fast64_t num);

