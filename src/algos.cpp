#include "algos.h"

#include <iostream>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>
#include <algorithm>
#include <iostream>
#include <random>
#include <exception>
#include <fmt/core.h>
#include <set>

/* naive implementation */
uint_least64_t CalculateHammingWeight(const char *data, int_fast64_t len)
{
    uint_least64_t weight = 0;
    for (int_fast64_t i = 0; i < len; i++)
        for (uint_fast8_t j = 1; j < sizeof(char); j++)
            weight += (data[i] & (1 << j)) >> j;
    
    return weight;
}

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

int_fast64_t ReduceGF2N(const GF2NGeneratorParameters &parameters, int_fast64_t polynomial)
{
    while(polynomial >= (1 << parameters.n))
    {
        int_fast64_t msb = GetMSB(polynomial);
        polynomial ^= parameters.ireduciblePolynomial << (msb - (parameters.n + 1));
    }

    return polynomial;
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

ElGamalData ElGamalEncrypt(const ElGamalPublicKey &pubKey, const ElGamalPrivateKey &privKey, int_fast64_t message, std::string *steps)
{
    int_fast64_t sharedSecret = ModExp(pubKey.y, privKey.k, pubKey.p);
    int_fast64_t y = ModExp(pubKey.q, privKey.k, pubKey.p);    
    int_fast64_t encrypted = (message * sharedSecret) % pubKey.p;

    if (steps)
    {
        std::string &str = *steps;
        str = fmt::format("sharedSecret = {} ^ {} mod {} = {}\n", pubKey.y, privKey.k, pubKey.p, sharedSecret);
        str += fmt::format("y = {} ^ {} mod {} = {}\n", pubKey.q, privKey.k, pubKey.p, y);
        str += fmt::format("encryptedMsg = ({} * {}) mod {} = {}", message, sharedSecret, pubKey.p, encrypted);
    }

    return { y, encrypted };
}

int_fast64_t ElGamalDecrypt(const ElGamalData &data, const ElGamalPrivateKey &privKey, std::string *steps)
{
    int_fast64_t sharedSecret = ModExp(data.y, privKey.k, privKey.p);
    int_fast64_t inverse = InverseMod(sharedSecret, privKey.p);
    int_fast64_t decrypted = (inverse * data.encData) % privKey.p;

    if (steps)
    {
        std::string &str = *steps;
        str = fmt::format("sharedSecret = {} ^ {} mod {} = {}\n", data.y, privKey.k, privKey.p, sharedSecret);
        str += fmt::format("sharedSecret^(-1) = {}^(-1) = {}\n", sharedSecret, inverse);
        str += fmt::format("decryptedMsg = ({} * {}) mod {} = {}", inverse, data.encData, privKey.p, decrypted);
    }

    return decrypted;
}

ElGamalPublicKey ElGamalDerivePublicKey(const ElGamalPrivateKey &privKey, std::string *steps)
{
    int_fast64_t publicPart = ModExp(privKey.q, privKey.k, privKey.p);
    if (steps)
        *steps = fmt::format("y = {}^{} mod {} = {}", privKey.q, privKey.k, privKey.p, publicPart);

    return ElGamalPublicKey{ privKey.p, privKey.q, publicPart };
}

bool ECAlignsOn(const ECCurve &curve, const ECPoint &p)
{
    return (p.y * p.y) % curve.p == PositiveMod((p.x * p.x * p.x) + curve.a * p.x + curve.b, curve.p);
}

ECPoint ECDouble(const ECCurve &curve, const ECPoint &p, std::string *steps)
{
    if (p.y == 0)
        throw std::invalid_argument("p.y must be highier than 0");

    //s = slope
    
    int_fast64_t s = ((3 * (p.x * p.x) + curve.a) * InverseMod(PositiveMod(2 * p.y, curve.p), curve.p)) % curve.p;
    
    ECPoint newPoint;
    newPoint.x = PositiveMod((s * s) - 2 * p.x, curve.p);
    newPoint.y = PositiveMod((s * (p.x - newPoint.x) - p.y), curve.p);

    if (steps)
    {
        std::string &str = *steps;
        str = fmt::format("s = (3 * {}^2 + {}) * ((2 * {})^-1) mod {} = {}\n", p.x, curve.a, p.y, curve.p, s);
        str += fmt::format("xr = {}^2 - 2 * {} mod {} = {}\n", s, p.x, curve.p, newPoint.x);
        str += fmt::format("yr = {} * ({} - {}) - {} mod {} = {}", s, p.x, newPoint.x, p.y, curve.p, newPoint.y);
    }

    return newPoint;
}

ECPoint ECSum(const ECCurve &curve, const ECPoint &p, const ECPoint &q, std::string *steps)
{
    /*Point doubling*/
    if (p == q)
        return ECDouble(curve, p, steps);
    //s = slope
    int_fast64_t s = PositiveMod(((p.y - q.y) * InverseMod(PositiveMod(p.x - q.x, curve.p), curve.p)), curve.p);
    if (s == 0)
        throw std::runtime_error("s equals 0, result point is point at infinity, (0, 0)");

    ECPoint newPoint;
    
    newPoint.x = PositiveMod(((s * s) - p.x - q.x), curve.p);

    /* 
     * y_-r(inverse y) = y_p + z
     * w(width) = x_r - x_p
     * z(height) = w * s = (x_r - x_p); 
     * y_-r(inverse y) = y_p + s (x_r - x_p)
     * y_r = -y_p + s * (x_p - x_r)
     */
    newPoint.y = PositiveMod((s * (p.x - newPoint.x) - p.y), curve.p);

    if (steps)
    {
        std::string &str = *steps;
        str = fmt::format("s = ({} - {}) * (({} - {})^-1) mod {} = {}\n", p.y, q.y, p.x, q.x, curve.p, s);
        str += fmt::format("xr = {}^2 - {} - {} mod {} = {}\n", s, p.x, q.x, curve.p, newPoint.x);
        str += fmt::format("yr = {} * ({} - {}) - {} mod {} = {}", s, p.x, newPoint.x, p.y, curve.p, newPoint.y);
    }

    return newPoint;
}

bool ECAlignsOnGF2N(const ECCurve &curve, const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps)
{
    int_fast64_t xPower      = p.x     == 0 ? 0 : ECPointGetPowerGF2N(generatorPoints, p.x); 
    int_fast64_t yPower      = p.y     == 0 ? 0 : ECPointGetPowerGF2N(generatorPoints, p.y);
    int_fast64_t curveAPower = curve.a == 0 ? 0 : ECPointGetPowerGF2N(generatorPoints, curve.a);

    int_fast64_t y2 = p.y == 0 ? 0 : generatorPoints.at((yPower * 2) % curve.p); //power multiplication
    int_fast64_t xy = (!p.x || !p.y) ? 0 : generatorPoints.at((xPower + yPower) % curve.p); //sum of powers
                                                                              
    int_fast64_t x3 = p.x == 0 ? 0 : generatorPoints.at((xPower * 3) % curve.p);
    int_fast64_t x2 = (!curve.a || !p.x) ? 0 : generatorPoints.at((curveAPower + xPower * 2) % curve.p);

    int_fast64_t leftSide = y2 ^ xy;
    int_fast64_t rightSide = x3 ^ x2 ^ curve.b;

    if (steps)
    {
        std::string yStr;
        std::string xStr;
        ECPointToStrGF2N(p, generatorPoints, xStr, yStr); 

        std::string curveAStr;
        std::string curveBStr;
        ECCurveToStrGF2N(curve, generatorPoints, curveAStr, curveBStr); 

        *steps = fmt::format("({})^2 + {} * {} = ({})^3 + {} * ({})^2 + {}\n", yStr, xStr, yStr, xStr, curveAStr, xStr, curveBStr);

        *steps += fmt::format("{} + {} = {} + {} + {}\n"
                , ConvertNumberToBinary(y2, parameters.n)
                , ConvertNumberToBinary(xy, parameters.n)
                , ConvertNumberToBinary(x3, parameters.n)
                , ConvertNumberToBinary(x2, parameters.n)
                , ConvertNumberToBinary(curve.b, parameters.n));
        *steps += fmt::format("{} = {}", ConvertNumberToBinary(leftSide, parameters.n), ConvertNumberToBinary(rightSide, parameters.n));
    }
    
    return leftSide == rightSide;
}

ECPoint ECDoubleGF2N(const ECCurve &curve, const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps)
{
    if (p.x == 0)
        throw std::runtime_error("2P = O, O = point at infinity");
   
    ECPoint newPoint;

    int_fast64_t s = 0;
    int_fast64_t sNum = p.x;
    int_fast64_t sNumSquare = 0;
    int_fast64_t xPower = ECPointGetPowerGF2N(generatorPoints, p.x);
    int_fast64_t xSquare = generatorPoints.at((xPower * 2) % curve.p); 

    if (p.y)
    {
        int_fast64_t yPower = ECPointGetPowerGF2N(generatorPoints, p.y);
        int_fast64_t pyDivPxGN = PositiveMod(yPower - xPower, curve.p);
        sNum ^= generatorPoints.at(pyDivPxGN);
    }

    if (sNum)
    {
        s = ECPointGetPowerGF2N(generatorPoints, sNum);
        sNumSquare = generatorPoints.at((s * 2) % curve.p);
    }

    int_fast64_t xr = sNumSquare ^ sNum ^ curve.a;
    int_fast64_t sOneSum = sNum ^ 1;
    int_fast64_t xrMultS = 0;
    if (xr && sOneSum)
        xrMultS = generatorPoints.at(PositiveMod(ECPointGetPowerGF2N(generatorPoints, xr) + ECPointGetPowerGF2N(generatorPoints, sOneSum), curve.p));

    int_fast64_t yr = xSquare ^ xrMultS;
    
    newPoint.x = xr;
    newPoint.y = yr;

    if (steps)
    {
        std::string pyStr;
        std::string pxStr;
        ECPointToStrGF2N(p, generatorPoints, pxStr, pyStr); 
        
        std::string curveAStr;
        std::string curveBStr;
        ECCurveToStrGF2N(curve, generatorPoints, curveAStr, curveBStr);

        std::string sStr = sNum ? fmt::format("g^{}", s) : "0";
        std::string newPointXStr;
        std::string newPointYStr;
        
        ECPointToStrGF2N(newPoint, generatorPoints, newPointXStr, newPointYStr);
        *steps = fmt::format("s = {} + ({} / {}) = {}\n", pxStr, pyStr, pxStr, sStr);
        *steps += fmt::format("xr = {} + {} + {} = {} = {}\n"
            , ConvertNumberToBinary(sNumSquare, parameters.n)
            , ConvertNumberToBinary(sNum, parameters.n)
            , ConvertNumberToBinary(curve.a, parameters.n)
            , ConvertNumberToBinary(xr, parameters.n)
            , newPointXStr
        );

        *steps += fmt::format("yr = ({})^2 + {} * ({} + 1) = "
            , pxStr
            , newPointXStr
            , sStr
        );

        *steps += fmt::format("{} + {} * {} = "
            , ConvertNumberToBinary(xSquare, parameters.n)
            , newPointXStr
            , sOneSum ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, sOneSum)) : "0" 
        );

        *steps += fmt::format("{} + {} = {} = {}"
            , ConvertNumberToBinary(xSquare, parameters.n)
            , ConvertNumberToBinary(xrMultS, parameters.n)
            , ConvertNumberToBinary(yr, parameters.n)
            , newPointYStr
        );

    }
    
    
    return newPoint;
}

ECPoint ECSumGF2N(const ECCurve &curve, const ECPoint &p, const ECPoint &q, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps)
{
    if (p == q)
        return ECDoubleGF2N(curve, p, generatorPoints, parameters, steps);
   
    ECPoint pInverse;
    pInverse.x = p.x;
    pInverse.y = p.x ^ p.y;

    if (pInverse == q) 
        throw std::runtime_error("Q equals to -P, P + -P = O, O = point at infinity");

    /*slope*/
    
    ECPoint newPoint;

    int_fast64_t pySum = p.y ^ q.y;
    int_fast64_t sNum = 0;
    int_fast64_t sSquareNum = 0;

    int_fast64_t s = PositiveMod(ECPointGetPowerGF2N(generatorPoints, pySum) - ECPointGetPowerGF2N(generatorPoints, p.x ^ q.x), curve.p);
    if (pySum)
    {
        sNum = generatorPoints.at(s);
        sSquareNum = generatorPoints.at((s * 2) % curve.p);
    }

    int_fast64_t xr = sSquareNum ^ sNum ^ p.x ^ q.x ^ curve.a;
    int_fast64_t pxXrSumNum = p.x ^ xr;
    int_fast64_t sMultPxXrSumNum = 0;
    if (sNum && pxXrSumNum)
        sMultPxXrSumNum = generatorPoints.at((s + ECPointGetPowerGF2N(generatorPoints, pxXrSumNum)) % curve.p);
        
    int_fast64_t yr = sMultPxXrSumNum ^ xr ^ p.y;
    
    newPoint.x = xr;
    newPoint.y = yr;
    
    if (steps)
    {
        std::string pyStr;
        std::string pxStr;
        ECPointToStrGF2N(p, generatorPoints, pxStr, pyStr); 
        
        std::string qyStr;
        std::string qxStr;
        ECPointToStrGF2N(q, generatorPoints, qxStr, qyStr); 

        std::string curveAStr;
        std::string curveBStr;
        ECCurveToStrGF2N(curve, generatorPoints, curveAStr, curveBStr);

        std::string sStr = sNum ? fmt::format("g^{}", s) : "0";
        
        std::string newPointXStr;
        std::string newPointYStr;
        ECPointToStrGF2N(newPoint, generatorPoints, newPointXStr, newPointYStr);

        *steps = fmt::format("s = ({} + {}) / ({} + {}) = {}\n", pyStr, qyStr, pxStr, qxStr, sStr);
        *steps += fmt::format("xr = {} + {} + {} + {} + {} = {} = {}\n"
                , ConvertNumberToBinary(sSquareNum, parameters.n)
                , ConvertNumberToBinary(sNum, parameters.n)
                , ConvertNumberToBinary(p.x, parameters.n)
                , ConvertNumberToBinary(q.x, parameters.n)
                , ConvertNumberToBinary(curve.a, parameters.n)
                , ConvertNumberToBinary(xr, parameters.n)
                , newPointXStr
        );

        *steps += fmt::format("yr = {} * ({} + {}) + {} + {} = " 
                , sStr
                , ConvertNumberToBinary(p.x, parameters.n)
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(p.y, parameters.n)
        );

        *steps += fmt::format("{} * {} + {} + {} = " 
                , sStr
                , ConvertNumberToBinary(pxXrSumNum, parameters.n)
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(p.y, parameters.n)
        );

        *steps += fmt::format("{} * {} + {} + {} = " 
                , sStr
                , ConvertNumberToBinary(pxXrSumNum, parameters.n)
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(p.y, parameters.n)
        );

        *steps += fmt::format("{} * {} + {} + {} = " 
                , sStr
                , pxXrSumNum ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, pxXrSumNum)) : "0"
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(p.y, parameters.n)
        );

        *steps += fmt::format("{} + {} + {} = {} = {}" 
                , ConvertNumberToBinary(sMultPxXrSumNum, parameters.n) 
                , ConvertNumberToBinary(xr, parameters.n)
                , ConvertNumberToBinary(p.y, parameters.n)
                , ConvertNumberToBinary(yr, parameters.n)
                , newPointYStr
        );

    }
    
    return newPoint;
}

ECPoint ECMultiplyGF2N(const ECCurve &curve, const ECPoint &p, int_fast64_t scalar, const std::vector<int_fast64_t> &generatorPoints, const GF2NGeneratorParameters &parameters, std::string *steps)
{
    ECPoint newPoint = p;
    std::string tempSteps;

    for (int_fast64_t i = 1; i < scalar; i++)
    {
        newPoint = ECSumGF2N(curve, newPoint, p, generatorPoints, parameters, &tempSteps);
        
        if (steps)
        {
            if (i + 1 == scalar)
                *steps += fmt::format("step {}:\n{}", i, tempSteps);
            else
                *steps += fmt::format("step {}:\n{}\n\n", i, tempSteps);    
        }
    }
    return newPoint;
}

int_fast64_t ECPointGetPowerGF2N(const std::vector<int_fast64_t> &generatorPoints, int_fast64_t num)
{
    auto pointIt = std::find(generatorPoints.cbegin(), generatorPoints.cend(), num);
    return std::distance(generatorPoints.cbegin(), pointIt);
}

void ECPointToStrGF2N(const ECPoint &p, const std::vector<int_fast64_t> &generatorPoints, std::string &xStr, std::string &yStr)
{
    xStr = p.x ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, p.x)) : "0";
    yStr = p.y ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, p.y)) : "0";
}

void ECCurveToStrGF2N(const ECCurve &curve, const std::vector<int_fast64_t> &generatorPoints, std::string &curveA, std::string &curveB)
{
    curveA = curve.a ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, curve.a)) : "0";
    curveB = curve.b ? fmt::format("g^{}", ECPointGetPowerGF2N(generatorPoints, curve.b)) : "0";
}

int_fast64_t RsaEncrypt(const RsaPublicKey &pubKey, int_fast64_t message, std::string *steps) 
{
    int_fast64_t encrypted = ModExp(message, pubKey.e, pubKey.n);

    if (steps)
        *steps = fmt::format("{}^{} mod {} = {}", message, pubKey.e, pubKey.n, encrypted);

    return encrypted;
}

int_fast64_t RsaDecrypt(const RsaPrivateKey &privKey, int_fast64_t message, std::string *steps)
{
    int_fast64_t decrypted = ModExp(message, privKey.d, privKey.n);

    if (steps)
        *steps = fmt::format("{}^{} mod {} = {}", message, privKey.d, privKey.n, decrypted);

    return decrypted;
}

RsaPrivateKey RsaDerivePrivateKeyFromModule(const RsaPublicKey &pubKey, std::string *steps)
{
    RsaPrivateKey privKey;
    privKey.n = pubKey.n;

    auto tuple = DoFermantFactorization(pubKey.n, steps);
    int_fast64_t p = std::get<0>(tuple);
    int_fast64_t q = std::get<1>(tuple);

    int_fast64_t phi = (p - 1) * (q - 1);

    privKey.d = InverseMod(pubKey.e, phi);

    if (steps)
    {
        *steps += fmt::format("p = {}\n", p);
        *steps += fmt::format("q = {}\n", q);
        *steps += fmt::format("phi = (p - 1) * (q - 1) = {} * {} = {}\n", p - 1, q - 1, phi);        
        *steps += fmt::format("d = {}^-1 mod {} = {}", pubKey.e, phi, privKey.d);
    }

    return privKey;
}

std::tuple<RsaPrivateKey, RsaPublicKey> RsaDeriveKeysFromPublicExponent(int_fast64_t p, int_fast64_t q, int_fast64_t e, std::string *steps) 
{
    RsaPrivateKey privateKey;
    RsaPublicKey publicKey;

    int_fast64_t n = p * q;
    int_fast64_t phi = (p - 1) * (q - 1);

    privateKey.n = n;
    privateKey.d = InverseMod(e, phi);

    publicKey.n = n;
    publicKey.e = e;

    if (steps)
    {
        *steps += fmt::format("p = {}\n", p);
        *steps += fmt::format("q = {}\n", q);
        *steps += fmt::format("phi = (p - 1) * (q - 1) = {} * {} = {}\n", p - 1, q - 1, phi);        
        *steps += fmt::format("d = {}^-1 mod {} = {}\n", e, phi, privateKey.d);
    }

    return std::make_tuple(privateKey, publicKey);
}

std::tuple<RsaPrivateKey, RsaPublicKey> RsaDeriveKeysFromPrivateExponent(int_fast64_t p, int_fast64_t q, int_fast64_t d, std::string *steps) 
{
    RsaPrivateKey privateKey;
    RsaPublicKey publicKey;

    int_fast64_t n = p * q;
    int_fast64_t phi = (p - 1) * (q - 1);

    privateKey.n = n;
    privateKey.d = d;

    publicKey.n = n;
    publicKey.e = InverseMod(d, phi);
    
    if (steps)
    {
        *steps += fmt::format("p = {}\n", p);
        *steps += fmt::format("q = {}\n", q);
        *steps += fmt::format("phi = (p - 1) * (q - 1) = {} * {} = {}\n", p - 1, q - 1, phi);        
        *steps += fmt::format("e = {}^-1 mod {} = {}\n", d, phi, publicKey.e);
    }

    return std::make_tuple(privateKey, publicKey);
}

ShamirSubject DoLagrangeInterpolation(const ShamirParameters &paramaters, int_fast64_t x, std::string *steps = nullptr)
{
    x = PositiveMod(x, paramaters.p);
    int_fast64_t outY = 0;

    int counter = 0;

    const std::vector<ShamirSubject> &subjects = paramaters.subjects;
    for (auto mainIt = subjects.begin(); mainIt != subjects.end(); mainIt++)
    {
        const ShamirSubject &iSubject = *mainIt;

        int_fast64_t numerator = 1;
        int_fast64_t denominator = 1;

        for (auto secondaryIt = subjects.begin(); secondaryIt != subjects.end(); secondaryIt++)
        {
            if (secondaryIt == mainIt)
                continue;

            const ShamirSubject &jSubject = *secondaryIt;
            numerator *= (x - jSubject.x);
            denominator *= (iSubject.x - jSubject.x);

            numerator %= paramaters.p;
            denominator %= paramaters.p;
            
            if (denominator == 0)
                throw std::runtime_error(fmt::format("{} - {} = 0, denominator equals zero", iSubject.x, jSubject.x));
        }

        if (steps)
            *steps += fmt::format("c{} = {} * {}^(-1) = ", counter, numerator, denominator);

        if(denominator < 0)
        {
            denominator *= -1;
            numerator *= -1;
        }

        denominator = InverseMod(denominator, paramaters.p);
        int_fast64_t ci = (numerator * denominator);
        int_fast64_t yici = (iSubject.y * ci);
        int_fast64_t ySum = outY + yici;

        if (steps)
        {
            *steps += fmt::format("{} * {} = {}\n", numerator, denominator, ci);
            *steps += fmt::format("y{}*c{} = {} * {} = {}\n", counter, counter, iSubject.y, ci, yici);
            if (mainIt == subjects.begin())
                *steps += fmt::format("y = {}\n\n", ySum);
            else
                *steps += fmt::format("y = {} + {} = {}\n\n", outY, yici, ySum);
        }

        outY = ySum;
        counter++;
    }

    int_fast64_t finalOutY = PositiveMod(outY, paramaters.p);
    
    if (steps)
        *steps += fmt::format("y = {} mod {} = {}", outY, paramaters.p, finalOutY);   

    return ShamirSubject{ x, PositiveMod(outY, paramaters.p) };
}

std::vector<ShamirSubject> GetShamirSubjects(const ShamirParameters &paramaters, int_fast64_t n, std::string *steps)
{
    std::vector<int_fast64_t> subjectsIndicies(paramaters.p - 1);
    std::iota(subjectsIndicies.begin(), subjectsIndicies.end(), 1);
    std::mt19937 generator(std::random_device{}());
    std::shuffle(subjectsIndicies.begin(), subjectsIndicies.end(), generator);

    std::vector<ShamirSubject> out;
    const std::vector<ShamirSubject> &subjects = paramaters.subjects;
    for (int_fast64_t i = 0; i < n; i++)
    {
        int_fast64_t subjectX = subjectsIndicies.at(i);

        auto it = std::find_if(subjects.begin(), subjects.end(), [&subjects, subjectX](const ShamirSubject &subject) { return subject.x == subjectX; });
        if (it != subjects.end())
        {
            out.push_back(*it);
            continue;
        }

        if (steps)
            *steps += fmt::format("<steps for x = {}>\n", subjectX);

        out.push_back(DoLagrangeInterpolation(paramaters, subjectX, steps));
        if (steps && i + 1 != n)
            *steps += "\n\n";
    }

    return out;
}

ShamirSubject DoShamirReconstruction(const ShamirParameters &paramaters, std::string *steps)
{
    return DoLagrangeInterpolation(paramaters, 0, steps);
}
