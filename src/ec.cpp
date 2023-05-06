#include "ec.h"

#include "generators.h"
#include "algos.h"
#include <algorithm>
#include <fmt/core.h>

/*GF(p) curves */

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

/*GF(2^n) curves */

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

/* utility functions */

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
