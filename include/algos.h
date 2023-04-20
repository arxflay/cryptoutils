#pragma once
#include <tuple>
#include <string>
#include <vector>

/*Additional algorithms*/
uint_least64_t CalculateHammingWeight(const char *data, size_t len);

/*Common algorithms */

std::tuple<int_fast64_t, int_fast64_t> ExtendedGCD(int_fast64_t a, int_fast64_t b);

int_fast64_t InverseMod(int_fast64_t num, int_fast64_t mod);
int_fast64_t GCD(int_fast64_t num, int_fast64_t mod);

int_fast64_t ModExp(int_fast64_t num, int_fast64_t power, int_fast64_t mod);
int_fast64_t Mod(int_fast64_t num, int_fast64_t mod);

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

/* ElGamal */

struct ElGamalPublicKey
{
    int_fast64_t p; //prime
    int_fast64_t q; //generator
    int_fast64_t y; //pubkey part
};

struct ElGamalPrivateKey
{
    int_fast64_t p; //prime
    int_fast64_t q; //generator  
    int_fast64_t k; //privkey
};

struct ElGamalData
{
    int_fast64_t y; //pubkey part
    int_fast64_t encData; //data
};

ElGamalData ElGamalEncrypt(const ElGamalPublicKey &pubKey, const ElGamalPrivateKey &privKey, int_fast64_t message, std::string *steps = nullptr);
int_fast64_t ElGamalDecrypt(const ElGamalData &data, const ElGamalPrivateKey &privKey, std::string *steps = nullptr);
ElGamalPublicKey ElGamalDerivePublicKey(const ElGamalPrivateKey &privKey, std::string *steps = nullptr);

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

bool ECAlignsOn(const ECCurve &curve, const ECPoint &p);
ECPoint ECDoubling(const ECCurve &curve, const ECPoint &p, std::string *steps = nullptr);
ECPoint ECSum(const ECCurve &curve, const ECPoint &p, const ECPoint &q, std::string *steps = nullptr);

