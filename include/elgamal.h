#pragma once
#include <string>
#include <cstdint>

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
