#include "rsa.h"
#include "algos.h"
#include <fmt/core.h>

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
