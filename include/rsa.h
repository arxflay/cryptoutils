#pragma once

#include <tuple>
#include <cstdint>
#include <string>

struct RsaPublicKey 
{
    int_fast64_t n; // Module
    int_fast64_t e; // Encryption exponent
};

struct RsaPrivateKey
{
    int_fast64_t n; // Module
    int_fast64_t d; // Decryption exponent
};

/* Encrypt message with public key. */
int_fast64_t RsaEncrypt(const RsaPublicKey &pubKey, int_fast64_t message, std::string *steps = nullptr);

/* Decrypt message with private key. */
int_fast64_t RsaDecrypt(const RsaPrivateKey &privKey, int_fast64_t message, std::string *steps = nullptr);

/* Derive private key from public key. */
RsaPrivateKey RsaDerivePrivateKeyFromModule(const RsaPublicKey &pubKey, std::string *steps = nullptr);

/* Derive private and public key from p q e. */
std::tuple<RsaPrivateKey, RsaPublicKey> RsaDeriveKeysFromPublicExponent(int_fast64_t p, int_fast64_t q, int_fast64_t e, std::string *steps = nullptr);

/* Derive private and public key from p q d. */
std::tuple<RsaPrivateKey, RsaPublicKey> RsaDeriveKeysFromPrivateExponent(int_fast64_t p, int_fast64_t q, int_fast64_t d, std::string *steps = nullptr);
