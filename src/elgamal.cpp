#include "elgamal.h"

#include "algos.h"
#include <fmt/core.h>

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
