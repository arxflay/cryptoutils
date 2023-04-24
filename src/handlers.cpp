#include "handlers.h"

#include "algos.h"
#include <array>
#include <string_view>
#include <iostream>
#include <map>
#include <functional>
#include <cstring>
#include <fmt/core.h>

int HandleExtgcd(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "enter number prime" << '\n';
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t prime = std::atol(argv[1]);

    auto [nsd, inverseMod] = ExtendedGCD(number, prime);
    
    std::cout << "nsd = " << nsd << "; inverseMod = " << (inverseMod >= 0 ? inverseMod : prime + inverseMod) << '\n';

    return 0;

}

int HandleModExp(int argc, const char **argv)
{
    if (argc < 3)
    {
        std::cout << "enter number power prime" << '\n';
        return -1;
    }

    int_fast64_t number = atol(argv[0]);
    int_fast64_t power = atol(argv[1]);
    int_fast64_t prime = atol(argv[2]);

    std::cout << number << "^" << power << " mod " << prime << " = " <<  ModExp(number, power, prime);

    return 0;
}


int HandleFermantFactorization(int argc, const char **argv)
{
    if (argc < 1)
    {
        std::cout << "enter number" << '\n';
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    auto [a, b] = DoFermantFactorization(number);
    std::cout << number << " = " << a << " * " << b << '\n';

    return 0;
}

int HandleRhoFactorization(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "enter number seed" << '\n';
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t seed = std::atol(argv[1]);
    std::string steps;
    auto [a, b] = DoRhoFactorization(number, seed, &steps);
    std::cout << steps << '\n' << number << " = " << a << " * " << b << '\n';

    return 0;
}

int HandleLhPeralt(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "enter examinedNumber numbers..." << '\n';
        return -1;
    }
    
    int_fast64_t examinedNumber = atol(argv[0]);
    
    std::vector<int_fast64_t> numbers(static_cast<size_t>(argc - 1));
    for (size_t i = 0; i < numbers.size(); i++)
        numbers[i] = atol(argv[i + 1]);
    
    std::string steps;

    LehmanPeraltResult result = LehmanPeralt(numbers, examinedNumber, &steps); 

    std::cout << steps << '\n';
    if (result.flags & LehmanPeraltFlags::NOT_PRIME)
        std::cout << "number is composite, not prime" << '\n';
    else
        std::cout << "number is " 
                  << ((result.flags & LehmanPeraltFlags::MINUS_ONE_OCCURED) ? "prime" : "non-prime")
                  << " with a chance " << result.chance << "%" << '\n';

    return 0;
}

int HandleElGamal(int argc, const char **argv)
{
    static const std::array<std::string_view, 3> commands{ "enc", "dec", "derivePubKey" };
    std::string_view desired(argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;


    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        std::cout << "enter enc/dec" << '\n';
        return -1;
    }

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;

    if (cmdIndex == 0)
    {
        if(argc < 6)
        {
            std::cout << beginMessage << "prime generator pubkey otherPartPrivKey(k) message" << '\n';
            return -1;
        }
        
        ElGamalPublicKey pubKey;
        pubKey.p = atol(argv[1]);
        pubKey.q = atol(argv[2]);
        pubKey.y = atol(argv[3]);
        
        ElGamalPrivateKey privKey;
        privKey.k = atol(argv[4]);
        
        int_fast64_t message = atol(argv[5]);

        ElGamalData encMessage = ElGamalEncrypt(pubKey, privKey, message, &steps);

        std::cout << steps << '\n' << fmt::format("EncMsg ({}, {})", encMessage.y, encMessage.encData) << '\n';  ;
    }
    else if (cmdIndex == 1)
    {
        if (argc < 5)
        {
            std::cout << beginMessage << "encY encMessage prime privKey" << '\n'; 
            return -1;
        }
        
        ElGamalData data;
        data.y = atol(argv[1]);
        data.encData = atol(argv[2]);
        
        ElGamalPrivateKey privKey;
        privKey.p = atol(argv[3]);
        privKey.k = atol(argv[4]);

        int_fast64_t message = ElGamalDecrypt(data, privKey, &steps);

        std::cout << steps << '\n' << "decrypted message = " << message << '\n';
    }
    else
    {
        if (argc < 4)
        {
            std::cout << beginMessage << "prime generator privKey" << '\n';
            return -1;
        }

        ElGamalPrivateKey privKey;
        privKey.p = atol(argv[1]);
        privKey.q = atol(argv[2]);
        privKey.k = atol(argv[3]);

        ElGamalPublicKey pubKey = ElGamalDerivePublicKey(privKey, &steps);

        std::cout << steps << '\n' << fmt::format("public key ({}, {}, {})", pubKey.p, pubKey.q, pubKey.y) << '\n';
    }

    return 0;
}

int HandleEcc(int argc, const char **argv)
{
    if (argc < 1)
    {
        std::cout << "enter GF(p)/GF(2^n)" << '\n';
        return -1;
    }

    std::string steps;

    if(std::strcmp(argv[0], "GF(p)") == 0)
    {
        if (argc < 2)
        {
            std::cout << "enter GF(p) sum/aligns" << '\n';
            return -1;
        }
        
        if (std::strcmp(argv[1], "sum") == 0)
        {
            if (argc < 9)
            {
                std::cout << "enter GF(p) sum curve_a curve_b curve_prime x0 y0 x1 y1" << '\n';
                return -1;
            }

            ECCurve curve;
            curve.a = atol(argv[2]);
            curve.b = atol(argv[3]);
            curve.p = atol(argv[4]);
            
            ECPoint p;
            p.x = atol(argv[5]);
            p.y = atol(argv[6]);

            ECPoint q;
            q.x = atol(argv[7]);
            q.y = atol(argv[8]);

            ECPoint r = ECSum(curve, p, q, &steps);

            std::cout << steps << '\n' << "EC sum = (" << r.x << ", " << r.y << ")\n";
        }
        else if (std::strcmp(argv[1], "aligns") == 0)
        {
            if (argc < 7)
            {
                std::cout << "enter GF(p) aligns curve_a curve_b curve_prime x y" << '\n';
                return -1;
            }

            ECCurve curve;
            curve.a = atol(argv[2]);
            curve.b = atol(argv[3]);
            curve.p = atol(argv[4]); 

            ECPoint p;
            p.x = atol(argv[5]);
            p.y = atol(argv[6]);

            std::cout << (ECAlignsOn(curve, p) ? "true" : "false") << '\n';
        }

    }

    return 0;
}

const std::map<std::string_view, UtilHandler> &GetUtilHandlers()
{
    static std::map<std::string_view, UtilHandler> handlers = {
        { "extgcd"  , &HandleExtgcd  },
        { "modexp"  , &HandleModExp  },
        { "fermant" , &HandleFermantFactorization },
        { "rhoalgo" , &HandleRhoFactorization },
        { "lhperalt", &HandleLhPeralt},
        { "elgamal" , &HandleElGamal },
        { "ecc"     , &HandleEcc }
    };

    return handlers;
}
