#include "handlers.h"

#include "algos.h"
#include "generators.h"
#include "ec.h"
#include "rsa.h"
#include "shamir.h"
#include "elgamal.h"

#include <array>
#include <string_view>
#include <iostream>
#include <map>
#include <functional>
#include <cstring>
#include <cmath>
#include <fmt/core.h>

int HandleExtgcd(int argc, const char **argv)
{
    fmt::print("[Extended Euclidian Algorithm]\n");
    if (argc < 2)
    {
        fmt::print("Enter number prime\n");
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t prime = std::atol(argv[1]);

    auto [nsd, inverseMod] = ExtendedGCD(number, prime);
    fmt::print("<result>\n");
    fmt::print("GCD = {}\ninverseMod of smallest = {}\n", nsd, (inverseMod >= 0 ? inverseMod : prime + inverseMod) );

    return 0;

}

int HandleModExp(int argc, const char **argv)
{
    fmt::print("[Fast Modular Exponention]\n");
    if (argc < 3)
    {
        fmt::print("Enter number power prime\n");
        return -1;
    }

    int_fast64_t number = atol(argv[0]);
    int_fast64_t power = atol(argv[1]);
    int_fast64_t prime = atol(argv[2]);

    fmt::print("{}^{} mod {} = {}",number, power, prime, ModExp(number, power, prime));

    return 0;
}


int HandleFermantFactorization(int argc, const char **argv)
{
    fmt::print("[Fermant Factorization]\n");
    if (argc < 1)
    {
        fmt::print("Enter number\n");
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    auto [a, b] = DoFermantFactorization(number);
    fmt::print("{} = {} * {}\n", number, a, b);
    return 0;
}

int HandleRhoFactorization(int argc, const char **argv)
{
    fmt::print("[Pollard Rho Factorization]\n");
    if (argc < 2)
    {
        fmt::print("Enter number seed\n");
        return -1;
    }
    
    fmt::print("<input>\n");
    fmt::print("-> number = {}\n", argv[0]);
    fmt::print("-> seed = {}\n", argv[1]);
    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t seed = std::atol(argv[1]);
    std::string steps;
    auto [a, b] = DoRhoFactorization(number, seed, &steps);
    fmt::print("<steps>\n");
    fmt::print("{}\n", steps);
    fmt::print("<result>\n");
    fmt::print("{} * {}\n", number, a, b);

    return 0;
}

int HandleLhPeralt(int argc, const char **argv)
{
    fmt::print("[Lehman Peralta Test]\n");
    if (argc < 2)
    {
        fmt::print("Enter examinedNumber numbers...\n");
        return -1;
    }
    
    fmt::print("<input>\n");
    int_fast64_t examinedNumber = atol(argv[0]);
    fmt::print("-> examinedNumber = {}\n", examinedNumber);
    
    fmt::print("-> numbers = [");
    std::vector<int_fast64_t> numbers(static_cast<size_t>(argc - 1));
    for (size_t i = 0; i < numbers.size(); i++)
    {
        numbers[i] = atol(argv[i + 1]);
        if (i + 1 == numbers.size())
            fmt::print("{}]\n", numbers[i]);
        else
            fmt::print("{}, ", numbers[i]);
    }

    std::string steps;
    LehmanPeraltResult result = LehmanPeralt(numbers, examinedNumber, &steps); 

    fmt::print("<steps>\n");
    fmt::print("{}\n", steps);
    fmt::print("<result>\n");
    if (result.flags & LehmanPeraltFlags::NOT_PRIME)
        fmt::print("Number is composite, not prime\n");
    else
    {
        std::string_view resultText(( (result.flags & LehmanPeraltFlags::MINUS_ONE_OCCURED) ? "prime" : "non-prime" ));
        fmt::print("Number is {} with chance {}%\n", resultText, result.chance);
    }

    return 0;
}

int HandleElGamal(int argc, const char **argv)
{
    fmt::print("[El Gamal]\n");
    static const std::array<std::string_view, 3> commands{ "enc", "dec", "derivePubKey" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter ");
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }

        fmt::print("\n");
        return -1;
    }

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;

    if (cmdIndex == 0)
    {
        if(argc < 6)
        {
            fmt::print("{}prime generator pubkey otherPartPrivKey(k) message\n", beginMessage);
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
        
        fmt::print("<steps>\n");
        fmt::print("Warning: skip Public key step if it isn't explicitly required in exercise\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Warning: skip Decipher Public key and Combined if it isn't explicitly required in exercise\n");
        fmt::print("Decipher Public key = {}\nEncrypted message = {}\nCombined ({}, {})\n", encMessage.y, encMessage.encData, encMessage.y, encMessage.encData);
    }
    else if (cmdIndex == 1)
    {
        if (argc < 5)
        {
            fmt::print("{}encY encMessage prime privKey \n", beginMessage);
            return -1;
        }
        
        ElGamalData data;
        data.y = atol(argv[1]);
        data.encData = atol(argv[2]);
        
        ElGamalPrivateKey privKey;
        privKey.p = atol(argv[3]);
        privKey.k = atol(argv[4]);

        int_fast64_t message = ElGamalDecrypt(data, privKey, &steps);
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Decrypted message = {}\n", message);
    }
    else
    {
        if (argc < 4)
        {
            fmt::print("{}prime generator privKey\n", beginMessage);
            return -1;
        }

        ElGamalPrivateKey privKey;
        privKey.p = atol(argv[1]);
        privKey.q = atol(argv[2]);
        privKey.k = atol(argv[3]);

        ElGamalPublicKey pubKey = ElGamalDerivePublicKey(privKey, &steps);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Public key ({}, {}, {})\n", pubKey.p, pubKey.q, pubKey.y);
    }

    return 0;
}

static inline int HandleEcGF2N(int argc, const char **argv)
{
    static const std::array<std::string_view, 3> commands{ "sum", "multiply", "alignson" };
    std::string_view desired(argc < 2 ? "" : argv[1]);
    
    std::string beginMessage = "Enter GF(2^n) ";

    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 2 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print(beginMessage);
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");

        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string steps;
    
    auto getGenPoints = [](const std::vector<int_fast64_t> &elements, std::string_view str) -> int_fast64_t {
        if (str.size() == 1)
            return 0;

        return (memcmp(str.data(), "g^", std::min(str.size(), static_cast<size_t>(2))) == 0) ? elements.at(atol(str.data() + 2)) : 0; 
    };

    beginMessage = fmt::format("{}{} ", beginMessage, desired);
    if (cmdIndex == 0)
    {
        if (argc < 11)
        {
            fmt::print("{}curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 x1_g^p/0 y1_g^p/0\n", beginMessage);
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<int_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
        ECCurve curve;
        curve.a = getGenPoints(generatorPoints, argv[2]);
        curve.b = getGenPoints(generatorPoints, argv[3]);
        curve.p = (int_fast64_t)std::pow(2, generatorParameters.n) - 1;

        ECPoint p;
        p.x = getGenPoints(generatorPoints, argv[7]);
        p.y = getGenPoints(generatorPoints, argv[8]);

        ECPoint q;
        q.x = getGenPoints(generatorPoints, argv[9]);
        q.y = getGenPoints(generatorPoints, argv[10]);
        
        try
        {
            ECPoint newPoint = ECSumGF2N(curve, p, q, generatorPoints, generatorParameters, &steps);
            std::string xStr;
            std::string yStr;
            ECPointToStrGF2N(newPoint, generatorPoints, xStr, yStr);
            
            fmt::print("<steps>\n");
            fmt::print("{}\n", steps);
            fmt::print("<result>\n");
            fmt::print("R = ({}, {})\n", xStr, yStr);
        }
        catch(const std::runtime_error &e)
        {
            fmt::print("<steps>\n");
            fmt::print("{}\n", steps);
            fmt::print("<result>\n");
            fmt::print("{}\n", e.what());
        }
    }
    else if (cmdIndex == 1)
    {
        if (argc < 10)
        {
            fmt::print("{}curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 scalar\n", beginMessage);
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<int_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
        ECCurve curve;
        curve.a = getGenPoints(generatorPoints, argv[2]);
        curve.b = getGenPoints(generatorPoints, argv[3]);
        curve.p = (int_fast64_t)std::pow(2, generatorParameters.n) - 1;

        ECPoint p;
        p.x = getGenPoints(generatorPoints, argv[7]);
        p.y = getGenPoints(generatorPoints, argv[8]);

        int_fast64_t scalar = atol(argv[9]);
        
        try
        {
            ECPoint newPoint = ECMultiplyGF2N(curve, p, scalar, generatorPoints, generatorParameters, &steps);
            
            std::string xStr;
            std::string yStr;
            ECPointToStrGF2N(newPoint, generatorPoints, xStr, yStr);
            
            fmt::print("<steps>\n");
            fmt::print("{}\n", steps);
            fmt::print("<result>\n");
            fmt::print("R = ({}, {})\n", xStr, yStr);
        }
        catch(const std::runtime_error &e)
        {
            fmt::print("<result>\n");
            fmt::print("{}\n", e.what());
        }
    }
    else if (cmdIndex == 2)
    {
        if (argc < 9)
        {
            fmt::print("{}curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0\n", beginMessage);
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<int_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
        ECCurve curve;
        curve.a = getGenPoints(generatorPoints, argv[2]);
        curve.b = getGenPoints(generatorPoints, argv[3]);
        curve.p = (int_fast64_t)std::pow(2, generatorParameters.n) - 1;

        ECPoint p;
        p.x = getGenPoints(generatorPoints, argv[7]);
        p.y = getGenPoints(generatorPoints, argv[8]);
        
        std::string xStr;
        std::string yStr;
        ECPointToStrGF2N(p, generatorPoints, xStr, yStr);

        std::string aStr;
        std::string bStr;
        ECCurveToStrGF2N(curve, generatorPoints, aStr, bStr);

        fmt::print("Does point ({}, {}) aligns on curve y^2 + xy = x^3 + {}*x^2 + {}?\n", xStr, yStr, aStr, bStr);
        bool alignsOn = ECAlignsOnGF2N(curve, p, generatorPoints, generatorParameters, &steps);

        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("{}\n", alignsOn ? "Yes, it does" : "No, it doesn't");
    }

    return 0;
}

static inline int HandleEcGFP(int argc, const char **argv)
{
    static const std::array<std::string_view, 2> commands{ "sum", "alignson" };
    std::string_view desired(argc < 2 ? "" : argv[1]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;
    
    std::string beginMessage = "Enter GF(p) ";
    if (argc < 2 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter GF(p) ");
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");

        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string steps;

    beginMessage = fmt::format("{}{} ", beginMessage, desired);
    if (cmdIndex == 0)
    {
        if (argc < 9)
        {
            fmt::print("{}curve_a curve_b curve_prime x0 y0 x1 y1\n", beginMessage);
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

        try
        {
            ECPoint r = ECSum(curve, p, q, &steps);
            fmt::print("<steps>\n");
            fmt::print("{}\n", steps);
            fmt::print("<result>\n");
            fmt::print("R = ({},{})\n", r.x, r.y);
        }
        catch (const std::runtime_error &e)
        {
            fmt::print("<steps>\n");
            fmt::print("{}\n", steps);
            fmt::print("<result>\n");
            fmt::print("{}\n", e.what());
        }
    }
    else if (cmdIndex == 1)
    {
        if (argc < 7)
        {
            fmt::print("{}curve_a curve_b curve_prime x y\n", beginMessage);
            return -1;
        }

        ECCurve curve;
        curve.a = atol(argv[2]);
        curve.b = atol(argv[3]);
        curve.p = atol(argv[4]); 

        ECPoint p;
        p.x = atol(argv[5]);
        p.y = atol(argv[6]);

        fmt::print("Does point ({}, {}) aligns on y^2 mod {} = x^3 + {}x + {} mod {}?\n", p.x, p.y, curve.p, curve.a, curve.b, curve.p);
        bool aligns = ECAlignsOn(curve, p, &steps);   
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("{}\n", aligns ? "Yes, it does" : "No, it doesn't");
    }

    return 0;
}

int HandleEc(int argc, const char **argv)
{
    fmt::print("[EC]\n");
    static const std::array<std::string_view, 2> commands{ "GF(p)", "GF(2^n)" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter ");
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");
        
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);

    if (cmdIndex == 0)
        return HandleEcGFP(argc, argv);
    else if(cmdIndex == 1)
        return HandleEcGF2N(argc, argv);

    return -1;
}

int HandleIsGenerator(int argc, const char **argv)
{
    fmt::print("[Generator Tools]\n");
    static const std::array<std::string_view, 2> commands{ "GF(p)", "GF(2^n)" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter ");
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;
    
    if (cmdIndex == 0)
    {
        if (argc < 3)
        {
            fmt::print("{}p numbers...", beginMessage);
            return -1;
        }

        int_fast64_t p = atol(argv[1]);

        for (int i = 2; i < argc; i++)
        {
            fmt::print ("Is {} a generator?", argv[i]);
            int_fast64_t number = atol(argv[i]);
            bool result = IsGeneratorGFP(number, p, &steps);
            fmt::print("{}\nnumber {} {}\n", steps, number, (result ? "is generator" : "is not generator"));
        }
    }
    else if (cmdIndex == 1)
    {
        if (argc < 4)
        {
            fmt::print("{} ireduciblePolynomial(binary) n polynomials(binary)...\n", beginMessage);
            return -1;
        }
        GF2NGeneratorParameters parameters;
        parameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[1]);
        parameters.n = atol(argv[2]);
        for (int i = 3; i < argc; i++)
        {
            fmt::print("Is polynomial {} a generator?\n", argv[i]);
            parameters.polynomial = ConvertBinaryToNumber(argv[i]);
            bool result = IsGeneratorGF2N(parameters, &steps);
            fmt::print("{}\npolynomial {} {}\n", steps, argv[i], (result ? "is generator" : "is not generator"));
        }
    }

    return 0;
}

int HandleRsa(int argc, const char **argv)
{
    fmt::print("[RSA Tools]\n");
    static const std::array<std::string_view, 5> commands{ "enc", "dec", "derivePrivKeyFromMod", "deriveKeysFromPubExp", "deriveKeysFromPrivExp" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    std::string beginMessage = "Enter ";
    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print(beginMessage);
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");
        return -1;
    } 
    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string steps;
    
    beginMessage = fmt::format("{}{} ", beginMessage, desired);
    if (cmdIndex == 0)
    {
        if (argc < 4)
        {
            fmt::print("{}message n e\n", beginMessage);
            return -1;
        }
        
        RsaPublicKey pubKey;
        pubKey.n = atol(argv[2]);
        pubKey.e = atol(argv[3]);

        int_fast64_t message = atol(argv[1]);
        int_fast64_t encrypted = RsaEncrypt(pubKey, message, &steps);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Encrypted message = {}\n", encrypted);
    }
    else if (cmdIndex == 1)
    {
        if (argc < 4)
        {
            fmt::print("{}encryptedMessage n d\n", beginMessage);
            return -1;
        }

        RsaPrivateKey privKey;
        privKey.n = atol(argv[2]);
        privKey.d = atol(argv[3]);        

        int_fast64_t encryptedMessage = atol(argv[1]);

        int_fast64_t decrypted = RsaDecrypt(privKey, encryptedMessage, &steps);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Decrypted message = {}\n", decrypted);
    }
    else if (cmdIndex == 2) 
    {
        if (argc < 3)
        {
            fmt::print("{}n e\n", beginMessage);
            return -1;
        }

        RsaPublicKey pubKey;
        pubKey.n = atol(argv[1]);
        pubKey.e = atol(argv[2]);

        RsaPrivateKey privKey = RsaDerivePrivateKeyFromModule(pubKey, &steps);
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Private key: (n, d) = ({}, {})\n", privKey.n, privKey.d);

    }
    else if (cmdIndex == 3) 
    {
        if (argc < 4)
        {
            fmt::print("{}p q e\n", beginMessage);
            return -1;
        }

        int_fast64_t p = atol(argv[1]);
        int_fast64_t q = atol(argv[2]);
        int_fast64_t e = atol(argv[3]);

        auto keys = RsaDeriveKeysFromPublicExponent(p, q, e, &steps); 
        RsaPrivateKey privKey = std::get<0>(keys);
        RsaPublicKey pubKey = std::get<1>(keys);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Private key: (n, d) = ({}, {})\nPublic key: (n, e) = ({}, {})\n", privKey.n, privKey.d, pubKey.n, pubKey.e);
    }
    else if (cmdIndex == 4) 
    {
        if (argc < 4)
        {
            fmt::print("{}p q d\n", beginMessage);
            return -1;
        }

        int_fast64_t p = atol(argv[1]);
        int_fast64_t q = atol(argv[2]);
        int_fast64_t d = atol(argv[3]);

        auto keys = RsaDeriveKeysFromPrivateExponent(p, q, d, &steps); 
        RsaPrivateKey privKey = std::get<0>(keys);
        RsaPublicKey pubKey = std::get<1>(keys);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("Private key: (n, d) = ({}, {})\nPublic key: (n, e) = ({}, {})\n", privKey.n, privKey.d, pubKey.n, pubKey.e);
    }

    return 0;
}

int HandleShamirProtocol(int argc, const char **argv)
{
    fmt::print("[Shamir protocol tools]\n");
    static const std::array<std::string_view, 2> commands{ "getSubjects", "reconstruction" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter ");
        for (size_t i = 0; i < commands.size(); i++)
        {
            if (i + 1 != commands.size())
                fmt::print("{}/", commands.at(i));
            else
                fmt::print("{}", commands.at(i));
        }
        fmt::print("\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;

    if (cmdIndex == 0)
    {
        if (argc < 5)
        {
            fmt::print("{}p K N S x y...\n", beginMessage);
            return -1;
        }

        ShamirParameters parameters;
        parameters.p = atol(argv[1]);
        int_fast64_t k = atol(argv[2]);

        if (k <= 1)
        {
            fmt::print("K has invalid value, {} <= 1\n", k);
            return -1;
        }

        int_fast64_t expectedArgc = (k - 1) * 2 + 5;
        
        if (argc != expectedArgc)
        {
            if (argc < expectedArgc)
                fmt::print("Invalid amount of x y points or point is incomplete\n");
            else
                fmt::print("too many args\n");

            fmt::print("expectedArgc {} but was {}\n", expectedArgc, argc);
            return -1;
        }
        
        int_fast64_t n = atol(argv[3]);
        if (n > parameters.p)
        {
            fmt::print("Amount of needed subjects N is highier than p, {} > {}\n", n, parameters.p);
            return -1;
        }

        int_fast64_t s = atol(argv[4]);
        
        parameters.subjects.push_back(ShamirSubject{ 0, s });

        for (int_fast64_t i = 5; i < argc; i += 2)
            parameters.subjects.push_back(ShamirSubject{ atol(argv[i]), atol(argv[i + 1]) });
        
        std::string steps;
        std::vector<ShamirSubject> subjects = GetShamirSubjects(parameters, n, &steps);
        
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("N Subjects: ");
        for (const ShamirSubject &subject : subjects)
            fmt::print("({}, {}) ", subject.x, subject.y);

        fmt::print("\n");
    }
    else if (cmdIndex == 1)
    {
        if (argc < 3)
        {
            fmt::print("{}p K x y...\n", beginMessage);
            return -1;
        }

        ShamirParameters parameters;
        parameters.p = atol(argv[1]);
        int_fast64_t k = atol(argv[2]);

        if (k <= 1)
        {
            fmt::print("K has invalid value, {} <= 1\n", k);
            return -1;
        }

        int_fast64_t expectedArgc = (k * 2) + 3;
        
        if (argc != expectedArgc)
        {
            if (argc < expectedArgc)
                fmt::print("Invalid amount of x y points or point is incomplete\n");
            else
                fmt::print("too many args\n");

            fmt::print("expectedArgc {} but was {}\n", expectedArgc, argc);
            return -1;
        }

        for (int_fast64_t i = 3; i < argc; i += 2)
            parameters.subjects.push_back(ShamirSubject{ atol(argv[i]), atol(argv[i + 1]) });
        
        ShamirSubject subject = DoShamirReconstruction(parameters, &steps);
        fmt::print("<steps>\n");
        fmt::print("{}\n", steps);
        fmt::print("<result>\n");
        fmt::print("S = {}\n", subject.y);
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
        { "ec"     , &HandleEc },
        { "isgenerator", &HandleIsGenerator },
        { "rsa", &HandleRsa },
        { "shamir_protocol", &HandleShamirProtocol }
    };

    return handlers;
}
