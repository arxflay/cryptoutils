#include "handlers.h"

#include "algos.h"
#include <array>
#include <string_view>
#include <iostream>
#include <map>
#include <functional>
#include <cstring>
#include <cmath>
#include <random>
#include <fmt/core.h>
#include <numeric>


int HandleExtgcd(int argc, const char **argv)
{
    if (argc < 2)
    {
        fmt::print("Enter number prime\n");
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t prime = std::atol(argv[1]);

    auto [nsd, inverseMod] = ExtendedGCD(number, prime);
    
    fmt::print(" nsd result = {}, Result of inverseMod = {} \n", nsd, (inverseMod >= 0 ? inverseMod : prime + inverseMod) );

    return 0;

}

int HandleModExp(int argc, const char **argv)
{
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
    if (argc < 2)
    {
        fmt::print("Enter number seed\n");
        return -1;
    }

    int_fast64_t number = std::atol(argv[0]);
    int_fast64_t seed = std::atol(argv[1]);
    std::string steps;
    auto [a, b] = DoRhoFactorization(number, seed, &steps);
    fmt::print("{}\n{} = {} * {}\n", steps, number, a, b);

    return 0;
}

int HandleLhPeralt(int argc, const char **argv)
{
    if (argc < 2)
    {
        fmt::print("Enter examinedNumber numbers...\n");
        return -1;
    }
    
    int_fast64_t examinedNumber = atol(argv[0]);
    
    std::vector<int_fast64_t> numbers(static_cast<size_t>(argc - 1));
    for (size_t i = 0; i < numbers.size(); i++)
        numbers[i] = atol(argv[i + 1]);
    
    std::string steps;

    LehmanPeraltResult result = LehmanPeralt(numbers, examinedNumber, &steps); 

    fmt::print("{}\n",steps);
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
    static const std::array<std::string_view, 3> commands{ "enc", "dec", "derivePubKey" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("Enter enc/dec\n");
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

        fmt::print("{}:\nEncMsg ({}, {})\n", steps, encMessage.y, encMessage.encData);
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
        fmt::print("{}\ndecrypted message = {}\n", steps, message);
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

        fmt::print("{}\npublic key ({}, {}, {})\n", steps, pubKey.p, pubKey.q, pubKey.y);
    }

    return 0;
}

static inline int HandleEcGF2N(int argc, const char **argv)
{
    static const std::array<std::string_view, 3> commands{ "sum", "multiply", "alignson" };
    std::string_view desired(argc < 2 ? "" : argv[1]);

    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 2 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter GF(2^n) sum/multiply/alignson\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string steps;
    
    auto getGenPoints = [](const std::vector<uint_fast64_t> &elements, std::string_view str) -> uint_fast64_t {
        if (str.size() == 1)
            return 0;

        return (memcmp(str.data(), "g^", std::min(str.size(), static_cast<size_t>(2))) == 0) ? elements.at(atol(str.data() + 2)) : 0; 
    };

    auto getPower = [](const std::vector<uint_fast64_t> &generatorPoints, uint_fast64_t num){
        auto it = std::find(generatorPoints.begin(), generatorPoints.end(), num);
        return std::distance(generatorPoints.begin(), std::find(generatorPoints.begin(), generatorPoints.end(), num));
    };

    if (cmdIndex == 0)
    {
        if (argc < 11)
        {
            fmt::print("Enter GF(2^n) sum curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 x1_g^p/0 y1_g^p/0\n");
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<uint_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
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
            ECPoint newPoint = ECSumGF2N(curve, p, q, generatorPoints, &steps);
            fmt::print("NewPoint = ({}, {})\n"
                    , newPoint.x == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, newPoint.x))
                    , newPoint.y == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, newPoint.y)));
        }
        catch(const std::runtime_error &e)
        {
            fmt::print("{}\n", e.what());
        }
    }
    else if (cmdIndex == 1)
    {
        if (argc < 10)
        {
            fmt::print("Enter GF(2^n) sum curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 scalar\n");
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<uint_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
        ECCurve curve;
        curve.a = getGenPoints(generatorPoints, argv[2]);
        curve.b = getGenPoints(generatorPoints, argv[3]);
        curve.p = (int_fast64_t)std::pow(2, generatorParameters.n) - 1;

        ECPoint p;
        p.x = getGenPoints(generatorPoints, argv[7]);
        p.y = getGenPoints(generatorPoints, argv[8]);

        uint_fast64_t scalar = atol(argv[9]);
        
        try
        {
            ECPoint newPoint = ECMultiplyGF2N(curve, p, scalar, generatorPoints, &steps);
            fmt::print("NewPoint = ({}, {})\n"
                    , newPoint.x == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, newPoint.x))
                    , newPoint.y == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, newPoint.y)));
        }
        catch(const std::runtime_error &e)
        {
            fmt::print("{}\n", e.what());
        }

    }
    else if (cmdIndex == 2)
    {
        if (argc < 9)
        {
            fmt::print("Enter GF(2^n) alignson curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0\n");
            return -1;
        }

        GF2NGeneratorParameters generatorParameters;
        generatorParameters.polynomial = ConvertBinaryToNumber(argv[4]);
        generatorParameters.ireduciblePolynomial = ConvertBinaryToNumber(argv[5]);
        generatorParameters.n = atol(argv[6]);
        std::vector<uint_fast64_t> generatorPoints = GetGF2NGeneratorElements(generatorParameters);
        
        ECCurve curve;
        curve.a = getGenPoints(generatorPoints, argv[2]);
        curve.b = getGenPoints(generatorPoints, argv[3]);
        curve.p = (int_fast64_t)std::pow(2, generatorParameters.n) - 1;

        ECPoint p;
        p.x = getGenPoints(generatorPoints, argv[7]);
        p.y = getGenPoints(generatorPoints, argv[8]);

        fmt::print("Does point ({}, {}) aligns on curve y^2 + xy = x^3 + {}*x^2 + {}?\n"
                , p.x == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, p.x))
                , p.y == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, p.y))
                , curve.a == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, curve.a))
                , curve.b == 0 ? "0" : fmt::format("g^{}", getPower(generatorPoints, curve.b))
        );

        fmt::print("{}\n{}\n", steps, (ECAlignsOnGF2N(curve, p, generatorPoints, generatorParameters, &steps)) ? "Yes" : "No");

    }

    return 0;
}

static inline int HandleEcGFP(int argc, const char **argv)
{
    static const std::array<std::string_view, 2> commands{ "sum", "alignson" };
    std::string_view desired(argc < 2 ? "" : argv[1]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 2 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter GF(p) sum/alignson\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string steps;

    if (cmdIndex == 0)
    {
        if (argc < 9)
        {
            fmt::print("Enter GF(p) sum curve_a curve_b curve_prime x0 y0 x1 y1\n");
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

        fmt::print("{}\nEC sum = ({},{})\n", steps, r.x, r.y);
    }
    else if (cmdIndex == 1)
    {
        if (argc < 7)
        {
            fmt::print("Enter GF(p) alignson curve_a curve_b curve_prime x y\n");
            return -1;
        }

        ECCurve curve;
        curve.a = atol(argv[2]);
        curve.b = atol(argv[3]);
        curve.p = atol(argv[4]); 

        ECPoint p;
        p.x = atol(argv[5]);
        p.y = atol(argv[6]);
        fmt::print("Does point ({}, {}) aligns on 'y^2 mod {} = x^3 + {}x + {} mod {}?\n", p.x, p.y, curve.p, curve.a, curve.b, curve.p);
        fmt::print("{}\n", (ECAlignsOn(curve, p) ? "Yes" : "No"));
    }

    return 0;
}

int HandleEc(int argc, const char **argv)
{   
    static const std::array<std::string_view, 2> commands{ "GF(p)", "GF(2^n)" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter GF(p)/GF(2^n)\n");
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
    static const std::array<std::string_view, 2> commands{ "GF(p)", "GF(2^n)" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter GF(p)/GF(2^n)\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;
    
    if (cmdIndex == 0)
    {
        if (argc < 3)
        {
            fmt::print("{} p numbers...", beginMessage);
            return -1;
        }

        uint_fast64_t p = atol(argv[1]);

        for (int i = 2; i < argc; i++)
        {
            fmt::print ("Is {} a generator?", argv[i]);
            uint_fast64_t number = atol(argv[i]);
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
    static const std::array<std::string_view, 5> commands{ "enc", "dec", "derivePrivKeyFromMod", "deriveKeysFromPubExp", "deriveKeysFromPrivExp" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter enc/dec/derivePrivKeyFromMod/deriveKeysFromPubExp/deriveKeysFromPrivExp\n");
        return -1;
    } 

    ptrdiff_t cmdIndex = std::distance(commands.cbegin(), it);
    std::string beginMessage(fmt::format("enter {} ", desired));
    std::string steps;
    
    if (cmdIndex == 0)
    {
        if (argc < 4)
        {
            fmt::print("{} message n e\n", beginMessage);
            return -1;
        }
        
        RsaPublicKey pubKey;
        pubKey.n = atol(argv[2]);
        pubKey.e = atol(argv[3]);

        uint_fast64_t message = atol(argv[1]);

        uint_fast64_t encrypted = RsaEncrypt(pubKey, message, &steps);
        fmt::print("{}\nEncrypted message = {}\n", steps, encrypted);
    }
    else if (cmdIndex == 1)
    {
        if (argc < 4)
        {
            fmt::print("{} encryptedMessage n d\n", beginMessage);
            return -1;
        }

        RsaPrivateKey privKey;
        privKey.n = atol(argv[2]);
        privKey.d = atol(argv[3]);        

        uint_fast64_t encryptedMessage = atol(argv[1]);

        uint_fast64_t decrypted = RsaDecrypt(privKey, encryptedMessage, &steps);
        fmt::print("{}\nDecrypted message = {}\n", steps, decrypted);
    }
    else if (cmdIndex == 2) 
    {
        if (argc < 3)
        {
            fmt::print("{} n e\n", beginMessage);
            return -1;
        }

        RsaPublicKey pubKey;
        pubKey.n = atol(argv[1]);
        pubKey.e = atol(argv[2]);

        RsaPrivateKey privKey = RsaDerivePrivateKeyFromModule(pubKey, &steps);
        fmt::print("{}\nPrivate key: (n, d) = ({}, {})\n", steps, privKey.n, privKey.d);
    }
    else if (cmdIndex == 3) 
    {
        if (argc < 4)
        {
            fmt::print("{} p q e\n", beginMessage);
            return -1;
        }

        uint_fast64_t p = atol(argv[1]);
        uint_fast64_t q = atol(argv[2]);
        uint_fast64_t e = atol(argv[3]);

        auto keys = RsaDeriveKeysFromPublicExponent(p, q, e, &steps); 
        RsaPrivateKey privKey = std::get<0>(keys);
        RsaPublicKey pubKey = std::get<1>(keys);

        fmt::print("{}\nPrivate key: (n, d) = ({}, {})\nPublic key: (n, e) = ({}, {})\n", steps, privKey.n, privKey.d, pubKey.n, pubKey.e);
    }
    else if (cmdIndex == 4) 
    {
        if (argc < 4)
        {
            fmt::print("{} p q d\n", beginMessage);
            return -1;
        }

        uint_fast64_t p = atol(argv[1]);
        uint_fast64_t q = atol(argv[2]);
        uint_fast64_t d = atol(argv[3]);

        auto keys = RsaDeriveKeysFromPrivateExponent(p, q, d, &steps); 
        RsaPrivateKey privKey = std::get<0>(keys);
        RsaPublicKey pubKey = std::get<1>(keys);

        fmt::print("{}\nPrivate key: (n, d) = ({}, {})\nPublic key: (n, e) = ({}, {})\n", steps, privKey.n, privKey.d, pubKey.n, pubKey.e);
    }

    return 0;
}

int HandleShamirProtocol(int argc, const char **argv)
{
    static const std::array<std::string_view, 2> commands{ "getSubjects", "reconstruction" };
    std::string_view desired(argc < 1 ? "" : argv[0]);
    
    auto containsCommand = [desired](std::string_view command){ return desired.compare(command) == 0; };
    decltype(commands)::const_iterator it;

    if (argc < 1 || (it = std::find_if(commands.cbegin(), commands.cend(), containsCommand)) == commands.cend())
    {
        fmt::print("enter getSubjects/reconstruction\n");
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

        std::vector<ShamirSubject> subjects = GetShamirSubjects(parameters);

        fmt::print("List of all points: ");
        for (const ShamirSubject &subj : subjects)
            fmt::print("({}, {}) ", subj.x, subj.y);
        fmt::print("\n");
        fmt::print("Picking {} random points from list...\n", n);

        std::vector<int_fast64_t> subjectsIndicies(subjects.size());
        std::iota(subjectsIndicies.begin(), subjectsIndicies.end(), 0);
        std::mt19937 generator(std::random_device{}());
        std::shuffle(subjectsIndicies.begin(), subjectsIndicies.end(), generator);

        fmt::print("Result: ");
        for (int_fast64_t i = 0; i < n; i++)
        {
            ShamirSubject &subject = subjects.at(subjectsIndicies[i]);
            fmt::print("({}, {}) ", subject.x, subject.y);
        }

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
        fmt::print("{}", steps);
        fmt::print("S = ({}, {})\n", subject.x, subject.y);
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
