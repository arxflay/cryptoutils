
**Documentation of usage of the program**

# ec - Elyptic Curve
"GF(p)"
    sum
        curve_a curve_b curve_prime x0 y0 x1 y1
        Example: 
            TODO

    alignson
        curve_a curve_b curve_prime x y
        Example: 
            TODO

"GF(2^n)"
    sum
        curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 x1_g^p/0 y1_g^p/0
        Example:     
            TODO     

    multiply
        curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0 scalar
        Example:
            TODO

    alignson
        curve_a_g^p/0 curve_b_g^p/0 polynomial(binary) ireduciblePolynomial(binary) n x0_g^p/0 y0_g^p/0
        Example:
            TODO
            

# elgamal - El Gamal
enc
    prime generator pubkey otherPartPrivKey(k) message 
    Example: 
        (p, z, b) = (607, 555, 7)
        m = 10
        k = 4
            *elgamal enc 607 555 7 4 10*
dec
    encY encMessage prime privKey
    Example:
        Bob uses the ElGamal cipher with prime number p and a multiplicative integer group generator modup p, z. His prime key is (p, z, a). Alice sent him a ciphertext encrypted with his public key s = (c, d). Decrypt the ciphertext s and find the OT m.
        (p, z, a) = (211, 75, 17) 
        s = (c, d) = (157, 141)
            *elgamal dec 157 141 211 17*

# extgcd - Extended Greatest Common Devisor
    number prime
    Example:
        *extgcd 12 8*

# fermant - Fermant Factorization
    number
    Example:
        Number 442931 is a module in the RSA algorithm. Find the decomposition into prime factors by Fermat by the method.
            **fermant 442931**

# isgenerator - Is Generator
    "GF(p)"
        p numbers...
        Example:
            TODO

    "GF(2^n)"
        ireduciblePolynomial(binary) n polynomials(binary)...
        Example:
            TODO
            
# lhperalt - Lehman Peralt Primality Test
    examinedNumber numbers...
    Example: 
        *lhperalt 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71*

# modexp - Modular Exponentiation
    number power prime
    Example: 
        *modexp 2 3 5*

# rhoalgo - Pollard Rho Algorithm
    number seed
    Example:
        *rhoalgo 442931 2*

# rsa - RSA
    enc
        message n e
        Example:
            *rsa enc 10 23701 1015*

    dec
        encryptedMessage n d
        Example:
            *rsa dec 4754 23701 20903*

    derivePrivKeyFromMod
        n e
        Example:
            *rsa derivePrivKeyFromMod 23701 1015*

    deriveKeysFromPubExp
        p q e
        Example:
            We have two primes p = 137, q = 173. We randomly chose the encryption exponent e = 1015. Determine the public and private key pair.
            p = 137
            q = 173
            e = 1015
                *rsa deriveKeysFromPubExp 137 173 1015*

    deriveKeysFromPrivExp
        p q d
        Example:
            p = 137
            q = 173
            d = 20903
                *rsa deriveKeysFromPrivExp 137 173 20903*
            
# shamir_protocol - Shamir Protocol
    getSubjects
        p K N S x y...
        Example:
        We have a prime p = 13 and a shared secret S = 10. We want to construct a shared secret using Shamir's threshold protocol, divided into N = 3 subjects, with the understanding that to compose we need the participation of K = 2 subjects. We randomly choose a pair in module p: A1=(x1,y1)=(4,8). Determine the other two pairs of numbers for the remaining two subjects.
            p = 13
            K = 2
            N = 3
            S = 10
            A0 = (x0, y0) = (5, 3)
            A1 = (x1, y1) = (4, 8)            
                *shamir_protocol getSubjects 13 2 3 10 4 8*

    reconstruction
        p K x y...
        Example:
            p = 13
            K = 2
            A1 = (x1, y1) = (11, 11)
            A2 = (x2, y2) = (1, 3)
                *shamir_protocol reconstruction 13 2 11 11 1 3*
            