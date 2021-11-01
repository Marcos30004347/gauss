#ifndef GALOIS_FIELD_H
#define GALOIS_FIELD_H

#include "Core/Polynomial/Polynomial.hpp"

namespace galoisField {

/**
 * @brief Compute the representation of u(x) in Z[x] over Zp[x]
 * 
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The representation of u(x) over Zp[x]
 */
ast::AST* gf(ast::AST* u, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Compute the representation of u in Z[x] over Zp[x]
 * 
 * @param u A polynomial in Z[x]
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The representation of u(x) over Zp[x]
 */
ast::AST* gf(ast::AST* u, Int p, bool symmetric);

/**
 * @brief Divide a polynomial a(x) by b(x) over Zp[x].
 * 
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return A list with que quotient and remainder over Zp[x]
 */
ast::AST* divPolyGf(ast::AST* a, ast::AST* b, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Compute the remainder of the division of a(x) by b(x) over Zp[x].
 * 
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The remainder of a(x) by b(x) over Zp[x]
 */
ast::AST* remPolyGf(ast::AST* a, ast::AST* b, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Compute the quotient of the division of a(x) by b(x) over Zp[x].
 * 
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The quotient of a(x) by b(x) over Zp[x]
 */
ast::AST* quoPolyGf(ast::AST* a, ast::AST* b, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Compute the monic form of f(x) over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return A list with the content of f(x) and the monic form of f(x) over Zp[x]
 */
ast::AST* monicPolyGf(ast::AST* f, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Computes the gcd of a(x) and b(x) over Zp[x]
 * 
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The gcd between a(x) and b(x) over Zp[x]
 */
ast::AST* gcdPolyGf(ast::AST* a, ast::AST* b, ast::AST* x, Int p, bool sym = true);


/**
 * @brief Add two polynomials over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The result of f(x) + g(x) over Zp[x]
 */
ast::AST* addPolyGf(ast::AST* f, ast::AST* g, ast::AST* x, Int p, bool sym = true);


/**
 * @brief Subtract two polynomials over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The result of f(x) - g(x) over Zp[x]
 */
ast::AST* subPolyGf(ast::AST* f, ast::AST* g, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Multiply two polynomials over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The result of f(x) * g(x) over Zp[x]
 */
ast::AST* mulPolyGf(ast::AST* f, ast::AST* g, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Compute f(x)^n in Zp[x]/g(x), that is  f(x)^n mod g(x) over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param n An integer
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @returnf(x)^n mod g(x) over Zp[x]
 */
ast::AST* powModPolyGf(ast::AST* f, ast::AST* g, ast::AST* x, Int n, Int p, bool sym = true);


/**
 * @brief Return a random polynomial over Zp[x]
 * 
 * @param d The degree of the polynomial
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return A random polynomial with degree 'd' over Zp[x]
 */
ast::AST* randPolyGf(Int d, ast::AST* x, Int p, bool sym = true);

/**
 * @brief Extended Euclidean Algorithm over Zp[x]
 * 
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return A list with [gcd, s, t] such as s*f + g*t = gcd over Zp[x]
 */
ast::AST* extendedEuclidGf(ast::AST* f, ast::AST* g, ast::AST* x, Int p, bool sym = true);


/**
 * @brief Computes the quotient of the division of s by t in Zp[x]
 * 
 * @param s An integer
 * @param t An integer
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The quotient of s / t in Zp[x]
 */
Int quoGf(Int s, Int t, Int p,  bool sym = true);

/**
 * @brief Computes the inverse of a in Zp
 * 
 * @param a A integer
 * @param p A prime number
 * @param sym true if Zp[x] is in symmetric representation 
 * 						or false otherwise, true by default
 * @return The inverse of a over Zp 
 */
Int inverseGf(Int a, Int p, bool symmetric = true); 

/**
 * @brief Computes the remainder of a divided by b
 * 
 * @param a A integer
 * @param b A integer
 * @param sym true if result should be in symmetric representation 
 * 						or false otherwise
 * @return The remainder of a / b
 */
Int mod(Int a, Int b, bool sym);

}

#endif
