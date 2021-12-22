#ifndef GALOIS_FIELD_H
#define GALOIS_FIELD_H

#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Algebra/List.hpp"

// TODO: optimize those functions using reference arguments and move semantics

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
//ast::Expr gf(ast::Expr u, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Compute the representation of u in Z[x] over Zp[x]
 *
 * @param u A polynomial in Z[x]
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The representation of u(x) over Zp[x]
 */
ast::Expr gf(ast::Expr u, Int p, bool symmetric);


/**
 * @brief Compute the representation of a polynomial expression in Z[x] over Zp[x]
 *
 * @param u A polynomial expression in Z[x]
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The representation of u(x) over Zp[x]
 */
ast::Expr gfPolyExpr(ast::Expr u, Int p, bool symmetric);


/**
 * @brief Compute the representation of u(x) in Z[x] over Zp[x] without expansion
 *
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The representation of u(x) over Zp[x]
 */
ast::Expr groundGf(ast::Expr u, Int s, bool symmetric);

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
ast::Expr divPolyGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Divide a polynomial expression a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return A list with que quotient and remainder over Zp[x]
 */
ast::Expr divPolyExprGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);



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
ast::Expr remPolyGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Compute the remainder of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The remainder of a(x) by b(x) over Zp[x]
 */
ast::Expr remPolyExprGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

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
ast::Expr quoPolyGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Compute the quotient of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The quotient of a(x) by b(x) over Zp[x]
 */
ast::Expr quoPolyExprGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

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
ast::Expr monicPolyGf(ast::Expr f, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Compute the monic form of f(x) over Zp[x]
 *
 * @param f A polynomial Expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return A list with the content of f(x) and the monic form of f(x) over Zp[x]
 */
ast::Expr monicPolyExprGf(ast::Expr f, ast::Expr x, Int p, bool sym = true);

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
ast::Expr gcdPolyGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Computes the gcd of a(x) and b(x) over Zp[x]
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return The gcd between a(x) and b(x) over Zp[x]
 */
ast::Expr gcdPolyExprGf(ast::Expr a, ast::Expr b, ast::Expr x, Int p, bool sym = true);


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
ast::Expr addPolyGf(ast::Expr f, ast::Expr g, ast::Expr x, Int p, bool sym = true);


/**
 * @brief Add two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by default
 * @return The result of f(...) + g(...) over Zp
 */
ast::Expr addPolyExprGf(ast::Expr f, ast::Expr g, Int p, bool sym = true);

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
ast::Expr subPolyGf(ast::Expr f, ast::Expr g, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Sub two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by default
 * @return The result of f(...) - g(...) over Zp
 */
ast::Expr subPolyExprGf(ast::Expr f, ast::Expr g, Int p, bool sym = true);

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
ast::Expr mulPolyGf(ast::Expr f, ast::Expr g, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Multiply two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by default
 * @return The result of f(...) * g(...) over Zp
 */
ast::Expr mulPolyExprGf(ast::Expr f, ast::Expr g, Int p, bool sym = true);

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
ast::Expr powModPolyGf(ast::Expr f, ast::Expr g, ast::Expr x, Int n, Int p, bool sym = true);


/**
 * @brief Compute f(x)^n in Zp[x]/g(x), that is  f(x)^n mod g(x) over Zp[x]
 *
 * @param f A polynomial expression in Zp[x]
 * @param g A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param n An integer
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @returnf(x)^n mod g(x) over Zp[x]
 */
ast::Expr powModPolyExprGf(ast::Expr f, ast::Expr g, ast::Expr x, Int n, Int p, bool sym = true);

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
ast::Expr randPolyGf(Int d, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Return a random polynomial expression over Zp[x]
 *
 * @param d The degree of the polynomial
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return A random polynomial with degree 'd' over Zp[x]
 */
ast::Expr randPolyExprGf(Int d, ast::Expr x, Int p, bool sym = true);

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
ast::Expr extendedEuclidGf(ast::Expr f, ast::Expr g, ast::Expr x, Int p, bool sym = true);

/**
 * @brief Extended Euclidean Algorithm over Zp[x]
 *
 * @param f A polynomial expression in Zp[x]
 * @param g A polynomial expression in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by default
 * @return A list with [gcd, s, t] such as s*f + g*t = gcd over Zp[x]
 */
ast::Expr extendedEuclidPolyExprGf(ast::Expr f, ast::Expr g, ast::Expr x, Int p, bool sym = true);


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
