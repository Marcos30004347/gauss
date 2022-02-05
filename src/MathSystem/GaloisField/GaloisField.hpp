#ifndef GALOIS_FIELD_H
#define GALOIS_FIELD_H

// #include "MathSystem/AST/AST.hpp"
// #include "MathSystem/Algebra/List.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"

namespace galoisField {

/**
 * @brief Compute the representation of u(x) in Z[x] over Zp[x]
 *
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The representation of u(x) over Zp[x]
 */
// alg::expr gf(alg::expr u, alg::expr x, Int p, bool sym = true);

/**
 * @brief Compute the representation of u in Z[x] over Zp[x]
 *
 * @param u A polynomial in Z[x]
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The representation of u(x) over Zp[x]
 */
alg::expr gf(alg::expr u, Int p, bool symmetric);

/**
 * @brief Compute the representation of a polynomial expression in Z[x] over
 * Zp[x]
 *
 * @param u A polynomial expression in Z[x]
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The representation of u(x) over Zp[x]
 */
alg::expr gfPolyExpr(alg::expr u, Int p, bool symmetric);

/**
 * @brief Compute the representation of u(x) in Z[x] over Zp[x] without
 * expansion
 *
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @param p A prime integer x
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The representation of u(x) over Zp[x]
 */
alg::expr groundGf(alg::expr u, Int s, bool symmetric);

/**
 * @brief Divide a polynomial a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A list with que quotient and remainder over Zp[x]
 */
alg::expr divPolyGf(alg::expr a, alg::expr b, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Divide a polynomial expression a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A list with que quotient and remainder over Zp[x]
 */
alg::expr divPolyExprGf(alg::expr a, alg::expr b, alg::expr L, Int p,
                        bool sym = true);

/**
 * @brief Compute the remainder of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param X The x symbol
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The remainder of a(x) by b(x) over Zp[x]
 */
alg::expr remPolyGf(alg::expr a, alg::expr b, alg::expr X, Int p,
                    bool sym = true);

/**
 * @brief Compute the remainder of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The remainder of a(x) by b(x) over Zp[x]
 */
alg::expr remPolyExprGf(alg::expr a, alg::expr b, alg::expr L, Int p,
                        bool sym = true);

/**
 * @brief Compute the quotient of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The quotient of a(x) by b(x) over Zp[x]
 */
alg::expr quoPolyGf(alg::expr a, alg::expr b, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Compute the quotient of the division of a(x) by b(x) over Zp[x].
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The quotient of a(x) by b(x) over Zp[x]
 */
alg::expr quoPolyExprGf(alg::expr a, alg::expr b, alg::expr L, Int p,
                        bool sym = true);

/**
 * @brief Compute the monic form of f(x) over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A list with the content of f(x) and the monic form of f(x) over Zp[x]
 */
alg::expr monicPolyGf(alg::expr f, alg::expr x, Int p, bool sym = true);

/**
 * @brief Compute the monic form of f(x) over Zp[x]
 *
 * @param f A polynomial Expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A list with the content of f(x) and the monic form of f(x) over Zp[x]
 */
alg::expr monicPolyExprGf(alg::expr f, alg::expr L, Int p, bool sym = true);

/**
 * @brief Computes the gcd of a(x) and b(x) over Zp[x]
 *
 * @param a A polynomial in Zp[x]
 * @param b A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The gcd between a(x) and b(x) over Zp[x]
 */
alg::expr gcdPolyGf(alg::expr a, alg::expr b, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Computes the gcd of a(x) and b(x) over Zp[x]
 *
 * @param a A polynomial expression in Zp[x]
 * @param b A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The gcd between a(x) and b(x) over Zp[x]
 */
alg::expr gcdPolyExprGf(alg::expr a, alg::expr b, alg::expr L, Int p,
                        bool sym = true);

/**
 * @brief Add two polynomials over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(x) + g(x) over Zp[x]
 */
alg::expr addPolyGf(alg::expr f, alg::expr g, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Add two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(...) + g(...) over Zp
 */
alg::expr addPolyExprGf(alg::expr f, alg::expr g, Int p, bool sym = true);

/**
 * @brief Subtract two polynomials over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(x) - g(x) over Zp[x]
 */
alg::expr subPolyGf(alg::expr f, alg::expr g, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Sub two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(...) - g(...) over Zp
 */
alg::expr subPolyExprGf(alg::expr f, alg::expr g, Int p, bool sym = true);

/**
 * @brief Multiply two polynomials over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(x) * g(x) over Zp[x]
 */
alg::expr mulPolyGf(alg::expr f, alg::expr g, alg::expr x, Int p,
                    bool sym = true);

/**
 * @brief Multiply two polynomials expressions over Zp
 *
 * @param f A polynomial expression in Zp
 * @param g A polynomial expression in Zp
 * @param p A prime integer
 * @param sym  true if Zp is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return The result of f(...) * g(...) over Zp
 */
alg::expr mulPolyExprGf(alg::expr f, alg::expr g, Int p, bool sym = true);

/**
 * @brief Compute f(x)^n in Zp[x]/g(x), that is  f(x)^n mod g(x) over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param n An integer
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @returnf(x)^n mod g(x) over Zp[x]
 */
alg::expr powModPolyGf(alg::expr f, alg::expr g, alg::expr x, Int n, Int p,
                       bool sym = true);

/**
 * @brief Compute f(x)^n in Zp[x]/g(x), that is  f(x)^n mod g(x) over Zp[x]
 *
 * @param f A polynomial expression in Zp[x]
 * @param g A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param n An integer
 * @param p A prime integer
 * @param sym  true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @returnf(x)^n mod g(x) over Zp[x]
 */
alg::expr powModPolyExprGf(alg::expr f, alg::expr g, alg::expr L, Int n, Int p,
                           bool sym = true);

/**
 * @brief Return a random polynomial over Zp[x]
 *
 * @param d The degree of the polynomial
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A random polynomial with degree 'd' over Zp[x]
 */
alg::expr randPolyGf(Int d, alg::expr x, Int p, bool sym = true);

/**
 * @brief Return a random polynomial expression over Zp[x]
 *
 * @param d The degree of the polynomial
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A random polynomial with degree 'd' over Zp[x]
 */
alg::expr randPolyExprGf(Int d, alg::expr L, Int p, bool sym = true);

/**
 * @brief Extended Euclidean Algorithm over Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param g A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
 * @return A list with [gcd, s, t] such as s*f + g*t = gcd over Zp[x]
 */
alg::expr extendedEuclidGf(alg::expr f, alg::expr g, alg::expr x, Int p,
                           bool sym = true);

/**
 * @brief Extended Euclidean Algorithm over Zp[x]
 *
 * @param f A polynomial expression in Zp[x]
 * @param g A polynomial expression in Zp[x]
 * @param L The list of symbols in a, this list can have at most one element
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * @return A list with [gcd, s, t] such as s*f + g*t = gcd over Zp[x]
 */
alg::expr extendedEuclidPolyExprGf(alg::expr f, alg::expr g, alg::expr L, Int p,
                                   bool sym = true);

/**
 * @brief Return the quotient of the poly expr u divided by the constant v on Zp
 *
 * @param u A poly expr
 * @param v a Integer
 * @param p a Integer
 * @return The poly expr u/v on Zp
 */
alg::expr groundQuoPolyExprGf(alg::expr u, Int v, Int p, bool symmetric = true);

/**
 * @brief Computes the quotient of the division of s by t in Zp[x]
 *
 * @param s An integer
 * @param t An integer
 * @param p A prime integer
 * @param sym true if Zp[x] is in symmetric representation
 * @return The quotient of s / t in Zp[x]
 */
Int quoGf(Int s, Int t, Int p, bool sym = true);

/**
 * @brief Computes the inverse of a in Zp
 *
 * @param a A integer
 * @param p A prime number
 * @param sym true if Zp[x] is in symmetric representation
 * 						or false otherwise, true by
 * default
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

} // namespace galoisField

#endif
