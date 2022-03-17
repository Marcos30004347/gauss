/**
 * @file      Gauss.hpp
 * @brief     Header of the main API of this libary.
 * @date      Sun Mar  6 20:15:51 2022
 * @author    Marcos Vinicius Moreira Santos.
 * @copyright BSD-3-Clause
 *
 * This module defines all the public API methods of the Gauss library.
 */

#include "Algebra/Expression.hpp"
#include "gauss/Algebra/Matrix.hpp"

#include <array>
#include <cstddef>

namespace gauss {

	typedef alg::expr expr;
	typedef alg::kind kind;

namespace algebra {


/**
 * @brief Create a number expression from a double type.
 *
 * @details Create a number expression from a double type.
 * The returned value will be a fraction if the given value
 * has a fractional part greather than the machine episilon
 * or a integer if the given value does not have a fractional
 * part.
 *
 * @param[in] v A double value.
 *
 * @return The closest fraction to the double value
 * considering the maximum denominator specified.
 */
expr numberFromDouble(double v);

/**
 * @brief Create a number expression from a C string value.
 *
 * @details Create a number expression from a C string value.
 * The string should represent a Integer.
 *
 * @param[in] v A C string corresponding to a Integer.
 *
 * @return A Integer expression.
 */
expr intFromString(const char *v);

/**
 * @brief Create a number expression from a long value.
 *
 * @details Create a number expression from a long value.
 *
 * @param[in] v long value.
 *
 * @return A Integer expression.
 */
expr intFromLong(long v);

/**
 * @brief Creates a symbol expression.
 *
 * @details Creates a symbol expression.
 *
 * @param[in] v A C string corresponding
 * to the symbol identifier.
 *
 * @return A symbol expression.
 */
expr symbol(const char *v);

/**
 * @brief Create a power expression a^e.
 * @param[in] a The base of the power.
 * @param[in] e The expoent of the power.
 * @return A power expression.
 */
expr pow(expr a, expr e);

/**
 * @brief Return the ith operand of a.
 * @param[in] a A expression with multiple operands.
 * @param[in] i A positive integer number.
 * @return The i'th operant of a.
 */
expr &getOperand(expr a, size_t i);

/**
 * @brief Set the i'th operand of a to b.
 * @param[in] a A multi operand expression.
 * @param[in] i A positive integer number.
 * @param[in] b A expression.
 */
void setOperand(expr &a, size_t i, expr b);

/**
 * @brief Return the kind of a expression.
 * @param[in] a A expression.
 * @return The kind of the expression.
 */
kind kindOf(expr a);

/**
 * @brief Verify if a expression is of one of the given types.
 * @param[in] a A expression.
 * @param[in] k A integer that can be constructed from kinds
 * with bitwise 'or'. ex: kind::INT | kind::FRAC.
 * @return True if the expression is of one of the given kinds,
 * False otherwise.
 */
bool is(expr a, int k);

/**
 * @brief A expression of kind root.
 * @param[in] a The radical of the root expression.
 * @param[in] b The index of the root expression.
 * @return A root expression.
 */
expr root(expr a, expr b);

/**
 * @brief A square root expression.
 * @param[in] a The radical of the expression.
 * @return A square root expression.
 */
expr sqrt(expr a);

/**
 * @brief Creates an expression of form a + b.
 *
 * @details Creates an expression of form a + b, this
 * function does not evaluate the addition, the result
 * can be computed by a reduction, that is 'reduce(add(a, b))'.
 *
 * @param[in] a An algebraic expression.
 * @param[in] b An algebraic expression.
 * @return A new expression with the form a + b.
 */
expr add(expr a, expr b);

/**
 * @brief Creates an expression of form a - b.
 *
 * @details Creates an expression of form a - b, this
 * function does not evaluate the subtraction, the result
 * can be computed by a reduction, that is 'reduce(sub(a, b))'.
 *
 * @param[in] a An algebraic expression.
 * @param[in] b An algebraic expression.
 * @return A new expression with the form a - b.
 */
expr sub(expr a, expr b);

/**
 * @brief Creates an expression of form a * b;
 *
 * @details Creates an expression of form a * b, this
 * function does not evaluate the subtraction, the result
 * can be computed by a reduction, that is 'reduce(mul(a, b))'.
 *
 * @param[in] a An algebraic expression.
 * @param[in] b An algebraic expression.
 * @return A new expression with the form a * b.
 */
expr mul(expr a, expr b);

/**
 * @brief Creates an expression of form \begin{equation}\frac{a}{b}\end{equation}.
 * \f$ x=2 \f$
 * @details Creates an expression of form a / b, this
 * function does not evaluate the subtraction, the result
 * can be computed by a reduction, that is 'reduce(div(a, b))'.
 *
 * @param[in] a An algebraic expression.
 * @param[in] b An algebraic expression.
 * @return A new expression with the form a / b.
 */
expr div(expr a, expr b);

/**
 * @brief Expand a expression.
 *
 * @details Expand and reduce an expression.
 *
 * @example expand(x(3x + 4)) = 3x^2 + 4x.
 *
 * @return A algebraic expression corresponding to the
 * expansion of the expression 'a'.
 */
expr expand(expr a);

/**
 * @brief Reduce an expression.
 *
 * @details Reduce a expression to the smallest possible form
 * not regarding algebraix equalities or expansions. That means
 * that it perform the operations of a given expression.
 *
 * @param[in] a An algebraic expression.
 *
 * @example 'reduce(3x + 4y^2 + 5x + (3x + 3y^2)) = 11x + 7y^2'
 *
 * @return The reduced form of 'a'
 */
expr reduce(expr a);

/**
 * @brief Return a expression corresponding to a
 * call of the logarithmic function on 'x' with
 * a given base.
 *
 * @param[in] x An expression.
 * @param[in] base An expression.
 *
 * @return a call to the logarithmic function on x.
 */
expr log(expr x, expr base);

/**
 * @brief Return a expression corresponding to a
 * call of the eponential function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the exponential function on x.
 */
expr exp(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the natural logarithmic function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the natural logarithmic function on x.
 */
expr ln(expr x);

/**
 * @brief Replaces x on u by v.
 * @param[in] u A expression.
 * @param[in] x A symbol.
 * @param[in] v A expression.
 * @return u with all occurences of x replaced by v.
 */
expr replace(expr u, expr x, expr v);

/**
 * @brief Replaces x on u by v and expand the resulting expression.
 * @param[in] u A expression.
 * @param[in] x A symbol.
 * @param[in] v A expression.
 * @return A new expression without x.
 */
expr eval(expr u, expr x, expr v);

/**
 * @brief Return all free variables of the expression.
 * @param[in] u A expression.
 * @return A set of symbols.
 */
expr freeVariables(expr u);

namespace trigonometry {

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic sine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic sine function on x.
 */
expr sinh(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic cosine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic cosine function on x.
 */
expr cosh(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic tangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic tangent function on x.
 */
expr tanh(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the cosine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the cosine function on x.
 */
expr cos(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the sine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the sine function on x.
 */
expr sin(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the tangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the tangent function on x.
 */
expr tan(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the cosecant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the cosecant function on x.
 */
expr csc(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the cotangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the cotangent function on x.
 */
expr cot(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the secant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the secant function on x.
 */
expr sec(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic cotangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic cotangent function on x.
 */
expr coth(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic secant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic secant function on x.
 */
expr sech(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the hyperbolic cosecant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the hyperbolic cosecant function on x.
 */
expr csch(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc cosine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc cosine function on x.
 */
expr arccos(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc sine function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc sine function on x.
 */
expr arcsin(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc tangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc tangent function on x.
 */
expr arctan(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc cotangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc cotangent function on x.
 */
expr arccot(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc secant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc secant function on x.
 */
expr arcsec(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc cosecant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc cosecant function on x.
 */
expr arccsc(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc hyperbolic cosecant function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc hyperbolic cosecant function on x.
 */
expr arccosh(expr x);

/**
 * @brief Return a expression corresponding to a
 * call of the arc hyperbolic tangent function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the arc hyperbolic tangent function on x.
 */
expr arctanh(expr x);

}; // namespace trigonometry

namespace linear {

/**
 * @brief Creates a zero matrix.
 * @param[in] l Number of rows of the matrix.
 * @param[in] c Number of columns of the matrix.
 * @return A matrix filled with zeros.
 */
expr matrix(unsigned l, unsigned c);

/**
 * @brief Creates a identity matrix.
 * @param[in] l Number of rows of the matrix.
 * @param[in] c Number of columns of the matrix.
 * @return A identity matrix.
 */
expr identity(unsigned l, unsigned c);

/**
 * @brief Get a element of the matrix.
 * @param[in] A The matrix expression.
 * @param[in] i The row of the element.
 * @param[in] j The column of the element.
 * @return A number in the A[i][j] position;
 */
expr matrixGet(expr A, unsigned i, unsigned j);

/**
 * @brief Set a element of the matrix.
 * @param[in] A The matrix expression.
 * @param[in] i The row.
 * @param[in] j The column.
 * @param[in] a A double value.
 */
void matrixSet(expr A, unsigned i, unsigned j, double a);

/**
 * @brief Computes the singular value decomposition of a matrix;
 * @param[in] A The matrix expression.
 * @return A list with the matrices [U, D, transpose(V)].
 */
expr svd(expr A);

/**
 * @brief Return the inverse of a given matrix.
 * @param[in] A The matrix expression.
 * @return The inverse of the matrix A.
 */
expr inverse(expr A);

/**
 * @brief Computes the determinant of a matrix.
 * @param[in] A The matrix expression.
 * @return The determinant of the matrix.
 */
expr det(expr A);

/**
 * @brief Computes the tranposed form of a matrix.
 * @param[in] A The matrix expression.
 * @return The transposed form of 'A'.
 */
expr transpose(expr A);

/**
 * @brief Solve a linear system A*x = b.
 * @param[in] A Matrix of coefficients.
 * @param[in] b Vector of solutions.
 * @return The vector x.
 */
expr solveLinear(expr A, expr b);

} // namespace linear

} // namespace algebra

namespace polynomial {

/**
 * @brief Return the greatest degree of f on x;
 * @param[in] f A expression.
 * @param[in] x A symbol.
 * @return The degree greatest of f on x.
 */
expr degreePoly(expr f, expr x);

/**
 * @brief Return the coefficient of f on x^d;
 * @param[in] f A expression.
 * @param[in] x A symbol.
 * @param[in] d A integer.
 * @return The coefficient of f on x^d.
 */
expr coefficientPoly(expr f, expr x,
                              expr d);

/**
 * @brief Return the greatest coefficient of f on x;
 * @param[in] f A expression.
 * @param[in] x A symbol.
 * @return The greatest coefficient of f on x.
 */
expr leadingCoefficientPoly(expr f, expr x);

/**
 * @brief Computes the roots of a univariate polynomial
 * @details Computes the roots of the polynomial using
 * the Jenkins and Traub Algorithm.
 * @param[in] a Univariate Polynomial
 * @return A list with the roots of the polynomial.
 */
expr rootsOfPoly(expr a);

/**
 * @brief Computes the the content and the factors of a
 * Multivariate Polynomial.
 *
 * @param[in] poly A polynomial expression
 * @return The factorized form of a polynomial expression.
 */
expr factorPoly(expr poly);

/**
 * @brief Computes the resultant of two Polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The resultant polynomial of a and b.
 */
expr resultantOfPoly(expr a, expr b);

/**
 * @brief Add two polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The polynomial resulting of the addition of a and b.
 */
expr addPoly(expr a, expr b);

/**
 * @brief Subtract two polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The polynomial resulting of the subtraction of a and b.
 */
expr subPoly(expr a, expr b);

/**
 * @brief Multiply two polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The polynomial resulting of the multiplication of a and b.
 */
expr mulPoly(expr a, expr b);

/**
 * @brief Divide two polynomial expressions.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The expression quotient(a, b) + remainder(a, b)
 */
expr divPoly(expr a, expr b);

/**
 * @brief Compute the quotient of the polynomial division of a and b.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The expression corresponding to the quotient of 'a / b'
 */
expr quoPoly(expr a, expr b);

/**
 * @brief Compute the remainder of the polynomial division of a and b.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The expression corresponding to the remainder of 'a / b'
 */
expr remPoly(expr a, expr b);

/**
 * @brief Compute the greathest commom divisor of two polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The greathest commom divisor of 'a' and 'b'.
 */
expr gcdPoly(expr a, expr b);

/**
 * @brief Compute the least commom multiple of two polynomials.
 * @param[in] a A polynomial expression.
 * @param[in] b A polynomial expression.
 * @return The least commom multiple of 'a' and 'b'.
 */
expr lcmPoly(expr a, expr b);

namespace finiteField {

/**
 * @brief Compute a mod p;
 * @param[in] a A polynomial expression.
 * @param[in] p A long long integer.
 * @return a mod p
 */
expr projectPolyFiniteField(expr a, long long p);

/**
 * @brief Add two polynomial on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return (a + b) mod p
 */
expr addPolyFiniteField(expr a, expr b, long long p);

/**
 * @brief Subtract two polynomial on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return (a - b) mod p
 */
expr subPolyFiniteField(expr a, expr b, long long p);

/**
 * @brief Multiply two polynomial on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return (a * b) mod p
 */
expr mulPolyFiniteField(expr a, expr b, long long p);

/**
 * @brief Divide two polynomial on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return (a / b) mod p
 */
expr divPolyFiniteField(expr a, expr b, long long p);

/**
 * @brief Compute the quotient of a/b on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return quotient(a, b) mod p
 */
expr quoPolyFiniteField(expr a, expr b, long long p);

/**
 * @brief Compute the remainder of a/b on the finite field of length 'p';
 * @param[in] a A polynomial expression
 * @param[in] b A polynomial expression
 * @param[in] p A integer.
 * @return remainder(a, b) mod p
 */
expr remPolyFiniteField(expr a, expr b, long long p);
} // namespace finiteField

} // namespace polynomial

namespace calculus {

/**
 * @brief Compute the derivative of a algebraic expression
 *
 * @details Computes the derivative of the expression 'a' on the variable 'x
 * using elementary calculus methods.
 *
 * @param[in] a An algebraic expression
 * @param[in] x A free variable of the expression 'a'
 *
 *@return The algebraic expression corresponding to the derivative of 'a' by 'x'
 */
expr derivative(expr a, expr x);

} // namespace calculus

/**
 * @brief Return a string corresponding to a given expression.
 *
 * @param a A algebraic expression.
 *
 * @return A human friendly string representation of a given expression.
 */
std::string toString(expr a);

/**
 * @brief Construct a latex representation of a given expression.
 *
 * @param[in] a A expression.
 *
 * @param[in] useFractions If true, print rational numbers as fractions.
 *
 * @param[in] maxDenominators This is the maximum denominator for a fraction
 * representing a number between [0, 1], bigger the number, bigger the precision
 * on representing double precision floating points. Because of the nature of
 * floating arithmetic, you may not always want this number as big as it can be.
 *
 * @return A string representing a expression on latex format.
 */
std::string toLatex(expr a, bool print_as_fractions,
                    unsigned long max_den);

} // namespace gauss
