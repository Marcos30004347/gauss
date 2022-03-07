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
#include "MathSystem/Algebra/Matrix.hpp"

namespace gauss {
namespace algebra {

typedef alg::expr expr;

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
 * @param[optional] max_den The maximum denominator of the
 * fraction returned.
 *
 * @return The closest fraction to the double value
 * considering the maximum denominator specified.
 */
expr number(double v, long max_den = 1000000000);

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
expr number(const char *v);

/**
 * @brief Create a number expression from a long value.
 *
 * @details Create a number expression from a long value.
 *
 * @param[in] v long value.
 *
 * @return A Integer expression.
 */
expr number(long long v);

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
 * @brief Creates an expression of form a / b;
 *
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
 * @example 'expand((3*x + 4) * x) = 3x^2 + 4x'.
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
 * call of the logarithmic function on 'x'.
 *
 * @param[in] x An expression.
 *
 * @return a call to the logarithmic function on x.
 */
expr log(expr x);

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
typedef alg::matrix mat;

mat addMat(mat a, mat b);

mat subMat(mat a, mat b);

mat mulMat(mat a, mat b);

mat mulMat(mat a, long b);

mat divMat(mat a, mat b);

mat divMat(mat a, long b);

mat svdMat(mat a);

mat invMat(mat a);

mat echelonFormMat(mat a);

mat transpMat(mat a);

mat solveSystem(mat a, mat b);

}; // namespace linear

}; // namespace algebra

namespace polynomial {

algebra::expr rootsOfPoly(algebra::expr a);

algebra::expr factorsOfPoly(algebra::expr poly);

algebra::expr resultantOfPoly(algebra::expr a, algebra::expr b);

algebra::expr addPoly(algebra::expr a, algebra::expr b);

algebra::expr subPoly(algebra::expr a, algebra::expr b);

algebra::expr mulPoly(algebra::expr a, algebra::expr b);

algebra::expr divPoly(algebra::expr a, algebra::expr b);

algebra::expr quoPoly(algebra::expr a, algebra::expr b);

algebra::expr remPoly(algebra::expr a, algebra::expr b);

algebra::expr gcdPoly(algebra::expr a, algebra::expr b);

algebra::expr lcmPoly(algebra::expr a, algebra::expr b);

namespace finiteField {
algebra::expr addPolyFF(algebra::expr a, algebra::expr b, Int p);

algebra::expr subPolyFF(algebra::expr a, algebra::expr b, Int p);

algebra::expr mulPolyFF(algebra::expr a, algebra::expr b, Int p);

algebra::expr divPolyFF(algebra::expr a, algebra::expr b, Int p);

algebra::expr quoPolyFF(algebra::expr a, algebra::expr b, Int p);

algebra::expr remPolyFF(algebra::expr a, algebra::expr b, Int p);
}; // namespace finiteField

}; // namespace polynomial

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
algebra::expr derivative(algebra::expr a, algebra::expr x);

}; // namespace calculus

/**
 * @brief Return a string corresponding to a given expression.
 *
 * @param a A algebraic expression.
 *
 * @return A human friendly string representation of a given expression.
 */
std::string toString(algebra::expr a);

/**
 * @brief Return a string corresponding to a given expression.
 *
 * @param a A algebraic expression.
 *
 * @return A string in latex representation of a given expression.
 */
std::string toLatex(algebra::expr a);

}; // namespace gauss
