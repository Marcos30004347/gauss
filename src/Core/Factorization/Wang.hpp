#ifndef FACTORIZATION_WANG_H
#define FACTORIZATION_WANG_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief This is a helper function that is used to test
 * 				evaluation point suitability.
 * 
 * @param G An integer that multiplies F[0]*...*F[i]
 * 
 * @param F The list of factors F[i] that will be tested.
 *
 * @param c The content of U[0](x)
 *
 * @param L A list containing the variablex x[i]
 * 
 * @param K The field K[L...] of F[i](L...)
 * 
 * @param d The set of integers d[i] defined as above,
 * 					this array must be pre-allocated before
 * 					the call.
 * 
 *  @return 'true' if all the d[i] where computed or 'false' in case of failure. 
 *
 */
bool nondivisors(long G, ast::AST* F, long c, ast::AST* L, ast::AST* K, long* d);

// bool getEvaluationPoint(ast::AST* F, ast::AST* L, ast::AST* K, long mod = 3);
ast::AST* factors(ast::AST* f, ast::AST* L, ast::AST* K);
ast::AST* factorsWang(ast::AST* f, ast::AST* L, ast::AST* K);
ast::AST* trialDivision(ast::AST* f, ast::AST* F, ast::AST* L, ast::AST* K);
ast::AST* groundLeadCoeff(ast::AST* f, ast::AST* L);

}

#endif
