#ifndef TYPECHECKER_HPP
#define TYPECHECKER_HPP

#include "AST.hpp"


union TypeContext { key key; };

struct TypeChecker;

TypeChecker* TypeCheckerCreate(ASTManager* manager);

void TypeCheckerDestroy(TypeChecker*, ASTManager*);

TypeContext TypeCheckerAddTypeBinding(TypeChecker* checker, TypeContext parent, ASTNodeKey term, ASTNodeKey type);

ASTNodeKey TypeCheckerGetTypeFromContext(TypeChecker* checker, TypeContext cxt, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateFunctionDeclarationType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateVariableDeclarationType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateCallType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateIntegerLiteralType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);
ASTNodeKey TypeCheckerDerivateBinaryOperationType(TypeChecker* checker, TypeContext ctx, ASTNodeKey term);

#endif
