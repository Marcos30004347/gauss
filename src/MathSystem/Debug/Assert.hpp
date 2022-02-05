#ifndef CORE_DEBUG_H
#define CORE_DEBUG_H

#include <iostream>
#include <stdarg.h>
		
#ifndef NDEBUG
#   define assert(Expr, ...) \
    __assert(#Expr, Expr, __FILE__, __LINE__, __VA_ARGS__)
#else
#   define assert(Expr, Msg, ...) ;
#endif

void __assert(const char* expr_str, bool expr, const char* file, int line, const char* msg, ...);

#endif
