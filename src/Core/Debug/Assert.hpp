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

void __assert(const char* expr_str, bool expr, const char* file, int line, const char* msg, ...)
{
	if (!expr)
	{
		va_list lst;
		va_start(lst, msg);

		std::cerr << "Assert failed: ";

		while(*msg != '\0') {
				if(*msg!='%')
				{
						putchar(*msg);
						msg++;
						continue;
				}
				msg++;
				switch(*msg)
				{
						case 's': fputs(va_arg(lst, char *), stdout); break;
						case 'c': putchar(va_arg(lst, int)); break;
				}
				msg++;
		}

		std::cerr	<< "\nExpected: '" << expr_str << "'\n"
		<< "Source: '" << file << "', line: '" << line << "'\n";
		abort();
	}
}

#endif
