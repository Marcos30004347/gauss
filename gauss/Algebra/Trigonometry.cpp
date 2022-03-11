#include "Trigonometry.hpp"

namespace alg {

namespace trig {

expr sinh(expr x) { return func_call("sinh", {x}); }

expr cosh(expr x) { return func_call("cosh", {x}); }

expr tanh(expr x) { return func_call("tanh", {x}); }

expr cos(expr x) { return func_call("cos", {x}); }

expr sin(expr x) { return func_call("sin", {x}); }

expr tan(expr x) { return func_call("tan", {x}); }

expr csc(expr x) { return func_call("csc", {x}); }

expr cot(expr x) { return func_call("cot", {x}); }

expr sec(expr x) { return func_call("sec", {x}); }

expr coth(expr x) { return func_call("coth", {x}); }

expr sech(expr x) { return func_call("sech", {x}); }

expr csch(expr x) { return func_call("csch", {x}); }

expr arccos(expr x) { return func_call("arccos", {x}); }

expr arcsin(expr x) { return func_call("arcsin", {x}); }

expr arctan(expr x) { return func_call("arctan", {x}); }

expr arccot(expr x) { return func_call("arccot", {x}); }

expr arcsec(expr x) { return func_call("arcsec", {x}); }

expr arccsc(expr x) { return func_call("arccsc", {x}); }

expr arccosh(expr x) { return func_call("arccosh", {x}); }

expr arctanh(expr x) { return func_call("arctanh", {x}); }

}

}
