
#include "Calculus.hpp"
#include "gauss/Algebra/Expression.hpp"

using namespace alg;

namespace calc {

// expr derivative(expr u, expr x) {
// 	return expr(
// 		Kind::Derivative,
// 		{ u, x }
// 	);
// }

expr derivateInverseTrig(expr u, expr x) {
  if (is(&u, kind::POW) && degree(u) == -1) {
    if (u[0].kind() == kind::FUNC) {
      if (u[0].funName() == "sin") {
        return reduce((1 / pow(1 - pow(u[0][0], 2), fraction(1,2))) *
                      derivate(u[0], x));
      }

      if (u[0].funName() == "cos") {
        return reduce(-1 * (1 / pow((1 - pow(u[0][0], 2)), fraction(1,2))) *
                      derivate(u[0], x));
      }

      if (u[0].funName() == "tan") {
        return reduce((1 / (1 + pow(u[0][0], 2))) * derivate(u[0], x));
      }

      if (u[0].funName() == "cot") {
        return reduce(-1 * (1 / (1 + pow(u[0][0], 2))) * derivate(u[0], x));
      }

      if (u[0].funName() == "sec") {
        return reduce(
            (1 / (abs(u[0][0]) * pow((pow(u[0][0], integer(2)) - 1), fraction(1,2)))) *
            derivate(u[0], x));
      }

      if (u[0].funName() == "csc") {
        return reduce(-1 *
                      (1 / (abs(u[0][0]) * pow((pow(u[0][0], 2) - 1), fraction(1,2)))) *
                      derivate(u[0], x));
      }
    }
  }

  return undefined();
}

expr derivatePow(expr u, expr x) {
  if (u.kind() == kind::POW) {
    expr v = base(u);
    expr w = degree(u);

    expr d = ((degree(u) * pow(base(u), degree(u) - 1) * derivate(v, x)) +
              (derivate(w, x) * pow(base(u), degree(u)) * ln(base(u))));

    return reduce(d);
  }

  return undefined();
}

expr derivateSumsAndSubs(expr u, expr x) {
  if (u.kind() == kind::ADD || u.kind() == kind::SUB) {
    expr dx = 0;

    for (unsigned int i = 0; i < u.size(); i++) {
      dx = dx + derivate(u[i], x);
    }

    return reduce(dx);
	}

  return undefined();
}

expr derivateMul(expr u, expr x) {
  if (u.kind() == kind::MUL) {
    expr v = u[0];
    expr w = reduce(u / v);

    return reduce((derivate(v, x) * w) + (v * derivate(w, x)));
  }

  return undefined();
}

expr derivateFuncs(expr u, expr x) {
  if (u.kind() == kind::FUNC) {

    if (u.funName() == "abs") {
      return abs(derivate(u[0], x));
    }

    if (u.funName() == "ln") {
      return reduce(1 / u[0]);
    }

    if (u.funName() == "log") {
      if (u.size() == 2) {
        return reduce(1 / (u[0] * ln(u[1])));
      } else {
        return reduce(1 / (u[0] * ln(2)));
      }
    }

    if (u.funName() == "exp") {
      return reduce(exp(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "tan") {
      return reduce(pow(sec(u[0]), 2) * derivate(u[0], x));
    }

    if (u.funName() == "sinh") {
      return reduce(cosh(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "cosh") {
      return reduce(sinh(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "tanh") {
      return reduce(pow(sech(u[0]), 2) * derivate(u[0], x));
    }

    if (u.funName() == "sec") {
      return reduce(sec(u[0]) * tan(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "csc") {
      return reduce(-1 * cot(u[0]) * csc(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "cot") {
      return reduce(-1 * pow(csc(u[0]), 2) * derivate(u[0], x));
    }

    if (u.funName() == "coth") {
      return reduce(-1 * pow(csch(u[0]), 2) * derivate(u[0], x));
    }

    if (u.funName() == "sech") {
      return reduce(-1 * tanh(u[0]) * sech(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "csch") {
      return reduce(-1 * coth(u[0]) * csch(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "sin") {
      return reduce(cos(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "cos") {
      return reduce(-1 * sin(u[0]) * derivate(u[0], x));
    }

    if (u.funName() == "arcsin") {
      return reduce((1 / pow(1 - pow(u[0], 2), fraction(1,2))) * derivate(u[0], x));
    }

    if (u.funName() == "arccos") {
      return reduce(-1 * (1 / pow((1 - pow(u[0], 2)), fraction(1,2))) *
                    derivate(u[0], x));
    }

    if (u.funName() == "arctan") {
      return reduce((1 / (1 + pow(u[0], 2))) * derivate(u[0], x));
    }

    if (u.funName() == "arccot") {
      return reduce(-1 * (1 / (1 + pow(u[0], 2))) * derivate(u[0], x));
    }

    if (u.funName() == "arcsec") {
      return reduce((1 / (abs(u[0]) * pow((pow(u[0], 2) - 1), fraction(1, 2)))) *
                    derivate(u[0], x));
    }

    if (u.funName() == "arccsc") {
      return reduce(-1 * (1 / (abs(u[0]) * pow((pow(u[0], 2) - 1), fraction(1, 2)))) *
                    derivate(u[0], x));
    }

    if (u.funName() == "arccosh") {
      return reduce((1 / pow((pow(u[0], 2) - 1), fraction(1,2))) * derivate(u[0], x));
    }

    if (u.funName() == "arctanh") {
      return reduce((1 / (1 - pow(u[0], 2))) * derivate(u[0], x));
    }
  }

  return undefined();
}

expr derivate(expr u, expr x) {
	expr dx;

  if (u == x) {
    return 1;
  }

	// TODO: add verification if u is a inverse trig or pow, or mul, or etc
	// before doing the differentiation

  dx = derivateInverseTrig(u, x);
  if (dx != undefined()) {
    return dx;
  }

  dx = derivatePow(u, x);
  if (dx != undefined()) {
    return dx;
  }

  dx = derivateSumsAndSubs(u, x);
  if (dx != undefined()) {
    return dx;
  }

  dx = derivateMul(u, x);
  if (dx != undefined()) {
    return dx;
  }

  dx = derivateFuncs(u, x);
  if (dx != undefined()) {
    return dx;
  }

  if (u.freeOf(x)) {
    return 0;
  }

	return diff(u, x);
}

} // namespace calculus
