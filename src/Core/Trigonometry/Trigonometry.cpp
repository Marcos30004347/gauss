#include "Trigonometry.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

namespace trigonometry {

AST* expandTrigRules(AST* A);
AST* contractTrigRules(AST* u);

AST* substituteTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->copy();

	AST* g = mapUnaryAST(u, substituteTrig);

	if(g->kind() == Kind::FunctionCall) {
		if(g->funName() == "tan") {
			AST* k = div(
				funCall("sin", { g->operand(0)->copy() }),
				funCall("cos", { g->operand(0)->copy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "cot") {
			AST* k = div(
				funCall("cos", { g->operand(0)->copy() }),
				funCall("sin", { g->operand(0)->copy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "sec") {
			AST* k = div(
				integer(1),
				funCall("cos", { g->operand(0)->copy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "csc") {
			AST* k = div(
				integer(1),
				funCall("sin", { g->operand(0)->copy() })
			);
			delete g;
			return k;
		}
	}

	return g;
}



signed long integer_fact(signed long i) {
	if(i == 1 || i == 0)
		return 1;

	return i * integer_fact(i - 1);
} 

AST* multipleAndlgeCos(AST* n, AST* theta) {
	AST* r = integer(0);

	for(int j=0; j <= n->value(); j++) {
		if(j%2 != 0)
			continue;
		
		signed long b = 
			integer_fact(n->value()) / 
			(integer_fact(j) * integer_fact(n->value() - j));
		
		AST* c_ = funCall("cos", {
			theta->copy()
		});

		AST* s_ = funCall("sin", {
			theta->copy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->copy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->copy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div(integer(j), integer(2))),
			integer(b),
			power(c, sub({ n->copy(), integer(j) })),
			power(s, integer(j))
		});

		r = add({r, e});

	}
	AST* d = reduceAST(r);
	delete r;
	return d;
}


AST* multipleAndlgeSin(AST* n, AST* theta) {
	AST* r = integer(0);

	for(int j=0; j <= n->value(); j++) {
		if(j%2 != 1)
			continue;
		
		signed long b = integer_fact(n->value())/(integer_fact(j) * integer_fact(n->value() - j));
		AST* c_ = funCall("cos", {
			theta->copy()
		});

		AST* s_ = funCall("sin", {
			theta->copy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->copy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->copy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div( sub({ integer(j), integer(1) }), integer(2))),
			integer(b),
			power(c, sub({ n->copy(), integer(j) })),
			power(s, integer(j))
		});
	
		r = add({r, e});
	}

	AST* expanded = reduceAST(r);

	delete r;

	return expanded;
}


AST* expandTrigRules(AST* A) {
	if(A->kind() == Kind::Addition) {
		AST* f = expandTrigRules(A->operand(0));
		
		AST* A__ = sub({
			A->copy(),
			A->operand(0)->copy()
		});
		
		AST* A_ = reduceAST(A__);

		AST* r = expandTrigRules(A_);

		AST* s_ = add({
			mul({
				f->operand(0)->copy(),
				r->operand(1)->copy(),
			}),
			mul({
				f->operand(1)->copy(),
				r->operand(0)->copy(),
			}),
		});

		AST* s = reduceAST(s_);

		AST* c_ = sub({
			mul({
				f->operand(1)->copy(),
				r->operand(1)->copy(),
			}),
			mul({
				f->operand(0)->copy(),
				r->operand(0)->copy(),
			}),
		});

		AST* c = reduceAST(s_);

		delete A_;
		delete A__;
		delete c_;
		delete s_;

		delete f;
		delete r;

		return list({s,c});
	}

	if(A->kind() == Kind::Multiplication) {
		AST* f = A->operand(0);
		if(f->kind() == Kind::Integer) {
			AST* k_ = div(
				A->copy(),
				f->copy()
			);

			AST* k = reduceAST(k_);
	
			AST* a = multipleAndlgeSin(f, k);
			AST* b = multipleAndlgeCos(f, k);
	

			delete k;
			delete k_;

			return list({a, b});
		}
	}

	return list({
		funCall("sin", { A->copy() }),
		funCall("cos", { A->copy() }),
	});
}

AST* expandTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->copy();

	AST* u_ = algebraicExpand(u);
	AST* v = mapUnaryAST(u_, expandTrig);
	delete u_;

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "sin"
	) {
		AST* a_ = expandTrigRules(v->operand(0));
		AST* a = reduceAST(a_);

		AST* r = a->operand(0)->copy();
		
		delete a;
		delete a_;
		delete v;

		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "cos"
	) {
		AST* a_ = expandTrigRules(v->operand(0));
		AST* a = reduceAST(a_);

		AST* r = a->operand(1)->copy();

		delete a;
		delete a_;
		delete v;

		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(isDivisionByZero(v)) {
		delete v;
		return undefined();
	}

	return v;
}

int floor(double x) {
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

AST* contractTrigPower(AST* u) {

	AST* b = u->operand(0);
	AST* n = u->operand(1);

	if(n->kind() == Kind::Integer && n->value() > 0) {
		if(
			b->kind() == Kind::FunctionCall && b->funName() == "cos"
		) {
			AST* theta = b->operand(0);
	
			if(n->value() % 2 == 0) {
				// n even
				AST* n_ = integer(n->value()/2);

				AST* p0 = div(
					div(
						integer(integer_fact(n->value())),
						integer(integer_fact(n_->value()) * integer_fact(n->value() - n_->value()))
					),
					power(integer(2), n->copy())
				);

				AST* p1 = div(
					integer(1),
					power(integer(2), sub({n->copy(), integer(1)}))
				);

				AST* p2 = integer(0);

				for(int j=0; j <= n_->value() - 1; j++) {
					AST* b = div(
						integer(integer_fact(n->value())),
						integer(integer_fact(j) * integer_fact(n->value() - j))
					);

					AST* c = funCall("cos", {
						mul({
							sub({
								n->copy(),
								mul({integer(2), integer(j)})
							}),
							theta->copy()
						})
					});

					p2 = add({ p2, mul({ b, c })});
				}

				delete n_;

				AST* r_ = add({
					p0,
					mul({p1, p2})
				});
				AST* r = reduceAST(r_);
	
				delete r_;

				return r;
			}


			if(n->value() % 2 == 1) {
				// n ood
				AST* n_ = integer(floor(n->value()/2.0f));

				AST* p1 = div(
					integer(1),
					power(integer(2), sub({n->copy(), integer(1)}))
				);

				AST* p2 = integer(0);
				for(int j=0; j<=n_->value() - 1; j++) {
					AST* b = div(
						integer(integer_fact(n->value())),
						integer(integer_fact(j) * integer_fact(n->value() - j))
					);
					AST* c = funCall("cos", {
						mul({
							sub({
								n->copy(),
								mul({integer(2), integer(j)})
							}),
							theta->copy()
						})
					});
					p2 = add({
						p2, mul({
							b, c
						})
					});
				}
				delete n_;

				AST* r_ = mul({p1, p2});

				AST* r = reduceAST(r_);
	
				delete r_;

				return r;
			}
		}
		if(
			b->kind() == Kind::FunctionCall && b->funName() == "sin" 
		) {
			
			AST* theta = b->operand(0);
	
			if(n->value() % 2 == 0) {
				// n even
				AST* n_ = integer(n->value()/2);

				AST* p0 = div(
					mul({
						power(integer(-1), n->copy()),
						div(
							integer(integer_fact(n->value())),
							integer(integer_fact(n->value()/2) * integer_fact(n->value() - (n->value()/2)))
						)
					}),
					power(integer(2), n->copy())
				);
	
				AST* p1 = div(
					power(integer(-1), integer(n->value()/2)),
					power(integer(2), sub({ n->copy(), integer(1) }))
				);

				AST* p2 = integer(0);
				for(int j=0; j<=n_->value() - 1; j++) {
					AST* a = power(
						integer(-1),
						integer(j)
					);
				
					AST* b = div(
						integer(integer_fact(n->value())),
						integer(integer_fact(j) * integer_fact(n->value() - j))
					);
				
					AST* c = funCall("cos", {
						mul({
							sub({
								n->copy(),
								mul({integer(2), integer(j)})
							}),
							theta->copy()
						})
					});
	
					p2 = add({
						p2, mul({
							a, b, c
						})
					});
				}
				delete n_;

				AST* r_ = add({
					p0,
					mul({p1, p2})
				});

				AST* r = reduceAST(r_);
	
				delete r_;

				return r;
			}

			if(n->value() % 2 == 1) {
				// n odd
				AST* n_ = integer(floor(n->value()/2.f));

				AST* p1 = div(
					power(integer(-1), integer((n->value() - 1)/2)),
					power(integer(2), sub({ n->copy(), integer(1) }))
				);
	
				AST* p2 = integer(0);
				for(int j=0; j<=n_->value() - 1; j++) {
					AST* a = power(
						integer(-1),
						integer(j)
					);
				
					AST* b = div(
						integer(integer_fact(n->value())),
						integer(integer_fact(j) * integer_fact(n->value() - j))
					);
				
					AST* c = funCall("sin", {
						mul({
							sub({
								n->copy(),
								mul({integer(2), integer(j)})
							}),
							theta->copy()
						})
					});
	
					p2 = add({
						p2, mul({
							b, a, c
						})
					});
				}

				delete n_;

				AST* r_ = mul({ p1, p2 });

				AST* r = reduceAST(r_);
	
				delete r_;

				return r;
			}
		}
	}
	return u->copy();
}

AST* separateSinCos(AST* u) {
	if(u->kind() == Kind::Multiplication) {
		AST* s = integer(1);
		AST* r = integer(1);

		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			AST* y = u->operand(i);
			
			if(
				(y->kind() == Kind::FunctionCall && y->funName() == "sin") ||
				(y->kind() == Kind::FunctionCall && y->funName() == "cos")
			) {
				s = mul({
					s,
					y->copy()
				});
			} else
			if(
				y->kind() == Kind::Power && 
				y->operand(0)->kind() == Kind::FunctionCall &&
				(
					y->operand(0)->funName() == "sin" ||
					y->operand(0)->funName() == "cos"
				) &&
				y->operand(1)->kind() == Kind::Integer &&
				y->operand(1)->value() > 0
			) {
				s = mul({
					s,
					y->copy()
				});
			} else {
				r = mul({
					r,
					y->copy()
				});
			}
		}
	
		AST* L = list({reduceAST(r), reduceAST(s)});
		
		delete r;
		delete s;
		
		return L;
	}

	if(
		(u->kind() == Kind::FunctionCall && u->funName() == "sin") ||
		(u->kind() == Kind::FunctionCall && u->funName() == "cos")
	) {
		return list({integer(1), u->copy()});
	}

	if(
		u->kind() == Kind::Power && 
		u->operand(0)->kind() == Kind::FunctionCall &&
		(
			u->operand(0)->funName() == "sin" ||
			u->operand(0)->funName() == "cos"
		) &&
		u->operand(1)->kind() == Kind::Integer &&
		u->operand(1)->value() > 0
	) {
		return list({integer(1), u->copy()});
	}

	return list({u->copy(), integer(1)});
}

AST* contractTrigProduct(AST* u) {

	if(u->kind() == Kind::Integer) {
		return u->copy();
	}

	if(u->numberOfOperands() == 1) {
		return u->operand(0)->copy();
	}

	if(u->numberOfOperands() == 2) {
		AST* A = u->operand(0);
		AST* B = u->operand(1);

		if(A->kind() == Kind::Power) {
			A = contractTrigPower(A);
			
			AST* C = mul({
				A->copy(),
				B->copy(),
			});
			AST* r = contractTrigRules(C);
		
			delete A;
			delete C;
		
			return r;
		}

		if(B->kind() == Kind::Power) {
			B = contractTrigPower(B);
			AST* C = mul({
				A->copy(),
				B->copy(),
			});
			AST* r = contractTrigRules(C);
		
			delete B;
			delete C;
		
			return r;
		}

		AST* theta 	= A->operand(0);
		AST* phi 		= B->operand(0);

		if(
			A->kind() == Kind::FunctionCall &&
			A->funName() == "sin" &&
			B->kind() == Kind::FunctionCall &&
			B->funName() == "sin"
		) {
			AST* t = sub({
				div(
					funCall("cos", { sub({theta->copy(), phi->copy()})}),
					integer(2)
				),
				div(
					funCall("cos", { add({theta->copy(), phi->copy()})}),
					integer(2)
				)
			});

			AST* r = reduceAST(t);

			delete t;
			return r;
		}

		if(
			A->kind() == Kind::FunctionCall &&
			A->funName() == "cos" &&
			B->kind() == Kind::FunctionCall &&
			B->funName() == "cos"
		) {
			AST* t = add({
				div(
					funCall("cos", { add({theta->copy(), phi->copy()})}),
					integer(2)
				),
				div(
					funCall("cos", { sub({theta->copy(), phi->copy()})}),
					integer(2)
				)
			});

			AST* r = reduceAST(t);

			delete t;
			return r;
		}

		if(
			A->kind() == Kind::FunctionCall &&
			A->funName() == "sin" &&
			B->kind() == Kind::FunctionCall &&
			B->funName() == "cos"
		) {
			AST* t = add({
				div(
					funCall("sin", { add({theta->copy(), phi->copy()})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({theta->copy(), phi->copy()})}),
					integer(2)
				)
			});

			AST* r = reduceAST(t);

			delete t;
			return r;
		}

		if(
			A->kind() == Kind::FunctionCall &&
			A->funName() == "cos" &&
			B->kind() == Kind::FunctionCall &&
			B->funName() == "sin"
		) {
			AST* t = add({
				div(
					funCall("sin", { add({theta->copy(), phi->copy()})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({phi->copy(), theta->copy()})}),
					integer(2)
				)
			});

			AST* r = reduceAST(t);

			delete t;
			return r;
		}
	}

	AST* A = u->operand(0);

	AST* k_ = div(
		u->copy(),
		A->copy()
	);

	AST* k = reduceAST(k_);

	delete k_;

	AST* B = contractTrigProduct(k);

	delete k;

	AST* r = mul({
		A->copy(),
		B->copy(),
	});

	AST* d = contractTrigRules(r);

	delete r;

	return d;
}

AST* contractTrigRules(AST* u) {

	AST* v = algebraicExpandRoot(u);

	if(v->kind() == Kind::Power) {
		AST* t = contractTrigPower(v);

		delete v;
		return t;
	}

	if(v->kind() == Kind::Multiplication) {
		AST* s = separateSinCos(v);
	

		AST* c = s->operand(0);
		AST* d = s->operand(1);
		
		if(d->kind() == Kind::Integer && d->value() == 1) {
			delete s;
			return v;
		}

		if(
			(d->kind() == Kind::FunctionCall && d->funName() == "sin") ||
			(d->kind() == Kind::FunctionCall && d->funName() == "cos")
		) {
			delete s;
			return v;
		}

		if(d->kind() == Kind::Power) {
			AST* k_ = mul({
				c->copy(),
				contractTrigPower(d)
			});

			AST* k = reduceAST(k_);
			AST* r = algebraicExpandRoot(k);
			
			delete s;
			delete v;
			delete k;
			delete k_;

			return r;
		}

		AST* k_ = mul({
			c->copy(),
			contractTrigProduct(d)
		});

		AST* k = reduceAST(k_);
		AST* r = algebraicExpandRoot(k);

		delete s;
		delete v;
		delete k;
		delete k_;

		return r;
	}

	if(v->kind() == Kind::Addition) {
		AST* s = integer(0);

		for(unsigned int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(
				y->kind() == Kind::Multiplication ||
				y->kind() == Kind::Power 
			) {

				s = add({
					s, contractTrigRules(y)
				});

			} else {
				s = add({
					s, y->copy()
				});
			}
		}

		AST* r = reduceAST(s);

		delete s;
		delete v;
		return r;
	}

	return v;
}

AST* contractTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->copy();

	AST* v_ = mapUnaryAST(u, contractTrig);
	AST* v = reduceAST(v_);
	delete v_;

	if(
		v->kind() == Kind::Multiplication ||
		v->kind() == Kind::Power
	) 
	{
		AST* t = contractTrigRules(v);
		delete v;
		return t;
	}

	return v;
}

}
