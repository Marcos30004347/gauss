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
	) return u->deepCopy();

	AST* g = mapUnaryAST(u, substituteTrig);

	if(g->kind() == Kind::FunctionCall) {
		
		if(g->funName() == "tan") {
			AST* k = div(
				funCall("sin", { g->operand(0)->deepCopy() }),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "cot") {
			AST* k = div(
				funCall("cos", { g->operand(0)->deepCopy() }),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "sec") {
			AST* k = div(
				integer(1),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "csc") {
			AST* k = div(
				integer(1),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}
	}

	return g;
}

bool isDivisionByZero(AST* k) {
	AST* d = denominator(k);
	
	if(d->kind() == Kind::Integer && d->value() == 0) {
		delete d;
		return true;
	}

	delete d;
	return false;
}

AST* expandExponentialRules(AST* A) {
	if(A->kind() == Kind::Addition) {
		AST* f = A->operand(0);
	
		AST* k = sub({
			A->deepCopy(),
			f->deepCopy()
		});

		AST* a_ = funCall(
			"exp", {
				f->deepCopy()
			}
		);

		AST* b_ = funCall(
			"exp", {
				reduceAST(k)
			}
		);

		AST* r_ =  mul({
			expandExponential(a_),
			expandExponential(b_),
		});
	
		AST* r = reduceAST(r_);
	
		delete k;
		delete a_;
		delete b_;
		delete r_;
	
		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(A->kind() == Kind::Multiplication) {
		AST* f = A->operand(0);

		if(f->kind() == Kind::Integer || f->kind() == Kind::Symbol) {
			AST* k = div(A->deepCopy(), f->deepCopy());
			
			AST* p_ = power(
				funCall(
					"exp", {
						reduceAST(k)
					}
				),
				f->deepCopy()
			);

			AST* p = reduceAST(p_);

			delete p_;
			delete k;

			if(isDivisionByZero(p)) {
				delete p;
				return undefined();
			}
	
			return p;
		}
	}

	return funCall("exp", { A->deepCopy() });
}

AST* expandExponential(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* u_ = algebraicExpand(u);
	
	AST* v = mapUnaryAST(u_, expandExponential);
	
	delete u_;

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "exp"
	) {
		AST* r = expandExponentialRules(v->operand(0));
		delete v;
		return r;
	}

	if(isDivisionByZero(v)) {
		delete v;
		return undefined();
	}

	return v;
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
			theta->deepCopy()
		});

		AST* s_ = funCall("sin", {
			theta->deepCopy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->deepCopy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->deepCopy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div(integer(j), integer(2))),
			integer(b),
			power(c, sub({ n->deepCopy(), integer(j) })),
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
			theta->deepCopy()
		});

		AST* s_ = funCall("sin", {
			theta->deepCopy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->deepCopy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->deepCopy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div( sub({ integer(j), integer(1) }), integer(2))),
			integer(b),
			power(c, sub({ n->deepCopy(), integer(j) })),
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
			A->deepCopy(),
			A->operand(0)->deepCopy()
		});
		
		AST* A_ = reduceAST(A__);

		AST* r = expandTrigRules(A_);

		AST* s_ = add({
			mul({
				f->operand(0)->deepCopy(),
				r->operand(1)->deepCopy(),
			}),
			mul({
				f->operand(1)->deepCopy(),
				r->operand(0)->deepCopy(),
			}),
		});

		AST* s = reduceAST(s_);

		AST* c_ = sub({
			mul({
				f->operand(1)->deepCopy(),
				r->operand(1)->deepCopy(),
			}),
			mul({
				f->operand(0)->deepCopy(),
				r->operand(0)->deepCopy(),
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
				A->deepCopy(),
				f->deepCopy()
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
		funCall("sin", { A->deepCopy() }),
		funCall("cos", { A->deepCopy() }),
	});
}

AST* expandTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* u_ = algebraicExpand(u);
	AST* v = mapUnaryAST(u_, expandExponential);
	delete u_;

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "sin"
	) {
		AST* a_ = expandTrigRules(v->operand(0));
		AST* a = reduceAST(a_);

		AST* r = a->operand(0)->deepCopy();
		
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

		AST* r = a->operand(1)->deepCopy();

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

// AST* expandMainProduct(AST* r, AST* s) {
// 	if(r->kind() == Kind::Addition) {
// 		AST* f = r->operand(0);

// 		AST* v = sub({
// 			r->deepCopy(),
// 			f->deepCopy(),
// 		});

// 		AST* k = reduceAST(v);
	
// 		AST* z = add({
// 			expandMainProduct(f, s),
// 			expandMainProduct(k, s),
// 		});

// 		AST* y = reduceAST(z);

// 		delete v;
// 		delete k;
// 		delete z;

// 		return y;
// 	}

// 	if(s->kind() == Kind::Addition) {
// 		return expandMainProduct(s, r);
// 	}

// 	AST* a = algebraicExpand(r);
// 	AST* b = algebraicExpand(s);

// 	if(a->kind() == Kind::Addition || b->kind() == Kind::Addition) {
// 		AST* t = expandMainProduct(a, b);
		
// 		delete a;
// 		delete b;
		
// 		return t;
// 	}

// 	AST* t = mul({ a, b });

// 	AST* k = reduceAST(t);

// 	delete t;

// 	return k;
// }

// AST* expandMainPower(AST* u, AST* n) {
// 	if(u->kind() == Kind::Addition) {
// 		AST* f = u->operand(0);

// 		AST* r_ = sub({ u->deepCopy(), f->deepCopy() });

// 		AST* r = reduceAST(r_);
		
// 		delete r_;

// 		AST* s = integer(0);
	
// 		for(int k_ = 0; k_ <= n->value(); k_++) {
// 			AST* k = integer(k_);

// 			AST* c_ = div(
// 				integer(integer_fact(n->value())),
// 				integer(integer_fact(k->value()) * integer_fact(n->value() - k->value()))
// 			);
	
// 			AST* c = reduceAST(c_);
	
	
// 			AST* z_ = mul({
// 				c->deepCopy(),
// 				power(
// 					f->deepCopy(),
// 					integer(n->value() - k->value())
// 				)
// 			});

// 			AST* z = reduceAST(z_);


// 			AST* t = expandMainPower(r, k);
		
// 			s = add({ s, expandMainProduct(z, t) });

// 			delete c;
// 			delete k;
// 			delete z;
// 			delete t;
// 			delete z_;
// 			delete c_;
// 		}
		
// 		delete r;

// 		return s;
// 	}
	
// 	AST* v = power(
// 		u->deepCopy(),
// 		n->deepCopy()
// 	);

// 	AST* reduced = reduceAST(v);
// 	delete v;

// 	return reduced;
// }

// AST* expandMainOperator(AST* u) {
// 	if(u->isTerminal())
// 		return reduceAST(u);

// 	AST* u_ = reduceAST(u);
	
// 	// Maybe this should be removed
// 	if(u_->kind() == Kind::Addition) {
// 		AST* v = u_->operand(0);
	
// 		AST* a = sub({
// 			u_->deepCopy(),
// 			v->deepCopy()
// 		});

// 		AST* k = reduceAST(a);
	
// 		AST* t = add({
// 			algebraicExpand(v),
// 			algebraicExpand(k)
// 		});

// 		delete u_;
// 		u_ = reduceAST(t);
		
// 		delete k;
// 		delete a;
// 		delete t;
// 	}

// 	if(u_->kind() == Kind::Multiplication) {

// 		AST* v = u_->operand(0);
// 		AST* e = div(
// 			u_->deepCopy(),
// 			v->deepCopy()
// 		);
		
// 		AST* t = reduceAST(e);

// 		AST* z = expandMainProduct(t, v);

// 		delete u_;
// 		u_ = reduceAST(z);

// 		delete t;
// 		delete e;
// 		delete z;

// 	}

// 	if(u_->kind() == Kind::Power) {

// 		AST* b = u_->operand(0)->deepCopy();
// 		AST* e = u_->operand(1)->deepCopy();

// 		if(e->kind() == Kind::Integer && e->value() >= 2) {
// 			AST* t = expandMainPower(b, e);

// 			delete u_;
// 			u_ = reduceAST(t);
// 			delete t;
// 		}
	
// 		if(e->kind() == Kind::Integer && e->value() <= -2) {
// 			AST* p_ = power(u_->deepCopy(), integer(-1));
// 			AST* p = reduceAST(p_);

// 			delete p_;
			
// 			AST* b_ = p->operand(0);
// 			AST* e_ = p->operand(1);
	
// 			AST* t = expandMainPower(b_, e_);
			
// 			delete p;
			
// 			delete u_;
// 			u_ = power(t, integer(-1));
// 		}

// 		delete b;
// 		delete e;
// 	}

	
// 	AST* k = reduceAST(u_);

// 	delete u_;

// 	return k;
// }


AST* contractExponentialRules(AST* u) {
	AST* v = algebraicExpandRoot(u);

	if(v->kind() == Kind::Power) {
		AST* b = v->operand(0);
		AST* s = v->operand(1);

		if(b->kind() == Kind::FunctionCall && b->funName() == "exp") {
			AST* p = mul({
				b->operand(0)->deepCopy(),
				s->deepCopy()
			});
	
			if(
				p->kind() == Kind::Multiplication || 
				p->kind() == Kind::Power
			) {
				AST* p_ = contractExponentialRules(p);
				delete p;
				p = p_;
				delete v;

				AST* r = funCall("exp", { p });
				return r;
			}

			delete p;
			return v;
		}
	}

	if(v->kind() == Kind::Multiplication) {
		AST* p = integer(1);
		AST* s = integer(0);

		for(int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(y->kind() == Kind::FunctionCall && y->funName() == "exp") {
				s = add({
					s, y->operand(0)->deepCopy()
				});
			} else {
				p = mul({
					p,
					y->deepCopy()
				});
			}
		}
	
		AST* r_;
	
		if(s->kind() == Kind::Integer && s->value() == 0) {
			delete s;
			r_ = p;
		} else {
			r_ = mul({ funCall("exp", {s}), p });
		}
	
		AST* r = reduceAST(r_);
		
		delete v;
		delete r_;
		
		return r;
	}

	if(v->kind() == Kind::Addition) {
		AST* s = integer(0);
		for(int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(
				y->kind() == Kind::Multiplication ||
				y->kind() == Kind::Power 
			) {
				s = add({
					s, contractExponentialRules(y)
				});
			} else {
				s = add({
					s,
					y->deepCopy()
				});
			}
		}

		AST* r = reduceAST(s);

		delete v;
		delete s;
		
		return r;
	}

	return v;
}

AST* contractExponential(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* v_ = mapUnaryAST(u, contractExponential);
	AST* v = reduceAST(v_);
	delete v_;

	if(
		v->kind() == Multiplication ||
		v->kind() == Power
	) {
		AST* t = contractExponentialRules(v);
		delete v;
		return t;
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
					power(integer(2), n->deepCopy())
				);

				AST* p1 = div(
					integer(1),
					power(integer(2), sub({n->deepCopy(), integer(1)}))
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
								n->deepCopy(),
								mul({integer(2), integer(j)})
							}),
							theta->deepCopy()
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
					power(integer(2), sub({n->deepCopy(), integer(1)}))
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
								n->deepCopy(),
								mul({integer(2), integer(j)})
							}),
							theta->deepCopy()
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
						power(integer(-1), n->deepCopy()),
						div(
							integer(integer_fact(n->value())),
							integer(integer_fact(n->value()/2) * integer_fact(n->value() - (n->value()/2)))
						)
					}),
					power(integer(2), n->deepCopy())
				);
	
				AST* p1 = div(
					power(integer(-1), integer(n->value()/2)),
					power(integer(2), sub({ n->deepCopy(), integer(1) }))
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
								n->deepCopy(),
								mul({integer(2), integer(j)})
							}),
							theta->deepCopy()
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
					power(integer(2), sub({ n->deepCopy(), integer(1) }))
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
								n->deepCopy(),
								mul({integer(2), integer(j)})
							}),
							theta->deepCopy()
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
	return u->deepCopy();
}

AST* separateSinCos(AST* u) {
	if(u->kind() == Kind::Multiplication) {
		AST* s = integer(1);
		AST* r = integer(1);

		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* y = u->operand(i);
			
			if(
				y->kind() == Kind::FunctionCall && y->funName() == "sin" ||
				y->kind() == Kind::FunctionCall && y->funName() == "cos"
			) {
				s = mul({
					s,
					y->deepCopy()
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
					y->deepCopy()
				});
			} else {
				r = mul({
					r,
					y->deepCopy()
				});
			}
		}
	
		AST* L = list({reduceAST(r), reduceAST(s)});
		
		delete r;
		delete s;
		
		return L;
	}

	if(
		u->kind() == Kind::FunctionCall && u->funName() == "sin" ||
		u->kind() == Kind::FunctionCall && u->funName() == "cos"
	) {
		return list({integer(1), u->deepCopy()});
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
		return list({integer(1), u->deepCopy()});
	}

	return list({u->deepCopy(), integer(1)});
}

AST* contractTrigProduct(AST* u) {

	if(u->kind() == Kind::Integer) {
		return u->deepCopy();
	}

	if(u->numberOfOperands() == 1) {
		return u->operand(0)->deepCopy();
	}

	if(u->numberOfOperands() == 2) {
		AST* A = u->operand(0);
		AST* B = u->operand(1);

		if(A->kind() == Kind::Power) {
			A = contractTrigPower(A);
			
			AST* C = mul({
				A->deepCopy(),
				B->deepCopy(),
			});
			AST* r = contractTrigRules(C);
		
			delete A;
			delete C;
		
			return r;
		}

		if(B->kind() == Kind::Power) {
			B = contractTrigPower(B);
			AST* C = mul({
				A->deepCopy(),
				B->deepCopy(),
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
					funCall("cos", { sub({theta->deepCopy(), phi->deepCopy()})}),
					integer(2)
				),
				div(
					funCall("cos", { add({theta->deepCopy(), phi->deepCopy()})}),
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
					funCall("cos", { add({theta->deepCopy(), phi->deepCopy()})}),
					integer(2)
				),
				div(
					funCall("cos", { sub({theta->deepCopy(), phi->deepCopy()})}),
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
					funCall("sin", { add({theta->deepCopy(), phi->deepCopy()})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({theta->deepCopy(), phi->deepCopy()})}),
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
					funCall("sin", { add({theta->deepCopy(), phi->deepCopy()})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({phi->deepCopy(), theta->deepCopy()})}),
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
		u->deepCopy(),
		A->deepCopy()
	);

	AST* k = reduceAST(k_);

	delete k_;

	AST* B = contractTrigProduct(k);

	delete k;

	AST* r = mul({
		A->deepCopy(),
		B->deepCopy(),
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
			d->kind() == Kind::FunctionCall && d->funName() == "sin" ||
			d->kind() == Kind::FunctionCall && d->funName() == "cos"
		) {
			delete s;
			return v;
		}

		if(d->kind() == Kind::Power) {
			AST* k_ = mul({
				c->deepCopy(),
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
			c->deepCopy(),
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
		// printf("asdasd\n");
		AST* s = integer(0);
		for(int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(
				y->kind() == Kind::Multiplication ||
				y->kind() == Kind::Power 
			) {

				s = add({
					s, contractTrigRules(y)
				});

				// printf("uuuuuu %s\n", s->toString().c_str());

	
			} else {
				s = add({
					s, y->deepCopy()
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
	) return u->deepCopy();

	AST* v_ = mapUnaryAST(u, contractTrig);
	AST* v = reduceAST(v_);
	delete v_;

	if(
		v->kind() == Kind::Multiplication ||
		v->kind() == Kind::Power
	) {
		AST* t = contractTrigRules(v);
		delete v;
		return t;
	}

	return v;
}

}
