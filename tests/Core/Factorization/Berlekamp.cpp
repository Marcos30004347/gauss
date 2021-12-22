#include "test.hpp"

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Polynomial/Polynomial.hpp"
//#include "Core/Polynomial/Zp.hpp"
#include "Core/Factorization/Berlekamp.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_form_berkelamp_basis_matrix() {
  Expr x = Expr("x");

  Expr p = power(x, 4) + power(x, 2) + x + 1;
  Expr t = power(x, 8) + power(x, 7) + power(x, 4) + power(x, 3) + x + 1;
  Expr k = power(x, 5) + 6 * power(x, 4) + 4 * power(x, 3) + 4 * power(x, 2) +
           4 * x + 6;

  Expr r = 2;
  Expr q = 3;
  Expr z = 7;

  Expr A = buildBerkelampMatrix(p, x, r, false);

  assert(A[0].kind() == Kind::List);
  assert(A.size() == 4);

  assert(A[0].size() == 4);
  assert(A[0][0] == 1);
  assert(A[0][1] == 0);
  assert(A[0][2] == 0);
  assert(A[0][3] == 0);

  assert(A[1].size() == 4);
  assert(A[1][0] == 0);
  assert(A[1][1] == 0);
  assert(A[1][2] == 1);
  assert(A[1][3] == 0);

  assert(A[2].size() == 4);
  assert(A[2][0] == 1);
  assert(A[2][1] == 1);
  assert(A[2][2] == 1);
  assert(A[2][3] == 0);

  assert(A[3].size() == 4);
  assert(A[3][0] == 1);
  assert(A[3][1] == 1);
  assert(A[3][2] == 0);
  assert(A[3][3] == 1);

  Expr Ab = buildBerlekampBasis(A, r, false);

  assert(Ab.kind() == Kind::List);
  assert(Ab.size() == 2);

  assert(Ab[0].kind() == Kind::List);
  assert(Ab[0].size() == 4);
  assert(Ab[0][0] == 1);
  assert(Ab[0][1] == 0);
  assert(Ab[0][2] == 0);
  assert(Ab[0][3] == 0);

  assert(Ab[1].kind() == Kind::List);
  assert(Ab[1].size() == 4);
  assert(Ab[1][0] == 0);
  assert(Ab[1][1] == 0);
  assert(Ab[1][2] == 1);
  assert(Ab[1][3] == 1);

  Expr B = buildBerkelampMatrix(t, x, q, false);

  assert(B[0].kind() == Kind::List);
  assert(B.size() == 8);

  assert(B[0].size() == 8);
  assert(B[0][0] == 1);
  assert(B[0][1] == 0);
  assert(B[0][2] == 0);
  assert(B[0][3] == 0);
  assert(B[0][4] == 0);
  assert(B[0][5] == 0);
  assert(B[0][6] == 0);
  assert(B[0][7] == 0);

  assert(B[1].size() == 8);
  assert(B[1][0] == 0);
  assert(B[1][1] == 0);
  assert(B[1][2] == 0);
  assert(B[1][3] == 1);
  assert(B[1][4] == 0);
  assert(B[1][5] == 0);
  assert(B[1][6] == 0);
  assert(B[1][7] == 0);

  assert(B[2].size() == 8);
  assert(B[2][0] == 0);
  assert(B[2][1] == 0);
  assert(B[2][2] == 0);
  assert(B[2][3] == 0);
  assert(B[2][4] == 0);
  assert(B[2][5] == 0);
  assert(B[2][6] == 1);
  assert(B[2][7] == 0);

  assert(B[3].size() == 8);
  assert(B[3][0] == 1);
  assert(B[3][1] == 0);
  assert(B[3][2] == 2);
  assert(B[3][3] == 1);
  assert(B[3][4] == 0);
  assert(B[3][5] == 2);
  assert(B[3][6] == 0);
  assert(B[3][7] == 1);

  assert(B[4].size() == 8);
  assert(B[4][0] == 0);
  assert(B[4][1] == 1);
  assert(B[4][2] == 0);
  assert(B[4][3] == 0);
  assert(B[4][4] == 1);
  assert(B[4][5] == 2);
  assert(B[4][6] == 0);
  assert(B[4][7] == 0);

  assert(B[5].size() == 8);
  assert(B[5][0] == 1);
  assert(B[5][1] == 1);
  assert(B[5][2] == 0);
  assert(B[5][3] == 1);
  assert(B[5][4] == 2);
  assert(B[5][5] == 0);
  assert(B[5][6] == 0);
  assert(B[5][7] == 2);

  assert(B[6].size() == 8);
  assert(B[6][0] == 1);
  assert(B[6][1] == 0);
  assert(B[6][2] == 0);
  assert(B[6][3] == 0);
  assert(B[6][4] == 1);
  assert(B[6][5] == 0);
  assert(B[6][6] == 2);
  assert(B[6][7] == 0);

  assert(B[7].size() == 8);
  assert(B[7][0] == 2);
  assert(B[7][1] == 0);
  assert(B[7][2] == 1);
  assert(B[7][3] == 0);
  assert(B[7][4] == 0);
  assert(B[7][5] == 1);
  assert(B[7][6] == 0);
  assert(B[7][7] == 0);

  Expr Bb = buildBerlekampBasis(B, q, false);

  assert(Bb.kind() == Kind::List);
  assert(Bb.size() == 3);

  assert(Bb[0].kind() == Kind::List);
  assert(Bb[0].size() == 8);
  assert(Bb[0][0] == 1);
  assert(Bb[0][1] == 0);
  assert(Bb[0][2] == 0);
  assert(Bb[0][3] == 0);
  assert(Bb[0][4] == 0);
  assert(Bb[0][5] == 0);
  assert(Bb[0][6] == 0);
  assert(Bb[0][7] == 0);

  assert(Bb[2].kind() == Kind::List);
  assert(Bb[2].size() == 8);
  assert(Bb[2][0] == 0);
  assert(Bb[2][1] == 0);
  assert(Bb[2][2] == 0);
  assert(Bb[2][3] == 1);
  assert(Bb[2][4] == 0);
  assert(Bb[2][5] == 0);
  assert(Bb[2][6] == 0);
  assert(Bb[2][7] == 1);

  assert(Bb[1].kind() == Kind::List);
  assert(Bb[1].size() == 8);
  assert(Bb[1][0] == 0);
  assert(Bb[1][1] == 2);
  assert(Bb[1][2] == 2);
  assert(Bb[1][3] == 1);
  assert(Bb[1][4] == 1);
  assert(Bb[1][5] == 1);
  assert(Bb[1][6] == 1);
  assert(Bb[1][7] == 0);

  Expr C = factorization::buildBerkelampMatrix(k, x, z, false);

  assert(C[0].kind() == Kind::List);
  assert(C.size() == 5);

  assert(C[0].size() == 5);
  assert(C[0][0] == 1);
  assert(C[0][1] == 0);
  assert(C[0][2] == 0);
  assert(C[0][3] == 0);
  assert(C[0][4] == 0);

  assert(C[1].size() == 5);
  assert(C[1][0] == 4);
  assert(C[1][1] == 6);
  assert(C[1][2] == 2);
  assert(C[1][3] == 4);
  assert(C[1][4] == 3);

  assert(C[2].size() == 5);
  assert(C[2][0] == 2);
  assert(C[2][1] == 3);
  assert(C[2][2] == 6);
  assert(C[2][3] == 1);
  assert(C[2][4] == 4);

  assert(C[3].size() == 5);
  assert(C[3][0] == 6);
  assert(C[3][1] == 3);
  assert(C[3][2] == 5);
  assert(C[3][3] == 3);
  assert(C[3][4] == 1);

  assert(C[4].size() == 5);
  assert(C[4][0] == 1);
  assert(C[4][1] == 5);
  assert(C[4][2] == 5);
  assert(C[4][3] == 6);
  assert(C[4][4] == 6);

  Expr Cb = buildBerlekampBasis(C, z, false);

  assert(Cb.kind() == Kind::List);
  assert(Cb.size() == 2);

  assert(Cb[0].kind() == Kind::List);
  assert(Cb[0].size() == 5);
  assert(Cb[0][0] == 1);
  assert(Cb[0][1] == 0);
  assert(Cb[0][2] == 0);
  assert(Cb[0][3] == 0);
  assert(Cb[0][4] == 0);

  assert(Cb[1].kind() == Kind::List);
  assert(Cb[1].size() == 5);
  assert(Cb[1][0] == 0);
  assert(Cb[1][1] == 5);
  assert(Cb[1][2] == 6);
  assert(Cb[1][3] == 5);
  assert(Cb[1][4] == 1);
}

void should_factorize_square_free_poly_with_berlekamp() {
  Expr x = Expr("x");

  Expr ux = power(x, 6) + power(x, 5) + power(x, 4) + power(x, 3) + 1;

  Expr p0 = 2;

  Expr F0 = berlekampFactors(ux, x, p0, false);

  Expr f00 = power(x, 2) + x + 1;

  Expr f01 = power(x, 4) + x + 1;

  assert(F0.kind() == Kind::Set);
  assert(F0.size() == 3);
  assert(F0[0].kind() == Kind::Integer);
  assert(F0[0] == 1);
  assert(F0[1] == (f00));
  assert(F0[2] == (f01));

  Expr vx = 3 * power(x, 5) + 4 * power(x, 4) + 5 * power(x, 3) +
            5 * power(x, 2) + 5 * x + 4;

  Expr p1 = 7;

  Expr F1 = berlekampFactors(vx, x, p1, false);

  Expr f10 = 3;

  Expr f11 = x + 1;

  Expr f12 = power(x, 4) + 5 * power(x, 3) + 6 * power(x, 2) + 5 * x + 6;

  assert(F1.kind() == Kind::Set);
  assert(F1.size() == 3);
  assert(F1[0] == (f10));
  assert(F1[1] == (f11));
  assert(F1[2] == (f12));

  Expr tx = add({power(x, 4), power(x, 2), x, 1});

  Expr F2 = berlekampFactors(tx, x, p0, false);

  Expr f20 = x + 1;

  Expr f21 = power(x, 2) + power(x, 3) + 1;

  assert(F2.kind() == Kind::Set);
  assert(F2.size() == 3);
  assert(F2[0].kind() == Kind::Integer);
  assert(F2[0] == 1);
  assert(F2[1] == (f20));
  assert(F2[2] == (f21));

  Expr p3 = 3;

  Expr gx = power(x, 8) + power(x, 7) + power(x, 4) + power(x, 3) + x + 1;

  Expr F3 = berlekampFactors(gx, x, p3, false);

  Expr f30 = x + 1;

  Expr f31 = x + 2;

  Expr f32 = power(x, 6) + power(x, 5) + power(x, 4) + power(x, 3) +
             2 * power(x, 2) + 2 * x + 2;

  assert(F3.kind() == Kind::Set);
  assert(F3.size() == 4);
  assert(F3[0].kind() == Kind::Integer);
  assert(F3[0] == 1);
  assert(F3[1] == (f30));
  assert(F3[2] == (f31));
  assert(F3[3] == (f32));
}

int main() {
  TEST(should_form_berkelamp_basis_matrix)
  TEST(should_factorize_square_free_poly_with_berlekamp)
  return 0;
}
