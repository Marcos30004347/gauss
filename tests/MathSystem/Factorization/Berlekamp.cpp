#include "MathSystem/Algebra/Expression.hpp"
#include "test.hpp"

#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/Factorization/Berlekamp.hpp"

using namespace alg;
using namespace polynomial;
using namespace factorization;

void should_form_berkelamp_basis_matrix() {
  expr x = expr("x");

  expr p = pow(x, 4) + pow(x, 2) + x + 1;
  expr t = pow(x, 8) + pow(x, 7) + pow(x, 4) + pow(x, 3) + x + 1;
  expr k = pow(x, 5) + 6 * pow(x, 4) + 4 * pow(x, 3) + 4 * pow(x, 2) +
           4 * x + 6;

  expr r = 2;
  expr q = 3;
  expr z = 7;

  expr A = buildBerkelampMatrix(p, x, r, false);

	assert(A[0].kind() == kind::LIST);

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
  expr Ab = buildBerlekampBasis(A, r, false);

  assert(Ab.kind() == kind::LIST);
  assert(Ab.size() == 2);

  assert(Ab[0].kind() == kind::LIST);
  assert(Ab[0].size() == 4);
  assert(Ab[0][0] == 1);
  assert(Ab[0][1] == 0);
  assert(Ab[0][2] == 0);
  assert(Ab[0][3] == 0);

  assert(Ab[1].kind() == kind::LIST);
  assert(Ab[1].size() == 4);
  assert(Ab[1][0] == 0);
  assert(Ab[1][1] == 0);
  assert(Ab[1][2] == 1);
  assert(Ab[1][3] == 1);

  expr B = buildBerkelampMatrix(t, x, q, false);

  assert(B[0].kind() == kind::LIST);
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

  expr Bb = buildBerlekampBasis(B, q, false);

  assert(Bb.kind() == kind::LIST);
  assert(Bb.size() == 3);

  assert(Bb[0].kind() == kind::LIST);
  assert(Bb[0].size() == 8);
  assert(Bb[0][0] == 1);
  assert(Bb[0][1] == 0);
  assert(Bb[0][2] == 0);
  assert(Bb[0][3] == 0);
  assert(Bb[0][4] == 0);
  assert(Bb[0][5] == 0);
  assert(Bb[0][6] == 0);
  assert(Bb[0][7] == 0);

  assert(Bb[2].kind() == kind::LIST);
  assert(Bb[2].size() == 8);
  assert(Bb[2][0] == 0);
  assert(Bb[2][1] == 0);
  assert(Bb[2][2] == 0);
  assert(Bb[2][3] == 1);
  assert(Bb[2][4] == 0);
  assert(Bb[2][5] == 0);
  assert(Bb[2][6] == 0);
  assert(Bb[2][7] == 1);

  assert(Bb[1].kind() == kind::LIST);
  assert(Bb[1].size() == 8);
  assert(Bb[1][0] == 0);
  assert(Bb[1][1] == 2);
  assert(Bb[1][2] == 2);
  assert(Bb[1][3] == 1);
  assert(Bb[1][4] == 1);
  assert(Bb[1][5] == 1);
  assert(Bb[1][6] == 1);
  assert(Bb[1][7] == 0);

  expr C = buildBerkelampMatrix(k, x, z, false);

  assert(C[0].kind() == kind::LIST);
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

  expr Cb = buildBerlekampBasis(C, z, false);

  assert(Cb.kind() == kind::LIST);
  assert(Cb.size() == 2);

  assert(Cb[0].kind() == kind::LIST);
  assert(Cb[0].size() == 5);
  assert(Cb[0][0] == 1);
  assert(Cb[0][1] == 0);
  assert(Cb[0][2] == 0);
  assert(Cb[0][3] == 0);
  assert(Cb[0][4] == 0);

  assert(Cb[1].kind() == kind::LIST);
  assert(Cb[1].size() == 5);
  assert(Cb[1][0] == 0);
  assert(Cb[1][1] == 5);
  assert(Cb[1][2] == 6);
  assert(Cb[1][3] == 5);
  assert(Cb[1][4] == 1);
}

void should_factorize_square_free_poly_with_berlekamp() {
  expr x = expr("x");

  expr ux = pow(x, 6) + pow(x, 5) + pow(x, 4) + pow(x, 3) + 1;

  expr p0 = 2;

  expr F0 = berlekampFactors(ux, x, p0, false);

  expr f00 = pow(x, 2) + x + 1;

  expr f01 = pow(x, 4) + x + 1;

  assert(F0.kind() == kind::LIST);
  assert(F0.size() == 3);
  assert(F0[0].kind() == kind::INT);
  assert(F0[0] == 1);
  assert(F0[1] == (f00));
  assert(F0[2] == (f01));

  expr vx = 3 * pow(x, 5) + 4 * pow(x, 4) + 5 * pow(x, 3) +
            5 * pow(x, 2) + 5 * x + 4;

  expr p1 = 7;

  expr F1 = berlekampFactors(vx, x, p1, false);

  expr f10 = 3;

  expr f11 = x + 1;

  expr f12 = pow(x, 4) + 5 * pow(x, 3) + 6 * pow(x, 2) + 5 * x + 6;

  assert(F1.kind() == kind::LIST);
  assert(F1.size() == 3);
  assert(F1[0] == (f10));
  assert(F1[1] == (f11));
  assert(F1[2] == (f12));

  expr tx = create(kind::ADD, {pow(x, 4), pow(x, 2), x, 1});

  expr F2 = berlekampFactors(tx, x, p0, false);

  expr f20 = x + 1;

  expr f21 = pow(x, 2) + pow(x, 3) + 1;

  assert(F2.kind() == kind::LIST);
  assert(F2.size() == 3);
  assert(F2[0].kind() == kind::INT);
  assert(F2[0] == 1);
  assert(F2[1] == (f20));
  assert(F2[2] == (f21));

  expr p3 = 3;

  expr gx = pow(x, 8) + pow(x, 7) + pow(x, 4) + pow(x, 3) + x + 1;

  expr F3 = berlekampFactors(gx, x, p3, false);

  expr f30 = x + 1;

  expr f31 = x + 2;

  expr f32 = pow(x, 6) + pow(x, 5) + pow(x, 4) + pow(x, 3) +
             2 * pow(x, 2) + 2 * x + 2;

  assert(F3.kind() == kind::LIST);
  assert(F3.size() == 4);
  assert(F3[0].kind() == kind::INT);
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
