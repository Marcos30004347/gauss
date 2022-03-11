#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "gauss/Algebra/Matrix.hpp"

#include "gauss/SVD/SVD.hpp"


using namespace alg;

void echelonFormTests()
{
	matrix A = matrix(3, 4, {
		1, 2, -1, -4,
		2, 3, -1, -11,
		-2, 0, -3, 22
	});
	matrix A_echelon = echelonForm(A);

	assert(A_echelon[0][0] == 1.f);
	assert(A_echelon[1][0] == 0.f);
	assert(A_echelon[2][0] == 0.f);

	assert(A_echelon[0][1] == 0.f);
	assert(A_echelon[1][1] == 1.f);
	assert(A_echelon[2][1] == 0.f);

	assert(A_echelon[0][2] == 0.f);
	assert(A_echelon[1][2] == 0.f);
	assert(A_echelon[2][2] == 1.f);

	assert(A_echelon[0][3] == -8.f);
	assert(A_echelon[1][3] == 1.f);
	assert(A_echelon[2][3] == -2.f);
}

void nullspaceTests()
{
	matrix A = matrix(3, 4, {
		1, 2, 3, 4,
		1, 3, 5, 6,
		2, 5, 8, 10
	});

	matrix A_space = nullspace(A);

	assert(A_space[0][0] == 1);
	assert(A_space[0][1] == -2);
	assert(A_space[0][2] == 1);
	assert(A_space[0][3] == 0);

	assert(A_space[1][0] == 0);
	assert(A_space[1][1] == -2);
	assert(A_space[1][2] == 0);
	assert(A_space[1][3] == 1);

	matrix B = matrix(2, 4, {
		-1, 1, 2, 4,
		2, 0, 1, -7,
	});

	matrix B_space = nullspace(B);

	assert(B_space[0][0] == -0.5);
	assert(B_space[0][1] == -2.5);
	assert(B_space[0][2] == 1);
	assert(B_space[0][3] == 0);

	assert(B_space[1][0] == 3.5);
	assert(B_space[1][1] == -0.5);
	assert(B_space[1][2] == 0);
	assert(B_space[1][3] == 1);

	matrix C = matrix(2, 2, {
		2, 1,
		1, 2,
	});

	matrix C_space = nullspace(C);

	assert(C_space[0][0] == 0);
	assert(C_space[0][1] == 0);
}

void svdTests()
{
	matrix at, a, res;

	a = matrix(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});

	std::tuple<matrix, matrix, matrix> SVD1 = svd(a);

	res = std::get<0>(SVD1) * std::get<1>(SVD1) * std::get<2>(SVD1);

	for(unsigned i=0; i< a.lines(); i++)
		for(unsigned j=0; j<a.columns(); j++) {
			assert(fabs(a[i][j] - res[i][j]) <= 2.22e-14);
		}

	at = transpose(a);

	std::tuple<matrix, matrix, matrix> SVD2 =  svd(at);

	res = std::get<0>(SVD2) * std::get<1>(SVD2) * std::get<2>(SVD2);

	for(unsigned i=0; i<at.lines(); i++)
		for(unsigned j=0; j<at.columns(); j++)
			assert(fabs(at[i][j] - res[i][j]) <= 2.22e-14);
}


int main()
{
	matrix A(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}, 2, 2);

	int v = 1;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			assert(A[i][j] == v++);

	matrix B(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(B[i][j] == j+1);


	matrix C(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	matrix D(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	matrix E = C + D;

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(E[i][j] == C[i][j] + D[i][j]);

	matrix F(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 2, 2);

	matrix G = C + F;

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(G[i][j] == C[i][j] + F[i][j]);

	matrix H = C*D;

	float data[] = {45, 90, 135, 180, 225, 270, 315, 360, 405};

	for(int i=0; i<9; i++)
	 	for(int j=0; j<9; j++)
	 		assert(H[i][j] == data[j]);

	matrix I = C*F;


  for(int i=0; i<9; i++) {
  	for(int j=0; j<9; j++) {
  		assert(I[i][j] == data[j]);
		}
	}

	matrix J(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	matrix K(9, 1, {
		1, 2, 3, 4, 5, 6, 7, 8 ,9
	}, 4, 1);

	matrix L = J*K;

	assert(L.lines() == 9);
	assert(L.columns() == 1);

	for(int i=0; i<9; i++)
		assert(L[i][0] == 285);

  matrix M = transpose(J);

  assert(M.lines() == J.columns());
  assert(M.columns() == J.lines());

	for(unsigned int i=0; i<M.lines(); i++)
	 	for(unsigned int j=0; j<M.columns(); j++)
  		assert(M[i][j] == J[j][i]);

	matrix N = transpose(L);
  assert(N.lines() == L.columns());
  assert(N.columns() == L.lines());

	for(unsigned int i=0; i<N.lines(); i++)
		for(unsigned int j=0; j<N.columns(); j++)
			assert(N[i][j] == L[j][i]);

	matrix P = N*L;

	assert(P.columns() == 1);
	assert(P.lines() == 1);
	assert(P[0][0] == 731025);

	matrix Q = M*3.f;

	for(unsigned int i=0; i<Q.lines(); i++)
  	for(unsigned int j=0; j<Q.columns(); j++)
  		assert(Q[i][j] == 3*M[i][j]);

	matrix W = Q/3.f;

	for(unsigned int i=0; i<Q.lines(); i++)
	 	for(unsigned int j=0; j<Q.columns(); j++)
	 		assert(W[i][j] == Q[i][j]/3.f);

	matrix X = matrix(1, 1, { 3 });

  matrix Z = Q/X;

	for(unsigned int i=0; i<Q.lines(); i++)
  	for(unsigned int j=0; j<Q.columns(); j++)
	 		assert(Z[i][j] == Q[i][j]/3.f);

	matrix Y = W*X;

	for(unsigned int i=0; i<Q.lines(); i++)
  	for(unsigned int j=0; j<Q.columns(); j++)
  		assert(Y[i][j] == W[i][j]*3.f);

	matrix T(3,3,{
		1,1,0,
		2,1,3,
		3,1,1
	});

	std::pair<matrix, matrix> T_LU = LUDecomposition(T);

	assert(T_LU.first[0][0] == 1);
	assert(T_LU.first[0][1] == 0);
	assert(T_LU.first[0][2] == 0);
	assert(T_LU.first[1][0] == 2);
	assert(T_LU.first[1][1] == -1);
	assert(T_LU.first[1][2] == 0);
	assert(T_LU.first[2][0] == 3);
	assert(T_LU.first[2][1] == -2);
	assert(T_LU.first[2][2] == -5);

	assert(T_LU.second[0][0] == 1);
	assert(T_LU.second[0][1] == 1);
	assert(T_LU.second[0][2] == 0);
	assert(T_LU.second[1][0] == 0);
	assert(T_LU.second[1][1] == 1);
	assert(T_LU.second[1][2] == -3);
	assert(T_LU.second[2][0] == 0);
	assert(T_LU.second[2][1] == 0);
	assert(T_LU.second[2][2] == 1);

	matrix S(2,2,{
		3,2,
		2,6,
	});

	// LU decomposition tests
	std::pair<matrix, matrix> System_Permutation = LUPDecomposition(S);

	matrix b(2, 1, {2,-8});

	matrix x = LUPSolve(System_Permutation.first, System_Permutation.second, b);

	assert(x[0][0] == 2);

	assert(x[1][0] == -2);

	matrix Inv = LUPInverse(System_Permutation.first, System_Permutation.second);

	matrix Res = S*Inv;

	assert(Res[0][0] - 1 < 0.00000001);
	assert(Res[0][1] - 0 < 0.00000001);
	assert(Res[1][0] - 0 < 0.00000001);
	assert(Res[1][1] - 1 < 0.00000001);

	double Det = LUPDeterminant(System_Permutation.first, System_Permutation.second);

	assert(Det == 14);


	echelonFormTests();

  nullspaceTests();

	svdTests();

	return 0;
}
