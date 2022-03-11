#include "Matrix.hpp"
#include "MathSystem/Algebra/Expression.hpp"
#include "MathSystem/SVD/SVD.hpp"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdio.h>
#include <string>
#include <tuple>

using namespace alg;

matrix::matrix() {
  this->m_data = nullptr;
  this->m_lines = 0;
  this->m_columns = 0;
  this->m_block_heigth = 0;
  this->m_block_width = 0;
  this->m_stored_lines = 0;
  this->m_stored_column = 0;
}

matrix::matrix(const matrix &other) {

  this->m_lines = other.m_lines;
  this->m_columns = other.m_columns;
  this->m_block_heigth = other.m_block_heigth;
  this->m_block_width = other.m_block_width;
  this->m_stored_lines = other.m_stored_lines;
  this->m_stored_column = other.m_stored_column;

  this->m_data = new double[this->storedLines() * this->storedColumns()];

  std::copy(other.m_data,
            other.m_data + (other.storedLines() * other.storedColumns()),
            this->m_data);
}

matrix *matrix::copy(matrix *o) {
  matrix *t = new matrix();

  t->m_lines = o->m_lines;
  t->m_columns = o->m_columns;
  t->m_block_heigth = o->m_block_heigth;
  t->m_block_width = o->m_block_width;
  t->m_stored_column = o->m_stored_column;
  t->m_stored_lines = o->m_stored_lines;

  t->m_data = new double[o->m_stored_lines * o->m_stored_column];

  std::copy(o->m_data, o->m_data + (o->storedLines() * o->storedColumns()),
            t->m_data);

  return t;
}

matrix::matrix(unsigned int l, unsigned int c, double *data)
    : matrix{l, c, 1, 1, data} {}

matrix::matrix(unsigned int l, unsigned int c) {
  m_lines = l;
  m_columns = c;
  m_block_heigth = 1;
  m_block_width = 1;
  this->m_stored_lines =
      std::ceil(this->lines() / (double)this->m_block_heigth) *
      this->m_block_heigth;
  this->m_stored_column =
      std::ceil(this->columns() / (double)this->blockWidth()) *
      this->m_block_width;
  this->m_data = new double[this->storedColumns() * this->storedLines()];
  std::fill(this->m_data,
            this->m_data + this->storedColumns() * this->storedLines(), 0);
}

// matrix::matrix(unsigned int l, unsigned int c,
//                std::initializer_list<double> data, unsigned int bx,
//                unsigned int by)
//     : matrix{l, c, bx, by, data} {

// }
// matrix::matrix(unsigned int l, unsigned int c, unsigned int bx, unsigned int
// by)
//     : matrix{l, c, bx, by} {

// }

matrix::matrix(unsigned int lines, unsigned int columns,
               std::initializer_list<double> data, unsigned int block_width,
               unsigned int block_heigth)
    : matrix{lines, columns, block_width, block_heigth, data} {}

matrix::matrix(unsigned int lines, unsigned int columns,
               unsigned int block_width, unsigned int block_heigth,
               std::initializer_list<double> data) {
  m_lines = lines;
  m_columns = columns;
  m_block_heigth = block_heigth;
  m_block_width = block_width;
  this->m_stored_lines =
      std::ceil(this->lines() / (double)this->m_block_heigth) *
      this->m_block_heigth;
  this->m_stored_column =
      std::ceil(this->columns() / (double)this->blockWidth()) *
      this->m_block_width;
  this->m_data = new double[this->storedColumns() * this->storedLines()];
  std::fill(this->m_data,
            this->m_data + this->storedColumns() * this->storedLines(), 0);

  // define extra space as identity
  for (unsigned int i = 0;
       i < std::min(this->storedLines() - this->lines(),
                    this->storedColumns() - this->columns());
       i++) {
    int block_y = (this->lines() + i) / this->m_block_heigth;
    int block_x = (this->columns() + i) / this->blockWidth();

    int y = (this->lines() + i) - block_y * this->m_block_heigth;
    int x = (this->columns() + i) - block_x * this->blockWidth();

    this->m_data[this->stride(block_y, block_x) + y * this->blockWidth() + x] =
        1.f;
  }

  for (unsigned int block_y = 0;
       block_y < this->storedLines() / this->m_block_heigth; block_y++) {
    for (unsigned int block_x = 0;
         block_x < this->storedColumns() / this->blockWidth(); block_x++) {
      // Current block width and current block heigth
      unsigned int row_margin = lines - block_y * this->m_block_heigth;
      unsigned int col_margin = columns - block_x * this->blockWidth();

      int width = std::min(m_block_width, col_margin);
      int heigth = std::min(m_block_heigth, row_margin);

      // where current block starts

      unsigned int idx = this->m_block_heigth * this->columns() * block_y +
                         block_x * this->blockWidth();

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          this->m_data[this->stride(block_y, block_x) + y * this->blockWidth() +
                       x] = *(data.begin() + idx++);
        }
        // this->storedLines() - this->lines() + m_block_width - width;
        idx += this->columns() - (block_x + 1) * this->blockWidth() +
               (block_x) * this->m_block_width + m_block_width - width;
      }
    }
  }
}

matrix::matrix(unsigned int lines, unsigned int columns,
               unsigned int block_width, unsigned int block_heigth,
               double *data) {
  m_lines = lines;
  m_columns = columns;
  m_block_heigth = block_heigth;
  m_block_width = block_width;

  this->m_stored_lines =
      std::ceil(this->lines() / (double)this->m_block_heigth) *
      this->m_block_heigth;
  this->m_stored_column =
      std::ceil(this->columns() / (double)this->blockWidth()) *
      this->m_block_width;

  this->m_data = new double[this->storedLines() * this->storedColumns()];
  std::fill(this->m_data,
            this->m_data + this->storedColumns() * this->storedLines(), 0);

  // define extra space as identity
  for (unsigned int i = 0;
       i < std::min(this->storedLines() - this->lines(),
                    this->storedColumns() - this->columns());
       i++) {
    int block_y = (this->lines() + i) / this->m_block_heigth;
    int block_x = (this->columns() + i) / this->blockWidth();

    int y = (this->lines() + i) - block_y * this->m_block_heigth;
    int x = (this->columns() + i) - block_x * this->blockWidth();

    this->m_data[this->stride(block_y, block_x) + y * this->blockWidth() + x] =
        1.f;
  }

  for (unsigned int block_y = 0;
       block_y < this->storedLines() / this->m_block_heigth; block_y++) {
    for (unsigned int block_x = 0;
         block_x < this->storedColumns() / this->blockWidth(); block_x++) {
      // Current block width and current block heigth
      unsigned int row_margin = lines - block_y * this->m_block_heigth;
      unsigned int col_margin = columns - block_x * this->blockWidth();

      int width = std::min(m_block_width, col_margin);
      int heigth = std::min(m_block_heigth, row_margin);

      // unsigned int idx = block_y*this->m_block_heigth*this->blockWidth() +
      // block_x*this->m_block_width;
      unsigned int idx = this->m_block_heigth * this->columns() * block_y +
                         block_x * this->blockWidth();

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          this->m_data[this->stride(block_y, block_x) + y * this->blockWidth() +
                       x] = data[idx];
        }
        idx += this->columns() - (block_x + 1) * this->blockWidth() +
               (block_x) * this->m_block_width + m_block_width - width;
      }
    }
  }
}

matrix::matrix(unsigned int lines, unsigned int columns,
               unsigned int block_width, unsigned int block_heigth) {
  m_lines = lines;
  m_columns = columns;
  m_block_heigth = block_heigth;
  m_block_width = block_width;

  this->m_stored_lines =
      std::ceil(this->lines() / (double)this->m_block_heigth) *
      this->m_block_heigth;
  this->m_stored_column =
      std::ceil(this->columns() / (double)this->blockWidth()) *
      this->m_block_width;

  this->m_data = new double[this->storedLines() * this->storedColumns()];
  std::fill(this->m_data,
            this->m_data + this->storedColumns() * this->storedLines(), 0);

  // define extra space as identity
  for (unsigned int i = 0;
       i < std::min(this->storedLines() - this->lines(),
                    this->storedColumns() - this->columns());
       i++) {
    int block_y = (this->lines() + i) / this->m_block_heigth;
    int block_x = (this->columns() + i) / this->blockWidth();

    int y = (this->lines() + i) - block_y * this->m_block_heigth;
    int x = (this->columns() + i) - block_x * this->blockWidth();

    this->m_data[this->stride(block_y, block_x) + y * this->blockWidth() + x] =
        1.f;
  }
}

matrix::~matrix() {
  if (this->m_data) {
    delete[] this->m_data;
  }
  this->m_data = 0;
}

unsigned int matrix::lines() const { return this->m_lines; }

unsigned int matrix::columns() const { return this->m_columns; }

unsigned int matrix::blockWidth() const { return this->m_block_width; }

unsigned int matrix::blockHeight() const { return this->m_block_heigth; }

// #include <string.h>

matrix &matrix::operator=(const matrix &other) {
  this->m_stored_column = other.m_stored_column;
  this->m_stored_lines = other.m_stored_lines;
  this->m_lines = other.m_lines;
  this->m_columns = other.m_columns;
  this->m_block_width = other.m_block_width;
  this->m_block_heigth = other.m_block_heigth;

  if (this->m_data)
    delete[] this->m_data;

  this->m_data = new double[this->storedLines() * this->storedColumns()];

  std::copy(other.m_data,
            other.m_data + (other.storedLines() * other.storedColumns()),
            this->m_data);

  return *this;
}

// matrix::matrix(double k): m_data{1, 1, 1, 1, {k}} {}
// matrix::matrix(): m_data{0, 0, 0, 0} {}
// matrix::matrix(unsigned int l, unsigned int c, double* data): m_data{l, c, 1,
// 1, data} {} matrix::matrix(unsigned int l, unsigned int c,
// std::initializer_list<double> data, unsigned int bx, unsigned int by):
// m_data{l, c, bx, by, data} {} matrix::matrix(unsigned int l, unsigned int c,
// unsigned int bx, unsigned int by): m_data{l, c, bx, by} {}

void add(matrix *C, const matrix *const A, const matrix *const B, bool A_T,
         bool B_T) {
  assert(A->columns() == B->columns());
  assert(A->lines() == B->lines());

  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      unsigned int row_margin = C->lines() - i;
      unsigned int col_margin = C->columns() - j;

      int width = std::min(C->blockWidth(), col_margin);
      int heigth = std::min(C->blockHeight(), row_margin);

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          int Ci = i + y;
          int Cj = j + x;

          C->set(
              Ci, Cj,
              A->get((1 - A_T) * Ci + A_T * Cj, (1 - A_T) * Cj + A_T * Ci) +
                  B->get((1 - B_T) * Ci + B_T * Cj, (1 - B_T) * Cj + B_T * Ci));
        }
      }
    }
  }
}

void sub(matrix *C, const matrix *const A, const matrix *const B, bool A_T,
         bool B_T) {
  assert(A->columns() == B->columns());
  assert(A->lines() == B->lines());

  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      unsigned int row_margin = C->lines() - i;
      unsigned int col_margin = C->columns() - j;

      int width = std::min(C->blockWidth(), col_margin);
      int heigth = std::min(C->blockHeight(), row_margin);

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          int Ci = i + y;
          int Cj = j + x;

          C->set(
              Ci, Cj,
              A->get((1 - A_T) * Ci + A_T * Cj, (1 - A_T) * Cj + A_T * Ci) -
                  B->get((1 - B_T) * Ci + B_T * Cj, (1 - B_T) * Cj + B_T * Ci));
        }
      }
    }
  }
}

void div(matrix *C, const matrix *const A, const matrix *const B, bool A_T,
         bool B_T) {
  assert(A->columns() == B->columns());
  assert(A->lines() == B->lines());

  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {

      unsigned int row_margin = C->lines() - i;
      unsigned int col_margin = C->columns() - j;

      int width = std::min(C->blockWidth(), col_margin);
      int heigth = std::min(C->blockHeight(), row_margin);

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          int Ci = i + y;
          int Cj = j + x;

          C->set(
              Ci, Cj,
              A->get((1 - A_T) * Ci + A_T * Cj, (1 - A_T) * Cj + A_T * Ci) /
                  B->get((1 - B_T) * Ci + B_T * Cj, (1 - B_T) * Cj + B_T * Ci));
        }
      }
    }
  }
}

void LUdecompose(matrix *L, matrix *U, const matrix *const A) {
  int i = 0, j = 0, k = 0;

  int n = A->lines();

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j < i) {
        L->set(j, i, 0);
      } else {
        L->set(j, i, A->get(j, i));
        for (k = 0; k < i; k++) {
          L->set(j, i, L->get(j, i) - L->get(j, k) * U->get(k, i));
        }
      }
    }
    for (j = 0; j < n; j++) {
      if (j < i) {
        U->set(i, j, 0);
      } else if (j == i) {
        U->set(i, j, 1);
      } else {
        U->set(i, j, A->get(i, j) / L->get(i, i));

        for (k = 0; k < i; k++) {
          U->set(i, j,
                 U->get(i, j) - ((L->get(i, k) * U->get(k, j)) / L->get(i, i)));
        }
      }
    }
  }
}

void pivot(matrix *A, int i, int k) {
  for (unsigned int j = 0; j < A->columns(); j++) {
    double tmp = A->get(i, j);

    A->set(i, j, A->get(k, j));
    A->set(k, j, tmp);
  }
}

int LUPdecompose(matrix *A, matrix *P) {
  int i, j, k, imax;
  double maxA, absA;
  int N = A->lines();

  for (i = 0; i <= N; i++) {
    P->set(i, 0, i);
  }

  for (i = 0; i < N; i++) {
    maxA = 0.f;
    imax = i;

    for (k = i; k < N; k++) {
      if ((absA = fabs(A->get(k, i))) > maxA) {
        maxA = absA;
        imax = k;
      }
    }

    if (maxA < 0.00001)
      return 0;

    if (imax != i) {
      j = P->get(i, 0);
      P->set(i, 0, P->get(imax, 0));
      P->set(imax, 0, j);
      pivot(A, i, imax);
      P->set(N, 0, P->get(N, 0) + 1);
    }

    for (j = i + 1; j < N; j++) {
      A->set(j, i, A->get(j, i) / A->get(i, i));
      for (k = i + 1; k < N; k++) {
        A->set(j, k, A->get(j, k) - A->get(j, i) * A->get(i, k));
      }
    }
  }
  return 1;
}

int _LUPSolve(const matrix *const A, const matrix *const P,
              const matrix *const b, matrix *x) {
  int N = A->lines();

  for (int i = 0; i < N; i++) {
    x->set(i, 0, b->get(P->get(i, 0), 0));
    for (int k = 0; k < i; k++) {
      x->set(i, 0, x->get(i, 0) - A->get(i, k) * x->get(k, 0));
    }
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++) {
      x->set(i, 0, x->get(i, 0) - A->get(i, k) * x->get(k, 0));
    }
    x->set(i, 0, x->get(i, 0) / A->get(i, i));
  }

  return 1;
}

int LUPInvet(const matrix *const A, const matrix *const P, matrix *A_Inv) {
  int N = A->lines();

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      A_Inv->set(i, j, static_cast<int>(P->get(i, 0)) == j ? 1.f : 0.f);

      for (int k = 0; k < i; k++) {
        A_Inv->set(i, j, A_Inv->get(i, j) - A->get(i, k) * A_Inv->get(k, j));
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int k = i + 1; k < N; k++) {
        A_Inv->set(i, j, A_Inv->get(i, j) - A->get(i, k) * A_Inv->get(k, j));
      }
      A_Inv->set(i, j, A_Inv->get(i, j) / A->get(i, i));
    }
  }

  return 1;
}

double _LUPDeterminant(const matrix *const A, const matrix *const P) {
  double det = A->get(0, 0);
  int N = A->lines();

  for (int i = 1; i < N; i++) {
    det *= A->get(i, i);
  }

  return (static_cast<int>(P->get(N, 0)) - N) % 2 == 0 ? det : -det;
}

void div(matrix *C, const matrix *const A, const float alpha) {
  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      unsigned int row_margin = C->lines() - i;
      unsigned int col_margin = C->columns() - j;

      int width = std::min(C->blockWidth(), col_margin);
      int heigth = std::min(C->blockHeight(), row_margin);

      for (int y = 0; y < heigth; y++) {
        for (int x = 0; x < width; x++) {
          int Ci = i + y;
          int Cj = j + x;

          C->set(Ci, Cj, A->get(Ci, Cj) / alpha);
        }
      }
    }
  }
}

void swapLines(matrix *A, int i, int k) {
  for (unsigned int j = 0; j < A->columns(); j++) {
    double tmp = A->get(i, j);

    A->set(i, j, A->get(k, j));
    A->set(k, j, tmp);
  }
}

unsigned int findPivo(matrix *M, unsigned int l, unsigned int c) {
  unsigned int best = l;
  double val = fabs(M->get(l, c));

  for (unsigned int r = l + 1; r < M->lines(); r++) {
    double curr = fabs(M->get(r, c));

    if (curr > val) {
      best = r;
      val = curr;
    }
  }

  return best;
}

void toEchelonForm(matrix *M) {
  unsigned int rows = M->lines();
  unsigned int cols = M->columns();

  unsigned int r = 0;

  for (unsigned int c = 0; c < cols; c++) {

    unsigned int pivo = findPivo(M, r, c);

    if (pivo != r) {
      swapLines(M, r, pivo);
    }

    double b = M->get(r, c);

    if (b > std::numeric_limits<double>::epsilon()) {
      for (unsigned int j = c; j < cols; j++) {
        double a = M->get(r, j);
        M->set(r, j, a / b);
      }
    }

    if (r > 0) {
      for (unsigned int l = 0; l < r; l++) {
        double k = M->get(l, c);
        for (unsigned int j = c; j < cols; j++) {
          double a = M->get(l, j);
          double b = M->get(r, j);

          M->set(l, j, a - b * k);
        }
      }
    }

    if (r < rows - 1) {

      for (unsigned int l = r + 1; l < rows; l++) {
        double k = M->get(l, c);
        for (unsigned int j = c; j < cols; j++) {
          double a = M->get(l, j);
          double b = M->get(r, j);

          M->set(l, j, a - b * k);
        }
      }
    }

    r += 1;

    if (r == rows)
      break;
  }
}

void nullSpace(matrix *M, matrix &ns) {

  toEchelonForm(M);

  unsigned int lead = 0;
  unsigned int columns = M->columns();

  while (lead < std::min(M->lines(), M->columns()) && M->get(lead, lead) == 1) {
    lead++;
  }

  unsigned int rank = columns - lead;

  if (rank == 0) {
    std::vector<double> zeros(columns, 0);
    ns = matrix(1, columns, 1, 1, zeros.data());
    return;
  }

  std::vector<double> basis(rank * columns, 0);

  for (long long i = 1; i <= rank; i++)
    basis[i * columns - rank + (i - 1)] = 1;

  for (long long i = 0; i < rank; i++)
    for (long long j = 0; j < M->lines(); j++)
      basis[i * columns + j] -= M->get(j, lead + i);

  ns = matrix(rank, columns, 1, 1, basis.data());
}

void _transpose(matrix *C, const matrix *const A) {
  // code for GPU multiplication
  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      // loops above will become one block
      for (unsigned int y = 0; y < C->blockHeight(); y++) {
        for (unsigned int x = 0; x < C->blockWidth(); x++) {
          // loops above will become one thread
          C->set(i + y, j + x, A->get(j + x, i + y));
        }
      }
    }
  }
}

void mul(matrix *C, const matrix *const A, const matrix *const B) {
  assert(A->columns() == B->lines());
  // code for GPU multiplication
  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      // loops above will become one block
      for (unsigned int y = 0; y < C->blockHeight(); y++) {
        for (unsigned int x = 0; x < C->blockWidth(); x++) {
          // loops above will become one thread
          double acc = 0;

          int Ci = i + y, Cj = j + x;
          for (unsigned int k = 0; k < A->columns();
               k += std::min(A->blockWidth(),
                             B->blockHeight()) /*C->blockWidth()*/) {
            for (unsigned int q = 0;
                 q < std::min(A->blockWidth(), B->blockHeight()); q++) {
              unsigned int Ai = (i + y), Aj = (k + q);
              unsigned int Bi = (k + q), Bj = (j + x);

              bool Ai_bounds = Ai < A->storedLines();
              bool Aj_bounds = Aj < A->storedColumns();
              bool Bi_bounds = Bi < B->storedLines();
              bool Bj_bounds = Bj < B->storedColumns();

              Ai = (A->storedLines() - 1) * !Ai_bounds + Ai * Ai_bounds;
              Aj = (A->storedColumns() - 1) * !Aj_bounds + Aj * Aj_bounds;
              Bi = (B->storedLines() - 1) * !Bi_bounds + Bi * Bi_bounds;
              Bj = (B->storedColumns() - 1) * !Bj_bounds + Bj * Bj_bounds;

              acc += A->get(Ai, Aj) * B->get(Bi, Bj);
            }
          }
          C->set(Ci, Cj, acc);
        }
      }
    }
  }

  for (unsigned int j = 0; j < C->storedColumns() - C->columns(); j++) {
    for (unsigned int i = 0; i < C->lines(); i++) {
      int block_y = i / C->blockHeight();
      int block_x = (C->columns() + j) / C->blockWidth();

      int y = i - block_y * C->blockHeight();
      int x = (C->columns() + j) - block_x * C->blockWidth();

      C->data()[C->stride(block_y, block_x) + y * C->blockWidth() + x] = 0.f;
    }
  }

  for (unsigned int i = 0; i < C->storedLines() - C->lines(); i++) {
    for (unsigned int j = 0; j < C->columns(); j++) {
      int block_y = (C->lines() + i) / C->blockHeight();
      int block_x = j / C->blockWidth();

      int y = (C->lines() + i) - block_y * C->blockHeight();
      int x = j - block_x * C->blockWidth();

      C->data()[C->stride(block_y, block_x) + y * C->blockWidth() + x] = 0.f;
    }
  }

  for (unsigned int i = 0; i < std::min(C->storedLines() - C->lines(),
                                        C->storedColumns() - C->columns());
       i++) {
    int block_y = (C->lines() + i) / C->blockHeight();
    int block_x = (C->columns() + i) / C->blockWidth();

    int y = (C->lines() + i) - block_y * C->blockHeight();
    int x = (C->columns() + i) - block_x * C->blockWidth();

    C->data()[C->stride(block_y, block_x) + y * C->blockWidth() + x] = 1.f;
  }
}

void mul(matrix *C, const matrix *const A, const double alpha) {
  // code for GPU multiplication
  for (unsigned int i = 0; i < C->lines(); i += C->blockHeight()) {
    for (unsigned int j = 0; j < C->columns(); j += C->blockWidth()) {
      // loops above will become one block
      for (unsigned int y = 0; y < C->blockHeight(); y++) {
        for (unsigned int x = 0; x < C->blockWidth(); x++) {
          // loops above will become one thread
          int Ci = i + y, Cj = j + x;
          C->set(Ci, Cj, A->get(Ci, Cj) * alpha);
        }
      }
    }
  }

  for (unsigned int i = 0; i < std::min(C->storedLines() - C->lines(),
                                        C->storedColumns() - C->columns());
       i++) {
    int block_y = (C->lines() + i) / C->blockHeight();
    int block_x = (C->columns() + i) / C->blockWidth();

    int y = (C->lines() + i) - block_y * C->blockHeight();
    int x = (C->columns() + i) - block_x * C->blockWidth();

    C->data()[C->stride(block_y, block_x) + y * C->blockWidth() + x] = 1.f;
  }

  return;
}

// matrix& matrix::operator=(const matrix& other)
// {

// 	this->m_data = other.m_data;
// 	return *this;
// }

matrix matrix::operator+(const matrix &other) {
  matrix C(this->lines(), this->columns(), this->blockHeight(),
           other.blockWidth());

  add(&C, this, &other, false, false);

  return C;
}

matrix matrix::operator-(const matrix &other) {
  matrix C(this->lines(), this->columns(), this->blockHeight(),
           other.blockWidth());

  sub(&C, this, &other, false, false);

  return C;
}

matrix matrix::operator*(const matrix &other) {
  if (other.columns() == 1 && other.lines() == 1) {
    matrix D(this->lines(), this->columns(), this->blockWidth(),
             this->blockHeight());
    mul(&D, this, other.get(0, 0));
    return D;
  }

  if (this->columns() == 1 && this->lines() == 1) {
    matrix D(other.lines(), other.columns(), other.blockWidth(),
             other.blockHeight());
    mul(&D, &other, this->get(0, 0));
    return D;
  }

  matrix C(this->lines(), other.columns(), this->blockWidth(),
           other.blockHeight());

  mul(&C, this, &other);

  return C;
}

matrix matrix::operator/(const matrix &other) {
  if (other.columns() == 1 && other.lines() == 1) {
    matrix D(this->lines(), this->columns(), this->blockWidth(),
             this->blockHeight());
    div(&D, this, other.get(0, 0));
    return D;
  }

  if (this->columns() == 1 && this->lines() == 1) {
    matrix D(other.lines(), other.columns(), other.blockWidth(),
             other.blockHeight());
    div(&D, &other, this->get(0, 0));
    return D;
  }

  matrix C(this->lines(), other.columns(), this->blockWidth(),
           other.blockHeight());
  div(&C, this, &other, false, false);
  return C;
}

matrix matrix::operator/(const double other) {
  matrix C(this->lines(), this->columns(), this->blockWidth(),
           this->blockHeight());
  div(&C, this, other);
  return C;
}

matrix matrix::operator*(const double other) {
  matrix C(this->lines(), this->columns(), this->blockWidth(),
           this->blockHeight());
  mul(&C, this, other);
  return C;
}

matrix::MatrixLineGetter::MatrixLineGetter(matrix *parent, unsigned int line) {
  this->line = line;
  this->parent = parent;
}

double &matrix::MatrixLineGetter::operator[](const unsigned int idx) {
  return parent->get(this->line, idx);
}

matrix::MatrixLineGetter matrix::operator[](const unsigned int idx) {
  return matrix::MatrixLineGetter(this, idx);
}

matrix alg::transpose(matrix &other) {
  matrix C(other.columns(), other.lines(), other.blockHeight(),
           other.blockWidth());
  _transpose(&C, &other);
  return C;
}

matrix alg::transpose(matrix &&other) {
  matrix C(other.columns(), other.lines(), other.blockHeight(),
           other.blockWidth());
  _transpose(&C, &other);
  return C;
}

matrix alg::transpose(matrix *other) {
  matrix C(other->columns(), other->lines(), other->blockHeight(),
           other->blockWidth());
  _transpose(&C, other);
  return C;
}

// std::pair<matrix, matrix> alg::LUDecomposition(const matrix* const A)
// {
// 	matrix L(A->lines(), A->columns(), A->blockWidth(), A->blockHeight());
// 	matrix U(A->lines(), A->columns(), A->blockWidth(), A->blockHeight());

// 	LUdecompose(&L, &U, A);

// 	return { L, U };
// }

std::pair<matrix, matrix> alg::LUDecomposition(const matrix &A) {
  matrix L(A.lines(), A.columns(), A.blockWidth(), A.blockHeight());
  matrix U(A.lines(), A.columns(), A.blockWidth(), A.blockHeight());
  LUdecompose(&L, &U, &A);
  return {L, U};
}

std::pair<matrix, matrix> alg::LUPDecomposition(const matrix &A) {
  matrix R = matrix(A);
  matrix P(A.lines() + 1, 1);

  LUPdecompose(&R, &P);

  return {R, P};
}

matrix alg::LUPSolve(const matrix &A, const matrix &P, const matrix &b) {
  matrix x(A.lines(), 1);

  _LUPSolve(&A, &P, &b, &x);

  return x;
}

matrix alg::LUPInverse(const matrix &A, const matrix &P) {
  matrix Inv(A.lines(), A.columns(), A.blockWidth(), A.blockHeight());

  LUPInvet(&A, &P, &Inv);

  return Inv;
}

double alg::LUPDeterminant(const matrix &A, const matrix &P) {
  return _LUPDeterminant(&A, &P);
}

matrix alg::echelonForm(matrix matrix) {
  toEchelonForm(&matrix);
  return matrix;
}

matrix alg::nullspace(matrix m) {
  matrix null;

  nullSpace(&m, null);

  return null;
}

matrix alg::diag(matrix &&diag) {
  matrix a(diag.lines(), diag.lines());

  for (unsigned int i = 0; i < diag.lines(); i++)
    a[i][i] = diag[i][0];
  return a;
}

matrix alg::diag(matrix &diag) {
  matrix a(diag.lines(), diag.lines());

  for (unsigned int i = 0; i < diag.lines(); i++)
    a[i][i] = diag[i][0];
  return a;
}

matrix alg::diag(double *diag, unsigned int n) {
  matrix a(n, n);
  for (unsigned int i = 0; i < n; i++)
    a[i][i] = diag[i];
  return a;
}

matrix alg::diag(double *diag, unsigned int m, unsigned int n) {
  matrix a(m, n);

  for (unsigned int i = 0; i < std::min(n, m); i++)
    a[i][i] = diag[i];

  return a;
}

matrix alg::bidiag(double *diag, double *sdiag, unsigned int m,
                   unsigned int n) {
  matrix a(m, n);

  for (unsigned int i = 0; i < std::min(n, m); i++) {
    a[i][i] = diag[i];

    if (i < std::min(n, m) - 1) {
      a[i][i + 1] = sdiag[i + 1];
    }
  }

  return a;
}

matrix diag(matrix &&diag, unsigned int m, unsigned int n) {
  assert(diag.lines() == n);
  matrix a(m, n);
  for (unsigned int i = 0; i < std::min(n, m); i++)
    a[i][i] = diag[i][0];
  return a;
}

matrix diag(matrix &diag, unsigned int m, unsigned int n) {
  assert(diag.lines() == n);

  matrix a(m, n);
  for (unsigned int i = 0; i < n; i++)
    a[i][i] = diag[i][0];
  return a;
}

// void printMatrixWithMargin(matrix* A, unsigned precision, double eps)
// {
// 	for(int i=0;i<A->storedLines(); i++)
// 	{
// 		for(int j=0; j<A->storedColumns();j++)
// 		{
// 			double val = fabs(A->get(i,j)) > eps ? A->get(i,j) :
// 0.0;

// 			std::cout << std::fixed << std::setprecision(precision)
// << std::setw(precision+8) << val << " ";
// 		}

// 		std::cout << std::endl;
// 	}
// }

// void printMatrixWithMargin(matrix& A, unsigned precision, double eps)
// {
// 	for(int i=0;i<A.storedLines(); i++)
// 	{
// 		for(int j=0; j<A.storedColumns();j++)
// 		{
// 			double val = fabs(A.get(i,j)) > eps ? A.get(i,j) : 0.0;

// 			std::cout << std::fixed << std::setprecision(precision)
// << std::setw(precision+8) <<val << " ";
// 		}
// 		std::cout << std::endl;
// 	}
// }

void alg::printMatrix(matrix &A) {
  for (unsigned i = 0; i < A.lines(); i++) {
    for (unsigned j = 0; j < A.columns(); j++) {
      double val = A.get(i, j);
      printf("%f ", val);
    }

    printf("\n");
  }
}
std::string alg::matrixToString(matrix *m) {
  std::string r = "Mat";

  r += std::to_string(m->lines());
  r += "x";
  r += std::to_string(m->columns());
  r += "[";

  for (unsigned i = 0; i < m->lines(); i++) {
    r += "[";
    for (unsigned j = 0; j < m->columns(); j++) {
      r += to_string(number(m->get(i, j)));
      // r += std::to_string(m->get(i,j));
      if (j < m->columns() - 1) {
        r += ", ";
      }
    }
    r += "]";
    if (i < m->lines() - 1)
      r += ", ";
  }

  r += "]";

  return r;
}

// void printSubMatrix(matrix& A, int p, int q, int r, int s, unsigned
// precision, double eps)
// {
// 	for(int i=p;i<q; i++)
// 	{
// 		for(int j=r; j<s;j++)
// 		{
// 			double val = fabs(A.get(i,j)) > eps ? A.get(i,j) : 0.0;

// 			std::cout << std::scientific <<
// std::setprecision(precision)
// << std::setw(precision+8) <<val << ", ";
// 		}
// 		std::cout << std::endl;
// 	}
// }

void alg::printMatrix(matrix &&A) {
  for (unsigned i = 0; i < A.lines(); i++) {
    for (unsigned j = 0; j < A.columns(); j++) {
      double val = A.get(i, j);
      printf("%f ", val);
    }

    printf("\n");
  }
}

matrix alg::identity(unsigned int m, unsigned int n) {
  matrix I(m, n);

  for (unsigned int i = 0; i < std::min(n, m); i++) {
    I[i][i] = 1.0;
  }

  return I;
}

matrix alg::inverse(matrix &A) {
  std::pair<matrix, matrix> B = LUPDecomposition(A);
  return LUPInverse(B.first, B.second);
}

matrix alg::solve(matrix &A, matrix &b) {
  std::pair<matrix, matrix> B = LUPDecomposition(A);
  return LUPSolve(B.first, B.second, b);
}

matrix *matrix::add_ptr(matrix *A, matrix *B) {
  assert(A->lines() == B->lines() && A->columns() == B->columns());

  matrix *C =
      new matrix(A->lines(), B->columns(), A->blockWidth(), B->blockHeight());
  ::add(C, A, B, false, false);
  return C;
}

matrix *matrix::mul_ptr(matrix *A, matrix *B) {
  matrix *C =
      new matrix(A->lines(), B->columns(), A->blockWidth(), B->blockHeight());

  ::mul(C, A, B);

  return C;
}

matrix *matrix::mul_ptr(matrix *A, double B) {
  matrix *C =
      new matrix(A->lines(), A->columns(), A->blockWidth(), A->blockHeight());
  ::mul(C, A, B);
  return C;
}

matrix *matrix::div_ptr(matrix *A, double B) {
  matrix *C =
      new matrix(A->lines(), A->columns(), A->blockWidth(), A->blockHeight());
  ::div(C, A, B);
  return C;
}

matrix *matrix::inv_ptr(matrix *A) {
  matrix *t = matrix::copy(A);
  matrix *P = new matrix(A->lines() + 1, 1);

  LUPdecompose(t, P);

  matrix *Inv =
      new matrix(A->lines(), A->columns(), A->blockWidth(), A->blockHeight());

  LUPInvet(t, P, Inv);

  delete t;
  delete P;

  return Inv;
}

double matrix::det_ptr(matrix *A) {
  matrix *t = matrix::copy(A);
  matrix *P = new matrix(A->lines() + 1, 1);

  LUPdecompose(t, P);

  double d = _LUPDeterminant(t, P);

  delete t;
  delete P;

  return d;
}

matrix *matrix::transp_ptr(matrix *A) {
  matrix *C =
      new matrix(A->columns(), A->lines(), A->blockHeight(), A->blockWidth());
  _transpose(C, A);
  return C;
}

matrix *matrix::solve_ptr(matrix *A, matrix *b) {
  matrix *t = matrix::copy(A);
  matrix *P = new matrix(A->lines() + 1, 1);

  LUPdecompose(t, P);

  matrix *x = new matrix(A->lines(), 1);

  _LUPSolve(t, P, b, x);

  return x;
}

std::tuple<matrix *, matrix *, matrix *> matrix::svd_ptr(matrix *a) {
  matrix A = *a;

  std::tuple<matrix, matrix, matrix> SVD = ::svd(A);

  return std::make_tuple(new matrix(std::get<0>(SVD)),
                         new matrix(std::get<1>(SVD)),
                         new matrix(std::get<2>(SVD)));
}
