#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <string>

namespace alg {

class matrix {

public:
  matrix();

	matrix(unsigned int lines, unsigned int columns);

  matrix(const matrix &other);

	matrix(unsigned int lines, unsigned int columns, unsigned int block_width,
         unsigned int block_heigth, std::initializer_list<double> data);

	matrix(unsigned int lines, unsigned int columns, unsigned int block_width,
         unsigned int block_heigth, double *data);

	matrix(unsigned int lines, unsigned int columns, unsigned int block_width,
         unsigned int block_heigth);

  ~matrix();

  unsigned int lines() const;
  unsigned int columns() const;

  unsigned int blockHeight() const;
  unsigned int blockWidth() const;

  inline double &get(int i, int j) const {
    int block_y = i / this->m_block_heigth;
    int block_x = j / this->m_block_width;

    int y = i - block_y * this->m_block_heigth;
    int x = j - block_x * this->m_block_width;

    return this
        ->m_data[this->stride(block_y, block_x) + y * this->m_block_width + x];
  }

  inline void set(int i, int j, double val) const {
    int block_y = i / this->m_block_heigth;
    int block_x = j / this->m_block_width;

    int y = i - block_y * this->m_block_heigth;
    int x = j - block_x * this->m_block_width;

    this->m_data[this->stride(block_y, block_x) + y * this->m_block_width + x] =
        val;
  }
	inline double* data() { return this->m_data; }
	inline unsigned int storedLines() const { return this->m_stored_lines; }
	inline unsigned int storedColumns() const { return this->m_stored_column; }

  inline unsigned int stride(unsigned int block_line,
                             unsigned int block_column) const {
    unsigned int blocks_per_block_line =
        this->m_stored_column / this->m_block_width;
    unsigned int block_size = this->m_block_heigth * this->m_block_width;
    return block_line * blocks_per_block_line * block_size +
           block_column * block_size;
  }

  matrix &operator=(const matrix &other);

  // matrix(double k);
  // matrix(const matrix& other);
  // matrix(const matrix& other);

  // matrix();
  matrix(unsigned int l, unsigned int c, double *data);
  matrix(unsigned int l, unsigned int c, std::initializer_list<double> data,
         unsigned int bx = 1, unsigned int by = 1);
  // matrix(unsigned int l, unsigned int c, unsigned int bx = 1, unsigned int by
  // = 1);

  // matrix& operator=(const matrix& other);
  matrix operator+(const matrix &other);
  matrix operator-(const matrix &other);
  matrix operator*(const matrix &other);
  matrix operator*(const double other);
  matrix operator/(const matrix &other);
  matrix operator/(const double other);

  class MatrixLineGetter {
    friend class matrix;
    matrix *parent;
    unsigned int line;
    MatrixLineGetter(matrix *parent, unsigned int line);

  public:
    double &operator[](const unsigned int idx);
  };

  MatrixLineGetter operator[](const unsigned int idx);

	static matrix* copy(matrix*);

	static matrix* add_ptr(matrix*, matrix*);
	static matrix* div_ptr(matrix*, double i);
	static matrix* mul_ptr(matrix*, matrix*);
	static matrix* mul_ptr(matrix*, double i);
	static matrix* inv_ptr(matrix*);
	static double det_ptr(matrix*);
	static matrix* transp_ptr(matrix*);
	static matrix *solve_ptr(matrix *A, matrix *b);
	static std::tuple<matrix*, matrix*, matrix*> svd_ptr(matrix* a);

private:
  double *m_data;

  unsigned int m_lines;
  unsigned int m_columns;
  unsigned int m_stored_lines;
  unsigned int m_stored_column;
  unsigned int m_block_heigth;
  unsigned int m_block_width;

};


matrix diag(double *diag, unsigned int m, unsigned int n);
matrix diag(matrix &diag, unsigned int m, unsigned int n);
matrix diag(matrix &&diag, unsigned int m, unsigned int n);

matrix diag(double *diag, unsigned int n);

matrix diag(matrix &diag);
matrix diag(matrix &&diag);

matrix bidiag(double *diag, double *sdiag, unsigned int m, unsigned int n);

matrix echelonForm(matrix matrix);
matrix nullspace(matrix matrix);

matrix transpose(matrix &matrix);
matrix transpose(matrix &&matrix);
matrix transpose(matrix *matrix);

matrix identity(unsigned int m, unsigned int n);

void identity(matrix &I);
void identity(matrix *I);

// std::pair<matrix, matrix> LUDecomposition(const matrix *const A);
std::pair<matrix, matrix> LUDecomposition(const matrix &A);

std::pair<matrix, matrix> LUPDecomposition(const matrix &A);

matrix LUPSolve(const matrix &A, const matrix &P, const matrix &b);
matrix LUPInverse(const matrix &A, const matrix &P);

double LUPDeterminant(const matrix &A, const matrix &P);

matrix inverse(matrix& A);
matrix solve(matrix& A, matrix& b);

void printMatrix(matrix &A);
void printMatrix(matrix &&A);

std::string matrixToString(matrix* m);
std::string matrixToLatex(matrix* m, bool fractions = false, long max_den = 10000);

} // namespace algebra

#endif
