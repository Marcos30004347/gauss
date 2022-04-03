#ifndef ALG_HPP
#define ALG_HPP

#include "Integer.hpp"
#include "Matrix.hpp"

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace alg {

enum kind {
  //ERROR = (1 << 0),
  FACT = (1 << 1),
  POW = (1 << 2),
  MUL = (1 << 3),
  ADD = (1 << 4),
  SUB = (1 << 5),
  DIV = (1 << 6),
  ROOT = (1 << 7),
  INF = (1 << 8),
  UNDEF = (1 << 9),
  SYM = (1 << 10),
  INT = (1 << 11),
  FRAC = (1 << 12),
  FAIL = (1 << 13),
  FUNC = (1 << 14),

  LIST = (1 << 15),
  SET = (1 << 16),

  MAT = (1 << 18),

  // UTILS
  CONST = INT | FRAC,
  SUMMABLE = MUL | POW | SYM | INF | FUNC,
  MULTIPLICABLE = POW | SYM | FUNC | ADD | INF | UNDEF | FAIL,
  NON_CONSTANT = SYM | FUNC | INF | UNDEF | FAIL,
  TERMINAL = FAIL | UNDEF | FAIL | INF | SYM | INT,
  ORDERED = POW | DIV | ROOT | FUNC,
};

#define SORTED_BIT 1
#define REDUCED_BIT 2
#define EXPANDED_BIT 3

enum info {
  UNKNOWN = (0 << 0),
  SORTED = (1 << SORTED_BIT),
  REDUCED = (1 << REDUCED_BIT),
  EXPANDED = (1 << EXPANDED_BIT)
};

struct list;
struct set;

struct expr {
  enum kind kind_of = kind::UNDEF;

  int expr_info = info::UNKNOWN;
  int sort_kind = kind::UNDEF;

  union {
    char *expr_sym;
    list *expr_list;
    Int *expr_int;
    set *expr_set;
    matrix *expr_mat;
  };

  std::vector<expr> expr_childs;

  expr();
  expr(enum kind k);
  expr(expr &&other);
  expr(const expr &other);
  expr(Int v);
  expr(int v);
  expr(long int v);
  expr(long long v);
  expr(std::string v);

  expr(list &);
  expr(list &&);

  expr(set &);
  expr(set &&);

  expr(enum kind k, std::initializer_list<expr> &&);

  ~expr();

  std::vector<expr> operands();

  bool freeOf(expr &);
  bool freeOf(expr &&);

  enum kind kind() const;

  void insert(const expr &b);
  void insert(expr &&b);

  void insert(const expr &b, size_t idx);
  void insert(expr &&b, size_t idx);

  void remove(list &l);
  void remove(list &&l);

  void remove(size_t idx);
  void remove();

  bool match(expr *);

  bool operator==(const expr &);
  bool operator==(expr &&);

  bool operator!=(const expr &);
  bool operator!=(expr &&);

  expr &operator=(const expr &);
  expr &operator=(expr &&);
  expr &operator[](size_t idx);
  expr &operator[](Int idx);

  expr operator+(const expr &);
  expr operator+(expr &&);
  expr operator-(const expr &);
  expr operator-(expr &&);
  expr operator*(const expr &);
  expr operator*(expr &&);
  expr operator/(const expr &);
  expr operator/(expr &&);

  expr &operator+=(const expr &);
  expr &operator+=(expr &&);
  expr &operator-=(const expr &);
  expr &operator-=(expr &&);

  expr operator+();
  expr operator-();

  friend expr pow(const expr &a, const expr &b);
  friend expr pow(expr &&a, expr &&b);
  friend expr pow(expr &&a, const expr &b);
  friend expr pow(const expr &a, expr &&b);

  friend expr sqrt(const expr &a, expr b);
  friend expr sqrt(expr &&a, expr b);

  friend expr fact(const expr &a);
  friend expr fact(expr &&a);

  friend expr undefined();
  friend expr fail();
  friend expr inf();

  friend expr expand(expr &);
  friend expr expand(expr &&);

  inline Int value() const { return *expr_int; }
  inline std::string identifier() const { return expr_sym; }

  size_t size();

  friend expr first(expr &a);
  friend expr rest(expr &a);

  // list only methods
  friend expr append(const expr &, const expr &);
  friend expr append(const expr &, expr &&);

  friend expr join(const expr &, const expr &);
  friend expr join(const expr &, expr &&);

  // set only methods
  friend expr difference(const expr &, const expr &);
  friend expr difference(const expr &, expr &&);

  friend expr unification(const expr &, const expr &);
  friend expr unification(const expr &, expr &&);

  friend expr intersection(const expr &, const expr &);
  friend expr intersection(const expr &, expr &&);

  friend bool exists(const expr &, expr &);

  inline std::string funName() { return this->expr_sym; }
};

expr pow(const expr &a, const expr &b);
expr pow(expr &&a, expr &&b);
expr pow(expr &&a, const expr &b);
expr pow(const expr &a, expr &&b);

expr sqrt(const expr &a, expr b);
expr sqrt(expr &&a, expr b);

expr operator*(Int i, expr &&other);
expr operator*(Int i, expr &other);
expr operator+(Int i, expr &&other);
expr operator+(Int i, expr &other);
expr operator-(Int i, expr &&other);
expr operator-(Int i, expr &other);
expr operator/(Int i, expr &&other);
expr operator/(Int i, expr &other);
expr operator*(int i, expr &&other);
expr operator*(int i, expr &other);
expr operator+(int i, expr &&other);
expr operator+(int i, expr &other);
expr operator-(int i, expr &&other);
expr operator-(int i, expr &other);
expr operator/(int i, expr &&other);
expr operator/(int i, expr &other);
expr operator*(long i, expr &&other);
expr operator*(long i, expr &other);
expr operator+(long i, expr &&other);
expr operator+(long i, expr &other);
expr operator-(long i, expr &&other);
expr operator-(long i, expr &other);
expr operator/(long i, expr &&other);
expr operator/(long i, expr &other);
expr operator*(long long i, expr &&other);
expr operator*(long long i, expr &other);
expr operator+(long long i, expr &&other);
expr operator+(long long i, expr &other);
expr operator-(long long i, expr &&other);
expr operator-(long long i, expr &other);
expr operator/(long long i, expr &&other);
expr operator/(long long i, expr &other);

expr replace(expr &a, expr &b, expr &c);
expr replace(expr &a, expr &&b, expr &&c);
expr replace(expr &a, expr &b, expr &&c);
expr replace(expr &a, expr &&b, expr &c);

expr map(expr &u, expr (*f)(expr &));
expr map(expr &u, expr &v, expr (*f)(expr &, expr &));

void sort(expr *a);

expr create(kind kind);

expr create(kind kind, std::initializer_list<expr> &&);

expr func_call(const char *id, std::initializer_list<expr> &&);

list freeVariables(expr &a);

// terminals
// expr error(const char *message);

// const char *error_message(expr);
// const char *error_message(expr *);

expr symbol(const char *id);
expr integer(Int value);
expr fraction(Int num, Int den);
expr inf();
expr fail();
expr undefined();

expr lcm(expr &a, expr &b);
expr gcd(expr &a, expr &b);

expr &base(expr &u);
expr &degree(expr &u);

expr &numerator(expr &u);
expr &denominator(expr &u);

expr binomial(Int n, std::vector<Int> &ks);
expr binomial(Int n, std::vector<Int> &&ks);

expr diff(expr f, expr dx);

expr abs(expr x);

expr exp(expr x);

expr log(expr x, expr base);

expr ln(expr x);

expr mat(unsigned int l, unsigned int c);
expr mat(unsigned int l, unsigned int c, std::initializer_list<double> d);

expr identity_matrix(unsigned int l, unsigned int c);

expr inverse_matrix(expr A);

expr svd_matrix(expr A);

expr solve_linear_system(expr A, expr b);

expr transpose_matrix(expr A);

expr determinant_matrix(expr A);

expr mat_get(expr &a, unsigned i, unsigned j);

void mat_set(expr &a, unsigned i, unsigned j, expr n);

Int rows(expr a);

Int columns(expr a);

void sort_vec(std::vector<expr> &a, kind k, long int l, long int r);

inline int is(const expr *a, int k) { return a->kind_of & k; }

inline kind kind_of(const expr *expr) { return expr->kind_of; }

inline char *get_id(expr *expr) { return expr->expr_sym; }

inline Int get_val(expr *expr) { return Int(*expr->expr_int); }

inline const char *get_func_id(expr *expr) { return expr->expr_sym; }

std::string to_latex(expr *a, bool fraction = false,
                     unsigned long max_den = 1000);
std::string to_latex(expr a, bool fraction = false,
                     unsigned long max_den = 1000);

std::string to_string(expr *a);
std::string to_string(expr &a);
std::string to_string(expr &&a);

void expand(expr *a);
// void expr_print(expr *a, int tabs = 0);

struct list {
  std::vector<expr> members;

  list(std::initializer_list<expr> &&);
  list(std::vector<expr> &&);
  list(std::vector<expr> &);

  list(const list &) = default;
  list(list &&) = default;

  void append(list &&);
  void append(list &);
  void append(list *);

  void insert(expr &&, size_t);
  void insert(const expr &, size_t);
  void insert(expr &&);
  void insert(const expr &);

  void remove(list &M);
  void remove(list &&M);

  list rest(size_t from = 1);

  void join(list &a);
  void join(list &&a);

  void remove(size_t i);
  void remove();

  bool match(list *other);

  inline size_t size() const { return members.size(); }
  inline expr &operator[](size_t i) { return members[i]; }
  inline expr &operator[](Int i) { return members[i.longValue()]; }

  void sortMembers();

  list &operator=(const list &a) = default;
  list &operator=(list &&a) = default;

  bool operator==(list &o);
  bool operator==(list &&o);

  bool operator!=(list &o);
  bool operator!=(list &&o);

  friend expr fist(list &l);
  friend list rest(list &, size_t from);

  friend list append(list &, list &);
  friend list append(list &, list &&);
  friend list append(list &, list *);
  friend list insert(list &, const expr &);
  friend list insert(list &, expr &&);

  friend list remove(list &, list &);
  friend list remove(list &, list &&);
  friend list remove(list &, size_t i);

  friend list join(list &, list &);
  friend list join(list &, list &&);

  friend std::string to_string(list &);
  friend std::string to_string(list *);

  friend std::string to_latex(list *a, bool fraction, unsigned long max_den);
};

list rest(list &, size_t from = 1);

struct set {
  std::vector<expr> members;

  set(std::initializer_list<expr> &&);
  set(std::vector<expr> &&);
  set(std::vector<expr> &);

  set(set &) = default;
  set(set &&) = default;

  inline size_t size() const { return members.size(); }
  inline expr &operator[](size_t i) { return members[i]; }
  inline expr &operator[](Int i) { return members[i.longValue()]; }

  bool insert(const expr &a);
  bool insert(expr &&a);

  void remove(expr &a);
  void remove(expr &&a);

  void remove(size_t idx);

  friend expr fist(set &l);
  friend set rest(set &, size_t from);

  friend set difference(set &L, set &M);
  friend set unification(set &L, set &M);
  friend set intersection(set &L, set &M);
  friend bool set_exists(set &L, expr &e);
  friend int search(set &L, expr &e);

  void sort();
  void sort(enum kind);

  long match(set *other);

  bool operator==(set &o);
  bool operator==(set &&o);

  bool operator!=(set &o);
  bool operator!=(set &&o);

  friend std::string to_string(set &);
  friend std::string to_string(set *);

  friend std::string to_latex(set *a, bool fraction, unsigned long max_den);

  set &operator=(const set &) = default;
  set &operator=(set &&) = default;
};

set rest(set &, size_t from = 1);

void trim(set *);

inline size_t size_of(const expr *expr) {
  if (is(expr, kind::SET)) {
    return expr->expr_set->size();
  }

  if (is(expr, kind::LIST)) {
    return expr->expr_list->size();
  }

  return expr->expr_childs.size();
}

inline expr *operand(expr *const a, size_t i) {
  if (is(a, kind::SET)) {
    return &a->expr_set->members[i];
  }

  if (is(a, kind::LIST)) {
    return &a->expr_list->members[i];
  }

  return &a->expr_childs[i];
}

expr fromDouble(double v, Int max_den = 999999999999999999);

void decimalToFraction(double input, Int maxden, Int &n, Int &d);

double doubleFromExpr(expr a);

} // namespace alg

#endif
