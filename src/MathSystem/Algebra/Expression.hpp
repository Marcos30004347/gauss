#ifndef ALG_HPP
#define ALG_HPP

#include "Integer.hpp"
#include "String.hpp"

#include <cstddef>
#include <initializer_list>
#include <vector>


namespace alg {

enum kind {
  FACT = (1 << 0),
  POW = (1 << 1),
  MUL = (1 << 2),
  ADD = (1 << 3),
  SUB = (1 << 4),
  DIV = (1 << 5),
  SQRT = (1 << 6),
  INF = (1 << 7),
  UNDEF = (1 << 9),
  SYM = (1 << 10),
  INT = (1 << 11),
  FRAC = (1 << 12),
  FAIL = (1 << 13),
  FUNC = (1 << 14),

  LIST = (1 << 15),
  SET = (1 << 16),

  // UTILS
  CONST = INT | FRAC,
  SUMMABLE = MUL | POW | SYM | INF,
  MULTIPLICABLE = POW | SYM | ADD | INF | UNDEF | FAIL,
  NON_CONSTANT = SYM | FUNC | INF | UNDEF | FAIL,
  TERMINAL = FAIL | UNDEF | FAIL | INF | SYM | INT,
  ORDERED = POW | DIV | SQRT | FUNC,
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
  expr(string v);

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

  void remove(list& l);
  void remove(list&& l);

  void remove(size_t idx);
  void remove();

	bool match(expr*);

  bool operator==(const expr &);
  bool operator==(expr &&);

  bool operator!=(const expr &);
  bool operator!=(expr &&);

  expr &operator=(const expr&);
  expr &operator=(expr&&);
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

  friend expr sqrt(const expr &a);
  friend expr sqrt(expr &&a);

  friend expr fact(const expr &a);
  friend expr fact(expr &&a);

  friend expr undefined();
  friend expr fail();
  friend expr inf();

  friend expr reduce(expr &);
  friend expr reduce(expr &&);

  friend expr expand(expr &);
  friend expr expand(expr &&);

  inline Int value() const { return *expr_int; }
  inline string identifier() const { return string(expr_sym); }

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

  friend bool exists(const expr &, const expr &);
  friend bool exists(const expr &, expr &&);

	inline string funName() { return this->expr_sym; }
};

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

expr func_call(const char* id, std::initializer_list<expr>&&);

// terminals
expr symbol(const char *id);
expr integer(Int value);
expr fraction(Int num, Int den);
expr inf();
expr fail();
expr undefined();

expr lcm(expr& a, expr& b);
expr gcd(expr& a, expr& b);

expr& base(expr& u);
expr& degree(expr& u);

expr& numerator(expr& u);
expr& denominator(expr& u);

expr binomial(Int n, std::vector<Int>& ks);
expr binomial(Int n, std::vector<Int>&& ks);

expr diff(expr f, expr dx);

expr abs(expr x);

expr sinh(expr x);
expr cosh(expr x);
expr tanh(expr x);
expr exp(expr x);
expr cos(expr x);
expr sin(expr x);
expr tan(expr x);
expr csc(expr x);
expr cot(expr x);
expr log(expr x);
expr ln(expr x);
expr sec(expr x);
expr coth(expr x);
expr sech(expr x);
expr csch(expr x);
expr arccos(expr x);
expr arcsin(expr x);
expr arctan(expr x);
expr arccot(expr x);
expr arcsec(expr x);
expr arccsc(expr x);
expr arccosh(expr x);
expr arctanh(expr x);

void sort_vec(std::vector<expr> &a, kind k, long int l, long int r);

inline int is(const expr *a, int k) { return a->kind_of & k; }

inline kind kind_of(const expr *expr) { return expr->kind_of; }

inline char *get_id(expr *expr) { return expr->expr_sym; }

inline Int get_val(expr *expr) { return Int(*expr->expr_int); }

inline char *get_func_id(expr *expr) { return expr->expr_sym; }

string to_string(expr *a);
string to_string(expr &a);
string to_string(expr &&a);

int compare(expr *a, expr *b, kind ctx);

void reduce(expr *a);

void expand(expr *a);
void expr_print(expr *a, int tabs = 0);

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

  bool match(list* other);

	inline size_t size() const { return members.size(); }
  inline expr &operator[](size_t i) { return members[i]; }
  inline expr &operator[](Int i) { return members[i.longValue()]; }

	void sortMembers();

  list &operator=(const list &a) = default;
  list &operator=(list &&a) = default;

	bool operator==(list& o) ;
	bool operator==(list&& o) ;

	bool operator!=(list& o) ;
	bool operator!=(list&& o) ;

  friend expr fist(list &l);
  friend list rest(list &, size_t from);

  friend list append(list &, list &);
  friend list append(list &, list &&);
  friend list append(list &, list*);
	friend list insert(list&, const expr&);
	friend list insert(list&, expr&&);

  friend list remove(list &, list &);
  friend list remove(list &, list &&);
	friend list remove(list&, size_t i);

  friend list join(list &, list &);
  friend list join(list &, list &&);

	friend string to_string(list&);

	friend string to_string(list&);
	friend string to_string(list*);
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

	bool insert(const expr& a);
	bool insert(expr&& a);

	void remove(expr& a);
	void remove(expr&& a);

	void remove(size_t idx);

  friend expr fist(set &l);
  friend set rest(set &, size_t from);

  friend set difference(set &L, set &M);
  friend set unification(set &L, set &M);
  friend set intersection(set &L, set &M);
  friend bool exists(set &L, expr &e);
  friend int search(set &L, expr &e);

	void sort();
	void sort(enum kind);

  long match(set* other);

	bool operator==(set& o) ;
	bool operator==(set&& o) ;

	bool operator!=(set& o) ;
	bool operator!=(set&& o) ;

	friend string to_string(set&);
	friend string to_string(set*);

	set& operator=(const set&) = default;
	set& operator=(set&&) = default;
};

set rest(set &, size_t from = 1);

void trim(set *);

inline size_t size_of(const expr *expr) {
	if(is(expr, kind::SET)) {
		return expr->expr_set->size();
	}

	if(is(expr, kind::LIST)) {
		return expr->expr_list->size();
	}

	return expr->expr_childs.size();
}

inline expr *operand(expr *const a, size_t i) {
	if(is(a, kind::SET)) {
		return &a->expr_set->members[i];
	}

	if(is(a, kind::LIST)) {
		return &a->expr_list->members[i];
	}

  return &a->expr_childs[i];
}


void toFraction(double input, unsigned long long maxden, unsigned long long &n, unsigned long long &d);

} // namespace alg

#endif
