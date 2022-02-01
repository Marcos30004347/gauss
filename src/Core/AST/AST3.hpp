#include "Integer.hpp"
#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace ast_teste {
enum kind {
  FACT = (1 << 0),
  POW = (1 << 1),
  MUL = (1 << 2),
  ADD = (1 << 3),
  SUB = (1 << 4),
  DIV = (1 << 5),
  SQRT = (1 << 6),
  INF = (1 << 7),
  NEG_INF = (1 << 8),
  UNDEF = (1 << 9),
  INT = (1 << 10),
  FRAC = (1 << 11),
  SYM = (1 << 12),
  FAIL = (1 << 13),
  FUNC = (1 << 14),

  LIST = (1 << 15),
  SET = (1 << 16),

  // UTILS
  CONST = INT | FRAC,
  SUMMABLE = MUL | POW | SYM,
  MULTIPLICABLE = POW | SYM | ADD,
  NON_CONSTANT = SYM | FUNC,
  TERMINAL = FAIL | UNDEF | FAIL | INF | NEG_INF | SYM | INT,
  ORDERED = POW | DIV | SQRT | FUNC
};

struct list;
struct set;

struct ast {
  enum kind kind_of = UNDEF;

  union {
    char *ast_sym;
    list *ast_list;
    Int *ast_int;
    set *ast_set;
  };

  std::vector<ast> ast_childs;

  ast();
  ast(enum kind k);
  ast(ast &&other);
  ast(const ast &other);
  ast(Int v);
  ast(int v);
  ast(long int v);
  ast(long long v);
  ast(std::string v);

  ast(list &);
  ast(list &&);

  ast(set &);
  ast(set &&);

  ast(enum kind k, std::initializer_list<ast> &&);

  ~ast();

  bool freeOf(ast &);
  bool freeOf(ast &&);

  std::string to_string() const;

  enum kind kind() const;

  void insert(const ast &b);
  void insert(ast &&b);

  void insert(const ast &b, size_t idx);
  void insert(ast &&b, size_t idx);

  void remove(list& l);
  void remove(list&& l);

  void remove(size_t idx);
  void remove();

  bool operator==(const ast &);
  bool operator==(ast &&);

  bool operator!=(const ast &);
  bool operator!=(ast &&);

  ast &operator=(const ast &);
  ast &operator=(ast &&);
  ast &operator[](size_t idx);

  ast operator+(const ast &);
  ast operator+(ast &&);
  ast operator-(const ast &);
  ast operator-(ast &&);
  ast operator*(const ast &);
  ast operator*(ast &&);
  ast operator/(const ast &);
  ast operator/(ast &&);

  ast &operator+=(const ast &);
  ast &operator+=(ast &&);
  ast &operator-=(const ast &);
  ast &operator-=(ast &&);

  ast operator+();
  ast operator-();

  friend ast pow(const ast &a, const ast &b);
  friend ast pow(ast &&a, ast &&b);
  friend ast pow(ast &&a, const ast &b);
  friend ast pow(const ast &a, ast &&b);

  friend ast sqrt(const ast &a);
  friend ast sqrt(ast &&a);

  friend ast fact(const ast &a);
  friend ast fact(ast &&a);

  friend ast undefined();
  friend ast fail();
  friend ast inf();

  friend ast reduce(ast &);
  friend ast reduce(ast &&);

  friend ast expand(ast &);
  friend ast expand(ast &&);

  inline Int value() const { return *ast_int; }
  inline std::string identifier() const { return std::string(ast_sym); }

  inline size_t size() { return ast_childs.size(); }

  friend ast first(ast &a);
  friend ast rest(ast &a);

  // list only methods
  friend ast append(const ast &, const ast &);
  friend ast append(const ast &, ast &&);

  friend ast join(const ast &, const ast &);
  friend ast join(const ast &, ast &&);

  // set only methods
  friend ast difference(const ast &, const ast &);
  friend ast difference(const ast &, ast &&);

  friend ast unnification(const ast &, const ast &);
  friend ast unnification(const ast &, ast &&);

  friend ast intersection(const ast &, const ast &);
  friend ast intersection(const ast &, ast &&);

  friend int exists(const ast &, const ast &);
  friend int exists(const ast &, ast &&);
};

ast operator*(Int i, ast &&other);
ast operator*(Int i, ast &other);
ast operator+(Int i, ast &&other);
ast operator+(Int i, ast &other);
ast operator-(Int i, ast &&other);
ast operator-(Int i, ast &other);
ast operator/(Int i, ast &&other);
ast operator/(Int i, ast &other);
ast operator*(int i, ast &&other);
ast operator*(int i, ast &other);
ast operator+(int i, ast &&other);
ast operator+(int i, ast &other);
ast operator-(int i, ast &&other);
ast operator-(int i, ast &other);
ast operator/(int i, ast &&other);
ast operator/(int i, ast &other);
ast operator*(long i, ast &&other);
ast operator*(long i, ast &other);
ast operator+(long i, ast &&other);
ast operator+(long i, ast &other);
ast operator-(long i, ast &&other);
ast operator-(long i, ast &other);
ast operator/(long i, ast &&other);
ast operator/(long i, ast &other);
ast operator*(long long i, ast &&other);
ast operator*(long long i, ast &other);
ast operator+(long long i, ast &&other);
ast operator+(long long i, ast &other);
ast operator-(long long i, ast &&other);
ast operator-(long long i, ast &other);
ast operator/(long long i, ast &&other);
ast operator/(long long i, ast &other);

ast replace(ast &a, ast &b, ast &c);
ast replace(ast &a, ast &&b, ast &&c);
ast replace(ast &a, ast &b, ast &&c);
ast replace(ast &a, ast &&b, ast &c);

ast map(ast &u, ast (*f)(ast &));
ast map(ast &u, ast &v, ast (*f)(ast &, ast &));

void sort(ast *a);

ast create(kind kind);

ast create(kind kind, std::initializer_list<ast> &&);

ast symbol(const char *id);

ast integer(Int value);

ast fraction(Int num, Int den);

inline ast *operand(ast *const a, size_t i) { return &a->ast_childs[i]; }

inline int is(const ast *a, int k) { return a->kind_of & k; }

inline size_t size_of(const ast *ast) { return ast->ast_childs.size(); }

inline kind kind_of(const ast *ast) { return ast->kind_of; }

inline char *get_id(ast *ast) { return ast->ast_sym; }

inline Int get_val(ast *ast) { return Int(*ast->ast_int); }

inline char *get_func_id(ast *ast) { return ast->ast_sym; }

std::string to_string(ast *a);

int compare(ast *a, ast *b, kind ctx);

void reduce(ast *a);

void expand(ast *a);

void ast_print(ast *a, int tabs = 0);

struct list {
  std::vector<ast> members;

  list(std::initializer_list<ast> &&);
  list(std::vector<ast> &&);
  list(std::vector<ast> &);

  list(const list &) = default;
  list(list &&) = default;

  void append(ast &&);
  void append(const ast &);

  void remove(list &M);
  void remove(list &&M);

  list rest(size_t from = 1);

  void join(list &a);
  void join(list &&a);

	void remove(size_t i);
	void remove();

  bool match(list* other);

	inline size_t size() const { return members.size(); }
  inline ast &operator[](size_t i) { return members[i]; }

  list &operator=(const list &a) = default;
  list &operator=(list &&a) = default;

	bool operator==(list& o) ;
	bool operator==(list&& o) ;

	bool operator!=(list& o) ;
	bool operator!=(list&& o) ;

  friend ast fist(list &l);
  friend list rest(list &, size_t from);

  friend list append(list &, const ast &);
  friend list append(list &, ast &&);

  friend list remove(list &, list &);
  friend list remove(list &, list &&);
	friend list remove(list&, size_t i);

  friend list join(list &, list &);
  friend list join(list &, list &&);

	friend std::string to_string(list&);

	friend std::string to_string(list&);
	friend std::string to_string(list*);
};

list rest(list &, size_t from = 1);

struct set {
  std::vector<ast> members;

  set(std::initializer_list<ast> &&);
  set(std::vector<ast> &&);
  set(std::vector<ast> &);

  set(set &) = default;
  set(set &&) = default;

  inline size_t size() const { return members.size(); }
  inline ast &operator[](size_t i) { return members[i]; }

  friend ast fist(set &l);
  friend set rest(set &, size_t from);

  friend set difference(set &L, set &M);
  friend set unnification(set &L, set &M);
  friend set intersection(set &L, set &M);
  friend int exists(set &L, ast &e);

  long match(set* other);

	bool operator==(set& o) ;
	bool operator==(set&& o) ;

	bool operator!=(set& o) ;
	bool operator!=(set&& o) ;

	friend std::string to_string(set&);
	friend std::string to_string(set*);
};

set rest(set &, size_t from = 1);

void trim(set *);

} // namespace ast_teste
