#include "Integer.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

Int::Int() {
	this->val = bint<30>::from<int>(0);
}

Int::Int(const Int &a) {
	this->val = a.val->copy();
}

Int::Int(Int &&a) {
  this->val = a.val;
  a.val = nullptr;
}

Int::Int(long int v) {
	this->val = bint<30>::from<long int>(v);
}

Int::Int(long long v) {
	this->val = bint<30>::from<long long>(v);
}

Int::Int(unsigned long long v) {
  this->val = bint<30>::from<unsigned long long>(v);
}

Int::Int(unsigned int v) {
	this->val = bint<30>::from<unsigned int>(v);
}

Int::Int(int v) {
	this->val = bint<30>::from<int>(v);
}

Int::Int(double b) {
	this->val = bint<30>::from(b);
}

Int::~Int() {
  if (this->val) delete this->val;
}

Int::Int(bint<30> *v) {
	this->val = v;
}

std::string Int::to_string() { return this->val->to_string(); }
